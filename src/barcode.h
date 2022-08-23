#ifndef BARCODE_H
#define BARCODE_H

#include <iostream>
#include <vector>
#include <fstream>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#ifdef PROFILE
#include "gperftools/profiler.h"
#endif

namespace bamstats {

  struct BarcodeConfig {
    bool hasInput;
    uint32_t targetham;
    uint32_t barlen;
    uint32_t enumall;
    float reqent;
    boost::filesystem::path outfile;
    boost::filesystem::path infile;
  };


  inline bool
  isDNA(std::string const& str) {
    for(uint32_t i = 0; i < str.size(); ++i) {
      if ((str[i] != 'A') && (str[i] != 'C') && (str[i] != 'G') && (str[i] != 'T')) return false;
    }
    return true;
  }

  template<typename TConfig>
  inline uint32_t
  hamming(TConfig const& c, std::string const& s1, std::string const& s2) {
    uint32_t ham = 0;
    for(uint32_t i = 0; ((i < s1.size()) && (ham < c.targetham)); ++i) {
      if (s1[i] != s2[i]) ++ham;
    }
    return ham;
  }


  template<typename TConfig>
  inline void
  loadbarcodes(TConfig const& c, std::vector<std::string>& barcodes, std::vector<std::string>& names) {
    faidx_t* fai = fai_load(c.infile.string().c_str());
    int32_t nseq = faidx_nseq(fai);
    names.resize(nseq);
    barcodes.resize(nseq);
    for(int32_t idx=0; idx < nseq; ++idx) {
      names[idx] = faidx_iseq(fai, idx);
      int32_t seqlen = -1;
      char* seq = faidx_fetch_seq(fai, names[idx].c_str(), 0, faidx_seq_len(fai, names[idx].c_str()), &seqlen);
      barcodes[idx] = boost::to_upper_copy(std::string(seq));
      free(seq);
    }
    fai_destroy(fai);
  }

  template<typename TConfig>
  inline void
  writebarcodes(TConfig const& c, std::vector<std::string> const& barcodes, std::vector<std::string> const& names, std::vector<bool> const& bl) {
    std::ofstream ofile(c.outfile.string().c_str());
    uint32_t count = 0;
    for(uint32_t idx=0; idx < barcodes.size(); ++idx) {
      if (bl[idx]) continue;
      if (names.empty()) ofile << ">Barcode" << count << std::endl;
      else ofile << ">" << names[idx] << std::endl;
      ofile << barcodes[idx] << std::endl;
      ++count;
    }
    ofile.close();
  }

  // Extend barcodes
  template<typename TConfig>
  inline void
  extendbarcodes(TConfig const& c, std::vector<std::string> const& barcodes) {
    std::ofstream ofile(c.outfile.string().c_str());
    uint32_t count = 0;
    for(uint32_t i1=0; i1 < barcodes.size(); ++i1) {
      for(uint32_t i2=0; i2 < barcodes.size(); ++i2) {
	ofile << ">Barcode" << count << std::endl;
	ofile << barcodes[i1] + barcodes[i2] << std::endl;
	++count;
      }
    }
    ofile.close();
  }
  
  template<typename TConfig>
  inline int32_t
  runBarcode(TConfig& c) {

#ifdef PROFILE
    ProfilerStart("barcode.prof");
#endif

    // Load or generate barcodes
    std::vector<char> alphabet = {'A', 'C', 'G', 'T'};
    std::vector<std::string> barcodes;
    std::vector<std::string> names;
    if (c.hasInput) loadbarcodes(c, barcodes, names);
    else {
      barcodes.push_back("A");
      barcodes.push_back("C");
      barcodes.push_back("G");
      barcodes.push_back("T");
      for(uint32_t i = 1; i < c.barlen; ++i) {
	uint32_t elst = barcodes.size();
	for(uint32_t k = 0; k < elst; ++k) {
	  std::string prefix = barcodes[k];
	  barcodes[k] += alphabet[i % 4];
	  for(uint32_t m = 1; m < 4; ++m) {
	    if (i <= c.enumall) barcodes.push_back(prefix + alphabet[(i + m) % 4]);
	    else if (entropy(prefix + alphabet[(i + m) % 4]) >= c.reqent) {
	      barcodes.push_back(prefix + alphabet[(i + m) % 4]);
	      break;
	    }
	  }
	}
	//std::cerr << i << ',' << barcodes.size() << std::endl;
      }
    }
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Loaded/Generated " << barcodes.size() << " barcode candidates" << std::endl;

    // Calculate entropy
    std::vector<float> entr(barcodes.size(), 0);
    std::vector<bool> bl(barcodes.size(), false);
    uint32_t blcount = 0;
    for(uint32_t i = 0; i < barcodes.size(); ++i) {
      entr[i] = entropy(barcodes[i]);
      if (entr[i] < c.reqent) {
	++blcount;
	bl[i] = true;
      }
    }
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Removed " << blcount << " barcodes because of low entropy" << std::endl;

    // Filter by hamming distance
    bool rewind = true;
    uint32_t iter = 0;
    std::vector<bool> validrow(barcodes.size(), false);
    while (rewind) {
      for(uint32_t i1=0; i1 < barcodes.size(); ++i1) {
	if ((validrow[i1]) || (bl[i1])) continue;
	validrow[i1] = true;
	for(uint32_t i2 = i1 + 1; i2 < barcodes.size(); ++i2) {
	  if (bl[i2]) continue;
	  if (hamming(c, barcodes[i1], barcodes[i2]) < c.targetham) {
	    rewind = false;
	    if (entr[i1] < entr[i2]) bl[i1] = true;
	    else bl[i2] = true;
	    ++blcount;
	    validrow[i1] = false;
	    break;
	  }
	}
      }
      if (!rewind) rewind = true;
      else rewind = false;
      ++iter;

      std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] Iteration " << iter << ": Removed " << blcount << " because of low entropy or hamming distance" << std::endl;
    }
    
    // Output
    writebarcodes(c, barcodes, names, bl);
    
#ifdef PROFILE
    ProfilerStop();
#endif

    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] ";
    std::cout << "Done." << std::endl;
    
    return 0;
  }
  

  int barcode(int argc, char **argv) {
    BarcodeConfig c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("target,t", boost::program_options::value<uint32_t>(&c.targetham)->default_value(3), "min. hamming distance")
      ("barlen,l", boost::program_options::value<uint32_t>(&c.barlen)->default_value(6), "barcode length")
      ("enumall,a", boost::program_options::value<uint32_t>(&c.enumall)->default_value(8), "enumerate all possible barcodes until this length")
      ("entropy,e", boost::program_options::value<float>(&c.reqent)->default_value(1.9), "min. barcode entropy")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("bar.fa"), "output FASTA file")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.infile), "input FASTA file")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if (vm.count("help")) {
      std::cout << "Usage: alfred " << argv[0] << " [OPTIONS]" << std::endl;
      std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] <barcodes.fasta>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    }

    if (vm.count("input-file")) c.hasInput = true;
    else c.hasInput = false;
    
    // Show cmd
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] ";
    std::cout << "alfred ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return runBarcode(c);
  }

}

#endif

