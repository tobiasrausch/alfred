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
    uint32_t targetham;
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

  inline uint32_t
  hamming(std::string const& s1, std::string const& s2) {
    uint32_t ham = 0;
    for(uint32_t i = 0; i < s1.size(); ++i) {
      if (s1[i] != s2[i]) ++ham;
    }
    return ham;
  }


  inline uint32_t
  medianHammingDist(std::vector<uint32_t> const& ham, uint32_t const nseq, uint32_t const idx) {
    std::vector<uint32_t> dist;
    for(uint32_t idx1 = 0; idx1 < idx; ++idx1) dist.push_back(ham[idx1 * nseq + idx]);
    for(uint32_t idx1 = idx + 1; idx1 < nseq; ++idx1) dist.push_back(ham[idx * nseq + idx1]);
    std::sort(dist.begin(), dist.end());
    return dist[(int) dist.size()/2];
  }

  template<typename TConfig>
  inline void
  writebarcodes(TConfig const& c, std::vector<bool> const& bl) {
    std::ofstream ofile(c.outfile.string().c_str());
    faidx_t* fai = fai_load(c.infile.string().c_str());
    int32_t nseq = faidx_nseq(fai);
    for(int32_t idx=0; idx < nseq; ++idx) {
      if (bl[idx]) continue;
      std::string tname(faidx_iseq(fai, idx));
      int32_t seqlen = -1;
      char* seq = faidx_fetch_seq(fai, tname.c_str(), 0, faidx_seq_len(fai, tname.c_str()), &seqlen);
      ofile << ">" << tname << std::endl;
      ofile << seq << std::endl;
      free(seq);
    }
    fai_destroy(fai);
    ofile.close();
  }
  
  
  template<typename TConfig>
  inline int32_t
  runBarcode(TConfig& c) {

#ifdef PROFILE
    ProfilerStart("barcode.prof");
#endif

    // Rebuild index file
    boost::filesystem::remove(c.infile.string() + ".fai");
    
    // Read barcodes
    faidx_t* fai = fai_load(c.infile.string().c_str());
    int32_t nseq = faidx_nseq(fai);
    std::vector<uint32_t> ham(nseq*nseq, std::numeric_limits<uint32_t>::max());
    std::vector<double> entr(nseq, 0);
    for(int32_t idx1=0; idx1 < nseq; ++idx1) {
      std::string tname1(faidx_iseq(fai, idx1));
      int32_t seqlen1 = -1;
      char* seq1 = faidx_fetch_seq(fai, tname1.c_str(), 0, faidx_seq_len(fai, tname1.c_str()), &seqlen1);
      std::string s1 = boost::to_upper_copy(std::string(seq1));
      if (!isDNA(s1)) {
	std::cerr << "Nucleotide is not [A, C, G, T]!" << std::endl;
	std::cerr << ">" << tname1 << std::endl;
	std::cerr << s1 << std::endl;
	return 1;
      }
      entr[idx1] = entropy(s1);
      for(int32_t idx2 = idx1 + 1; idx2 < nseq; ++idx2) {
	std::string tname2(faidx_iseq(fai, idx2));
	int32_t seqlen2 = -1;
	char* seq2 = faidx_fetch_seq(fai, tname2.c_str(), 0, faidx_seq_len(fai, tname2.c_str()), &seqlen2);
	std::string s2 = boost::to_upper_copy(std::string(seq2));
	if (!isDNA(s2)) {
	  std::cerr << "Nucleotide is not [A, C, G, T]!" << std::endl;
	  std::cerr << ">" << tname2 << std::endl;
	  std::cerr << s2 << std::endl;
	  return 1;
	}
	if (seqlen1 != seqlen2) {
	  std::cerr << "Barcodes have different length!" << std::endl;
	  std::cerr << ">" << tname1 << std::endl;
	  std::cerr << seq1 << std::endl;
	  std::cerr << ">" << tname2 << std::endl;
	  std::cerr << seq2 << std::endl;
	  return 1;
	}
	//std::cerr << idx1 << ',' << idx2 << std::endl;
	//std::cerr << s1 << std::endl;
	//std::cerr << s2 << std::endl;
	ham[idx1 * nseq + idx2] = hamming(s1, s2);
	free(seq2);
      }
      free(seq1);
    }
    fai_destroy(fai);

    // Blacklist barcodes
    std::vector<bool> bl(nseq, false);
    bool rewind = true;
    while (rewind) {
      for(int32_t idx1=0; (idx1 < nseq) && (rewind); ++idx1) {
	if (bl[idx1]) continue;
	for(int32_t idx2 = idx1 + 1; (idx2 < nseq) && (rewind); ++idx2) {
	  if (bl[idx2]) continue;
	  if (ham[idx1 * nseq + idx2] < c.targetham) {
	    rewind = false;
	    if (entr[idx1] < entr[idx2]) bl[idx1] = true;
	    else bl[idx2] = true;
	  }
	}
      }
      if (!rewind) rewind = true;
      else rewind = false;
    }

    // Output
    writebarcodes(c, bl);
    
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
    if ((vm.count("help")) || (!vm.count("input-file"))) {
      std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] -o out.fasta <input.fasta>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    } 
    
    // Show cmd
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] ";
    std::cout << "alfred ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return runBarcode(c);
  }

}

#endif

