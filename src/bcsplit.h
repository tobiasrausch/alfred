#ifndef BCSPLIT_H
#define BCSPLIT_H

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

#include <htslib/sam.h>
#include <htslib/faidx.h>

#ifdef PROFILE
#include "gperftools/profiler.h"
#endif

namespace bamstats {

  struct BcsplitConfig {
    uint32_t hamming;
    uint32_t ncount;
    uint32_t umilen;
    uint32_t barlen;
    std::string outprefix;
    std::string pattern;
    boost::filesystem::path barcodes;
    boost::filesystem::path fastqfile;
    boost::filesystem::path idxfile;
  };

  template<typename TConfig>
  inline bool
  hammingExt(TConfig const& c, std::map<std::string, uint32_t>& barcode_fmap) {
    typedef typename std::map<std::string, uint32_t> TBarFileMap;

    // Debug
    //for(typename TBarFileMap::const_iterator it = barcode_fmap.begin(); it != barcode_fmap.end(); ++it) std::cerr << it->first << '\t' << it->second << std::endl;

    // Alphabet
    typedef std::set<char> TAlphabet;
    char tmp[] = {'A', 'C', 'G', 'T'};
    TAlphabet alphabet(tmp, tmp + sizeof(tmp) / sizeof(tmp[0]));
    
    // Extend barcode set
    for(uint32_t i = 0; i < c.hamming; ++i) {
      // Copy original barcodes
      TBarFileMap copy;
      for(typename TBarFileMap::const_iterator it = barcode_fmap.begin(); it != barcode_fmap.end(); ++it) copy.insert(std::make_pair(it->first, it->second));

      // Mutate each barcode
      for(typename TBarFileMap::const_iterator it = copy.begin(); it != copy.end(); ++it) {
	for(uint32_t k = 0; k < it->first.size(); ++k) {
	  // Mutate each position
	  for(typename TAlphabet::const_iterator alit = alphabet.begin(); alit != alphabet.end(); ++alit) {
	    std::string bar = it->first;
	    bar[k] = *alit;
	    if (barcode_fmap.find(bar) != barcode_fmap.end()) {
	      if (barcode_fmap[bar] != it->second) {
		std::cerr << "Hamming extension causes non-unique barcodes!" << std::endl;
		std::cerr << bar << '\t' << barcode_fmap[bar] << '\t' << it->second << std::endl;
		return false;
	      }
	    } else barcode_fmap.insert(std::make_pair(bar, it->second));
	  }
	}
      }
    }

    // Debug
    //for(typename TBarFileMap::const_iterator it = barcode_fmap.begin(); it != barcode_fmap.end(); ++it) std::cerr << it->first << '\t' << it->second << std::endl;

    // All fine
    return true;
  }


  template<typename TConfig>
  inline bool
  allowNs(TConfig const& c, std::map<std::string, uint32_t>& barcode_fmap) {
    typedef typename std::map<std::string, uint32_t> TBarFileMap;

    // Debug
    //for(typename TBarFileMap::const_iterator it = barcode_fmap.begin(); it != barcode_fmap.end(); ++it) std::cerr << it->first << '\t' << it->second << std::endl;

    // Extend barcode set
    for(uint32_t i = 0; i < c.ncount; ++i) {
      // Copy original barcodes
      TBarFileMap copy;
      for(typename TBarFileMap::const_iterator it = barcode_fmap.begin(); it != barcode_fmap.end(); ++it) copy.insert(std::make_pair(it->first, it->second));

      // Mutate each barcode
      for(typename TBarFileMap::const_iterator it = copy.begin(); it != copy.end(); ++it) {
	for(uint32_t k = 0; k < it->first.size(); ++k) {
	  std::string bar = it->first;
	  bar[k] = 'N';
	  if (barcode_fmap.find(bar) != barcode_fmap.end()) {
	    if (barcode_fmap[bar] != it->second) {
	      std::cerr << "N extension causes non-unique barcodes!" << std::endl;
	      std::cerr << bar << '\t' << barcode_fmap[bar] << '\t' << it->second << std::endl;
	      return false;
	    }
	  } else barcode_fmap.insert(std::make_pair(bar, it->second));
	}
      }
    }

    // Debug
    //for(typename TBarFileMap::const_iterator it = barcode_fmap.begin(); it != barcode_fmap.end(); ++it) std::cerr << it->first << '\t' << it->second << std::endl;

    // All fine
    return true;
  }

  template<typename TConfig>
  inline bool
  loadBarcodes(TConfig const& c, std::map<std::string, uint32_t>& barcode_fmap, std::vector<std::string>& ids) {
    std::set<std::string> idset;
    
    std::ifstream barFile;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    if (is_gz(c.barcodes)) {
      barFile.open(c.barcodes.string().c_str(), std::ios_base::in | std::ios_base::binary);
      dataIn.push(boost::iostreams::gzip_decompressor(), 16*1024);
    } else {
      barFile.open(c.barcodes.string().c_str(), std::ios_base::in);
    }
    dataIn.push(barFile);

    // Read barcodes
    std::istream instream(&dataIn);
    std::string gline;
    while(std::getline(instream, gline)) {
      typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
      boost::char_separator<char> sep("\t, ");
      Tokenizer tokens(gline, sep);
      Tokenizer::iterator tokIter = tokens.begin();
      if (tokIter != tokens.end()) {
	if (*tokIter == "#") continue;
	std::string barcode = *tokIter;
	if (c.barlen != barcode.size()) {
	  std::cerr << "Barcode length does not match pattern!" << std::endl;
	  return false;
	}
	if (barcode_fmap.find(barcode) != barcode_fmap.end()) {
	  std::cerr << "Duplicate barcodes: " << barcode << std::endl;
	  return false;
	}
	++tokIter;
	if (tokIter != tokens.end()) {
	  std::string idname = *tokIter;
	  if (idset.find(idname) != idset.end()) {
	    std::cerr << "Duplicate barcode IDs: " << idname << std::endl;
	    return false;
	  }
	  barcode_fmap.insert(std::make_pair(barcode, ids.size()));
	  ids.push_back(idname);
	  idset.insert(idname);
	}
      }
    }
    dataIn.pop();
    if (is_gz(c.barcodes)) dataIn.pop();
    barFile.close();

    // All fine
    return true;
  }    
  
  template<typename TConfig>
  inline int32_t
  runBcsplit(TConfig& c) {

#ifdef PROFILE
    ProfilerStart("barcode.prof");
#endif

    // Get barcode and umi pattern
    std::vector<bool> barpos(c.pattern.size(), false);
    std::vector<bool> umipos(c.pattern.size(), false);
    c.barlen = 0;
    c.umilen = 0;
    for(uint32_t i = 0; i < c.pattern.size(); ++i) {
      if (c.pattern[i] == 'B') {
	barpos[i] = true;
	++c.barlen;
      } else if (c.pattern[i] == 'U') {
	umipos[i] = true;
	++c.umilen;
      }
    }

    // Open barcode file
    typedef std::map<std::string, uint32_t> TBarFileMap;
    TBarFileMap barcode_fmap;
    std::vector<std::string> ids;
    if (!loadBarcodes(c, barcode_fmap, ids)) return -1;

    // Extend by hamming distance
    if (!hammingExt(c, barcode_fmap)) return -1;

    // Extend using Ns
    if (!allowNs(c, barcode_fmap)) return -1;
    
    // Open output files
    std::vector<boost::iostreams::filtering_ostream> dataOut(ids.size());
    for(uint32_t i = 0; i < dataOut.size(); ++i) {
      dataOut[i].push(boost::iostreams::gzip_compressor());
      std::string filename = c.outprefix + "." + ids[i] + ".fq.gz";
      dataOut[i].push(boost::iostreams::file_sink(filename, std::ios_base::out | std::ios_base::binary));
    }
    
    // Open index file
    std::ifstream idxfile;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIdx;
    if (is_gz(c.idxfile)) {
      idxfile.open(c.idxfile.string().c_str(), std::ios_base::in | std::ios_base::binary);
      dataIdx.push(boost::iostreams::gzip_decompressor(), 16*1024);
    } else idxfile.open(c.idxfile.string().c_str(), std::ios_base::in);
    dataIdx.push(idxfile);
    std::istream idxstream(&dataIdx);
    std::string idxline;
    std::string idxseq;
    
    // Open FASTQ
    std::ifstream fqfile;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    if (is_gz(c.fastqfile)) {
      fqfile.open(c.fastqfile.string().c_str(), std::ios_base::in | std::ios_base::binary);
      dataIn.push(boost::iostreams::gzip_decompressor(), 16*1024);
    } else fqfile.open(c.fastqfile.string().c_str(), std::ios_base::in);
    dataIn.push(fqfile);
    std::istream instream(&dataIn);
    std::string gline;
    uint64_t lnum = 0;
    std::string qname;
    std::string seq;

    // Parse jointly
    std::string barloc(c.barlen, 'N');
    std::string umiloc(c.umilen, 'N');
    uint64_t assigned = 0;
    uint64_t unassigned = 0;
    while(std::getline(instream, gline)) {
      if (!std::getline(idxstream, idxline)) {
	std::cerr << "Index FASTQ and read FASTQ differ in length!" << std::endl;
	return -1;
      }
      if (lnum % 4 == 0) qname = gline;
      else if (lnum % 4 == 1) {
	seq = gline;
	idxseq = idxline;
      }
      else if (lnum % 4 == 3) {
	if (c.pattern.size() != idxseq.size()) {
	  std::cerr << "Pattern size is not equal to index read size!" << std::endl;
	  return -1;
	}
	uint32_t bari = 0;
	uint32_t umii = 0;
	for(uint32_t i = 0; i < idxseq.size(); ++i) {
	  if (barpos[i]) barloc[bari++] = idxseq[i];
	  else if (umipos[i]) umiloc[umii++] = idxseq[i];
	}
	TBarFileMap::const_iterator it = barcode_fmap.find(barloc);
	if (it != barcode_fmap.end()) {
	  ++assigned;
	  dataOut[it->second] << qname << " umi:" << umiloc << " barcode:" << barloc << std::endl;
	  dataOut[it->second] << seq << std::endl;
	  dataOut[it->second] << "+" << std::endl;
	  dataOut[it->second] << gline << std::endl;
	} else ++unassigned;
      }
      ++lnum;
    }
    // Clean-up
    dataIn.pop();
    if (is_gz(c.fastqfile)) dataIn.pop();
    fqfile.close();
    dataIdx.pop();
    if (is_gz(c.idxfile)) dataIdx.pop();
    idxfile.close();	
    for(uint32_t i = 0; i < dataOut.size(); ++i) {
      dataOut[i].pop();
      dataOut[i].pop();
    }

    // Statistics
    double percass = 100 * (double) assigned / ((double) (assigned + unassigned));
    std::cerr << "Assigned: " << assigned << ", Unassigned: " << unassigned << ", PercentageAssigned: " << percass << "%" << std::endl;
    
#ifdef PROFILE
    ProfilerStop();
#endif    


    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] ";
    std::cout << "Done." << std::endl;
    
    return 0;
  }
  

  int bcsplit(int argc, char **argv) {
    BcsplitConfig c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("hamming,a", boost::program_options::value<uint32_t>(&c.hamming)->default_value(0), "max. hamming distance to barcode")
      ("ncount,n", boost::program_options::value<uint32_t>(&c.ncount)->default_value(0), "max. number of Ns per barcode")
      ("barcodes,b", boost::program_options::value<boost::filesystem::path>(&c.barcodes), "barcode file [barcode, id]")
      ("idxfile,i", boost::program_options::value<boost::filesystem::path>(&c.idxfile), "input index FASTQ file")
      ("pattern,p", boost::program_options::value<std::string>(&c.pattern)->default_value("BBBBBBBBUUUUUU"), "index read pattern [U: UMI, B: Barcode, N: Ignore]")
      ("outprefix,o", boost::program_options::value<std::string>(&c.outprefix)->default_value("split"), "output prefix")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.fastqfile), "input FASTQ file")
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
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("barcodes")) || (!vm.count("idxfile"))) {
      std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] -i index.fq.gz -b bar.tsv reads.fq.gz" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    }

    // Show cmd
    std::cout << '[' << boost::posix_time::to_simple_string(boost::posix_time::second_clock::local_time()) << "] ";
    std::cout << "alfred ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return runBcsplit(c);
  }

}

#endif

