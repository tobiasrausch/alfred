#ifndef PWALIGN_H
#define PWALIGN_H

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <vector>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <math.h>
#include <stdio.h>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "util.h"
#include "align.h"
#include "needle.h"
#include "gotoh.h"
#include "swneedle.h"
#include "swgotoh.h"

namespace bamstats {


struct PWAlignConsensus {
  bool seq1endsfree;
  bool seq2endsfree;
  bool localAlignment;
  int32_t gapopen;
  int32_t gapext;
  int32_t match;
  int32_t mismatch;
  std::string format;
  boost::filesystem::path alignment;
  std::vector<boost::filesystem::path> inputfiles;
};


int pwalign(int argc, char **argv) {
  PWAlignConsensus c;

  // Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("gapopen,g", boost::program_options::value<int32_t>(&c.gapopen)->default_value(-10), "gap open")
    ("gapext,e", boost::program_options::value<int32_t>(&c.gapext)->default_value(-1), "gap extension")
    ("match,m", boost::program_options::value<int32_t>(&c.match)->default_value(5), "match")
    ("mismatch,n", boost::program_options::value<int32_t>(&c.mismatch)->default_value(-4), "mismatch")
    ("endsfree1,p", "leading/trailing gaps free for seq1")
    ("endsfree2,q", "leading/trailing gaps free for seq2")
    ("local,l", "local alignment")
    ;

  boost::program_options::options_description otp("Output options");
  otp.add_options()
    ("format,f", boost::program_options::value<std::string>(&c.format)->default_value("h"), "output format [v|h]")
    ("alignment,a", boost::program_options::value<boost::filesystem::path>(&c.alignment)->default_value("al.fa.gz"), "vertical/horizontal alignment")
    ;
  
  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value< std::vector<boost::filesystem::path> >(&c.inputfiles), "input fasta file")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(otp).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic).add(otp);
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);

  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file")) || (c.inputfiles.size() != 2)) {
    std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] <seq1.fasta> <seq2.fasta>" << std::endl;
    std::cout << visible_options << "\n";
    return 1;
  }

  // Flags
  if (vm.count("endsfree1")) c.seq1endsfree = true;
  else c.seq1endsfree = false;
  if (vm.count("endsfree2")) c.seq2endsfree = true;
  else c.seq2endsfree = false;
  if (vm.count("local")) c.localAlignment = true;
  else c.localAlignment = false;
  
  // Check input files
  for(unsigned int file_c = 0; file_c < c.inputfiles.size(); ++file_c) {
    if (!(boost::filesystem::exists(c.inputfiles[file_c]) && boost::filesystem::is_regular_file(c.inputfiles[file_c]) && boost::filesystem::file_size(c.inputfiles[file_c]))) {
      std::cerr << "Input fasta file is missing: " << c.inputfiles[file_c].string() << std::endl;
      return 1;
    }
  }
  
  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  std::cout << "alfred ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Load FASTA sequences
  std::string faname1;
  std::string seq1;
  loadSingleFasta(c.inputfiles[0].string(), faname1, seq1);
  std::cout << "Sequence1: " << faname1 << ", Length: " << seq1.size() << std::endl;
  std::string faname2;
  std::string seq2;
  loadSingleFasta(c.inputfiles[1].string(), faname2, seq2);
  std::cout << "Sequence2: " << faname2 << ", Length: " << seq2.size() << std::endl;

  // Alignment
  typedef boost::multi_array<char, 2> TAlign;
  TAlign align;
  DnaScore<int> sc(c.match, c.mismatch, c.gapopen, c.gapext);
  int32_t alScore = 0;
  if (c.localAlignment) {
    AlignConfig<false, false> alignconf;
    alScore = swGotoh(seq1, seq2, align, alignconf, sc);
  } else {
    if (c.seq1endsfree) {
      if (c.seq2endsfree) {
	AlignConfig<true, true> alignconf;
	alScore = gotoh(seq1, seq2, align, alignconf, sc);
      } else {
	AlignConfig<true, false> alignconf;
	alScore = gotoh(seq1, seq2, align, alignconf, sc);
      }
    } else {
      if (c.seq2endsfree) {
	AlignConfig<false, true> alignconf;
	alScore = gotoh(seq1, seq2, align, alignconf, sc);
      } else {
	AlignConfig<false, false> alignconf;
	alScore = gotoh(seq1, seq2, align, alignconf, sc);
      }
    }
  }
  std::cout << "Alignment score: " << alScore << std::endl;

  // Output
  if (c.format == "h") {
    boost::iostreams::filtering_ostream rcfile;
    rcfile.push(boost::iostreams::gzip_compressor());
    rcfile.push(boost::iostreams::file_sink(c.alignment.c_str(), std::ios_base::out | std::ios_base::binary));
    typedef typename TAlign::index TAIndex;
    for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {
      if (i == 0) rcfile << ">" << faname1 << std::endl;
      else rcfile << ">" << faname2 << std::endl;
      for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
	rcfile << align[i][j];
      }
      rcfile << std::endl;
    }
    rcfile.pop();
  } else {
    boost::iostreams::filtering_ostream rcfile;
    rcfile.push(boost::iostreams::gzip_compressor());
    rcfile.push(boost::iostreams::file_sink(c.alignment.c_str(), std::ios_base::out | std::ios_base::binary));
    typedef typename TAlign::index TAIndex;
    for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
      for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {
	rcfile << align[i][j];
      }
      rcfile << std::endl;
    }
    rcfile.pop();
  }
  
  // Done
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
  
  return 0;
}


}

#endif
