#ifndef PWEDIT_H
#define PWEDIT_H

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

#include "edlib.h"
#include "util.h"

namespace bamstats {


  struct PWEditConfig {
    bool reverseComplement;
    std::string mode;
    std::string format;
    boost::filesystem::path alignment;
    std::vector<boost::filesystem::path> inputfiles;
  };

  inline void
  printAlignmentPretty(std::string const& query, std::string const& target, EdlibAlignMode const modeCode, EdlibAlignResult& align) {
    int32_t tIdx = -1;
    int32_t qIdx = -1;
    if (modeCode == EDLIB_MODE_HW) {
        tIdx = align.endLocations[0];
        for (int32_t i = 0; i < align.alignmentLength; i++) {
            if (align.alignment[i] != EDLIB_EDOP_INSERT) tIdx--;
        }
    }
    std::cout << std::endl;
    for (int start = 0; start < align.alignmentLength; start += 50) {
      std::cout << "T: ";
      int32_t startTIdx = -1;
      for (int32_t j = start; ((j < start + 50) && (j < align.alignmentLength)); ++j) {
	if (align.alignment[j] == EDLIB_EDOP_INSERT) std::cout << "-";
	else std::cout << target[++tIdx];
	if (j == start) startTIdx = tIdx;
      }
      std::cout << " (" << std::max(startTIdx, 0) << " - " << tIdx << ")" << std::endl;

      // match / mismatch
      std::cout << ("   ");
      for (int32_t j = start; j < start + 50 && j < align.alignmentLength; j++) {
	if (align.alignment[j] == EDLIB_EDOP_MATCH) std::cout <<  "|";
	else std::cout << " ";
      }
      std::cout << std::endl;

      // query
      std::cout << "Q: ";
      int32_t startQIdx = qIdx;
      for (int32_t j = start; j < start + 50 && j < align.alignmentLength; j++) {
	if (align.alignment[j] == EDLIB_EDOP_DELETE) std::cout << "-";
	else std::cout << query[++qIdx];
	if (j == start) startQIdx = qIdx;
      }
      std::cout << " ("<< std::max(startQIdx, 0) << " - " << qIdx << ")" << std::endl;
      std::cout << std::endl;
    }
  }

  template<typename TConfig>
  inline void
  printAlignmentHorizontal(TConfig const& c, std::string const& seqI, std::string const& tname, std::string const& seqJ, std::string const& qname, EdlibAlignMode const modeCode, EdlibAlignResult& cigar) {
    boost::iostreams::filtering_ostream rcfile;
    rcfile.push(boost::iostreams::gzip_compressor());
    rcfile.push(boost::iostreams::file_sink(c.alignment.c_str(), std::ios_base::out | std::ios_base::binary));
    rcfile << ">" << tname << std::endl;
      
    int32_t tIdx = -1;
    int32_t qIdx = -1;
    uint32_t missingEnd = 0;
    uint32_t missingStart = 0;
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      tIdx = cigar.endLocations[0];
      if (tIdx < (int32_t) seqJ.size()) missingEnd = seqJ.size() - tIdx - 1;
      for (int32_t i = 0; i < cigar.alignmentLength; i++) {
	if (cigar.alignment[i] != EDLIB_EDOP_INSERT) tIdx--;
      }
      if (tIdx >= 0) missingStart = tIdx + 1;
    }
    // infix alignment, fix start
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      if (missingStart) {
	for (uint32_t j = 0; j < missingStart; ++j) rcfile << '-';
      }
    }
    // seqI
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_DELETE) rcfile << '-';
      else rcfile << seqI[++qIdx];
    }
    // infix alignment, fix end
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      if (missingEnd) {
	for (uint32_t j = 0; j < missingEnd; ++j) rcfile << '-';
      }
    }
    rcfile << std::endl;

    // Query
    rcfile << ">" << qname << std::endl;

    // infix alignment, fix start
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      if (missingStart) {
	for (uint32_t j = 0; j < missingStart; ++j) rcfile << seqJ[j];
      }
    }
    // seqJ
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_INSERT) rcfile << '-';
      else rcfile << seqJ[++tIdx];
    }
    // infix alignment, fix end
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      if (missingEnd) {
	for (uint32_t j = 0; j < missingEnd; ++j) rcfile << seqJ[++tIdx];
      }
    }
    rcfile << std::endl;
    rcfile.pop();
  }

  template<typename TConfig>
  inline void
  printAlignmentVertical(TConfig const& c, std::string const& seqI, std::string const& seqJ, EdlibAlignMode const modeCode, EdlibAlignResult& cigar) {
    std::string qa;
    std::string ta;
    
    int32_t tIdx = -1;
    int32_t qIdx = -1;
    uint32_t missingEnd = 0;
    uint32_t missingStart = 0;
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      tIdx = cigar.endLocations[0];
      if (tIdx < (int32_t) seqJ.size()) missingEnd = seqJ.size() - tIdx - 1;
      for (int32_t i = 0; i < cigar.alignmentLength; i++) {
	if (cigar.alignment[i] != EDLIB_EDOP_INSERT) tIdx--;
      }
      if (tIdx >= 0) missingStart = tIdx + 1;
    }
    // infix alignment, fix start
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      if (missingStart) {
	for (uint32_t j = 0; j < missingStart; ++j) qa += '-';
      }
    }
    // seqI
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_DELETE) qa += '-';
      else qa += seqI[++qIdx];
    }
    // infix alignment, fix end
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      if (missingEnd) {
	for (uint32_t j = 0; j < missingEnd; ++j) qa += '-';
      }
    }

    // infix alignment, fix start
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      if (missingStart) {
	for (uint32_t j = 0; j < missingStart; ++j) ta += seqJ[j];
      }
    }
    // seqJ
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_INSERT) ta += '-';
      else ta += seqJ[++tIdx];
    }
    // infix alignment, fix end
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      if (missingEnd) {
	for (uint32_t j = 0; j < missingEnd; ++j) ta +=  seqJ[++tIdx];
      }
    }

    // Output
    boost::iostreams::filtering_ostream rcfile;
    rcfile.push(boost::iostreams::gzip_compressor());
    rcfile.push(boost::iostreams::file_sink(c.alignment.c_str(), std::ios_base::out | std::ios_base::binary));
    for(uint32_t i = 0; (i < qa.size()) && (i < ta.size()); ++i) {
      rcfile << qa[i] << ta[i] << std::endl;
    }
    rcfile << std::endl;
    rcfile.pop();
  }
  

  int pwedit(int argc, char **argv) {
    PWEditConfig c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("mode,m", boost::program_options::value<std::string>(&c.mode)->default_value("infix"), "alignment mode [global|prefix|infix]")
      ("revcomp,r", "reverse complement query")
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
      std::cerr << "Usage: alfred " << argv[0] << " [OPTIONS] <target.fasta> <query.fasta>" << std::endl;
      std::cerr << visible_options << "\n";
      return 1;
    }
    
    // Flags
    if (vm.count("revcomp")) c.reverseComplement = true;
    else c.reverseComplement = false;
    
    // Check input files
    for(unsigned int file_c = 0; file_c < c.inputfiles.size(); ++file_c) {
      if (!(boost::filesystem::exists(c.inputfiles[file_c]) && boost::filesystem::is_regular_file(c.inputfiles[file_c]) && boost::filesystem::file_size(c.inputfiles[file_c]))) {
	std::cerr << "Input fasta file is missing: " << c.inputfiles[file_c].string() << std::endl;
	return 1;
      }
    }
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cerr << "alfred ";
    for(int i=0; i<argc; ++i) { std::cerr << argv[i] << ' '; }
    std::cerr << std::endl;
    
    // Load FASTA sequences
    std::string tname;
    std::string target;
    if (!loadSingleFasta(c.inputfiles[0].string(), tname, target, false)) return 1;
    std::cerr << "Target: " << tname << ", Length: " << target.size() << std::endl;
    std::string qname;
    std::string query;
    if (!loadSingleFasta(c.inputfiles[1].string(), qname, query, false)) return 1;
    std::cerr << "Query: " << qname << ", Length: " << query.size() << std::endl;

    // Reverse complement
    if (c.reverseComplement) reverseComplement(query);
    
    // Alignment
    EdlibAlignMode alnMode = EDLIB_MODE_HW;
    if (c.mode == "global") alnMode = EDLIB_MODE_NW;
    else if (c.mode == "prefix") alnMode = EDLIB_MODE_SHW;
    EdlibAlignResult aln = edlibAlign(query.c_str(), query.size(), target.c_str(), target.size(), edlibNewAlignConfig(-1, alnMode, EDLIB_TASK_PATH, NULL, 0));
    std::cout << "Alignment score: " << aln.editDistance << std::endl;
    std::cout << "Normalized score (by query length): " << ((double) aln.editDistance / (double) query.size()) << std::endl;
    printAlignmentPretty(query, target, alnMode, aln);
    
    // Output
    if (c.format == "h") {
      printAlignmentHorizontal(c, query, qname, target, tname, alnMode, aln);
    } else {
      printAlignmentVertical(c, query, target, alnMode, aln);
    }

    // Clean-up
    edlibFreeAlignResult(aln);

    // Done
    now = boost::posix_time::second_clock::local_time();
    std::cerr << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    
    return 0;
  }


}

#endif
