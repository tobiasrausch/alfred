/*
============================================================================
Alfred: BAM alignment statistics
============================================================================
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
============================================================================
Contact: Tobias Rausch (rausch@embl.de)
============================================================================
*/

#ifndef REPLISEQ_H
#define REPLISEQ_H

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

#include "repliseqproc.h"

namespace bamstats {

  struct RepliSeqConfig {
    uint16_t minq;
    int32_t wsize;
    int32_t step;
    std::string outprefix;
    boost::filesystem::path genome;
    std::vector<boost::filesystem::path> files;
  };


  int repliseq(int argc, char **argv) {
    
#ifdef PROFILE
    ProfilerStart("repliseq.prof");
#endif

    RepliSeqConfig c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("qual,q", boost::program_options::value<uint16_t>(&c.minq)->default_value(1), "min. mapping quality")
      ("window,w", boost::program_options::value<int32_t>(&c.wsize)->default_value(50000), "sliding window size")
      ("step,s", boost::program_options::value<int32_t>(&c.step)->default_value(1000), "window offset (step size)")
      ("reference,r", boost::program_options::value<boost::filesystem::path>(&c.genome), "reference fasta file (required)")
      ("outprefix,o", boost::program_options::value<std::string>(&c.outprefix)->default_value("pref"), "output file prefix")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<std::vector<boost::filesystem::path> >(&c.files), "input bam files")
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
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("reference"))) {
      std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] -r <ref.fa> -o outprefix <g1b.bam> <s1.bam> <s2.bam> <s3.bam> <s4.bam> <g2.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    } 
    
    // Check genome
    if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
      std::cerr << "Input reference file is missing: " << c.genome.string() << std::endl;
      return 1;
    } else {
      faidx_t* fai = fai_load(c.genome.string().c_str());
      if (fai == NULL) {
	if (fai_build(c.genome.string().c_str()) == -1) {
	  std::cerr << "Fail to open genome fai index for " << c.genome.string() << std::endl;
	  return 1;
	} else fai = fai_load(c.genome.string().c_str());
      }
      fai_destroy(fai);
    }
    
    // Check bam files
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      if (!(boost::filesystem::exists(c.files[file_c]) && boost::filesystem::is_regular_file(c.files[file_c]) && boost::filesystem::file_size(c.files[file_c]))) {
	std::cerr << "Alignment file is missing: " << c.files[file_c].string() << std::endl;
	return 1;
      }
      samFile* samfile = sam_open(c.files[file_c].string().c_str(), "r");
      if (samfile == NULL) {
	std::cerr << "Fail to open file " << c.files[file_c].string() << std::endl;
	return 1;
      }
      hts_idx_t* idx = sam_index_load(samfile, c.files[file_c].string().c_str());
      if (idx == NULL) {
	std::cerr << "Fail to open index for " << c.files[file_c].string() << std::endl;
	return 1;
      }
      bam_hdr_t* hdr = sam_hdr_read(samfile);
      if (hdr == NULL) {
	std::cerr << "Fail to open header for " << c.files[file_c].string() << std::endl;
	return 1;
      }
      faidx_t* fai = fai_load(c.genome.string().c_str());
      for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
	std::string tname(hdr->target_name[refIndex]);
	if (!faidx_has_seq(fai, tname.c_str())) {
	  std::cerr << "BAM file chromosome " << hdr->target_name[refIndex] << " is NOT present in your reference file " << c.genome.string() << std::endl;
	  return 1;
	}
      }
      fai_destroy(fai);
      bam_hdr_destroy(hdr);
      hts_idx_destroy(idx);
      sam_close(samfile);
    }
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "alfred ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return repliseqRun(c);
  }

}

#endif
