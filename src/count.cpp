/*
============================================================================
Alfred: BAM alignment statistics
============================================================================
Copyright (C) 2017 Tobias Rausch

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

#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <fstream>

#define BOOST_DISABLE_ASSERTS
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

#include "count.h"
#include "util.h"
#include "gtf.h"
#include "version.h"

using namespace bamstats;

struct Config {
  std::string sampleName;
  std::string outprefix;
  boost::filesystem::path gtfFile;
  boost::filesystem::path bamFile;
};


int main(int argc, char **argv) {

#ifdef PROFILE
  ProfilerStart("pbBamStats.prof");
#endif

  Config c;

  // Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("gtf,g", boost::program_options::value<boost::filesystem::path>(&c.gtfFile), "gtf file (required)")
    ("outprefix,o", boost::program_options::value<std::string>(&c.outprefix)->default_value("count"), "output file prefix")
    ;

  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bamFile), "input bam file")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic);

  // Only help/license/warranty/version information
  if ((argc < 2) || (std::string(argv[1]) == "help") || (std::string(argv[1]) == "--help") || (std::string(argv[1]) == "-h") || (std::string(argv[1]) == "-?")) {
    printTitle("Alfred");
    std::cout << "Usage: " << argv[0] << " [OPTIONS] -r <ref.fa> <aligned.bam>" << std::endl;
    std::cout << visible_options << "\n";
    return 0;
  } else if ((std::string(argv[1]) == "version") || (std::string(argv[1]) == "--version") || (std::string(argv[1]) == "--version-only") || (std::string(argv[1]) == "-v")) {
    std::cout << "Alfred version: v" << alfredVersionNumber << std::endl;
    return 0;
  } else if ((std::string(argv[1]) == "warranty") || (std::string(argv[1]) == "--warranty") || (std::string(argv[1]) == "-w")) {
    displayWarranty();
    return 0;
  } else if ((std::string(argv[1]) == "license") || (std::string(argv[1]) == "--license") || (std::string(argv[1]) == "-l")) {
    gplV3();
    return 0;
  }

  // Parse command-line
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);

  // Check command line arguments
  if ((!vm.count("input-file")) || (!vm.count("gtf"))) {
    printTitle("Alfred");
    std::cout << "Usage: " << argv[0] << " [OPTIONS] -g <hg19.gtf> <aligned.bam>" << std::endl;
    std::cout << visible_options << "\n";
    return 1;
  }

  // Check bam file
  if (!(boost::filesystem::exists(c.bamFile) && boost::filesystem::is_regular_file(c.bamFile) && boost::filesystem::file_size(c.bamFile))) {
    std::cerr << "Alignment file is missing: " << c.bamFile.string() << std::endl;
    return 1;
  } else {
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    if (samfile == NULL) {
      std::cerr << "Fail to open file " << c.bamFile.string() << std::endl;
      return 1;
    }
    hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
    if (idx == NULL) {
      if (bam_index_build(c.bamFile.string().c_str(), 0) != 0) {
	std::cerr << "Fail to open index for " << c.bamFile.string() << std::endl;
	return 1;
      }
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    
    // Get sample name
    std::string sampleName;
    if (!getSMTag(std::string(hdr->text), c.bamFile.stem().string(), sampleName)) {
      std::cerr << "Only one sample (@RG:SM) is allowed per input BAM file " << c.bamFile.string() << std::endl;
      return 1;
    } else c.sampleName = sampleName;
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
  }
  
  // Check region file
  if (!(boost::filesystem::exists(c.gtfFile) && boost::filesystem::is_regular_file(c.gtfFile) && boost::filesystem::file_size(c.gtfFile))) {
    std::cerr << "Input gtf file is missing: " << c.gtfFile.string() << std::endl;
    return 1;
  }

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  return countRun(c);
}
