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

#include "bamstats.h"
#include "util.h"
#include "version.h"

using namespace bamstats;

struct Config {
  bool hasRegionFile;
  std::string sampleName;
  std::string outprefix;
  boost::filesystem::path genome;
  boost::filesystem::path regionFile;
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
    ("reference,r", boost::program_options::value<boost::filesystem::path>(&c.genome), "reference fasta file (required)")
    ("bed,b", boost::program_options::value<boost::filesystem::path>(&c.regionFile), "bed file with regions to analyze (optional)")
    ("outprefix,o", boost::program_options::value<std::string>(&c.outprefix)->default_value("stats"), "output file prefix")
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
  if ((!vm.count("input-file")) || (!vm.count("reference"))) {
    printTitle("Alfred");
    std::cout << "Usage: " << argv[0] << " [OPTIONS] -r <ref.fa> <aligned.bam>" << std::endl;
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

  // Check bam file
  if (!(boost::filesystem::exists(c.bamFile) && boost::filesystem::is_regular_file(c.bamFile) && boost::filesystem::file_size(c.bamFile))) {
    std::cerr << "Alignment file is missing: " << c.bamFile.string() << std::endl;
    return 1;
  }
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
  faidx_t* fai = fai_load(c.genome.string().c_str());
  for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
    std::string tname(hdr->target_name[refIndex]);
    if (!faidx_has_seq(fai, tname.c_str())) {
      std::cerr << "BAM file chromosome " << hdr->target_name[refIndex] << " is NOT present in your reference file " << c.genome.string() << std::endl;
      return 1;
    }
  }
  fai_destroy(fai);
  std::string sampleName;
  if (!getSMTag(std::string(hdr->text), c.bamFile.stem().string(), sampleName)) {
    std::cerr << "Only one sample (@RG:SM) is allowed per input BAM file " << c.bamFile.string() << std::endl;
    return 1;
  } else c.sampleName = sampleName;
  bam_hdr_destroy(hdr);
  hts_idx_destroy(idx);
  sam_close(samfile);
  
  // Check region file
  if (vm.count("bed")) {
    if (!(boost::filesystem::exists(c.regionFile) && boost::filesystem::is_regular_file(c.regionFile) && boost::filesystem::file_size(c.regionFile))) {
      std::cerr << "Input region file in bed format is missing: " << c.regionFile.string() << std::endl;
      return 1;
    }
    std::string oldChr;
    faidx_t* fai = fai_load(c.genome.string().c_str());
    std::ifstream interval_file(c.regionFile.string().c_str(), std::ifstream::in);
    if (interval_file.is_open()) {
      while (interval_file.good()) {
	std::string intervalLine;
	getline(interval_file, intervalLine);
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t,;");
	Tokenizer tokens(intervalLine, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter!=tokens.end()) {
	  std::string chrName=*tokIter++;
	  if (chrName.compare(oldChr) != 0) {
	    oldChr = chrName;
	    if (!faidx_has_seq(fai, chrName.c_str())) {
	      std::cerr << "Chromosome from bed file " << chrName << " is NOT present in your reference file " << c.genome.string() << std::endl;
	      return 1;
	    }
	  }
	}
      }
      interval_file.close();
    }
    fai_destroy(fai);
    c.hasRegionFile = true;
  } else c.hasRegionFile = false;

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  return bamStatsRun(c); 
}
