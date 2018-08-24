/*
============================================================================
Alfred: BAM alignment statistics
============================================================================
Copyright (C) 2017-2018 Tobias Rausch

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

#ifndef SCAN_H
#define SCAN_H

#include <limits>

#include <boost/icl/split_interval_map.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/progress.hpp>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "version.h"
#include "util.h"
#include "motif.h"


namespace bamstats
{

  struct ScanConfig {
    std::map<std::string, int32_t> nchr;
    boost::filesystem::path infile;
    boost::filesystem::path motifFile;
    boost::filesystem::path outfile;
  };


  template<typename TConfig>
  inline int32_t
  scanRun(TConfig const& c) {

#ifdef PROFILE
    ProfilerStart("alfred.prof");
#endif

    std::vector<Pfm> pfms;
    if (!parseJaspar(c, pfms)) {
      std::cerr << "Motif file cannot be parsed!" << std::endl;
      return 1;
    }
    std::cout << pfms.size() << std::endl;
      
    
    // Done
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    
#ifdef PROFILE
    ProfilerStop();
#endif

    return 0;
  }


  int scan(int argc, char **argv) {
    ScanConfig c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("anno.bed"), "output file")
      ;

    boost::program_options::options_description motifopt("Motif options");
    motifopt.add_options()
      ("motif,m", boost::program_options::value<boost::filesystem::path>(&c.motifFile), "motif file in jaspar or raw format")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.infile), "input file")
      ;

    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);

    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(motifopt).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(motifopt);

    // Parse command-line
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("motif"))) {
      std::cout << std::endl;
      std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] -m <motif.jaspar> <peaks.bed>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    }

    // Input BED file
    if (!(boost::filesystem::exists(c.infile) && boost::filesystem::is_regular_file(c.infile) && boost::filesystem::file_size(c.infile))) {
      std::cerr << "Input BED file is missing." << std::endl;
      return 1;
    } else {
      std::string oldChr = "";
      typedef std::set<std::string> TChrSet;
      TChrSet chrSet;
      std::ifstream chrFile(c.infile.string().c_str(), std::ifstream::in);
      if (chrFile.is_open()) {
	while (chrFile.good()) {
	  std::string chrFromFile;
	  getline(chrFile, chrFromFile);
	  typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	  boost::char_separator<char> sep(" \t,;");
	  Tokenizer tokens(chrFromFile, sep);
	  Tokenizer::iterator tokIter = tokens.begin();
	  if (tokIter!=tokens.end()) {
	    std::string chrName = *tokIter++;
	    if (chrName != oldChr) chrSet.insert(chrName);
	  }
	}
	chrFile.close();
      }
      int32_t refIndex = 0;
      for(TChrSet::iterator itc = chrSet.begin(); itc != chrSet.end(); ++itc, ++refIndex) c.nchr.insert(std::make_pair(*itc, refIndex));
    }
    
    // Check region file
    if (!(boost::filesystem::exists(c.motifFile) && boost::filesystem::is_regular_file(c.motifFile) && boost::filesystem::file_size(c.motifFile))) {
      std::cerr << "Input motif file is missing." << std::endl;
      return 1;
    }
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "alfred motif ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;

    return scanRun(c);
  }
  


  
}

#endif
