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
    typedef std::map<std::string, uint32_t> TChrMap;
    TChrMap nchr;
    boost::filesystem::path infile;
    boost::filesystem::path motifFile;
    boost::filesystem::path genome;
    boost::filesystem::path outfile;
  };


  template<typename TConfig>
  inline int32_t
  scanRun(TConfig const& c) {

#ifdef PROFILE
    ProfilerStart("alfred.prof");
#endif

    // Load motifs
    std::vector<Pwm> pwms;
    if (!parseJasparPwm(c, pwms)) {
      std::cerr << "Motif file cannot be parsed!" << std::endl;
      return 1;
    }

    // Motif finding
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Motif finding" << std::endl;
    boost::progress_display show_progress(c.nchr.size());
    
    faidx_t* fai = fai_load(c.genome.string().c_str());
    char* seq = NULL;
    for(uint32_t refIndex=0; refIndex < c.nchr.size(); ++refIndex) {
      ++show_progress;

      // Load chromosome sequence
      std::string tname = "NA";
      for(typename ScanConfig::TChrMap::const_iterator itChr = c.nchr.begin(); itChr != c.nchr.end(); ++itChr) {
	if (refIndex == itChr->second) {
	  tname = itChr->first;
	}
      }
      int32_t seqlen = -1;
      seq = faidx_fetch_seq(fai, tname.c_str(), 0, faidx_seq_len(fai, tname.c_str()) + 1, &seqlen);
      std::vector<uint8_t> ref(seqlen, 0);
      for(int32_t i = 0; i < seqlen; ++i) {
	if ((seq[i] == 'A') || (seq[i] == 'a')) ref[i] = 0;
	else if ((seq[i] == 'C') || (seq[i] == 'c')) ref[i] = 1;
	else if ((seq[i] == 'G') || (seq[i] == 'g')) ref[i] = 2;
	else if ((seq[i] == 'T') || (seq[i] == 't')) ref[i] = 3;
	else ref[i] = 4;
      }
      if (seq != NULL) free(seq);

      
      // Score PWMs
      typedef std::pair<int32_t, int32_t> TPosScore;
      typedef std::vector<TPosScore> TChrPosScore;
      typedef std::vector<TChrPosScore> TStrandChrPosScore;
      typedef std::vector<TStrandChrPosScore> TMotifHits;
      TMotifHits mh(pwms.size(), TStrandChrPosScore());
      typedef std::vector<int32_t> TMotifCounts;
      TMotifCounts mc(pwms.size(), 0);
      for(uint32_t i = 0; i<pwms.size(); ++i) mh[i].resize(2, TChrPosScore());      
      for(uint32_t i = 0; i<pwms.size(); ++i) {
	typedef std::pair<int32_t, int32_t> TFwdRevHits;
	TFwdRevHits hits = scorePwm(ref, pwms[i], 0.8, mh);
	//std::cout << i << "," << hits.first << "," << hits.second << std::endl;
      }
    }
    fai_destroy(fai);
    
    // Done
    now = boost::posix_time::second_clock::local_time();
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
      ("reference,r", boost::program_options::value<boost::filesystem::path>(&c.genome), "reference fasta file (required)")
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
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("motif")) || (!vm.count("reference"))) {
      std::cout << std::endl;
      std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] -r <genome.fa> -m <motif.jaspar> <peaks.bed>" << std::endl;
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
    
    // Input BED file
    if (!(boost::filesystem::exists(c.infile) && boost::filesystem::is_regular_file(c.infile) && boost::filesystem::file_size(c.infile))) {
      std::cerr << "Input BED file is missing." << std::endl;
      return 1;
    } else {
      std::string oldChr = "";
      faidx_t* fai = fai_load(c.genome.string().c_str());
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
	    if (chrName.compare(oldChr) != 0) {
	      oldChr = chrName;
	      if (!faidx_has_seq(fai, chrName.c_str())) {
		std::cerr << "Chromosome from bed file " << chrName << " is NOT present in your reference file " << c.genome.string() << std::endl;
		return 1;
	      } else {
		chrSet.insert(chrName);
	      }
	    }
	  }
	}
	chrFile.close();
      }
      uint32_t refIndex = 0;
      for(TChrSet::iterator itc = chrSet.begin(); itc != chrSet.end(); ++itc, ++refIndex) c.nchr.insert(std::make_pair(*itc, refIndex));
      fai_destroy(fai);
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
