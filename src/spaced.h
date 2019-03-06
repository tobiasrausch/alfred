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

#ifndef SPACED_H
#define SPACED_H

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
#include "gtf.h"
#include "gff3.h"
#include "bed.h"
#include "motif.h"

namespace bamstats
{

  struct SpacedConfig {
    int32_t low;
    int32_t high;
    std::string motif1;
    std::string motif2;
    boost::filesystem::path infile;
    boost::filesystem::path outfile;
  };

  template<typename TConfig>
  inline int32_t
  spacedMotifRun(TConfig const& c) {
    int32_t motiflen1 = -1;
    int32_t motiflen2 =	-1;

    // Motif hits
    typedef std::vector<uint32_t> TMotifHit;
    typedef std::vector<TMotifHit> TGenomicMotifHit;
    typedef std::vector<TGenomicMotifHit> TStrandGenomicHit;
    TStrandGenomicHit mo1(2, TGenomicMotifHit());
    TStrandGenomicHit mo2(2, TGenomicMotifHit());

    // Chromosome map
    typedef std::map<std::string, uint32_t> TChrMap;
    TChrMap chrMap;
    uint32_t numChr = 0;
    std::vector<std::string> revMap;
    
    // Parse motif hits
    std::ifstream file(c.infile.string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string gline;
    while(std::getline(instream, gline)) {
      typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
      boost::char_separator<char> sep(" \t");
      Tokenizer tokens(gline, sep);
      Tokenizer::iterator tokIter = tokens.begin();
      if (tokIter != tokens.end()) {
	std::string chrName = *tokIter++;
	std::string pos = *tokIter;
	if (pos == "start") continue; // Header
	int32_t start = boost::lexical_cast<int32_t>(*tokIter++);
	int32_t end = boost::lexical_cast<int32_t>(*tokIter++);
	std::string id = *tokIter++;
	if ((id == c.motif1) || (id == c.motif2)) {
	  std::string strandStr = *tokIter;
	  char strand = strandStr[0];
	  uint32_t refIndex = numChr;
	  TChrMap::const_iterator it = chrMap.find(chrName);
	  if (it == chrMap.end()) {
	    chrMap.insert(std::make_pair(chrName, numChr));
	    revMap.push_back(chrName);
	    ++numChr;
	    mo1[0].resize(numChr, TMotifHit());
	    mo1[1].resize(numChr, TMotifHit());
	    mo2[0].resize(numChr, TMotifHit());
	    mo2[1].resize(numChr, TMotifHit());
	  } else refIndex = it->second;
	  if (id == c.motif1) {
	    int32_t len1 = (end - start) + 1;
	    if (motiflen1 == -1) motiflen1 = len1;
	    else if (motiflen1 != len1) {
	      std::cerr << "Warning: Motif hits have different lengths!" << std::endl;
	    }
	    if (strand == '+') mo1[0][refIndex].push_back(start);
	    else mo1[1][refIndex].push_back(start);
	  } else {
	    int32_t len2 = (end - start) + 1;
	    if (motiflen2 == -1) motiflen2 = len2;
	    else if (motiflen2 != len2) {
	      std::cerr	<< "Warning: Motif hits have different lengths!" << std::endl;
	    }
	    if (strand == '+') mo2[0][refIndex].push_back(start);
	    else mo2[1][refIndex].push_back(start);
	  }
	}
      }
    }
    dataIn.pop();

    // Motifs
    for(uint32_t strand = 0; strand < 2; ++strand) {
      char strandlabel = '+';
      if (strand) strandlabel = '-';
      for(uint32_t refIndex = 0; refIndex < revMap.size(); ++refIndex) {
	std::sort(mo1[strand][refIndex].begin(), mo1[strand][refIndex].end());
	std::sort(mo2[strand][refIndex].begin(), mo2[strand][refIndex].end());
	for(uint32_t i = 0; i < mo1[strand][refIndex].size(); ++i) {
	  for(uint32_t j = 0; j < mo2[strand][refIndex].size(); ++j) {
	    if (strand) {
	      if (mo2[strand][refIndex][j] + c.high < mo1[strand][refIndex][i]) continue;
	      if (mo2[strand][refIndex][j] + c.low > mo1[strand][refIndex][i]) break;
	      std::cout << revMap[refIndex] << ',' << mo1[strand][refIndex][i] << ',' << mo1[strand][refIndex][i] + motiflen1 - 1 << ',' << c.motif1 << ',' << strandlabel << ',' << revMap[refIndex] << ',' << mo2[strand][refIndex][j] << ',' << mo2[strand][refIndex][j] + motiflen2 - 1 << ',' << c.motif2 << ',' << strandlabel << std::endl;
	    } else {
	      if (mo1[strand][refIndex][i] + c.high < mo2[strand][refIndex][j]) break;
	      if (mo1[strand][refIndex][i] + c.low > mo2[strand][refIndex][j]) continue;
	      std::cout << revMap[refIndex] << ',' << mo1[strand][refIndex][i] << ',' << mo1[strand][refIndex][i] + motiflen1 - 1 << ',' << c.motif1 << ',' << strandlabel << ',' << revMap[refIndex] << ',' << mo2[strand][refIndex][j] << ',' << mo2[strand][refIndex][j] + motiflen2 - 1 << ',' << c.motif2 << ',' << strandlabel << std::endl;
	    }
	  }
	}
      }
    }
    
    
    
    return 0;
  }


  int spaced(int argc, char **argv) {
    SpacedConfig c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("motif1,m", boost::program_options::value<std::string>(&c.motif1)->default_value("Heptamer"), "motif1 name")
      ("motif2,n", boost::program_options::value<std::string>(&c.motif2)->default_value("Nonamer"), "motif2 name")
      ("spacer-low,l", boost::program_options::value<int32_t>(&c.low)->default_value(18), "min. spacer length (left-most to left-most)")
      ("spacer-high,h", boost::program_options::value<int32_t>(&c.high)->default_value(20), "max. spacer length (left-most to left-most)")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("joined.bed"), "joined motif hits")
      ;

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.infile), "input file")
      ;

    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);

    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);

    // Parse command-line
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file"))) {
      std::cout << std::endl;
      std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] <motif.hits.gz>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    }

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "alfred ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;

    return spacedMotifRun(c);
  }
  


  
}

#endif
