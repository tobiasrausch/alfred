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

#ifndef BED_H
#define BED_H

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/algorithm/string.hpp>

#include <htslib/sam.h>

#include "util.h"

namespace bamstats
{

  template<typename TConfig, typename TGenomicRegions, typename TGeneIds, typename TProteinCoding>
  inline int32_t
  parseBEDAll(TConfig const& c, TGenomicRegions& overlappingRegions, TGeneIds& geneIds, TProteinCoding& pCoding) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "BED feature parsing" << std::endl;
    
    // Check gzip
    if (!is_gz(c.bedFile)) {
      std::cerr << "BED file is not gzipped!" << std::endl;
      return 0;
    }

    // Map IDs to integer
    typedef std::map<std::string, int32_t> TIdMap;
    TIdMap idMap;

    // Keep track of unique exon IDs
    int32_t eid = 0;

    // Parse BED
    std::ifstream file(c.bedFile.string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string gline;
    while(std::getline(instream, gline)) {
      if ((gline.size()) && (gline[0] == '#')) continue;
      typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
      boost::char_separator<char> sep(" \t,;");
      Tokenizer tokens(gline, sep);
      Tokenizer::iterator tokIter = tokens.begin();
      if (tokIter==tokens.end()) {
	std::cerr << "Empty line in BED file!" << std::endl;
	return 0;
      }
      std::string chrName=*tokIter++;
      if (c.nchr.find(chrName) == c.nchr.end()) continue;
      int32_t chrid = c.nchr.find(chrName)->second;
      if (tokIter == tokens.end()) {
	std::cerr << "Corrupted BED file!" << std::endl;
	return 0;
      }
      int32_t start = boost::lexical_cast<int32_t>(*tokIter++);
      int32_t end = boost::lexical_cast<int32_t>(*tokIter++);
      if (tokIter == tokens.end()) {
	std::cerr << "Name is missing in BED file!" << std::endl;
	return 0;
      }
      std::string val = *tokIter++;
      char strand = '*';
      std::string biotype = "NA";
      if (tokIter != tokens.end()) {
	++tokIter; // skip score
	strand = boost::lexical_cast<char>(*tokIter++);
	biotype = *tokIter++;
      }
      int32_t idval = geneIds.size();
      typename TIdMap::const_iterator idIter = idMap.find(val);
      if (idIter == idMap.end()) {
	idMap.insert(std::make_pair(val, idval));
	geneIds.push_back(val);
	if (biotype == "protein_coding") pCoding.push_back(true);
	else pCoding.push_back(false);
      } else idval = idIter->second;
      // BED is 0-based and right-open, no need to convert
      if (start > end) {
	std::cerr << "Feature start is greater than feature end!" << std::endl;
	return 0;
      }
      //std::cerr << geneIds[idval] << "\t" << start << "\t" << end << std::endl;
      _insertInterval(overlappingRegions[chrid], start, end, strand, idval, eid++);
    }
    return geneIds.size();
  }   

  template<typename TConfig, typename TGenomicRegions, typename TGeneIds>
  inline int32_t
  parseBEDAll(TConfig const& c, TGenomicRegions& overlappingRegions, TGeneIds& geneIds) {
    std::vector<bool> pCoding;
    return parseBEDAll(c, overlappingRegions, geneIds, pCoding);
  }
  

  template<typename TConfig, typename TGenomicRegions, typename TGeneIds, typename TProteinCoding>
  inline int32_t
  parseBED(TConfig const& c, TGenomicRegions& gRegions, TGeneIds& geneIds, TProteinCoding& pCoding) {
    typedef typename TGenomicRegions::value_type TChromosomeRegions;

    // Overlapping intervals for each label
    TGenomicRegions overlappingRegions;
    overlappingRegions.resize(gRegions.size(), TChromosomeRegions());
    parseBEDAll(c, overlappingRegions, geneIds, pCoding);

    // Make intervals non-overlapping for each label
    for(uint32_t refIndex = 0; refIndex < overlappingRegions.size(); ++refIndex) {
      // Sort by ID
      std::sort(overlappingRegions[refIndex].begin(), overlappingRegions[refIndex].end(), SortIntervalLabel<IntervalLabel>());
      int32_t runningId = -1;
      char runningStrand = '*';
      typedef boost::icl::interval_set<uint32_t> TIdIntervals;
      typedef typename TIdIntervals::interval_type TIVal;
      TIdIntervals idIntervals;
      for(uint32_t i = 0; i < overlappingRegions[refIndex].size(); ++i) {
	if (overlappingRegions[refIndex][i].lid != runningId) {
	  for(typename TIdIntervals::iterator it = idIntervals.begin(); it != idIntervals.end(); ++it) {
	    //std::cerr << "merged\t" << geneIds[runningId] << "\t" << it->lower() << "\t" << it->upper() << std::endl;  
	    gRegions[refIndex].push_back(IntervalLabel(it->lower(), it->upper(), runningStrand, runningId));
	  }
	  idIntervals.clear();
	  runningId = overlappingRegions[refIndex][i].lid;
	  runningStrand = overlappingRegions[refIndex][i].strand;
	}
	idIntervals.insert(TIVal::right_open(overlappingRegions[refIndex][i].start, overlappingRegions[refIndex][i].end));
      }
      // Process last id
      for(typename TIdIntervals::iterator it = idIntervals.begin(); it != idIntervals.end(); ++it) gRegions[refIndex].push_back(IntervalLabel(it->lower(), it->upper(), runningStrand, runningId));
    }
    return geneIds.size();
  }

  template<typename TConfig, typename TGenomicRegions, typename TGeneIds>
  inline int32_t
  parseBED(TConfig const& c, TGenomicRegions& gRegions, TGeneIds& geneIds) {
    std::vector<bool> pCoding;
    return parseBED(c, gRegions, geneIds, pCoding);
  }
  

}

#endif
