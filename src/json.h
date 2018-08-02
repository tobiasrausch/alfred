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

#ifndef JSON_H
#define JSON_H

#include <boost/progress.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include "tsv.h"

namespace bamstats
{
  

  template<typename TConfig, typename TRGMap>
  inline void
  qcJsonOut(TConfig const& c, bam_hdr_t const* hdr, TRGMap const& rgMap, BedCounts const& be, ReferenceFeatures const& rf) {
    boost::iostreams::filtering_ostream rfile;
    rfile.push(boost::iostreams::gzip_compressor());
    rfile.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    
    rfile << "{ \"data\": [" << std::endl;
    rfile << "{" << std::endl;
    rfile << "\"sample\": \"" << c.sampleName << "\"," << std::endl;
    rfile << "\"rg\": [" << std::endl;

    // All read-groups
    for(typename TRGMap::const_iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
      if (itRg != rgMap.begin()) rfile << ", ";
      rfile << "{" << std::endl;
      rfile << "\"readGroup\": \"" << itRg->first << "\"," << std::endl;
      rfile << "\"metrics\": [" << std::endl;

      // Base content
      {
	rfile << "{" << std::endl;
	rfile << "\"name\": \"Base content distribution\"," << std::endl;
	rfile << "\"pos\": [";
	uint32_t lastValidBQIdx = _lastNonZeroIdxACGTN(itRg->second.rc);
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  if (i > 0) rfile << ", ";
	  rfile << i;
	}
	rfile << "]," << std::endl;
	rfile << "\"A\": [" << std::endl;
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  if (i > 0) rfile << ", ";
	  rfile << itRg->second.rc.aCount[i];
	}
	rfile << "]," << std::endl;
	rfile << "\"C\": [" << std::endl;
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  if (i > 0) rfile << ", ";
	  rfile << itRg->second.rc.cCount[i];
	}
	rfile << "]," << std::endl;
	rfile << "\"G\": [" << std::endl;
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  if (i > 0) rfile << ", ";
	  rfile << itRg->second.rc.gCount[i];
	}
	rfile << "]," << std::endl;
	rfile << "\"T\": [" << std::endl;
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  if (i > 0) rfile << ", ";
	  rfile << itRg->second.rc.tCount[i];
	}
	rfile << "]," << std::endl;
	rfile << "\"N\": [" << std::endl;
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  if (i > 0) rfile << ", ";
	  rfile << itRg->second.rc.nCount[i];
	}
	rfile << "]" << std::endl;
	rfile << "}," << std::endl;
      }
	
      // Read-length
      {
	rfile << "{" << std::endl;
	rfile << "\"name\": \"Read length distribution\"," << std::endl;
	rfile << "\"length\": [";
	uint32_t lastValidRL = _lastNonZeroIdx(itRg->second.rc.lRc);
	for(uint32_t i = 0; i <= lastValidRL; ++i) {
	  if (i > 0) rfile << ", ";
	  rfile << i;
	}
	rfile << "]," << std::endl;
	rfile << "\"count\": [" << std::endl;
	for(uint32_t i = 0; i <= lastValidRL; ++i) {
	  if (i > 0) rfile << ", ";
	  rfile << itRg->second.rc.lRc[i];
	}
	rfile << "]" << std::endl;
	rfile << "}," << std::endl;
      }
      
      // Mean Base Quality
      {
	rfile << "{" << std::endl;
	rfile << "\"name\": \"Mean base quality distribution\"," << std::endl;
	rfile << "\"pos\": [";
	uint32_t lastValidBQIdx = _lastNonZeroIdxACGTN(itRg->second.rc);
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  if (i > 0) rfile << ", ";
	  rfile << i;
	}
	rfile << "]," << std::endl;
	rfile << "\"qual\": [" << std::endl;
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  if (i > 0) rfile << ", ";
	  uint64_t bcount = itRg->second.rc.aCount[i] + itRg->second.rc.cCount[i] + itRg->second.rc.gCount[i] + itRg->second.rc.tCount[i] + itRg->second.rc.nCount[i];
	  if (bcount > 0) rfile << (double) (itRg->second.rc.bqCount[i]) / (double) (bcount);
	  else rfile << 0;
	}
	rfile << "]" << std::endl;
	rfile << "}" << std::endl;
      }
      

      
      rfile << "]" << std::endl;
      rfile << "}";
    }
    rfile << "]" << std::endl;
    rfile << "}" << std::endl;
    rfile << "]}" << std::endl;
    rfile.pop();
  }
 
}

#endif
