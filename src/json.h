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
    std::string filename = c.outfile.string();
    if (c.format == "both") filename += ".json.gz";
    
    boost::iostreams::filtering_ostream rfile;
    rfile.push(boost::iostreams::gzip_compressor());
    rfile.push(boost::iostreams::file_sink(filename.c_str(), std::ios_base::out | std::ios_base::binary));
    
    rfile << "{ \"data\": [" << std::endl;
    // Reference information
    rfile << "{" << std::endl;
    rfile << "\"chromosomes\": [";
    bool firstVal = true;
    for(uint32_t i = 0; i < rf.chrGC.size(); ++i) {
      if (rf.chrGC[i].ncount + rf.chrGC[i].gccount > 0) {
	double total = hdr->target_len[i] - rf.chrGC[i].ncount;
	if (total > 0) {
	  double frac = (double) rf.chrGC[i].gccount / total;
	  if (!firstVal) rfile << ",";
	  else firstVal = false;
	  rfile << "{" << std::endl;
	  rfile << "\"name\": \"" << hdr->target_name[i] << "\"," << std::endl;
	  rfile << "\"size\": " << hdr->target_len[i] << "," << std::endl;
	  rfile << "\"ncount\": " << rf.chrGC[i].ncount << "," << std::endl;
	  rfile << "\"gccount\": " << rf.chrGC[i].gccount << "," << std::endl;
	  rfile << "\"gcfrac\": " << frac << std::endl;
	  rfile << "}" << std::endl;
	}
      }
    }
    rfile << "]" << std::endl;
    rfile << "}," << std::endl;

    // Sample information
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
	rfile << "}," << std::endl;
      }

      // Mapping quality histogram
      {
	rfile << "{" << std::endl;
	rfile << "\"name\": \"Mapping quality distribution\"," << std::endl;
	rfile << "\"pos\": [";
	uint32_t lastValidMQ = _lastNonZeroIdx(itRg->second.qc.qcount);
	for(uint32_t i = 0; i <= lastValidMQ; ++i) {
	  if (i > 0) rfile << ", ";
	  rfile << i;
	}
	rfile << "]," << std::endl;
	rfile << "\"qual\": [" << std::endl;
	for(uint32_t i = 0; i <= lastValidMQ; ++i) {
	  if (i > 0) rfile << ", ";
	  rfile << itRg->second.qc.qcount[i];
	}
	rfile << "]" << std::endl;
	rfile << "}," << std::endl;
      }

      // Coverage Histogram
      {
	rfile << "{" << std::endl;
	rfile << "\"name\": \"Coverage histogram\"," << std::endl;
	rfile << "\"coverage\": [";
	uint32_t lastValidCO = _lastNonZeroIdx(itRg->second.bc.bpWithCoverage);
	for(uint32_t i = 0; i <= lastValidCO; ++i) {
	  if (i > 0) rfile << ", ";
	  rfile << i;
	}
	rfile << "]," << std::endl;
	rfile << "\"count\": [" << std::endl;
	for(uint32_t i = 0; i <= lastValidCO; ++i) {
	  if (i > 0) rfile << ", ";
	  rfile << itRg->second.bc.bpWithCoverage[i];
	}
	rfile << "]" << std::endl;
	rfile << "}," << std::endl;
      }

      // Insert Size Histogram
      {
	rfile << "{" << std::endl;
	rfile << "\"name\": \"Insert size histogram\"," << std::endl;
	rfile << "\"insertSize\": [";
	uint32_t lastValidIS = _lastNonZeroIdxISize(itRg->second.pc);
	for(uint32_t i = 0; i <= lastValidIS; ++i) {
	  if (i > 0) rfile << ", ";
	  rfile << i;
	}
	rfile << "]," << std::endl;
	rfile << "\"fPlus\": [" << std::endl;
	for(uint32_t i = 0; i <= lastValidIS; ++i) {
	  if (i > 0) rfile << ", ";
	  rfile << itRg->second.pc.fPlus[i];
	}
	rfile << "]," << std::endl;
	rfile << "\"fMinus\": [" << std::endl;
	for(uint32_t i = 0; i <= lastValidIS; ++i) {
	  if (i > 0) rfile << ", ";
	  rfile << itRg->second.pc.fMinus[i];
	}
	rfile << "]," << std::endl;
	rfile << "\"rPlus\": [" << std::endl;
	for(uint32_t i = 0; i <= lastValidIS; ++i) {
	  if (i > 0) rfile << ", ";
	  rfile << itRg->second.pc.rPlus[i];
	}
	rfile << "]," << std::endl;
	rfile << "\"rMinus\": [" << std::endl;
	for(uint32_t i = 0; i <= lastValidIS; ++i) {
	  if (i > 0) rfile << ", ";
	  rfile << itRg->second.pc.rMinus[i];
	}
	rfile << "]" << std::endl;
	rfile << "}," << std::endl;
      }

      // Mapping statistics by chromosome
      {
	rfile << "{" << std::endl;
	rfile << "\"name\": \"Chromosome mapping statistics\"," << std::endl;
	rfile << "\"chromosomes\": [";
	uint64_t totalMappedChr = 0;
	for(uint32_t i = 0; i < itRg->second.rc.mappedchr.size(); ++i) totalMappedChr += itRg->second.rc.mappedchr[i];
	for(uint32_t i = 0; i < itRg->second.rc.mappedchr.size(); ++i) {
	  if (i > 0) rfile << ", ";
	  double frac = 0;
	  if (totalMappedChr > 0) frac = (double) itRg->second.rc.mappedchr[i] / (double) totalMappedChr;
	  double expect = (double) (hdr->target_len[i] - rf.chrGC[i].ncount) / (double) (rf.referencebp - rf.ncount);
	  double obsexprat = frac / expect;
	  rfile << "{" << std::endl;
	  rfile << "\"name\": \"" << hdr->target_name[i] << "\"," << std::endl;
	  rfile << "\"size\": " << hdr->target_len[i] << "," << std::endl;
	  rfile << "\"mapped\": " << itRg->second.rc.mappedchr[i] << "," << std::endl;
	  rfile << "\"fracTotal\": " << frac << "," << std::endl;
	  rfile << "\"observedExpected\": " << obsexprat << std::endl;
	  rfile << "}" << std::endl;
	}
	rfile << "]}" << std::endl;
      }

      // Bed specific data
      if (c.hasRegionFile) {
	// On target rate
	{
	  rfile << ",{" << std::endl;
	  rfile << "\"name\": \"On-target rate\"," << std::endl;
	  rfile << "\"targetExtension\": [";
	  uint64_t alignedbases = itRg->second.bc.matchCount + itRg->second.bc.mismatchCount;
	  typename BedCounts::TOnTargetMap::const_iterator itOT = be.onTarget.find(itRg->first);
	  for(uint32_t i = 0; i < itOT->second.size(); ++i) {
	    if (i > 0) rfile << ", ";
	    rfile << i * be.stepsize;
	  }
	  rfile << "]," << std::endl;
	  rfile << "\"fractionOnTarget\": [" << std::endl;
	  for(uint32_t i = 0; i < itOT->second.size(); ++i) {
	    if (i > 0) rfile << ", ";
	    rfile << (double) itOT->second[i] / (double) alignedbases;
	  }
	  rfile << "]" << std::endl;
	  rfile << "}" << std::endl;
	}

	// Avg. target coverage
	{
	  rfile << ",{" << std::endl;
	  rfile << "\"name\": \"Targets above coverage threshold\"," << std::endl;
	  std::vector<double> fracAboveCov;
	  uint32_t maxBedCov = _lastCoverageLevel(be, rf, hdr->n_targets, fracAboveCov);
	  rfile << "\"coverageLevel\": [";
	  for(uint32_t i = 0; i < maxBedCov; ++i) {
	    if (i > 0) rfile << ", ";
	    rfile << i;
	  }
	  rfile << "]," << std::endl;
	  rfile << "\"fractionAboveCoverage\": [" << std::endl;
	  for(uint32_t i = 0; i < maxBedCov; ++i) {
	    if (i > 0) rfile << ", ";
	    rfile << fracAboveCov[i];
	  }
	  rfile << "]" << std::endl;
	  rfile << "}" << std::endl;
	}
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
