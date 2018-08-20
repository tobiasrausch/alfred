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

    // Sample information
    rfile << "{\"samples\": [{";
    rfile << "\"id\": \"" << c.sampleName << "\",";

    // Summary Table
    rfile << "\"summary\": ";
    {
      rfile << "{\"id\": \"summaryTable\",";
      rfile << "\"title\": \"Summary Statistics\",";
      rfile << "\"data\": {\"columns\": [\"Sample\", \"Library\", \"#QCFail\", \"QCFailFraction\", \"#DuplicateMarked\", \"DuplicateFraction\", \"#Unmapped\", \"UnmappedFraction\", \"#Mapped\", \"MappedFraction\"], \"rows\": ["; 
      for(typename TRGMap::const_iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
	uint64_t totalReadCount = _totalReadCount(itRg);
	uint64_t mappedCount = itRg->second.rc.mapped1 + itRg->second.rc.mapped2;
	if (itRg != rgMap.begin()) rfile << ",";
	rfile << "[";
	rfile << "\"" << c.sampleName << "\"" << ",";
	rfile << "\"" << itRg->first << "\"" << ",";
	rfile << "\"" << itRg->second.rc.qcfail << "\"" << ",";
	rfile << "\"" << (double) itRg->second.rc.qcfail / (double) totalReadCount << "\"" << ",";
	rfile << "\"" << itRg->second.rc.dup << "\"" << ",";
	rfile << "\"" << (double) itRg->second.rc.dup / (double) totalReadCount << "\"" << ",";
	rfile << "\"" << itRg->second.rc.unmap << "\"" << ",";
	rfile << "\"" << (double) itRg->second.rc.unmap / (double) totalReadCount << "\"" << ",";
	rfile << "\"" << mappedCount << "\"" << ",";
	rfile << "\"" << (double) mappedCount / (double) totalReadCount << "\"";
	rfile << "]";
      }
      rfile << "]},";
      rfile << "\"type\": \"table\"}";
    }
    rfile << ",";
    
    // Read-group information
    rfile << "\"readGroups\": [";

    // All read-groups
    for(typename TRGMap::const_iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
      if (itRg != rgMap.begin()) rfile << ", ";
      rfile << "{";
      rfile << "\"id\": \"" << itRg->first << "\",";
      rfile << "\"metrics\": [";

      // Base content
      {
	rfile << "{\"id\": \"baseContent\",";
	rfile << "\"title\": \"Base content distribution\",";
	rfile << "\"x\": {\"data\": [{\"values\": [";
	uint32_t lastValidBQIdx = _lastNonZeroIdxACGTN(itRg->second.rc);
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << i;
	}
	rfile << "]}], \"axis\": {\"title\": \"Position in read\"}},";
	rfile << "\"y\": {\"data\": [";
	rfile << "{\"values\": [";
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << itRg->second.rc.aCount[i];
	}
	rfile << "], \"title\": \"A\"},";
	rfile << "{\"values\": [";
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << itRg->second.rc.cCount[i];
	}
	rfile << "], \"title\": \"C\"},";
	rfile << "{\"values\": [";
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << itRg->second.rc.gCount[i];
	}
	rfile << "], \"title\": \"G\"},";
	rfile << "{\"values\": [";
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << itRg->second.rc.tCount[i];
	}
	rfile << "], \"title\": \"T\"},";
	rfile << "{\"values\": [";
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << itRg->second.rc.nCount[i];
	}
	rfile << "], \"title\": \"N\"}], \"axis\": {\"title\": \"Count\"}}, \"type\": \"line\"}";
      }
	
      // Read-length
      {
	rfile << ",{\"id\": \"readLength\",";
	rfile << "\"title\": \"Read length distribution\",";
	rfile << "\"x\": {\"data\": [{\"values\": [";
	uint32_t lastValidRL = _lastNonZeroIdx(itRg->second.rc.lRc);
	for(uint32_t i = 0; i <= lastValidRL; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << i;
	}
	rfile << "]}], \"axis\": {\"title\": \"Read length\"}},";
	rfile << "\"y\": {\"data\": [{\"values\": [";
	for(uint32_t i = 0; i <= lastValidRL; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << itRg->second.rc.lRc[i];
	}
	rfile << "]}], \"axis\": {\"title\": \"Count\"}}, \"type\": \"line\"}";
      }
      
      // Mean Base Quality
      {
	rfile << ",{\"id\": \"baseQuality\",";	
	rfile << "\"title\": \"Mean base quality distribution\",";
	rfile << "\"x\": {\"data\": [{\"values\": [";
	uint32_t lastValidBQIdx = _lastNonZeroIdxACGTN(itRg->second.rc);
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << i;
	}
	rfile << "]}], \"axis\": {\"title\": \"Read position\"}},";
	rfile << "\"y\": {\"data\": [{\"values\": [";
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  if (i > 0) rfile << ",";
	  uint64_t bcount = itRg->second.rc.aCount[i] + itRg->second.rc.cCount[i] + itRg->second.rc.gCount[i] + itRg->second.rc.tCount[i] + itRg->second.rc.nCount[i];
	  if (bcount > 0) rfile << (double) (itRg->second.rc.bqCount[i]) / (double) (bcount);
	  else rfile << 0;
	}
	rfile << "]}], \"axis\": {\"title\": \"Average base quality\"}}, \"type\": \"line\"}";
      }

      // Mapping quality histogram
      {
	rfile << ",{\"id\": \"mappingQuality\", \"title\": \"Mapping quality distribution\",";
	rfile << "\"x\": {\"data\": [{\"values\": [";
	uint32_t lastValidMQ = _lastNonZeroIdx(itRg->second.qc.qcount);
	for(uint32_t i = 0; i <= lastValidMQ; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << i;
	}
	rfile << "]}], \"axis\": {\"title\": \"Mapping Quality\"}},";
	rfile << "\"y\": {\"data\": [{\"values\": [";
	for(uint32_t i = 0; i <= lastValidMQ; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << itRg->second.qc.qcount[i];
	}
	rfile << "]}], \"axis\": {\"title\": \"Count\"}}, \"type\": \"line\"}";
      }

      // Coverage Histogram
      {
	rfile << ",{\"id\": \"coverageHistogram\", \"title\": \"Coverage histogram\",";
	rfile << "\"x\": {\"data\": [{\"values\": [";
	uint32_t lastValidCO = _lastNonZeroIdx(itRg->second.bc.bpWithCoverage);
	for(uint32_t i = 0; i <= lastValidCO; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << i;
	}
	rfile << "]}], \"axis\": {\"title\": \"Coverage\", \"range\": [0,60]}},";
	rfile << "\"y\": {\"data\": [{\"values\": [";
	for(uint32_t i = 0; i <= lastValidCO; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << itRg->second.bc.bpWithCoverage[i];
	}
	rfile << "]}], \"axis\": {\"title\": \"Count\"}}, \"type\": \"line\"}";
      }

      // Insert Size Histogram
      {
	rfile << ",{\"id\": \"insertSize\", \"title\": \"Insert size histogram\",";
	rfile << "\"x\": {\"data\": [{\"values\": [";
	uint32_t lastValidIS = _lastNonZeroIdxISize(itRg->second.pc);
	for(uint32_t i = 0; i <= lastValidIS; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << i;
	}
	rfile << "]}], \"axis\": {\"title\": \"Insert Size\", \"range\": [0,1000]}},";
	rfile << "\"y\": {\"data\": [";
	rfile << "{\"values\": [";
	for(uint32_t i = 0; i <= lastValidIS; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << itRg->second.pc.fPlus[i];
	}
	rfile << "], \"title\": \"F+\"},";
	rfile << "{\"values\": [";
	for(uint32_t i = 0; i <= lastValidIS; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << itRg->second.pc.fMinus[i];
	}
	rfile << "], \"title\": \"F-\"},";
	rfile << "{\"values\": [";
	for(uint32_t i = 0; i <= lastValidIS; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << itRg->second.pc.rPlus[i];
	}
	rfile << "], \"title\": \"R+\"},";
	rfile << "{\"values\": [";
	for(uint32_t i = 0; i <= lastValidIS; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << itRg->second.pc.rMinus[i];
	}
	rfile << "], \"title\": \"R-\"}], \"axis\": {\"title\": \"Count\"}}, \"type\": \"line\"}";
      }

      // Bed specific data
      if (c.hasRegionFile) {
	// On target rate
	{
	  rfile << ",{\"id\": \"onTarget\",";
	  rfile << "\"title\": \"On-target rate\",";
	  rfile << "\"x\": {\"data\": [{\"values\": [";
	  uint64_t alignedbases = itRg->second.bc.matchCount + itRg->second.bc.mismatchCount;
	  typename BedCounts::TOnTargetMap::const_iterator itOT = be.onTarget.find(itRg->first);
	  for(uint32_t i = 0; i < itOT->second.size(); ++i) {
	    if (i > 0) rfile << ",";
	    rfile << i * be.stepsize;
	  }
	  rfile << "]}], \"axis\": {\"title\": \"Target Extension\"}},";
	  rfile << "\"y\": {\"data\": [{\"values\": [";
	  for(uint32_t i = 0; i < itOT->second.size(); ++i) {
	    if (i > 0) rfile << ",";
	    rfile << (double) itOT->second[i] / (double) alignedbases;
	  }
	  rfile << "]}], \"axis\": {\"title\": \"Fraction on target\"}}, \"type\": \"line\"}";
	}

	// Avg. target coverage
	{
	  rfile << ",{\"id\": \"targetCoverage\",";
	  rfile << "\"title\": \"Targets above coverage threshold\",";
	  rfile << "\"x\": {\"data\": [{\"values\": [";
	  std::vector<double> fracAboveCov;
	  uint32_t maxBedCov = _lastCoverageLevel(be, rf, hdr->n_targets, fracAboveCov);
	  for(uint32_t i = 0; i < maxBedCov; ++i) {
	    if (i > 0) rfile << ",";
	    rfile << i;
	  }
	  rfile << "]}], \"axis\": {\"title\": \"Coverage Level\"}},";
	  rfile << "\"y\": {\"data\": [{\"values\": [";
	  for(uint32_t i = 0; i < maxBedCov; ++i) {
	    if (i > 0) rfile << ", ";
	    rfile << fracAboveCov[i];
	  }
	  rfile << "]}], \"axis\": {\"title\": \"Fraction above coverage\"}}, \"type\": \"line\"}";
	}
      }

      // InDel Size
      {
	rfile << ",{\"id\": \"indelSize\", \"title\": \"InDel Size\",";
	rfile << "\"x\": {\"data\": [{\"values\": [";
	uint32_t lastSize = itRg->second.bc.delSize.size();
	if (itRg->second.bc.insSize.size() > lastSize) lastSize = itRg->second.bc.insSize.size();
	for(uint32_t i = 0; i < lastSize; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << i;
	}
	rfile << "]}], \"axis\": {\"title\": \"InDel Size\", \"range\": [0,10]}},";
	rfile << "\"y\": {\"data\": [";
	rfile << "{\"values\": [";
	for(uint32_t i = 0; i < lastSize; ++i) {
	  if (i > 0) rfile << ",";
	  if (i < itRg->second.bc.delSize.size()) rfile << itRg->second.bc.delSize[i];
	  else rfile << '0';
	}
	rfile << "], \"title\": \"Deletion\"},";
	rfile << "{\"values\": [";
	for(uint32_t i = 0; i < lastSize; ++i) {
	  if (i > 0) rfile << ",";
	  if (i < itRg->second.bc.insSize.size()) rfile << itRg->second.bc.insSize[i];
	  else rfile << '0';
	}
	rfile << "], \"title\": \"Insertion\"}], \"axis\": {\"title\": \"Count\"}}, \"type\": \"line\"}";
      }


      // GC Content
      {
	rfile << ",{\"id\": \"gcContent\", \"title\": \"GC content\",";
	rfile << "\"x\": {\"data\": [{\"values\": [";
	double refTotal = 0;
	double sampleTotal = 0;
	for(uint32_t i = 0; i < 101; ++i) {
	  refTotal += rf.refGcContent[i];
	  sampleTotal += itRg->second.rc.gcContent[i];
	}
	for(uint32_t i = 0; i < 101; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << (double) i / (double) 100;
	}
	rfile << "]}], \"axis\": {\"title\": \"GC fraction\"}},";
	rfile << "\"y\": {\"data\": [";
	rfile << "{\"values\": [";
	for(uint32_t i = 0; i < 101; ++i) {
	  if (i > 0) rfile << ",";
	  double frac = 0;
	  if (refTotal > 0) frac = (double) (rf.refGcContent[i]) / refTotal;
	  rfile << frac;
	}
	rfile << "], \"title\": \"Reference\"},";
	rfile << "{\"values\": [";
	for(uint32_t i = 0; i < 101; ++i) {
	  if (i > 0) rfile << ",";
	  double frac = 0;
	  if (sampleTotal > 0) frac = (double) itRg->second.rc.gcContent[i] / sampleTotal;
	  rfile << frac;
	}
	rfile << "], \"title\": \"Sample\"}], \"axis\": {\"title\": \"Normalized Fraction\"}}, \"type\": \"line\"}";
      }

      // Homopolymer InDel Context
      {
	rfile << ",{\"id\": \"homIndelContext\", \"title\": \"InDel Context\",";
	rfile << "\"x\": {\"data\": [{\"values\": [\"A\",\"C\",\"G\",\"T\",\"N\",\"None\"]";
	rfile << "}], \"axis\": {\"title\": \"Homopolymer\"}},";
	rfile << "\"y\": {\"data\": [";
	rfile << "{\"values\": [";
	for(uint32_t i = 0; i < itRg->second.bc.delHomACGTN.size(); ++i) {
	  if (i > 0) rfile << ",";
	  rfile << itRg->second.bc.delHomACGTN[i];
	}
	rfile << "], \"title\": \"Deletion\"},";
	rfile << "{\"values\": [";	
	for(uint32_t i = 0; i < itRg->second.bc.insHomACGTN.size(); ++i) {
	  if (i > 0) rfile << ",";
	  rfile << itRg->second.bc.insHomACGTN[i];
	}
	rfile << "], \"title\": \"Insertion\"}], \"axis\": {\"title\": \"Count\"}}, \"type\": \"bar\", \"options\": {\"layout\": \"group\"}}";
      }
      
      // Mapping statistics by chromosome
      {
	rfile << ",{\"id\": \"mappingByChromosome\",";
	rfile << "\"title\": \"Mapping statistics by chromosome\",";
	rfile << "\"data\": {\"columns\": [\"Chr\", \"Size\", \"#N\", \"#GC\", \"GC-fraction\", \"mapped\", \"fracTotal\", \"observedExpected\"], \"rows\": [";
	uint64_t totalMappedChr = 0;
	for(uint32_t i = 0; i < itRg->second.rc.mappedchr.size(); ++i) totalMappedChr += itRg->second.rc.mappedchr[i];
	for(uint32_t i = 0; i < itRg->second.rc.mappedchr.size(); ++i) {
	  if (i > 0) rfile << ",";
	  rfile << "[";
	  double frac = 0;
	  if (totalMappedChr > 0) frac = (double) itRg->second.rc.mappedchr[i] / (double) totalMappedChr;
	  double expect = (double) (hdr->target_len[i] - rf.chrGC[i].ncount) / (double) (rf.referencebp - rf.ncount);
	  double obsexprat = frac / expect;
	  rfile << "\"" << hdr->target_name[i] << "\",";
	  rfile << hdr->target_len[i] << ",";
	  rfile << rf.chrGC[i].ncount << ",";
	  rfile << rf.chrGC[i].gccount << ",";
	  double totalBases = hdr->target_len[i] - rf.chrGC[i].ncount;
	  double gcfrac = 0;
	  if (totalBases > 0) gcfrac = (double) rf.chrGC[i].gccount / totalBases;
	  rfile << gcfrac << ",";
	  rfile << itRg->second.rc.mappedchr[i] << ",";
	  rfile << frac << ",";
	  rfile << obsexprat;
	  rfile << "]";
	}
	rfile << "]},";
	rfile << "\"type\": \"table\"}";
      }

      
      rfile << "]}";
    }
    rfile << "]}]}" << std::endl;
    rfile.pop();
  }
 
}

#endif
