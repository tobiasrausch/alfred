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
      rfile << "\"data\": {\"columns\": [\"Sample\", \"Library\", \"#QCFail\", \"QCFailFraction\", \"#DuplicateMarked\", \"DuplicateFraction\", \"#Unmapped\", \"UnmappedFraction\", \"#Mapped\", \"MappedFraction\", \"#MappedRead1\", \"#MappedRead2\", \"RatioMapped2vsMapped1\", \"#MappedForward\", \"MappedForwardFraction\", \"#MappedReverse\", \"MappedReverseFraction\", \"#SecondaryAlignments\", \"SecondaryAlignmentFraction\", \"#SupplementaryAlignments\", \"SupplementaryAlignmentFraction\", \"#SplicedAlignments\", \"SplicedAlignmentFraction\", ";
      rfile << "\"#Pairs\", \"#MappedPairs\", \"MappedPairsFraction\", \"#MappedSameChr\", \"MappedSameChrFraction\", \"#MappedProperPair\", \"MappedProperFraction\", ";
      rfile << "\"#ReferenceBp\", \"#ReferenceNs\", \"#AlignedBases\", \"#MatchedBases\", \"MatchRate\", \"#MismatchedBases\", \"MismatchRate\", \"#DeletionsCigarD\", \"DeletionRate\", \"HomopolymerContextDel\", \"#InsertionsCigarI\", \"InsertionRate\", \"HomopolymerContextIns\", \"#SoftClippedBases\", \"SoftClipRate\", \"#HardClippedBases\", \"HardClipRate\", \"ErrorRate\", ";
      rfile << "\"MedianReadLength\", \"DefaultLibraryLayout\", \"MedianInsertSize\", \"MedianCoverage\", \"SDCoverage\", \"CoveredBp\", \"FractionCovered\", \"BpCov1ToCovNRatio\", \"BpCov1ToCov2Ratio\", \"MedianMAPQ\"";
      if (c.hasRegionFile) {
	rfile << ",";
	rfile << "\"#TotalBedBp\", \"#AlignedBasesInBed\", \"FractionInBed\", \"EnrichmentOverBed\"";
      }
      if (c.isMitagged) {
	rfile << ",";
	rfile << "\"#MItagged\", \"FractionMItagged\", \"#UMIs\"";
      }
      if (c.isHaplotagged) {
	rfile << ",";
	rfile << "\"#HaploTagged\", \"FractionHaploTagged\", \"#PhasedBlocks\", \"N50PhasedBlockLength\"";
      }
      rfile << "], \"rows\": ["; 
      for(typename TRGMap::const_iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
	uint64_t totalReadCount = _totalReadCount(itRg);
	uint64_t mappedCount = itRg->second.rc.mapped1 + itRg->second.rc.mapped2;
	if (itRg != rgMap.begin()) rfile << ",";
	rfile << "[";
	rfile << "\"" << c.sampleName << "\"" << ",";
	rfile << "\"" << itRg->first << "\"" << ",";
	rfile << itRg->second.rc.qcfail << ",";
	rfile << (double) itRg->second.rc.qcfail / (double) totalReadCount << ",";
	rfile << itRg->second.rc.dup << ",";
	rfile << (double) itRg->second.rc.dup / (double) totalReadCount << ",";
	rfile << itRg->second.rc.unmap << ",";
	rfile << (double) itRg->second.rc.unmap / (double) totalReadCount << ",";
	rfile << mappedCount << ",";
	rfile << (double) mappedCount / (double) totalReadCount << ",";
	rfile << itRg->second.rc.mapped1 << ",";
	rfile << itRg->second.rc.mapped2 << ",";
	rfile << (double) itRg->second.rc.mapped2 / (double) itRg->second.rc.mapped1 << ",";
	rfile << itRg->second.rc.forward << ",";
	rfile << (double) itRg->second.rc.forward / (double) mappedCount << ",";
	rfile << itRg->second.rc.reverse << ",";
	rfile << (double) itRg->second.rc.reverse / (double) mappedCount << ",";
	rfile << itRg->second.rc.secondary << ",";
	rfile << (double) itRg->second.rc.secondary / (double) mappedCount << ",";
	rfile << itRg->second.rc.supplementary << ",";
	rfile << (double) itRg->second.rc.supplementary / (double) mappedCount << ",";
	rfile << itRg->second.rc.spliced << ",";
	rfile << (double) itRg->second.rc.spliced / (double) mappedCount << ",";

	// Paired counts
	int64_t paired = itRg->second.pc.paired / 2;
	int64_t mapped = itRg->second.pc.mapped / 2;
	int64_t mappedSameChr = itRg->second.pc.mappedSameChr / 2;
	int64_t mappedProper = itRg->second.pc.mappedProper / 2;
	double mappedpairedfrac = 0;
	if (paired > 0) mappedpairedfrac = (double) mapped / (double) paired;
	double mappedpairedchrfrac = 0;
	if (paired > 0) mappedpairedchrfrac = (double) mappedSameChr / (double) paired;
	double mappedproperfrac = 0;
	if (paired > 0) mappedproperfrac = (double) mappedProper / (double) paired;
	rfile << paired << ",";
	rfile << mapped << ",";
	rfile << mappedpairedfrac << ",";
	rfile << mappedSameChr << ",";
	rfile << mappedpairedchrfrac << ",";
	rfile << mappedProper << ",";
	rfile << mappedproperfrac << ",";

	// Homopolymer Context of InDels
	double insFrac = _homopolymerIndel(itRg->second.bc.insHomACGTN);
	double delFrac = _homopolymerIndel(itRg->second.bc.delHomACGTN);

	// Error rates
	uint64_t alignedbases = itRg->second.bc.matchCount + itRg->second.bc.mismatchCount;
	double errRate = (double) (itRg->second.bc.mismatchCount + itRg->second.bc.delCount + itRg->second.bc.insCount + itRg->second.bc.softClipCount + itRg->second.bc.hardClipCount) / (double) alignedbases;
	rfile << rf.referencebp << ",";
	rfile << rf.ncount << ",";
	rfile << alignedbases << ",";
	rfile << itRg->second.bc.matchCount << ",";
	rfile << (double) itRg->second.bc.matchCount / (double) alignedbases << ",";
	rfile << itRg->second.bc.mismatchCount << ",";
	rfile << (double) itRg->second.bc.mismatchCount / (double) alignedbases << ",";
	rfile << itRg->second.bc.delCount << ",";
	rfile << (double) itRg->second.bc.delCount / (double) alignedbases << ",";
	rfile << delFrac << ",";
	rfile << itRg->second.bc.insCount << ",";
	rfile << (double) itRg->second.bc.insCount / (double) alignedbases << ",";
	rfile << insFrac << ",";
	rfile << itRg->second.bc.softClipCount << ",";
	rfile << (double) itRg->second.bc.softClipCount / (double) alignedbases << ",";
	rfile << itRg->second.bc.hardClipCount << ",";
	rfile << (double) itRg->second.bc.hardClipCount / (double) alignedbases << ",";
	rfile << errRate << ",";

	// Median coverage, read length, etc.
	int32_t deflayout = _defLayout(itRg);
	int32_t medISize = _medISize(itRg, deflayout);

	// Standardized SD of genomic coverage
	double ssdcov = 1000 * sdFromHistogram(itRg->second.bc.bpWithCoverage) / std::sqrt((double) mappedCount);
	double fraccovbp = (double) itRg->second.bc.nd / (double) (rf.referencebp - rf.ncount);
	double pbc1 = (double) itRg->second.bc.n1 / (double) itRg->second.bc.nd;
	double pbc2 = (double) itRg->second.bc.n1 / (double) itRg->second.bc.n2;

	rfile << medianFromHistogram(itRg->second.rc.lRc) << ",";
	rfile << deflayout << ",";
	rfile << medISize << ",";
	rfile <<  medianFromHistogram(itRg->second.bc.bpWithCoverage) << ",";
	rfile << ssdcov << ",";
	rfile << itRg->second.bc.nd << ",";
	rfile << fraccovbp << ",";
	rfile << pbc1 << ",";
	rfile << pbc2 << ",";
	rfile << medianFromHistogram(itRg->second.qc.qcount);

	// Bed metrics
	if (c.hasRegionFile) {
	  uint64_t nonN = rf.referencebp - rf.ncount;
	  typename BedCounts::TOnTargetMap::const_iterator itOT = be.onTarget.find(itRg->first);
	  uint64_t alignedBedBases = itOT->second[0];
	  double fractioninbed = (double) alignedBedBases / (double) alignedbases;
	  double enrichment = fractioninbed / ((double) rf.totalBedSize / (double) nonN);
	  rfile << ",";
	  rfile << rf.totalBedSize << ",";
	  rfile << alignedBedBases << ",";
	  rfile << fractioninbed << ",";
	  rfile << enrichment;
	}
	if (c.isMitagged) {
	  rfile << ",";
	  rfile << itRg->second.rc.mitagged << ",";
	  rfile << (double) itRg->second.rc.mitagged / (double) totalReadCount << ",";
	  rfile << itRg->second.rc.umi.count();
	}
	if (c.isHaplotagged) {
	  int32_t n50ps = n50PhasedBlockLength(itRg->second.rc.brange);
	  rfile << ",";
	  rfile << itRg->second.rc.haplotagged << ",";
	  rfile << (double) itRg->second.rc.haplotagged / (double) totalReadCount << ",";
	  rfile << phasedBlocks(itRg->second.rc.brange) << ",";
	  rfile << n50ps;
	}
	rfile << "]";
      }
      rfile << "]},";
      rfile << "\"type\": \"table\"}";
    }
    rfile << ",";
    rfile << std::endl;
    
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
	  uint64_t bcount = itRg->second.rc.aCount[i] + itRg->second.rc.cCount[i] + itRg->second.rc.gCount[i] + itRg->second.rc.tCount[i] + itRg->second.rc.nCount[i];
	  if (bcount > 0) {
	    if (i > 0) rfile << ",";
	    rfile << (double) itRg->second.rc.aCount[i] / (double) bcount;
	  }
	}
	rfile << "], \"title\": \"A\"},";
	rfile << "{\"values\": [";
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  uint64_t bcount = itRg->second.rc.aCount[i] + itRg->second.rc.cCount[i] + itRg->second.rc.gCount[i] + itRg->second.rc.tCount[i] + itRg->second.rc.nCount[i];
	  if (bcount > 0) {
	    if (i > 0) rfile << ",";
	    rfile << (double) itRg->second.rc.cCount[i] / (double) bcount;
	  }
	}
	rfile << "], \"title\": \"C\"},";
	rfile << "{\"values\": [";
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  uint64_t bcount = itRg->second.rc.aCount[i] + itRg->second.rc.cCount[i] + itRg->second.rc.gCount[i] + itRg->second.rc.tCount[i] + itRg->second.rc.nCount[i];
	  if (bcount > 0) {
	    if (i > 0) rfile << ",";
	    rfile << (double) itRg->second.rc.gCount[i] / (double) bcount;
	  }
	}
	rfile << "], \"title\": \"G\"},";
	rfile << "{\"values\": [";
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  uint64_t bcount = itRg->second.rc.aCount[i] + itRg->second.rc.cCount[i] + itRg->second.rc.gCount[i] + itRg->second.rc.tCount[i] + itRg->second.rc.nCount[i];
	  if (bcount > 0) {
	    if (i > 0) rfile << ",";
	    rfile << (double) itRg->second.rc.tCount[i] / (double) bcount;
	  }
	}
	rfile << "], \"title\": \"T\"},";
	rfile << "{\"values\": [";
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  uint64_t bcount = itRg->second.rc.aCount[i] + itRg->second.rc.cCount[i] + itRg->second.rc.gCount[i] + itRg->second.rc.tCount[i] + itRg->second.rc.nCount[i];
	  if (bcount > 0) {
	    if (i > 0) rfile << ",";
	    rfile << (double) itRg->second.rc.nCount[i] / (double) bcount;
	  }
	}
	rfile << "], \"title\": \"N\"}], \"axis\": {\"title\": \"Base Fraction\"}}, \"type\": \"line\"}";
      }
	
      // Read-length
      {
	uint32_t lastValidRL = _lastNonZeroIdx(itRg->second.rc.lRc);
	float lastFrac = _lastPercentage(itRg->second.rc.lRc, itRg->second.rc.maxReadLength);
	rfile << ",{\"id\": \"readLength\",";
	rfile << "\"title\": \"Read length distribution\",";
	if (lastFrac > 0) {
	  rfile << "\"subtitle\": \"" << lastFrac << "% of all reads >= " << itRg->second.rc.maxReadLength << "bp\",";
	}
	rfile << "\"x\": {\"data\": [{\"values\": [";
	for(uint32_t i = 0; i < lastValidRL; ++i) {
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
	uint32_t lastValidCO = _lastNonZeroIdx(itRg->second.bc.bpWithCoverage);
	float lastFrac = _lastPercentage(itRg->second.bc.bpWithCoverage, itRg->second.bc.maxCoverage);
	rfile << ",{\"id\": \"coverageHistogram\", \"title\": \"Coverage histogram\",";
	if (lastFrac > 0) {
	  rfile << "\"subtitle\": \"" << lastFrac << "% of all bases with >= " << itRg->second.bc.maxCoverage << "x coverage\",";
	}
	rfile << "\"x\": {\"data\": [{\"values\": [";
	for(uint32_t i = 0; i < lastValidCO; ++i) {
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
	uint32_t lastValidIS = _lastNonZeroIdxISize(itRg->second.pc);
	// Only output for PE data
	if (lastValidIS > 0) {
	  rfile << ",{\"id\": \"insertSize\", \"title\": \"Insert size histogram\",";
	  rfile << "\"x\": {\"data\": [{\"values\": [";
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
	  rfile << "]}], \"axis\": {\"title\": \"Fraction on target\", \"range\": [0,1]}}, \"type\": \"line\"}";
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
	  rfile << "]}], \"axis\": {\"title\": \"Fraction above coverage\", \"range\": [0,1]}}, \"type\": \"line\"}";
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
	double beTotal = 0;
	for(uint32_t i = 0; i < 102; ++i) {
	  refTotal += rf.refGcContent[i];
	  sampleTotal += itRg->second.rc.gcContent[i];
	  if (c.hasRegionFile) beTotal += be.bedGcContent[i];
	}
	for(uint32_t i = 0; i < 102; ++i) {
	  if (i > 0) rfile << ",";
	  rfile << (double) i / (double) 101;
	}
	rfile << "]}], \"axis\": {\"title\": \"GC fraction\"}},";
	rfile << "\"y\": {\"data\": [";
	rfile << "{\"values\": [";
	for(uint32_t i = 0; i < 102; ++i) {
	  if (i > 0) rfile << ",";
	  double frac = 0;
	  if (refTotal > 0) frac = (double) (rf.refGcContent[i]) / refTotal;
	  rfile << frac;
	}
	rfile << "], \"title\": \"Reference\"},";
	if (c.hasRegionFile) {
	  rfile << "{\"values\": [";
	  for(uint32_t i = 0; i < 102; ++i) {
	    if (i > 0) rfile << ",";
	    double frac = 0;
	    if (beTotal > 0) frac = (double) be.bedGcContent[i] / beTotal;
	    rfile << frac;
	  }
	  rfile << "], \"title\": \"Target\"},";
	}
	rfile << "{\"values\": [";
	for(uint32_t i = 0; i < 102; ++i) {
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
	bool firstVal = true;
	for(uint32_t i = 0; i < itRg->second.rc.mappedchr.size(); ++i) totalMappedChr += itRg->second.rc.mappedchr[i];
	for(uint32_t i = 0; i < itRg->second.rc.mappedchr.size(); ++i) {
	  if (hdr->target_len[i] > c.minChrLen) {
	    if (!firstVal) rfile << ",";
	    else firstVal = false;
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
