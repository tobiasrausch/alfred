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

#include <nlohmann/json.hpp>

#include "tsv.h"

namespace bamstats
{
  

  template<typename TConfig, typename TRGMap>
  inline void
  qcJsonOut(TConfig const& c, bam_hdr_t const* hdr, TRGMap const& rgMap, BedCounts const& be, ReferenceFeatures const& rf) {
    boost::iostreams::filtering_ostream rfile;
    rfile.push(boost::iostreams::gzip_compressor());
    rfile.push(boost::iostreams::file_sink(c.jsonout.string().c_str(), std::ios_base::out | std::ios_base::binary));

    // Sample information
    rfile << "{\"samples\": [";
    nlohmann::json sp;
    sp["id"] = c.sampleName;

    // Summary Table
    {
      nlohmann::json j;
      j["id"] = "summaryTable";
      j["title"] = "Summary Statistics";
      j["type"] = "table";
      
      // Columns
      j["data"]["columns"] = {"Sample", "Library", "#QCFail", "QCFailFraction", "#DuplicateMarked", "DuplicateFraction", "#Unmapped", "UnmappedFraction", "#Mapped", "MappedFraction", "#MappedRead1", "#MappedRead2", "RatioMapped2vsMapped1", "#MappedForward", "MappedForwardFraction", "#MappedReverse", "MappedReverseFraction", "#SecondaryAlignments", "SecondaryAlignmentFraction", "#SupplementaryAlignments", "SupplementaryAlignmentFraction", "#SplicedAlignments", "SplicedAlignmentFraction", "#Pairs", "#MappedPairs", "MappedPairsFraction", "#MappedSameChr", "MappedSameChrFraction", "#MappedProperPair", "MappedProperFraction", "#ReferenceBp", "#ReferenceNs", "#AlignedBases", "#MatchedBases", "MatchRate", "#MismatchedBases", "MismatchRate", "#DeletionsCigarD", "DeletionRate", "HomopolymerContextDel", "#InsertionsCigarI", "InsertionRate", "HomopolymerContextIns", "#SoftClippedBases", "SoftClipRate", "#HardClippedBases", "HardClipRate", "ErrorRate", "MedianReadLength", "DefaultLibraryLayout", "MedianInsertSize", "MedianCoverage", "SDCoverage", "CoveredBp", "FractionCovered", "BpCov1ToCovNRatio", "BpCov1ToCov2Ratio", "MedianMAPQ"};
      if (c.hasRegionFile) j["data"]["columns"].insert(j["data"]["columns"].end(), {"#TotalBedBp", "#AlignedBasesInBed", "FractionInBed", "EnrichmentOverBed"});
      if (c.isMitagged) j["data"]["columns"].insert(j["data"]["columns"].end(), {"#MItagged", "FractionMItagged", "#UMIs"});
      if (c.isHaplotagged) j["data"]["columns"].insert(j["data"]["columns"].end(), {"#HaploTagged", "FractionHaploTagged", "#PhasedBlocks", "N50PhasedBlockLength"});

      // Rows
      j["data"]["rows"] = nlohmann::json::array();
      for(typename TRGMap::const_iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
	uint64_t totalReadCount = _totalReadCount(itRg);
	uint64_t mappedCount = itRg->second.rc.mapped1 + itRg->second.rc.mapped2;
	
	nlohmann::json row = nlohmann::json::array();
	row.push_back(c.sampleName);
	row.push_back(itRg->first);
	row.push_back(itRg->second.rc.qcfail);
	if (totalReadCount > 0) row.push_back(itRg->second.rc.qcfail / (double) totalReadCount);
	else row.push_back(nullptr);
	row.push_back(itRg->second.rc.dup);
	if (totalReadCount > 0) row.push_back((double) itRg->second.rc.dup / (double) totalReadCount);
	else row.push_back(nullptr);
	row.push_back(itRg->second.rc.unmap);
	if (totalReadCount > 0) row.push_back((double) itRg->second.rc.unmap / (double) totalReadCount);
	else row.push_back(nullptr);
	row.push_back(mappedCount);
	if (totalReadCount > 0) row.push_back((double) mappedCount / (double) totalReadCount);
	else row.push_back(nullptr);
	row.push_back(itRg->second.rc.mapped1);
	row.push_back(itRg->second.rc.mapped2);
	if (itRg->second.rc.mapped1 > 0) row.push_back((double) itRg->second.rc.mapped2 / (double) itRg->second.rc.mapped1);
	else row.push_back(nullptr);
	row.push_back(itRg->second.rc.forward);
	if (mappedCount > 0) row.push_back((double) itRg->second.rc.forward / (double) mappedCount);
	else row.push_back(nullptr);
	row.push_back(itRg->second.rc.reverse);
	if (mappedCount > 0) row.push_back((double) itRg->second.rc.reverse / (double) mappedCount);
	else row.push_back(nullptr);
	row.push_back(itRg->second.rc.secondary);
	if (mappedCount > 0) row.push_back((double) itRg->second.rc.secondary / (double) mappedCount);
	else row.push_back(nullptr);
	row.push_back(itRg->second.rc.supplementary);
	if (mappedCount > 0) row.push_back((double) itRg->second.rc.supplementary / (double) mappedCount);
	else row.push_back(nullptr);
	row.push_back(itRg->second.rc.spliced);
	if (mappedCount > 0) row.push_back((double) itRg->second.rc.spliced / (double) mappedCount);
	else row.push_back(nullptr);

	// Paired counts
	int64_t paired = itRg->second.pc.paired / 2;
	int64_t mapped = itRg->second.pc.mapped / 2;
	int64_t mappedSameChr = itRg->second.pc.mappedSameChr / 2;
	int64_t mappedProper = itRg->second.pc.mappedProper / 2;
	row.push_back(paired);
	row.push_back(mapped);
	if (paired > 0) row.push_back((double) mapped / (double) paired);
	else row.push_back(nullptr);
	row.push_back(mappedSameChr);
	if (paired > 0) row.push_back((double) mappedSameChr / (double) paired);
	else row.push_back(nullptr);
	row.push_back(mappedProper);
	if (paired > 0) row.push_back((double) mappedProper / (double) paired);
	else row.push_back(nullptr);
	
	// Homopolymer Context of InDels
	double insFrac = _homopolymerIndel(itRg->second.bc.insHomACGTN);
	double delFrac = _homopolymerIndel(itRg->second.bc.delHomACGTN);

	// Error rates
	uint64_t alignedbases = itRg->second.bc.matchCount + itRg->second.bc.mismatchCount;
	row.push_back(rf.referencebp);
	row.push_back(rf.ncount);
	row.push_back(alignedbases);
	row.push_back(itRg->second.bc.matchCount);
	if (alignedbases > 0) row.push_back((double) itRg->second.bc.matchCount / (double) alignedbases);
	else row.push_back(nullptr);
	row.push_back(itRg->second.bc.mismatchCount);
	if (alignedbases > 0) row.push_back((double) itRg->second.bc.mismatchCount / (double) alignedbases);
	else row.push_back(nullptr);
	row.push_back(itRg->second.bc.delCount);
	if (alignedbases > 0) row.push_back((double) itRg->second.bc.delCount / (double) alignedbases);
	else row.push_back(nullptr);
	row.push_back(delFrac);
	row.push_back(itRg->second.bc.insCount);
	if (alignedbases > 0) row.push_back((double) itRg->second.bc.insCount / (double) alignedbases);
	else row.push_back(nullptr);
	row.push_back(insFrac);
	row.push_back(itRg->second.bc.softClipCount);
	if (alignedbases > 0) row.push_back((double) itRg->second.bc.softClipCount / (double) alignedbases);
	else row.push_back(nullptr);
	row.push_back(itRg->second.bc.hardClipCount);
	if (alignedbases > 0) row.push_back((double) itRg->second.bc.hardClipCount / (double) alignedbases);
	else row.push_back(nullptr);
	if (alignedbases > 0) row.push_back((double) (itRg->second.bc.mismatchCount + itRg->second.bc.delCount + itRg->second.bc.insCount + itRg->second.bc.softClipCount + itRg->second.bc.hardClipCount) / (double) alignedbases);
	else row.push_back(nullptr);
	
	// Median coverage, read length, standardized SD of genomic coverage, ...
	int32_t deflayout = _defLayout(itRg);
	if (itRg->second.pc.paired) {
	  std::string rlstr = boost::lexical_cast<std::string>(medianFromHistogram(itRg->second.rc.lRc[0])) + ":" + boost::lexical_cast<std::string>(medianFromHistogram(itRg->second.rc.lRc[1]));
	  row.push_back(rlstr);
	  row.push_back(_defLayoutToString(deflayout));
	  row.push_back(_medISize(itRg, deflayout));
	} else {
	  row.push_back(medianFromHistogram(itRg->second.rc.lRc[0]));
	  row.push_back(nullptr);
	  row.push_back(nullptr);
	}
	row.push_back(medianFromHistogram(itRg->second.bc.bpWithCoverage));
	if (mappedCount > 0) row.push_back((double) (1000.0 * sdFromHistogram(itRg->second.bc.bpWithCoverage) / std::sqrt((double) mappedCount)));
	else row.push_back(nullptr);
	row.push_back(itRg->second.bc.nd);
	if ((rf.referencebp - rf.ncount) > 0) row.push_back((double) itRg->second.bc.nd / (double) (rf.referencebp - rf.ncount));
	else row.push_back(nullptr);
	if (itRg->second.bc.nd > 0) row.push_back((double) itRg->second.bc.n1 / (double) itRg->second.bc.nd);
	else row.push_back(nullptr);
	if (itRg->second.bc.n2 > 0) row.push_back((double) itRg->second.bc.n1 / (double) itRg->second.bc.n2);
	else row.push_back(nullptr);
	row.push_back(medianFromHistogram(itRg->second.qc.qcount));
	
	// Bed metrics
	if (c.hasRegionFile) {
	  uint64_t nonN = rf.referencebp - rf.ncount;
	  typename BedCounts::TOnTargetMap::const_iterator itOT = be.onTarget.find(itRg->first);
	  uint64_t alignedBedBases = itOT->second[0];
	  double fractioninbed = (double) alignedBedBases / (double) alignedbases;
	  double bedfraction = ((double) rf.totalBedSize / (double) nonN);
	  row.push_back(rf.totalBedSize);
	  row.push_back(alignedBedBases);
	  if (alignedbases > 0) row.push_back((double) alignedBedBases / (double) alignedbases);
	  else row.push_back(nullptr);
	  if (bedfraction > 0) row.push_back(fractioninbed / bedfraction);
	  else row.push_back(nullptr);
	}
	if (c.isMitagged) {
	  row.push_back(itRg->second.rc.mitagged);
	  if (totalReadCount > 0) row.push_back((double) itRg->second.rc.mitagged / (double) totalReadCount);
	  else row.push_back(nullptr);
	  row.push_back(itRg->second.rc.umi.count());
	}
	if (c.isHaplotagged) {
	  row.push_back(itRg->second.rc.haplotagged);
	  if (totalReadCount > 0) row.push_back((double) itRg->second.rc.haplotagged / (double) totalReadCount);
	  else row.push_back(nullptr);
	  row.push_back(phasedBlocks(itRg->second.rc.brange));
	  row.push_back(n50PhasedBlockLength(itRg->second.rc.brange));
	}
	j["data"]["rows"].push_back(row);
      }
      sp["summary"] = j;
    }
    
    // Read-group information
    sp["readGroups"] = nlohmann::json::array();

    // All read-groups
    for(typename TRGMap::const_iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
      nlohmann::json rg;
      rg["id"] = itRg->first;
      rg["metrics"] = nlohmann::json::array();

      // Base content read1
      {
	nlohmann::json j;
	// Paired-end?
	if (itRg->second.pc.paired) {
	  j["id"] = "baseContentRead1";
	  j["title"] = "Base content distribution read1";
	} else {
	  j["id"] = "baseContent";
	  j["title"] = "Base content distribution";
	}
	uint32_t lastValidBQIdx = _lastNonZeroIdxACGTN(itRg->second.rc, 0);
	j["x"]["data"] = nlohmann::json::array();
	j["x"]["axis"]["title"] = "Position in read";
	{
	  nlohmann::json axisX;
	  nlohmann::json valx = nlohmann::json::array();
	  for(uint32_t i = 0; i <= lastValidBQIdx; ++i) valx.push_back(i);
	  axisX["values"] = valx;
	  j["x"]["data"].push_back(axisX);
	}
	j["y"]["data"] = nlohmann::json::array();
	{
	  nlohmann::json axisY;
	  nlohmann::json valy = nlohmann::json::array();
	  for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	    uint64_t bcount = itRg->second.rc.aCount[0][i] + itRg->second.rc.cCount[0][i] + itRg->second.rc.gCount[0][i] + itRg->second.rc.tCount[0][i] + itRg->second.rc.nCount[0][i];
	    if (bcount > 0) valy.push_back((double) itRg->second.rc.aCount[0][i] / (double) bcount);
	    else valy.push_back(nullptr);
	  }
	  axisY["values"] = valy;
	  axisY["title"] = "A";
	  j["y"]["data"].push_back(axisY);
	}
	{
	  nlohmann::json axisY;
	  nlohmann::json valy = nlohmann::json::array();
	  for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	    uint64_t bcount = itRg->second.rc.aCount[0][i] + itRg->second.rc.cCount[0][i] + itRg->second.rc.gCount[0][i] + itRg->second.rc.tCount[0][i] + itRg->second.rc.nCount[0][i];
	    if (bcount > 0) valy.push_back((double) itRg->second.rc.cCount[0][i] / (double) bcount);
	    else valy.push_back(nullptr);
	  }
	  axisY["values"] = valy;
	  axisY["title"] = "C";
	  j["y"]["data"].push_back(axisY);
	}
	{
	  nlohmann::json axisY;
	  nlohmann::json valy = nlohmann::json::array();
	  for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	    uint64_t bcount = itRg->second.rc.aCount[0][i] + itRg->second.rc.cCount[0][i] + itRg->second.rc.gCount[0][i] + itRg->second.rc.tCount[0][i] + itRg->second.rc.nCount[0][i];
	    if (bcount > 0) valy.push_back((double) itRg->second.rc.gCount[0][i] / (double) bcount);
	    else valy.push_back(nullptr);
	  }
	  axisY["values"] = valy;
	  axisY["title"] = "G";
	  j["y"]["data"].push_back(axisY);
	}
	{
	  nlohmann::json axisY;
	  nlohmann::json valy = nlohmann::json::array();
	  for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	    uint64_t bcount = itRg->second.rc.aCount[0][i] + itRg->second.rc.cCount[0][i] + itRg->second.rc.gCount[0][i] + itRg->second.rc.tCount[0][i] + itRg->second.rc.nCount[0][i];
	    if (bcount > 0) valy.push_back((double) itRg->second.rc.tCount[0][i] / (double) bcount);
	    else valy.push_back(nullptr);
	  }
	  axisY["values"] = valy;
	  axisY["title"] = "T";
	  j["y"]["data"].push_back(axisY);
	}
	{
	  nlohmann::json axisY;
	  nlohmann::json valy = nlohmann::json::array();
	  for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	    uint64_t bcount = itRg->second.rc.aCount[0][i] + itRg->second.rc.cCount[0][i] + itRg->second.rc.gCount[0][i] + itRg->second.rc.tCount[0][i] + itRg->second.rc.nCount[0][i];
	    if (bcount > 0) valy.push_back((double) itRg->second.rc.nCount[0][i] / (double) bcount);
	    else valy.push_back(nullptr);
	  }
	  axisY["values"] = valy;
	  axisY["title"] = "N";
	  j["y"]["data"].push_back(axisY);
	}
	j["y"]["axis"]["title"] = "Base Fraction";
	j["type"] = "line";
	rg["metrics"].push_back(j);
      }

      // Paired-end?
      if (itRg->second.pc.paired) {
	// Base content read2
	{
	  nlohmann::json j;
	  j["id"] = "baseContentRead2";
	  j["title"] = "Base content distribution read2";
	  uint32_t lastValidBQIdx = _lastNonZeroIdxACGTN(itRg->second.rc, 1);
	  j["x"]["data"] = nlohmann::json::array();
	  j["x"]["axis"]["title"] = "Position in read";
	  {
	    nlohmann::json axisX;
	    nlohmann::json valx = nlohmann::json::array();
	    for(uint32_t i = 0; i <= lastValidBQIdx; ++i) valx.push_back(i);
	    axisX["values"] = valx;
	    j["x"]["data"].push_back(axisX);
	  }
	  j["y"]["data"] = nlohmann::json::array();
	  {
	    nlohmann::json axisY;
	    nlohmann::json valy = nlohmann::json::array();
	    for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	      uint64_t bcount = itRg->second.rc.aCount[1][i] + itRg->second.rc.cCount[1][i] + itRg->second.rc.gCount[1][i] + itRg->second.rc.tCount[1][i] + itRg->second.rc.nCount[1][i];
	      if (bcount > 0) valy.push_back((double) itRg->second.rc.aCount[1][i] / (double) bcount);
	      else valy.push_back(nullptr);
	    }
	    axisY["values"] = valy;
	    axisY["title"] = "A";
	    j["y"]["data"].push_back(axisY);
	  }
	  {
	    nlohmann::json axisY;
	    nlohmann::json valy = nlohmann::json::array();
	    for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	      uint64_t bcount = itRg->second.rc.aCount[1][i] + itRg->second.rc.cCount[1][i] + itRg->second.rc.gCount[1][i] + itRg->second.rc.tCount[1][i] + itRg->second.rc.nCount[1][i];
	      if (bcount > 0) valy.push_back((double) itRg->second.rc.cCount[1][i] / (double) bcount);
	      else valy.push_back(nullptr);
	    }
	    axisY["values"] = valy;
	    axisY["title"] = "C";
	    j["y"]["data"].push_back(axisY);
	  }
	  {
	    nlohmann::json axisY;
	    nlohmann::json valy = nlohmann::json::array();
	    for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	      uint64_t bcount = itRg->second.rc.aCount[1][i] + itRg->second.rc.cCount[1][i] + itRg->second.rc.gCount[1][i] + itRg->second.rc.tCount[1][i] + itRg->second.rc.nCount[1][i];
	      if (bcount > 0) valy.push_back((double) itRg->second.rc.gCount[1][i] / (double) bcount);
	      else valy.push_back(nullptr);
	    }
	    axisY["values"] = valy;
	    axisY["title"] = "G";
	    j["y"]["data"].push_back(axisY);
	  }
	  {
	    nlohmann::json axisY;
	    nlohmann::json valy = nlohmann::json::array();
	    for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	      uint64_t bcount = itRg->second.rc.aCount[1][i] + itRg->second.rc.cCount[1][i] + itRg->second.rc.gCount[1][i] + itRg->second.rc.tCount[1][i] + itRg->second.rc.nCount[1][i];
	      if (bcount > 0) valy.push_back((double) itRg->second.rc.tCount[1][i] / (double) bcount);
	      else valy.push_back(nullptr);
	    }
	    axisY["values"] = valy;
	    axisY["title"] = "T";
	    j["y"]["data"].push_back(axisY);
	  }
	  {
	    nlohmann::json axisY;
	    nlohmann::json valy = nlohmann::json::array();
	    for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	      uint64_t bcount = itRg->second.rc.aCount[1][i] + itRg->second.rc.cCount[1][i] + itRg->second.rc.gCount[1][i] + itRg->second.rc.tCount[1][i] + itRg->second.rc.nCount[1][i];
	      if (bcount > 0) valy.push_back((double) itRg->second.rc.nCount[1][i] / (double) bcount);
	      else valy.push_back(nullptr);
	    }
	    axisY["values"] = valy;
	    axisY["title"] = "N";
	    j["y"]["data"].push_back(axisY);
	  }
	  j["y"]["axis"]["title"] = "Base Fraction";
	  j["type"] = "line";
	  rg["metrics"].push_back(j);
	}
      }
	
      // Read-length
      {
	float lastFrac1 = _lastPercentage(itRg->second.rc.lRc[0], itRg->second.rc.maxReadLength);
	float lastFrac2 = _lastPercentage(itRg->second.rc.lRc[1], itRg->second.rc.maxReadLength);
	nlohmann::json j;
	j["id"] = "readLength";
	j["title"] = "Read length distribution";
	if (itRg->second.pc.paired) {
	  if ((lastFrac1 > 0) || (lastFrac2 > 0)) {
	    std::string rlstr = boost::lexical_cast<std::string>(lastFrac1) + "%, " + boost::lexical_cast<std::string>(lastFrac2) + "% of all reads >= " + boost::lexical_cast<std::string>(itRg->second.rc.maxReadLength) + "bp";
	    j["subtitle"] = rlstr;
	  }
	} else {
	  if (lastFrac1 > 0) {
	    std::string rlstr = boost::lexical_cast<std::string>(lastFrac1) + "% of all reads >= " + boost::lexical_cast<std::string>(itRg->second.rc.maxReadLength) + "bp";
	    j["subtitle"] = rlstr;
	  }
	}
	j["x"]["data"] = nlohmann::json::array();
	j["x"]["axis"]["title"] = "Read length";
	{
	  nlohmann::json axisX;
	  nlohmann::json valx = nlohmann::json::array();
	  uint32_t lastValidRL = _lastNonZeroIdx(itRg->second.rc.lRc[0], itRg->second.rc.maxReadLength);
	  for(uint32_t i = 0; i <= lastValidRL; ++i) valx.push_back(i);
	  axisX["values"] = valx;
	  j["x"]["data"].push_back(axisX);
	}
	// Paired-end?
	if (itRg->second.pc.paired) {
	  nlohmann::json axisX;
	  nlohmann::json valx = nlohmann::json::array();
	  uint32_t lastValidRL = _lastNonZeroIdx(itRg->second.rc.lRc[1], itRg->second.rc.maxReadLength);
	  for(uint32_t i = 0; i <= lastValidRL; ++i) valx.push_back(i);
	  axisX["values"] = valx;
	  j["x"]["data"].push_back(axisX);
	}
	j["y"]["data"] = nlohmann::json::array();
	{
	  nlohmann::json axisY;
	  nlohmann::json valy = nlohmann::json::array();
	  uint32_t lastValidRL = _lastNonZeroIdx(itRg->second.rc.lRc[0], itRg->second.rc.maxReadLength);
	  for(uint32_t i = 0; i <= lastValidRL; ++i) valy.push_back(itRg->second.rc.lRc[0][i]);
	  axisY["values"] = valy;
	  if (itRg->second.pc.paired) axisY["title"] = "Read1";
	  j["y"]["data"].push_back(axisY);
	}
	// Paired-end?
	if (itRg->second.pc.paired) {
	  nlohmann::json axisY;
	  nlohmann::json valy = nlohmann::json::array();
	  uint32_t lastValidRL = _lastNonZeroIdx(itRg->second.rc.lRc[1], itRg->second.rc.maxReadLength);
	  for(uint32_t i = 0; i <= lastValidRL; ++i) valy.push_back(itRg->second.rc.lRc[1][i]);
	  axisY["values"] = valy;
	  axisY["title"] = "Read2";
	  j["y"]["data"].push_back(axisY);
	}
	j["y"]["axis"]["title"] = "Count";
	j["type"] = "bar";
	j["options"]["layout"] = "group";
	rg["metrics"].push_back(j);
      }
      
      // Mean Base Quality
      {
	nlohmann::json j;
	j["id"] = "baseQuality";
	j["title"] = "Mean base quality distribution";
	j["x"]["data"] = nlohmann::json::array();
	j["x"]["axis"]["title"] = "Read position";
	{
	  nlohmann::json axisX;
	  nlohmann::json valx = nlohmann::json::array();
	  uint32_t lastValidBQIdx = _lastNonZeroIdxACGTN(itRg->second.rc, 0);
	  for(uint32_t i = 0; i <= lastValidBQIdx; ++i) valx.push_back(i);
	  axisX["values"] = valx;
	  j["x"]["data"].push_back(axisX);
	}	  
	// Paired-end?
	if (itRg->second.pc.paired) {
	  // Multiple x-axis (for different read length)
	  nlohmann::json axisX;
	  nlohmann::json valx = nlohmann::json::array();
	  uint32_t lastValidBQIdx = _lastNonZeroIdxACGTN(itRg->second.rc, 1);
	  for(uint32_t i = 0; i <= lastValidBQIdx; ++i) valx.push_back(i);
	  axisX["values"] = valx;
	  j["x"]["data"].push_back(axisX);
	}
	j["y"]["data"] = nlohmann::json::array();
	{
	  nlohmann::json axisY;
	  nlohmann::json valy = nlohmann::json::array();
	  uint32_t lastValidBQIdx = _lastNonZeroIdxACGTN(itRg->second.rc, 0);
	  for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	    uint64_t bcount = itRg->second.rc.aCount[0][i] + itRg->second.rc.cCount[0][i] + itRg->second.rc.gCount[0][i] + itRg->second.rc.tCount[0][i] + itRg->second.rc.nCount[0][i];
	    if (bcount > 0) valy.push_back((double) (itRg->second.rc.bqCount[0][i]) / (double) (bcount));
	    else valy.push_back(nullptr);
	  }
	  axisY["values"] = valy;
	  if (itRg->second.pc.paired) axisY["title"] = "Read1";
	  j["y"]["data"].push_back(axisY);
	}
	// Paired-end?
	if (itRg->second.pc.paired) {
	  nlohmann::json axisY;
	  nlohmann::json valy = nlohmann::json::array();
	  uint32_t lastValidBQIdx = _lastNonZeroIdxACGTN(itRg->second.rc, 1);
	  for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	    uint64_t bcount = itRg->second.rc.aCount[1][i] + itRg->second.rc.cCount[1][i] + itRg->second.rc.gCount[1][i] + itRg->second.rc.tCount[1][i] + itRg->second.rc.nCount[1][i];
	    if (bcount > 0) valy.push_back((double) (itRg->second.rc.bqCount[1][i]) / (double) (bcount));
	    else valy.push_back(nullptr);
	  }
	  axisY["values"] = valy;
	  axisY["title"] = "Read2";
	  j["y"]["data"].push_back(axisY);
	}
	j["y"]["axis"]["title"] = "Average base quality";
	j["type"] = "line";
	rg["metrics"].push_back(j);
      }

      // Mapping quality histogram
      {
	nlohmann::json j;
	j["id"] = "mappingQuality";
	j["title"] = "Mapping quality distribution";
	j["x"]["data"] = nlohmann::json::array();
	j["x"]["axis"]["title"] = "Mapping Quality";
	uint32_t lastValidMQ = _lastNonZeroIdx(itRg->second.qc.qcount);
	{
	  nlohmann::json axisX;
	  nlohmann::json valx = nlohmann::json::array();
	  for(uint32_t i = 0; i <= lastValidMQ; ++i) valx.push_back(i);
	  axisX["values"] = valx;
	  j["x"]["data"].push_back(axisX);
	}
	j["y"]["data"] = nlohmann::json::array();
	{
	  nlohmann::json axisY;
	  nlohmann::json valy = nlohmann::json::array();
	  for(uint32_t i = 0; i <= lastValidMQ; ++i) valy.push_back(itRg->second.qc.qcount[i]);
	  axisY["values"] = valy;
	  j["y"]["data"].push_back(axisY);
	}
	j["y"]["axis"]["title"] = "Count";
	j["type"] = "bar";
	j["options"]["layout"] = "group";
	rg["metrics"].push_back(j);
      }

      // Coverage Histogram
      {
	uint32_t lastValidCO = _lastNonZeroIdx(itRg->second.bc.bpWithCoverage, itRg->second.bc.maxCoverage);
	float lastFrac = _lastPercentage(itRg->second.bc.bpWithCoverage, itRg->second.bc.maxCoverage);
	nlohmann::json j;
	j["id"] = "coverageHistogram";
	j["title"] = "Coverage histogram";
	if (lastFrac > 0) {
	  std::string rlstr = boost::lexical_cast<std::string>(lastFrac) + "% of all bases with >= " + boost::lexical_cast<std::string>(itRg->second.bc.maxCoverage) + "x coverage";
	  j["subtitle"] = rlstr;
	}
	j["x"]["data"] = nlohmann::json::array();
	j["x"]["axis"]["title"] = "Coverage";
	j["x"]["axis"]["range"] = {1, 60};
	{
	  nlohmann::json axisX;
	  nlohmann::json valx = nlohmann::json::array();
	  for(uint32_t i = 0; i <= lastValidCO; ++i) valx.push_back(i);
	  axisX["values"] = valx;
	  j["x"]["data"].push_back(axisX);
	}
	j["y"]["data"] = nlohmann::json::array();
	{
	  nlohmann::json axisY;
	  nlohmann::json valy = nlohmann::json::array();
	  for(uint32_t i = 0; i <= lastValidCO; ++i) valy.push_back(itRg->second.bc.bpWithCoverage[i]);
	  axisY["values"] = valy;
	  j["y"]["data"].push_back(axisY);
	}
	j["y"]["axis"]["title"] = "Count";
	j["type"] = "line";
	rg["metrics"].push_back(j);
      }

      // Insert Size Histogram
      {
	uint32_t lastValidIS = _lastNonZeroIdxISize(itRg->second.pc, itRg->second.pc.maxInsertSize);	
	// Only output for PE data
	if (lastValidIS > 0) {
	  int32_t deflayout = _defLayout(itRg);
	  float lastFrac = 0;
	  if (deflayout == 0) lastFrac = _lastPercentage(itRg->second.pc.fPlus, itRg->second.pc.maxInsertSize);
	  else if (deflayout == 1) lastFrac = _lastPercentage(itRg->second.pc.fMinus, itRg->second.pc.maxInsertSize);
	  else if (deflayout == 2) lastFrac = _lastPercentage(itRg->second.pc.rPlus, itRg->second.pc.maxInsertSize);
	  else if (deflayout == 3) lastFrac = _lastPercentage(itRg->second.pc.rMinus, itRg->second.pc.maxInsertSize);
	  else lastFrac = 0;
	  nlohmann::json j;
	  j["id"] = "insertSize";
	  j["title"] = "Insert size histogram";
	  j["x"]["data"] = nlohmann::json::array();
	  j["x"]["axis"]["title"] = "Insert Size";
	  j["x"]["axis"]["range"] = {0, 1000};
	  if (lastFrac > 0) {
	    std::string rlstr = boost::lexical_cast<std::string>(lastFrac) + "% of all ";
	    rlstr += _defLayoutToString(deflayout);
	    rlstr += " pairs span >= " + boost::lexical_cast<std::string>(itRg->second.pc.maxInsertSize) + "bp";
	    j["subtitle"] = rlstr;
	  }
	  {
	    nlohmann::json axisX;
	    nlohmann::json valx = nlohmann::json::array();
	    for(uint32_t i = 0; i <= lastValidIS; ++i) valx.push_back(i);
	    axisX["values"] = valx;
	    j["x"]["data"].push_back(axisX);
	  }
	  j["y"]["data"] = nlohmann::json::array();
	  {
	    nlohmann::json axisY;
	    nlohmann::json valy = nlohmann::json::array();
	    for(uint32_t i = 0; i <= lastValidIS; ++i) valy.push_back(itRg->second.pc.fPlus[i]);
	    axisY["values"] = valy;
	    axisY["title"] = "F+";
	    j["y"]["data"].push_back(axisY);
	  }
	  {
	    nlohmann::json axisY;
	    nlohmann::json valy = nlohmann::json::array();
	    for(uint32_t i = 0; i <= lastValidIS; ++i) valy.push_back(itRg->second.pc.fMinus[i]);
	    axisY["values"] = valy;
	    axisY["title"] = "F-";
	    j["y"]["data"].push_back(axisY);
	  }
	  {
	    nlohmann::json axisY;
	    nlohmann::json valy = nlohmann::json::array();
	    for(uint32_t i = 0; i <= lastValidIS; ++i) valy.push_back(itRg->second.pc.rPlus[i]);
	    axisY["values"] = valy;
	    axisY["title"] = "R+";
	    j["y"]["data"].push_back(axisY);
	  }
	  {
	    nlohmann::json axisY;
	    nlohmann::json valy = nlohmann::json::array();
	    for(uint32_t i = 0; i <= lastValidIS; ++i) valy.push_back(itRg->second.pc.rMinus[i]);
	    axisY["values"] = valy;
	    axisY["title"] = "R-";
	    j["y"]["data"].push_back(axisY);
	  }
	  j["y"]["axis"]["title"] = "Count";
	  j["type"] = "line";
	  rg["metrics"].push_back(j);
	}
      }

      // Bed specific data
      if (c.hasRegionFile) {
	// On target rate
	{
	  nlohmann::json j;
	  j["id"] = "onTarget";
	  j["title"] = "On-target rate";
	  j["x"]["data"] = nlohmann::json::array();
	  j["x"]["axis"]["title"] = "Target Extension";
	  uint64_t alignedbases = itRg->second.bc.matchCount + itRg->second.bc.mismatchCount;
	  typename BedCounts::TOnTargetMap::const_iterator itOT = be.onTarget.find(itRg->first);
	  {
	    nlohmann::json axisX;
	    nlohmann::json valx = nlohmann::json::array();
	    for(uint32_t i = 0; i < itOT->second.size(); ++i) valx.push_back(i * be.stepsize);
	    axisX["values"] = valx;
	    j["x"]["data"].push_back(axisX);
	  }
	  j["y"]["data"] = nlohmann::json::array();
	  {
	    nlohmann::json axisY;
	    nlohmann::json valy = nlohmann::json::array();
	    for(uint32_t i = 0; i < itOT->second.size(); ++i) {
	      if (alignedbases > 0) valy.push_back((double) itOT->second[i] / (double) alignedbases);
	      else valy.push_back(nullptr);
	    }
	    axisY["values"] = valy;
	    j["y"]["data"].push_back(axisY);
	  }
	  j["y"]["axis"]["title"] = "Fraction on target";
	  j["y"]["axis"]["range"] = {0, 1};
	  j["type"] = "line";
	  rg["metrics"].push_back(j);
	}

	// Avg. target coverage
	{
	  nlohmann::json j;
	  j["id"] = "targetCoverage";
	  j["title"] = "Targets above coverage threshold";
	  j["x"]["data"] = nlohmann::json::array();
	  j["x"]["axis"]["title"] = "Coverage Level";
	  std::vector<double> fracAboveCov;
	  uint32_t maxBedCov = _lastCoverageLevel(be, rf, hdr->n_targets, fracAboveCov);
	  {
	    nlohmann::json axisX;
	    nlohmann::json valx = nlohmann::json::array();
	    for(uint32_t i = 0; i < maxBedCov; ++i) valx.push_back(i);
	    axisX["values"] = valx;
	    j["x"]["data"].push_back(axisX);
	  }
	  j["y"]["data"] = nlohmann::json::array();
	  {
	    nlohmann::json axisY;
	    nlohmann::json valy = nlohmann::json::array();
	    for(uint32_t i = 0; i < maxBedCov; ++i) valy.push_back(fracAboveCov[i]);
	    axisY["values"] = valy;
	    j["y"]["data"].push_back(axisY);
	  }
	  j["y"]["axis"]["title"] = "Fraction above coverage";
	  j["y"]["axis"]["range"] = {0, 1};
	  j["type"] = "line";
	  rg["metrics"].push_back(j);
	}
      }

      // InDel Size
      {
	nlohmann::json j;
	j["id"] = "indelSize";
	j["title"] = "InDel Size";
	j["x"]["data"] = nlohmann::json::array();
	j["x"]["axis"]["title"] = "InDel Size";
	j["x"]["axis"]["range"] = {0, 10};
	uint32_t lastSize = itRg->second.bc.delSize.size();
	if (itRg->second.bc.insSize.size() > lastSize) lastSize = itRg->second.bc.insSize.size();
	{
	  nlohmann::json axisX;
	  nlohmann::json valx = nlohmann::json::array();
	  for(uint32_t i = 0; i < lastSize; ++i) valx.push_back(i);
	  axisX["values"] = valx;
	  j["x"]["data"].push_back(axisX);
	}
	j["y"]["data"] = nlohmann::json::array();
	{
	  nlohmann::json axisY;
	  nlohmann::json valy = nlohmann::json::array();
	  for(uint32_t i = 0; i < lastSize; ++i) {
	    if (i < itRg->second.bc.delSize.size()) valy.push_back(itRg->second.bc.delSize[i]);
	    else valy.push_back(0);
	  }
	  axisY["values"] = valy;
	  axisY["title"] = "Deletion";
	  j["y"]["data"].push_back(axisY);
	}
	{
	  nlohmann::json axisY;
	  nlohmann::json valy = nlohmann::json::array();
	  for(uint32_t i = 0; i < lastSize; ++i) {
	    if (i < itRg->second.bc.insSize.size()) valy.push_back(itRg->second.bc.insSize[i]);
	    else valy.push_back(0);
	  }
	  axisY["values"] = valy;
	  axisY["title"] = "Insertion";
	  j["y"]["data"].push_back(axisY);
	}
	j["y"]["axis"]["title"] = "Count";
	j["type"] = "bar";
	j["options"]["layout"] = "group";
	rg["metrics"].push_back(j);
      }

      // GC Content
      {
	nlohmann::json j;
	j["id"] = "gcContent";
	j["title"] = "GC content";
	j["x"]["data"] = nlohmann::json::array();
	j["x"]["axis"]["title"] = "GC fraction";
	double refTotal = 0;
	double sampleTotal = 0;
	double beTotal = 0;
	for(uint32_t i = 0; i < 102; ++i) {
	  refTotal += rf.refGcContent[i];
	  sampleTotal += itRg->second.rc.gcContent[i];
	  if (c.hasRegionFile) beTotal += be.bedGcContent[i];
	}
	{
	  nlohmann::json axisX;
	  nlohmann::json valx = nlohmann::json::array();
	  for(uint32_t i = 0; i < 102; ++i) valx.push_back((double) i / (double) 101);
	  axisX["values"] = valx;
	  j["x"]["data"].push_back(axisX);
	}
	j["y"]["data"] = nlohmann::json::array();
	{
	  nlohmann::json axisY;
	  nlohmann::json valy = nlohmann::json::array();
	  for(uint32_t i = 0; i < 102; ++i) {
	    if (refTotal > 0) valy.push_back((double) (rf.refGcContent[i]) / (double) refTotal);
	    else valy.push_back(nullptr);
	  }
	  axisY["values"] = valy;
	  axisY["title"] = "Reference";
	  j["y"]["data"].push_back(axisY);
	}
	if (c.hasRegionFile) {
	  nlohmann::json axisY;
	  nlohmann::json valy = nlohmann::json::array();
	  for(uint32_t i = 0; i < 102; ++i) {
	    if (beTotal > 0) valy.push_back((double) be.bedGcContent[i] / (double) beTotal);
	    else valy.push_back(nullptr);
	  }
	  axisY["values"] = valy;
	  axisY["title"] = "Target";
	  j["y"]["data"].push_back(axisY);
	}
	{
	  nlohmann::json axisY;
	  nlohmann::json valy = nlohmann::json::array();
	  for(uint32_t i = 0; i < 102; ++i) {
	    if (sampleTotal > 0) valy.push_back((double) itRg->second.rc.gcContent[i] / (double) sampleTotal);
	    else valy.push_back(nullptr);
	  }
	  axisY["values"] = valy;
	  axisY["title"] = "Sample";
	  j["y"]["data"].push_back(axisY);
	}
	j["y"]["axis"]["title"] = "Normalized Fraction";
	j["type"] = "line";
	rg["metrics"].push_back(j);
      }

      // Homopolymer InDel Context
      {
	nlohmann::json j;
	j["id"] = "homIndelContext";
	j["title"] = "InDel Context";
	j["x"]["data"] = nlohmann::json::array();
	nlohmann::json axisX;
	axisX["values"] = {"A", "C", "G", "T", "N", "None"};
	j["x"]["data"].push_back(axisX);
	j["x"]["axis"]["title"] = "Homopolymer";
	{
	  nlohmann::json axisY;
	  nlohmann::json valy = nlohmann::json::array();	
	  for(uint32_t i = 0; i < itRg->second.bc.delHomACGTN.size(); ++i) valy.push_back(itRg->second.bc.delHomACGTN[i]);
	  axisY["values"] = valy;
	  axisY["title"] = "Deletion";
	  j["y"]["data"].push_back(axisY);
	}
	{
	  nlohmann::json axisY;
	  nlohmann::json valy = nlohmann::json::array();
	  for(uint32_t i = 0; i < itRg->second.bc.insHomACGTN.size(); ++i) valy.push_back(itRg->second.bc.insHomACGTN[i]);
	  axisY["values"] = valy;
	  axisY["title"] = "Insertion";
	  j["y"]["data"].push_back(axisY);
	}
	j["y"]["axis"]["title"] = "Count";
	j["type"] = "bar";
	j["options"]["layout"] = "group";
	rg["metrics"].push_back(j);
      }
      
      // Mapping statistics by chromosome
      {
	nlohmann::json j;
	j["id"] = "mappingByChromosome";
	j["title"] = "Mapping statistics by chromosome";
	j["data"]["columns"] = {"Chr", "Size", "#N", "#GC", "GC-fraction", "mapped", "fracTotal", "observedExpected"};
	uint64_t totalMappedChr = 0;
	for(uint32_t i = 0; i < itRg->second.rc.mappedchr.size(); ++i) totalMappedChr += itRg->second.rc.mappedchr[i];
	j["data"]["rows"] = nlohmann::json::array();
	for(uint32_t i = 0; i < itRg->second.rc.mappedchr.size(); ++i) {
	  if (hdr->target_len[i] > c.minChrLen) {	    	
	    nlohmann::json row = nlohmann::json::array();
	    row.push_back(hdr->target_name[i]);
	    row.push_back(hdr->target_len[i]);
	    row.push_back(rf.chrGC[i].ncount);
	    row.push_back(rf.chrGC[i].gccount);
	    double frac = 0;
	    if (totalMappedChr > 0) frac = (double) itRg->second.rc.mappedchr[i] / (double) totalMappedChr;
	    double expect = 0;
	    if ((rf.referencebp - rf.ncount) > 0) expect = (double) (hdr->target_len[i] - rf.chrGC[i].ncount) / (double) (rf.referencebp - rf.ncount);
	    double totalBases = hdr->target_len[i] - rf.chrGC[i].ncount;
	    if (totalBases > 0) row.push_back((double) rf.chrGC[i].gccount / totalBases);
	    else row.push_back(nullptr);
	    row.push_back(itRg->second.rc.mappedchr[i]);
	    if (totalMappedChr > 0) row.push_back(frac);
	    else row.push_back(nullptr);
	    if ((totalMappedChr > 0) && (expect > 0)) row.push_back(frac / expect);
	    else row.push_back(nullptr);
	    j["data"]["rows"].push_back(row);
	  }
	}
	j["type"] = "table";
	rg["metrics"].push_back(j);
      }
      sp["readGroups"].push_back(rg);
    }
    rfile << sp.dump();
    rfile << "]}" << std::endl;
    rfile.pop();
  }
 
}

#endif
