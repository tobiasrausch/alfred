#ifndef TSV_H
#define TSV_H

#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include "qcstruct.h"

namespace bamstats
{

  template<typename TVector>
  inline double
  _homopolymerIndel(TVector const& vec) {
    double vecTotal = 0;
    for(uint32_t i = 0; i < vec.size(); ++i) vecTotal += vec[i];
    double vecFrac = 0;
    if (vecTotal > 0) {
      for(uint32_t i = 0; i < vec.size() - 1; ++i) vecFrac += vec[i];
      vecFrac /= vecTotal;
    }
    return vecFrac;
  }

  template<typename TRgMapIterator>
  inline int32_t
  _medISize(TRgMapIterator const& itRg, int32_t const deflayout) {
    int32_t medISize = 0;
    switch(deflayout) {
    case 0:
      medISize = medianFromHistogram(itRg->second.pc.fPlus);
      break;
    case 1:
      medISize = medianFromHistogram(itRg->second.pc.fMinus);
      break;
    case 2:
      medISize = medianFromHistogram(itRg->second.pc.rPlus);
      break;
    case 3:
      medISize = medianFromHistogram(itRg->second.pc.rMinus);
      break;
    default:
      break;
    }
    return medISize;
  }

  template<typename TVector>
  inline uint32_t
  _lastNonZeroIdx(TVector const& vec, uint32_t const lastIdx) {
    uint32_t lastNonZeroIdx = 0;
    for(uint32_t i = 0; ((i < vec.size()) && (i < lastIdx)); ++i)
      if (vec[i] > 0) lastNonZeroIdx = i;
    return lastNonZeroIdx;
  }

  template<typename TVector>
  inline uint32_t
  _lastNonZeroIdx(TVector const& vec) {
    return _lastNonZeroIdx(vec, vec.size());
  }

  template<typename TVector>
  inline float
  _lastPercentage(TVector const& vec, uint32_t const lastIdx) {
    float cumsum = 0;
    for(uint32_t i = 0; i <= lastIdx; ++i) cumsum += vec[i];
    return ((float) (vec[lastIdx]) * 100.0) / cumsum;
  }

  template<typename TVector>
  inline float
  _lastPercentage(TVector const& vec) {
    return _lastPercentage(vec, vec.size() - 1);
  }

  template<typename TVector>
  inline uint32_t
  _lastNonZeroIdxACGTN(TVector const& vec, uint32_t const read12) {
    uint32_t lastNonZeroIdx = 0;
    for(uint32_t i = 0; i < vec.nCount[read12].size(); ++i) {
      uint64_t bcount = vec.aCount[read12][i] + vec.cCount[read12][i] + vec.gCount[read12][i] + vec.tCount[read12][i] + vec.nCount[read12][i];
      if (bcount > 0) lastNonZeroIdx = i;
    }
    return lastNonZeroIdx;
  }

  template<typename TVector>
  inline uint32_t
  _lastNonZeroIdxISize(TVector const& vec, uint32_t const lastIdx) {
    uint32_t maxFPlus = _lastNonZeroIdx(vec.fPlus, lastIdx);
    uint32_t maxFMinus = _lastNonZeroIdx(vec.fMinus, lastIdx);
    uint32_t maxRPlus = _lastNonZeroIdx(vec.rPlus, lastIdx);
    uint32_t maxRMinus = _lastNonZeroIdx(vec.rMinus, lastIdx);
    uint32_t lastValidISIdx = maxFPlus;
    if (maxFMinus > lastValidISIdx) lastValidISIdx = maxFMinus;
    if (maxRPlus > lastValidISIdx) lastValidISIdx = maxRPlus;
    if (maxRMinus > lastValidISIdx) lastValidISIdx = maxRMinus;
    return lastValidISIdx;
  }

  template<typename TVector>
  inline uint32_t
  _lastNonZeroIdxISize(TVector const& vec) {
    return _lastNonZeroIdxISize(vec, vec.fPlus.size());
  }
  
  template<typename TVector>
  inline uint32_t
  _lastCoverageLevel(BedCounts const& be, ReferenceFeatures const& rf, uint32_t const nchr, TVector& fracAboveCov) {
    typedef typename TVector::value_type TValue;
    
    uint32_t lastLevel = 0;
    for(uint32_t level = 0; level<200000; ++level) {
      uint32_t aboveLevel = 0;
      uint32_t totalTargets = 0;
      for(uint32_t refIndex = 0; refIndex < nchr; ++refIndex) {
	for(typename BedCounts::TRgBpMap::const_iterator itChr = be.gCov[refIndex].begin(); itChr != be.gCov[refIndex].end(); ++itChr) {
	  for(uint32_t i = 0; i < rf.gRegions[refIndex].size(); ++i) {
	    ++totalTargets;
	    if (itChr->second[i] > level) ++aboveLevel;
	  }
	}
      }
      TValue frac = (TValue) aboveLevel / (TValue) totalTargets;
      fracAboveCov.push_back(frac);
      lastLevel = level + 1;
      if (frac < 0.01) break;
    }
    return lastLevel;
  }

  template<typename TRgMapIterator>
  inline uint64_t
  _totalReadCount(TRgMapIterator const& itRg) {
    return itRg->second.rc.qcfail + itRg->second.rc.dup + itRg->second.rc.unmap + itRg->second.rc.mapped1 + itRg->second.rc.mapped2;
  }

  template<typename TRgMapIterator>
  inline int32_t
  _defLayout(TRgMapIterator const& itRg) {
    int32_t deflayout = 0;
    int32_t maxcount = itRg->second.pc.orient[0];
    for(int32_t i = 1; i<4; ++i) {
      if (itRg->second.pc.orient[i] > maxcount) {
	maxcount = itRg->second.pc.orient[i];
	deflayout = i;
      }
    }
    return deflayout;
  }

  inline std::string
  _defLayoutToString(int32_t const dl) {
    if (dl == 0) return "F+";
    else if (dl == 1) return "F-";
    else if (dl == 2) return "R+";
    else if (dl == 3) return "R-";
    else return "NA";
  }
  
  template<typename TConfig, typename TRGMap>
  inline void
  qcTsvOut(TConfig const& c, bam_hdr_t const* hdr, TRGMap const& rgMap, BedCounts const& be, ReferenceFeatures const& rf) {
    // Outfile
    boost::iostreams::filtering_ostream rcfile;
    rcfile.push(boost::iostreams::gzip_compressor());
    rcfile.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));

    // Output header
    rcfile << "# This file was produced by alfred v" << alfredVersionNumber << "." << std::endl;
    rcfile << "# This file contains alignment statistics for all reads." << std::endl;

    // Output metrics
    rcfile << "# Alignment summary metrics (ME)." << std::endl;
    rcfile << "# Use `zgrep ^ME <outfile> | cut -f 2-` to extract this part." << std::endl;
    rcfile << "ME\tSample\tLibrary\t#QCFail\tQCFailFraction\t#DuplicateMarked\tDuplicateFraction\t#Unmapped\tUnmappedFraction\t#Mapped\tMappedFraction\t#MappedRead1\t#MappedRead2\tRatioMapped2vsMapped1\t#MappedForward\tMappedForwardFraction\t#MappedReverse\tMappedReverseFraction\t#SecondaryAlignments\tSecondaryAlignmentFraction\t#SupplementaryAlignments\tSupplementaryAlignmentFraction\t#SplicedAlignments\tSplicedAlignmentFraction" << "\t";
    rcfile << "#Pairs\t#MappedPairs\tMappedPairsFraction\t#MappedSameChr\tMappedSameChrFraction\t#MappedProperPair\tMappedProperFraction" << "\t";
    rcfile << "#ReferenceBp\t#ReferenceNs\t#AlignedBases\t#MatchedBases\tMatchRate\t#MismatchedBases\tMismatchRate\t#DeletionsCigarD\tDeletionRate\tHomopolymerContextDel\t#InsertionsCigarI\tInsertionRate\tHomopolymerContextIns\t#SoftClipped\tSoftClipRate\t#HardClipped\tHardClipRate\tErrorRate" << "\t";
    rcfile << "MedianReadLength\tDefaultLibraryLayout\tMedianInsertSize\tMedianCoverage\tSDCoverage\tCoveredBp\tFractionCovered\tBpCov1ToCovNRatio\tBpCov1ToCov2Ratio\tMedianMAPQ";
    if (c.hasRegionFile) rcfile << "\t#TotalBedBp\t#AlignedBasesInBed\tFractionInBed\tEnrichmentOverBed";
    if (c.isMitagged) rcfile << "\t#MItagged\tFractionMItagged\t#UMIs";
    if (c.isHaplotagged) rcfile << "\t#HaploTagged\tFractionHaploTagged\t#PhasedBlocks\tN50PhasedBlockLength"; 
    rcfile << std::endl;
    for(typename TRGMap::const_iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
      // Read counts
      uint64_t totalReadCount = _totalReadCount(itRg);
      uint64_t mappedCount = itRg->second.rc.mapped1 + itRg->second.rc.mapped2;
      rcfile << "ME\t" << c.sampleName << "\t" << itRg->first << "\t" << itRg->second.rc.qcfail << "\t" << (double) itRg->second.rc.qcfail / (double) totalReadCount << "\t" << itRg->second.rc.dup << "\t" << (double) itRg->second.rc.dup / (double) totalReadCount << "\t" << itRg->second.rc.unmap << "\t" << (double) itRg->second.rc.unmap / (double) totalReadCount << "\t" << mappedCount << "\t" << (double) mappedCount / (double) totalReadCount << "\t" << itRg->second.rc.mapped1 << "\t" << itRg->second.rc.mapped2 << "\t" << (double) itRg->second.rc.mapped2 / (double) itRg->second.rc.mapped1 << "\t" << itRg->second.rc.forward << "\t" << (double) itRg->second.rc.forward / (double) mappedCount << "\t" << itRg->second.rc.reverse << "\t" << (double) itRg->second.rc.reverse / (double) mappedCount << "\t" << itRg->second.rc.secondary << "\t" << (double) itRg->second.rc.secondary / (double) mappedCount << "\t" << itRg->second.rc.supplementary << "\t" << (double) itRg->second.rc.supplementary / (double) mappedCount << "\t" << itRg->second.rc.spliced << "\t" << (double) itRg->second.rc.spliced / (double) mappedCount << "\t";

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
      rcfile << paired << "\t" << mapped << "\t" << mappedpairedfrac << "\t" << mappedSameChr << "\t" << mappedpairedchrfrac << "\t" << mappedProper << "\t" << mappedproperfrac << "\t";

      // Homopolymer Context of InDels
      double insFrac = _homopolymerIndel(itRg->second.bc.insHomACGTN);
      double delFrac = _homopolymerIndel(itRg->second.bc.delHomACGTN);
      
      // Error rates
      uint64_t alignedbases = itRg->second.bc.matchCount + itRg->second.bc.mismatchCount + itRg->second.bc.delCount + itRg->second.bc.insCount;
      double errRate = (double) (itRg->second.bc.mismatchCount + itRg->second.bc.delCount + itRg->second.bc.insCount) / (double) alignedbases;
      rcfile << rf.referencebp << "\t" << rf.ncount << "\t" << alignedbases << "\t" << itRg->second.bc.matchCount << "\t" << (double) itRg->second.bc.matchCount / (double) alignedbases << "\t" << itRg->second.bc.mismatchCount << "\t" << (double) itRg->second.bc.mismatchCount / (double) alignedbases << "\t" << itRg->second.bc.delCount << "\t" << (double) itRg->second.bc.delCount / (double) alignedbases << "\t" << delFrac << "\t" << itRg->second.bc.insCount << "\t" << (double) itRg->second.bc.insCount / (double) alignedbases << "\t" << insFrac << "\t" << itRg->second.bc.softClipCount << "\t" << (double) itRg->second.bc.softClipCount / (double) mappedCount << "\t" << itRg->second.bc.hardClipCount << "\t" << (double) itRg->second.bc.hardClipCount / (double) mappedCount << "\t" << errRate  << "\t";

      // Median coverage, read length, etc.
      int32_t deflayout = _defLayout(itRg);
      int32_t medISize = _medISize(itRg, deflayout);

      // Standardized SD of genomic coverage
      double ssdcov = 1000 * sdFromHistogram(itRg->second.bc.bpWithCoverage) / std::sqrt((double) mappedCount);
      double fraccovbp = (double) itRg->second.bc.nd / (double) (rf.referencebp - rf.ncount);
      double pbc1 = (double) itRg->second.bc.n1 / (double) itRg->second.bc.nd;
      double pbc2 = (double) itRg->second.bc.n1 / (double) itRg->second.bc.n2;

      if (itRg->second.pc.paired) rcfile << medianFromHistogram(itRg->second.rc.lRc[0]) << ":" << medianFromHistogram(itRg->second.rc.lRc[1]);
      else rcfile << medianFromHistogram(itRg->second.rc.lRc[0]);
      // or mean coverage: meanFromHistogram()
      rcfile << "\t" << deflayout << "\t" << medISize << "\t";
      if (c.meanCoverage) rcfile << meanFromHistogram(itRg->second.bc.bpWithCoverage) << "\t";
      else rcfile << medianFromHistogram(itRg->second.bc.bpWithCoverage) << "\t";
      rcfile << ssdcov << "\t" << itRg->second.bc.nd << "\t" << fraccovbp << "\t" << pbc1 << "\t" << pbc2 << "\t" << medianFromHistogram(itRg->second.qc.qcount);
      
      // Bed metrics
      if (c.hasRegionFile) {
	uint64_t nonN = rf.referencebp - rf.ncount;
	typename BedCounts::TOnTargetMap::const_iterator itOT = be.onTarget.find(itRg->first);
	uint64_t alignedBedBases = itOT->second[0];
	double fractioninbed = (double) alignedBedBases / (double) alignedbases;
	double enrichment = fractioninbed / ((double) rf.totalBedSize / (double) nonN);
	rcfile << "\t" << rf.totalBedSize << "\t" << alignedBedBases << "\t" << fractioninbed << "\t" << enrichment;
      }
      if (c.isMitagged) {
	rcfile << "\t" << itRg->second.rc.mitagged << "\t" << (double) itRg->second.rc.mitagged / (double) totalReadCount << "\t" << itRg->second.rc.umi.count();
      }
      if (c.isHaplotagged) {
	int32_t n50ps = n50PhasedBlockLength(itRg->second.rc.brange);
	rcfile << "\t" << itRg->second.rc.haplotagged << "\t" << (double) itRg->second.rc.haplotagged / (double) totalReadCount << "\t" << phasedBlocks(itRg->second.rc.brange) << "\t" << n50ps;
      }
      rcfile << std::endl;
    }

    // Output read length histogram
    rcfile << "# Read length distribution (RL)." << std::endl;
    rcfile << "# Use `zgrep ^RL <outfile> | cut -f 2-` to extract this part." << std::endl;
    rcfile << "RL\tSample\tReadlength\tCount\tFraction\tLibrary\tRead" << std::endl;
    for(typename TRGMap::const_iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
      uint32_t lastValidRL = _lastNonZeroIdx(itRg->second.rc.lRc[0]);
      double total = 0;
      for(uint32_t i = 0; i <= lastValidRL; ++i) total += itRg->second.rc.lRc[0][i];
      for(uint32_t i = 0; i <= lastValidRL; ++i) {
	double frac = 0;
	if (total > 0) frac = (double) itRg->second.rc.lRc[0][i] / total;
	rcfile << "RL\t" << c.sampleName << "\t" << i << "\t" << itRg->second.rc.lRc[0][i] << "\t" << frac << "\t" << itRg->first << "\tRead1" << std::endl;
      }
      if (itRg->second.pc.paired) {
	total = 0;
	lastValidRL = _lastNonZeroIdx(itRg->second.rc.lRc[1]);
	for(uint32_t i = 0; i <= lastValidRL; ++i) total += itRg->second.rc.lRc[1][i];
	for(uint32_t i = 0; i <= lastValidRL; ++i) {
	  double frac = 0;
	  if (total > 0) frac = (double) itRg->second.rc.lRc[1][i] / total;
	  rcfile << "RL\t" << c.sampleName << "\t" << i << "\t" << itRg->second.rc.lRc[1][i] << "\t" << frac << "\t" << itRg->first << "\tRead2" << std::endl;
	}
      }
    }

    // Output mean base quality
    rcfile << "# Mean base quality (BQ)." << std::endl;
    rcfile << "# Use `zgrep ^BQ <outfile> | cut -f 2-` to extract this part." << std::endl;
    rcfile << "BQ\tSample\tPosition\tBaseQual\tLibrary\tRead" << std::endl;
    for(typename TRGMap::const_iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
      uint32_t lastValidBQIdx = _lastNonZeroIdxACGTN(itRg->second.rc, 0);
      for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	uint64_t bcount = itRg->second.rc.aCount[0][i] + itRg->second.rc.cCount[0][i] + itRg->second.rc.gCount[0][i] + itRg->second.rc.tCount[0][i] + itRg->second.rc.nCount[0][i];
	if (bcount > 0) rcfile << "BQ\t" << c.sampleName << "\t" << i << "\t" << (double) (itRg->second.rc.bqCount[0][i]) / (double) (bcount) << "\t" << itRg->first << "\tRead1" << std::endl;
	else rcfile << "BQ\t" << c.sampleName << "\t" << i << "\tNA\t" << itRg->first << "\tRead1" << std::endl;
      }
      if (itRg->second.pc.paired) {
	lastValidBQIdx = _lastNonZeroIdxACGTN(itRg->second.rc, 1);
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  uint64_t bcount = itRg->second.rc.aCount[1][i] + itRg->second.rc.cCount[1][i] + itRg->second.rc.gCount[1][i] + itRg->second.rc.tCount[1][i] + itRg->second.rc.nCount[1][i];
	  if (bcount > 0) rcfile << "BQ\t" << c.sampleName << "\t" << i << "\t" << (double) (itRg->second.rc.bqCount[1][i]) / (double) (bcount) << "\t" << itRg->first << "\tRead2" << std::endl;
	  else rcfile << "BQ\t" << c.sampleName << "\t" << i << "\tNA\t" << itRg->first << "\tRead2" << std::endl;
	}
      }
    }

    // Output per base ACGTN content
    rcfile << "# Base content (BC)." << std::endl;
    rcfile << "# Use `zgrep ^BC <outfile> | cut -f 2-` to extract this part." << std::endl;
    rcfile << "BC\tSample\tPosition\tBase\tCount\tFraction\tLibrary\tRead" << std::endl;
    for(typename TRGMap::const_iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
      uint32_t lastValidBQIdx = _lastNonZeroIdxACGTN(itRg->second.rc, 0);
      for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	uint64_t bcount = itRg->second.rc.aCount[0][i] + itRg->second.rc.cCount[0][i] + itRg->second.rc.gCount[0][i] + itRg->second.rc.tCount[0][i] + itRg->second.rc.nCount[0][i];
	if (bcount > 0) {
	  rcfile << "BC\t" << c.sampleName << "\t" << i << "\tA\t" << itRg->second.rc.aCount[0][i] << "\t" << (double) itRg->second.rc.aCount[0][i] / (double) bcount << "\t" << itRg->first << "\tRead1" << std::endl;
	  rcfile << "BC\t" << c.sampleName << "\t" << i << "\tC\t" << itRg->second.rc.cCount[0][i] << "\t" << (double) itRg->second.rc.cCount[0][i] / (double) bcount << "\t" << itRg->first << "\tRead1" << std::endl;
	  rcfile << "BC\t" << c.sampleName << "\t" << i << "\tG\t" << itRg->second.rc.gCount[0][i] << "\t" << (double) itRg->second.rc.gCount[0][i] / (double) bcount << "\t" << itRg->first << "\tRead1" << std::endl;
	  rcfile << "BC\t" << c.sampleName << "\t" << i << "\tT\t" << itRg->second.rc.tCount[0][i] << "\t" << (double) itRg->second.rc.tCount[0][i] / (double) bcount << "\t" << itRg->first << "\tRead1" << std::endl;
	  rcfile << "BC\t" << c.sampleName << "\t" << i << "\tN\t" << itRg->second.rc.nCount[0][i] << "\t" << (double) itRg->second.rc.nCount[0][i] / (double) bcount << "\t" << itRg->first << "\tRead1" << std::endl;
	} else {
	  rcfile << "BC\t" << c.sampleName << "\t" << i << "\tA\t" << itRg->second.rc.aCount[0][i] << "\tNA\t" << itRg->first << "\tRead1" << std::endl;
	  rcfile << "BC\t" << c.sampleName << "\t" << i << "\tC\t" << itRg->second.rc.cCount[0][i] << "\tNA\t" << itRg->first << "\tRead1" << std::endl;
	  rcfile << "BC\t" << c.sampleName << "\t" << i << "\tG\t" << itRg->second.rc.gCount[0][i] << "\tNA\t" << itRg->first << "\tRead1" << std::endl;
	  rcfile << "BC\t" << c.sampleName << "\t" << i << "\tT\t" << itRg->second.rc.tCount[0][i] << "\tNA\t" << itRg->first << "\tRead1" << std::endl;
	  rcfile << "BC\t" << c.sampleName << "\t" << i << "\tN\t" << itRg->second.rc.nCount[0][i] << "\tNA\t" << itRg->first << "\tRead1" << std::endl;
	}
      }
      if (itRg->second.pc.paired) {
	lastValidBQIdx = _lastNonZeroIdxACGTN(itRg->second.rc, 1);
	for(uint32_t i = 0; i <= lastValidBQIdx; ++i) {
	  uint64_t bcount = itRg->second.rc.aCount[1][i] + itRg->second.rc.cCount[1][i] + itRg->second.rc.gCount[1][i] + itRg->second.rc.tCount[1][i] + itRg->second.rc.nCount[1][i];
	  if (bcount > 0) {
	    rcfile << "BC\t" << c.sampleName << "\t" << i << "\tA\t" << itRg->second.rc.aCount[1][i] << "\t" << (double) itRg->second.rc.aCount[1][i] / (double) bcount << "\t" << itRg->first << "\tRead2" << std::endl;
	    rcfile << "BC\t" << c.sampleName << "\t" << i << "\tC\t" << itRg->second.rc.cCount[1][i] << "\t" << (double) itRg->second.rc.cCount[1][i] / (double) bcount << "\t" << itRg->first << "\tRead2" << std::endl;
	    rcfile << "BC\t" << c.sampleName << "\t" << i << "\tG\t" << itRg->second.rc.gCount[1][i] << "\t" << (double) itRg->second.rc.gCount[1][i] / (double) bcount << "\t" << itRg->first << "\tRead2" << std::endl;
	    rcfile << "BC\t" << c.sampleName << "\t" << i << "\tT\t" << itRg->second.rc.tCount[1][i] << "\t" << (double) itRg->second.rc.tCount[1][i] / (double) bcount << "\t" << itRg->first << "\tRead2" << std::endl;
	    rcfile << "BC\t" << c.sampleName << "\t" << i << "\tN\t" << itRg->second.rc.nCount[1][i] << "\t" << (double) itRg->second.rc.nCount[1][i] / (double) bcount << "\t" << itRg->first << "\tRead2" << std::endl;
	  } else {
	    rcfile << "BC\t" << c.sampleName << "\t" << i << "\tA\t" << itRg->second.rc.aCount[1][i] << "\tNA\t" << itRg->first << "\tRead2" << std::endl;
	    rcfile << "BC\t" << c.sampleName << "\t" << i << "\tC\t" << itRg->second.rc.cCount[1][i] << "\tNA\t" << itRg->first << "\tRead2" << std::endl;
	    rcfile << "BC\t" << c.sampleName << "\t" << i << "\tG\t" << itRg->second.rc.gCount[1][i] << "\tNA\t" << itRg->first << "\tRead2" << std::endl;
	    rcfile << "BC\t" << c.sampleName << "\t" << i << "\tT\t" << itRg->second.rc.tCount[1][i] << "\tNA\t" << itRg->first << "\tRead2" << std::endl;
	    rcfile << "BC\t" << c.sampleName << "\t" << i << "\tN\t" << itRg->second.rc.nCount[1][i] << "\tNA\t" << itRg->first << "\tRead2" << std::endl;
	  }
	}
      }
    }


    // Output mapping quality histogram
    rcfile << "# Mapping quality histogram (MQ)." << std::endl;
    rcfile << "# Use `zgrep ^MQ <outfile> | cut -f 2-` to extract this part." << std::endl;
    rcfile << "MQ\tSample\tMappingQuality\tCount\tFraction\tLibrary" << std::endl;
    for(typename TRGMap::const_iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
      uint32_t lastValidMQ = 0;
      for(uint32_t i = 0; i < itRg->second.qc.qcount.size(); ++i)
	if (itRg->second.qc.qcount[i] > 0) lastValidMQ = i;
      double total = 0;
      for(uint32_t i = 0; i <= lastValidMQ; ++i) total += itRg->second.qc.qcount[i];
      for(uint32_t i = 0; i <= lastValidMQ; ++i) {
	double frac = 0;
	if (total > 0) frac = (double) itRg->second.qc.qcount[i] / total;
	rcfile << "MQ\t" << c.sampleName << "\t" << i << "\t" << itRg->second.qc.qcount[i] << "\t" << frac << "\t" << itRg->first << std::endl;
      }
    }    

    // Output coverage histograms
    rcfile << "# Coverage histogram (CO)." << std::endl;
    rcfile << "# Use `zgrep ^CO <outfile> | cut -f 2-` to extract this part." << std::endl;
    rcfile << "CO\tSample\tCoverage\tCount\tQuantile\tLibrary" << std::endl;
    for(typename TRGMap::const_iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
      uint32_t lastValidCO = 0;
      for(uint32_t i = 0; i < itRg->second.bc.bpWithCoverage.size(); ++i)
	if (itRg->second.bc.bpWithCoverage[i] > 0) lastValidCO = i;
      // Ignore last bucket that collects all higher coverage bases
      uint64_t totalCO = 0;
      for(uint32_t i = 0; i < itRg->second.bc.bpWithCoverage.size() - 1; ++i) totalCO += itRg->second.bc.bpWithCoverage[i];
      uint64_t cumsum = 0;
      for(uint32_t i = 0; i <= lastValidCO; ++i) {
	double quant = 0;
	if (totalCO > 0) quant = (double) cumsum / (double) totalCO;
	rcfile << "CO\t" << c.sampleName << "\t" << i << "\t" << itRg->second.bc.bpWithCoverage[i] << "\t" << quant << "\t" << itRg->first << std::endl;
	cumsum += itRg->second.bc.bpWithCoverage[i];
      }
    }

    // Output mapping statistics by chromosome
    rcfile << "# Chromosome mapping statistics (CM)." << std::endl;
    rcfile << "# Use `zgrep ^CM <outfile> | cut -f 2-` to extract this part." << std::endl;
    rcfile << "CM\tSample\tLibrary\tChrom\tSize\tMapped\tMappedFraction\tObsExpRatio" << std::endl;
    for(typename TRGMap::const_iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
      uint64_t totalMappedChr = 0;
      for(uint32_t i = 0; i < itRg->second.rc.mappedchr.size(); ++i) totalMappedChr += itRg->second.rc.mappedchr[i];
      for(uint32_t i = 0; i < itRg->second.rc.mappedchr.size(); ++i) {
	double frac = 0;
	if (totalMappedChr > 0) frac = (double) itRg->second.rc.mappedchr[i] / (double) totalMappedChr;
	double expect = (double) (hdr->target_len[i] - rf.chrGC[i].ncount) / (double) (rf.referencebp - rf.ncount);
	double obsexprat = frac / expect;
	rcfile << "CM\t" << c.sampleName << "\t" << itRg->first << "\t" << hdr->target_name[i] << "\t" << hdr->target_len[i] << "\t" << itRg->second.rc.mappedchr[i] << "\t" << frac << "\t" << obsexprat << std::endl;
      }
    }
    
    // Output insert size histograms
    rcfile << "# Insert size histogram (IS)." << std::endl;
    rcfile << "# Use `zgrep ^IS <outfile> | cut -f 2-` to extract this part." << std::endl;
    rcfile << "IS\tSample\tInsertSize\tCount\tLayout\tQuantile\tLibrary" << std::endl;
    for(typename TRGMap::const_iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
      uint32_t lastValidISIdx = _lastNonZeroIdxISize(itRg->second.pc);
      // Ignore last bucket that collects all other pairs
      uint64_t totalFR = 0;
      for(uint32_t i = 0; i < itRg->second.pc.fPlus.size() - 1; ++i) totalFR += itRg->second.pc.fPlus[i] + itRg->second.pc.fMinus[i] + itRg->second.pc.rPlus[i] + itRg->second.pc.rMinus[i];
      uint64_t cumsum = 0;
      for(uint32_t i = 0; i <= lastValidISIdx; ++i) {
	double quant = 0;
	if (totalFR > 0) quant = (double) cumsum / (double) totalFR;
	rcfile << "IS\t" << c.sampleName << "\t" << i << "\t" << itRg->second.pc.fPlus[i] << "\tF+\t" << quant << "\t" << itRg->first << std::endl;
	rcfile << "IS\t" << c.sampleName << "\t" << i << "\t" << itRg->second.pc.fMinus[i] << "\tF-\t" << quant << "\t" << itRg->first << std::endl;
	rcfile << "IS\t" << c.sampleName << "\t" << i << "\t" << itRg->second.pc.rPlus[i] << "\tR+\t" << quant << "\t" << itRg->first << std::endl;
	rcfile << "IS\t" << c.sampleName << "\t" << i << "\t" << itRg->second.pc.rMinus[i] << "\tR-\t" << quant << "\t" << itRg->first << std::endl;		
	cumsum += itRg->second.pc.fPlus[i] + itRg->second.pc.fMinus[i] + itRg->second.pc.rPlus[i] + itRg->second.pc.rMinus[i];
      }
    }

    // Homopolymer InDel context
    rcfile << "# InDel context (IC)." << std::endl;
    rcfile << "# Use `zgrep ^IC <outfile> | cut -f 2-` to extract this part." << std::endl;
    rcfile << "IC\tSample\tLibrary\tInDel\tHomopolymer\tCount\tFraction" << std::endl;
    for(typename TRGMap::const_iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
      double total = 0;
      for(uint32_t i = 0; i < itRg->second.bc.delHomACGTN.size(); ++i) total += itRg->second.bc.delHomACGTN[i];
      for(uint32_t i = 0; i < itRg->second.bc.delHomACGTN.size(); ++i) {
	double frac = 0;
	if (total > 0) frac = (double) itRg->second.bc.delHomACGTN[i] / total;
	rcfile << "IC\t" << c.sampleName << "\t" << itRg->first << "\tDEL\t";
	if (i == 0) rcfile << 'A';
	else if (i == 1) rcfile << 'C';
	else if (i == 2) rcfile << 'G';
	else if (i == 3) rcfile << 'T';
	else if (i == 4) rcfile << 'N';
	else rcfile << "None";
	rcfile << "\t" << itRg->second.bc.delHomACGTN[i] << "\t" << frac << std::endl;
      }
    }
    for(typename TRGMap::const_iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
      double total = 0;
      for(uint32_t i = 0; i < itRg->second.bc.insHomACGTN.size(); ++i) total += itRg->second.bc.insHomACGTN[i];
      for(uint32_t i = 0; i < itRg->second.bc.insHomACGTN.size(); ++i) {
	double frac = 0;
	if (total > 0) frac = (double) itRg->second.bc.insHomACGTN[i] / total;
	rcfile << "IC\t" << c.sampleName << "\t" << itRg->first << "\tINS\t";
	if (i == 0) rcfile << 'A';
	else if (i == 1) rcfile << 'C';
	else if (i == 2) rcfile << 'G';
	else if (i == 3) rcfile << 'T';
	else if (i == 4) rcfile << 'N';
	else rcfile << "None";
	rcfile << "\t" << itRg->second.bc.insHomACGTN[i] << "\t" << frac << std::endl;
      }
    }

    // InDel size
    rcfile << "# InDel size (IZ)." << std::endl;
    rcfile << "# Use `zgrep ^IZ <outfile> | cut -f 2-` to extract this part." << std::endl;
    rcfile << "IZ\tSample\tLibrary\tInDel\tSize\tCount" << std::endl;
    for(typename TRGMap::const_iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
      for(uint32_t i = 1; i < itRg->second.bc.delSize.size(); ++i) 
	rcfile << "IZ\t" << c.sampleName << "\t" << itRg->first << "\tDEL\t" << i << "\t" << itRg->second.bc.delSize[i] << std::endl;
      for(uint32_t i = 1; i < itRg->second.bc.insSize.size(); ++i) 
	rcfile << "IZ\t" << c.sampleName << "\t" << itRg->first << "\tINS\t" << i << "\t" << itRg->second.bc.insSize[i] << std::endl;
    }

    // Reference GC content
    rcfile << "# Chromosome GC-content (CG)." << std::endl;
    rcfile << "# Use `zgrep ^CG <outfile> | cut -f 2-` to extract this part." << std::endl;
    rcfile << "CG\tChromosome\tSize\tnumN\tnumGC\tGCfraction" << std::endl;
    for(uint32_t i = 0; i < rf.chrGC.size(); ++i) {
      // Only chromosomes with mapped data
      if (rf.chrGC[i].ncount + rf.chrGC[i].gccount > 0) {
	double total = hdr->target_len[i] - rf.chrGC[i].ncount;
	if (total > 0) {
	  double frac = (double) rf.chrGC[i].gccount / total;
	  rcfile << "CG\t" << hdr->target_name[i] << "\t" << hdr->target_len[i] << "\t" << rf.chrGC[i].ncount << "\t" << rf.chrGC[i].gccount << "\t" << frac << std::endl;
	}
      }
    }

    // GC-content
    rcfile << "# GC-content (GC)." << std::endl;
    rcfile << "# Use `zgrep ^GC <outfile> | cut -f 2-` to extract this part." << std::endl;
    rcfile << "GC\tSample\tLibrary\tGCcontent\tfractionOfReads" << std::endl;
    {
      double total = 0;
      for(uint32_t i = 0; i < 102; ++i) total += rf.refGcContent[i];
      for(uint32_t i = 0; i < 102; ++i) {
	double frac = 0;
	if (total > 0) frac = (double) rf.refGcContent[i] / total;
	rcfile << "GC\tReference\tReference\t" << (double) i / (double) 101 << "\t" << frac << std::endl;
      }
    }
    if (c.hasRegionFile) {
      double total = 0;
      for(uint32_t i = 0; i < 102; ++i) total += be.bedGcContent[i];
      for(uint32_t i = 0; i < 102; ++i) {
	double frac = 0;
	if (total > 0) frac = (double) be.bedGcContent[i] / total;
	rcfile << "GC\tTarget\tTarget\t" << (double) i / (double) 101 << "\t" << frac << std::endl;
      }
    }
    for(typename TRGMap::const_iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
      double total = 0;
      for(uint32_t i = 0; i < 102; ++i) total += itRg->second.rc.gcContent[i];
      for(uint32_t i = 0; i < 102; ++i) {
	double frac = 0;
	if (total > 0) frac = (double) itRg->second.rc.gcContent[i] / total;
	rcfile << "GC\t" << c.sampleName << "\t" << itRg->first << "\t" << (double) i / (double) 101 << "\t" << frac << std::endl;
      }
    }

    if (c.hasRegionFile) {
      // Output avg. bed coverage
      rcfile << "# Avg. target coverage (TC)." << std::endl;
      rcfile << "# Use `zgrep ^TC <outfile> | cut -f 2-` to extract this part." << std::endl;
      rcfile << "TC\tSample\tLibrary\tChr\tStart\tEnd\tAvgCov" << std::endl;
      for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
	for(typename BedCounts::TRgBpMap::const_iterator itChr = be.gCov[refIndex].begin(); itChr != be.gCov[refIndex].end(); ++itChr) {
	  for(uint32_t i = 0; i < rf.gRegions[refIndex].size(); ++i) {
	    rcfile << "TC\t" << c.sampleName << "\t" << itChr->first << "\t" << hdr->target_name[refIndex] << "\t" << rf.gRegions[refIndex][i].start << "\t" << rf.gRegions[refIndex][i].end << "\t" << itChr->second[i] << std::endl;
	  }
	}
      }
      
      // Output on-target rate
      rcfile << "# On target rate (OT)." << std::endl;
      rcfile << "# Use `zgrep ^OT <outfile> | cut -f 2-` to extract this part." << std::endl;
      rcfile << "OT\tSample\tLibrary\tExtension\tOnTarget" << std::endl;
      for(typename TRGMap::const_iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
	uint64_t alignedbases = itRg->second.bc.matchCount + itRg->second.bc.mismatchCount;
	typename BedCounts::TOnTargetMap::const_iterator itOT = be.onTarget.find(itRg->first);
	for(uint32_t k = 0; k < itOT->second.size(); ++k) {
	  rcfile << "OT\t" << c.sampleName << "\t" << itRg->first << "\t" << k * be.stepsize << "\t" << (double) itOT->second[k] / (double) alignedbases << std::endl;
	}
      }
    }

    if (c.isHaplotagged) {
      // Output phased block length histogram
      rcfile << "# Phased block length (PS)." << std::endl;
      rcfile << "# Use `zgrep ^PS <outfile> | cut -f 2-` to extract this part." << std::endl;
      rcfile << "PS\tSample\tChr\tStart\tEnd\tPSid\tSize\tLibrary" << std::endl;
      for(typename TRGMap::const_iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
	for(int32_t refIndex = 0; refIndex < (int32_t) itRg->second.rc.brange.size(); ++refIndex) {
	  for(typename ReadCounts::TBlockRange::const_iterator itBR = itRg->second.rc.brange[refIndex].begin(); itBR != itRg->second.rc.brange[refIndex].end(); ++itBR) {
	    if (itBR->second.first < itBR->second.second) {
	      rcfile << "PS\t" << c.sampleName << "\t" <<  hdr->target_name[refIndex] << "\t" << itBR->second.first << "\t" << itBR->second.second << "\t" << itBR->first << "\t" << (itBR->second.second - itBR->second.first) << "\t" << itRg->first << std::endl;
	    }
	  }
	}
      }
    }
  }  
 
}

#endif
