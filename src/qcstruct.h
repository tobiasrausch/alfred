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

#ifndef QCSTRUCT_H
#define QCSTRUCT_H

#include <limits>

#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/progress.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "util.h"

namespace bamstats
{

  struct ChrGC {
    uint32_t ncount;
    uint32_t gccount;
  };

  struct ReferenceFeatures {
    typedef std::vector<uint64_t> TGCContent;
    typedef std::vector<Interval> TChromosomeRegions;
    typedef std::vector<TChromosomeRegions> TGenomicRegions;
    
    uint64_t referencebp;
    uint64_t ncount;
    uint32_t totalBedSize;
    uint32_t nchr;
    std::vector<ChrGC> chrGC;
    TGenomicRegions gRegions;
    TGCContent refGcContent;
    

    explicit ReferenceFeatures(uint32_t const nc) : referencebp(0), ncount(0), totalBedSize(0), nchr(nc) {
      chrGC.resize(nc, ChrGC());
      gRegions.resize(nc, TChromosomeRegions());
      refGcContent.resize(102, 0);
    }
  };

  
  struct BaseCounts {
    typedef uint32_t TCountType;
    typedef std::vector<TCountType> TCoverageBp;
    
    typedef uint16_t TMaxCoverage;
    typedef std::vector<TMaxCoverage> TBpCoverage;

    uint32_t maxCoverage;
    uint32_t maxIndelSize;
    uint64_t n1;
    uint64_t n2;
    uint64_t nd;
    uint64_t matchCount;
    uint64_t mismatchCount;
    uint64_t delCount;
    uint64_t insCount;
    uint64_t softClipCount;
    uint64_t hardClipCount;
    std::vector<uint32_t> delHomACGTN;  // A:0, C:1, G:2, T:3, N:4, none:5
    std::vector<uint32_t> insHomACGTN;  // A:0, C:1, G:2, T:3, N:4, none:5
    std::vector<uint32_t> delSize;
    std::vector<uint32_t> insSize;
    TCoverageBp bpWithCoverage;
    TBpCoverage cov;

    BaseCounts() : maxCoverage(std::numeric_limits<TMaxCoverage>::max()), maxIndelSize(50), n1(0), n2(0), nd(0), matchCount(0), mismatchCount(0), delCount(0), insCount(0), softClipCount(0), hardClipCount(0) {
      delHomACGTN.resize(6, 0);
      insHomACGTN.resize(6, 0);
      bpWithCoverage.resize(maxCoverage + 1, 0);
      delSize.resize(maxIndelSize + 1, 0);
      insSize.resize(maxIndelSize + 1, 0);
      cov.clear();
    }
  };

  struct ReadCounts {
    typedef uint16_t TMaxReadLength;
    typedef uint32_t TCountType;
    typedef std::vector<TCountType> TLengthReadCount;
    typedef std::vector<uint64_t> TBaseQualitySum;
    typedef ReferenceFeatures::TGCContent TGCContent;
    typedef boost::dynamic_bitset<> TBitSet;
    typedef std::pair<int32_t, int32_t> TStartEndPair;
    typedef std::map<int32_t, TStartEndPair> TBlockRange;
    typedef std::vector<TBlockRange> TGenomicBlockRange;
    typedef std::vector<uint64_t> TMappedChr;
    
    int32_t maxReadLength;
    int32_t maxUMI;
    int64_t secondary;
    int64_t qcfail;
    int64_t dup;
    int64_t supplementary;
    int64_t unmap;
    int64_t forward;
    int64_t reverse;
    int64_t spliced;
    int64_t mapped1;
    int64_t mapped2;
    int64_t haplotagged;
    int64_t mitagged;
    TMappedChr mappedchr;
    TLengthReadCount lRc;
    TLengthReadCount nCount;
    TLengthReadCount aCount;
    TLengthReadCount cCount;
    TLengthReadCount gCount;
    TLengthReadCount tCount;
    TBaseQualitySum bqCount;
    TGCContent gcContent;
    TBitSet umi;
    TGenomicBlockRange brange;

    ReadCounts(uint32_t const n_targets) : maxReadLength(std::numeric_limits<TMaxReadLength>::max()), maxUMI(10000000), secondary(0), qcfail(0), dup(0), supplementary(0), unmap(0), forward(0), reverse(0), spliced(0), mapped1(0), mapped2(0), haplotagged(0), mitagged(0) {
      mappedchr.resize(n_targets, 0);
      lRc.resize(maxReadLength + 1, 0);
      aCount.resize(maxReadLength + 1, 0);
      cCount.resize(maxReadLength + 1, 0);
      gCount.resize(maxReadLength + 1, 0);
      tCount.resize(maxReadLength + 1, 0);
      nCount.resize(maxReadLength + 1, 0);
      bqCount.resize(maxReadLength + 1, 0);
      gcContent.resize(102, 0);
    }
  };

  struct PairCounts {
    typedef uint16_t TMaxInsertSize;
    typedef uint32_t TCountType;
    typedef std::vector<TCountType> TISizePairCount;
    int32_t maxInsertSize;
    int64_t paired;
    int64_t mapped;
    int64_t mappedSameChr;
    int64_t mappedProper;
    int64_t orient[4];
    int64_t totalISizeCount;
    TISizePairCount fPlus;
    TISizePairCount rPlus;
    TISizePairCount fMinus;
    TISizePairCount rMinus;
    
    
    PairCounts() : maxInsertSize(std::numeric_limits<TMaxInsertSize>::max()), paired(0), mapped(0), mappedSameChr(0), mappedProper(0), totalISizeCount(0) {
      orient[0] = 0;
      orient[1] = 0;
      orient[2] = 0;
      orient[3] = 0;
      fPlus.resize(maxInsertSize + 1, 0);
      rPlus.resize(maxInsertSize + 1, 0);
      fMinus.resize(maxInsertSize + 1, 0);
      rMinus.resize(maxInsertSize + 1, 0);
    }
  };

  
  struct QualCounts {
    typedef uint8_t TMaxQuality;
    typedef uint32_t TCountType;
    typedef std::vector<TCountType> TQualCount;
    int32_t maxQuality;
    TQualCount qcount;

    QualCounts() : maxQuality(std::numeric_limits<TMaxQuality>::max()) {
      qcount.resize(maxQuality + 1, 0);
    }
  };
    
  
  struct ReadGroupStats {
    BaseCounts bc;
    ReadCounts rc;
    PairCounts pc;
    QualCounts qc;
    
  ReadGroupStats(uint32_t const n_targets) : bc(BaseCounts()), rc(ReadCounts(n_targets)), pc(PairCounts()), qc(QualCounts()) {}
  };


  struct BedCounts {
    typedef double TAvgCov;
    typedef std::vector<TAvgCov> TBpCov;
    typedef boost::unordered_map<std::string, TBpCov> TRgBpMap;
    typedef std::vector<TRgBpMap> TGenomicBp;

    typedef std::vector<int64_t> TOnTargetBp;
    typedef boost::unordered_map<std::string, TOnTargetBp> TOnTargetMap;
    typedef std::vector<uint64_t> TGCContent;
    
    int32_t stepsize;
    int32_t onTSize;
    TGenomicBp gCov;
    TOnTargetMap onTarget;
    TGCContent bedGcContent;
    
    BedCounts(int32_t nchr, int32_t s, int32_t vs) : stepsize(s), onTSize(vs) {
      gCov.resize(nchr, TRgBpMap());
      bedGcContent.resize(102, 0);
    }
  };
  
}

#endif
