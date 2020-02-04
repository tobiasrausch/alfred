#ifndef TENX_H
#define TENX_H

#include <boost/unordered_map.hpp>
#include <boost/algorithm/string.hpp>
#include <htslib/sam.h>


namespace bamstats
{

  template<typename TGenomicBlockRanges>
  inline int32_t
  phasedBlocks(TGenomicBlockRanges const& brange) {
    int32_t count = 0;
    for(int32_t refIndex = 0; refIndex < (int32_t) brange.size(); ++refIndex) count += brange[refIndex].size();
    return count;
  }
  
  template<typename TGenomicBlockRanges>
  inline int32_t
  n50PhasedBlockLength(TGenomicBlockRanges const& brange) {
    typedef typename TGenomicBlockRanges::value_type TBlockRange;
    typedef std::vector<int32_t> TSizes;
    TSizes sz;
    int64_t totalSize = 0;
    for(int32_t refIndex = 0; refIndex < (int32_t) brange.size(); ++refIndex) {
      for(typename TBlockRange::const_iterator itBR = brange[refIndex].begin(); itBR != brange[refIndex].end(); ++itBR) {
	if (itBR->second.first < itBR->second.second) {
	  sz.push_back(itBR->second.second - itBR->second.first);
	  totalSize += (itBR->second.second - itBR->second.first);
	}
	else std::cerr << "Warning: Phased block start after phased block end!" << std::endl;
      }
    }
    std::sort(sz.begin(), sz.end(), std::greater<int32_t>());
    totalSize /= 2;
    int64_t cumSize = 0;
    for(int32_t i = 0; i < (int32_t) sz.size(); ++i) {
      cumSize += sz[i];
      if (cumSize > totalSize) return sz[i];
    }
    return 0;
  }

}

#endif
