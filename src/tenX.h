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
    for(int32_t refIndex = 0; refIndex < (int32_t) brange.size(); ++refIndex) {
      for(int32_t i = 0; i < (int32_t) brange[refIndex].size(); ++i) {
	// Ignore "unused" phased block identifiers
	if ((brange[refIndex][i].first != -1) && (brange[refIndex][i].second != -1)) ++count;
      }
    }
    return count;
  }
  
  template<typename TGenomicBlockRanges>
  inline int32_t
  n50PhasedBlockLength(TGenomicBlockRanges const& brange) {
    typedef std::vector<int32_t> TSizes;
    TSizes sz;
    int64_t totalSize = 0;
    for(int32_t refIndex = 0; refIndex < (int32_t) brange.size(); ++refIndex) {
      for(int32_t i = 0; i < (int32_t) brange[refIndex].size(); ++i) {
	if (brange[refIndex][i].first < brange[refIndex][i].second) {
	  sz.push_back(brange[refIndex][i].second - brange[refIndex][i].first);
	  totalSize += (brange[refIndex][i].second - brange[refIndex][i].first);
	}
	else {
	  if ((brange[refIndex][i].first != -1) && (brange[refIndex][i].second != -1)) {
	    std::cerr << "Warning: Phased block start after phased block end!" << std::endl;
	  }
	}
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
