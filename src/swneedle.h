/*
============================================================================
Alfred
============================================================================
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

#ifndef SWNEEDLE_H
#define SWNEEDLE_H

#define BOOST_DISABLE_ASSERTS
#include <boost/multi_array.hpp>
#include <iostream>
#include "align.h"

namespace bamstats
{

  template<typename TAlign1, typename TAlign2, typename TAlign, typename TAlignConfig, typename TScoreObject>
  inline int
  swNeedle(TAlign1 const& a1, TAlign2 const& a2, TAlign& align, TAlignConfig const& ac, TScoreObject const& sc)
  {
    typedef typename TScoreObject::TValue TScoreValue;
    TScoreValue maxVal = 0;
    std::size_t rowMax = 0;
    std::size_t colMax = 0;
    
    // DP Matrix
    std::size_t m = _size(a1, 1);
    std::size_t n = _size(a2, 1);
    std::vector<TScoreValue> s(n+1, 0);
    TScoreValue prevsub = 0;

    // Trace Matrix
    std::size_t mf = n+1;
    typedef boost::dynamic_bitset<> TBitSet;
    TBitSet bit3( (m+1) * (n+1), false);
    TBitSet bit4( (m+1) * (n+1), false);
    
    // Create profile
    typedef boost::multi_array<double, 2> TProfile;
    TProfile p1;
    TProfile p2;
    if ((_size(a1, 0) != 1) || (_size(a2, 0) != 1)) {
      _createProfile(a1, p1);
      _createProfile(a2, p2);
    }

    // DP
    for(std::size_t row = 0; row <= m; ++row) {
      for(std::size_t col = 0; col <= n; ++col) {
	// Initialization
	if ((row == 0) || (col == 0)) {
	  s[0] = 0;
	  prevsub = 0;
	  bit3[row * mf + col] = true;
	  bit4[row * mf + col] = true;
	} else {
	  // Recursion
	  TScoreValue prevprevsub = prevsub;
	  prevsub = s[col];
	  s[col] = std::max(std::max(std::max(prevprevsub + _score(a1, a2, p1, p2, row-1, col-1, sc), prevsub + _verticalGap(ac, col, n, sc.ge)), s[col-1] + _horizontalGap(ac, row, m, sc.ge)), 0);

	  // Trace
	  if (s[col] == 0) {
	    bit3[row * mf + col] = true;
	    bit4[row * mf + col] = true;
	  }
	  else if (s[col] ==  s[col-1] + _horizontalGap(ac, row, m, sc.ge)) bit3[row * mf + col] = true;
	  else if (s[col] == prevsub + _verticalGap(ac, col, n, sc.ge)) bit4[row * mf + col] = true;

	  // Highest Score
	  if (s[col] > maxVal) {
	    maxVal = s[col];
	    rowMax = row;
	    colMax = col;
	  }
	}
      }
    }
	
    // Trace-back using pointers
    std::size_t row = rowMax;
    std::size_t col = colMax;
    typedef std::vector<char> TTrace;
    TTrace trace;
    while ((row>0) || (col>0)) {
      if ((bit3[row * mf + col]) && (bit4[row * mf + col])) {
	break;
      } else if (bit3[row * mf + col]) {
	--col;
	trace.push_back('h');
      } else if (bit4[row * mf + col]) {
	--row;
	trace.push_back('v');
      } else {
	--row;
	--col;
	trace.push_back('s');
      }
    }
    
    // Create alignment
    _createLocalAlignment(trace, a1, a2, align, row, col);

    // Score
    return maxVal;
  }

  template<typename TAlign1, typename TAlign2, typename TAlign, typename TAlignConfig>
  inline int
  swNeedle(TAlign1 const& a1, TAlign2 const& a2, TAlign& align, TAlignConfig const& ac)
  {
    DnaScore<int> dnasc;
    return swNeedle(a1, a2, align, ac, dnasc);
  }

  template<typename TAlign1, typename TAlign2, typename TAlign>
  inline int
  swNeedle(TAlign1 const& a1, TAlign2 const& a2, TAlign& align)
  {
    AlignConfig<false, false> ac;
    return swNeedle(a1, a2, align, ac);
  }

}

#endif
