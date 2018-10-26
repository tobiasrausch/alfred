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

#ifndef SWGOTOH_H
#define SWGOTOH_H

#include <boost/dynamic_bitset.hpp>

#include <iostream>
#include "align.h"

namespace bamstats
{

  template<typename TAlign1, typename TAlign2, typename TAlign, typename TAlignConfig, typename TScoreObject>
  inline int
  swGotoh(TAlign1 const& a1, TAlign2 const& a2, TAlign& align, TAlignConfig const& ac, TScoreObject const& sc)
  {
    typedef typename TScoreObject::TValue TScoreValue;

    // DP variables
    std::size_t m = _size(a1, 1);
    std::size_t n = _size(a2, 1);
    std::vector<TScoreValue> s(n+1, 0);
    std::vector<TScoreValue> v(n+1, 0);
    TScoreValue newhoz = 0;
    TScoreValue prevsub = 0;
    TScoreValue maxVal = 0;
    std::size_t rowMax = 0;
    std::size_t colMax = 0;
    
    // Trace Matrix
    std::size_t mf = n+1;
    typedef boost::dynamic_bitset<> TBitSet;
    TBitSet bit1( (m+1) * (n+1), false);
    TBitSet bit2( (m+1) * (n+1), false);
    TBitSet bit3( (m+1) * (n+1), false);
    TBitSet bit4( (m+1) * (n+1), false);

    // Create profile
    typedef boost::multi_array<float, 2> TProfile;
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
	  v[0] = -sc.inf;
	  newhoz = -sc.inf;
	  prevsub = 0;
	  bit1[row * mf + col] = true;
          bit2[row * mf + col] = true;
	  bit3[row * mf + col] = true;
	  bit4[row * mf + col] = true;
	} else {
	  // Recursion
	  TScoreValue prevhoz = newhoz;
	  TScoreValue prevver = v[col];
	  TScoreValue prevprevsub = prevsub;
	  prevsub = s[col];
	  newhoz = std::max(s[col-1] + _horizontalGap(ac, row, m, sc.go + sc.ge), prevhoz + _horizontalGap(ac, row, m, sc.ge));
	  v[col] = std::max(prevsub + _verticalGap(ac, col, n, sc.go + sc.ge), prevver + _verticalGap(ac, col, n, sc.ge));
	  s[col] = std::max(std::max(std::max(prevprevsub + _score(a1, a2, p1, p2, row-1, col-1, sc), newhoz), v[col]), 0);

	  // Trace
	  if (s[col] == 0) {
	    bit3[row * mf + col] = true;
	    bit4[row * mf + col] = true;
	  } else {
	    if (s[col] == newhoz) bit3[row * mf + col] = true;
	    else if (s[col] == v[col]) bit4[row * mf + col] = true;
	    if (newhoz != prevhoz + _horizontalGap(ac, row, m, sc.ge)) bit1[row * mf + col] = true;
	    if (v[col] != prevver + _verticalGap(ac, col, n, sc.ge)) bit2[row * mf + col] = true;
	  }

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
    char lastMatrix = 's';
    typedef std::vector<char> TTrace;
    TTrace btr;
    while ((row>0) || (col>0)) {
      if ((bit3[row * mf + col]) && (bit4[row * mf + col])) {
	break;
      } else {
	if (lastMatrix == 's') {
	  if (bit3[row * mf + col]) lastMatrix = 'h';
	  else if (bit4[row * mf + col]) lastMatrix = 'v';
	  else {
	    --row;
	    --col;
	    btr.push_back('s');
	  }
	} else if (lastMatrix == 'h') {
	  if (bit1[row * mf + col]) lastMatrix = 's';
	  --col;
	  btr.push_back('h');
	} else if (lastMatrix == 'v') {
	  if (bit2[row * mf + col]) lastMatrix = 's';
	  --row;
	  btr.push_back('v');
	}
      }
    }
      

    // Create alignment
    _createLocalAlignment(btr, a1, a2, align, row, col);

    // Score
    return maxVal;
  }

  template<typename TAlign1, typename TAlign2, typename TAlign, typename TAlignConfig>
  inline int
  swGotoh(TAlign1 const& a1, TAlign2 const& a2, TAlign& align, TAlignConfig const& ac) 
  {
    DnaScore<int> dnasc;    
    return swGotoh(a1, a2, align, ac, dnasc);
  }

  template<typename TAlign1, typename TAlign2, typename TAlign>
  inline int
  swGotoh(TAlign1 const& a1, TAlign2 const& a2, TAlign& align)
  {
    AlignConfig<false, false> ac;
    return swGotoh(a1, a2, align, ac);
  }

}

#endif
