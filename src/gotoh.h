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

#ifndef GOTOH_H
#define GOTOH_H

#include <boost/dynamic_bitset.hpp>

#include <iostream>
#include "align.h"

namespace bamstats
{
  
  template<typename TAlign1, typename TAlign2, typename TAlign, typename TAlignConfig, typename TScoreObject>
  inline int
    gotoh(TAlign1 const& a1, TAlign2 const& a2, TAlign& align, TAlignConfig const& ac, TScoreObject const& sc)
  {
    typedef typename TScoreObject::TValue TScoreValue;

    // DP variables
    std::size_t m = _size(a1, 1);
    std::size_t n = _size(a2, 1);
    std::vector<TScoreValue> s(n+1, 0);
    std::vector<TScoreValue> v(n+1, 0);
    TScoreValue newhoz = 0;
    TScoreValue prevsub = 0;
    
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
	if ((row == 0) && (col == 0)) {
	  s[0] = 0;
	  v[0] = -sc.inf;
	  newhoz = -sc.inf;
	  bit1[0 * mf + 0] = true;
	  bit2[0 * mf + 0] = true;
	  bit3[0 * mf + 0] = false;
	  bit4[0 * mf + 0] = false;
	  //tracept[0][0] = 'a';
	} else if (row == 0) {
	  v[col] = -sc.inf;
	  s[col] = _horizontalGap(ac, 0, m, sc.go + col * sc.ge);
	  newhoz = _horizontalGap(ac, 0, m, sc.go + col * sc.ge);
	  bit1[0 * mf + col] = false;
	  bit2[0 * mf + col] = false;
	  bit3[0 * mf + col] = true;
	  bit4[0 * mf + col] = false;
	  //tracept[0][col] = 'h';
	} else if (col == 0) {
	  newhoz = -sc.inf;
	  s[0] = _verticalGap(ac, 0, n, sc.go + row * sc.ge);
	  prevsub = s[0];
	  v[0] = _verticalGap(ac, 0, n, sc.go + row * sc.ge);
	  bit1[row * mf + 0] = false;
	  bit2[row * mf + 0] = false;
	  bit3[row * mf + 0] = false;
	  bit4[row * mf + 0] = true;
	  //tracept[row][0] = 'l';
	} else {
	  // Recursion
	  TScoreValue prevhoz = newhoz;
	  TScoreValue prevver = v[col];
	  TScoreValue prevprevsub = prevsub;
	  prevsub = s[col];
	  newhoz = std::max(s[col-1] + _horizontalGap(ac, row, m, sc.go + sc.ge), prevhoz + _horizontalGap(ac, row, m, sc.ge));
	  v[col] = std::max(prevsub + _verticalGap(ac, col, n, sc.go + sc.ge), prevver + _verticalGap(ac, col, n, sc.ge));
	  s[col] = std::max(std::max(prevprevsub + _score(a1, a2, p1, p2, row-1, col-1, sc), newhoz), v[col]);

	  // Trace
	  char subval = 's';
	  char hozval = 'h';
	  char verval = 'v'; 
	  if (s[col] == newhoz) subval = 'h';
	  else if (s[col] == v[col]) subval = 'v';
	  if (newhoz != prevhoz + _horizontalGap(ac, row, m, sc.ge)) hozval = 's';
	  if (v[col] != prevver + _verticalGap(ac, col, n, sc.ge)) verval = 's';
	  if (subval == 's') {
	    if (hozval == 's') {
	      if (verval == 's') {
		bit1[row * mf + col] = true;
		bit2[row * mf + col] = true;
		bit3[row * mf + col] = false;
		bit4[row * mf + col] = false;
		//tracept[row][col] = 'a';
	      } else {
		bit1[row * mf + col] = true;
		bit2[row * mf + col] = false;
		bit3[row * mf + col] = false;
		bit4[row * mf + col] = false;
		//tracept[row][col] = 'b';
	      }
	    } else {
	      if (verval == 's') {
		bit1[row * mf + col] = false;
		bit2[row * mf + col] = true;
		bit3[row * mf + col] = false;
		bit4[row * mf + col] = false;
		//tracept[row][col] = 'c';
	      } else {
		bit1[row * mf + col] = false;
		bit2[row * mf + col] = false;
		bit3[row * mf + col] = false;
		bit4[row * mf + col] = false;
		//tracept[row][col] = 'd';
	      }
	    }
	  } else if (subval == 'h') {
	    if (hozval == 's') {
	      if (verval == 's') {
		bit1[row * mf + col] = true;
		bit2[row * mf + col] = true;
		bit3[row * mf + col] = true;
		bit4[row * mf + col] = false;
		//tracept[row][col] = 'e';
	      } else {
		bit1[row * mf + col] = true;
		bit2[row * mf + col] = false;
		bit3[row * mf + col] = true;
		bit4[row * mf + col] = false;
		//tracept[row][col] = 'f';
	      }
	    } else {
	      if (verval == 's') {
		bit1[row * mf + col] = false;
		bit2[row * mf + col] = true;
		bit3[row * mf + col] = true;
		bit4[row * mf + col] = false;
		//tracept[row][col] = 'g';
	      } else {
		bit1[row * mf + col] = false;
		bit2[row * mf + col] = false;
		bit3[row * mf + col] = true;
		bit4[row * mf + col] = false;
		//tracept[row][col] = 'h';
	      }
	    }
	  } else {
	    if (hozval == 's') {
	      if (verval == 's') {
		bit1[row * mf + col] = true;
		bit2[row * mf + col] = true;
		bit3[row * mf + col] = false;
		bit4[row * mf + col] = true;
		//tracept[row][col] = 'i';
	      } else {
		bit1[row * mf + col] = true;
		bit2[row * mf + col] = false;
		bit3[row * mf + col] = false;
		bit4[row * mf + col] = true;
		//tracept[row][col] = 'j';
	      }
	    } else {
	      if (verval == 's') {
		bit1[row * mf + col] = false;
		bit2[row * mf + col] = true;
		bit3[row * mf + col] = false;
		bit4[row * mf + col] = true;
		//tracept[row][col] = 'k';
	    } else {
		bit1[row * mf + col] = false;
		bit2[row * mf + col] = false;
		bit3[row * mf + col] = false;
		bit4[row * mf + col] = true;
		//tracept[row][col] = 'l';
	      }
	    }
	  }
	}
      }
    }


    // Trace-back using pointers
    std::size_t row = m;
    std::size_t col = n;
    char lastMatrix = 's';
    typedef std::vector<char> TTrace;
    TTrace btr;
    while ((row>0) || (col>0)) {
      if (lastMatrix == 's') {
	if (bit3[row * mf + col]) lastMatrix = 'h';
	else if (bit4[row * mf + col]) lastMatrix = 'v';
	//if ((tracept[row][col] == 'e') || (tracept[row][col] == 'f') || (tracept[row][col] == 'g') || (tracept[row][col] == 'h')) lastMatrix = 'h';
	//else if ((tracept[row][col] == 'i') || (tracept[row][col] == 'j') || (tracept[row][col] == 'k') || (tracept[row][col] == 'l')) lastMatrix = 'v';
	else {
	  --row;
	  --col;
	  btr.push_back('s');
	}
      } else if (lastMatrix == 'h') {
	if (bit1[row * mf + col]) lastMatrix = 's';
	//if ((tracept[row][col] == 'a') || (tracept[row][col] == 'b') || (tracept[row][col] == 'e') || (tracept[row][col] == 'f') || (tracept[row][col] == 'i') || (tracept[row][col] == 'j')) lastMatrix = 's';
	--col;
	btr.push_back('h');
      } else if (lastMatrix == 'v') {
	if (bit2[row * mf + col]) lastMatrix = 's';
	//if ((tracept[row][col] == 'a') || (tracept[row][col] == 'c') || (tracept[row][col] == 'e') || (tracept[row][col] == 'g') || (tracept[row][col] == 'i') || (tracept[row][col] == 'k')) lastMatrix = 's';
	--row;
	btr.push_back('v');
      }
    }

    // Create alignment
    _createAlignment(btr, a1, a2, align);

    // Score
    return s[n];
  }

  template<typename TAlign1, typename TAlign2, typename TAlign, typename TAlignConfig>
  inline int
  gotoh(TAlign1 const& a1, TAlign2 const& a2, TAlign& align, TAlignConfig const& ac) 
  {
    DnaScore<int> dnasc;    
    return gotoh(a1, a2, align, ac, dnasc);
  }

  template<typename TAlign1, typename TAlign2, typename TAlign>
  inline int
  gotoh(TAlign1 const& a1, TAlign2 const& a2, TAlign& align)
  {
    AlignConfig<false, false> ac;
    return gotoh(a1, a2, align, ac);
  }

}

#endif
