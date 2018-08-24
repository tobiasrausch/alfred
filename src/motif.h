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

#ifndef MOTIF_H
#define MOTIF_H

#include <limits>

#include <boost/multi_array.hpp>
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

  struct Pfm {
    typedef boost::multi_array<int32_t, 2> T2DArray;
    T2DArray matrix;

    int32_t id;
    std::string matrixId;
    std::string symbol;
  };

  struct Pwm {
    typedef boost::multi_array<double, 2> T2DArray;
    T2DArray matrix;

    int32_t id;
    std::string matrixId;
    std::string symbol;
  };

  inline double
  _minScore(Pwm const& pwm) {
    double ms = 0;
    for(uint32_t j = 0; j < pwm.matrix.shape()[1]; ++j) {
      double minJ = pwm.matrix[0][j];
      for(uint32_t i = 1; i < 4; ++i) {
	if (pwm.matrix[i][j] < minJ) minJ = pwm.matrix[i][j];
      }
      ms += minJ;
    }
    return ms;
  }

  inline double
  _maxScore(Pwm const& pwm) {
    double ms = 0;
    for(uint32_t j = 0; j < pwm.matrix.shape()[1]; ++j) {
      double maxJ = pwm.matrix[0][j];
      for(uint32_t i = 1; i < 4; ++i) {
	if (pwm.matrix[i][j] > maxJ) maxJ = pwm.matrix[i][j];
      }
      ms += maxJ;
    }
    return ms;
  }

  inline void
  scale(Pwm& pwm) {
    double minsc = _minScore(pwm);
    double maxsc = _maxScore(pwm);
    double cols = pwm.matrix.shape()[1];
    for(uint32_t i = 0; i < 4; ++i) {
      for(uint32_t j = 0; j < pwm.matrix.shape()[1]; ++j) {
	pwm.matrix[i][j] = ((pwm.matrix[i][j] - minsc / cols) / (maxsc - minsc));
      }
    }
  }
  
  inline void
  convert(Pfm const& pfm, Pwm& pwm, std::vector<double> const& bg, double const pc) {
    pwm.matrixId = pfm.matrixId;
    pwm.symbol = pfm.symbol;
    pwm.id = pfm.id;
    pwm.matrix.resize(boost::extents[4][pfm.matrix.shape()[1]]);
    double totalBg = 0;
    for(uint32_t i = 0; i < 4; ++i) totalBg += bg[i];
    for(uint32_t j = 0; j < pwm.matrix.shape()[1]; ++j) {
      int32_t total = 0;
      for(uint32_t i = 0; i < 4; ++i) total += pfm.matrix[i][j];
      for(uint32_t i = 0; i < 4; ++i) {
	pwm.matrix[i][j] = ((double) pfm.matrix[i][j] + bg[i] * pc) / ((double) total + totalBg * pc);
	pwm.matrix[i][j] = std::log(pwm.matrix[i][j] / (bg[i] / totalBg)) / log(2);
      }
    }
  }

  inline void
  convert(Pfm const& pfm, Pwm& pwm, double const pc) {
    std::vector<double> bg;
    bg.push_back(0.25);
    bg.push_back(0.25);
    bg.push_back(0.25);
    bg.push_back(0.25);
    convert(pfm, pwm, bg, pc);
  }

  inline void
  convert(Pfm const& pfm, Pwm& pwm) {
    convert(pfm, pwm, 0.8);
  }
  
  template<typename TPositionMatrix>
  inline void
  revComp(TPositionMatrix const& pfm, TPositionMatrix& out) {
    out.matrixId = pfm.matrixId;
    out.symbol = pfm.symbol;
    out.id = pfm.id;
    out.matrix.resize(boost::extents[4][pfm.matrix.shape()[1]]);
    for(uint32_t i = 0; i < 4; ++i) {
      uint32_t r = out.matrix.shape()[1] - 1;
      for(uint32_t j = 0; j < out.matrix.shape()[1]; ++j, --r) {
	out.matrix[(3-i)][r] = pfm.matrix[i][j];
      }
    }
  }

  template<typename TMotifHits>
  inline std::pair<int32_t, int32_t>
  scorePwm(std::vector<uint8_t> const& ref, Pwm const& pwm, double const fraction, TMotifHits& mh) {
    int32_t hitsFwd = 0;
    int32_t hitsRev = 0;
    int32_t motiflen = pwm.matrix.shape()[1];
    double maxscore = _maxScore(pwm);
    double minscore = _minScore(pwm);
    double threshold = minscore + (maxscore - minscore) * fraction;
    Pwm pwmRev;
    revComp(pwm, pwmRev);
    
    // Parse sequence
    double scoreFwd = 0;
    double scoreRev = 0;
    typedef std::vector<uint8_t> TNucMap;
    TNucMap::const_iterator n = ref.begin();
    int32_t k = 0;
    for(TNucMap::const_iterator itNuc = ref.begin(); itNuc < ref.end() - motiflen + 1; ++itNuc) {
      scoreFwd = 0;
      scoreRev = 0;
      n = itNuc;
      for(k = 0; k < motiflen; ++k, ++n) {
	if (*n < 4) {
	  scoreFwd += pwm.matrix[*n][k];
	  scoreRev += pwmRev.matrix[*n][k];
	} else break;
      }
      if (scoreFwd > threshold) ++hitsFwd;
      if (scoreRev > threshold) ++hitsRev;
    }
    return std::make_pair(hitsFwd, hitsRev);
  }

  template<typename TConfig>
  inline bool
  parseJasparPfm(TConfig const& c, std::vector<Pfm>& pfms) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "JASPAR parsing" << std::endl;

    // Check gzip
    if (!is_gz(c.motifFile)) {
      std::cerr << "JASPAR file is not gzipped!" << std::endl;
      return false;
    }

    // Parse JASPAR
    std::ifstream file(c.motifFile.string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string gline;
    int32_t id = 0;
    int32_t acgt = 0;
    while(std::getline(instream, gline)) {
      // Header
      if ((gline.size()) && (gline[0] == '>')) {
	id = pfms.size();
	pfms.resize(id+1);
	gline = gline.substr(1);
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t");
	Tokenizer tokens(gline, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter != tokens.end()) {
	  pfms[id].id = id;
	  pfms[id].matrixId = *tokIter++;
	  pfms[id].symbol = "NA";
	  if (tokIter != tokens.end()) {
	    pfms[id].symbol = *tokIter++;
	  }
	}
	acgt = 0;
      } else {
	if ((gline.size()) && ((gline[0] == 'A') || (gline[0] == 'C') || (gline[0] == 'G') || (gline[0] == 'T'))) {
	  // JASPAR format
	  typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	  boost::char_separator<char> sep("[");
	  Tokenizer tokens(gline, sep);
	  Tokenizer::iterator tokIter = tokens.begin();
	  if ((tokIter!=tokens.end()) && (++tokIter != tokens.end())) {
	    gline = *tokIter;
	    boost::char_separator<char> sep2("]");
	    Tokenizer tokens2(gline, sep2);
	    Tokenizer::iterator tokIter2 = tokens2.begin();
	    if (tokIter2 != tokens2.end()) {
	      gline = *tokIter2;
	    } else {
	      std::cerr << "JASPAR cannot be parsed!" << std::endl;
	      return false;
	    }
	  } else {
	    std::cerr << "JASPAR cannot be parsed!" << std::endl;
	    return false;
	  }
	}

	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t");
	Tokenizer tokens(gline, sep);
	if (acgt == 0) { 
	  int32_t lenMotif = 0;
	  for(Tokenizer::iterator tokIter = tokens.begin(); tokIter!=tokens.end(); ++tokIter) ++lenMotif;
	  pfms[id].matrix.resize(boost::extents[4][lenMotif]);
	}
	uint32_t col = 0;
	for(Tokenizer::iterator tokIter = tokens.begin(); tokIter!=tokens.end(); ++tokIter, ++col) pfms[id].matrix[acgt][col] = boost::lexical_cast<int32_t>(*tokIter);

	// Debug code
	//if (acgt == 3) {
	//std::cout << ">" << pfms[id].matrixId << ',' << pfms[id].symbol << std::endl;
	//for(uint32_t i = 0; i < pfms[id].matrix.shape()[0]; ++i) {
	// for(uint32_t j = 0; j < pfms[id].matrix.shape()[1]; ++j) {
	//    std::cerr << pfms[id].matrix[i][j] << ',';
	//  }
	//  std::cerr << std::endl;
	//}
	//}
	
	++acgt;
      }

    }
    dataIn.pop();
    return true;
  }

  template<typename TConfig>
  inline bool
  parseJasparPwm(TConfig const& c, std::vector<Pwm>& pwms) {
    std::vector<Pfm> pfms;
    if (!parseJasparPfm(c, pfms)) return false;
    pwms.resize(pfms.size(), Pwm());
    for(uint32_t i = 0; i < pfms.size(); ++i) {
      convert(pfms[i], pwms[i]);
    }
    return true;
  }
  
}

#endif
