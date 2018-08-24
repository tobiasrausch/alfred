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
    T2DArray pfm;

    std::string matrixId;
    std::string symbol;
  };


  template<typename TConfig>
  inline bool
  parseJaspar(TConfig const& c, std::vector<Pfm>& pfms) {
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
	  pfms[id].pfm.resize(boost::extents[4][lenMotif]);
	}
	uint32_t col = 0;
	for(Tokenizer::iterator tokIter = tokens.begin(); tokIter!=tokens.end(); ++tokIter, ++col) pfms[id].pfm[acgt][col] = boost::lexical_cast<int32_t>(*tokIter);

	// Debug code
	//if (acgt == 3) {
	//std::cout << ">" << id << "," << pfms[id].matrixId << ',' << pfms[id].symbol << std::endl;
	//for(uint32_t i = 0; i < pfms[id].pfm.shape()[0]; ++i) {
	//for(uint32_t j = 0; j < pfms[id].pfm.shape()[1]; ++j) {
	//std::cerr << pfms[id].pfm[i][j] << ',';
	//}
	//std::cerr << std::endl;
	//}
	//}

	
	++acgt;
      }

    }
    dataIn.pop();
    return true;
  }
  
}

#endif
