/*
============================================================================
Alfred: BAM alignment statistics
============================================================================
Copyright (C) 2017 Tobias Rausch

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

#ifndef GTF_H
#define GTF_H

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/algorithm/string.hpp>

#include <htslib/sam.h>

#include "util.h"

namespace bamstats
{


  struct IntervalLabel {
    int32_t start;
    int32_t end;
    char strand;
    int32_t lid;

    IntervalLabel(int32_t s, int32_t e, char t, int32_t l) : start(s), end(e), strand(t), lid(l) {}
  };


  template<typename TGenomicRegions, typename TGeneIds>
  inline bool
    parseGTF(bam_hdr_t* hdr, boost::filesystem::path const& gtf, std::string const& feature, std::string const& attribute, TGenomicRegions& gRegions, TGeneIds& geneIds) {
    if (!is_gz(gtf)) {
      std::cerr << "GTF file is not gzipped!" << std::endl;
      return false;
    }
    typedef std::map<std::string, int32_t> TIdMap;
    TIdMap idMap;
    std::ifstream file(gtf.string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string gline;
    while(std::getline(instream, gline)) {
      if ((gline.size()) && (gline[0] == '#')) continue;
      typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
      boost::char_separator<char> sep("\t");
      Tokenizer tokens(gline, sep);
      Tokenizer::iterator tokIter = tokens.begin();
      if (tokIter==tokens.end()) {
	std::cerr << "Empty line in GTF file!" << std::endl;
	return false;
      }
      std::string chrName=*tokIter++;
      int32_t chrid = bam_name2id(hdr, chrName.c_str());
      if (chrid < 0) continue;
      if (tokIter == tokens.end()) {
	std::cerr << "Corrupted GTF file!" << std::endl;
	return false;
      }
      ++tokIter;
      if (tokIter == tokens.end()) {
	std::cerr << "Corrupted GTF file!" << std::endl;
	return false;
      }
      std::string ft = *tokIter++;
      if (ft == feature) {
	if (tokIter != tokens.end()) {
	  int32_t start = boost::lexical_cast<int32_t>(*tokIter++);
	  int32_t end = boost::lexical_cast<int32_t>(*tokIter++);
	  ++tokIter; // score
	  if (tokIter == tokens.end()) {
	    std::cerr << "Corrupted GTF file!" << std::endl;
	    return false;
	  }
	  char strand = boost::lexical_cast<char>(*tokIter++);
	  ++tokIter; // frame
	  std::string attr = *tokIter;
	  boost::char_separator<char> sepAttr(";");
	  Tokenizer attrTokens(attr, sepAttr);
	  for(Tokenizer::iterator attrIter = attrTokens.begin(); attrIter != attrTokens.end(); ++attrIter) {
	    std::string keyval = *attrIter;
	    boost::trim(keyval);
	    boost::char_separator<char> sepKeyVal(" ");
	    Tokenizer kvTokens(keyval, sepKeyVal);
	    Tokenizer::iterator kvTokensIt = kvTokens.begin();
	    std::string key = *kvTokensIt++;
	    if (key == attribute) {
	      std::string val = *kvTokensIt;
	      int32_t idval = geneIds.size();
	      typename TIdMap::const_iterator idIter = idMap.find(val);
	      if (idIter == idMap.end()) {
		idMap.insert(std::make_pair(val, idval));
		geneIds.push_back(val);
	      } else idval = idIter->second;
	      gRegions[chrid].push_back(IntervalLabel(start, end, strand, idval));
	    }
	  }
	}
      }
    }
    return true;
  }

}

#endif
