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

#ifndef UTIL_H
#define UTIL_H

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include <htslib/sam.h>


namespace bamstats
{

  struct Interval {
    int32_t start;
    int32_t end;
    
    Interval(int32_t s, int32_t e) : start(s), end(e) {}
  };
  
  inline bool is_gz(boost::filesystem::path const& f) {
    std::ifstream in(f.string().c_str());
    if (!in) return false;
    in.close();
    
    boost::iostreams::filtering_istream gzin;
    gzin.push(boost::iostreams::gzip_decompressor());
    gzin.push(boost::iostreams::file_source(f.string().c_str()), std::ios_base::in | std::ios_base::binary);
    char c;
    try {
      gzin >> c;
    } catch (boost::iostreams::gzip_error& e) {
      gzin.pop();
      return false;
    }
    gzin.pop();
    return true;
  }
      
    

  // F+ 0
  // F- 1
  // R+ 2
  // R- 3
  inline uint8_t layout(bam1_t const* rec) {
    if (rec->core.flag & BAM_FREAD1) {
      if (!(rec->core.flag & BAM_FREVERSE)) {
	if (!(rec->core.flag & BAM_FMREVERSE)) return (rec->core.pos < rec->core.mpos) ? 0 : 1;
	else return (rec->core.pos < rec->core.mpos) ? 2 : 3;
      } else {
	if (!(rec->core.flag & BAM_FMREVERSE)) return (rec->core.pos > rec->core.mpos) ? 2 : 3;
	else return (rec->core.pos > rec->core.mpos) ? 0 : 1;
      }
    } else {
      if (!(rec->core.flag & BAM_FREVERSE)) {
	if (!(rec->core.flag & BAM_FMREVERSE)) return (rec->core.pos < rec->core.mpos) ? 1 : 0;
	else return (rec->core.pos < rec->core.mpos) ? 2 : 3;
      } else {
	if (!(rec->core.flag & BAM_FMREVERSE)) return (rec->core.pos > rec->core.mpos) ? 2 : 3;
	else return (rec->core.pos > rec->core.mpos) ? 1 : 0;
      }
    }
  }

  
  inline uint32_t alignmentLength(bam1_t const* rec) {
    uint32_t* cigar = bam_get_cigar(rec);
    uint32_t alen = 0;
    for (uint32_t i = 0; i < rec->core.n_cigar; ++i)
      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CDEL)) alen += bam_cigar_oplen(cigar[i]);
    return alen;
  }

  inline uint32_t
  lastAlignedPosition(bam1_t const* rec) {
    return rec->core.pos + alignmentLength(rec);
  }

  inline uint32_t halfAlignmentLength(bam1_t* rec) {
    return (alignmentLength(rec) / 2);
  }

  inline void
  reverseComplement(std::string& sequence) 
  {
    std::string rev = boost::to_upper_copy(std::string(sequence.rbegin(), sequence.rend()));
    std::size_t i = 0;
    for(std::string::iterator revIt = rev.begin(); revIt != rev.end(); ++revIt, ++i) {
      switch (*revIt) {
      case 'A': sequence[i]='T'; break;
      case 'C': sequence[i]='G'; break;
      case 'G': sequence[i]='C'; break;
      case 'T': sequence[i]='A'; break;
      case 'N': sequence[i]='N'; break;
      default: break;
      }
    }
  }

  inline bool
  getSMTag(std::string const& header, std::string const& fileName, std::string& sampleName) {
    std::set<std::string> smIdentifiers;
    std::string delimiters("\n");
    typedef std::vector<std::string> TStrParts;
    TStrParts lines;
    boost::split(lines, header, boost::is_any_of(delimiters));
    TStrParts::const_iterator itH = lines.begin();
    TStrParts::const_iterator itHEnd = lines.end();
    bool rgPresent = false;
    for(;itH!=itHEnd; ++itH) {
      if (itH->find("@RG")==0) {
	std::string delim("\t ");
	TStrParts keyval;
	boost::split(keyval, *itH, boost::is_any_of(delim));
	TStrParts::const_iterator itKV = keyval.begin();
	TStrParts::const_iterator itKVEnd = keyval.end();
	for(;itKV != itKVEnd; ++itKV) {
	  size_t sp = itKV->find(":");
	  if (sp != std::string::npos) {
	    std::string field = itKV->substr(0, sp);
	    if (field == "SM") {
	      rgPresent = true;
	      std::string rgSM = itKV->substr(sp+1);
	      smIdentifiers.insert(rgSM);
	    }
	  }
	}
      }
    }
    if (!rgPresent) {
      sampleName = fileName;
      return true;
    } else if (smIdentifiers.size() == 1) {
      sampleName = *(smIdentifiers.begin());
      return true;
    } else {
      sampleName = "";
      return false;
    }
  }

  inline void
  getRGs(std::string const& header, std::set<std::string>& rgs) {
    // Get read groups
    std::string delimiters("\n");
    typedef std::vector<std::string> TStrParts;
    TStrParts lines;
    boost::split(lines, header, boost::is_any_of(delimiters));
    TStrParts::const_iterator itH = lines.begin();
    TStrParts::const_iterator itHEnd = lines.end();
    bool rgPresent = false;
    for(;itH!=itHEnd; ++itH) {
      if (itH->find("@RG")==0) {
	std::string delim("\t ");
	TStrParts keyval;
	boost::split(keyval, *itH, boost::is_any_of(delim));
	TStrParts::const_iterator itKV = keyval.begin();
	TStrParts::const_iterator itKVEnd = keyval.end();
	for(;itKV != itKVEnd; ++itKV) {
	  size_t sp = itKV->find(":");
	  if (sp != std::string::npos) {
	    std::string field = itKV->substr(0, sp);
	    if (field == "ID") {
	      rgPresent = true;
	      std::string rgID = itKV->substr(sp+1);
	      rgs.insert(rgID);
	    }
	  }
	}
      }
    }
    if (!rgPresent) rgs.insert("DefaultLib");
  }


  template<typename TVector>
  inline int32_t
  medianFromHistogram(TVector const& vec) {
    int64_t tc = 0;
    for(typename TVector::const_iterator it = vec.begin(); it != vec.end(); ++it) tc += *it;
    int64_t medind = tc / 2;
    tc = 0;
    for(int32_t i = 0; i < (int32_t) vec.size(); ++i) {
      tc += (int64_t) vec[i];
      if (tc >= medind) return i;
    }
    return 0;
  }

  template<typename TVector>
  inline double
  meanFromHistogram(TVector const& vec) {
    int64_t tc = 0;
    for(typename TVector::const_iterator it = vec.begin(); it != vec.end(); ++it) tc += *it;
    int64_t mean = 0;
    for(int32_t i = 0; i < (int32_t) vec.size(); ++i) mean += (int64_t) (vec[i]) * (int64_t) (i);
    return (double) mean / (double) tc;
  }

  template<typename TVector>
  inline double
  sdFromHistogram(TVector const& vec) {
    double mu = meanFromHistogram(vec);
    int64_t tc = 0;
    for(typename TVector::const_iterator it = vec.begin(); it != vec.end(); ++it) tc += *it;
    double sd = 0;
    for(int32_t i = 0; i < (int32_t) vec.size(); ++i) sd += (double) (vec[i]) * ((double) (i) - mu) * ((double) (i) - mu);
    return std::sqrt((double) sd / (double) tc);
  }


}

#endif
