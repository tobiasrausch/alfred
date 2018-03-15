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
#include <boost/math/distributions/binomial.hpp>

#include <htslib/sam.h>


namespace bamstats
{

  struct Interval {
    int32_t start;
    int32_t end;
    
    Interval(int32_t s, int32_t e) : start(s), end(e) {}
  };


  inline double
  binomTest(uint32_t x, uint32_t n, double p) {
    boost::math::binomial binomialdist(n, p);
    double cutoff = pdf(binomialdist, x);
    double pval = 0.0;
    for(uint32_t k = 0; k <= n; ++k) {
      double p = pdf(binomialdist, k);
      if (p <= cutoff) pval +=p;
    }
    return pval;
  }


  
  inline unsigned hash_string(const char *s) {
    unsigned h = 37;
    while (*s) {
      h = (h * 54059) ^ (s[0] * 76963);
      s++;
    }
    return h;
  }


  struct IntervalLabel {
    int32_t start;
    int32_t end;
    char strand;
    int32_t lid;

    explicit IntervalLabel(int32_t s) : start(s), end(s+1), strand('*'), lid(-1) {}
    IntervalLabel(int32_t s, int32_t e, char t, int32_t l) : start(s), end(e), strand(t), lid(l) {}
  };

  template<typename TRecord>
  struct SortIntervalLabel : public std::binary_function<TRecord, TRecord, bool> {
    inline bool operator()(TRecord const& s1, TRecord const& s2) const {
      return s1.lid < s2.lid;
    }
  };

  template<typename TRecord>
  struct SortIntervalStart : public std::binary_function<TRecord, TRecord, bool> {
    inline bool operator()(TRecord const& s1, TRecord const& s2) const {
      return s1.start < s2.start;
    }
  };

  
  inline std::size_t hash_pair(bam1_t* rec) {
    std::size_t seed = hash_string(bam_get_qname(rec));
    boost::hash_combine(seed, rec->core.tid);
    boost::hash_combine(seed, rec->core.pos);
    boost::hash_combine(seed, rec->core.mtid);
    boost::hash_combine(seed, rec->core.mpos);
    return seed;
  }

  inline int32_t
  homopolymerContext(std::string const& s, int32_t idx, int32_t homlen) {
    for(int32_t i = std::max(0, idx - (homlen - 1)); i <= (idx + 1); ++i) {
      if (i + homlen <= (int32_t) s.size()) {
	bool hompoly = true;
	for(int32_t k = i + 1; k < i + homlen; ++k) {
	  if (s[k] != s[i]) {
	    hompoly = false;
	    break;
	  }
	}
	if (hompoly) {
	  if (s[i] == 'A') return 0;
	  else if (s[i] == 'C') return 1;
	  else if (s[i] == 'G') return 2;
	  else if (s[i] == 'T') return 3;
	  else if (s[i] == 'N') return 4;
	}
      }
    }
    return 5; // None
  }

  inline std::size_t hash_pair_mate(bam1_t* rec) {
    std::size_t seed = hash_string(bam_get_qname(rec));
    boost::hash_combine(seed, rec->core.mtid);
    boost::hash_combine(seed, rec->core.mpos);
    boost::hash_combine(seed, rec->core.tid);
    boost::hash_combine(seed, rec->core.pos);
    return seed;
  }

  inline bool is_gff3(boost::filesystem::path const& f) {
    std::ifstream in(f.string().c_str());
    if (!in) return false;
    in.close();

    std::ifstream file(f.string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string gline;
    std::getline(instream, gline);
    bool gff = false;
    if ((gline.size()) && (gline == "##gff-version 3")) gff = true;
    file.close();
    return gff;
  }
    
  
  inline bool is_gz(boost::filesystem::path const& f) {
    std::ifstream in(f.string().c_str());
    if (!in) return false;
    in.close();
    
    boost::iostreams::filtering_istream gzin;
    gzin.push(boost::iostreams::gzip_decompressor());
    gzin.push(boost::iostreams::file_source(f.string().c_str()), std::ios_base::in | std::ios_base::binary);
    try {
      char c;
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
      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CDEL) || (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP)) alen += bam_cigar_oplen(cigar[i]);
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
