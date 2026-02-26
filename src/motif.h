#ifndef MOTIF_H
#define MOTIF_H

#include <limits>

#include <boost/multi_array.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
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

    std::string matrixId;
    std::string symbol;
  };

  struct Pwm {
    typedef boost::multi_array<double, 2> T2DArray;
    T2DArray matrix;

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

  inline std::string
  _maxSimpleMotif(Pwm const& pwm) {
    std::string motif;
    for(uint32_t j = 0; j < pwm.matrix.shape()[1]; ++j) {
      double maxJ = pwm.matrix[0][j];
      uint32_t maxI = 0;
      for(uint32_t i = 1; i < 4; ++i) {
	if (pwm.matrix[i][j] > maxJ) {
	  maxJ = pwm.matrix[i][j];
	  maxI = i;
	}
      }
      if (maxI == 0) motif += 'A';
      else if (maxI == 1) motif += 'C';
      else if (maxI == 2) motif += 'G';
      else if (maxI == 3) motif += 'T';
    }
    return motif;
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
    out.matrix.resize(boost::extents[4][pfm.matrix.shape()[1]]);
    for(uint32_t i = 0; i < 4; ++i) {
      uint32_t r = out.matrix.shape()[1] - 1;
      for(uint32_t j = 0; j < out.matrix.shape()[1]; ++j, --r) {
	out.matrix[(3-i)][r] = pfm.matrix[i][j];
      }
    }
  }

  template<typename TConfig, typename TBitSet, typename TMotifHits>
  inline void
  scorePwm(TConfig const& c, char const* seq, TBitSet const& evalPos, Pwm const& inpwm, std::string const& tname, TMotifHits& mh, boost::iostreams::filtering_ostream& dataOut) {
    Pwm pwmFwd(inpwm);
    scale(pwmFwd);
    Pwm pwmRev;
    revComp(inpwm, pwmRev);
    scale(pwmRev);
    int32_t motiflen = pwmFwd.matrix.shape()[1];
    int32_t lastHit = -(motiflen + 1);
    //std::cerr << _minScore(pwmFwd) << ',' << _maxScore(pwmFwd) << std::endl;
    //std::cerr << _minScore(pwmRev) << ',' << _maxScore(pwmRev) << std::endl;
    
    // Parse sequence
    for(uint32_t pos = 0; pos < evalPos.size() - motiflen + 1; ++pos) {
      if (evalPos[pos]) {
	std::string ref = boost::to_upper_copy(std::string(seq + pos, seq + pos + motiflen));
	double scoreFwd = 0;
	double scoreRev = 0;
	int32_t k = 0;
	for(; k < motiflen; ++k) {
	  if (evalPos[pos+k]) {
	    int32_t n = 4;
	    if (ref[k] == 'A') n = 0;
	    else if (ref[k] == 'C') n = 1;
	    else if (ref[k] == 'G') n = 2;
	    else if (ref[k] == 'T') n = 3;
	    if (n < 4) {
	      scoreFwd += pwmFwd.matrix[n][k];
	      scoreRev += pwmRev.matrix[n][k];
	    } else break;
	  } else break;
	}
	if ((k == motiflen) && ((scoreFwd > c.motifScoreQuantile) || (scoreRev > c.motifScoreQuantile))) {
	  //std::cerr << "Genom:" << ref << "," << ref << std::endl;
	  //std::cerr << "Query:" << _maxSimpleMotif(pwmFwd) << "," << _maxSimpleMotif(pwmRev) << std::endl;
	  if ((c.overlappingHits) || (lastHit + motiflen < (int32_t) pos)) {
	    mh.push_back(pos);
	    lastHit = pos;
	    if (c.motifPosOut) {
	      if (scoreFwd > c.motifScoreQuantile) {
		dataOut << tname << "\t" << (pos + 1) << "\t" << pos + motiflen << "\t" << inpwm.symbol << "\t+\t" << scoreFwd << std::endl;
	      }
	      if (scoreRev > c.motifScoreQuantile) {
		dataOut << tname << "\t" << (pos + 1) << "\t" << pos + motiflen << "\t" << inpwm.symbol << "\t-\t" << scoreRev << std::endl;
	      }
	    }
	  }
	}
      }
    }
  }

  template<typename TConfig>
  inline bool
  parseJasparPfm(TConfig const& c, std::vector<Pfm>& pfms) {
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
    int32_t acgt = 0;
    int32_t id = 0;
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

  template<typename TConfig, typename TGenomicRegions, typename TMotifIds>
  inline int32_t
  parseJasparAll(TConfig const& c, TGenomicRegions& overlappingRegions, TMotifIds& motifIds) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Motif file parsing" << std::endl;

    // Generate PWMs
    std::vector<Pwm> pwms;
    parseJasparPwm(c, pwms);

    // Motif search
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Motif search" << std::endl;

    // Output motif positions
    boost::iostreams::filtering_ostream dataOut;
    if (c.motifPosOut) {
      dataOut.push(boost::iostreams::gzip_compressor());
      dataOut.push(boost::iostreams::file_sink(c.outpos.string(), std::ios_base::out | std::ios_base::binary));
      dataOut << "chr\tstart\tend\tid\tstrand\tquantile" << std::endl;
    }
    
    // Iterate chromosomes
    faidx_t* fai = fai_load(c.genome.string().c_str());
    char* seq = NULL;
    for(int32_t refIndex=0; refIndex < (int32_t) c.nchr.size(); ++refIndex) {

      // Chromosome name and length
      std::string tname = "NA";
      for(typename TConfig::TChrMap::const_iterator itChr = c.nchr.begin(); itChr != c.nchr.end(); ++itChr) {
	if (refIndex == itChr->second) {
	  tname = itChr->first;
	}
      }
      int32_t seqlen = faidx_seq_len(fai, tname.c_str());

      // Pre-process bed file so we can speed-up motif search
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet evalPos(seqlen, false);
      std::ifstream chrFile(c.infile.string().c_str(), std::ifstream::in);
      if (chrFile.is_open()) {
	while (chrFile.good()) {
	  std::string chrFromFile;
	  getline(chrFile, chrFromFile);
	  typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	  boost::char_separator<char> sep(" \t,;");
	  Tokenizer tokens(chrFromFile, sep);
	  Tokenizer::iterator tokIter = tokens.begin();
	  if (tokIter!=tokens.end()) {
	    std::string chrName = *tokIter++;
	    if (c.nchr.find(chrName)->second != refIndex) continue;
	    int32_t start = boost::lexical_cast<int32_t>(*tokIter++);
	    int32_t end = boost::lexical_cast<int32_t>(*tokIter++);
	    std::string name = "NA";
	    if (start >= end) continue;  // Bed has right-open intervals
	    int32_t realstart = std::max(0, start - c.maxDistance);
            int32_t realend = std::min(seqlen, end + c.maxDistance);
	    for(int32_t i = realstart; i < realend; ++i) evalPos[i] = true;
	  }
	}
	chrFile.close();
      }

      // Anything to annotate on this chromosome?
      if (evalPos.count()) {
	seqlen = -1;
	seq = faidx_fetch_seq(fai, tname.c_str(), 0, faidx_seq_len(fai, tname.c_str()) + 1, &seqlen);

	// Blacklist Ns
	for(int32_t i = 0; i < seqlen; ++i) {
	  if ((seq[i] == 'n') || (seq[i] == 'N')) evalPos[i] = false;
	}
	
	// Score PWMs
	for(uint32_t i = 0; i<pwms.size(); ++i) {
	  typedef std::vector<int32_t> TMotifHits;
	  TMotifHits mh;
	  scorePwm(c, seq, evalPos, pwms[i], tname, mh, dataOut);

	  int32_t motiflen = pwms[i].matrix.shape()[1];
	  for(uint32_t hit = 0; hit < mh.size(); ++hit) {
	    _insertInterval(overlappingRegions[refIndex], mh[hit], mh[hit] + motiflen, '*', i, 0);
	  }
	}

	// Clean-up
	if (seq != NULL) free(seq);
      }
    }
    // Close gzipped motif positions
    if (c.motifPosOut) dataOut.pop();
    
    // Assign Motif Ids
    for(uint32_t i = 0; i<pwms.size(); ++i) motifIds.push_back(pwms[i].symbol);

    return motifIds.size();
  }

  
  template<typename TConfig, typename TGenomicRegions, typename TMotifIds>
  inline int32_t
  parseJaspar(TConfig const& c, TGenomicRegions& gRegions, TMotifIds& motifIds) {
    typedef typename TGenomicRegions::value_type TChromosomeRegions;

    // Overlapping intervals for each label
    TGenomicRegions overlappingRegions;
    overlappingRegions.resize(gRegions.size(), TChromosomeRegions());
    parseJasparAll(c, overlappingRegions, motifIds);
    
    // Make intervals non-overlapping for each label
    for(uint32_t refIndex = 0; refIndex < overlappingRegions.size(); ++refIndex) {
      // Sort by ID
      std::sort(overlappingRegions[refIndex].begin(), overlappingRegions[refIndex].end(), SortIntervalLabel<IntervalLabel>());
      int32_t runningId = -1;
      char runningStrand = '*';
      typedef boost::icl::interval_set<uint32_t> TIdIntervals;
      typedef typename TIdIntervals::interval_type TIVal;
      TIdIntervals idIntervals;
      for(uint32_t i = 0; i < overlappingRegions[refIndex].size(); ++i) {
	if (overlappingRegions[refIndex][i].lid != runningId) {
	  for(typename TIdIntervals::iterator it = idIntervals.begin(); it != idIntervals.end(); ++it) {
	    gRegions[refIndex].push_back(IntervalLabel(it->lower(), it->upper(), runningStrand, runningId));
	  }
	  idIntervals.clear();
	  runningId = overlappingRegions[refIndex][i].lid;
	  runningStrand = overlappingRegions[refIndex][i].strand;
	}
	idIntervals.insert(TIVal::right_open(overlappingRegions[refIndex][i].start, overlappingRegions[refIndex][i].end));
      }
      // Process last id
      for(typename TIdIntervals::iterator it = idIntervals.begin(); it != idIntervals.end(); ++it) gRegions[refIndex].push_back(IntervalLabel(it->lower(), it->upper(), runningStrand, runningId));
    }

    return motifIds.size();
  }
  
}

#endif
