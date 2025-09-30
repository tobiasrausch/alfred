#ifndef CONSEDLIB_H
#define CONSEDLIB_H

#include <iostream>

#include "util.h"
#include "msa.h"
#include "edlib.h"

namespace bamstats
{

  inline void
  printAlignmentPretty(std::string const& query, std::string const& target, EdlibAlignMode const modeCode, EdlibAlignResult const& align) {
    int32_t tIdx = -1;
    int32_t qIdx = -1;
    if (modeCode == EDLIB_MODE_HW) {
        tIdx = align.endLocations[0];
        for (int32_t i = 0; i < align.alignmentLength; i++) {
            if (align.alignment[i] != EDLIB_EDOP_INSERT) tIdx--;
        }
    }
    std::cerr << std::endl;
    for (int start = 0; start < align.alignmentLength; start += 50) {
      std::cerr << "T: ";
      int32_t startTIdx = -1;
      for (int32_t j = start; ((j < start + 50) && (j < align.alignmentLength)); ++j) {
	if (align.alignment[j] == EDLIB_EDOP_INSERT) std::cerr << "-";
	else std::cerr << target[++tIdx];
	if (j == start) startTIdx = tIdx;
      }
      std::cerr << " (" << std::max(startTIdx, 0) << " - " << tIdx << ")" << std::endl;

      // match / mismatch
      std::cerr << ("   ");
      for (int32_t j = start; j < start + 50 && j < align.alignmentLength; j++) {
	if (align.alignment[j] == EDLIB_EDOP_MATCH) std::cerr <<  "|";
	else std::cerr << " ";
      }
      std::cerr << std::endl;

      // query
      std::cerr << "Q: ";
      int32_t startQIdx = qIdx;
      for (int32_t j = start; j < start + 50 && j < align.alignmentLength; j++) {
	if (align.alignment[j] == EDLIB_EDOP_DELETE) std::cerr << "-";
	else std::cerr << query[++qIdx];
	if (j == start) startQIdx = qIdx;
      }
      std::cerr << " ("<< std::max(startQIdx, 0) << " - " << qIdx << ")" << std::endl;
      std::cerr << std::endl;
    }
  }

  inline void
  printAlignment(std::string const& seqI, std::string const& seqJ, EdlibAlignMode const modeCode, EdlibAlignResult const& cigar) {
    int32_t tIdx = -1;
    int32_t qIdx = -1;
    uint32_t missingEnd = 0;
    uint32_t missingStart = 0;
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      tIdx = cigar.endLocations[0];
      if (tIdx < (int32_t) seqJ.size()) missingEnd = seqJ.size() - tIdx - 1;
      for (int32_t i = 0; i < cigar.alignmentLength; i++) {
	if (cigar.alignment[i] != EDLIB_EDOP_INSERT) tIdx--;
      }
      if (tIdx >= 0) missingStart = tIdx + 1;
      if (missingStart) {
	for (uint32_t j = 0; j < missingStart; ++j) std::cerr << '-';
      }
    }
    // seqI
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_DELETE) std::cerr << '-';
      else std::cerr << seqI[++qIdx];
    }
    // infix alignment, fix end
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      if (missingEnd) {
	for (uint32_t j = 0; j < missingEnd; ++j) std::cerr << '-';
      }
    }
    std::cerr << std::endl;
    // infix alignment, fix start
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      if (missingStart) {
	for (uint32_t j = 0; j < missingStart; ++j) std::cerr << seqJ[j];
      }
    }
    // seqJ
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_INSERT) std::cerr << '-';
      else std::cerr << seqJ[++tIdx];
    }
    // infix alignment, fix end
    if ((modeCode == EDLIB_MODE_HW) || (modeCode == EDLIB_MODE_SHW)) {
      if (missingEnd) {
	for (uint32_t j = 0; j < missingEnd; ++j) std::cerr << seqJ[++tIdx];
      }
    }
    std::cerr << std::endl;
  }

  template<typename TAlign>
  inline void
  convertAlignment(std::string const& query, TAlign& align, EdlibAlignMode const modeCode, EdlibAlignResult const& cigar) {
    // Input alignment
    TAlign alignIn;
    alignIn.resize(boost::extents[align.shape()[0]][align.shape()[1]]);
    for(uint32_t i = 0; i < align.shape()[0]; ++i) {
      for(uint32_t j = 0; j < align.shape()[1]; ++j) {
	alignIn[i][j] = align[i][j];
      }
    }
	
    // Create new alignment
    uint32_t seqPos = alignIn.shape()[0];
    int32_t tIdx = -1;
    int32_t qIdx = -1;
    uint32_t missingEnd = 0;
    uint32_t missingStart = 0;
    if (modeCode == EDLIB_MODE_HW) {
      tIdx = cigar.endLocations[0];
      if (tIdx < (int32_t) alignIn.shape()[1]) missingEnd = alignIn.shape()[1] - tIdx - 1;
      for (int32_t i = 0; i < cigar.alignmentLength; i++) {
	if (cigar.alignment[i] != EDLIB_EDOP_INSERT) tIdx--;
      }
      if (tIdx >= 0) missingStart = tIdx + 1;
    }
    align.resize(boost::extents[alignIn.shape()[0]+1][missingStart + cigar.alignmentLength + missingEnd]);

    // infix alignment, fix start
    if (modeCode == EDLIB_MODE_HW) {
      if (missingStart) {
	for (uint32_t j = 0; j < missingStart; ++j) {
	  for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j] = alignIn[seqIdx][j];
	  align[seqPos][j] = '-';
	}
      }
    }
    
    // target
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_INSERT) {
	for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j + missingStart] = '-';
      } else {
	++tIdx;
	for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j + missingStart] = alignIn[seqIdx][tIdx];
      }
    }

    // query
    for (int32_t j = 0; j < cigar.alignmentLength; ++j) {
      if (cigar.alignment[j] == EDLIB_EDOP_DELETE) align[seqPos][j + missingStart] = '-';
      else align[seqPos][j + missingStart] = query[++qIdx];
    }

    // infix alignment, fix end
    if (modeCode == EDLIB_MODE_HW) {
      if (missingEnd) {
	for (uint32_t j = cigar.alignmentLength + missingStart; j < cigar.alignmentLength + missingStart + missingEnd; ++j) {
	  ++tIdx;
	  for(uint32_t seqIdx = 0; seqIdx < seqPos; ++seqIdx) align[seqIdx][j] = alignIn[seqIdx][tIdx];
	  align[seqPos][j] = '-';
	}
      }
    }
  }

  template<typename TAlign>
  inline void
  consensusEdlib(TAlign const& align, std::string& cons) {
    typedef typename TAlign::index TAIndex;

    cons.resize(align.shape()[1]);
    for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
      std::vector<int32_t> count(5, 0); // ACGT-
      for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {
	if ((align[i][j] == 'A') || (align[i][j] == 'a')) ++count[0];
	else if ((align[i][j] == 'C') || (align[i][j] == 'c')) ++count[1];
	else if ((align[i][j] == 'G') || (align[i][j] == 'g')) ++count[2];
	else if ((align[i][j] == 'T') || (align[i][j] == 't')) ++count[3];
	else ++count[4];
      }
      uint32_t maxIdx = 0;
      uint32_t sndIdx = 1;
      if (count[maxIdx] < count[sndIdx]) {
	maxIdx = 1;
	sndIdx = 0;
      }
      for(uint32_t i = 2; i<5; ++i) {
	if (count[i] > count[maxIdx]) {
	  sndIdx = maxIdx;
	  maxIdx = i;
	}
	else if (count[i] > count[sndIdx]) {
	  sndIdx = i;
	}
      }
      if (2 * count[sndIdx] < count[maxIdx]) {
	switch (maxIdx) {
	case 0: cons[j] = 'A'; break;
	case 1: cons[j] = 'C'; break;
	case 2: cons[j] = 'G'; break;
	case 3: cons[j] = 'T'; break;
	default: cons[j] = '-'; break;
	}
      } else {
	uint32_t k1 = maxIdx;
	uint32_t k2 = sndIdx;
	if (k1 > k2) {
	  k1 = sndIdx;
	  k2 = maxIdx;
	}
	// ACGT-
	if ((k1 == 0) && (k2 == 1)) cons[j] = 'M';
	else if ((k1 == 0) && (k2 == 2)) cons[j] = 'R';
	else if ((k1 == 0) && (k2 == 3)) cons[j] = 'W';
	else if ((k1 == 0) && (k2 == 4)) cons[j] = 'B';
	else if ((k1 == 1) && (k2 == 2)) cons[j] = 'S';
	else if ((k1 == 1) && (k2 == 3)) cons[j] = 'Y';
	else if ((k1 == 1) && (k2 == 4)) cons[j] = 'D';
	else if ((k1 == 2) && (k2 == 3)) cons[j] = 'K';
	else if ((k1 == 2) && (k2 == 4)) cons[j] = 'E';
	else if ((k1 == 3) && (k2 == 4)) cons[j] = 'F';
	else cons[j] = '-';
      }
    }
  }

  
  template<typename TConfig, typename TSplitReadSet>
  inline int
  msaEdlib(TConfig const& c, TSplitReadSet& sps, std::string& cs) {
    float pidth = 0.05; // Expected percent mis-identity threshold
    int32_t minOverlap = 50;   // Min. overlap length

    // Overlap matrices
    std::vector<float> edit(sps.size() * sps.size(), 1);  //1: no overlap (n), <1: overlap, i<j: fwd-fwd (ff), i>j: fwd-rev (fr)
    std::vector<int32_t> loc(sps.size() * sps.size(), 0);  //negative: prefix (p), positive: suffix (s), 0: containment (c)
    
    // Compute overlap matrix
    for(uint32_t i = 0; i < sps.size(); ++i) {
      int32_t iSize = sps[i].size();
      std::string iRev = sps[i];
      reverseComplement(iRev);
      for(uint32_t j = 0; j < sps.size(); ++j) {
	if (i == j) {
	  edit[i * sps.size() + j] = 0;
	  loc[i * sps.size() + j] = 'i';
	}
	int32_t jSize = sps[j].size();
	// Forward - Forward cases
	if (i < j) {
	  // Containment
	  // --------->
	  //    -->
	  if (edit[i * sps.size() + j] == 1) {
	    if (iSize < jSize) {
	      EdlibAlignResult align = edlibAlign(sps[i].c_str(), iSize, sps[j].c_str(), jSize, edlibNewAlignConfig((int) (pidth * iSize), EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
	      if (align.editDistance >= 0) {
		edit[i * sps.size() + j] = (float) align.editDistance / (float) iSize;
		loc[i * sps.size() + j] = 0;
	      }
	      edlibFreeAlignResult(align);
	    } else {
	      EdlibAlignResult align = edlibAlign(sps[j].c_str(), jSize, sps[i].c_str(), iSize, edlibNewAlignConfig((int) (pidth * jSize), EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
	      if (align.editDistance >= 0) {
		edit[i * sps.size() + j] = (float) align.editDistance / (float) jSize;
		loc[i * sps.size() + j] = 0;
	      }
	      edlibFreeAlignResult(align);
	    }
	  }
	  // Suffix
	  //  ----->
	  //     ----->
	  if (edit[i * sps.size() + j] == 1) {
	    for(int32_t osize = std::min(iSize, jSize); osize >= minOverlap; osize -= minOverlap) {
	      std::string suffix = sps[i].substr(iSize - osize);
	      std::string prefix = sps[j].substr(0, osize);
	      EdlibAlignResult align = edlibAlign(suffix.c_str(), osize, prefix.c_str(), osize, edlibNewAlignConfig((int) (pidth * osize), EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
	      if (align.editDistance >= 0) {
		edit[i * sps.size() + j] = (float) align.editDistance / (float) osize;
		loc[i * sps.size() + j] = osize;
		edlibFreeAlignResult(align);
		break;
	      }
	      edlibFreeAlignResult(align);
	    }
	  }
	  // Prefix
	  //    ----->
	  // ----->
	  if (edit[i * sps.size() + j] == 1) {
	    for(int32_t osize = std::min(iSize, jSize); osize >= minOverlap; osize -= minOverlap) {
	      std::string prefix = sps[i].substr(0, osize);
	      std::string suffix = sps[j].substr(jSize - osize);
	      EdlibAlignResult align = edlibAlign(prefix.c_str(), osize, suffix.c_str(), osize, edlibNewAlignConfig((int) (pidth * osize), EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
	      if (align.editDistance >= 0) {
		edit[i * sps.size() + j] = (float) align.editDistance / (float) osize;
		loc[i * sps.size() + j] = -1 * osize;
		edlibFreeAlignResult(align);
		break;
	      }
	      edlibFreeAlignResult(align);
	    }
	  }
	} else { // Forward-Reverse cases
	  // --------->
	  //    <--
	  if (iSize < jSize) {
	    EdlibAlignResult align = edlibAlign(iRev.c_str(), iSize, sps[j].c_str(), jSize, edlibNewAlignConfig((int) (pidth * iSize), EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
	    if (align.editDistance >= 0) {
	      edit[i * sps.size() + j] = (float) align.editDistance / (float) iSize;
	      loc[i * sps.size() + j] = 0;
	    }
	    edlibFreeAlignResult(align);
	  } else {
	    EdlibAlignResult align = edlibAlign(sps[j].c_str(), jSize, iRev.c_str(), iSize, edlibNewAlignConfig((int) (pidth * jSize), EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0));
	    if (align.editDistance >= 0) {
	      edit[i * sps.size() + j] = (float) align.editDistance / (float) jSize;
	      loc[i * sps.size() + j] = 0;
	    }
	    edlibFreeAlignResult(align);
	  }
	  //  <-----
	  //     ----->
	  if (edit[i * sps.size() + j] == 1) {
	    for(int32_t osize = std::min(iSize, jSize); osize >= minOverlap; osize -= minOverlap) {
	      std::string suffix = iRev.substr(iSize - osize);
	      std::string prefix = sps[j].substr(0, osize);
	      EdlibAlignResult align = edlibAlign(suffix.c_str(), osize, prefix.c_str(), osize, edlibNewAlignConfig((int) (pidth * osize), EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
	      if (align.editDistance >= 0) {
		edit[i * sps.size() + j] = (float) align.editDistance / (float) osize;
		loc[i * sps.size() + j] = osize;
		edlibFreeAlignResult(align);
		break;
	      }
	      edlibFreeAlignResult(align);
	    }
	  }
	  //    <-----
	  // ----->
	  if (edit[i * sps.size() + j] == 1) {
	    for(int32_t osize = std::min(iSize, jSize); osize >= minOverlap; osize -= minOverlap) {
	      std::string prefix = iRev.substr(0, osize);
	      std::string suffix = sps[j].substr(jSize - osize);
	      EdlibAlignResult align = edlibAlign(prefix.c_str(), osize, suffix.c_str(), osize, edlibNewAlignConfig((int) (pidth * osize), EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
	      if (align.editDistance >= 0) {
		edit[i * sps.size() + j] = (float) align.editDistance / (float) osize;
		loc[i * sps.size() + j] = -1 * osize;
		edlibFreeAlignResult(align);
		break;
	      }
	      edlibFreeAlignResult(align);
	    }
	  }
	}
      }
    }
    // Debug
    std::cerr << "Containment percent identity" << std::endl;
    for(uint32_t i = 0; i < sps.size(); ++i) {      
      for(uint32_t j = 0; j < sps.size(); ++j) {
	if (j != 0) std::cerr << "\t";
	if (i == j) std::cerr << edit[i * sps.size() + j] << " (ff,i," <<  sps[i].size() << ')' << '\t';
	else {
	  std::string adir = "fr";
	  if (i < j) adir = "ff";
	  if (edit[i * sps.size() + j] == 1) std::cerr << edit[i * sps.size() + j] << " (" << adir << ",n,0)";
	  else {
	    if (loc[i * sps.size() + j] == 0) std::cerr << edit[i * sps.size() + j] << " (" << adir << ",c," << std::min(sps[i].size(), sps[j].size()) << ')';
	    else if (loc[i * sps.size() + j] < 0) std::cerr << edit[i * sps.size() + j] << " (" << adir << ",p," << -1 * loc[i * sps.size() + j] << ')';
	    else std::cerr << edit[i * sps.size() + j] << " (" << adir << ",s," << loc[i * sps.size() + j] << ')';
	  }
	}
      }
      std::cerr << std::endl;
    }

    // Find best sequence to start alignment
    uint32_t bestIdx = 0;
    float bestVal = sps[0].size();
    for(uint32_t i = 0; i < sps.size(); ++i) {
      std::vector<float> dist(sps.size());
      for(uint32_t j = 0; j < sps.size(); ++j) dist[j] = std::min(edit[i * sps.size() + j], edit[j * sps.size() + i]);
      std::sort(dist.begin(), dist.end());
      if (dist[1] >= 1) std::cerr << "Sequence index " << i << " has no alignment!" << std::endl;
      if (dist[sps.size()/2] < bestVal) {
	bestVal = dist[sps.size()/2];
	bestIdx = i;
      }
    }
    
    // Order according to best sequence
    std::vector<uint32_t> selectedIdx;
    if (selectedIdx.empty()) {
      std::vector<std::pair<float, int32_t> > qscores;
      qscores.push_back(std::make_pair(0, bestIdx));
      for(uint32_t j = 0; j < sps.size(); ++j) {
	if (j != bestIdx) qscores.push_back(std::make_pair(std::min(edit[bestIdx * sps.size() + j], edit[j * sps.size() + bestIdx]), j));
      }
      std::sort(qscores.begin(), qscores.end());    
      for(uint32_t i = 0; i < qscores.size(); ++i) {
	// Debug
	std::cerr << qscores[i].first << ',' << qscores[i].second << std::endl;
	selectedIdx.push_back(qscores[i].second);
      }
    }
    exit(-1);
    
    // Extended IUPAC code
    EdlibEqualityPair additionalEqualities[20] = {{'M', 'A'}, {'M', 'C'}, {'R', 'A'}, {'R', 'G'}, {'W', 'A'}, {'W', 'T'}, {'B', 'A'}, {'B', '-'}, {'S', 'C'}, {'S', 'G'}, {'Y', 'C'}, {'Y', 'T'}, {'D', 'C'}, {'D', '-'}, {'K', 'G'}, {'K', 'T'}, {'E', 'G'}, {'E', '-'}, {'F', 'T'}, {'F', '-'}};

    // Incrementally align sequences    
    typedef boost::multi_array<char, 2> TAlign;
    TAlign align;
    align.resize(boost::extents[1][sps[selectedIdx[0]].size()]);
    uint32_t ind = 0;
    for(typename std::string::const_iterator str = sps[selectedIdx[0]].begin(); str != sps[selectedIdx[0]].end(); ++str) align[0][ind++] = *str;
    for(uint32_t i = 1; i < selectedIdx.size(); ++i) {
      // Convert to consensus
      std::string alignStr;
      consensusEdlib(align, alignStr);
      // Debug MSA
      //std::cerr << "Progressive MSA: " << i << '(' << align.shape()[0] << ':' << align.shape()[1] << ')' << std::endl;
      //for(uint32_t i = 0; i<align.shape()[0]; ++i) {
      //for(uint32_t j = 0; j<align.shape()[1]; ++j) std::cerr << align[i][j];
      //std::cerr << std::endl;
      //}
      //std::cerr << "Consensus: " << std::endl;
      //std::cerr << alignStr << std::endl;
      //std::cerr << "ToBeAligned: " << sps[selectedIdx[i]] << std::endl;
      // Compute alignment
      EdlibAlignResult cigar = edlibAlign(sps[selectedIdx[i]].c_str(), sps[selectedIdx[i]].size(), alignStr.c_str(), alignStr.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, additionalEqualities, 20));
      convertAlignment(sps[selectedIdx[i]], align, EDLIB_MODE_HW, cigar);
      edlibFreeAlignResult(cigar);
    }
    
    // Debug MSA
    //std::cerr << "Output MSA " << '(' << align.shape()[0] << ':' << align.shape()[1] << ')' << std::endl;
    //for(uint32_t i = 0; i<align.shape()[0]; ++i) {
    //for(uint32_t j = 0; j<align.shape()[1]; ++j) std::cerr << align[i][j];
    //std::cerr << std::endl;
    //}

    // Consensus
    std::string gapped;
    consensus(c, align, gapped, cs);
    //std::cerr << "Consensus:" << std::endl;
    //std::cerr << gapped << std::endl;

    // OutputAlignment
    msaAlignOut(c, align, gapped);
    
    // Return split-read support
    return align.shape()[0];
  }


}

#endif
