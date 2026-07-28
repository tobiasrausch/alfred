#ifndef ANCHORED_H
#define ANCHORED_H

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <set>
#include <array>
#include <algorithm>
#include <cstdint>
#include <cctype>
#include <cmath>
#include <limits>

#include <boost/multi_array.hpp>
#include <htslib/sam.h>

#include "msa.h"
#include "consedlib.h"

namespace bamstats {

  struct AlignedRead {
    bool anchored;
    bool rev;
    bool distantSA;
    uint16_t mapq;
    int32_t refStart;
    int32_t leftClip;
    int32_t rightClip;
    std::vector<uint32_t> cigar;
    std::string seq;
    std::string name;

    AlignedRead() : anchored(false), rev(false), distantSA(false), mapq(0), refStart(-1), leftClip(0), rightClip(0) {}
    explicit AlignedRead(std::string const& s) : anchored(false), rev(false), distantSA(false), mapq(0), refStart(-1), leftClip(0), rightClip(0), seq(s) {}
  };

  // Candidate allele
  struct AlleleCand {
    std::vector<int32_t> rows;
    int32_t indel;
    bool isRef;
  };
  
  typedef std::map<int32_t, std::map<int32_t, std::string> > TInsClusters;

  // Extract plain sequence
  inline void
  readSequences(std::vector<AlignedRead> const& areads, std::vector<std::string>& rs) {
    rs.clear();
    rs.reserve(areads.size());
    for(std::vector<AlignedRead>::const_iterator it = areads.begin(); it != areads.end(); ++it) rs.push_back(it->seq);
  }

  // Reference span
  inline int32_t
  cigarRefLength(std::vector<uint32_t> const& cigar) {
    int32_t rlen = 0;
    for(uint32_t k = 0; k < cigar.size(); ++k) {
      int op = bam_cigar_op(cigar[k]);
      if ((op == BAM_CMATCH) || (op == BAM_CEQUAL) || (op == BAM_CDIFF) || (op == BAM_CDEL) || (op == BAM_CREF_SKIP)) rlen += bam_cigar_oplen(cigar[k]);
    }
    return rlen;
  }

  // Progressive edlib MSA
  inline void
  insertionMSA(std::vector<std::string> const& seqs, std::vector<std::string>& rows) {
    // Extended IUPAC code
    EdlibEqualityPair additionalEqualities[20] = {{'M', 'A'}, {'M', 'C'}, {'R', 'A'}, {'R', 'G'}, {'W', 'A'}, {'W', 'T'}, {'B', 'A'}, {'B', '-'}, {'S', 'C'}, {'S', 'G'}, {'Y', 'C'}, {'Y', 'T'}, {'D', 'C'}, {'D', '-'}, {'K', 'G'}, {'K', 'T'}, {'E', 'G'}, {'E', '-'}, {'F', 'T'}, {'F', '-'}};

    // Seed with the first sequence
    boost::multi_array<char, 2> align;
    align.resize(boost::extents[1][seqs[0].size()]);
    for(uint32_t j = 0; j < seqs[0].size(); ++j) align[0][j] = seqs[0][j];

    // Incrementally add the remaining sequences
    for(uint32_t i = 1; i < seqs.size(); ++i) {
      std::string cons;
      consensusEdlib(align, cons);
      EdlibAlignResult cigar = edlibAlign(seqs[i].c_str(), seqs[i].size(), cons.c_str(), cons.size(), edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, additionalEqualities, 20));
      convertAlignment(seqs[i], align, EDLIB_MODE_NW, cigar);
      edlibFreeAlignResult(cigar);
    }

    // gapped strings
    rows.assign(align.shape()[0], std::string(align.shape()[1], '-'));
    for(uint32_t i = 0; i < align.shape()[0]; ++i)
      for(uint32_t j = 0; j < align.shape()[1]; ++j) rows[i][j] = align[i][j];
  }

  // Reference-anchored MSA
  inline void
  collectAnchored(std::vector<AlignedRead> const& areads, std::vector<AlignedRead const*>& reads) {
    reads.clear();
    for(std::vector<AlignedRead>::const_iterator it = areads.begin(); it != areads.end(); ++it)
      if ((it->anchored) && (!it->cigar.empty())) reads.push_back(&(*it));
  }

  inline int
  baseToIdx(char c) {
    switch(c) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    }
    return -1;
  }

  // Trim soft-clips
  inline std::string
  trimClipToInsertion(std::string const& clip, std::string const& downstreamRef) {
    if (downstreamRef.empty() || clip.empty()) return clip;
    std::string probe = downstreamRef.substr(0, std::min((std::size_t) 300, downstreamRef.size()));
    EdlibAlignResult r = edlibAlign(probe.c_str(), probe.size(), clip.c_str(), clip.size(), edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    std::string ins = clip;
    if ((r.editDistance >= 0) && (r.editDistance <= (int) (0.15 * probe.size()))) {
      uint32_t s = infixStart(r);
      if (s <= clip.size()) ins = clip.substr(0, s);
    }
    edlibFreeAlignResult(r);
    return ins;
  }


  template<typename TConfig, typename TAlign>
  inline bool
  buildAnchoredMSA(TConfig const& c, std::vector<AlignedRead const*> const& reads, TAlign& align, TInsClusters* clustersOut = NULL, int32_t* refMinOut = NULL) {
    if (reads.empty()) return false;

    // Reference span [refMin, refMax)
    int32_t refMin = std::numeric_limits<int32_t>::max();
    int32_t refMax = 0;
    for(uint32_t i = 0; i < reads.size(); ++i) {
      int32_t rs = reads[i]->refStart;
      int32_t re = rs + cigarRefLength(reads[i]->cigar);
      if (rs < refMin) refMin = rs;
      if (re > refMax) refMax = re;
    }
    int32_t span = refMax - refMin;  // number of reference positions
    if (span <= 0) return false;
    if (refMinOut) *refMinOut = refMin;

    // gather insertions
    struct InsItem {
      int32_t pos;
      int32_t readIdx;
      std::string seq;
    };
    
    struct ClipItem {
      int32_t pos;
      int32_t readIdx;
      bool leading;
      std::string seq;
    };
    
    std::vector<InsItem> items;
    std::vector<ClipItem> clips;
    std::vector<std::array<int32_t, 4> > refCnt(span, std::array<int32_t, 4>{{0, 0, 0, 0}});  // A,C,G,T per ref pos
    const int32_t minClip = c.minClip;
    const int32_t clipCap = c.clipCap;
    for(uint32_t i = 0; i < reads.size(); ++i) {
      int32_t rp = reads[i]->refStart;
      int32_t sp = 0;
      std::string const& seq = reads[i]->seq;
      std::vector<uint32_t> const& cig = reads[i]->cigar;
      for(uint32_t k = 0; k < cig.size(); ++k) {
	int op = bam_cigar_op(cig[k]);
	int ol = bam_cigar_oplen(cig[k]);
	if ((op == BAM_CMATCH) || (op == BAM_CEQUAL) || (op == BAM_CDIFF)) {
	  for(int b = 0; b < ol; ++b) {
	    int32_t idx = rp - refMin;
	    int bi = baseToIdx(seq[sp]);
	    if ((bi >= 0) && (idx >= 0) && (idx < span)) ++refCnt[idx][bi];
	    ++rp; ++sp;
	  }
	}
	else if ((op == BAM_CDEL) || (op == BAM_CREF_SKIP)) {
	  rp += ol;
	}
	else if (op == BAM_CINS) {
	  int32_t idx = rp - refMin;
	  if ((idx >= 0) && (idx <= span)) {
	    InsItem it;
	    it.pos = idx;
	    it.readIdx = (int32_t) i;
	    it.seq = seq.substr(sp, ol);
	    items.push_back(it);
	  }
	  sp += ol;
	} else if (op == BAM_CSOFT_CLIP) {
	  int32_t idx = rp - refMin;
	  bool leading = (sp == 0);
	  if ((ol >= minClip) && (idx >= 0) && (idx <= span)) {
	    int32_t take = std::min(ol, clipCap);
	    std::string cs = leading ? seq.substr(sp + ol - take, take) : seq.substr(sp, take);
	    ClipItem ci; ci.pos = idx; ci.readIdx = (int32_t) i; ci.leading = leading; ci.seq = cs;
	    clips.push_back(ci);
	  }
	  sp += ol;
	}
      }
    }

    // cluster insertions
    const int32_t insWindow = c.insWindow;
    std::sort(items.begin(), items.end(), [](InsItem const& a, InsItem const& b) { return a.pos < b.pos; });
    std::string refCons(span, 'N');
    for(int32_t p = 0; p < span; ++p) {
      int best = -1, bc = 0;
      for(int b = 0; b < 4; ++b) {
	if (refCnt[p][b] > bc) {
	  bc = refCnt[p][b];
	  best = b;
	}
      }
      if (best >= 0) refCons[p] = "ACGT"[best];
    }

    // I-op clusters
    struct Cluster {
      int32_t anchorIdx;
      int32_t posLo;
      int32_t posHi;
      int32_t ltarget;
      std::map<int32_t, std::string> byRead;
    };
    
    std::vector<Cluster> clusters;
    size_t a = 0;
    while (a < items.size()) {
      size_t b = a + 1;
      int32_t clusterMax = items[a].pos;
      while ((b < items.size()) && (items[b].pos - clusterMax <= insWindow)) {
	clusterMax = items[b].pos;
	++b;
      }
      Cluster cl;
      cl.posLo = items[a].pos;
      cl.posHi = clusterMax;
      cl.ltarget = 0;
      // Anchor the block
      int32_t bestLen = -1;
      int32_t bestPos = items[a].pos;
      for(size_t k = a; k < b; ++k) {
	if ((int32_t) items[k].seq.size() > bestLen) {
	  bestLen = (int32_t) items[k].seq.size();
	  bestPos = items[k].pos;
	}
	if (items[k].seq.size() > cl.byRead[items[k].readIdx].size()) cl.byRead[items[k].readIdx] = items[k].seq;
      }
      cl.anchorIdx = bestPos;
      for(std::map<int32_t, std::string>::const_iterator it = cl.byRead.begin(); it != cl.byRead.end(); ++it) {
	if ((int32_t) it->second.size() > cl.ltarget) cl.ltarget = (int32_t) it->second.size();
      }
      clusters.push_back(cl);
      a = b;
    }

    // soft-clips
    const int32_t minInsForClip = c.minInsForClip;
    std::vector<char> clipUsed(clips.size(), 0);
    int32_t nClipAdd = 0;
    for(uint32_t ci = 0; ci < clips.size(); ++ci) {
      for(uint32_t cl = 0; cl < clusters.size(); ++cl) {
	if (clusters[cl].ltarget < minInsForClip) continue;
	if ((clips[ci].pos >= clusters[cl].posLo - insWindow) && (clips[ci].pos <= clusters[cl].posHi + insWindow)) {
	  if (clusters[cl].byRead.find(clips[ci].readIdx) == clusters[cl].byRead.end()) {
	    int32_t L = clusters[cl].ltarget + insWindow;
	    std::string s = clips[ci].seq;
	    if ((int32_t) s.size() > L) s = clips[ci].leading ? s.substr(s.size() - L) : s.substr(0, L);
	    clusters[cl].byRead[clips[ci].readIdx] = s;
	    ++nClipAdd;
	  }
	  clipUsed[ci] = 1;
	  break;
	}
      }
    }

    // new clusters
    const int32_t minClipSupport = c.minClipSupport;
    int32_t nSeed = 0;
    {
      std::vector<ClipItem> un;
      for(uint32_t ci = 0; ci < clips.size(); ++ci) if ((!clipUsed[ci]) && (!clips[ci].leading)) un.push_back(clips[ci]);
      std::sort(un.begin(), un.end(), [](ClipItem const& a, ClipItem const& b) { return a.pos < b.pos; });
      size_t u = 0;
      while (u < un.size()) {
	size_t v = u + 1;
	int32_t cmax = un[u].pos;
	while ((v < un.size()) && (un[v].pos - cmax <= insWindow)) {
	  cmax = un[v].pos;
	  ++v;
	}
	std::map<int32_t, std::string> byRead;
	for(size_t k = u; k < v; ++k) {
	  int32_t idx = un[k].pos;
	  std::string ds = (idx < span) ? refCons.substr(idx, std::min((std::size_t)(span - idx), un[k].seq.size())) : std::string();
	  std::string ins = (ds.size() >= 40) ? trimClipToInsertion(un[k].seq, ds) : un[k].seq;
	  if (((int32_t) ins.size() >= minClip) && (ins.size() > byRead[un[k].readIdx].size())) byRead[un[k].readIdx] = ins;
	}
	if ((int32_t) byRead.size() >= minClipSupport) {
	  Cluster cl;
	  cl.anchorIdx = un[u].pos;
	  cl.posLo = un[u].pos;
	  cl.posHi = cmax;
	  cl.ltarget = 0;
	  cl.byRead = byRead;
	  clusters.push_back(cl);
	  ++nSeed;
	}
	u = v;
      }
    }

    std::map<int32_t, std::map<int32_t, std::string> > merged;
    for(uint32_t cl = 0; cl < clusters.size(); ++cl) {
      for(std::map<int32_t, std::string>::const_iterator it = clusters[cl].byRead.begin(); it != clusters[cl].byRead.end(); ++it) {
	std::string& s = merged[clusters[cl].anchorIdx][it->first];
	if (it->second.size() > s.size()) s = it->second;
      }
    }
    if (clustersOut) *clustersOut = merged;
    std::map<int32_t, std::vector<std::string> > insBlock;
    std::map<int32_t, std::vector<int32_t> > insRows;
    std::vector<int32_t> insW(span + 1, 0);
    int32_t nMulti = 0;
    for(std::map<int32_t, std::map<int32_t, std::string> >::const_iterator mit = merged.begin(); mit != merged.end(); ++mit) {
      std::vector<std::string> seqs;
      std::vector<int32_t> rows;
      for(std::map<int32_t, std::string>::const_iterator it = mit->second.begin(); it != mit->second.end(); ++it) {
	rows.push_back(it->first);
	seqs.push_back(it->second);
      }
      std::vector<std::string> blk;
      insertionMSA(seqs, blk);
      insW[mit->first] = blk.empty() ? 0 : (int32_t) blk[0].size();
      insBlock[mit->first] = blk;
      insRows[mit->first] = rows;
      if (rows.size() >= 2) ++nMulti;
    }

    std::vector<int64_t> colOff(span + 2, 0);
    int64_t total = 0;
    for(int32_t p = 0; p <= span; ++p) {
      colOff[p] = total;
      total += insW[p];
      if (p < span) total += 1;
    }
    colOff[span + 1] = total;
    int64_t ncol = total;

    // alignment matrix
    align.resize(boost::extents[reads.size()][ncol]);
    for(uint32_t i = 0; i < reads.size(); ++i)
      for(int64_t j = 0; j < ncol; ++j) align[i][j] = '-';
    
    // place reference-matching bases
    for(uint32_t i = 0; i < reads.size(); ++i) {
      int32_t rp = reads[i]->refStart;
      int32_t sp = 0;
      std::string const& seq = reads[i]->seq;
      std::vector<uint32_t> const& cig = reads[i]->cigar;
      for(uint32_t k = 0; k < cig.size(); ++k) {
	int op = bam_cigar_op(cig[k]);
	int ol = bam_cigar_oplen(cig[k]);
	if ((op == BAM_CMATCH) || (op == BAM_CEQUAL) || (op == BAM_CDIFF)) {
	  for(int b = 0; b < ol; ++b) {
	    int32_t idx = rp - refMin;
	    int64_t col = colOff[idx] + insW[idx];
	    align[i][col] = seq[sp];
	    ++rp;
	    ++sp;
	  }
	} else if ((op == BAM_CDEL) || (op == BAM_CREF_SKIP)) {
	  rp += ol;
	} else if (op == BAM_CINS) {
	  sp += ol;
	}
	else if (op == BAM_CSOFT_CLIP) {
	  sp += ol;
	}
      }
    }

    // splice the aligned insertion blocks
    for(std::map<int32_t, std::vector<std::string> >::const_iterator it = insBlock.begin(); it != insBlock.end(); ++it) {
      int32_t idx = it->first;
      std::vector<std::string> const& blk = it->second;
      std::vector<int32_t> const& rows = insRows[idx];
      int64_t base = colOff[idx];
      for(uint32_t r = 0; r < rows.size(); ++r) {
	int32_t readIdx = rows[r];
	for(int64_t col = 0; col < (int64_t) blk[r].size(); ++col) align[readIdx][base + col] = blk[r][col];
      }
    }
    return true;
  }

  // Per-read divergence
  template<typename TAlign>
  inline void
  perReadDivergence(TAlign const& align, std::string const& gapped, std::vector<int32_t>& cov, std::vector<double>& div) {
    typedef typename TAlign::index TAIndex;
    uint32_t nrow = align.shape()[0];
    uint32_t ncol = align.shape()[1];
    std::vector<int32_t> mism(nrow, 0);
    cov.assign(nrow, 0);
    for(TAIndex j = 0; j < (TAIndex) ncol; ++j) {
      char g = gapped[j];
      if (g == '-') continue;
      char gu = (char) std::toupper(g);
      for(TAIndex i = 0; i < (TAIndex) nrow; ++i) {
	char a = align[i][j];
	if (a == '-') continue;
	++cov[i];
	if ((char) std::toupper(a) != gu) ++mism[i];
      }
    }
    div.assign(nrow, 0.0);
    for(uint32_t i = 0; i < nrow; ++i) if (cov[i] > 0) div[i] = (double) mism[i] / (double) cov[i];
  }

  // Write MSA
  template<typename TConfig, typename TAlign>
  inline void
  writeAlignment(TConfig const& c, TAlign const& align, std::string const& gapped, std::vector<std::string> const& names, std::string const& path) {
    typedef typename TAlign::index TAIndex;
    boost::iostreams::filtering_ostream rcfile;
    rcfile.push(boost::iostreams::gzip_compressor());
    rcfile.push(boost::iostreams::file_sink(path, std::ios_base::out | std::ios_base::binary));
    if (c.outformat == "h") {
      for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) {
	if (((std::size_t) i < names.size()) && (!names[i].empty())) rcfile << '>' << names[i] << std::endl;
	else rcfile << ">Read" << i << std::endl;
	for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) rcfile << align[i][j];
	rcfile << std::endl;
      }
    } else {
      for(TAIndex j = 0; j < (TAIndex) align.shape()[1]; ++j) {
	for(TAIndex i = 0; i < (TAIndex) align.shape()[0]; ++i) rcfile << align[i][j];
	rcfile << '|' << gapped[j] << std::endl;
      }
    }
    rcfile.pop();
    rcfile.pop();
  }

  // Insert allele number
  inline std::string
  alleleAlignPath(std::string const& base, int32_t k) {
    std::size_t slash = base.find_last_of('/');
    std::size_t start = (slash == std::string::npos) ? 0 : slash + 1;
    std::size_t dot = base.find('.', start);
    if (dot == std::string::npos) return base + "." + std::to_string(k);
    return base.substr(0, dot) + "." + std::to_string(k) + base.substr(dot);
  }

  // Read names
  inline void
  rowNames(std::vector<AlignedRead const*> const& reads, std::vector<std::string>& names) {
    names.clear();
    names.reserve(reads.size());
    for(uint32_t i = 0; i < reads.size(); ++i) names.push_back(reads[i]->name);
  }

  // Drop divergence outliers
  template<typename TConfig>
  inline void
  filterAnchoredReads(TConfig const& c, std::vector<AlignedRead const*>& reads) {
    typedef boost::multi_array<char, 2> TAlign;
    if (reads.empty()) return;
    TAlign align;
    if (!buildAnchoredMSA(c, reads, align)) return;
    std::string gapped, cs;
    consensus(c, align, gapped, cs);
    std::vector<int32_t> cov;
    std::vector<double> div;
    perReadDivergence(align, gapped, cov, div);

    const int32_t minCovStat = 50;
    std::vector<double> ds;
    for(uint32_t i = 0; i < div.size(); ++i) if (cov[i] >= minCovStat) ds.push_back(div[i]);
    double thresh = c.maxdiv;
    if (ds.size() >= 4) {
      std::sort(ds.begin(), ds.end());
      double median = ds[ds.size()/2];
      std::vector<double> ad(ds.size());
      for(uint32_t i = 0; i < ds.size(); ++i) ad[i] = std::fabs(ds[i] - median);
      std::sort(ad.begin(), ad.end());
      double madThresh = median + 5.0 * 1.4826 * ad[ad.size()/2];
      if (madThresh > thresh) thresh = madThresh;
    }
    std::vector<AlignedRead const*> keep;
    for(uint32_t i = 0; i < reads.size(); ++i) {
      bool divOutlier = (cov[i] >= minCovStat) && (div[i] > thresh);
      int32_t readLen = (int32_t) reads[i]->seq.size();
      double clipFrac = (readLen > 0) ? (double) std::max(reads[i]->leftClip, reads[i]->rightClip) / (double) readLen : 0.0;
      bool chimeric = (reads[i]->distantSA) && (clipFrac >= c.maxclip);
      if ((!divOutlier) && (!chimeric)) keep.push_back(reads[i]);
    }
    reads.swap(keep);
  }

  template<typename TConfig>
  inline int
  msaAnchored(TConfig const& c, std::vector<AlignedRead> const& areads, std::string& cs) {
    typedef boost::multi_array<char, 2> TAlign;
    std::vector<AlignedRead const*> reads;
    collectAnchored(areads, reads);
    filterAnchoredReads(c, reads);
    TAlign align;
    if (!buildAnchoredMSA(c, reads, align)) {
      std::cerr << "No reference-anchored reads for consensus!" << std::endl;
      return 0;
    }
    std::string gapped;
    consensus(c, align, gapped, cs);

    std::vector<std::string> names;
    rowNames(reads, names);
    writeAlignment(c, align, gapped, names, c.alignment.string());
    return align.shape()[0];
  }

  // Allele consensus struct
  struct AlleleConsensus {
    std::string name;
    std::string cons;
    int32_t support;
    int32_t indel;
  };

  inline int
  ufFind(std::vector<int>& uf, int x) {
    while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; }
    return x;
  }

  // signed indel length per read
  inline int32_t
  netIndelAt(TInsClusters const& clusters, int32_t refMin, AlignedRead const* rd, int32_t readIdx, int32_t lo, int32_t hi) {
    int32_t insLen = 0;
    int32_t aLo = lo - refMin, aHi = hi - refMin;
    for(TInsClusters::const_iterator it = clusters.begin(); it != clusters.end(); ++it) {
      if ((it->first < aLo) || (it->first > aHi)) continue;
      std::map<int32_t, std::string>::const_iterator jt = it->second.find(readIdx);
      if (jt != it->second.end()) insLen = std::max(insLen, (int32_t) jt->second.size());
    }
    int32_t delLen = 0;
    int32_t rp = rd->refStart;
    for(uint32_t k = 0; k < rd->cigar.size(); ++k) {
      int op = bam_cigar_op(rd->cigar[k]);
      int ol = bam_cigar_oplen(rd->cigar[k]);
      if ((op == BAM_CMATCH) || (op == BAM_CEQUAL) || (op == BAM_CDIFF)) rp += ol;
      else if ((op == BAM_CDEL) || (op == BAM_CREF_SKIP)) {
	if ((rp >= lo) && (rp <= hi)) delLen += ol;
	rp += ol;
      }
    }
    return insLen - delLen;
  }

  // Group reads into alleles
  template<typename TConfig>
  inline void
  alleleGroupsAt(TConfig const& c, std::vector<AlignedRead const*> const& reads, TInsClusters const& clusters, int32_t refMin, int32_t lo, int32_t hi, std::vector<AlleleCand>& cands) {
    cands.clear();
    const int32_t insThresh = c.minInsForClip;
    int32_t mid = (lo + hi) / 2;

    const int32_t W = c.insWindow;

    // Classify spanning reads
    std::vector<int32_t> carrierIdx;
    std::vector<int32_t> carrierNet;
    std::vector<char> carrierOpen;
    std::vector<int32_t> refReads;
    for(uint32_t i = 0; i < reads.size(); ++i) {
      int32_t rs = reads[i]->refStart;
      int32_t re = rs + cigarRefLength(reads[i]->cigar);
      if (!((rs <= mid) && (re > mid))) continue;
      int32_t net = netIndelAt(clusters, refMin, reads[i], (int32_t) i, lo, hi);
      if (std::abs(net) >= insThresh) {
        bool open = false;
        if (net > 0) {
          if ((reads[i]->rightClip >= c.minClip) && (re >= lo - W) && (re <= hi + W)) open = true;
          if ((reads[i]->leftClip  >= c.minClip) && (rs >= lo - W) && (rs <= hi + W)) open = true;
        }
        carrierIdx.push_back((int32_t) i);
        carrierNet.push_back(net);
        carrierOpen.push_back(open ? 1 : 0);
      }
      else refReads.push_back((int32_t) i);
    }

    // Single-linkage cluster by length
    int nn = (int) carrierNet.size();
    std::vector<int> closed;
    std::vector<int> open;
    for(int i = 0; i < nn; ++i) {
      if (carrierOpen[i]) open.push_back(i);
      else closed.push_back(i);
    }
    std::vector<int> uf(nn);
    for(int i = 0; i < nn; ++i) uf[i] = i;
    for(uint32_t x = 0; x < closed.size(); ++x)
      for(uint32_t y = x + 1; y < closed.size(); ++y) {
        int i = closed[x];
	int j = closed[y];
        int32_t mx = std::max(std::abs(carrierNet[i]), std::abs(carrierNet[j]));
        if (std::abs(carrierNet[i] - carrierNet[j]) <= (int32_t)(c.alleleDist * mx)) uf[ufFind(uf, i)] = ufFind(uf, j);
      }
    std::map<int, std::vector<int> > groups;
    for(uint32_t x = 0; x < closed.size(); ++x) groups[ufFind(uf, closed[x])].push_back(closed[x]);

    // Per closed group representative length (median)
    struct Grp { std::vector<int> members; int32_t len; };
    std::vector<Grp> grps;
    for(std::map<int, std::vector<int> >::const_iterator g = groups.begin(); g != groups.end(); ++g) {
      Grp gr; gr.members = g->second;
      std::vector<int32_t> nets;
      for(uint32_t t = 0; t < g->second.size(); ++t) nets.push_back(carrierNet[g->second[t]]);
      std::sort(nets.begin(), nets.end());
      gr.len = nets[nets.size() / 2];
      grps.push_back(gr);
    }

    // Soft-clipped reads
    std::vector<int> leftover;
    for(uint32_t o = 0; o < open.size(); ++o) {
      int oc = open[o];
      int best = -1; int32_t bestLen = std::numeric_limits<int32_t>::max();
      for(uint32_t gi = 0; gi < grps.size(); ++gi) {
        if ((grps[gi].len >= carrierNet[oc]) && (grps[gi].len < bestLen)) {
	  bestLen = grps[gi].len;
	  best = (int) gi;
	}
      }
      if (best >= 0) grps[best].members.push_back(oc);
      else leftover.push_back(oc);
    }
    if (!leftover.empty()) {
      std::vector<int> uf2(nn);
      for(int i = 0; i < nn; ++i) uf2[i] = i;
      for(uint32_t x = 0; x < leftover.size(); ++x)
        for(uint32_t y = x + 1; y < leftover.size(); ++y) {
          int i = leftover[x], j = leftover[y];
          int32_t mx = std::max(std::abs(carrierNet[i]), std::abs(carrierNet[j]));
          if (std::abs(carrierNet[i] - carrierNet[j]) <= (int32_t)(c.alleleDist * mx)) uf2[ufFind(uf2, i)] = ufFind(uf2, j);
        }
      std::map<int, std::vector<int> > lg;
      for(uint32_t x = 0; x < leftover.size(); ++x) lg[ufFind(uf2, leftover[x])].push_back(leftover[x]);
      for(std::map<int, std::vector<int> >::const_iterator g = lg.begin(); g != lg.end(); ++g) {
        Grp gr; gr.members = g->second; gr.len = 0;
        std::vector<int32_t> nets;
        for(uint32_t t = 0; t < g->second.size(); ++t) nets.push_back(carrierNet[g->second[t]]);
        std::sort(nets.begin(), nets.end());
        gr.len = nets[nets.size() / 2];
        grps.push_back(gr);
      }
    }

    // Allele groups, drop singletons 
    for(uint32_t gi = 0; gi < grps.size(); ++gi) {
      if ((int32_t) grps[gi].members.size() < c.minAlleleSupport) continue;
      AlleleCand cd; cd.isRef = false;
      for(uint32_t t = 0; t < grps[gi].members.size(); ++t) cd.rows.push_back(carrierIdx[grps[gi].members[t]]);
      cd.indel = grps[gi].len;
      cands.push_back(cd);
    }
    if ((int32_t) refReads.size() >= c.minAlleleSupport) {
      AlleleCand cd;
      cd.isRef = true;
      cd.indel = 0;
      cd.rows = refReads;
      cands.push_back(cd);
    }
  }

  // Reference position of InDel
  template<typename TConfig>
  inline int32_t
  locusTruePos(TConfig const& c, std::vector<AlignedRead const*> const& reads, TInsClusters const& clusters, int32_t refMin, int32_t lo, int32_t hi) {
    const int32_t insThresh = c.minInsForClip;
    int32_t aLo = lo - refMin;
    int32_t aHi = hi - refMin;
    std::vector<int32_t> positions;
    for(TInsClusters::const_iterator it = clusters.begin(); it != clusters.end(); ++it) {
      if ((it->first < aLo) || (it->first > aHi)) continue;
      for(std::map<int32_t, std::string>::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt) {
        int32_t tlen = (int32_t) jt->second.size();
        if (tlen < insThresh) continue;
        AlignedRead const* rd = reads[jt->first];
        int32_t rp = rd->refStart, best = -1, bestd = std::numeric_limits<int32_t>::max();
        for(uint32_t k = 0; k < rd->cigar.size(); ++k) {
          int op = bam_cigar_op(rd->cigar[k]);
          int ol = bam_cigar_oplen(rd->cigar[k]);
          if ((op == BAM_CMATCH) || (op == BAM_CEQUAL) || (op == BAM_CDIFF)) rp += ol;
          else if ((op == BAM_CDEL) || (op == BAM_CREF_SKIP)) rp += ol;
          else if (op == BAM_CINS) {
            if (ol >= insThresh) {
	      int d = std::abs(ol - tlen);
	      if (d < bestd) {
		bestd = d;
		best = rp;
	      }
	    }
          }
        }
        if (best >= 0) positions.push_back(best);
      }
    }
    for(uint32_t i = 0; i < reads.size(); ++i) {
      int32_t rp = reads[i]->refStart;
      for(uint32_t k = 0; k < reads[i]->cigar.size(); ++k) {
        int op = bam_cigar_op(reads[i]->cigar[k]);
        int ol = bam_cigar_oplen(reads[i]->cigar[k]);
        if ((op == BAM_CMATCH) || (op == BAM_CEQUAL) || (op == BAM_CDIFF)) rp += ol;
        else if ((op == BAM_CDEL) || (op == BAM_CREF_SKIP)) {
	  if ((ol >= insThresh) && (rp >= lo) && (rp <= hi)) positions.push_back(rp);
	  rp += ol;
	}
      }
    }
    if (positions.empty()) return (lo + hi) / 2;
    std::sort(positions.begin(), positions.end());
    return positions[positions.size() / 2];
  }

  // Target position
  inline int32_t
  parseTargetPos(std::string const& position) {
    std::size_t colon = position.find_last_of(':');
    if (colon == std::string::npos) return -1;
    std::string ps = position.substr(colon + 1);
    if (ps.empty() || (ps.find_first_not_of("0123456789") != std::string::npos)) return -1;
    return (int32_t) std::atol(ps.c_str());
  }

  // Split reads into N alleles
  template<typename TConfig>
  inline void
  msaAlleles(TConfig const& c, std::vector<AlignedRead> const& areads, std::vector<AlleleConsensus>& out) {
    typedef boost::multi_array<char, 2> TAlign;
    out.clear();
    std::vector<AlignedRead const*> reads;
    collectAnchored(areads, reads);
    filterAnchoredReads(c, reads);
    if (reads.empty()) return;

    // Anchor MSA
    TAlign align;
    TInsClusters clusters;
    int32_t refMin = 0;
    if (!buildAnchoredMSA(c, reads, align, &clusters, &refMin)) return;

    // Reference span
    int32_t refMax = 0;
    for(uint32_t i = 0; i < reads.size(); ++i) refMax = std::max(refMax, reads[i]->refStart + cigarRefLength(reads[i]->cigar));

    const int32_t insThresh = c.minInsForClip; 

    // Candidates
    std::vector<int32_t> evPos;
    for(TInsClusters::const_iterator it = clusters.begin(); it != clusters.end(); ++it) {
      for(std::map<int32_t, std::string>::const_iterator jt = it->second.begin(); jt != it->second.end(); ++jt) {
        if ((int32_t) jt->second.size() >= insThresh) {
	  evPos.push_back(refMin + it->first); break;
	}
      }
    }
    for(uint32_t i = 0; i < reads.size(); ++i) {
      int32_t rp = reads[i]->refStart;
      for(uint32_t k = 0; k < reads[i]->cigar.size(); ++k) {
        int op = bam_cigar_op(reads[i]->cigar[k]);
        int ol = bam_cigar_oplen(reads[i]->cigar[k]);
        if ((op == BAM_CMATCH) || (op == BAM_CEQUAL) || (op == BAM_CDIFF)) rp += ol;
        else if ((op == BAM_CDEL) || (op == BAM_CREF_SKIP)) {
	  if (ol >= insThresh) evPos.push_back(rp);
	  rp += ol;
	}
      }
    }
    std::sort(evPos.begin(), evPos.end());

    // Candidate loci [lo, hi]
    std::vector<std::pair<int32_t, int32_t> > loci;
    for(size_t a = 0; a < evPos.size(); ) {
      size_t b = a + 1;
      int32_t hi = evPos[a];
      while ((b < evPos.size()) && (evPos[b] - hi <= c.insWindow)) { hi = evPos[b]; ++b; }
      loci.push_back(std::make_pair(evPos[a], hi));
      a = b;
    }

    // Target position
    int32_t targetPos = parseTargetPos(c.position);

    // Find primary InDel locus
    bool primFound = false;
    int32_t primScore = 0;
    int32_t primCarriers = 0;
    int32_t primDist = std::numeric_limits<int32_t>::max();
    int32_t primTruePos = -1;
    std::vector<AlleleCand> primCands;
    for(size_t L = 0; L < loci.size(); ++L) {
      int32_t truePos = locusTruePos(c, reads, clusters, refMin, loci[L].first, loci[L].second);
      if ((targetPos >= 0) && (std::abs(truePos - targetPos) > c.alleleWindow)) continue;
      std::vector<AlleleCand> cands;
      alleleGroupsAt(c, reads, clusters, refMin, loci[L].first, loci[L].second, cands);
      int32_t score = (int32_t) cands.size();
      if (score < 1) continue;
      int32_t carriers = 0;
      for(uint32_t t = 0; t < cands.size(); ++t) if (!cands[t].isRef) carriers += (int32_t) cands[t].rows.size();
      int32_t dist = (targetPos >= 0) ? std::abs(truePos - targetPos) : 0;
      bool better = (score > primScore) || ((score == primScore) && ((dist < primDist) || ((dist == primDist) && (carriers > primCarriers))));
      if (better) {
	primFound = true;
	primScore = score;
	primDist = dist;
	primCarriers = carriers;
	primTruePos = truePos;
	primCands.swap(cands);
      }
    }

    // No indel found
    if (!primFound) {
      std::string gapped, cs;
      consensus(c, align, gapped, cs);
      AlleleConsensus a;
      a.name = "Allele0_ref_n" + std::to_string(reads.size());
      a.cons = cs;
      a.support = (int32_t) reads.size();
      a.indel = 0;
      out.push_back(a);
      return;
    }

    // Sort by support
    std::vector<AlleleCand>& cands = primCands;
    std::sort(cands.begin(), cands.end(), [](AlleleCand const& a, AlleleCand const& b) { return a.rows.size() > b.rows.size(); });
    if ((c.nAlleles > 0) && ((int32_t) cands.size() > c.nAlleles)) cands.resize(c.nAlleles);

    // Per-allele consensus
    for(uint32_t k = 0; k < cands.size(); ++k) {
      std::vector<AlignedRead const*> sub;
      for(uint32_t t = 0; t < cands[k].rows.size(); ++t) sub.push_back(reads[cands[k].rows[t]]);
      TAlign aal;
      if (!buildAnchoredMSA(c, sub, aal)) continue;
      std::string g2, cs2;
      consensus(c, aal, g2, cs2);
      AlleleConsensus ac;
      ac.support = (int32_t) cands[k].rows.size();
      ac.indel = cands[k].indel;
      ac.cons = cs2;
      std::string tag;
      if (cands[k].isRef) tag = "ref";
      else if (cands[k].indel >= 0) tag = "ins" + std::to_string(cands[k].indel);
      else tag = "del" + std::to_string(-cands[k].indel);
      ac.name = "Allele" + std::to_string(k) + "_" + tag + "_n" + std::to_string(ac.support);
      out.push_back(ac);

      // Per-allele alignment file
      std::vector<std::string> names;
      rowNames(sub, names);
      writeAlignment(c, aal, g2, names, alleleAlignPath(c.alignment.string(), (int32_t) k));

      // Read membership
      std::cout << ac.name << " reads:";
      for(uint32_t t = 0; t < names.size(); ++t) std::cout << ' ' << names[t];
      std::cout << std::endl;
    }
    std::cout << "Alleles: " << out.size() << " (primary indel near " << primTruePos << ", " << primCarriers << " carriers)" << std::endl;
  }

}

#endif
