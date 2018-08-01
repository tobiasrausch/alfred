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

#ifndef BAMSTATS_H
#define BAMSTATS_H

#include <limits>

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

#include "tenX.h"
#include "util.h"
#include "json.h"
#include "tsv.h"
#include "qcstruct.h"

namespace bamstats
{

  template<typename TChromosomeRegions, typename TBpCoverage, typename TBedCounts>
  inline void
  _summarizeBedCoverage(TChromosomeRegions const& chrRegions, TBpCoverage const& cov, int32_t refIndex, std::string const& rg, TBedCounts& be) {
    typename BedCounts::TRgBpMap::iterator itChr = be.gCov[refIndex].find(rg);
    typename BedCounts::TOnTargetMap::iterator itOT = be.onTarget.find(rg);
    for(int32_t s = 0; s < (int32_t) itOT->second.size(); ++s) {
      // Avoid over-counting
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet used(cov.size());
      for(uint32_t i = 0; i<chrRegions.size(); ++i) {
	int64_t avgCov = 0;
	int32_t rStart = std::max(0, chrRegions[i].start - s * be.stepsize);
	int32_t rEnd = std::min((int32_t) cov.size(), chrRegions[i].end + s * be.stepsize);
	for(int32_t k = rStart; k < rEnd; ++k) {
	  if (!used[k]) {
	    avgCov += cov[k];
	    used[k] = 1;
	  }
	}
	itOT->second[s] += avgCov;
	if (s == 0) {
	  if (chrRegions[i].start < chrRegions[i].end) itChr->second[i] = (typename BedCounts::TAvgCov) ( (double) avgCov / (double) (chrRegions[i].end - chrRegions[i].start));
	  else itChr->second[i] = (typename BedCounts::TAvgCov) (0);
	}
      }
    }
  }

  template<typename TConfig>
  inline int32_t
  bamStatsRun(TConfig& c) {
    // Load bam file
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // One list of regions for every chromosome
    typedef std::vector<Interval> TChromosomeRegions;
    typedef std::vector<TChromosomeRegions> TGenomicRegions;
    TGenomicRegions gRegions;
    gRegions.resize(hdr->n_targets, TChromosomeRegions());

    // Parse regions from BED file or create one region per chromosome
    int32_t totalBedSize = 0;
    if (c.hasRegionFile) {
      if (is_gz(c.regionFile)) {
	std::ifstream file(c.regionFile.string().c_str(), std::ios_base::in | std::ios_base::binary);
	boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
	dataIn.push(boost::iostreams::gzip_decompressor());
	dataIn.push(file);
	std::istream instream(&dataIn);
	std::string line;
	while(std::getline(instream, line)) {
	  typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	  boost::char_separator<char> sep(" \t,;");
	  Tokenizer tokens(line, sep);
	  Tokenizer::iterator tokIter = tokens.begin();
	  std::string chrName = *tokIter++;
	  // Map chromosome names to the bam header chromosome IDs
	  int32_t chrid = bam_name2id(hdr, chrName.c_str());
	  // Valid ID?
	  if (chrid >= 0) {
	    if (tokIter!=tokens.end()) {
	      int32_t start = boost::lexical_cast<int32_t>(*tokIter++);
	      int32_t end = boost::lexical_cast<int32_t>(*tokIter++);
	      gRegions[chrid].push_back(Interval(start, end));
	    }
	  }
	}
	dataIn.pop();
      } else {
	// Parse regions from BED file
	std::ifstream bedFile(c.regionFile.string().c_str(), std::ifstream::in);
	if (bedFile.is_open()) {
	  while (bedFile.good()) {
	    std::string line;
	    getline(bedFile, line);
	    typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	    boost::char_separator<char> sep(" \t,;");
	    Tokenizer tokens(line, sep);
	    Tokenizer::iterator tokIter = tokens.begin();
	    std::string chrName = *tokIter++;
	    // Map chromosome names to the bam header chromosome IDs
	    int32_t chrid = bam_name2id(hdr, chrName.c_str());
	    // Valid ID?
	    if (chrid >= 0) {
	      if (tokIter!=tokens.end()) {
		int32_t start = boost::lexical_cast<int32_t>(*tokIter++);
		int32_t end = boost::lexical_cast<int32_t>(*tokIter++);
		gRegions[chrid].push_back(Interval(start, end));
	      }
	    }
	  }
	  bedFile.close();
	}
      }

      // Get total bed size
      for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
	typedef boost::dynamic_bitset<> TBitSet;
	TBitSet bedcovered(hdr->target_len[refIndex]);
	for(uint32_t i = 0; i < gRegions[refIndex].size(); ++i)
	  for(int32_t k = gRegions[refIndex][i].start; (k < gRegions[refIndex][i].end) && (k < (int32_t) hdr->target_len[refIndex]); ++k) bedcovered[k] = 1;
	totalBedSize += bedcovered.count();
      }
    }
    
    // Debug code
    //uint32_t rIndex = 0;
    //for(TGenomicRegions::const_iterator itG = gRegions.begin(); itG != gRegions.end(); ++itG, ++rIndex) {
    //for(TChromosomeRegions::const_iterator itC = itG->begin(); itC != itG->end(); ++itC) {
    //std::cout << rIndex << ',' << hdr->target_name[rIndex] << ',' << itC->start << ',' << itC->end << std::endl;
    //}
    //}

    // BED file statistics
    BedCounts be(hdr->n_targets, 25, 20);
    
    // Read group statistics
    typedef std::set<std::string> TRgSet;
    TRgSet rgs;
    getRGs(std::string(hdr->text), rgs);
    if (c.ignoreRG) rgs.insert("DefaultLib");
    if ((c.singleRG) && (rgs.find(c.rgname) == rgs.end())) {
	std::cerr << "Read group is not present in BAM file: " << c.rgname << std::endl;
	return 1;
    }
    typedef boost::unordered_map<std::string, ReadGroupStats> TRGMap;
    TRGMap rgMap;
    for(typename TRgSet::const_iterator itRg = rgs.begin(); itRg != rgs.end(); ++itRg) {
      if (((c.ignoreRG) && (*itRg == "DefaultLib")) || ((c.singleRG) && (*itRg == c.rgname)) || ((!c.ignoreRG) && (!c.singleRG))) {
	rgMap.insert(std::make_pair(*itRg, ReadGroupStats(hdr->n_targets)));
	for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
	  typename BedCounts::TRgBpMap::iterator itChr = be.gCov[refIndex].insert(std::make_pair(*itRg, typename BedCounts::TBpCov())).first;
	  itChr->second.resize(gRegions[refIndex].size());
	  typename BedCounts::TOnTargetMap::iterator itOT = be.onTarget.insert(std::make_pair(*itRg, typename BedCounts::TOnTargetBp())).first;
	  itOT->second.resize(be.onTSize, 0);
	}
      }
    }

    // Parse reference and BAM file
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "BAM file parsing" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );

    // GC- and N-content
    typedef boost::dynamic_bitset<> TBitSet;
    TBitSet nrun;
    TBitSet gcref;
    uint64_t referencebp = 0;
    uint64_t ncount = 0;
    int32_t gcRunnerIdx = 0;
    int32_t gcRunnerCount = 0;
    std::vector<ChrGC> chrGC(hdr->n_targets, ChrGC());
    ReadCounts::TGCContent refGcContent(101, 0);

    // Find N90 chromosome length
    uint32_t minChrLen = 0;
    {
      std::vector<uint32_t> chrlen(hdr->n_targets, 0);
      uint64_t genomelen = 0;
      for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
	chrlen[refIndex] = hdr->target_len[refIndex];
	genomelen += hdr->target_len[refIndex];
      }
      std::sort(chrlen.begin(), chrlen.end(), std::greater<uint32_t>());
      uint64_t cumsum = 0;
      for(uint32_t i = 0; i < chrlen.size(); ++i) {
	cumsum += chrlen[i];
	if (cumsum > genomelen * 0.9) {
	  minChrLen = chrlen[i];
	  break;
	}
      }
    }
    
    // Parse genome
    int32_t refIndex = -1;
    char* seq = NULL;
    faidx_t* fai = fai_load(c.genome.string().c_str());
    bam1_t* rec = bam_init1();
    while (sam_read1(samfile, hdr, rec) >= 0) {
      // New chromosome?
      if ((!(rec->core.flag & BAM_FUNMAP)) && (rec->core.tid != refIndex)) {
	++show_progress;
	
	// Summarize bp-level coverage
	if (refIndex != -1) {
	  for(typename TRGMap::iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
	    if ((c.hasRegionFile) && (!gRegions[refIndex].empty())) _summarizeBedCoverage(gRegions[refIndex], itRg->second.bc.cov, refIndex, itRg->first, be);
	    for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	      if (itRg->second.bc.cov[i] >= 1) {
		++itRg->second.bc.nd;
		if (itRg->second.bc.cov[i] == 1) ++itRg->second.bc.n1;
		if (itRg->second.bc.cov[i] == 2) ++itRg->second.bc.n2;
	      }
	      if (!nrun[i]) ++itRg->second.bc.bpWithCoverage[itRg->second.bc.cov[i]];
	    }
	    itRg->second.bc.cov.clear();
	  }
	  if (seq != NULL) free(seq);
	}
	refIndex = rec->core.tid;
	
	// Load chromosome
	int32_t seqlen = -1;
	std::string tname(hdr->target_name[refIndex]);
	seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);

	// Set N-mask
	nrun.clear();
	nrun.resize(hdr->target_len[refIndex], false);
	gcref.clear();
	gcref.resize(hdr->target_len[refIndex], false);
	referencebp += hdr->target_len[refIndex];
	for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	  if ((seq[i] == 'c') || (seq[i] == 'C') || (seq[i] == 'g') || (seq[i] == 'G')) gcref[i] = 1;
	  if ((seq[i] == 'n') || (seq[i] == 'N')) {
	    nrun[i] = 1;
	    ++ncount;
	  }
	}
	// Reference GC
	chrGC[refIndex].ncount = nrun.count();
	chrGC[refIndex].gccount = gcref.count();
	if ((hdr->target_len[refIndex] > 100) && (hdr->target_len[refIndex] >= minChrLen)) {
	  uint32_t nsum = 0;
	  uint32_t gcsum = 0;
	  for(uint32_t pos = 0; pos < hdr->target_len[refIndex] - 100; ++pos) {
	    if (pos == 0) {
	      for(uint32_t i = 0; i < 100; ++i) {
		nsum += nrun[i];
		gcsum += gcref[i];
	      }
	    } else {
	      nsum -= nrun[pos - 1];
	      gcsum -= gcref[pos - 1];
	      nsum += nrun[pos + 99];
	      gcsum += gcref[pos + 99];
	    }
	    if (!nsum) ++refGcContent[gcsum];
	  }
	}
	
	// Resize coverage vectors
	for(typename TRGMap::iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) itRg->second.bc.cov.resize(hdr->target_len[refIndex], 0);
      }
      
      // Get the library information
      std::string rG = "DefaultLib";
      if (!c.ignoreRG) {
	uint8_t *rgptr = bam_aux_get(rec, "RG");
	if (rgptr) {
	  char* rg = (char*) (rgptr + 1);
	  rG = std::string(rg);
	}
	if ((c.singleRG) && (rG != c.rgname)) continue;
      }
      typename TRGMap::iterator itRg = rgMap.find(rG);
      if (itRg == rgMap.end()) {
	std::cerr << "Missing read group: " << rG << std::endl;
	return 1;
      }

      // Alignments behind the reference end
      if ((!(rec->core.flag & BAM_FUNMAP)) && (((rec->core.pos >= (int32_t) hdr->target_len[refIndex]) || (lastAlignedPosition(rec) > hdr->target_len[refIndex])))) {
	std::cerr << "Alignment is past the reference end: " << hdr->target_name[refIndex] << ':' << rec->core.pos << std::endl;
	continue;
      }

      // Paired counts
      if (rec->core.flag & BAM_FPAIRED) {
	++itRg->second.pc.paired;
	if (!((rec->core.flag & BAM_FUNMAP) || (rec->core.flag & BAM_FMUNMAP))) {
	  ++itRg->second.pc.mapped;
	  if (rec->core.tid == rec->core.mtid) {
	    ++itRg->second.pc.mappedSameChr;
	    if (rec->core.flag & BAM_FPROPER_PAIR) ++itRg->second.pc.mappedProper;
	  }
	  if (rec->core.pos > rec->core.mpos) {
	    ++itRg->second.pc.totalISizeCount;
	    int32_t outerISize = rec->core.pos - rec->core.mpos + alignmentLength(rec);
	    switch(layout(rec)) {
	    case 0:
	      ++itRg->second.pc.orient[0];
	      if (outerISize < itRg->second.pc.maxInsertSize) ++itRg->second.pc.fPlus[outerISize];
	      else ++itRg->second.pc.fPlus[itRg->second.pc.maxInsertSize];
	      break;
	    case 1:
	      if (outerISize < itRg->second.pc.maxInsertSize) ++itRg->second.pc.fMinus[outerISize];
	      else ++itRg->second.pc.fMinus[itRg->second.pc.maxInsertSize];
	      break;
	    case 2:
	      ++itRg->second.pc.orient[2];
	      if (outerISize < itRg->second.pc.maxInsertSize) ++itRg->second.pc.rPlus[outerISize];
	      else ++itRg->second.pc.rPlus[itRg->second.pc.maxInsertSize];
	      break;
	    case 3:
	      ++itRg->second.pc.orient[3];
	      if (outerISize < itRg->second.pc.maxInsertSize) ++itRg->second.pc.rMinus[outerISize];
	      else ++itRg->second.pc.rMinus[itRg->second.pc.maxInsertSize];
	      break;
	    default:
	      break;
	    }
	  }
	}
      }

      // Read counts
      if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) {
	if (rec->core.flag & BAM_FSECONDARY) ++itRg->second.rc.secondary;
	if (rec->core.flag & BAM_FQCFAIL) ++itRg->second.rc.qcfail;
	if (rec->core.flag & BAM_FDUP) ++itRg->second.rc.dup;
	if (rec->core.flag & BAM_FSUPPLEMENTARY) ++itRg->second.rc.supplementary;
	if (rec->core.flag & BAM_FUNMAP) {
	  ++itRg->second.rc.unmap;
	} else {
	  ++itRg->second.rc.mappedchr[refIndex];
	}
	continue;
      }
      ++itRg->second.qc.qcount[(int32_t) rec->core.qual];
      ++itRg->second.rc.mappedchr[refIndex];
      if (rec->core.flag & BAM_FREAD2) ++itRg->second.rc.mapped2;
      else ++itRg->second.rc.mapped1;
      if (rec->core.flag & BAM_FREVERSE) ++itRg->second.rc.reverse;
      else ++itRg->second.rc.forward;
      if (rec->core.l_qseq < itRg->second.rc.maxReadLength) ++itRg->second.rc.lRc[rec->core.l_qseq];
      else ++itRg->second.rc.lRc[itRg->second.rc.maxReadLength];

      // Fetch molecule identifier
      uint8_t* miptr = bam_aux_get(rec, "MI");
      if (miptr) {
	c.isMitagged = true;
	++itRg->second.rc.mitagged;
	int32_t mitag = bam_aux2i(miptr);
	if ((mitag>=0) && (mitag < itRg->second.rc.maxUMI)) {
	  // Lazy resize
	  if (itRg->second.rc.umi.empty()) itRg->second.rc.umi.resize(itRg->second.rc.maxUMI, false);
	  itRg->second.rc.umi[mitag] = true;
	}
      }
      
      // Fetch haplotype tag
      uint8_t* hpptr = bam_aux_get(rec, "HP");
      if (hpptr) {
	c.isHaplotagged = true;
	++itRg->second.rc.haplotagged;

	// If no phased block assume all in one phased block
	int32_t psId = 0;
	uint8_t* psptr = bam_aux_get(rec, "PS");
	if (psptr) psId = bam_aux2i(psptr);
	if ((int32_t) itRg->second.rc.brange.size() <= refIndex) itRg->second.rc.brange.resize(refIndex + 1, ReadCounts::TBlockRange());
	if (itRg->second.rc.brange[refIndex].find(psId) == itRg->second.rc.brange[refIndex].end()) {
	  itRg->second.rc.brange[refIndex].insert(std::make_pair(psId, std::make_pair(rec->core.pos, lastAlignedPosition(rec))));
	} else {
	  itRg->second.rc.brange[refIndex][psId].first = std::min(rec->core.pos, itRg->second.rc.brange[refIndex][psId].first);
	  itRg->second.rc.brange[refIndex][psId].second = std::max((int32_t) lastAlignedPosition(rec), itRg->second.rc.brange[refIndex][psId].second);
	}
      }
      
      // Get the read sequence
      typedef std::vector<uint8_t> TQuality;
      TQuality quality;
      quality.resize(rec->core.l_qseq);
      std::string sequence;
      sequence.resize(rec->core.l_qseq);
      uint8_t* seqptr = bam_get_seq(rec);
      uint8_t* qualptr = bam_get_qual(rec);
      for (int32_t i = 0; i < rec->core.l_qseq; ++i) {
	quality[i] = qualptr[i];
	sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	//char c = 33 + quality[i];
	int32_t relpos = i;
	if (rec->core.flag & BAM_FREVERSE) {
	  relpos = rec->core.l_qseq - i - 1;
	  if (relpos < itRg->second.rc.maxReadLength) {
	    itRg->second.rc.bqCount[relpos] += (uint64_t) quality[i];
	    if ((sequence[i] == 'N') || (sequence[i] == 'n')) ++itRg->second.rc.nCount[relpos];
	    else if ((sequence[i] == 'A') || (sequence[i] == 'a')) ++itRg->second.rc.tCount[relpos];
	    else if ((sequence[i] == 'C') || (sequence[i] == 'c')) ++itRg->second.rc.gCount[relpos];
	    else if ((sequence[i] == 'G') || (sequence[i] == 'g')) ++itRg->second.rc.cCount[relpos];
	    else if ((sequence[i] == 'T') || (sequence[i] == 't')) ++itRg->second.rc.aCount[relpos];
	  } else {
	    itRg->second.rc.bqCount[itRg->second.rc.maxReadLength] += (uint64_t) quality[i];
	    if ((sequence[i] == 'N') || (sequence[i] == 'n')) ++itRg->second.rc.nCount[itRg->second.rc.maxReadLength];
	    else if ((sequence[i] == 'A') || (sequence[i] == 'a')) ++itRg->second.rc.tCount[itRg->second.rc.maxReadLength];
	    else if ((sequence[i] == 'C') || (sequence[i] == 'c')) ++itRg->second.rc.gCount[itRg->second.rc.maxReadLength];
	    else if ((sequence[i] == 'G') || (sequence[i] == 'g')) ++itRg->second.rc.cCount[itRg->second.rc.maxReadLength];
	    else if ((sequence[i] == 'T') || (sequence[i] == 't')) ++itRg->second.rc.aCount[itRg->second.rc.maxReadLength];
	  }
	} else {
	  if (relpos < itRg->second.rc.maxReadLength) {
	    itRg->second.rc.bqCount[relpos] += (uint64_t) quality[i];
	    if ((sequence[i] == 'N') || (sequence[i] == 'n')) ++itRg->second.rc.nCount[relpos];
	    else if ((sequence[i] == 'A') || (sequence[i] == 'a')) ++itRg->second.rc.aCount[relpos];
	    else if ((sequence[i] == 'C') || (sequence[i] == 'c')) ++itRg->second.rc.cCount[relpos];
	    else if ((sequence[i] == 'G') || (sequence[i] == 'g')) ++itRg->second.rc.gCount[relpos];
	    else if ((sequence[i] == 'T') || (sequence[i] == 't')) ++itRg->second.rc.tCount[relpos];
	  } else {
	    itRg->second.rc.bqCount[itRg->second.rc.maxReadLength] += (uint64_t) quality[i];
	    if ((sequence[i] == 'N') || (sequence[i] == 'n')) ++itRg->second.rc.nCount[itRg->second.rc.maxReadLength];
	    else if ((sequence[i] == 'A') || (sequence[i] == 'a')) ++itRg->second.rc.aCount[itRg->second.rc.maxReadLength];
	    else if ((sequence[i] == 'C') || (sequence[i] == 'c')) ++itRg->second.rc.cCount[itRg->second.rc.maxReadLength];
	    else if ((sequence[i] == 'G') || (sequence[i] == 'g')) ++itRg->second.rc.gCount[itRg->second.rc.maxReadLength];
	    else if ((sequence[i] == 'T') || (sequence[i] == 't')) ++itRg->second.rc.tCount[itRg->second.rc.maxReadLength];
	  }
	}
      }

      // Sequence GC content
      if (sequence.size() > 100) {
	uint32_t offset = (sequence.size() - 100) / 2;
	int32_t gccont = 0;
	for(uint32_t i = offset; i < (offset + 100); ++i) {
	  if ((sequence[i] == 'c') || (sequence[i] == 'C') || (sequence[i] == 'g') || (sequence[i] == 'G')) ++gccont;
	}
	++itRg->second.rc.gcContent[gccont];
      } else {
	// Most likely some tag-counting application, ignore the first 20bp
	if (!(rec->core.flag & BAM_FREVERSE)) {
	  for(uint32_t i = 20; i < sequence.size(); ++i, ++gcRunnerIdx) {
	    if ((sequence[i] == 'c') || (sequence[i] == 'C') || (sequence[i] == 'g') || (sequence[i] == 'G')) ++gcRunnerCount;
	    if (gcRunnerIdx == 100) {
	      ++itRg->second.rc.gcContent[gcRunnerCount];
	      gcRunnerIdx = 0;
	      gcRunnerCount = 0;
	    }
	  }
	}
      }
    
      // Get the reference slice
      std::string refslice = boost::to_upper_copy(std::string(seq + rec->core.pos, seq + lastAlignedPosition(rec)));
      
      // Debug 
      //std::cout << matchCount << ',' << mismatchCount << ',' << delCount << ',' << insCount << ',' << softClipCount << ',' << hardClipCount << std::endl;
      //std::cout << refslice << std::endl;
      //std::cout << sequence << std::endl;

      uint32_t rp = 0; // reference pointer
      uint32_t sp = 0; // sequence pointer
      
      // Parse the CIGAR
      uint32_t* cigar = bam_get_cigar(rec);
      bool spliced = false;
      for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	  // match or mismatch
	  for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]);++k) {
	    if (sequence[sp] == refslice[rp]) ++itRg->second.bc.matchCount;
	    else ++itRg->second.bc.mismatchCount;
	    // Count bp-level coverage
	    if (itRg->second.bc.cov[rec->core.pos + rp] < itRg->second.bc.maxCoverage) ++itRg->second.bc.cov[rec->core.pos + rp];
	    ++sp;
	    ++rp;
	  }
	} else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	  ++itRg->second.bc.delCount;
	  ++itRg->second.bc.delHomACGTN[homopolymerContext(sequence, sp, 3)];
	  if (bam_cigar_oplen(cigar[i]) < itRg->second.bc.maxIndelSize) ++itRg->second.bc.delSize[bam_cigar_oplen(cigar[i])];
	  else ++itRg->second.bc.delSize[itRg->second.bc.maxIndelSize];
	  rp += bam_cigar_oplen(cigar[i]);
	} else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	  ++itRg->second.bc.insCount;
	  ++itRg->second.bc.insHomACGTN[homopolymerContext(sequence, sp, 3)];
	  if (bam_cigar_oplen(cigar[i]) < itRg->second.bc.maxIndelSize) ++itRg->second.bc.insSize[bam_cigar_oplen(cigar[i])];
	  else ++itRg->second.bc.insSize[itRg->second.bc.maxIndelSize];
	  sp += bam_cigar_oplen(cigar[i]);
	} else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	  ++itRg->second.bc.softClipCount;
	  sp += bam_cigar_oplen(cigar[i]);
	} else if(bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	  ++itRg->second.bc.hardClipCount;
	} else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	  if (!spliced) {
	    ++itRg->second.rc.spliced;
	    spliced = true;
	  }
	  rp += bam_cigar_oplen(cigar[i]);
	} else {
	  std::cerr << "Unknown Cigar options" << std::endl;
	  return 1;
	}
      }
    }
    // Summarize bp-level coverage
    if (refIndex != -1) {
      for(typename TRGMap::iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
	if ((c.hasRegionFile) && (!gRegions[refIndex].empty())) _summarizeBedCoverage(gRegions[refIndex], itRg->second.bc.cov, refIndex, itRg->first, be);
	for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) {
	  if (itRg->second.bc.cov[i] >= 1) {
	    ++itRg->second.bc.nd;
	    if (itRg->second.bc.cov[i] == 1) ++itRg->second.bc.n1;
	    if (itRg->second.bc.cov[i] == 2) ++itRg->second.bc.n2;
	  }
	  if (!nrun[i]) ++itRg->second.bc.bpWithCoverage[itRg->second.bc.cov[i]];
	}
	itRg->second.bc.cov.clear();
      }
      if (seq != NULL) free(seq);
    }

    // Output
    qcTsvOut(c, hdr, rgMap, be, gRegions, refGcContent, chrGC, referencebp, ncount, totalBedSize);
    
    // clean-up
    bam_destroy1(rec);
    fai_destroy(fai);
    bam_hdr_destroy(hdr);
    sam_close(samfile);
    
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    
#ifdef PROFILE
    ProfilerStop();
#endif


    return 0;
  }

}

#endif
