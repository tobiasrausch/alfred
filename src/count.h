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

#ifndef BAMSTATS_H
#define BAMSTATS_H

#include <limits>

#include <boost/icl/split_interval_map.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/progress.hpp>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "util.h"
#include "gtf.h"


namespace bamstats
{

  template<typename TConfig>
  inline int32_t
  countRun(TConfig const& c) {
    // Load bam file
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Parse GTF file
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "GTF file parsing" << std::endl;
    typedef std::vector<IntervalLabel> TChromosomeRegions;
    typedef std::vector<TChromosomeRegions> TGenomicRegions;
    TGenomicRegions gRegions;
    gRegions.resize(hdr->n_targets, TChromosomeRegions());
    typedef std::vector<std::string> TGeneIds;
    TGeneIds geneIds;
    int32_t tf = parseGTF(hdr, c.gtfFile, "exon", "gene_id", gRegions, geneIds);
    if (tf == 0) {
      std::cerr << "Error parsing GTF file!" << std::endl;
      return 1;
    }
    
    // Parse BAM file
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "BAM file parsing" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );

    // Pair qualities and features
    typedef boost::unordered_map<std::size_t, uint8_t> TQualities;
    TQualities qualities;
    typedef boost::unordered_map<std::size_t, int32_t> TFeatures;
    TFeatures features;

    // Feature counter
    typedef std::vector<int32_t> TFeatureCounter;
    TFeatureCounter fc(tf, 0);

    // Iterate chromosomes
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      ++show_progress;
      if (gRegions[refIndex].empty()) continue;

      // Sort by position
      std::sort(gRegions[refIndex].begin(), gRegions[refIndex].end(), SortIntervalStart<IntervalLabel>());

      // Flag ambiguous positions
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet featureBitMap(hdr->target_len[refIndex]);
      TBitSet ambiguousBitMap(hdr->target_len[refIndex]);
      for(uint32_t i = 0; i < gRegions[refIndex].size(); ++i) {
	for(int32_t k = gRegions[refIndex][i].start; k < gRegions[refIndex][i].end; ++k) {
	  if (featureBitMap[k]) ambiguousBitMap[k] = 1;
	  else featureBitMap[k] = 1;
	}
      }

      // Count reads
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP | BAM_FMUNMAP)) continue;
	if (rec->core.tid != rec->core.mtid) continue;
	if (!(rec->core.flag & BAM_FPAIRED) || (rec->core.qual < c.minQual)) continue;

	// Get read sequence
	std::string sequence;
	sequence.resize(rec->core.l_qseq);
	uint8_t* seqptr = bam_get_seq(rec);
	for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

	// Parse CIGAR
	uint32_t* cigar = bam_get_cigar(rec);
	int32_t gp = rec->core.pos; // Genomic position
	int32_t sp = 0; // Sequence position
	int32_t featurepos = -1; // Genomic position where the read intersects a feature
	bool ambiguous = false;
	for (std::size_t i = 0; ((i < rec->core.n_cigar) && (!ambiguous)); ++i) {
	  if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) sp += bam_cigar_oplen(cigar[i]);
	  else if (bam_cigar_op(cigar[i]) == BAM_CINS) sp += bam_cigar_oplen(cigar[i]);
	  else if (bam_cigar_op(cigar[i]) == BAM_CDEL) gp += bam_cigar_oplen(cigar[i]);
	  else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) gp += bam_cigar_oplen(cigar[i]);
	  else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	    //Nop
	  } else if (bam_cigar_op(cigar[i]) == BAM_CMATCH) {
	    for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]); ++k, ++sp, ++gp) {
	      if (ambiguousBitMap[gp]) {
		ambiguous = true;
		break;
	      }
	      if (featureBitMap[gp]) featurepos = gp;
	    }
	  } else {
	    std::cerr << "Unknown Cigar options" << std::endl;
	    return 1;
	  }
	}
	if (ambiguous) continue;  // Ambiguous Read
	if (featurepos == -1) continue; // No feature

	// Find feature
	int32_t featureid = -1;
	TChromosomeRegions::const_iterator vIt = std::lower_bound(gRegions[refIndex].begin(), gRegions[refIndex].end(), IntervalLabel(rec->core.pos), SortIntervalStart<IntervalLabel>());
	TChromosomeRegions::const_iterator prevIt = vIt;
	TChromosomeRegions::const_iterator vItEnd = std::lower_bound(gRegions[refIndex].begin(), gRegions[refIndex].end(), IntervalLabel(featurepos), SortIntervalStart<IntervalLabel>());
	for(;((vIt!=vItEnd) && (featureid == -1)); ++vIt) {
	  if ((vIt->start <= featurepos) && (vIt->end > featurepos)) featureid = vIt->lid;
	}
	if (featureid == -1) {
	  // Backtrack because interval start might be smaller than rec->core.pos
	  while (featureid == -1) {
	    if ((prevIt->start <= featurepos) && (prevIt->end > featurepos)) featureid = vIt->lid;
	    if (prevIt == gRegions[refIndex].begin()) break;
	    else --prevIt;
	  }
	  if (featureid == -1) {
	    std::cerr << "Fatal error: corresponding read feature not found!" << std::endl;
	    return 1;
	  }
	}

	// Check pair
	if (rec->core.pos < rec->core.mpos) {
	  // First read
	  std::size_t hv = hash_pair(rec);
	  qualities[hv] = rec->core.qual;
	  features[hv] = featureid;
	} else {
	  // Second read
	  std::size_t hv=hash_pair_mate(rec);
	  uint8_t pairQuality = std::min((uint8_t) qualities[hv], (uint8_t) rec->core.qual);
	  int32_t featuremate = features[hv];
	  qualities[hv] = 0;
	  features[hv] = -1;

	  // Pair quality
	  if (pairQuality < c.minQual) continue;

	  // Feature disagreement
	  if (featureid != featuremate) continue;

	  // Hurray, we finally have a valid pair
	  ++fc[featureid];
	}
      }
      // Clean-up
      bam_destroy1(rec);
      hts_itr_destroy(iter);
      qualities.clear();
      features.clear();
    }
	  
    // clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);

    // Output count table
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Output count table" << std::endl;
    std::ofstream fcfile(c.outfile.string().c_str());
    fcfile << "gene\t" << c.sampleName << std::endl;
    for(uint32_t idval = 0; idval < geneIds.size(); ++idval) fcfile << geneIds[idval] << "\t" << fc[idval] << std::endl;
    fcfile.close();
    
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    
#ifdef PROFILE
    ProfilerStop();
#endif


    return 0;
  }

}

#endif
