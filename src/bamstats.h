#ifndef BAMSTATS_H
#define BAMSTATS_H

#include <limits>

#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/progress.hpp>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "util.h"

namespace bamstats
{

  struct Interval {
    int32_t start;
    int32_t end;
    
    Interval(int32_t s, int32_t e) : start(s), end(e) {}
  };

  struct ReadGroupStats {
    typedef std::vector<uint64_t> TCoverageBp;
    typedef uint8_t TBpCovInt;
    typedef std::vector<TBpCovInt> TBpCoverage;
    
    uint32_t maxBpCov;
    uint64_t matchCount;
    uint64_t mismatchCount;
    uint64_t delCount;
    uint64_t insCount;
    uint64_t softClipCount;
    uint64_t hardClipCount;
    TCoverageBp bpWithCoverage;
    TBpCoverage cov;
    
    ReadGroupStats() : maxBpCov(std::numeric_limits<TBpCovInt>::max()), matchCount(0), mismatchCount(0), delCount(0), insCount(0), softClipCount(0), hardClipCount(0) {
      bpWithCoverage.resize(maxBpCov + 1, 0);
      cov.clear();
    }
  };

  
  template<typename TConfig>
  inline int32_t
  bamStatsRun(TConfig const& c) {
    // Load bam file
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // One list of regions for every chromosome
    typedef std::vector<Interval> TChromosomeRegions;
    typedef std::vector<TChromosomeRegions> TGenomicRegions;
    TGenomicRegions gRegions;
    gRegions.resize(hdr->n_targets, TChromosomeRegions());

    // Parse regions from BED file or create one region per chromosome
    if (c.hasRegionFile) {
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
      }
    } else {
      // Make one region for every chromosome
      for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) gRegions[refIndex].push_back(Interval(0, hdr->target_len[refIndex]));
    }
    
    // Debug code
    //uint32_t rIndex = 0;
    //for(TGenomicRegions::const_iterator itG = gRegions.begin(); itG != gRegions.end(); ++itG, ++rIndex) {
    //for(TChromosomeRegions::const_iterator itC = itG->begin(); itC != itG->end(); ++itC) {
    //std::cout << rIndex << ',' << hdr->target_name[rIndex] << ',' << itC->start << ',' << itC->end << std::endl;
    //}
    //}

    // Read group statistics
    typedef std::set<std::string> TRgSet;
    TRgSet rgs;
    getRGs(std::string(hdr->text), rgs);
    typedef boost::unordered_map<std::string, ReadGroupStats> TRGMap;
    TRGMap rgMap;
    for(typename TRgSet::const_iterator itRg = rgs.begin(); itRg != rgs.end(); ++itRg) rgMap.insert(std::make_pair(*itRg, ReadGroupStats()));

    // Total reference base pairs
    uint64_t referencebp = 0;


    // Parse reference and BAM file
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "BAM file parsing" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );

    // Parse genome
    faidx_t* fai = fai_load(c.genome.string().c_str());
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      ++show_progress;

      // Load chromosome
      char* seq = NULL;
      int32_t seqlen = -1;
      std::string tname(hdr->target_name[refIndex]);
      seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);
      referencebp += hdr->target_len[refIndex];

      // Resize coverage vectors
      for(typename TRGMap::iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) itRg->second.cov.resize(hdr->target_len[refIndex], 0);

      // Iterate chromosome
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	// We may want to count these (secondary alignments, duplicates, supplementary alignments)
	if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;

	// Get the library information
	std::string rG = "DefaultLib";
	uint8_t *rgptr = bam_aux_get(rec, "RG");
	if (rgptr) {
	  char* rg = (char*) (rgptr + 1);
	  rG = std::string(rg);
	}
	typename TRGMap::iterator itRg = rgMap.find(rG);
	if (itRg == rgMap.end()) {
	  std::cerr << "Missing read group: " << rG << std::endl;
	  return 1;
	}
	
	// Get the read sequence
	std::string sequence;
	sequence.resize(rec->core.l_qseq);
	uint8_t* seqptr = bam_get_seq(rec);
	for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	
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
	for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	  if (bam_cigar_op(cigar[i]) == BAM_CMATCH) {
	    // match or mismatch
	    for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]);++k) {
	      if (sequence[sp] == refslice[rp]) ++itRg->second.matchCount;
	      else ++itRg->second.mismatchCount;
	      // Count bp-level coverage
	      if (itRg->second.cov[rec->core.pos + rp] < itRg->second.maxBpCov) ++itRg->second.cov[rec->core.pos + rp];
	      ++sp;
	      ++rp;
	    }
	  } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	    ++itRg->second.delCount;
	    rp += bam_cigar_oplen(cigar[i]);
	  } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	    ++itRg->second.insCount;
	    sp += bam_cigar_oplen(cigar[i]);
	  } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	    ++itRg->second.softClipCount;
	    sp += bam_cigar_oplen(cigar[i]);
	  } else if(bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	    ++itRg->second.hardClipCount;
	  } else {
	    std::cerr << "Unknown Cigar options" << std::endl;
	    return 1;
	  }
	}
      }
      // clean-up
      bam_destroy1(rec);
      hts_itr_destroy(iter);

      // Summarize bp-level coverage
      for(typename TRGMap::iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
	for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i) ++itRg->second.bpWithCoverage[itRg->second.cov[i]];
	itRg->second.cov.clear();
      }
      if (seq != NULL) free(seq);
    }
    
    // clean-up
    fai_destroy(fai);
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
    
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    
#ifdef PROFILE
    ProfilerStop();
#endif

    // Output statistics
    std::string statFileName = c.outprefix + "." + c.sampleName + ".tsv";
    std::ofstream ofile(statFileName.c_str());
    for(typename TRGMap::iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
      uint64_t alignedbases = itRg->second.matchCount + itRg->second.mismatchCount;
      ofile << "Sample\tLibrary\t#ReferenceBases\t#AlignedBases\t#Coverage\t#MatchedBases\tMatchRate\t#MismatchedBases\tMismatchRate\t#DeletionsCigarD\tDeletionRate\t#InsertionsCigarI\tInsertionRate\t#SoftClippedBases\tSoftClipRate\t#HardClippedBases\tHardClipRate\tErrorRate" << std::endl;
    ofile << c.sampleName << "\t" << itRg->first << "\t" << referencebp << "\t" << alignedbases << "\t" << (double) alignedbases / (double) referencebp << "\t" << itRg->second.matchCount << "\t" << (double) itRg->second.matchCount / (double) alignedbases << "\t" << itRg->second.mismatchCount << "\t" << (double) itRg->second.mismatchCount / (double) alignedbases << "\t" << itRg->second.delCount << "\t" << (double) itRg->second.delCount / (double) alignedbases << "\t" << itRg->second.insCount << "\t" << (double) itRg->second.insCount / (double) alignedbases << "\t" << itRg->second.softClipCount << "\t" << (double) itRg->second.softClipCount / (double) alignedbases << "\t" << itRg->second.hardClipCount << "\t" << (double) itRg->second.hardClipCount / (double) alignedbases << "\t" << (double) (itRg->second.mismatchCount + itRg->second.delCount + itRg->second.insCount + itRg->second.softClipCount + itRg->second.hardClipCount) / (double) alignedbases  << std::endl;
    }
    ofile.close();
    
    // Output coverage histograms
    for(typename TRGMap::iterator itRg = rgMap.begin(); itRg != rgMap.end(); ++itRg) {
      std::string statFileName = c.outprefix + "." + c.sampleName + "." + itRg->first + ".coverage.tsv";
      std::ofstream cfile(statFileName.c_str());
      cfile << "sample\tcoverage\tcount" << std::endl;
      for(uint32_t i = 0; i < itRg->second.bpWithCoverage.size(); ++i) cfile << c.sampleName << "\t" << i << "\t" << itRg->second.bpWithCoverage[i] << std::endl;
      cfile.close();
    }
    
    return 0;
  }

}

#endif
