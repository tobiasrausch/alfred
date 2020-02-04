#ifndef REPLISEQPROC_H
#define REPLISEQPROC_H

#include <limits>

#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/progress.hpp>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "util.h"

namespace bamstats
{

  template<typename TConfig>
  inline int32_t
  repliseqRun(TConfig const& c) {
    // Open file handles
    typedef std::vector<samFile*> TSamFile;
    typedef std::vector<hts_idx_t*> TIndex;
    TSamFile samfile(c.files.size());
    TIndex idx(c.files.size());
    for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) {
      samfile[file_c] = sam_open(c.files[file_c].string().c_str(), "r");
      idx[file_c] = sam_index_load(samfile[file_c], c.files[file_c].string().c_str());
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile[0]);

    // Genomic counts
    typedef std::vector<int32_t> TBinCount;
    typedef std::vector<TBinCount> TGenomicCount;
    typedef std::vector<TGenomicCount> TFileCounts;
    TFileCounts fc(c.files.size(), TGenomicCount());
    for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) fc[file_c].resize(hdr->n_targets, TBinCount());
    std::vector<int32_t> totalByFile(c.files.size(), 0);
    
    // Parse reference and BAM file
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "BAM file parsing" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );

    // Parse genome
    faidx_t* fai = fai_load(c.genome.string().c_str());
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      ++show_progress;

      // Fetch sequence
      char* seq = NULL;
      int32_t seqlen = -1;
      std::string tname(hdr->target_name[refIndex]);
      seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);

      // Fetch Ns
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet nrun(hdr->target_len[refIndex]);
      for(uint32_t i = 0; i < hdr->target_len[refIndex]; ++i)
	if ((seq[i] == 'n') || (seq[i] == 'N')) nrun[i] = 1;
      
      
      for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
	// Set up fragment counter
	typedef uint8_t TCountType;
	int32_t maxCount = std::numeric_limits<TCountType>::max();
	typedef std::vector<TCountType> TChrCounts;
	TChrCounts cc(hdr->target_len[refIndex], 0);

	// Iterate bam
	hts_itr_t* iter = sam_itr_queryi(idx[file_c], refIndex, 0, hdr->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	while (sam_itr_next(samfile[file_c], iter, rec) >= 0) {
	  if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	  if (rec->core.qual < c.minq) continue;

	  int32_t midPoint = 0;
	  if (rec->core.flag & BAM_FPAIRED) {
	    if ((rec->core.flag & BAM_FMUNMAP) || (rec->core.pos < rec->core.mpos) || (rec->core.tid != rec->core.mtid)) continue;
	    int32_t outerISize = rec->core.pos - rec->core.mpos + rec->core.l_qseq;
	    if (outerISize < 1000) midPoint = rec->core.pos + outerISize / 2;
	    else midPoint = rec->core.pos + halfAlignmentLength(rec);
	  } else {
	    midPoint = rec->core.pos + halfAlignmentLength(rec);
	  }
	  if ((midPoint >= 0) && (midPoint < (int32_t) hdr->target_len[refIndex])) {
	    ++totalByFile[file_c];
	    if (cc[midPoint] < maxCount) ++cc[midPoint];
	  }
	}
	bam_destroy1(rec);
	hts_itr_destroy(iter);

	// Summarize counts
	for(int32_t i = 0; (i + c.wsize) < (int32_t) hdr->target_len[refIndex]; i = i + c.step) {
	  int32_t sumf = 0;
	  int32_t nsum = 0;
	  for(int32_t k = i; k < i + c.wsize; ++k) {
	    sumf += cc[k];
	    nsum += nrun[k];
	  }
	  // Blacklist windows with Ns
	  if (!nsum) fc[file_c][refIndex].push_back(sumf);
	  else fc[file_c][refIndex].push_back(0);
	}
      }
      if (seq != NULL) free(seq);
    }
    fai_destroy(fai);

    // Median normalize counts
    std::sort(totalByFile.begin(), totalByFile.end());
    int32_t med = totalByFile[totalByFile.size() / 2];
    for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) {
      double corf = (double) med / (double) totalByFile[file_c];
      for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
	for(uint32_t k = 0; k < fc[file_c][refIndex].size(); ++k) {
	  fc[file_c][refIndex][k] = (int32_t) (fc[file_c][refIndex][k] * corf);
	}
      }
    }

    // Percent normalized values
    typedef std::vector<double> TWindows;
    typedef std::vector<TWindows> TGenomicWindows;
    TGenomicWindows gw(hdr->n_targets, TWindows());
    TWindows all;
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) gw[refIndex].resize(fc[0][refIndex].size(), 0);
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      for(uint32_t k = 0; k < fc[0][refIndex].size(); ++k) {
	for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) gw[refIndex][k] += fc[file_c][refIndex][k];
	all.push_back(gw[refIndex][k]);
      }
    }
    std::sort(all.begin(), all.end());
    double medrep = all[all.size() / 2];
    all.clear();
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      for(uint32_t k = 0; k < fc[0][refIndex].size(); ++k) {
	double corf = 1.0;
	if (gw[refIndex][k] != 0) corf = medrep / (double) gw[refIndex][k];
	for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) fc[file_c][refIndex][k] = (int32_t) (fc[file_c][refIndex][k] * corf);
      }
    }

    // Replication track
    std::vector<double> magicformula;
    magicformula.push_back(0.917);
    magicformula.push_back(0.75);
    magicformula.push_back(0.583);
    magicformula.push_back(0.417);
    magicformula.push_back(0.25);
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      for(uint32_t k = 0; k < gw[refIndex].size(); ++k) {
	gw[refIndex][k] = 0;
	for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) {
	  if (file_c < magicformula.size()) {
	    gw[refIndex][k] += magicformula[file_c] * fc[file_c][refIndex][k];
	  }
	}
      }
    }

    // Moving avg. smoothing
    int32_t smoothw = 75;
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      typename TGenomicWindows::value_type tmpgw(gw[refIndex].size(), 0);
      for(int32_t k = 0; k < (int32_t) gw[refIndex].size(); ++k) {
	double mavg = 0;
	int32_t cavg = 0;
	int32_t ks = std::max(0, k-smoothw);
	int32_t ke = std::min((int32_t) gw[refIndex].size(), k+smoothw);
	for(int32_t ki = ks; ki<ke; ++ki) {
	  if (gw[refIndex][ki] != 0) {
	    mavg += gw[refIndex][ki];
	    ++cavg;
	  }
	}
	tmpgw[k] = mavg / (double) cavg;
      }
      double maxVal = 0;
      double minVal = 40000000;
      for(int32_t k = 0; k < (int32_t) gw[refIndex].size(); ++k) {
	if (gw[refIndex][k] != 0) {
	  gw[refIndex][k] = tmpgw[k];
	  if (gw[refIndex][k] < minVal) minVal = gw[refIndex][k];
	  if (gw[refIndex][k] > maxVal) maxVal = gw[refIndex][k];
	} else gw[refIndex][k] = -1;
      }
      // Normalize to [0,1]
      for(int32_t k = 0; k < (int32_t) gw[refIndex].size(); ++k) {
	if (gw[refIndex][k] != -1) gw[refIndex][k] = (gw[refIndex][k] - minVal) / (maxVal - minVal);
      }
    }

    // Output profile
    std::string statFileName = c.outprefix + ".profile.tsv";
    std::ofstream pfile(statFileName.c_str());
    pfile << "chr\tpos";
    for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) pfile << "\t" << c.files[file_c].stem().string();
    pfile << std::endl;
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      if (!fc[0][refIndex].empty()) {
	for(uint32_t k = 0; k < fc[0][refIndex].size(); ++k) {
	  pfile << hdr->target_name[refIndex] << '\t' << k * c.step + c.wsize / 2;
	  for(uint32_t file_c = 0; file_c < c.files.size(); ++file_c) pfile << '\t' << fc[file_c][refIndex][k];
	  pfile << std::endl;
	}
      }
    }
    pfile.close();

    // Output replication timing (higher values correspond to earlier replication)
    statFileName = c.outprefix + ".reptime.tsv";
    std::ofstream rfile(statFileName.c_str());
    rfile << "chr\tpos\treptime" << std::endl;
    for(int32_t refIndex = 0; refIndex < hdr->n_targets; ++refIndex) {
      for(uint32_t k = 0; k < gw[refIndex].size(); ++k) {
	rfile << hdr->target_name[refIndex] << '\t' << k * c.step + c.wsize / 2 << '\t' << gw[refIndex][k] << std::endl;
      }
    }
    rfile.close();
    
    // clean-up
    bam_hdr_destroy(hdr);
    for(unsigned int file_c = 0; file_c < c.files.size(); ++file_c) {
      hts_idx_destroy(idx[file_c]);
      sam_close(samfile[file_c]);
    }
    
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    
#ifdef PROFILE
    ProfilerStop();
#endif


    return 0;
  }

}

#endif
