/*
============================================================================
Alfred
============================================================================
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

#ifndef TRACKS_H
#define TRACKS_H

#include <limits>

#include <boost/icl/split_interval_map.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/progress.hpp>

#include <htslib/sam.h>

#include "version.h"
#include "util.h"


namespace bamstats
{

  struct TrackConfig {
    bool wiggleFormat;
    uint16_t covtype;   // 0: sequencing coverage, 1: spanning coverage, 2: footprints
    uint16_t minQual;
    uint32_t normalize;
    float resolution;
    std::string sampleName;
    std::string format;
    boost::filesystem::path bamFile;
    boost::filesystem::path outfile;
  };

  struct Track {
    uint32_t start;
    uint32_t end;
    double score;
    Track(uint32_t s, uint32_t e, double sc) : start(s), end(e), score(sc) {}
  };
  
  template<typename TConfig>
  inline int32_t
  create_tracks(TConfig const& c) {
    // Load bam file
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Pair qualities and features
    typedef boost::unordered_map<std::size_t, uint8_t> TQualities;
    TQualities qualities;

    // Normalize read-counts
    double normFactor = 1;
    if (c.normalize) {
      uint64_t totalPairs = 0;
      boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
      std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Total read count normalization" << std::endl;
      boost::progress_display show_progress( hdr->n_targets );
      for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
	++show_progress;

	hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	int32_t lastAlignedPos = 0;
	std::set<std::size_t> lastAlignedPosReads;
	while (sam_itr_next(samfile, iter, rec) >= 0) {
	  if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	  if (rec->core.qual < c.minQual) continue;
	  if (rec->core.flag & BAM_FPAIRED) {
	    if ((rec->core.flag & BAM_FMUNMAP) || (rec->core.tid != rec->core.mtid)) continue;
	  
	    // Clean-up the read store for identical alignment positions
	    if (rec->core.pos > lastAlignedPos) {
	      lastAlignedPosReads.clear();
	      lastAlignedPos = rec->core.pos;
	    }
	
	    if ((rec->core.pos < rec->core.mpos) || ((rec->core.pos == rec->core.mpos) && (lastAlignedPosReads.find(hash_string(bam_get_qname(rec))) == lastAlignedPosReads.end()))) {
	      // First read
	      lastAlignedPosReads.insert(hash_string(bam_get_qname(rec)));
	      std::size_t hv = hash_pair(rec);
	      qualities[hv] = rec->core.qual;
	    } else {
	      // Second read
	      std::size_t hv = hash_pair_mate(rec);
	      if (qualities.find(hv) == qualities.end()) continue; // Mate discarded
	      uint8_t pairQuality = std::min((uint8_t) qualities[hv], (uint8_t) rec->core.qual);
	      qualities[hv] = 0;
	      
	      // Pair quality
	      if (pairQuality < c.minQual) continue; // Low quality pair
	      
	  
	      // Get bases
	      uint32_t* cigar = bam_get_cigar(rec);
	      for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
		if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
		  totalPairs += bam_cigar_oplen(cigar[i]);
		}
	      }
	    }
	  } else {
	    // Single end
	    uint32_t* cigar = bam_get_cigar(rec);
	    for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	      if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
		totalPairs += bam_cigar_oplen(cigar[i]);
	      }
	    }
	  }
	}
	// Clean-up
	bam_destroy1(rec);
	hts_itr_destroy(iter);
	qualities.clear();
      }
      // Normalize to 100bp paired-end reads
      normFactor = ((double) ((uint64_t) (c.normalize)) / (double) totalPairs) * 100 * 2;
    }
    
    // Open output file
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    if (c.format == "bedgraph") {
      // bedgraph
      dataOut << "track type=bedGraph name=\"" << c.sampleName << "\" description=\"" << c.sampleName << "\" visibility=full color=44,162,95" << std::endl;
    } else if (c.format == "wiggle") {
      // wiggle
      dataOut << "track type=wiggle_0 name=\"" << c.sampleName << "\" description=\"" << c.sampleName << "\" visibility=full color=44,162,95" << std::endl;  
    } else if (c.format == "raw") {
      dataOut << "chr\tpos\t" << c.sampleName << std::endl;
    } else {
      // bed
      dataOut << "chr\tstart\tend\tid\t" << c.sampleName << std::endl;
    }

    // Iterate chromosomes
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "BAM file parsing" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      ++show_progress;
      if (c.wiggleFormat) dataOut << "fixedStep chrom=" << hdr->target_name[refIndex] << " start=1 step=1" << std::endl;

      // Find valid pairs
      std::set<std::size_t> validPairs;
      {
	hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	int32_t lastAlignedPos = 0;
	std::set<std::size_t> lastAlignedPosReads;
	while (sam_itr_next(samfile, iter, rec) >= 0) {
	  if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	  if (rec->core.qual < c.minQual) continue;
	  if (rec->core.flag & BAM_FPAIRED) {
	    if ((rec->core.flag & BAM_FMUNMAP) || (rec->core.tid != rec->core.mtid)) continue;

	    // Clean-up the read store for identical alignment positions
	    if (rec->core.pos > lastAlignedPos) {
	      lastAlignedPosReads.clear();
	      lastAlignedPos = rec->core.pos;
	    }
	    
	    if ((rec->core.pos < rec->core.mpos) || ((rec->core.pos == rec->core.mpos) && (lastAlignedPosReads.find(hash_string(bam_get_qname(rec))) == lastAlignedPosReads.end()))) {
	      // First read
	      lastAlignedPosReads.insert(hash_string(bam_get_qname(rec)));
	      std::size_t hv = hash_pair(rec);
	      qualities[hv] = rec->core.qual;
	    } else {
	      // Second read
	      std::size_t hv = hash_pair_mate(rec);
	      if (qualities.find(hv) == qualities.end()) continue; // Mate discarded
	      uint8_t pairQuality = std::min((uint8_t) qualities[hv], (uint8_t) rec->core.qual);
	      qualities[hv] = 0;

	      // Pair quality
	      if (pairQuality < c.minQual) continue; // Low quality pair

	      // Insert valid pair
	      validPairs.insert(hash_pair_mate(rec));
	    }
	  } else {
	    // Single-end
	    validPairs.insert(hash_read(rec));
	  }
	}
	// Clean-up
	bam_destroy1(rec);
	hts_itr_destroy(iter);
	qualities.clear();
      }

      // Create Coverage track
      typedef uint16_t TCount;
      uint32_t maxCoverage = std::numeric_limits<TCount>::max();
      typedef std::vector<TCount> TCoverage;
      TCoverage cov;
      TCoverage scov;
      if (c.covtype == 0) cov.resize(hdr->target_len[refIndex], 0);
      else {
	cov.resize(hdr->target_len[refIndex], 0);
	scov.resize(hdr->target_len[refIndex], 0);
      }
      if (validPairs.size()) {
	hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
	bam1_t* rec = bam_init1();
	int32_t lastAlignedPos = 0;
	std::set<std::size_t> lastAlignedPosReads;
	while (sam_itr_next(samfile, iter, rec) >= 0) {
	  if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	  if (rec->core.qual < c.minQual) continue;
	  if (rec->core.flag & BAM_FPAIRED) {
	    if ((rec->core.flag & BAM_FMUNMAP) || (rec->core.tid != rec->core.mtid)) continue;

	    // Clean-up the read store for identical alignment positions
	    if (rec->core.pos > lastAlignedPos) {
	      lastAlignedPosReads.clear();
	      lastAlignedPos = rec->core.pos;
	    }

	    if (c.covtype > 0) {
	      // Spanning coverage
	      if ((rec->core.pos < rec->core.mpos) || ((rec->core.pos == rec->core.mpos) && (lastAlignedPosReads.find(hash_string(bam_get_qname(rec))) == lastAlignedPosReads.end()))) {
		// Do nothing
	      } else {
		std::size_t hv = hash_pair_mate(rec);
		if ((validPairs.find(hv) != validPairs.end()) && (layout(rec) == 2)) {
		  int32_t pStart = rec->core.mpos + 50;
		  int32_t pEnd = rec->core.pos - 50;
		  if (pStart < pEnd) {
		    for(int32_t i = pStart; i < pEnd; ++i) {
		      if (scov[i] < maxCoverage) ++scov[i];
		    }
		  }
		}
	      }
	    }
	    if ((c.covtype == 0) || (c.covtype == 2)) {
	      // Sequence coverage
	      std::size_t hv = 0;
	      if ((rec->core.pos < rec->core.mpos) || ((rec->core.pos == rec->core.mpos) && (lastAlignedPosReads.find(hash_string(bam_get_qname(rec))) == lastAlignedPosReads.end()))) hv = hash_pair(rec);
	      else hv = hash_pair_mate(rec);
	      if (validPairs.find(hv) != validPairs.end()) {

		// Reference pointer
		uint32_t rp = rec->core.pos;
	    
		// Parse the CIGAR
		uint32_t* cigar = bam_get_cigar(rec);
		for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
		  if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
		    // match or mismatch
		    for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]);++k) {
		      if (cov[rp] < maxCoverage) ++cov[rp];
		      ++rp;
		    }
		  }
		  else if (bam_cigar_op(cigar[i]) == BAM_CDEL) rp += bam_cigar_oplen(cigar[i]);
		  else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) rp += bam_cigar_oplen(cigar[i]);
		}
	      }
	    }
	  } else {
	    // Single end, only sequencing coverage is possible
	    if (c.covtype == 0) {
	      std::size_t hv = hash_read(rec);
	      if (validPairs.find(hv) != validPairs.end()) {
		
		// Reference pointer
		uint32_t rp = rec->core.pos;
		
		// Parse the CIGAR
		uint32_t* cigar = bam_get_cigar(rec);
		for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
		  if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
		    // match or mismatch
		    for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]);++k) {
		      if (cov[rp] < maxCoverage) ++cov[rp];
		      ++rp;
		    }
		  }
		  else if (bam_cigar_op(cigar[i]) == BAM_CDEL) rp += bam_cigar_oplen(cigar[i]);
		  else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) rp += bam_cigar_oplen(cigar[i]);
		}
	      }
	    }
	  }
	}
	// Clean-up
	bam_destroy1(rec);
	hts_itr_destroy(iter);

	// Coverage mangling
	if (c.covtype == 1) {
	  // Spanning coverage
	  for(uint32_t i = 0; i < scov.size(); ++i) cov[i] = scov[i];
	} else if (c.covtype == 2) {
	  // slope window
	  uint32_t slopewin = 85;
	  TCoverage covStore(slopewin, 0);	  
	  for(uint32_t i = 0; i < scov.size(); ++i) {
	    TCount oldCovVal = covStore[i % 50];
	    covStore[i % 50] = cov[i];
	    if ((i < slopewin) || (i + slopewin >= scov.size())) cov[i] = 0;
	    else {
	      int32_t spanPeakLeft = (int32_t) scov[i]  - (int32_t) scov[i - slopewin];
	      int32_t spanPeakRight = (int32_t) scov[i]  - (int32_t) scov[i + slopewin];
	      int32_t spanPeak = (spanPeakLeft + spanPeakRight) / 2;
	      int32_t covPeakLeft = (int32_t) cov[i]  - (int32_t) oldCovVal;
	      int32_t covPeakRight = (int32_t) cov[i]  - (int32_t) cov[i + slopewin];
	      int32_t covPeak = (covPeakLeft + covPeakRight) / 2;
	      int32_t footprint = spanPeak - covPeak;
	      if (spanPeak < footprint) footprint = spanPeak;
	      if (footprint < 0) cov[i] = 0;
	      else if (footprint < (int32_t) maxCoverage) cov[i] = (TCount) footprint;
	      else cov[i] = maxCoverage;
	    }
	  }
	}

	// Coverage track
	if (c.wiggleFormat) {
	  if (c.normalize) {
	    for(uint32_t i = 0; i < cov.size(); ++i) dataOut << normFactor * cov[i] << std::endl;
	  } else {
	    for(uint32_t i = 0; i < cov.size(); ++i) dataOut << cov[i] << std::endl;
	  }
	} else if (c.format == "raw") {
	  if (c.normalize) {
	    for(uint32_t i = 0; i < cov.size(); ++i) dataOut << hdr->target_name[refIndex] << '\t' << (i+1) << '\t' << normFactor * cov[i] << std::endl;
	  } else {
	    for(uint32_t i = 0; i < cov.size(); ++i) dataOut << hdr->target_name[refIndex] << '\t' << (i+1) << '\t' << cov[i] << std::endl;
	  }
	} else {
	  typedef std::list<Track> TrackLine;
	  TrackLine tl;
	  uint32_t wb = 0;
	  uint32_t we = 0;
	  double wval = cov[0];
	  for(uint32_t i = 1; i<cov.size(); ++i) {
	    if (cov[i] == wval) ++we;
	    else {
	      tl.push_back(Track(wb, we+1, normFactor * wval));
	      wb = i;
	      we = i;
	      wval = cov[i];
	    }
	  }
	  tl.push_back(Track(wb, we+1, normFactor * wval));
	  
	  // Reduce file size
	  if ((c.resolution > 0) && (c.resolution < 1)) {
	    double red = 1;
	    uint32_t origs = tl.size();
	    while ((tl.size() > 1) && (red > c.resolution)) {
	      TrackLine::iterator idx = tl.begin();
	      TrackLine::iterator idxNext = tl.begin();
	      ++idxNext;
	      std::vector<double> errs;
	      for(;idxNext != tl.end(); ++idx, ++idxNext) {
		uint32_t w1 = idx->end - idx->start;
		uint32_t w2 = idxNext->end - idxNext->start;
		double nwavg = (w1 * idx->score + w2 * idxNext->score) / (w1 + w2);
		double nerr = w1 * ((idx->score - nwavg) * (idx->score - nwavg));
		nerr += w2 * ((idxNext->score - nwavg) * (idxNext->score - nwavg));
		errs.push_back(nerr);
	      }
	      std::sort(errs.begin(), errs.end());
	      uint32_t bpidx = (red - c.resolution) * tl.size();
	      if (bpidx > 0) bpidx = bpidx - 1;
	      double thres = errs[bpidx];
	      idx = tl.begin();
	      idxNext = tl.begin();
	      ++idxNext;
	      while(idxNext != tl.end()) {
		uint32_t w1 = idx->end - idx->start;
		uint32_t w2 = idxNext->end - idxNext->start;
		double nwavg = (w1 * idx->score + w2 * idxNext->score) / (w1 + w2);
		double nerr = w1 * ((idx->score - nwavg) * (idx->score - nwavg));
		nerr += w2 * ((idxNext->score - nwavg) * (idxNext->score - nwavg));
		if (nerr <= thres) {
		  ++idxNext;
		  uint32_t oldst = idx->start;
		  tl.erase(idx++);
		  idx->start = oldst;
		  idx->score = nwavg;
		} else {
		  ++idxNext;
		  ++idx;
		}
	      }
	      red = (double) tl.size() / (double) origs;
	    }
	  }
	  if (c.format == "bedgraph") {
	    for(TrackLine::iterator idx = tl.begin(); idx != tl.end(); ++idx) dataOut << hdr->target_name[refIndex] << "\t" << idx->start << "\t" << idx->end << "\t" << idx->score << std::endl;
	  } else {
	    for(TrackLine::iterator idx = tl.begin(); idx != tl.end(); ++idx) dataOut << hdr->target_name[refIndex] << "\t" << idx->start << "\t" << idx->end << "\t" << hdr->target_name[refIndex] << ":" << idx->start << "-" << idx->end << "\t" << idx->score << std::endl;
	  }
	}
      }
    }
    
    // clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
    dataOut.pop();
    
    return 0;
  }


  int tracks(int argc, char **argv) {
    TrackConfig c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("map-qual,m", boost::program_options::value<uint16_t>(&c.minQual)->default_value(10), "min. mapping quality")
      ("normalize,n", boost::program_options::value<uint32_t>(&c.normalize)->default_value(30000000), "#pairs to normalize to (0: no normalization)")
      ("covtype,c", boost::program_options::value<uint16_t>(&c.covtype)->default_value(0), "coverage type (0: sequencing coverage, 1: spanning coverage, 2: footprinting)")
      ;

    boost::program_options::options_description resolution("Resolution options (bedgraph/bed format)");
    resolution.add_options()
      ("resolution,r", boost::program_options::value<float>(&c.resolution)->default_value(0.2), "fractional resolution ]0,1]")
      ;      

    boost::program_options::options_description otp("Output options");
    otp.add_options()
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("track.gz"), "track file")
      ("format,f", boost::program_options::value<std::string>(&c.format)->default_value("bedgraph"), "output format [bedgraph|bed|wiggle|raw]")
      ;

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bamFile), "input bam file")
      ;

    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);

    // Set the visibility
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(resolution).add(otp).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(resolution).add(otp);

    // Parse command-line
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file"))) {
      std::cout << std::endl;
      std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] <aligned.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    }

    // Wiggle format at 1bp resolution
    c.wiggleFormat = false;
    if (vm.count("format")) {
      if (c.format == "wiggle") {
	std::cerr << "Warning: Wiggle format triggers single-basepair resolution coverage track!" << std::endl;
	c.wiggleFormat = true;
      }
    }

    // Coverage type
    if (c.covtype > 2) c.covtype = 0;

    // Check bam file
    if (!(boost::filesystem::exists(c.bamFile) && boost::filesystem::is_regular_file(c.bamFile) && boost::filesystem::file_size(c.bamFile))) {
      std::cerr << "Alignment file is missing: " << c.bamFile.string() << std::endl;
      return 1;
    } else {
      samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
      if (samfile == NULL) {
	std::cerr << "Fail to open file " << c.bamFile.string() << std::endl;
	return 1;
      }
      hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
      if (idx == NULL) {
	if (bam_index_build(c.bamFile.string().c_str(), 0) != 0) {
	  std::cerr << "Fail to open index for " << c.bamFile.string() << std::endl;
	  return 1;
	}
      }
      bam_hdr_t* hdr = sam_hdr_read(samfile);
      if (hdr == NULL) {
	std::cerr << "Fail to open header for " << c.bamFile.string() << std::endl;
	return 1;
      }

      // Get sample name
      std::string sampleName;
      if (!getSMTag(std::string(hdr->text), c.bamFile.stem().string(), sampleName)) {
	std::cerr << "Only one sample (@RG:SM) is allowed per input BAM file " << c.bamFile.string() << std::endl;
	return 1;
      } else c.sampleName = sampleName;

      // Clean-up
      bam_hdr_destroy(hdr);
      hts_idx_destroy(idx);
      sam_close(samfile);
    }

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "alfred ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;

    return create_tracks(c);
  }

  
}

#endif
