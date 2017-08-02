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

#ifndef COUNT_H
#define COUNT_H

#include <limits>

#include <boost/icl/split_interval_map.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/progress.hpp>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "version.h"
#include "util.h"
#include "gtf.h"


namespace bamstats
{

  struct CountConfig {
    bool stranded;
    unsigned short minQual;
    std::string sampleName;
    std::string idname;
    std::string feature;
    boost::filesystem::path gtfFile;
    boost::filesystem::path bamFile;
    boost::filesystem::path outfile;
  };


  
  template<typename TConfig>
  inline int32_t
  countRun(TConfig const& c) {

#ifdef PROFILE
    ProfilerStart("alfred.prof");
#endif

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
    int32_t tf = parseGTF(hdr, c.gtfFile, c.feature, c.idname, gRegions, geneIds);
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

      // Flag feature positions
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet featureBitMap(hdr->target_len[refIndex]);
      for(uint32_t i = 0; i < gRegions[refIndex].size(); ++i)
	for(int32_t k = gRegions[refIndex][i].start; k < gRegions[refIndex][i].end; ++k) featureBitMap[k] = 1;

      // Count reads
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* rec = bam_init1();
      int32_t lastAlignedPos = 0;
      std::set<std::size_t> lastAlignedPosReads;
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if ((rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP | BAM_FMUNMAP)) || (rec->core.tid != rec->core.mtid) || (!(rec->core.flag & BAM_FPAIRED))) continue; 

	// Clean-up the read store for identical alignment positions
	if (rec->core.pos > lastAlignedPos) {
	  lastAlignedPosReads.clear();
	  lastAlignedPos = rec->core.pos;
	}

	// Get read sequence
	std::string sequence;
	sequence.resize(rec->core.l_qseq);
	uint8_t* seqptr = bam_get_seq(rec);
	for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

	// Parse CIGAR
	uint32_t* cigar = bam_get_cigar(rec);
	int32_t gp = rec->core.pos; // Genomic position
	int32_t sp = 0; // Sequence position
	typedef std::vector<int32_t> TFeaturePos;
	TFeaturePos featurepos;
	for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	  if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) sp += bam_cigar_oplen(cigar[i]);
	  else if (bam_cigar_op(cigar[i]) == BAM_CINS) sp += bam_cigar_oplen(cigar[i]);
	  else if (bam_cigar_op(cigar[i]) == BAM_CDEL) gp += bam_cigar_oplen(cigar[i]);
	  else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) gp += bam_cigar_oplen(cigar[i]);
	  else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	    //Nop
	  } else if (bam_cigar_op(cigar[i]) == BAM_CMATCH) {
	    for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]); ++k, ++sp, ++gp)
	      if (featureBitMap[gp]) featurepos.push_back(gp);
	  } else {
	    std::cerr << "Unknown Cigar options" << std::endl;
	    return 1;
	  }
	}

	// Find feature
	bool ambiguous = false;
	int32_t featureid = -1;  // No feature by default
	if (!featurepos.empty()) {
	  int32_t fpfirst = featurepos[0];
	  int32_t fplast = featurepos[featurepos.size()-1];
	  for(TChromosomeRegions::const_iterator vIt = gRegions[refIndex].begin(); vIt != gRegions[refIndex].end(); ++vIt) {
	    if (vIt->end <= fpfirst) continue;
	    if (vIt->start > fplast) break; // Sorted intervals so we can stop searching
	    for(TFeaturePos::const_iterator fIt = featurepos.begin(); fIt != featurepos.end(); ++fIt) {
	      if ((vIt->start <= *fIt) && (vIt->end > *fIt) && (featureid != vIt->lid)) {
		if (c.stranded) {
		  if (rec->core.flag & BAM_FREAD1) {
		    if (rec->core.flag & BAM_FREVERSE) {
		      if (vIt->strand != '-') continue;
		    } else {
		      if (vIt->strand != '+') continue;
		    }
		  } else {
		    if (rec->core.flag & BAM_FREVERSE) {
		      if (vIt->strand != '+') continue;
		    } else {
		      if (vIt->strand != '-') continue;
		    }
		  }
		}
		if (featureid == -1) featureid = vIt->lid;
		else {
		  ambiguous = true;
		  break;
		}
	      }
	    }
	  }
	}
	if (ambiguous) continue; // Ambiguous read

	// First or Second Read?
	if ((rec->core.pos < rec->core.mpos) || ((rec->core.pos == rec->core.mpos) && (lastAlignedPosReads.find(hash_string(bam_get_qname(rec))) == lastAlignedPosReads.end()))) {
	  // First read
	  lastAlignedPosReads.insert(hash_string(bam_get_qname(rec)));
	  std::size_t hv = hash_pair(rec);
	  qualities[hv] = rec->core.qual;
	  features[hv] = featureid;
	} else {
	  // Second read
	  std::size_t hv = hash_pair_mate(rec);
	  if (qualities.find(hv) == qualities.end()) continue; // Mate discarded
	  uint8_t pairQuality = std::min((uint8_t) qualities[hv], (uint8_t) rec->core.qual);
	  int32_t featuremate = features[hv];
	  qualities[hv] = 0;
	  features[hv] = -1;

	  // Pair quality
	  if (pairQuality < c.minQual) continue; // Low quality pair

	  // Check feature agreement
	  if ((featureid == -1) && (featuremate == -1)) continue; // No feature
	  else if ((featureid == -1) && (featuremate != -1)) featureid = featuremate;
	  else if ((featureid != -1) && (featuremate == -1)) featuremate = featureid;
	  else {
	    // Both reads have a feature assignment
	    if (featureid != featuremate) continue; // Feature disagreement
	  }

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


  int count(int argc, char **argv) {
    CountConfig c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("map-qual,m", boost::program_options::value<unsigned short>(&c.minQual)->default_value(10), "min. mapping quality")
      ("stranded,s", "strand-specific counting")
      ("gtf,g", boost::program_options::value<boost::filesystem::path>(&c.gtfFile), "gtf file (required)")
      ("id,i", boost::program_options::value<std::string>(&c.idname)->default_value("gene_id"), "GTF attribute")
      ("feature,f", boost::program_options::value<std::string>(&c.feature)->default_value("exon"), "GTF feature")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("gene.count"), "output file")
      ;

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bamFile), "input bam file")
      ;

    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);

    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);

    // Parse command-line
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("gtf"))) {
      printTitle("Alfred");
      std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] -g <hg19.gtf.gz> <aligned.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    }

    // Strand-specific counting
    if (vm.count("stranded")) c.stranded = true;
    else c.stranded = false;
    
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

      // Get sample name
      std::string sampleName;
      if (!getSMTag(std::string(hdr->text), c.bamFile.stem().string(), sampleName)) {
	std::cerr << "Only one sample (@RG:SM) is allowed per input BAM file " << c.bamFile.string() << std::endl;
	return 1;
      } else c.sampleName = sampleName;
      bam_hdr_destroy(hdr);
      hts_idx_destroy(idx);
      sam_close(samfile);
    }

    // Check region file
    if (!(boost::filesystem::exists(c.gtfFile) && boost::filesystem::is_regular_file(c.gtfFile) && boost::filesystem::file_size(c.gtfFile))) {
      std::cerr << "Input gtf file is missing: " << c.gtfFile.string() << std::endl;
      return 1;
    }

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "alfred ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;

    return countRun(c);
  }
  


  
}

#endif
