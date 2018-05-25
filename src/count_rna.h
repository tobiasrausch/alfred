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

#ifndef COUNT_RNA_H
#define COUNT_RNA_H

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
#include "gff3.h"
#include "bed.h"


namespace bamstats
{

  struct CountRNAConfig {
    uint8_t inputFileFormat;   // 0 = gtf, 1 = bed, 2 = gff3
    uint8_t inputBamFormat; // 0 = bam, 1 = bed
    uint16_t stranded;  // 0 = unstranded, 1 = stranded, 2 = stranded (opposite)
    uint16_t minQual;
    std::map<std::string, int32_t> nchr;
    std::string sampleName;
    std::string idname;
    std::string feature;
    boost::filesystem::path gtfFile;
    boost::filesystem::path bedFile;
    boost::filesystem::path bamFile;
    boost::filesystem::path outfile;
  };

  template<typename TConfig, typename TGenomicRegions, typename TFeatureCounter>
  inline int32_t
  bed_counter(TConfig const& c, TGenomicRegions& gRegions, TFeatureCounter& fc) {
    typedef typename TGenomicRegions::value_type TChromosomeRegions;

    // Parse BED file
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "BED file parsing" << std::endl;
    boost::progress_display show_progress(c.nchr.size());

    // Iterate chromosomese
    for(int32_t refIndex=0; refIndex < (int32_t) c.nchr.size(); ++refIndex) {
      ++show_progress;
      if (gRegions[refIndex].empty()) continue;

      // Sort by position
      std::sort(gRegions[refIndex].begin(), gRegions[refIndex].end(), SortIntervalStart<IntervalLabel>());

      // Flag feature positions
      typedef boost::dynamic_bitset<> TBitSet;
      TBitSet featureBitMap(250000000);
      for(uint32_t i = 0; i < gRegions[refIndex].size(); ++i)
	for(int32_t k = gRegions[refIndex][i].start; k < gRegions[refIndex][i].end; ++k) featureBitMap[k] = 1;

      // Count hits
      std::ifstream chrFile(c.bamFile.string().c_str(), std::ifstream::in);
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
	    char strand = '*';
	    if (tokIter != tokens.end()) {
	      ++tokIter; // skip name
	      if (tokIter != tokens.end()) {
		++tokIter; // skip score
		if (tokIter != tokens.end()) {
		  strand = boost::lexical_cast<char>(*tokIter++);
		}
	      }
	    }
	    if (start >= end) continue;  // Bed has right-open intervals
	    typedef std::vector<int32_t> TFeaturePos;
	    TFeaturePos featurepos;
	    for(int32_t i = start; i<end; ++i)
	      if (featureBitMap[i]) featurepos.push_back(i);

	    // Find feature
	    bool ambiguous = false;
	    int32_t featureid = -1;  // No feature by default
	    if (!featurepos.empty()) {
	      int32_t fpfirst = featurepos[0];
	      int32_t fplast = featurepos[featurepos.size()-1];
	      for(typename TChromosomeRegions::const_iterator vIt = gRegions[refIndex].begin(); vIt != gRegions[refIndex].end(); ++vIt) {
		if (vIt->end <= fpfirst) continue;
		if (vIt->start > fplast) break; // Sorted intervals so we can stop searching
		for(TFeaturePos::const_iterator fIt = featurepos.begin(); fIt != featurepos.end(); ++fIt) {
		  if ((vIt->start <= *fIt) && (vIt->end > *fIt) && (featureid != vIt->lid)) {
		    if (c.stranded) {
		      if (c.stranded == 1) {
			if (vIt->strand != strand) continue;
		      } else {
			if (vIt->strand == strand) continue;
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

	    // Check feature agreement
	    if (featureid == -1) continue; // No feature
	    ++fc[featureid];
	  }
	}
	chrFile.close();
      }
    }
    return 0;
  }


  template<typename TConfig, typename TGenomicRegions, typename TFeatureCounter>
  inline int32_t
  bam_counter(TConfig const& c, TGenomicRegions& gRegions, TFeatureCounter& fc) {
    typedef typename TGenomicRegions::value_type TChromosomeRegions;
    
    // Load bam file
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Parse BAM file
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "BAM file parsing" << std::endl;
    boost::progress_display show_progress( hdr->n_targets );

    // Pair qualities and features
    typedef boost::unordered_map<std::size_t, uint8_t> TQualities;
    TQualities qualities;
    typedef boost::unordered_map<std::size_t, int32_t> TFeatures;
    TFeatures features;

    // Iterate chromosomes
    for(int32_t refIndex=0; refIndex < (int32_t) hdr->n_targets; ++refIndex) {
      ++show_progress;
      if (gRegions[refIndex].empty()) continue;

      // Sort by position
      std::sort(gRegions[refIndex].begin(), gRegions[refIndex].end(), SortIntervalStart<IntervalLabel>());
      int32_t maxFeatureLength = 0;
      for(uint32_t i = 0; i < gRegions[refIndex].size(); ++i) {
	if ((gRegions[refIndex][i].end - gRegions[refIndex][i].start) > maxFeatureLength) {
	  maxFeatureLength = gRegions[refIndex][i].end - gRegions[refIndex][i].start;
	}
      }

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
	  typename TChromosomeRegions::const_iterator vIt = std::lower_bound(gRegions[refIndex].begin(), gRegions[refIndex].end(), IntervalLabel(std::max(0, fpfirst - maxFeatureLength)), SortIntervalStart<IntervalLabel>());
	  for(; vIt != gRegions[refIndex].end(); ++vIt) {
	    if (vIt->end <= fpfirst) continue;
	    if (vIt->start > fplast) break; // Sorted intervals so we can stop searching
	    for(TFeaturePos::const_iterator fIt = featurepos.begin(); fIt != featurepos.end(); ++fIt) {
	      if ((vIt->start <= *fIt) && (vIt->end > *fIt) && (featureid != vIt->lid)) {
		if (!_strandOkay(rec, vIt->strand, c.stranded)) continue;
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
    return 0;
  }

  
  template<typename TConfig>
  inline int32_t
  countRNARun(TConfig const& c) {

#ifdef PROFILE
    ProfilerStart("alfred.prof");
#endif

    // Parse GTF file
    typedef std::vector<IntervalLabel> TChromosomeRegions;
    typedef std::vector<TChromosomeRegions> TGenomicRegions;
    TGenomicRegions gRegions;
    gRegions.resize(c.nchr.size(), TChromosomeRegions());
    typedef std::vector<std::string> TGeneIds;
    TGeneIds geneIds;
    int32_t tf = 0;
    if (c.inputFileFormat == 0) tf = parseGTF(c, gRegions, geneIds);
    else if (c.inputFileFormat == 1) tf = parseBED(c, gRegions, geneIds);
    else if (c.inputFileFormat == 2) tf = parseGFF3(c, gRegions, geneIds);
    if (tf == 0) {
      std::cerr << "Error parsing GTF/GFF3/BED file!" << std::endl;
      return 1;
    }

    // Feature counter
    typedef std::vector<int32_t> TFeatureCounter;
    TFeatureCounter fc(tf, 0);
    int32_t retparse = 1;
    if (c.inputBamFormat == 0) retparse = bam_counter(c, gRegions, fc);
    else if (c.inputBamFormat == 1) retparse = bed_counter(c, gRegions, fc);
    if (retparse != 0) {
      std::cerr << "Error feature counting!" << std::endl;
      return 1;
    }

    // Output count table
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
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


  int count_rna(int argc, char **argv) {
    CountRNAConfig c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("map-qual,m", boost::program_options::value<uint16_t>(&c.minQual)->default_value(10), "min. mapping quality")
      ("stranded,s", boost::program_options::value<uint16_t>(&c.stranded)->default_value(0), "strand-specific counting (0: unstranded, 1: stranded, 2: reverse stranded)")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("gene.count"), "output file")
      ;

    boost::program_options::options_description gtfopt("GTF/GFF3 input file options");
    gtfopt.add_options()
      ("gtf,g", boost::program_options::value<boost::filesystem::path>(&c.gtfFile), "gtf/gff3 file")
      ("id,i", boost::program_options::value<std::string>(&c.idname)->default_value("gene_id"), "gtf/gff3 attribute")
      ("feature,f", boost::program_options::value<std::string>(&c.feature)->default_value("exon"), "gtf/gff3 feature")
      ;
    
    boost::program_options::options_description bedopt("BED input file options, columns chr, start, end, name [, score, strand]");
    bedopt.add_options()
      ("bed,b", boost::program_options::value<boost::filesystem::path>(&c.bedFile), "bed file")
      ;

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bamFile), "input bam file")
      ;

    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);

    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(gtfopt).add(bedopt).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(gtfopt).add(bedopt);

    // Parse command-line
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || ((!vm.count("gtf")) && (!vm.count("bed")))) {
      std::cout << std::endl;
      std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] -g <hg19.gtf.gz> <aligned.bam>" << std::endl;
      std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] -b <hg19.bed.gz> <aligned.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    }

    // Check bam file
    if (!(boost::filesystem::exists(c.bamFile) && boost::filesystem::is_regular_file(c.bamFile) && boost::filesystem::file_size(c.bamFile))) {
      std::cerr << "Alignment file is missing: " << c.bamFile.string() << std::endl;
      return 1;
    } else {
      if ((c.bamFile.string().length() > 3) && (c.bamFile.string().substr(c.bamFile.string().length() - 3) == "bed")) {
	c.inputBamFormat = 1;
	c.sampleName = c.bamFile.stem().string();
	std::string oldChr = "";
	typedef std::set<std::string> TChrSet;
	TChrSet chrSet;
	std::ifstream chrFile(c.bamFile.string().c_str(), std::ifstream::in);
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
	      if (chrName != oldChr) chrSet.insert(chrName);
	    }
	  }
	  chrFile.close();
	}
	int32_t refIndex = 0;
	for(TChrSet::iterator itc = chrSet.begin(); itc != chrSet.end(); ++itc, ++refIndex) c.nchr.insert(std::make_pair(*itc, refIndex));
      } else {
	c.inputBamFormat = 0;
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
	for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) c.nchr.insert(std::make_pair(hdr->target_name[refIndex], refIndex));
	
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
    }

    // Check region file
    if (!(boost::filesystem::exists(c.gtfFile) && boost::filesystem::is_regular_file(c.gtfFile) && boost::filesystem::file_size(c.gtfFile))) {
      if (!(boost::filesystem::exists(c.bedFile) && boost::filesystem::is_regular_file(c.bedFile) && boost::filesystem::file_size(c.bedFile))) {
	std::cerr << "Input gtf/bed file is missing." << std::endl;
	return 1;
      } else c.inputFileFormat = 1;
    } else {
      if (is_gff3(c.gtfFile)) c.inputFileFormat = 2;
      else c.inputFileFormat = 0;
    }

    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "alfred ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;

    return countRNARun(c);
  }
  


  
}

#endif
