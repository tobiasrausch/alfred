#ifndef TELMOTIF_H
#define TELMOTIF_H

#include <iostream>
#include <vector>
#include <fstream>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#ifdef PROFILE
#include "gperftools/profiler.h"
#endif

namespace bamstats {

  struct TelMotifConfig {
    uint16_t minSeqQual;
    uint16_t wsize;
    boost::filesystem::path outfile;
    boost::filesystem::path genome;
    boost::filesystem::path bam;
  };


  
  template<typename TConfig>
  inline int32_t
  runTelMotif(TConfig& c) {

    std::vector<std::string> motifs = { "CCCTCACCCTAACCCTCA", "CCCTGACCCTGACCCCAA", "CCCCAACCCTAACCCTCA", "CCCCAACCCTAACCCTGA", "TCAGGGTTAGGGTTGGGG", "TGAGGGTGAGGGTCAGGG", "CCCTAACCCTGACCCTAA", "TTGGGGTTGGGGTTGGGG", "CCCTAACCCTCACCCTGA", "TTGGGGTCAGGGTTGGGG", "CCCTCACCCCAACCCCAA", "TTAGGGTCAGGGTTAGGG", "CCCTGACCCTAACCCTAA", "TGAGGGTCAGGGTTAGGG", "CCCCAACCCTCACCCTAA", "CCCTCACCCCAACCCTCA", "CCCCAACCCTCACCCTGA", "TCAGGGTCAGGGTTGGGG", "TGAGGGTCAGGGTTGGGG", "TTGGGGTTAGGGTCAGGG", "TTGGGGTTAGGGTTAGGG", "CCCTCACCCTCACCCTCA", "CCCTGACCCTAACCCCAA", "CCCCAACCCTGACCCTCA", "TTAGGGTTGGGGTCAGGG", "CCCTGACCCTGACCCTAA", "CCCCAACCCTAACCCTAA", "TGAGGGTTGGGGTTGGGG", "CCCTAACCCTAACCCTGA", "CCCTGACCCTCACCCTAA", "TTAGGGTTGGGGTTGGGG", "TCAGGGTTGGGGTGAGGG", "CCCTCACCCTAACCCCAA", "CCCCAACCCCAACCCTCA", "TCAGGGTTGGGGTTAGGG", "TCAGGGTTAGGGTCAGGG", "TTAGGGTGAGGGTTGGGG", "TTGGGGTTAGGGTTGGGG", "CCCCAACCCTCACCCCAA", "TTAGGGTGAGGGTGAGGG", "CCCTGACCCCAACCCTCA", "CCCTCACCCTAACCCTAA", "CCCTGACCCCAACCCCAA", "TCAGGGTCAGGGTTAGGG", "TGAGGGTTGGGGTGAGGG", "CCCTCACCCTCACCCCAA", "CCCTAACCCTGACCCTGA", "CCCCAACCCTGACCCTGA", "CCCTAACCCTAACCCTAA", "CCCTGACCCTAACCCTGA", "CCCTCACCCCAACCCTAA", "CCCCAACCCTGACCCCAA", "TGAGGGTTGGGGTCAGGG", "CCCTAACCCCAACCCTGA", "TCAGGGTGAGGGTCAGGG", "TTGGGGTTAGGGTGAGGG", "TTAGGGTGAGGGTCAGGG", "TGAGGGTTAGGGTTGGGG", "TTGGGGTGAGGGTTGGGG", "TGAGGGTTGGGGTTAGGG", "TTAGGGTTGGGGTGAGGG", "CCCCAACCCTGACCCTAA", "TTAGGGTCAGGGTGAGGG", "CCCTGACCCTCACCCTCA", "TGAGGGTTAGGGTTAGGG", "TTAGGGTTGGGGTTAGGG", "TCAGGGTCAGGGTGAGGG", "CCCCAACCCTCACCCTCA", "TTGGGGTGAGGGTTAGGG", "TCAGGGTTGGGGTCAGGG", "TTGGGGTTGGGGTTAGGG", "TTAGGGTTAGGGTTGGGG", "CCCTAACCCTCACCCTCA", "CCCTCACCCCAACCCTGA", "TTAGGGTTAGGGTGAGGG", "TTGGGGTCAGGGTGAGGG", "TTGGGGTTGGGGTCAGGG", "CCCTCACCCTAACCCTGA", "TTGGGGTCAGGGTCAGGG", "TCAGGGTGAGGGTGAGGG", "CCCTGACCCTGACCCTGA", "TTGGGGTCAGGGTTAGGG", "CCCTGACCCCAACCCTAA", "TTAGGGTTAGGGTTAGGG", "TGAGGGTCAGGGTCAGGG", "CCCTCACCCTGACCCTGA", "TTAGGGTCAGGGTTGGGG", "CCCTAACCCTAACCCTCA", "CCCTAACCCCAACCCTAA", "CCCTAACCCTGACCCCAA", "TCAGGGTCAGGGTCAGGG", "CCCTCACCCTCACCCTGA", "TTAGGGTGAGGGTTAGGG", "CCCCAACCCCAACCCTGA", "CCCTAACCCTCACCCCAA", "CCCTGACCCTCACCCTGA", "TTGGGGTGAGGGTGAGGG", "TTGGGGTGAGGGTCAGGG", "TCAGGGTTGGGGTTGGGG", "TGAGGGTTAGGGTGAGGG", "TCAGGGTTAGGGTTAGGG", "TGAGGGTCAGGGTGAGGG", "CCCCAACCCTAACCCCAA", "TGAGGGTGAGGGTTAGGG", "CCCTAACCCTCACCCTAA", "TGAGGGTGAGGGTGAGGG", "TTAGGGTCAGGGTCAGGG", "CCCTAACCCTAACCCCAA", "CCCTAACCCTGACCCTCA", "TCAGGGTTAGGGTGAGGG", "CCCTAACCCCAACCCCAA", "CCCTGACCCTAACCCTCA", "TCAGGGTGAGGGTTGGGG", "TTAGGGTTAGGGTCAGGG", "CCCTGACCCCAACCCTGA", "CCCCAACCCCAACCCCAA", "CCCTCACCCTCACCCTAA", "CCCTCACCCTGACCCCAA", "TGAGGGTGAGGGTTGGGG", "CCCTCACCCTGACCCTAA", "TCAGGGTGAGGGTTAGGG", "CCCTAACCCCAACCCTCA", "TGAGGGTTAGGGTCAGGG", "CCCTGACCCTGACCCTCA", "CCCCAACCCCAACCCTAA", "TTGGGGTTGGGGTGAGGG", "CCCTGACCCTCACCCCAA", "CCCTCACCCTGACCCTCA" };
    int32_t seqMotifSize = motifs[0].size();
    
#ifdef PROFILE
    ProfilerStart("telmotif.prof");
#endif

    // Open file handles
    samFile* samfile = sam_open(c.bam.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.bam.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Initialize counters
    typedef std::vector<uint16_t> TChrCount;
    typedef std::vector<TChrCount> TGenomeCount;
    TGenomeCount gcount(hdr->n_targets, TChrCount());
    for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
      gcount[refIndex].resize((int32_t) (hdr->target_len[refIndex] / c.wsize) + 1, 0);
    }

    // Parse BAM
    bam1_t* rec = bam_init1();
    while (sam_read1(samfile, hdr, rec) >= 0) {
      if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
      int32_t refidx = 0;
      int32_t anchor = 0;
      if (rec->core.flag & BAM_FUNMAP) {
	if ((rec->core.flag & BAM_FMUNMAP) || (rec->core.mtid < 0)) continue; // Both unmapped
	anchor = rec->core.mpos;
	refidx = rec->core.mtid;
      } else {
	anchor = rec->core.pos;
	refidx = rec->core.tid;
      }
	
      // Load sequence and quality
      typedef std::vector<uint8_t> TQuality;
      TQuality quality(rec->core.l_qseq);
      std::string sequence(rec->core.l_qseq, 'N');
      uint8_t* seqptr = bam_get_seq(rec);
      uint8_t* qualptr = bam_get_qual(rec);
      for (int32_t i = 0; i < rec->core.l_qseq; ++i) {
        quality[i] = qualptr[i];
        sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
      }

      // Search motifs, count only first hit
      for(uint32_t i = 0; i < motifs.size(); ++i) {
	std::size_t pos = sequence.find(motifs[i]);
	if (pos != std::string::npos) {
	  // Sufficient quality?
	  int32_t avgqual = 0;
	  for(uint32_t k = pos; ((k < pos + seqMotifSize) && (k < quality.size())); ++k) avgqual += (int32_t) quality[k];
	  avgqual /= seqMotifSize;
	  if (avgqual > c.minSeqQual) {
	    int32_t hitpos = anchor + pos;
	    ++gcount[refidx][(int) (hitpos/c.wsize)];
	    break; // Max. one hit per sequence
	  }
	}
      }
    }
    bam_destroy1(rec);


    // Output
    std::ofstream ofile(c.outfile.string().c_str());
    for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
      for(uint32_t i = 0; i < gcount[refIndex].size(); ++i) {
	if (gcount[refIndex][i] > 0) {
	  ofile << hdr->target_name[refIndex] << '\t' << i * c.wsize << '\t' << (i+1) * c.wsize << '\t' << gcount[refIndex][i] << std::endl;
	}
      }
    }
    ofile.close();
    
    // Clean-up
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
    
#ifdef PROFILE
    ProfilerStop();
#endif

    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "Done." << std::endl;
    
    return 0;
  }
  

  int telmotif(int argc, char **argv) {
    TelMotifConfig c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("quality,q", boost::program_options::value<uint16_t>(&c.minSeqQual)->default_value(20), "min. sequence quality")
      ("wsize,w", boost::program_options::value<uint16_t>(&c.wsize)->default_value(1000), "window size")
      ("reference,r", boost::program_options::value<boost::filesystem::path>(&c.genome), "reference fasta file (required)")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("neotelomere.bed"), "output file")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bam), "input bam file")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);
    
    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);
    
    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("reference"))) {
      std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] -r <ref.fa> -o outprefix <input.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    } 
    
    // Check genome
    if (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome))) {
      std::cerr << "Input reference file is missing: " << c.genome.string() << std::endl;
      return 1;
    } else {
      faidx_t* fai = fai_load(c.genome.string().c_str());
      if (fai == NULL) {
	if (fai_build(c.genome.string().c_str()) == -1) {
	  std::cerr << "Fail to open genome fai index for " << c.genome.string() << std::endl;
	  return 1;
	} else fai = fai_load(c.genome.string().c_str());
      }
      fai_destroy(fai);
    }
    
    // Check bam files
    if (vm.count("input-file")) {
      if (!(boost::filesystem::exists(c.bam) && boost::filesystem::is_regular_file(c.bam) && boost::filesystem::file_size(c.bam))) {
	std::cerr << "Alignment file is missing: " << c.bam.string() << std::endl;
	return 1;
      }
      samFile* samfile = sam_open(c.bam.string().c_str(), "r");
      if (samfile == NULL) {
	std::cerr << "Fail to open file " << c.bam.string() << std::endl;
	return 1;
      }
      hts_idx_t* idx = sam_index_load(samfile, c.bam.string().c_str());
      if (idx == NULL) {
	std::cerr << "Fail to open index for " << c.bam.string() << std::endl;
	return 1;
      }
      bam_hdr_t* hdr = sam_hdr_read(samfile);
      if (hdr == NULL) {
	std::cerr << "Fail to open header for " << c.bam.string() << std::endl;
	return 1;
      }
      faidx_t* fai = fai_load(c.genome.string().c_str());
      for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
	std::string tname(hdr->target_name[refIndex]);
	if (!faidx_has_seq(fai, tname.c_str())) {
	  std::cerr << "BAM file chromosome " << hdr->target_name[refIndex] << " is NOT present in your reference file " << c.genome.string() << std::endl;
	  return 1;
	}
      }
      fai_destroy(fai);
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
    
    return runTelMotif(c);
  }

}

#endif

