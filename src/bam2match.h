#ifndef BAM2MATCH_H
#define BAM2MATCH_H

#include <iostream>
#include <vector>
#include <fstream>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "util.h"
#include "variants.h"

namespace bamstats
{

  struct Bam2MatchConfig {
    uint16_t minMapQual;
    boost::filesystem::path genome;
    boost::filesystem::path outfile;
    boost::filesystem::path bamfile;
  };

  template<typename TConfig>
  inline int32_t
  bam2MatchRun(TConfig const& c) {

#ifdef PROFILE
    ProfilerStart("alfred.prof");
#endif

    // Load bam files
    samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.bamfile.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Data out
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    dataOut << "chr\trefstart\trefend\tread\treadstart\treadend\tdirection" << std::endl;
    
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Extract matches" << std::endl;
    boost::progress_display show_progress(hdr->n_targets);

    // Parse BAM alignments
    int32_t refIndex = -1;
    char* seq = NULL;
    faidx_t* fai = fai_load(c.genome.string().c_str());
    bam1_t* rec = bam_init1();
    while (sam_read1(samfile, hdr, rec) >= 0) {
      if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP | BAM_FSECONDARY)) continue;
      if ((rec->core.qual < c.minMapQual) || (rec->core.tid<0)) continue;

      // New chromosome?
      if (rec->core.tid != refIndex) {
	++show_progress;
	if (refIndex != -1) {
	  if (seq != NULL) free(seq);
	}
	refIndex = rec->core.tid;
	// Load chromosome
	int32_t seqlen = -1;
	std::string tname(hdr->target_name[refIndex]);
	seq = faidx_fetch_seq(fai, tname.c_str(), 0, hdr->target_len[refIndex], &seqlen);
      }
      
      // Get read sequence
      std::string sequence;
      sequence.resize(rec->core.l_qseq);
      uint8_t* seqptr = bam_get_seq(rec);
      for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	  
      // Parse CIGAR
      uint32_t* cigar = bam_get_cigar(rec);
      int32_t gp = rec->core.pos; // Genomic position
      int32_t gpStart = -1; //Match start
      int32_t gpEnd = -1; //Match end
      int32_t sp = 0; // Sequence position
      int32_t seqStart = -1;  // Match start
      int32_t seqEnd = -1; // Match end
      for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	  if (seqStart == -1) {
	    seqStart = sp;
	    gpStart = gp;
	  }
	  gp += bam_cigar_oplen(cigar[i]);
	  sp += bam_cigar_oplen(cigar[i]);
	  seqEnd = sp;
	  gpEnd = gp;
	} else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	  if (seqStart == -1) {
	    seqStart = sp;
	    gpStart = gp;
	  }
	  sp += bam_cigar_oplen(cigar[i]);
	  seqEnd = sp;
	  gpEnd = gp;
	} else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	  if (seqStart == -1) {
	    seqStart = sp;
	    gpStart = gp;
	  }
	  gp += bam_cigar_oplen(cigar[i]);
	  seqEnd = sp;
	  gpEnd = gp;
	} else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	  sp += bam_cigar_oplen(cigar[i]);
	} else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	  gp += bam_cigar_oplen(cigar[i]);
	} else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	  sp += bam_cigar_oplen(cigar[i]);
	} else {
	  std::cerr << "Unknown Cigar options" << std::endl;
	  return 1;
	}
      }
      std::string dir = "fwd";
      if (rec->core.flag & BAM_FREVERSE) {
	dir = "rev";
	int32_t seqTmp = seqStart;
	seqStart = sp - seqEnd;
	seqEnd = sp - seqTmp;
      }
      dataOut << hdr->target_name[refIndex] << '\t' << gpStart << '\t' << gpEnd << '\t' << bam_get_qname(rec) << '\t' << seqStart << '\t' << seqEnd << '\t' << dir << std::endl;
    }
    // Close output file
    dataOut.pop();
    dataOut.pop();
    
    // Clean-up
    bam_destroy1(rec);
    if (seq != NULL) free(seq);
    fai_destroy(fai);
      
    // Close bam
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
    
    // End
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;

#ifdef PROFILE
    ProfilerStop();
#endif
    
    return 0;
  }

  int bam2match(int argc, char **argv) {
    Bam2MatchConfig c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("map-qual,m", boost::program_options::value<unsigned short>(&c.minMapQual)->default_value(0), "min. mapping quality")
      ("reference,r", boost::program_options::value<boost::filesystem::path>(&c.genome), "reference fasta file")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("match.gz"), "gzipped output file")
      ;

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bamfile), "input bam file")
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
      std::cout << std::endl;
      std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] -r <ref.fa> <contig.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    }

    // Check input BAM file
    if (vm.count("input-file")) {
      if (!(boost::filesystem::exists(c.bamfile) && boost::filesystem::is_regular_file(c.bamfile) && boost::filesystem::file_size(c.bamfile))) {
	std::cerr << "Input BAM file is missing: " << c.bamfile.string() << std::endl;
	return 1;
      }
      samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
      if (samfile == NULL) {
	std::cerr << "Fail to open file " << c.bamfile.string() << std::endl;
	return 1;
      }
      hts_idx_t* idx = sam_index_load(samfile, c.bamfile.string().c_str());
      if (idx == NULL) {
	std::cerr << "Fail to open index for " << c.bamfile.string() << std::endl;
	return 1;
      }
      bam_hdr_t* hdr = sam_hdr_read(samfile);
      if (hdr == NULL) {
	std::cerr << "Fail to open header for " << c.bamfile.string() << std::endl;
	return 1;
      }
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
    
    return bam2MatchRun(c);
}



}

#endif
