#ifndef CONSENSUS_H
#define CONSENSUS_H

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem.hpp>
#include <boost/tokenizer.hpp>
#include <iostream>
#include <vector>
#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <math.h>
#include <stdio.h>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "util.h"
#include "align.h"
#include "msa.h"


namespace bamstats {


  struct ConfigConsensus {
    bool secondary;
    bool trimreads;
    bool hasGenome;
    uint16_t minMapQual;
    uint32_t window;
    int32_t gapopen;
    int32_t gapext;
    int32_t match;
    int32_t mismatch;
    float fractionCalled;
    std::string position;
    std::string seqtype;
    std::string outformat;
    DnaScore<int32_t> aliscore;
    boost::filesystem::path alignment;
    boost::filesystem::path genome;
    boost::filesystem::path consensus;
    boost::filesystem::path inputfile;
  };


  template<typename TConfig>
  inline void
  _loadFastaReads(TConfig const& c, std::vector<std::string>& rs) {
    // Load FASTA/FASTQ
    rs.clear();
    std::ifstream fqfile;
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    if (is_gz(c.inputfile)) {
      fqfile.open(c.inputfile.string().c_str(), std::ios_base::in | std::ios_base::binary);
      dataIn.push(boost::iostreams::gzip_decompressor(), 16*1024);
    } else fqfile.open(c.inputfile.string().c_str(), std::ios_base::in);
    dataIn.push(fqfile);
    std::istream instream(&dataIn);
    std::string gline;
    uint64_t lnum = 0;
    std::string qname;
    bool validRec = true;
    while(std::getline(instream, gline)) {
      if (lnum % 2 == 0) {
	// FASTA or FASTQ
        if ((gline[0] == '>') || (gline[0] == '@')) {
          validRec = true;
          qname = gline.substr(1);
          qname = qname.substr(0, qname.find(' '));
          qname = qname.substr(0, qname.find('\t'));
        } else validRec = false;
      } else if (lnum % 2 == 1) {
        if (validRec) {
	  std::cout << "Read name: " << qname << ", Length: " << gline.size() << std::endl;	  
	  rs.push_back(gline);
        }
      }
      ++lnum;
    }
    // Clean-up
    dataIn.pop();
    if (is_gz(c.inputfile)) dataIn.pop();
    fqfile.close();
  }

  inline bool
  findReadTrim(bam1_t const* rec, int32_t const start, int32_t const end, int32_t& leftPos, int32_t& rightPos) {
    int32_t rp = rec->core.pos; // reference pointer
    int32_t sp = 0; // sequence pointer
    leftPos = -1;
    rightPos = -1;

    // Parse the CIGAR
    if (start < end) {
      uint32_t* cigar = bam_get_cigar(rec);
      for (std::size_t i = 0; i < rec->core.n_cigar; ++i) {
	if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
	  // match or mismatch
	  for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]);++k) {
	    ++sp;
	    ++rp;
	    if ((leftPos == -1) && (rp >= start)) leftPos = sp;
	    if ((rightPos == -1) && (rp >= end)) {
	      rightPos = sp;
	      return true;
	    }
	  }
	} else {
	  if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
	    rp += bam_cigar_oplen(cigar[i]);
	  } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
	    sp += bam_cigar_oplen(cigar[i]);
	  } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
	    sp += bam_cigar_oplen(cigar[i]);
	  } else if(bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
	    // Nothing
	  } else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) {
	    rp += bam_cigar_oplen(cigar[i]);
	  } else {
	    std::cerr << "Unknown Cigar options" << std::endl;
	    return 1;
	  }
	  if ((leftPos == -1) && (rp >= start)) leftPos = sp;
	  if ((rightPos == -1) && (rp >= end)) {
	    rightPos = sp;
	    return true;
	  }
	}
      }
    }
    return false;
  }
 
  template<typename TConfig>
  inline bool
  _loadBamReads(TConfig const& c, std::vector<std::string>& rs) {
    if (!(boost::filesystem::exists(c.inputfile) && boost::filesystem::is_regular_file(c.inputfile) && boost::filesystem::file_size(c.inputfile))) {
      std::cerr << "Alignment file is missing: " << c.inputfile.string() << std::endl;
      return false;
    }
    samFile* samfile = sam_open(c.inputfile.string().c_str(), "r");
    if (samfile == NULL) {
      std::cerr << "Fail to open file " << c.inputfile.string() << std::endl;
      return false;
    }
    if (c.hasGenome) hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.inputfile.string().c_str());
    if (idx == NULL) {
      std::cerr << "Fail to open index for " << c.inputfile.string() << std::endl;
      return false;
    }
    bam_hdr_t* hdr = sam_hdr_read(samfile);
    if (hdr == NULL) {
      std::cerr << "Fail to open header for " << c.inputfile.string() << std::endl;
      return false;
    }
    
    // Parse position
    typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
    boost::char_separator<char> sep(":");
    Tokenizer tokens(c.position, sep);
    Tokenizer::iterator tokIter = tokens.begin();
    bool posError = true;
    int32_t refIndex = -1;
    int32_t pos = -1;
    if (tokIter!=tokens.end()) {
      std::string chrName = *tokIter++;
      refIndex = bam_name2id(hdr, chrName.c_str());
      if (refIndex >= 0) {
	pos = boost::lexical_cast<int32_t>(*tokIter++);
	if ((pos >= 0) && (pos < (int32_t) hdr->target_len[refIndex])) posError = false;
      }
    }
    if (posError) {
      std::cerr << "Position needs to be in the format chr:pos" << std::endl;
      return false;
    }
    
    // Fetch reads
    std::set<unsigned> read_set;
    typedef boost::unordered_map<unsigned, bool> TMissingReads;
    TMissingReads missing_reads;
    {
      std::cout << "Primary alignments" << std::endl;
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, pos, pos+1);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
	if (rec->core.qual < c.minMapQual) continue;

	// Secondary Alignments
	if (rec->core.flag & BAM_FSECONDARY) {
	  if (c.secondary) {
	    // No sequence information
	    unsigned seed = hash_string(bam_get_qname(rec));
	    if (read_set.find(seed) == read_set.end()) missing_reads[seed] = (rec->core.flag & BAM_FREVERSE);
	  }
	} else {
	  // Overlaps a minimal window?
	  if (rec->core.pos + (int32_t) c.window <= pos) {
	    if (rec->core.pos + alignmentLength(rec) >= pos + c.window) {
	      unsigned seed = hash_string(bam_get_qname(rec));
	      if (read_set.find(seed) == read_set.end()) {
		// Any sequence information?
		if (rec->core.l_qseq > 1) {
		  std::string sequence;
		  sequence.resize(rec->core.l_qseq);
		  uint8_t* seqptr = bam_get_seq(rec);
		  for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
		  std::cout << "Read name: " << bam_get_qname(rec) << ", Length: " << rec->core.l_qseq << std::endl;
		  if (c.trimreads) {
		    int32_t leftPos = -1;
		    int32_t rightPos = -1;
		    bool success = findReadTrim(rec, pos - c.window, pos+c.window, leftPos, rightPos);
		    if (success) rs.push_back(sequence.substr(leftPos, rightPos - leftPos));
		  } else rs.push_back(sequence);
		  read_set.insert(seed);
		} else {
		  missing_reads[seed] = (rec->core.flag & BAM_FREVERSE);
		}
	      }
	    }
	  }
	}
      }
      // Clean-up
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }
    
    // Any missing reads?
    if (!missing_reads.empty()) {
      std::cout << "Secondary alignments" << std::endl;
      int32_t regstart = std::max((int32_t) pos - 100000, 0);
      int32_t regend = std::min((int32_t) pos + 100000, (int32_t) hdr->target_len[refIndex]);
      hts_itr_t* iter = sam_itr_queryi(idx, refIndex, regstart, regend);
      bam1_t* rec = bam_init1();
      while (sam_itr_next(samfile, iter, rec) >= 0) {
	if (rec->core.flag & (BAM_FQCFAIL | BAM_FDUP | BAM_FUNMAP)) continue;
	if (rec->core.qual < c.minMapQual) continue;
	unsigned seed = hash_string(bam_get_qname(rec));
	if (missing_reads.find(seed) != missing_reads.end()) {
	  if (read_set.find(seed) == read_set.end()) {
	    // Any sequence information?
	    if (rec->core.l_qseq > 1) {
	      std::string sequence;
	      sequence.resize(rec->core.l_qseq);
	      uint8_t* seqptr = bam_get_seq(rec);
	      for (int i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];
	      std::cout << "Read name: " << bam_get_qname(rec) << ", Length: " << rec->core.l_qseq << std::endl;
	      
	      // Check alignment direction
	      if ( (rec->core.flag & BAM_FREVERSE) == (missing_reads[seed]) ) {
		if (c.trimreads) {
		  int32_t leftPos = -1;
		  int32_t rightPos = -1;
		  bool success = findReadTrim(rec, pos - c.window, pos+c.window, leftPos, rightPos);
		  if (success) rs.push_back(sequence.substr(leftPos, rightPos - leftPos));
		} else rs.push_back(sequence);
	      } else {
		if (c.trimreads) {
		  // Tricky, todo
		} else {
		  reverseComplement(sequence);
		  rs.push_back(sequence);
		}
	      }
	      read_set.insert(seed);
	    }
	  }
	}
      }
      // Clean-up
      bam_destroy1(rec);
      hts_itr_destroy(iter);
    }
    std::cout << "Number of reads: " << rs.size() << std::endl;
    
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);
    return true;
  }

  
  int consensus(int argc, char **argv) {
    ConfigConsensus c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("called,d", boost::program_options::value<float>(&c.fractionCalled)->default_value(0.5), "fraction of reads required for consensus")
      ("seqtype,t", boost::program_options::value<std::string>(&c.seqtype)->default_value("ill"), "seq. type [ill|ont|pacbio|custom]")
      ;
    
    boost::program_options::options_description bamopt("BAM input options");
    bamopt.add_options()
      ("mapqual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(10), "min. mapping quality")
      ("position,p", boost::program_options::value<std::string>(&c.position)->default_value("chr4:500500"), "position to generate consensus")
      ("window,w", boost::program_options::value<uint32_t>(&c.window), "window around pos that reads need to span")
      ("genome,z", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome [req. for CRAM decoding]")
      ("secondary,s", "consider secondary alignments")
      ("trimreads,r", "trim reads to window")
      ;
    
    boost::program_options::options_description alignment("Alignment scoring options for 'custom' sequencing type");
    alignment.add_options()
      ("gapopen,g", boost::program_options::value<int32_t>(&c.gapopen)->default_value(-10), "gap open")
      ("gapext,e", boost::program_options::value<int32_t>(&c.gapext)->default_value(-1), "gap extension")
      ("match,m", boost::program_options::value<int32_t>(&c.match)->default_value(5), "match")
      ("mismatch,n", boost::program_options::value<int32_t>(&c.mismatch)->default_value(-4), "mismatch")
      ;
    
    boost::program_options::options_description otp("Output options");
    otp.add_options()
      ("outformat,u", boost::program_options::value<std::string>(&c.outformat)->default_value("v"), "output format [v|h]")
      ("alignment,a", boost::program_options::value<boost::filesystem::path>(&c.alignment)->default_value("al.fa.gz"), "vertical/horizontal alignment")
      ("consensus,c", boost::program_options::value<boost::filesystem::path>(&c.consensus)->default_value("cs.fa.gz"), "consensus")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.inputfile), "input bam/fasta file")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);

    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(bamopt).add(alignment).add(otp).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(bamopt).add(alignment).add(otp);
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file"))) {
      std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] <input.bam|reads.fa|reads.fq>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    }
    
    // Secondary alignments
    if (vm.count("secondary")) c.secondary = true;
    else c.secondary = false;

    // Genome file
    if (vm.count("genome")) c.hasGenome = true;
    else c.hasGenome = false;

    // Trim reads?
    if (vm.count("trimreads")) c.trimreads = true;
    else c.trimreads = false;
    
    // Set alignment score
    if (c.seqtype == "ill") {
      c.aliscore = DnaScore<int>(5, -4, -10, -1);
      if (!vm.count("window")) c.window = 5;
    } else if (c.seqtype == "ont") {
      c.aliscore = DnaScore<int>(3, -2, -3, -1);
      if (!vm.count("window")) c.window = 250;
    } else if (c.seqtype == "pacbio") {
      c.aliscore = DnaScore<int>(3, -2, -3, -1);
      if (!vm.count("window")) c.window = 250;
    } else {
      c.aliscore = DnaScore<int32_t>(c.match, c.mismatch, c.gapopen, c.gapext);
      if (!vm.count("window")) c.window = 5;
    }
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    std::cout << "alfred ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;

    // Some status information
    std::cout << "Sequencing type: " << c.seqtype << std::endl;
    std::cout << "Alignment scoring (match: " << c.aliscore.match << ", mismatch: " << c.aliscore.mismatch << ", gapopen: " << c.aliscore.go << ", gapext: " << c.aliscore.ge << ")" << std::endl;
    std::cout << "Window: " << c.window << std::endl;

    // Load reads for consensus
    typedef std::vector<std::string> TReads;
    TReads rs;
    if (vm.count("input-file")) {
      int32_t iftype = inputType(c.inputfile.string());
      if (iftype == 0) { // BAM
	bool rtval = _loadBamReads(c, rs);
	if (!rtval) return 1;
      } else if ((iftype == 1) || (iftype == 2)) { // FASTA or FASTQ
	_loadFastaReads(c, rs);
      } else {
	std::cerr << "Unknown input file format! " << c.inputfile.string() << std::endl;
	return 1;
      }
    }

    // Any reads?
    if (rs.empty()) {
      std::cerr << "No reads for consensus found!" << std::endl;
      return 1;
    }

    // Generate Consensus
    std::string consensus;
    msa(c, rs, consensus);
    
    // Output consensus
    boost::iostreams::filtering_ostream rcfile;
    rcfile.push(boost::iostreams::gzip_compressor());
    rcfile.push(boost::iostreams::file_sink(c.consensus.c_str(), std::ios_base::out | std::ios_base::binary));
    rcfile << ">Consensus" << std::endl;
    rcfile << consensus << std::endl;
    rcfile.pop();
    
    // Done
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
  
    return 0;
  }


}

#endif
