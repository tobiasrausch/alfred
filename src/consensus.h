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
#include "anchored.h"
#include "consedlib.h"

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
    float maxdiv;
    float maxclip;
    int32_t insWindow;
    int32_t chimDist;
    int32_t minInsForClip;
    int32_t minClipSupport;
    int32_t minClip;
    int32_t clipCap;
    int32_t nAlleles;
    int32_t minAlleleSupport;
    int32_t alleleWindow;
    float alleleDist;
    std::string position;
    std::string seqtype;
    std::string consmode;
    std::string outformat;
    DnaScore<int32_t> aliscore;
    boost::filesystem::path alignment;
    boost::filesystem::path genome;
    boost::filesystem::path consensus;
    boost::filesystem::path inputfile;
  };


  template<typename TConfig>
  inline void
  _loadFastaReads(TConfig const& c, std::vector<AlignedRead>& areads) {
    // Load FASTA/FASTQ
    areads.clear();
    gzFile fp = gzopen(c.inputfile.string().c_str(), "r");
    if (fp) {
      kseq_t* kseq = kseq_init(fp);
      int l;
      while((l = kseq_read(kseq)) >= 0) {
	//kseq->qual.s are the qualities
	areads.push_back(AlignedRead(std::string(kseq->seq.s)));
	areads.back().name = kseq->name.s;
	std::cout << "Read name: " << kseq->name.s << ", Length: " << areads.back().seq.size() << std::endl;
      }
      kseq_destroy(kseq);
      gzclose(fp);
    }
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
 
  // supplementary alignments
  inline bool
  _hasDistantSA(bam1_t const* rec, bam_hdr_t const* hdr, int32_t chimDist) {
    uint8_t* saTag = bam_aux_get(rec, "SA");
    if (saTag == NULL) return false;
    char* saStr = bam_aux2Z(saTag);
    if (saStr == NULL) return false;
    std::string primChrom = hdr->target_name[rec->core.tid];
    int32_t primPos = rec->core.pos;
    std::string s(saStr);
    std::size_t start = 0;
    while (start < s.size()) {
      std::size_t semi = s.find(';', start);
      std::string entry = (semi == std::string::npos) ? s.substr(start) : s.substr(start, semi - start);
      if (!entry.empty()) {
	std::size_t c1 = entry.find(',');
	if (c1 != std::string::npos) {
	  std::string rname = entry.substr(0, c1);
	  std::size_t c2 = entry.find(',', c1 + 1);
	  int32_t saPos = std::atoi(entry.substr(c1 + 1, (c2 == std::string::npos) ? std::string::npos : c2 - (c1 + 1)).c_str());
	  if (rname != primChrom) return true;
	  if (std::abs(saPos - primPos) > chimDist) return true;
	}
      }
      if (semi == std::string::npos) break;
      start = semi + 1;
    }
    return false;
  }

  inline void
  _pushAligned(std::vector<AlignedRead>& areads, bam1_t const* rec, bam_hdr_t const* hdr, int32_t chimDist, std::string const& seq, bool anchored) {
    AlignedRead ar;
    ar.seq = seq;
    ar.name = bam_get_qname(rec);
    ar.mapq = rec->core.qual;
    ar.rev = (rec->core.flag & BAM_FREVERSE);
    ar.distantSA = _hasDistantSA(rec, hdr, chimDist);
    if (anchored) {
      ar.anchored = true;
      ar.refStart = rec->core.pos;
      uint32_t* cig = bam_get_cigar(rec);
      ar.cigar.assign(cig, cig + rec->core.n_cigar);
      if (rec->core.n_cigar) {
	if (bam_cigar_op(cig[0]) == BAM_CSOFT_CLIP) ar.leftClip = bam_cigar_oplen(cig[0]);
	if (bam_cigar_op(cig[rec->core.n_cigar - 1]) == BAM_CSOFT_CLIP) ar.rightClip = bam_cigar_oplen(cig[rec->core.n_cigar - 1]);
      }
    }
    areads.push_back(ar);
  }

  template<typename TConfig>
  inline bool
  _loadBamReads(TConfig const& c, std::vector<AlignedRead>& areads) {
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
    const int32_t chimDist = c.chimDist;
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
		    if (success) _pushAligned(areads, rec, hdr, chimDist, sequence.substr(leftPos, rightPos - leftPos), false);
		  } else _pushAligned(areads, rec, hdr, chimDist, sequence, true);
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
		  if (success) _pushAligned(areads, rec, hdr, chimDist, sequence.substr(leftPos, rightPos - leftPos), false);
		} else _pushAligned(areads, rec, hdr, chimDist, sequence, true);
	      } else {
		if (c.trimreads) {
		  // Nop
		} else {
		  reverseComplement(sequence);
		  _pushAligned(areads, rec, hdr, chimDist, sequence, false);
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
    std::cout << "Number of reads: " << areads.size() << std::endl;
    
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
      ("mode,b", boost::program_options::value<std::string>(&c.consmode)->default_value("dp"), "msa algorithm [dp|ed|anchored], dp: dynamic programming (slow), ed: edit distance (fast), anchored: reference-anchored using BAM")
      ("called,d", boost::program_options::value<float>(&c.fractionCalled)->default_value(0.5), "fraction of reads required for consensus")
      ("maxdiv,x", boost::program_options::value<float>(&c.maxdiv)->default_value(0.25), "max per-read divergence (anchored mode)")
      ("maxclip,y", boost::program_options::value<float>(&c.maxclip)->default_value(0.3), "max soft-clip fraction with SA (anchored mode)")
      ("seqtype,t", boost::program_options::value<std::string>(&c.seqtype)->default_value("ill"), "seq. type [ill|ont|pacbio|custom]")
      ;
    
    boost::program_options::options_description bamopt("BAM input options");
    bamopt.add_options()
      ("mapqual,q", boost::program_options::value<uint16_t>(&c.minMapQual)->default_value(10), "min. mapping quality")
      ("position,p", boost::program_options::value<std::string>(&c.position)->default_value("chr4:500500"), "position to generate consensus")
      ("window,w", boost::program_options::value<uint32_t>(&c.window), "window around pos that reads need to span")
      ("genome,z", boost::program_options::value<boost::filesystem::path>(&c.genome), "genome")
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
    
    boost::program_options::options_description anchored("Anchored mode options");
    anchored.add_options()
      ("inswindow", boost::program_options::value<int32_t>(&c.insWindow)->default_value(50), "window to cluster insertions/soft-clips")
      ("mininsclip", boost::program_options::value<int32_t>(&c.minInsForClip)->default_value(10), "min. insertion length")
      ("minclip", boost::program_options::value<int32_t>(&c.minClip)->default_value(20), "min. soft-clip length")
      ("clipsupport", boost::program_options::value<int32_t>(&c.minClipSupport)->default_value(2), "min. clipped reads")
      ("clipcap", boost::program_options::value<int32_t>(&c.clipCap)->default_value(8000), "max. soft-clip length")
      ("chimdist", boost::program_options::value<int32_t>(&c.chimDist)->default_value(50000), "chimera SA distance")
      ("alleles", boost::program_options::value<int32_t>(&c.nAlleles)->default_value(0), "max. number of alleles (0=auto)")
      ("allelewindow", boost::program_options::value<int32_t>(&c.alleleWindow)->default_value(200), "window for local allele calling")
      ("allelesupport", boost::program_options::value<int32_t>(&c.minAlleleSupport)->default_value(2), "min. reads per allele")
      ("alleledist", boost::program_options::value<float>(&c.alleleDist)->default_value(0.2), "max. relative insertion-length difference")
      ;

    boost::program_options::options_description otp("Output options");
    otp.add_options()
      ("outformat,u", boost::program_options::value<std::string>(&c.outformat)->default_value("h"), "output format [v|h]")
      ("alignment,a", boost::program_options::value<boost::filesystem::path>(&c.alignment)->default_value("al.fa.gz"), "vertical/horizontal alignment (per-allele al.<k>.fa.gz in allele mode; -u h rows are read names)")
      ("consensus,c", boost::program_options::value<boost::filesystem::path>(&c.consensus)->default_value("cs.fa.gz"), "consensus")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.inputfile), "input bam/fasta file")
      ;
    
    boost::program_options::positional_options_description pos_args;
    pos_args.add("input-file", -1);

    boost::program_options::options_description cmdline_options;
    cmdline_options.add(generic).add(bamopt).add(alignment).add(anchored).add(otp).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(bamopt).add(alignment).add(anchored).add(otp);
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

    // Load reads for consensus (keeps reference anchoring from the BAM)
    std::vector<AlignedRead> areads;
    if (vm.count("input-file")) {
      int32_t iftype = inputType(c.inputfile.string());
      if (iftype == 0) { // BAM
	bool rtval = _loadBamReads(c, areads);
	if (!rtval) return 1;
      } else if ((iftype == 1) || (iftype == 2)) { // FASTA or FASTQ
	_loadFastaReads(c, areads);
      } else {
	std::cerr << "Unknown input file format! " << c.inputfile.string() << std::endl;
	return 1;
      }
    }

    // Any reads?
    if (areads.empty()) {
      std::cerr << "No reads for consensus found!" << std::endl;
      return 1;
    }

    // Plain sequences for the global consensus modes
    std::vector<std::string> rs;
    readSequences(areads, rs);

    // De novo allele splitting (anchored mode only)
    bool alleleMode = ((c.consmode == "anchored") && (c.nAlleles != 1));

    // Generate Consensus
    std::string consensus;
    std::vector<AlleleConsensus> alleles;
    if (alleleMode) msaAlleles(c, areads, alleles);
    else if (c.consmode == "dp") msa(c, rs, consensus);
    else if (c.consmode == "anchored") msaAnchored(c, areads, consensus);
    else msaEdlib(c, rs, consensus);

    // Output consensus
    boost::iostreams::filtering_ostream rcfile;
    rcfile.push(boost::iostreams::gzip_compressor());
    rcfile.push(boost::iostreams::file_sink(c.consensus.string(), std::ios_base::out | std::ios_base::binary));
    if (alleleMode) {
      for(uint32_t i = 0; i < alleles.size(); ++i) rcfile << '>' << alleles[i].name << std::endl << alleles[i].cons << std::endl;
    } else {
      rcfile << ">Consensus" << std::endl;
      rcfile << consensus << std::endl;
    }
    rcfile.pop();
    
    // Done
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
  
    return 0;
  }


}

#endif
