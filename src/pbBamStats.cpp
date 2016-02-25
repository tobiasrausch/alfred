
#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <fstream>
#include <zlib.h>

#define BOOST_DISABLE_ASSERTS
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <htslib/kseq.h>
#include <htslib/sam.h>

#ifdef PROFILE
#include "gperftools/profiler.h"
#endif


KSEQ_INIT(gzFile, gzread)

struct Config {
  bool hasRegionFile;
  boost::filesystem::path referenceFile;
  boost::filesystem::path regionFile;
  boost::filesystem::path bamFile;
};


struct Interval {
  int32_t start;
  int32_t end;
  
  Interval(int32_t s, int32_t e) : start(s), end(e) {}
};


inline uint32_t
lastAlignedPosition(bam1_t const* rec) {
  uint32_t* cigar = bam_get_cigar(rec);
  uint32_t alen = 0;
  for (std::size_t i = 0; i < rec->core.n_cigar; ++i)
    if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CDEL)) alen += bam_cigar_oplen(cigar[i]);
  return rec->core.pos + alen;
}

int main(int argc, char **argv) {

#ifdef PROFILE
  ProfilerStart("pbBamStats.prof");
#endif

  Config c;

  // Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("reference,r", boost::program_options::value<boost::filesystem::path>(&c.referenceFile), "reference fasta file (required)")
    ("bed,b", boost::program_options::value<boost::filesystem::path>(&c.regionFile), "BED file with regions to analyze (optional)")
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
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);

  // Check command line arguments
  if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("reference"))) {
    std::cout << "Usage: " << argv[0] << " [OPTIONS] -r <ref.fa> <pacbio.bam>" << std::endl;
    std::cout << visible_options << "\n";
    return 1;
  } 

  // Check reference
  if (vm.count("reference")) {
    if (!(boost::filesystem::exists(c.referenceFile) && boost::filesystem::is_regular_file(c.referenceFile) && boost::filesystem::file_size(c.referenceFile))) {
      std::cerr << "Input reference file is missing: " << c.referenceFile.string() << std::endl;
      return 1;
    }
  }
  
  // Check region file
  if (vm.count("bed")) {
    if (!(boost::filesystem::exists(c.regionFile) && boost::filesystem::is_regular_file(c.regionFile) && boost::filesystem::file_size(c.regionFile))) {
      std::cerr << "Input BED region file is missing: " << c.regionFile.string() << std::endl;
      return 1;
    }
    c.hasRegionFile = true;
  } else c.hasRegionFile = false;

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  // Load bam file
  samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
  if (samfile == NULL) {
    std::cerr << "Fail to open file " << c.bamFile.string() << std::endl;
    return 1;
  }

  // Load bam index
  hts_idx_t* idx = sam_index_load(samfile, c.bamFile.string().c_str());
  if (idx == NULL) {
    std::cerr << "Fail to open index for " << c.bamFile.string() << std::endl;
    return 1;
  }

  // Load bam header
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

  // Counters
  uint32_t matchCount = 0;
  uint32_t mismatchCount = 0;
  uint32_t delCount = 0;
  uint32_t insCount = 0;
  uint32_t softClipCount = 0;
  uint32_t hardClipCount = 0;

  // Parse reference and BAM file
  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "BAM file parsing" << std::endl;
  boost::progress_display show_progress( hdr->n_targets );
  kseq_t *seq;
  int l;
  gzFile fp = gzopen(c.referenceFile.string().c_str(), "r");
  seq = kseq_init(fp);
  while ((l = kseq_read(seq)) >= 0) {
    // Find BAM reference index for the given sequence name from the FASTA file
    for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
      if (std::string(seq->name.s) == std::string(hdr->target_name[refIndex])) {
	++show_progress;

	// Parse all alignments in regions of that chromosome
	for(TChromosomeRegions::const_iterator itC = gRegions[refIndex].begin(); itC != gRegions[refIndex].end(); ++itC) {
	  hts_itr_t* iter = sam_itr_queryi(idx, refIndex, itC->start, itC->end);
	  bam1_t* rec = bam_init1();
	  while (sam_itr_next(samfile, iter, rec) >= 0) {
	    // We may want to count these (secondary alignments, duplicates, supplementary alignments)
	    if (rec->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;

	    // Get the read sequence
	    std::string sequence;
	    sequence.resize(rec->core.l_qseq);
	    uint8_t* seqptr = bam_get_seq(rec);
	    for (int32_t i = 0; i < rec->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

	    // Get the reference slice
	    std::string refslice = boost::to_upper_copy(std::string(seq->seq.s + rec->core.pos, seq->seq.s + lastAlignedPosition(rec)));
	    
	    // Debug 
	    //std::cout << matchCount << ',' << mismatchCount << ',' << delCount << ',' << insCount << ',' << softClipCount << ',' << hardClipCount << std::endl;
	    //std::cout << refslice << std::endl;
	    //std::cout << sequence << std::endl;

	    uint32_t rp = 0; // reference pointer
	    uint32_t sp = 0; // sequence pointer

	    // Finally the fun part, parse the CIGAR
	    uint32_t* cigar = bam_get_cigar(rec);
	    for (std::size_t i = 0; i < rec->core.n_cigar; ++i)
	      if (bam_cigar_op(cigar[i]) == BAM_CMATCH) {
		// match or mismatch
		for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]);++k) {
		  if (sequence[sp] == refslice[rp]) ++matchCount;
		  else ++mismatchCount;
		  ++sp;
		  ++rp;
		}
	      } else if (bam_cigar_op(cigar[i]) == BAM_CDEL) {
		delCount += bam_cigar_oplen(cigar[i]);
		rp += bam_cigar_oplen(cigar[i]);
	      } else if (bam_cigar_op(cigar[i]) == BAM_CINS) {
		insCount += bam_cigar_oplen(cigar[i]);
		sp += bam_cigar_oplen(cigar[i]);
	      } else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) {
		softClipCount += bam_cigar_oplen(cigar[i]);
		sp += bam_cigar_oplen(cigar[i]);
	      } else if(bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
		hardClipCount += bam_cigar_oplen(cigar[i]);
	      } else {
		std::cerr << "Unknown Cigar options" << std::endl;
		return 1;
	      }
	  }
	  // clean-up
	  bam_destroy1(rec);
	  hts_itr_destroy(iter);
	}
      }        
    }
  }
  // clea-up for reference
  kseq_destroy(seq);
  gzclose(fp);

  // clean-up for bam
  bam_hdr_destroy(hdr);
  hts_idx_destroy(idx);
  sam_close(samfile);

  now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;

  // Output statistics

#ifdef PROFILE
  ProfilerStop();
#endif

  // Counts
  std::cout << "Counts" << std::endl;
  std::cout << "#Matched bases=" << matchCount << std::endl;
  std::cout << "#Mismatched bases=" << mismatchCount << std::endl;
  std::cout << "#Deleted bases=" << delCount << std::endl;
  std::cout << "#Inserted bases=" << insCount << std::endl;
  std::cout << "#SoftClipped bases=" << softClipCount << std::endl;
  std::cout << "#HardClipped bases=" << hardClipCount << std::endl;


  return 0;
}
