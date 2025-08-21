#ifndef QC_H
#define QC_H

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

#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "bamstats.h"
#include "util.h"
#include "version.h"

namespace bamstats
{

struct ConfigQC {
  bool hasRegionFile;
  bool hasJsonOut;
  bool hasTsvOut;
  bool ignoreRG;
  bool singleRG;
  bool isHaplotagged;
  bool isMitagged;
  bool secondary;
  bool supplementary;
  bool meanCoverage;
  float nXChrLen;
  uint32_t minChrLen;
  std::string rgname;
  std::string sampleName;
  boost::filesystem::path jsonout;
  boost::filesystem::path outfile;
  boost::filesystem::path genome;
  boost::filesystem::path regionFile;
  boost::filesystem::path bamFile;
};


int qc(int argc, char **argv) {
  ConfigQC c;
  c.isHaplotagged = false;
  c.isMitagged = false;
  c.minChrLen = 10000000;
  std::string sampleName;
  
  // Parameter
  boost::program_options::options_description generic("Generic options");
  generic.add_options()
    ("help,?", "show help message")
    ("reference,r", boost::program_options::value<boost::filesystem::path>(&c.genome), "reference fasta file (required)")
    ("bed,b", boost::program_options::value<boost::filesystem::path>(&c.regionFile), "bed file with target regions (optional)")
    ("name,a", boost::program_options::value<std::string>(&sampleName), "sample name (optional, otherwise SM tag is used)")
    ("jsonout,j", boost::program_options::value<boost::filesystem::path>(&c.jsonout), "gzipped json output file")
    ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile), "gzipped tsv output file")
    ("meancov,m", "report mean coverage as float instead of median integer coverage")
    ("secondary,s", "evaluate secondary alignments")
    ("supplementary,u", "evaluate supplementary alignments") 
    ;

  boost::program_options::options_description rgopt("Read-group options");
  rgopt.add_options()
    ("rg,g", boost::program_options::value<std::string>(&c.rgname), "only analyze this read group (optional)")
    ("ignore,i", "ignore read-groups")
    ;

  boost::program_options::options_description hidden("Hidden options");
  hidden.add_options()
    ("nxchrlen,n", boost::program_options::value<float>(&c.nXChrLen)->default_value(0.95), "N95 chromosome length to trim mapping table [0,1]")
    ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bamFile), "input bam file")
    ;

  boost::program_options::positional_options_description pos_args;
  pos_args.add("input-file", -1);

  boost::program_options::options_description cmdline_options;
  cmdline_options.add(generic).add(rgopt).add(hidden);
  boost::program_options::options_description visible_options;
  visible_options.add(generic).add(rgopt);

  // Parse command-line
  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
  boost::program_options::notify(vm);

  // Check command line arguments
  if ((vm.count("help")) || ((!vm.count("outfile")) && (!vm.count("jsonout"))) || (!vm.count("input-file")) || (!vm.count("reference"))) {
    std::cout << std::endl;
    std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] -r <ref.fa> -j <qc.json.gz> <aligned.bam>" << std::endl;
    std::cout << visible_options << "\n";
    return 1;
  }

  // Tsv output file
  if (vm.count("outfile")) c.hasTsvOut = true;
  else c.hasTsvOut = false;

  // Json output file
  if (vm.count("jsonout")) c.hasJsonOut = true;
  else c.hasJsonOut = false;

  // Secondary alignments
  if (vm.count("secondary")) c.secondary = true;
  else c.secondary = false;

  // Mean coverage
  if (vm.count("meancov")) c.meanCoverage = true;
  else c.meanCoverage = false;
  
  // Supplementary alignments
  if (vm.count("supplementary")) c.supplementary = true;
  else c.supplementary = false;
  
  // Check N95
  if (c.nXChrLen > 1) c.nXChrLen = 1;
  else if (c.nXChrLen < 0) c.nXChrLen = 0;

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

  // Check bam file
  if (!(boost::filesystem::exists(c.bamFile) && boost::filesystem::is_regular_file(c.bamFile) && boost::filesystem::file_size(c.bamFile))) {
    std::cerr << "Alignment file is missing: " << c.bamFile.string() << std::endl;
    return 1;
  }
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
  faidx_t* fai = fai_load(c.genome.string().c_str());
  for(int32_t refIndex=0; refIndex < hdr->n_targets; ++refIndex) {
    std::string tname(hdr->target_name[refIndex]);
    if (!faidx_has_seq(fai, tname.c_str())) {
      std::cerr << "BAM file chromosome " << hdr->target_name[refIndex] << " is NOT present in your reference file " << c.genome.string() << std::endl;
      return 1;
    }
  }
  fai_destroy(fai);
  if (!vm.count("name")) {
    if (!getSMTag(std::string(hdr->text), c.bamFile.stem().string(), sampleName)) {
      c.sampleName = c.bamFile.stem().string();
      std::cerr << "Warning: Multiple read-groups (@RG). Sample name is set to " << c.sampleName << std::endl;
    } else c.sampleName = sampleName;
  } else {
    if (sampleName.size()) c.sampleName = sampleName;
    else {
      std::cerr << "Empty sample name is not allowed!" << std::endl;
      return 1;
    }
  }
  bam_hdr_destroy(hdr);
  hts_idx_destroy(idx);
  sam_close(samfile);

  // Ignore read groups
  c.ignoreRG = true;
  c.singleRG = false;
  if (!vm.count("ignore")) {
    // Single RG
    if (vm.count("rg")) {
      c.ignoreRG = false;
      c.singleRG = true;
    } else {
      // Check number of RGs
      int32_t rgc = countRGs(c);
      if (rgc > 3) {
	std::cerr << "Warning: BAM file has more than 3 RGs. Please run read-groups separately!" << std::endl;
	std::cerr << "Warning: Program continues but ignores read-groups." << std::endl;
      } else {
	c.ignoreRG = false;
	c.singleRG = false;
      }
    }
  }

  
  // Check region file
  if (vm.count("bed")) {
    if (!(boost::filesystem::exists(c.regionFile) && boost::filesystem::is_regular_file(c.regionFile) && boost::filesystem::file_size(c.regionFile))) {
      std::cerr << "Input region file in bed format is missing: " << c.regionFile.string() << std::endl;
      return 1;
    }
    std::string oldChr;
    faidx_t* fai = fai_load(c.genome.string().c_str());
    if (is_gz(c.regionFile)) {
      std::ifstream file(c.regionFile.string().c_str(), std::ios_base::in | std::ios_base::binary);
      boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
      dataIn.push(boost::iostreams::gzip_decompressor());
      dataIn.push(file);
      std::istream instream(&dataIn);
      std::string intervalLine;
      while(std::getline(instream, intervalLine)) {
	typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	boost::char_separator<char> sep(" \t,;");
	Tokenizer tokens(intervalLine, sep);
	Tokenizer::iterator tokIter = tokens.begin();
	if (tokIter!=tokens.end()) {
	  std::string chrName=*tokIter++;
	  if (chrName.compare(oldChr) != 0) {
	    oldChr = chrName;
	    if (!faidx_has_seq(fai, chrName.c_str())) {
	      std::cerr << "Chromosome from bed file " << chrName << " is NOT present in your reference file " << c.genome.string() << std::endl;
	      return 1;
	    }
	  }
	}
      }
      dataIn.pop();
    } else {
      std::ifstream interval_file(c.regionFile.string().c_str(), std::ifstream::in);
      if (interval_file.is_open()) {
	while (interval_file.good()) {
	  std::string intervalLine;
	  getline(interval_file, intervalLine);
	  typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
	  boost::char_separator<char> sep(" \t,;");
	  Tokenizer tokens(intervalLine, sep);
	  Tokenizer::iterator tokIter = tokens.begin();
	  if (tokIter!=tokens.end()) {
	    std::string chrName=*tokIter++;
	    if (chrName.compare(oldChr) != 0) {
	      oldChr = chrName;
	      if (!faidx_has_seq(fai, chrName.c_str())) {
		std::cerr << "Chromosome from bed file " << chrName << " is NOT present in your reference file " << c.genome.string() << std::endl;
		return 1;
	      }
	    }
	  }
	}
	interval_file.close();
      }
    }
    fai_destroy(fai);
    c.hasRegionFile = true;
  } else c.hasRegionFile = false;

  // Show cmd
  boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
  std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
  std::cout << "alfred ";
  for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
  std::cout << std::endl;

  return bamStatsRun(c); 
}

}

#endif
