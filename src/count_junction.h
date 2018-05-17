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

#ifndef COUNT_JUNCTION_H
#define COUNT_JUNCTION_H

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

  struct CountJunctionConfig {
    unsigned short minQual;
    uint8_t inputFileFormat;   // 0 = gtf, 1 = bed, 2 = gff3
    std::map<std::string, int32_t> nchr;
    std::string sampleName;
    std::string idname;
    std::string feature;
    boost::filesystem::path gtfFile;
    boost::filesystem::path bedFile;
    boost::filesystem::path bamFile;
    boost::filesystem::path outintra;
    boost::filesystem::path outinter;
  };

  template<typename TConfig>
  inline int32_t
  countJunctionRun(TConfig const& c) {

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
    if (c.inputFileFormat == 0) tf = parseGTFAll(c, gRegions, geneIds);
    else if (c.inputFileFormat == 1) tf = parseBEDAll(c, gRegions, geneIds);
    else if (c.inputFileFormat == 2) tf = parseGFF3All(c, gRegions, geneIds);
    if (tf == 0) {
      std::cerr << "Error parsing GTF/GFF3/BED file!" << std::endl;
      return 1;
    }

    /*
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
    */
    
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    
#ifdef PROFILE
    ProfilerStop();
#endif

    return 0;
  }


  int count_junction(int argc, char **argv) {
    CountJunctionConfig c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("map-qual,m", boost::program_options::value<unsigned short>(&c.minQual)->default_value(10), "min. mapping quality")
      ("outintra,o", boost::program_options::value<boost::filesystem::path>(&c.outintra)->default_value("intra.tsv"), "intra-gene exon-exon junction reads")
      ("outinter,p", boost::program_options::value<boost::filesystem::path>(&c.outinter)->default_value("inter.tsv"), "inter-gene exon-exon junction reads")
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
      std::cout << visible_options << "\n";
      return 1;
    }

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

    return countJunctionRun(c);
  }
  


  
}

#endif
