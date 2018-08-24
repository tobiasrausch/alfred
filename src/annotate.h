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

#ifndef ANNOTATE_H
#define ANNOTATE_H

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

  struct AnnotateConfig {
    uint8_t inputFileFormat;   // 0 = gtf, 1 = bed, 2 = gff3
    int32_t maxDistance;
    std::map<std::string, int32_t> nchr;
    std::string idname;
    std::string feature;
    boost::filesystem::path gtfFile;
    boost::filesystem::path bedFile;
    boost::filesystem::path infile;
    boost::filesystem::path outfile;
  };


  template<typename TConfig, typename TGenomicRegions, typename TGeneIds>
  inline int32_t
  bed_anno(TConfig const& c, TGenomicRegions& gRegions, TGeneIds& geneIds) {
    typedef typename TGenomicRegions::value_type TChromosomeRegions;

    // Parse BED file
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "BED file parsing" << std::endl;
    boost::progress_display show_progress(c.nchr.size());

    // Distance vector
    std::vector<int32_t> dist(geneIds.size(), 0);
    
    // Open output file
    std::ofstream ofile(c.outfile.string().c_str());
    ofile << "chrom\tstart\tend\tid\tfeature\tdistance" << std::endl;
    
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

      // Annotate intervals
      std::ifstream chrFile(c.infile.string().c_str(), std::ifstream::in);
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
	    std::string name = "NA";
	    if (tokIter != tokens.end()) name = *tokIter++;
	    if (start >= end) continue;  // Bed has right-open intervals
	    typedef std::vector<int32_t> TFeaturePos;
	    TFeaturePos featurepos;
	    int32_t realstart = std::max(0, start - c.maxDistance);
	    int32_t realend = std::min(250000000, end + c.maxDistance);
	    for(int32_t i = realstart; i<realend; ++i)
	      if (featureBitMap[i]) featurepos.push_back(i);

	    // Find feature
	    typedef std::set<int32_t> TFeatureIds;
	    TFeatureIds featureid;  // No feature by default
	    if (!featurepos.empty()) {
	      int32_t fpfirst = featurepos[0];
	      int32_t fplast = featurepos[featurepos.size()-1];
	      for(typename TChromosomeRegions::const_iterator vIt = gRegions[refIndex].begin(); vIt != gRegions[refIndex].end(); ++vIt) {
		if (vIt->end <= fpfirst) continue;
		if (vIt->start > fplast) break; // Sorted intervals so we can stop searching
		for(TFeaturePos::const_iterator fIt = featurepos.begin(); fIt != featurepos.end(); ++fIt) {
		  if ((vIt->start <= *fIt) && (vIt->end > *fIt)) {
		    featureid.insert(vIt->lid);

		    // Get distance
		    int32_t locdist = 0;
		    if (vIt->end < start) { locdist = vIt->end - start; }
		    if (end < vIt->start) { locdist = vIt->start - end; }
		    dist[vIt->lid] = locdist;
		    break;
		  }
		}
	      }
	    }

	    // Output overlapping features
	    ofile << chrName << "\t" << start << "\t" << end << "\t" << name << "\t";
	    // Feature names
	    if (featureid.empty()) {
	      ofile << "NA";
	    } else {
	      bool firstF = true;
	      for(typename TFeatureIds::const_iterator itF = featureid.begin(); itF != featureid.end(); ++itF) {
		if (!firstF) ofile << ',';
		else firstF = false;
		ofile << geneIds[*itF];
	      }
	    }
	    ofile << "\t";
	    // Feature distances
	    if (featureid.empty()) {
	      ofile << "NA";
	    } else {
	      bool firstF = true;
	      for(typename TFeatureIds::const_iterator itF = featureid.begin(); itF != featureid.end(); ++itF) {
		if (!firstF) ofile << ',';
		else firstF = false;
		ofile << dist[*itF];
	      }
	    }
	    ofile << std::endl;
	  }
	}
	chrFile.close();
      }
    }
    ofile.close();
    
    return 0;
  }


  
  template<typename TConfig>
  inline int32_t
  annotateRun(TConfig const& c) {

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
      std::cerr << "Please check that the chromosome names agree (chr1 versus 1) between input and annotation file." << std::endl;
      return 1;
    }

    // Feature annotation
    int32_t retparse = bed_anno(c, gRegions, geneIds);
    if (retparse != 0) {
      std::cerr << "Error in BED annotation!" << std::endl;
      return 1;
    }
    
    // Done
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    
#ifdef PROFILE
    ProfilerStop();
#endif

    return 0;
  }


  int annotate(int argc, char **argv) {
    AnnotateConfig c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("distance,d", boost::program_options::value<int32_t>(&c.maxDistance)->default_value(0), "max. distance (0: overlapping features only)")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("anno.bed"), "output file")
      ;

    boost::program_options::options_description gtfopt("GTF/GFF3 annotation file options");
    gtfopt.add_options()
      ("gtf,g", boost::program_options::value<boost::filesystem::path>(&c.gtfFile), "gtf/gff3 file")
      ("id,i", boost::program_options::value<std::string>(&c.idname)->default_value("gene_name"), "gtf/gff3 attribute")
      ("feature,f", boost::program_options::value<std::string>(&c.feature)->default_value("gene"), "gtf/gff3 feature")
      ;
    
    boost::program_options::options_description bedopt("BED annotation file options, columns chr, start, end, name");
    bedopt.add_options()
      ("bed,b", boost::program_options::value<boost::filesystem::path>(&c.bedFile), "bed file")
      ;

    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.infile), "input file")
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
      std::cout << "Usage: alfred annotate " << argv[0] << " [OPTIONS] -g <hg19.gtf.gz> <peaks.bed>" << std::endl;
      std::cout << "Usage: alfred annotate " << argv[0] << " [OPTIONS] -b <hg19.bed.gz> <peaks.bed>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    }

    // Input BED file
    if (!(boost::filesystem::exists(c.infile) && boost::filesystem::is_regular_file(c.infile) && boost::filesystem::file_size(c.infile))) {
      std::cerr << "Input BED file is missing." << std::endl;
      return 1;
    } else {
      std::string oldChr = "";
      typedef std::set<std::string> TChrSet;
      TChrSet chrSet;
      std::ifstream chrFile(c.infile.string().c_str(), std::ifstream::in);
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

    return annotateRun(c);
  }
  


  
}

#endif
