#ifndef ANNOTATE_H
#define ANNOTATE_H

#include <limits>

#include <boost/icl/split_interval_map.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "version.h"
#include "util.h"
#include "gtf.h"
#include "gff3.h"
#include "bed.h"
#include "motif.h"

namespace bamstats
{

  struct AnnotateConfig {
    typedef std::map<std::string, int32_t> TChrMap;
    bool motifPosOut;
    bool nearest;
    bool overlappingHits;
    uint8_t inputFileFormat;   // 0 = gtf, 1 = bed, 2 = gff3, 3 = motif file
    int32_t maxDistance;
    float motifScoreQuantile;
    TChrMap nchr;
    std::string idname;
    std::string feature;
    boost::filesystem::path motifFile;
    boost::filesystem::path genome;
    boost::filesystem::path gtfFile;
    boost::filesystem::path bedFile;
    boost::filesystem::path infile;
    boost::filesystem::path outpos;
    boost::filesystem::path outgene;
    boost::filesystem::path outfile;
  };


  template<typename TConfig, typename TGenomicRegions, typename TGeneIds>
  inline int32_t
  bed_anno(TConfig const& c, TGenomicRegions& gRegions, TGeneIds& geneIds) {
    typedef typename TGenomicRegions::value_type TChromosomeRegions;

    // Parse BED file
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "BED file parsing" << std::endl;

    // Distance vector
    std::vector<int32_t> dist(geneIds.size(), 0);
    
    // Open output file
    std::ofstream ofile(c.outfile.string().c_str());
    ofile << "chrom\tstart\tend\tid\tfeature\tdistance" << std::endl;

    // Peak count and names
    std::vector<std::string> peakNames;

    // Gene-level summary
    typedef std::pair<int32_t, int32_t> TDistPeak;
    typedef std::vector<TDistPeak> TPeaksPerGene;
    std::vector<TPeaksPerGene> geneView(geneIds.size(), TPeaksPerGene());
    
    // Iterate chromosomese
    for(int32_t refIndex=0; refIndex < (int32_t) c.nchr.size(); ++refIndex) {

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
	    {
	      std::string name = "Interval" + boost::lexical_cast<std::string>(peakNames.size());
	      if (tokIter != tokens.end()) name = *tokIter++;
	      peakNames.push_back(name);
	    }
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
	    int32_t bestDistance = 250000000;
	    int32_t peakId = peakNames.size() - 1;
	    ofile << chrName << "\t" << start << "\t" << end << "\t" << peakNames[peakId] << "\t";
	    // Feature names
	    if (featureid.empty()) {
	      ofile << "NA";
	    } else {
	      for(typename TFeatureIds::const_iterator itF = featureid.begin(); itF != featureid.end(); ++itF) geneView[*itF].push_back(std::make_pair(dist[*itF], peakId));
	      if (c.nearest) {
		int32_t bestFeature = 0;
		for(typename TFeatureIds::const_iterator itF = featureid.begin(); itF != featureid.end(); ++itF) {
		  if (std::abs(dist[*itF]) < std::abs(bestDistance)) {
		    bestDistance = dist[*itF];
		    bestFeature = *itF;
		  }
		}
		ofile << geneIds[bestFeature];
	      } else {
		bool firstF = true;
		for(typename TFeatureIds::const_iterator itF = featureid.begin(); itF != featureid.end(); ++itF) {
		  if (!firstF) ofile << ',';
		  else firstF = false;
		  ofile << geneIds[*itF];
		}
	      }
	    }
	    ofile << "\t";
	    // Feature distances
	    if (featureid.empty()) {
	      ofile << "NA";
	    } else {
	      if (c.nearest) {
		ofile << bestDistance;
	      } else {
		bool firstF = true;
		for(typename TFeatureIds::const_iterator itF = featureid.begin(); itF != featureid.end(); ++itF) {
		  if (!firstF) ofile << ',';
		  else firstF = false;
		  ofile << dist[*itF];
		}
	      }
	    }
	    ofile << std::endl;
	  }
	}
	chrFile.close();
      }
    }
    ofile.close();

    // Output gene-level view
    std::ofstream gfile(c.outgene.string().c_str());
    gfile << "gene\tpeak\tdistance" << std::endl;
    for(uint32_t i = 0; i < geneIds.size(); ++i) {
      if (!geneView[i].empty()) {
	if (c.nearest) {
	  int32_t bestDistance = 250000000;
	  int32_t bestIdx = -1;
	  for(typename TPeaksPerGene::const_iterator itDP = geneView[i].begin(); itDP != geneView[i].end(); ++itDP) {
	    if (std::abs(itDP->first) < std::abs(bestDistance)) {
	      bestIdx = itDP->second;
	      bestDistance = itDP->first;
	    }
	  }
	  gfile << geneIds[i] << "\t" << peakNames[bestIdx] << "\t" << bestDistance << std::endl;
	} else {
	  gfile << geneIds[i] << "\t";
	  std::sort(geneView[i].begin(), geneView[i].end());
	  bool firstF = true;
	  for(typename TPeaksPerGene::const_iterator itDP = geneView[i].begin(); itDP != geneView[i].end(); ++itDP) {
	    if (!firstF) gfile << ',';
	    else firstF = false;
	    gfile << peakNames[itDP->second];
	  }
	  gfile << "\t";
	  firstF = true;
	  for(typename TPeaksPerGene::const_iterator itDP = geneView[i].begin(); itDP != geneView[i].end(); ++itDP) {
	    if (!firstF) gfile << ',';
	    else firstF = false;
	    gfile << itDP->first;
	  }
	  gfile << std::endl;
	}
      }
    }
    gfile.close();

    // Done
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
    else if (c.inputFileFormat == 3) tf = parseJaspar(c, gRegions, geneIds);
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
      ("outgene,u", boost::program_options::value<boost::filesystem::path>(&c.outgene)->default_value("gene.bed"), "gene/motif-level output")
      ("outfile,o", boost::program_options::value<boost::filesystem::path>(&c.outfile)->default_value("anno.bed"), "annotated peaks output")
      ("nearest,n", "nearest feature only")
      ;

    boost::program_options::options_description gtfopt("GTF/GFF3 annotation file options");
    gtfopt.add_options()
      ("gtf,g", boost::program_options::value<boost::filesystem::path>(&c.gtfFile), "gtf/gff3 file")
      ("id,i", boost::program_options::value<std::string>(&c.idname)->default_value("gene_name"), "gtf/gff3 attribute")
      ("feature,f", boost::program_options::value<std::string>(&c.feature)->default_value("gene"), "gtf/gff3 feature")
      ;

    boost::program_options::options_description motifopt("Motif annotation file options");
    motifopt.add_options()
      ("motif,m", boost::program_options::value<boost::filesystem::path>(&c.motifFile), "motif file in jaspar or raw format")
      ("reference,r", boost::program_options::value<boost::filesystem::path>(&c.genome), "reference file")
      ("quantile,q", boost::program_options::value<float>(&c.motifScoreQuantile)->default_value(0.95), "motif quantile score [0,1]")
      ("position,p", boost::program_options::value<boost::filesystem::path>(&c.outpos), "gzipped output file of motif hits")
      ("exclude,x", "exclude overlapping hits of the same motif")
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
    cmdline_options.add(generic).add(gtfopt).add(bedopt).add(motifopt).add(hidden);
    boost::program_options::options_description visible_options;
    visible_options.add(generic).add(gtfopt).add(bedopt).add(motifopt);

    // Parse command-line
    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::command_line_parser(argc, argv).options(cmdline_options).positional(pos_args).run(), vm);
    boost::program_options::notify(vm);

    // Check command line arguments
    if ((vm.count("help")) || (!vm.count("input-file")) || ((!vm.count("gtf")) && (!vm.count("bed")) && (!vm.count("motif")))) {
      std::cout << std::endl;
      std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] -g <hg19.gtf.gz> <peaks.bed>" << std::endl;
      std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] -b <hg19.bed.gz> <peaks.bed>" << std::endl;
      std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] -m <motif.jaspar.gz> -r <genome.fa> <peaks.bed>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    }

    // Motif position
    if (vm.count("position")) c.motifPosOut = true;
    else c.motifPosOut = false;

    // Nearest feature only
    if (vm.count("nearest")) c.nearest = true;
    else c.nearest = false;

    // Overlapping motif hits
    if (vm.count("exclude")) c.overlappingHits = false;
    else c.overlappingHits = true;
    
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
	if (!(boost::filesystem::exists(c.motifFile) && boost::filesystem::is_regular_file(c.motifFile) && boost::filesystem::file_size(c.motifFile))) {
	  std::cerr << "Input gtf/bed/motif annotation file is missing." << std::endl;
	  return 1;
	} else {
	  c.inputFileFormat = 3;
	  if ((!(vm.count("reference")))  || (!(boost::filesystem::exists(c.genome) && boost::filesystem::is_regular_file(c.genome) && boost::filesystem::file_size(c.genome)))) {
	    std::cerr << "Motif annotation requires a reference genome file." << std::endl;
	    return 1;
	  } else {
	    faidx_t* fai = fai_load(c.genome.string().c_str());
	    if (fai == NULL) {
	      if (fai_build(c.genome.string().c_str()) == -1) {
		std::cerr << "Fail to open genome fai index for " << c.genome.string() << std::endl;
		return 1;
	      } else fai = fai_load(c.genome.string().c_str());
	    }
	    // Check that all chromosomes in input file are present
	    for(typename AnnotateConfig::TChrMap::const_iterator itC = c.nchr.begin(); itC != c.nchr.end(); ++itC) {
	      std::string chrName(itC->first);
	      if (!faidx_has_seq(fai, chrName.c_str())) {
		 std::cerr << "Chromosome from bed file " << chrName << " is NOT present in your reference file " << c.genome.string() << std::endl;
		 return 1;
	      }
	    }
	    fai_destroy(fai);
	  }
	}
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
