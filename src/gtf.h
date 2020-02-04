#ifndef GTF_H
#define GTF_H

#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/algorithm/string.hpp>

#include <htslib/sam.h>

#include "util.h"

namespace bamstats
{

  template<typename TConfig, typename TGenomicRegions, typename TGeneIds, typename TProteinCoding>
  inline int32_t
  parseGTFAll(TConfig const& c, TGenomicRegions& overlappingRegions, TGeneIds& geneIds, TProteinCoding& pCoding) {
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "GTF feature parsing" << std::endl;

    // Check gzip
    if (!is_gz(c.gtfFile)) {
      std::cerr << "GTF file is not gzipped!" << std::endl;
      return 0;
    }

    // Map IDs to integer
    typedef std::map<std::string, int32_t> TIdMap;
    TIdMap idMap;

    // Keep track of unique exon IDs
    int32_t eid = 0;

    // Parse GTF
    std::ifstream file(c.gtfFile.string().c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> dataIn;
    dataIn.push(boost::iostreams::gzip_decompressor());
    dataIn.push(file);
    std::istream instream(&dataIn);
    std::string gline;
    while(std::getline(instream, gline)) {
      if ((gline.size()) && (gline[0] == '#')) continue;
      typedef boost::tokenizer< boost::char_separator<char> > Tokenizer;
      boost::char_separator<char> sep("\t");
      Tokenizer tokens(gline, sep);
      Tokenizer::iterator tokIter = tokens.begin();
      if (tokIter==tokens.end()) {
	std::cerr << "Empty line in GTF file!" << std::endl;
	return 0;
      }
      std::string chrName=*tokIter++;
      if (c.nchr.find(chrName) == c.nchr.end()) continue;
      int32_t chrid = c.nchr.find(chrName)->second;      
      if (tokIter == tokens.end()) {
	std::cerr << "Corrupted GTF file!" << std::endl;
	return 0;
      }
      ++tokIter;
      if (tokIter == tokens.end()) {
	std::cerr << "Corrupted GTF file!" << std::endl;
	return 0;
      }
      std::string ft = *tokIter++;
      if (ft == c.feature) {
	if (tokIter != tokens.end()) {
	  int32_t start = boost::lexical_cast<int32_t>(*tokIter++);
	  int32_t end = boost::lexical_cast<int32_t>(*tokIter++);
	  ++tokIter; // score
	  if (tokIter == tokens.end()) {
	    std::cerr << "Corrupted GTF file!" << std::endl;
	    return 0;
	  }
	  char strand = boost::lexical_cast<char>(*tokIter++);
	  ++tokIter; // frame
	  std::string attr = *tokIter;
	  boost::char_separator<char> sepAttr(";");
	  Tokenizer attrTokens(attr, sepAttr);
	  for(Tokenizer::iterator attrIter = attrTokens.begin(); attrIter != attrTokens.end(); ++attrIter) {
	    std::string keyval = *attrIter;
	    boost::trim(keyval);
	    boost::char_separator<char> sepKeyVal(" ");
	    Tokenizer kvTokens(keyval, sepKeyVal);
	    Tokenizer::iterator kvTokensIt = kvTokens.begin();
	    std::string key = *kvTokensIt++;
	    if (key == c.idname) {
	      std::string val = *kvTokensIt;
	      if (val.size() >= 3) val = val.substr(1, val.size()-2); // Trim off the bloody "
	      int32_t idval = geneIds.size();
	      typename TIdMap::const_iterator idIter = idMap.find(val);
	      if (idIter == idMap.end()) {
		idMap.insert(std::make_pair(val, idval));
		geneIds.push_back(val);
		// Protein Coding?
		bool pCode = false;
		for(Tokenizer::iterator arIter = attrTokens.begin(); arIter != attrTokens.end(); ++arIter) {
		  std::string kvl = *arIter;
		  boost::trim(kvl);
		  boost::char_separator<char> sKV2(" ");
		  Tokenizer kvT2(kvl, sKV2);
		  Tokenizer::iterator kvT2It = kvT2.begin();
		  std::string procod = *kvT2It++;
		  if (procod == "gene_biotype") {
		    std::string gbio = *kvT2It;
		    if (gbio.size() >= 3) gbio = gbio.substr(1, gbio.size()-2);
		    if (gbio == "protein_coding") pCode = true;
		  }
		}
		pCoding.push_back(pCode);
	      } else idval = idIter->second;
	      // Convert to 0-based and right-open
	      if (start == 0) {
		std::cerr << "GTF is 1-based format!" << std::endl;
		return 0;
	      }
	      if (start > end) {
		std::cerr << "Feature start is greater than feature end!" << std::endl;
		return 0;
	      }
	      //std::cerr << geneIds[idval] << "\t" << start << "\t" << end << std::endl;
	      _insertInterval(overlappingRegions[chrid], start - 1, end, strand, idval, eid++);
	    }
	  }
	}
      }
    }
    return geneIds.size();
  }


  template<typename TConfig, typename TGenomicRegions, typename TGeneIds>
  inline int32_t
  parseGTFAll(TConfig const& c, TGenomicRegions& overlappingRegions, TGeneIds& geneIds) {
    std::vector<bool> pCoding;
    return parseGTFAll(c, overlappingRegions, geneIds, pCoding);
  }
  
  template<typename TConfig, typename TGenomicRegions, typename TGeneIds, typename TProteinCoding>
  inline int32_t
  parseGTF(TConfig const& c, TGenomicRegions& gRegions, TGeneIds& geneIds, TProteinCoding& pCoding) {
    typedef typename TGenomicRegions::value_type TChromosomeRegions;

    // Overlapping intervals for each label
    TGenomicRegions overlappingRegions;
    overlappingRegions.resize(gRegions.size(), TChromosomeRegions());
    parseGTFAll(c, overlappingRegions, geneIds, pCoding);
    
    // Make intervals non-overlapping for each label
    for(uint32_t refIndex = 0; refIndex < overlappingRegions.size(); ++refIndex) {
      // Sort by ID
      std::sort(overlappingRegions[refIndex].begin(), overlappingRegions[refIndex].end(), SortIntervalLabel<IntervalLabel>());
      int32_t runningId = -1;
      char runningStrand = '*';
      typedef boost::icl::interval_set<uint32_t> TIdIntervals;
      typedef typename TIdIntervals::interval_type TIVal;
      TIdIntervals idIntervals;
      for(uint32_t i = 0; i < overlappingRegions[refIndex].size(); ++i) {
	if (overlappingRegions[refIndex][i].lid != runningId) {
	  for(typename TIdIntervals::iterator it = idIntervals.begin(); it != idIntervals.end(); ++it) {
	    //std::cerr << "merged\t" << geneIds[runningId] << "\t" << it->lower() << "\t" << it->upper() << std::endl;  
	    gRegions[refIndex].push_back(IntervalLabel(it->lower(), it->upper(), runningStrand, runningId));
	  }
	  idIntervals.clear();
	  runningId = overlappingRegions[refIndex][i].lid;
	  runningStrand = overlappingRegions[refIndex][i].strand;
	}
	idIntervals.insert(TIVal::right_open(overlappingRegions[refIndex][i].start, overlappingRegions[refIndex][i].end));
      }
      // Process last id
      for(typename TIdIntervals::iterator it = idIntervals.begin(); it != idIntervals.end(); ++it) gRegions[refIndex].push_back(IntervalLabel(it->lower(), it->upper(), runningStrand, runningId));
    }
    
    return geneIds.size();
  }

  template<typename TConfig, typename TGenomicRegions, typename TGeneIds>
  inline int32_t
  parseGTF(TConfig const& c, TGenomicRegions& gRegions, TGeneIds& geneIds) {
    std::vector<bool> pCoding;
    return parseGTF(c, gRegions, geneIds, pCoding);
  }
}

#endif
