/*
============================================================================
Alfred: BAM alignment statistics
============================================================================
Copyright (C) 2017 Tobias Rausch

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

#ifndef BAMSTATS_H
#define BAMSTATS_H

#include <limits>

#include <boost/dynamic_bitset.hpp>
#include <boost/unordered_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/progress.hpp>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "util.h"
#include "gtf.h"


namespace bamstats
{

  template<typename TConfig>
  inline int32_t
  countRun(TConfig const& c) {
    // Load bam file
    samFile* samfile = sam_open(c.bamFile.string().c_str(), "r");
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Parse GTF file
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "GTF file parsing" << std::endl;

    typedef std::vector<IntervalLabel> TChromosomeRegions;
    typedef std::vector<TChromosomeRegions> TGenomicRegions;
    TGenomicRegions gRegions;
    gRegions.resize(hdr->n_targets, TChromosomeRegions());
    typedef std::vector<std::string> TGeneIds;
    TGeneIds geneIds;
    parseGTF(hdr, c.gtfFile, "exon", "gene_id", gRegions, geneIds);

    
    // Debug code
    //uint32_t rIndex = 0;
    //for(TGenomicRegions::const_iterator itG = gRegions.begin(); itG != gRegions.end(); ++itG, ++rIndex) {
    //for(TChromosomeRegions::const_iterator itC = itG->begin(); itC != itG->end(); ++itC) {
    //std::cout << rIndex << ',' << hdr->target_name[rIndex] << ',' << itC->start << ',' << itC->end << std::endl;
    //}
    //}

    // Parse BAM file
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "BAM file parsing" << std::endl;
    boost::progress_display sp( hdr->n_targets );


    // clean-up
    bam_hdr_destroy(hdr);
    sam_close(samfile);
    
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    
#ifdef PROFILE
    ProfilerStop();
#endif


    return 0;
  }

}

#endif
