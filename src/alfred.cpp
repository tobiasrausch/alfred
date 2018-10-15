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

#define _SECURE_SCL 0
#define _SCL_SECURE_NO_WARNINGS
#include <iostream>
#include <vector>
#include <fstream>

#define BOOST_DISABLE_ASSERTS
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/tokenizer.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#ifdef PROFILE
#include "gperftools/profiler.h"
#endif


#include "version.h"
#include "util.h"
#include "bamstats.h"
#include "count_rna.h"
#include "count_dna.h"
#include "count_junction.h"
#include "annotate.h"
#include "tracks.h"
#include "split.h"
#include "ase.h"
#include "qc.h"
#include "consensus.h"
#include "pwalign.h"

using namespace bamstats;


inline void
asciiArt() {
  std::cout << "           _  __              _ " << std::endl;
  std::cout << "     /\\   | |/ _|            | |" << std::endl;
  std::cout << "    /  \\  | | |_ _ __ ___  __| |" << std::endl;
  std::cout << "   / /\\ \\ | |  _| '__/ _ \\/ _` |" << std::endl;
  std::cout << "  / ____ \\| | | | | |  __/ (_| |" << std::endl;
  std::cout << " /_/    \\_\\_|_| |_|  \\___|\\__,_|" << std::endl;
  std::cout << std::endl;
}

inline void
displayUsage() {
  std::cout << "Usage: alfred <command> <arguments>" << std::endl;
  std::cout << std::endl;
  std::cout << "Commands:" << std::endl;
  std::cout << std::endl;
  std::cout << "    qc           alignment quality control" << std::endl;
  std::cout << "    count_dna    counting DNA reads in windows" << std::endl;
  std::cout << "    count_rna    counting RNA reads in features" << std::endl;
  std::cout << "    count_jct    counting RNA split-reads at exon junctions" << std::endl;
  std::cout << "    tracks       create browser tracks" << std::endl;
  std::cout << "    annotate     annotate peaks" << std::endl;
  std::cout << "    split        split BAM into haplotypes" << std::endl;
  std::cout << "    consensus    consensus computation for error-prone reads" << std::endl;
  std::cout << "    pwalign      pairwise alignment using dynamic programming" << std::endl;
  std::cout << "    ase          allele-specific expression" << std::endl;
  std::cout << std::endl;
  std::cout << std::endl;
}


int main(int argc, char **argv) {
  if (argc < 2) {
    asciiArt();
    printTitle("Alfred");
    displayUsage();
    return 0;
  }
  
  if ((std::string(argv[1]) == "version") || (std::string(argv[1]) == "--version") || (std::string(argv[1]) == "--version-only") || (std::string(argv[1]) == "-v")) {
    std::cout << "Alfred version: v" << alfredVersionNumber << std::endl;
    return 0;
  }
  else if ((std::string(argv[1]) == "help") || (std::string(argv[1]) == "--help") || (std::string(argv[1]) == "-h") || (std::string(argv[1]) == "-?")) {
    printTitle("Alfred");
    displayUsage();
    return 0;
  }
  else if ((std::string(argv[1]) == "warranty") || (std::string(argv[1]) == "--warranty") || (std::string(argv[1]) == "-w")) {
    displayWarranty();
    return 0;
  }
  else if ((std::string(argv[1]) == "license") || (std::string(argv[1]) == "--license") || (std::string(argv[1]) == "-l")) {
    gplV3();
    return 0;
  }
  else if ((std::string(argv[1]) == "qc")) {
    return qc(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "count_rna")) {
    return count_rna(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "count_dna")) {
    return count_dna(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "count_jct")) {
    return count_junction(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "tracks")) {
    return tracks(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "annotate")) {
    return annotate(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "split")) {
    return split(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "consensus")) {
    return consensus(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "pwalign")) {
    return pwalign(argc-1,argv+1);
  }
  else if ((std::string(argv[1]) == "ase")) {
    return ase(argc-1,argv+1);
  }
  std::cerr << "Unrecognized command " << std::string(argv[1]) << std::endl;
  return 1;
}
