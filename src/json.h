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

#ifndef JSON_H
#define JSON_H

#include <boost/progress.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

namespace bamstats
{
  

  template<typename TConfig, typename THeader, typename TSampleWC>
  inline void
  scJsonOut(TConfig const& c, THeader const& hdr, TSampleWC const& sWC) {
    boost::iostreams::filtering_ostream rfile;
    rfile.push(boost::iostreams::gzip_compressor());
    rfile.push(boost::iostreams::file_sink(c.outfile.string().c_str(), std::ios_base::out | std::ios_base::binary));
    rfile << "{" << std::endl;
    rfile << "\"data\": [" << std::endl;
    for(uint32_t i = 0; i < c.sampleName.size(); ++i) {
      if (i > 0) rfile << "," << std::endl;
      rfile << "{" << std::endl;
      rfile << "\"sample\": \"" << c.sampleName[i] << "\"," << std::endl;
      rfile << "\"coverages\": [" << std::endl;
      for (int32_t refIndex = 0; refIndex<hdr[i]->n_targets; ++refIndex) {
	if (!sWC[0][refIndex].size()) continue;
	if (refIndex > 0) rfile << "," << std::endl;
	rfile << "{" << std::endl;
	rfile << "\"chromosome\": \"" << hdr[i]->target_name[refIndex] << "\"," << std::endl;
	rfile << "\"positions\": [" << std::endl;
	int32_t pos = 0;
	for (uint32_t k = 0; k < sWC[i][refIndex].size(); ++k) {
	  if (k > 0) rfile << ", ";
	  rfile << pos;
	  pos += c.window;
	}
	rfile << "]," << std::endl;
	rfile << "\"counts\": [" << std::endl;
	rfile << "{" << std::endl;
	rfile << "\"label\": \"Watson\"," << std::endl;
	rfile << "\"values\": [" << std::endl;
	for (uint32_t k = 0; k < sWC[i][refIndex].size(); ++k) {
	  if (k > 0) rfile << ", ";
	  rfile << sWC[i][refIndex][k].first;
	}
	rfile << "]" << std::endl;
	rfile << "}," << std::endl;
	rfile << "{" <<	std::endl;
	rfile << "\"label\": \"Crick\"," << std::endl;
	rfile << "\"values\": [" << std::endl;
	for (uint32_t k = 0; k < sWC[i][refIndex].size(); ++k) {
	  if (k > 0) rfile << ", ";
	  rfile << sWC[i][refIndex][k].second;
	}
	rfile << "]" << std::endl;
	rfile << "}" << std::endl;
	rfile << "]" << std::endl;
	rfile << "}";
      }
      rfile << "]" << std::endl;
      rfile << "}";
    }
    rfile << "]" << std::endl;
    rfile << std::endl;
    rfile << "}" << std::endl;
    rfile.pop();
  }
 
}

#endif
