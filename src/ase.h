#ifndef ASE_H
#define ASE_H

#include <iostream>
#include <vector>
#include <fstream>

#include <boost/math/distributions/binomial.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/stream_buffer.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/progress.hpp>
#include <htslib/faidx.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>

#include "util.h"
#include "variants.h"

namespace bamstats {


  struct AseConfig {
    bool isPhased;
    bool outputAll;
    unsigned short minMapQual;
    unsigned short minBaseQual;
    std::string sample;
    boost::filesystem::path as;
    boost::filesystem::path genome;
    boost::filesystem::path bamfile;
    boost::filesystem::path vcffile;
  };

  template<typename TConfig>
  inline int32_t
  aseRun(TConfig& c) {

#ifdef PROFILE
    ProfilerStart("alfred.prof");
#endif
    
    // Load bam files
    samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
    hts_set_fai_filename(samfile, c.genome.string().c_str());
    hts_idx_t* idx = sam_index_load(samfile, c.bamfile.string().c_str());
    bam_hdr_t* hdr = sam_hdr_read(samfile);

    // Load bcf file
    htsFile* ibcffile = bcf_open(c.vcffile.string().c_str(), "r");
    hts_idx_t* bcfidx = bcf_index_load(c.vcffile.string().c_str());
    bcf_hdr_t* bcfhdr = bcf_hdr_read(ibcffile);
    
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] " << "Assign reads to haplotypes" << std::endl;
    boost::progress_display show_progress(hdr->n_targets);

    // Allele support file
    boost::iostreams::filtering_ostream dataOut;
    dataOut.push(boost::iostreams::gzip_compressor());
    dataOut.push(boost::iostreams::file_sink(c.as.string().c_str(), std::ios_base::out | std::ios_base::binary));
    dataOut << "chr\tpos\tid\tref\talt\tdepth\trefsupport\taltsupport\tgt\taf\tpvalue" << std::endl;
  
    // Assign reads to SNPs
    faidx_t* fai = fai_load(c.genome.string().c_str());
    for (int refIndex = 0; refIndex<hdr->n_targets; ++refIndex) {
      std::string chrName(hdr->target_name[refIndex]);
      ++show_progress;

      // Load het. markers
      typedef std::vector<BiallelicVariant> TPhasedVariants;
      TPhasedVariants pv;
      if (!_loadVariants(ibcffile, bcfidx, bcfhdr, c.sample, chrName, pv)) continue;
      if (pv.empty()) continue;

      // Sort variants
      std::sort(pv.begin(), pv.end(), SortVariants<BiallelicVariant>());

      // Load reference
      int32_t seqlen = -1;
      char* seq = NULL;
      seq = faidx_fetch_seq(fai, chrName.c_str(), 0, hdr->target_len[refIndex], &seqlen);
      
      // Annotate REF and ALT support
      typedef std::vector<uint32_t> TAlleleSupport;
      TAlleleSupport ref(pv.size(), 0);
      TAlleleSupport alt(pv.size(), 0);
      hts_itr_t* itr = sam_itr_queryi(idx, refIndex, 0, hdr->target_len[refIndex]);
      bam1_t* r = bam_init1();
      while (sam_itr_next(samfile, itr, r) >= 0) {
	if (r->core.flag & (BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY | BAM_FUNMAP)) continue;
	if ((r->core.qual < c.minMapQual) || (r->core.tid<0)) continue;
	if ((r->core.flag & BAM_FPAIRED) && (r->core.flag & BAM_FMUNMAP)) continue;

	// Fetch contained variants
	TPhasedVariants::const_iterator vIt = std::lower_bound(pv.begin(), pv.end(), BiallelicVariant(r->core.pos), SortVariants<BiallelicVariant>());
	TPhasedVariants::const_iterator vItEnd = std::upper_bound(pv.begin(), pv.end(), BiallelicVariant(lastAlignedPosition(r)), SortVariants<BiallelicVariant>());
	if (vIt != vItEnd) {
	  // Get read sequence
	  std::string sequence;
	  sequence.resize(r->core.l_qseq);
	  uint8_t* seqptr = bam_get_seq(r);
	  for (int32_t i = 0; i < r->core.l_qseq; ++i) sequence[i] = "=ACMGRSVTWYHKDBN"[bam_seqi(seqptr, i)];

	  // Get base qualities
	  typedef std::vector<uint8_t> TQuality;
	  TQuality quality;
	  quality.resize(r->core.l_qseq);
	  uint8_t* qualptr = bam_get_qual(r);
	  for (int i = 0; i < r->core.l_qseq; ++i) quality[i] = qualptr[i];
	  
	  // Parse CIGAR
	  uint32_t* cigar = bam_get_cigar(r);
	  for(;vIt != vItEnd; ++vIt) {
	    int32_t gp = r->core.pos; // Genomic position
	    int32_t sp = 0; // Sequence position
	    bool varFound = false;
	    for (std::size_t i = 0; ((i < r->core.n_cigar) && (!varFound)); ++i) {
	      if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) sp += bam_cigar_oplen(cigar[i]);
	      else if (bam_cigar_op(cigar[i]) == BAM_CINS) sp += bam_cigar_oplen(cigar[i]);
	      else if (bam_cigar_op(cigar[i]) == BAM_CSOFT_CLIP) sp += bam_cigar_oplen(cigar[i]);
	      else if (bam_cigar_op(cigar[i]) == BAM_CDEL) gp += bam_cigar_oplen(cigar[i]);
	      else if (bam_cigar_op(cigar[i]) == BAM_CREF_SKIP) gp += bam_cigar_oplen(cigar[i]);
	      else if (bam_cigar_op(cigar[i]) == BAM_CHARD_CLIP) {
		//Nop
	      } else if ((bam_cigar_op(cigar[i]) == BAM_CMATCH) || (bam_cigar_op(cigar[i]) == BAM_CEQUAL) || (bam_cigar_op(cigar[i]) == BAM_CDIFF)) {
		if (gp + (int32_t) bam_cigar_oplen(cigar[i]) < vIt->pos) {
		  gp += bam_cigar_oplen(cigar[i]);
		  sp += bam_cigar_oplen(cigar[i]);
		} else {
		  for(std::size_t k = 0; k<bam_cigar_oplen(cigar[i]); ++k, ++sp, ++gp) {
		    if (gp == vIt->pos) {
		      varFound = true;
		      if (quality[sp] >= c.minBaseQual) {
			// Check REF allele
			if (vIt->ref == std::string(seq + gp, seq + gp + vIt->ref.size())) {
			  // Check ALT allele
			  if ((sp + vIt->alt.size() < sequence.size()) && (sp + vIt->ref.size() < sequence.size())) {
			    if (vIt->ref.size() == vIt->alt.size()) {
			      // SNP
			      if ((sequence.substr(sp, vIt->alt.size()) == vIt->alt) && (sequence.substr(sp, vIt->ref.size()) != vIt->ref)) {
				++alt[vIt-pv.begin()];
			      } else if ((sequence.substr(sp, vIt->alt.size()) != vIt->alt) && (sequence.substr(sp, vIt->ref.size()) == vIt->ref)) {
				++ref[vIt-pv.begin()];
			      }
			    }
			  } else if (vIt->ref.size() < vIt->alt.size()) {
			    // Insertion
			    int32_t diff = vIt->alt.size() - vIt->ref.size();
			    std::string refProbe = vIt->ref + std::string(seq + gp + vIt->ref.size(), seq + gp + vIt->ref.size() + diff);
			    if ((sequence.substr(sp, vIt->alt.size()) == vIt->alt) && (sequence.substr(sp, vIt->alt.size()) != refProbe)) {
			      ++alt[vIt-pv.begin()];
			    } else if ((sequence.substr(sp, vIt->alt.size()) != vIt->alt) && (sequence.substr(sp, vIt->alt.size()) == refProbe)) {
			      ++ref[vIt-pv.begin()];
			    }
			  } else {
			    // Deletion
			    int32_t diff = vIt->ref.size() - vIt->alt.size();
			    std::string altProbe = vIt->alt + std::string(seq + gp + vIt->ref.size(), seq + gp + vIt->ref.size() + diff);
			    if ((sequence.substr(sp, vIt->ref.size()) == altProbe) && (sequence.substr(sp, vIt->ref.size()) != vIt->ref)) {
			      ++alt[vIt-pv.begin()];
			    } else if ((sequence.substr(sp, vIt->ref.size()) != altProbe) && (sequence.substr(sp, vIt->ref.size()) == vIt->ref)) {
			      ++ref[vIt-pv.begin()];
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	      else {
		std::cerr << "Unknown Cigar options" << std::endl;
		return 1;
	      }
	    }
	  }
	}
      }
      bam_destroy1(r);
      hts_itr_destroy(itr);
      if (seqlen) free(seq);

      // Output (phased) allele support
      hts_itr_t* itervcf = bcf_itr_querys(bcfidx, bcfhdr, chrName.c_str());
      if (itervcf != NULL) {
	bcf1_t* recvcf = bcf_init1();
	for (uint32_t i = 0; i<pv.size(); ++i) {
	  // Fetch variant annotation from VCF
	  int32_t itrRet = 0;
	  do {
	    itrRet = bcf_itr_next(ibcffile, itervcf, recvcf);
	    if (itrRet >= 0) {
	      bcf_unpack(recvcf, BCF_UN_SHR);
	      std::vector<std::string> alleles;
	      for(std::size_t k = 0; k<recvcf->n_allele; ++k) alleles.push_back(std::string(recvcf->d.allele[k]));
	      if ((recvcf->pos == pv[i].pos) && (pv[i].ref == alleles[0]) && (pv[i].alt == alleles[1])) break;
	    } else {
	      std::cerr << "Error: Variant not found! " << chrName << ":" << (pv[i].pos + 1) << std::endl;
	      return 1;
	    }
	  } while (itrRet >= 0);
	  uint32_t totalcov = ref[i] + alt[i];
	  std::string hapstr = "0/1";
	  if (c.isPhased) {
	    if (pv[i].hap) hapstr = "1|0";
	    else hapstr = "0|1";
	  }
	  if (totalcov > 0) {
	    double h1af = 0;
	    double vaf = (double) alt[i] / (double) totalcov;
	    if (pv[i].hap) h1af = (double) alt[i] / (double) totalcov;
	    else h1af = (double) ref[i] / (double) totalcov;
	    double pval = binomTest(alt[i], totalcov, 0.5);
	    dataOut << chrName << "\t" << (pv[i].pos + 1) << "\t" << recvcf->d.id << "\t" << pv[i].ref << "\t" << pv[i].alt << "\t" << totalcov << "\t" << ref[i] << "\t" << alt[i] << "\t" << hapstr << "\t";
	    if (c.isPhased) dataOut << h1af << "\t";
	    else dataOut << vaf << "\t";
	    dataOut << pval << std::endl;
	  } else {
	    if (c.outputAll) {
	      // No coverage
	      dataOut << chrName << "\t" << (pv[i].pos + 1) << "\t" << recvcf->d.id << "\t" << pv[i].ref << "\t" << pv[i].alt << "\t" << totalcov << "\t" << ref[i] << "\t" << alt[i] << "\t" << hapstr << "\tNA\tNA" << std::endl;
	    }
	  }
	}
	bcf_destroy(recvcf);
	hts_itr_destroy(itervcf);
      }
    }
    fai_destroy(fai);

    // Close bam
    bam_hdr_destroy(hdr);
    hts_idx_destroy(idx);
    sam_close(samfile);

    // Close output allele file
    dataOut.pop();

    // Close BCF
    bcf_hdr_destroy(bcfhdr);
    hts_idx_destroy(bcfidx);
    bcf_close(ibcffile);
    
    // End
    now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] Done." << std::endl;
    
#ifdef PROFILE
    ProfilerStop();
#endif
    
    
    return 0;
  }

  int ase(int argc, char **argv) {
    AseConfig c;

    // Parameter
    boost::program_options::options_description generic("Generic options");
    generic.add_options()
      ("help,?", "show help message")
      ("map-qual,m", boost::program_options::value<unsigned short>(&c.minMapQual)->default_value(10), "min. mapping quality")
      ("base-qual,b", boost::program_options::value<unsigned short>(&c.minBaseQual)->default_value(10), "min. base quality")
      ("reference,r", boost::program_options::value<boost::filesystem::path>(&c.genome), "reference fasta file")
      ("sample,s", boost::program_options::value<std::string>(&c.sample)->default_value("NA12878"), "sample name")
      ("ase,a", boost::program_options::value<boost::filesystem::path>(&c.as)->default_value("as.tsv.gz"), "allele-specific output file")
      ("vcffile,v", boost::program_options::value<boost::filesystem::path>(&c.vcffile), "input (phased) BCF file")
      ("phased,p", "BCF file is phased and BAM is haplo-tagged")
      ("full,f", "output all het. input SNPs")
      ;
    
    boost::program_options::options_description hidden("Hidden options");
    hidden.add_options()
      ("input-file", boost::program_options::value<boost::filesystem::path>(&c.bamfile), "input bam file")
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
    if ((vm.count("help")) || (!vm.count("input-file")) || (!vm.count("reference")) || (!vm.count("vcffile"))) {
      std::cout << "Usage: alfred " << argv[0] << " [OPTIONS] -r <ref.fa> -s NA12878 -v <snps.bcf> -a <ase.tsv> <input.bam>" << std::endl;
      std::cout << visible_options << "\n";
      return 1;
    }

    // Phased running mode?
    if (!vm.count("phased")) c.isPhased = false;
    else c.isPhased = true;

    // Output all input het. SNPs
    if (!vm.count("full")) c.outputAll = false;
    else c.outputAll = true;
    
    // Check input BAM file
    if (vm.count("input-file")) {
      if (!(boost::filesystem::exists(c.bamfile) && boost::filesystem::is_regular_file(c.bamfile) && boost::filesystem::file_size(c.bamfile))) {
	std::cerr << "Input BAM file is missing: " << c.bamfile.string() << std::endl;
	return 1;
      }
      samFile* samfile = sam_open(c.bamfile.string().c_str(), "r");
      if (samfile == NULL) {
      std::cerr << "Fail to open file " << c.bamfile.string() << std::endl;
      return 1;
      }
      hts_idx_t* idx = sam_index_load(samfile, c.bamfile.string().c_str());
      if (idx == NULL) {
	std::cerr << "Fail to open index for " << c.bamfile.string() << std::endl;
	return 1;
      }
      bam_hdr_t* hdr = sam_hdr_read(samfile);
      if (hdr == NULL) {
	std::cerr << "Fail to open header for " << c.bamfile.string() << std::endl;
	return 1;
      }
      bam_hdr_destroy(hdr);
      hts_idx_destroy(idx);
      sam_close(samfile);
    }
    
    // Check VCF/BCF file
    if (vm.count("vcffile")) {
      if (!(boost::filesystem::exists(c.vcffile) && boost::filesystem::is_regular_file(c.vcffile) && boost::filesystem::file_size(c.vcffile))) {
	std::cerr << "Input SNP VCF/BCF file is missing: " << c.vcffile.string() << std::endl;
	return 1;
      }
      htsFile* ifile = bcf_open(c.vcffile.string().c_str(), "r");
      if (ifile == NULL) {
	std::cerr << "Fail to open file " << c.vcffile.string() << std::endl;
	return 1;
      }
      hts_idx_t* bcfidx = bcf_index_load(c.vcffile.string().c_str());
      if (bcfidx == NULL) {
	std::cerr << "Fail to open index file for " << c.vcffile.string() << std::endl;
	return 1;
      }
      bcf_hdr_t* hdr = bcf_hdr_read(ifile);
      if (hdr == NULL) {
	std::cerr << "Fail to open header for " << c.vcffile.string() << std::endl;
	return 1;
      }
      bcf_hdr_destroy(hdr);
    hts_idx_destroy(bcfidx);
    bcf_close(ifile);
    }
    
    // Show cmd
    boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
    std::cout << '[' << boost::posix_time::to_simple_string(now) << "] ";
    for(int i=0; i<argc; ++i) { std::cout << argv[i] << ' '; }
    std::cout << std::endl;
    
    return aseRun(c);
  }


}

#endif
