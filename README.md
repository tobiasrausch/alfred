bamStats installation (using recursive clone)
------------------------------------------

`git clone --recursive https://github.com/tobiasrausch/bamStats.git`

`cd bamStats/`

`make all`

bamStats installation (using package manager)
-------------------------------------------

Install Linux packages:

`apt-get update`

`apt-get install -y build-essential g++ git cmake zlib1g-dev ant libbz2-dev libboost-date-time-dev libboost-program-options-dev libboost-system-dev libboost-filesystem-dev libboost-iostreams-dev`

Install htslib:

`git clone https://github.com/samtools/htslib.git`

`cd htslib && make && make lib-static && cd ..`

Set environment variables for boost libraries and htslib:

`export BOOST_ROOT=/usr`

`export SEQTK_ROOT=<htslib_path>`

Build bamStats (no recursive clone required):

`git clone https://github.com/tobiasrausch/bamStats.git`

`cd bamStats/ && touch .htslib .boost && make all && cd ..`


Running bamStats
------------------

`./src/bamStats -r <ref.fa> -o outprefix <align.bam>`

`Rscript R/isize.R outprefix.isize.tsv`

`Rscript R/coverage.R outprefix.coverage.tsv`
	


