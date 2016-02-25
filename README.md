pbBamStats installation (using recursive clone)
------------------------------------------

`git clone --recursive https://github.com/tobiasrausch/pbBamStats.git`

`cd pbBamStats/`

`make all`

pbBamStats installation (using package manager)
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

Build pbBamStats (no recursive clone required):

`git clone https://github.com/tobiasrausch/pbBamStats.git`

`cd pbBamStats/ && touch .htslib .boost && make all && cd ..`


Running pbBamStats
------------------

`./src/pbBamStats -r example/ref.fa example/align.bam`


