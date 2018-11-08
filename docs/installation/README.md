# Installation

Alfred is available as a [Bioconda package](https://anaconda.org/bioconda/alfred), a pre-compiled statically linked binary from the [Alfred github release page](https://github.com/tobiasrausch/alfred/releases/) or as a minimal [Docker container](https://hub.docker.com/r/trausch/alfred/). All code is open-source and hosted on [Alfred's GitHub page](https://github.com/tobiasrausch/alfred).

## Installation from Source

To build Alfred from source you need some build essentials and the Boost libraries, i.e. for Ubuntu:

```bash
apt install \
    build-essential g++ \
    cmake \
    git-all \
    liblzma-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev
```

Once you have installed these system libraries you can compile and link Alfred.

```bash
git clone --recursive https://github.com/tobiasrausch/alfred.git
cd alfred/
make all
make install
./bin/alfred -h
```

## Non-default Boost installation directory

Delly requires Boost and you can install Boost locally using

```bash
git clone --recursive https://github.com/boostorg/boost.git
cd boost/
./bootstrap.sh --prefix=`pwd` --without-icu --with-libraries=iostreams,filesystem,system,program_options,date_time
./b2
./b2 headers
```

You can then specify a non-default Boost installation directory (i.e., /opt/boost below) using

```bash
make CXXFLAGS=-I/opt/boost LDFLAGS=-L/opt/boost/stage/lib all 
```

To clean the local Boost installation

```bash
./b2 --clean-all
```
