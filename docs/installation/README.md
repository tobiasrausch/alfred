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

## Non-default boost installation directory

You can set a non-default boost installation directory using


