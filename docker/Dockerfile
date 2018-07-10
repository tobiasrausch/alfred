# use the ubuntu base image
FROM ubuntu:14.04

MAINTAINER Tobias Rausch rausch@embl.de

# install required packages
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    g++ \
    gfortran \
    git \
    hdf5-tools \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libbz2-dev \
    libhdf5-dev \
    libncurses-dev \
    liblzma-dev \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# set environment
ENV BOOST_ROOT /usr
ENV EBROOTHTSLIB /opt/htslib

# install alfred
RUN cd /opt \
    && git clone https://github.com/samtools/htslib.git \
    && cd /opt/htslib \
    && make \
    && make lib-static \
    && make install
RUN cd /opt \
    && git clone https://github.com/tobiasrausch/alfred.git \
    && cd /opt/alfred/ \
    && make STATIC=1 all \
    && make install

# Multi-stage build
FROM alpine:latest
RUN mkdir -p /opt/alfred/bin
WORKDIR /opt/alfred/bin
COPY --from=0 /opt/alfred/bin/alfred .

# Workdir
WORKDIR /root/

# Add Alfred to PATH
ENV PATH="/opt/alfred/bin:${PATH}"

# by default /bin/sh is executed
CMD ["/bin/sh"]
