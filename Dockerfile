# use the ubuntu base image
FROM ubuntu:22.04

MAINTAINER Tobias Rausch rausch@embl.de

# install required packages
RUN apt-get update && apt-get install -y \
    autoconf \
    build-essential \
    cmake \
    g++ \
    gfortran \
    git \
    libcurl4-gnutls-dev \
    hdf5-tools \
    libboost-date-time-dev \
    libboost-program-options-dev \
    libboost-system-dev \
    libboost-filesystem-dev \
    libboost-iostreams-dev \
    libbz2-dev \
    libdeflate-dev \
    libhdf5-dev \
    libncurses-dev \
    liblzma-dev \
    zlib1g-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# set environment
ENV BOOST_ROOT /usr

# install alfred
RUN cd /opt \
    && git clone --recursive https://github.com/tobiasrausch/alfred.git \
    && cd /opt/alfred/ \
    && make STATIC=1 all \
    && make install

# Multi-stage build
FROM alpine:latest
RUN apk add --no-cache bash
RUN mkdir -p /opt/alfred/bin
WORKDIR /opt/alfred/bin
COPY --from=0 /opt/alfred/bin/alfred .

# Workdir
WORKDIR /home

# Add alfred to PATH
ENV PATH="/opt/alfred/bin:${PATH}"

# by default /bin/sh is executed
CMD ["/bin/bash"]
