SHELL := /bin/bash

# Targets
TARGETS = .conda .channels .alfred .test
PBASE=$(shell pwd)

all:  	$(TARGETS)

.conda:
	wget 'https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh' && bash Miniconda3-latest-Linux-x86_64.sh -b -p ${PBASE}/bin && rm -f Miniconda3-latest-Linux-x86_64.sh && touch .conda

.channels: .conda
	export PATH=${PBASE}/bin/bin:${PATH} && conda config --add channels defaults && conda config --add channels conda-forge && conda config --add channels bioconda && touch .channels

.alfred: .conda .channels
	export PATH=${PBASE}/bin/bin:${PATH} && conda install alfred && touch .alfred

.test: .conda .channels .alfred
	export PATH=${PBASE}/bin/bin:${PATH} && alfred --help && touch .test

clean:
	export PATH=${PBASE}/bin/bin:${PATH} && conda remove alfred && rm .alfred

distclean:
	rm -rf $(TARGETS) $(TARGETS:=.o) bin/

.PHONY: distclean clean all
