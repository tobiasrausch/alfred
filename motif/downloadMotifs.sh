#!/bin/bash

curl --output jaspar.zip 'http://jaspar.genereg.net/download/CORE/JASPAR2018_CORE_vertebrates_non-redundant_pfms_jaspar.zip'
unzip jaspar.zip
rm jaspar.zip
cat *.jaspar | gzip -c > jaspar.gz
rm *.jaspar
