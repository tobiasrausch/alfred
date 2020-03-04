#!/bin/bash

# Raw motifs
curl --output jaspar.zip 'http://jaspar.genereg.net/download/CORE/JASPAR2020_CORE_vertebrates_non-redundant_pfms_jaspar.zip'
unzip jaspar.zip
rm jaspar.zip
cat *.jaspar | gzip -c > jaspar.gz
rm *.jaspar

# Clustered motifs
curl --output JASPAR_2020_matrix_clustering_vertebrates_cluster_root_motifs.tf http://folk.uio.no/jamondra/JASPAR_2020_clusters/vertebrates/interactive_trees/JASPAR_2020_matrix_clustering_vertebrates_cluster_root_motifs.tf
python ./convert.py -f JASPAR_2020_matrix_clustering_vertebrates_cluster_root_motifs.tf  > jaspar.cluster
rm JASPAR_2020_matrix_clustering_vertebrates_cluster_root_motifs.tf
gzip jaspar.cluster
