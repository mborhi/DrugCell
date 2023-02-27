#!/bin/bash
data_dir="./../../classification_data/"
pathway_relations=$data_dir"ReactomeRelations.txt"
pathway_genes=$data_dir"ReactomePathways.gmt"

source activate pytorch3drugcell

python -u ./reactome.py -output_dir=$data_dir -pathway_relations=$pathway_relations -pathway_genes=$pathway_genes