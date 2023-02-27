#!/bin/bash
data_dir="./../../classification_data/"
sample2mutations=$data_dir"sample2mutations_ILC_IDC.txt"
gene2ind=$data_dir"reactome-genes2ind.txt"

source activate pytorch3drugcell

python -u ./write_data.py -output_dir=$data_dir -sample2mutations=$sample2mutations -gene2ind=$gene2ind