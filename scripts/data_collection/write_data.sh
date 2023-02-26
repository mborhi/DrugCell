#!/bin/bash
data_dir="./../../classification_data/"
sample2mutations=$data_dir"sample2mutations_ILC_IDC.txt"
sample2mutations=$data_dir"gene2ind.txt"

source activate pytorch3drugcell

python -u ./write_data.py -output_dir=$data_dir -sample2mutations=$combined_file -gene2ind