#!/bin/sh
data_dir="./../../classification_data/"
combined_file=$data_dir"combined_study_clinical_ILC_IDC.tsv"

# source activate pytorch3drugcell

python -u ./clean_data.py -output_dir=$data_dir -combined_file=$combined_file