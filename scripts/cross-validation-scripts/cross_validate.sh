#!/bin/bash
collection="Stomach_Esophageal"
pathway_files="../../classification_data/"
inputdir="../../datasets/"
collection_data="../../datasets/"$collection"/"
#[datasets]/[collection]/[collection]_
datadir=$inputdir$collection"/"
gene2idfile=$collection_data$collection"_top_varied_gene2ind.txt"
cell2idfile=$collection_data$collection"_sample2ind.txt" 
ontfile=$pathway_files"ReactomePathwaysStd.txt"
splitsdir=$inputdir$collection"/splits_txt/"

traindatafile=$collection_data$collection"_sample2type.txt" # not used
valdatafile=$collection_data$collection"_sample2type.txt" # not used

mutationfile=$collection_data$collection"_sample2binary_top_varied_gene_mutation.txt"

cudaid=1

modeldir=$collection"_Cross_Validation"
# mkdir $modeldir
mkdir $collection

source activate pytorch3drugcell

python -u ../../code/validation.py \
    -onto $ontfile \
    -splitsdir $splitsdir \
    -gene2id $gene2idfile \
    -cell2id $cell2idfile \
    -train $traindatafile \
    -test $valdatafile \
    -modeldir $modeldir \
    -basedir $collection \
    -cuda $cudaid \
    -genotype $mutationfile \
    -genotype_hiddens 6 \
    -final_hiddens 2 \
    -epoch 300 \
    -batchsize 3000 \
    > $collection"_cross_validation.log"