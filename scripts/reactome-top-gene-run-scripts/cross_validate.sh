#!/bin/bash
inputdir="../../classification_data/"
datadir=$inputdir"cross_validation_data/"
gene2idfile=$inputdir"top_varied_gene2ind.txt"
cell2idfile=$inputdir"reactome-sample2ind.txt" 
traindatafile=$inputdir"sample2type_test_ILC_IDC.txt"
valdatafile=$inputdir"sample2type_train_ILC_IDC.txt"
ontfile=$inputdir"ReactomePathwaysStd.txt"

mutationfile=$inputdir"top-gene-reactome-sample2binary_mutation.txt"

cudaid=1

modeldir=Reactome_Top_Genes_Cross_Validation
mkdir $modeldir

source activate pytorch3drugcell

python -u ../../code/validation.py -onto $ontfile -gene2id $gene2idfile -cell2id $cell2idfile -train $traindatafile -test $valdatafile -model $modeldir -cuda $cudaid -genotype $mutationfile -genotype_hiddens 6 -final_hiddens 2 -epoch 300 -batchsize 3000 > cross_validation.log