#!/bin/bash
inputdir="../classification_data/"
gene2idfile=$inputdir"gene2ind.txt" #"reactome-genes2ind.txt" 
cell2idfile=$inputdir"reactome-sample2ind.txt"
testdatafile=$inputdir"sample2type_test_ILC_IDC.txt"

mutationfile=$inputdir"sample2binary_mutation.txt" # "reactome-sample2binary_mutation.txt"

modelfile="./ReactomeModelTest/model_final.pt"

resultdir="Reactome_Classification_Result_sample"
hiddendir="Reactome_Classification_Hidden_sample"

cudaid=0

mkdir $resultdir
mkdir $hiddendir

source activate pytorch3drugcell

python -u ../code/predict_drugcell.py -gene2id $gene2idfile -cell2id $cell2idfile -genotype $mutationfile -hidden $hiddendir -result $resultdir -predict $testdatafile -load $modelfile -cuda $cudaid > test_sample.log
