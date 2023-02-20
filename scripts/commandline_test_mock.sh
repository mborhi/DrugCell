#!/bin/bash
inputdir="../custom_data/"
gene2idfile=$inputdir"gene2ind.txt"
cell2idfile=$inputdir"cell2ind.txt"
testdatafile=$inputdir"mock_test_data.txt"

mutationfile=$inputdir"cell2mutation.txt"

modelfile="./ModelTest/model_final.pt"

resultdir="Result_sample"
hiddendir="Hidden_sample"

cudaid=3

# if [$cudaid = ""]; then
# 	cudaid=0
# fi

mkdir $resultdir
mkdir $hiddendir

source activate pytorch3drugcell

python -u ../code/predict_drugcell.py -gene2id $gene2idfile -cell2id $cell2idfile -genotype $mutationfile -hidden $hiddendir -result $resultdir -predict $testdatafile -load $modelfile -cuda $cudaid > test_sample.log
