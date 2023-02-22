#!/bin/bash
inputdir="../custom_data/"
gene2idfile=$inputdir"gene2ind.txt"
cell2idfile=$inputdir"cell2ind.txt"
# drug2idfile=$inputdir"drug2ind.txt"
traindatafile=$inputdir"mock_train_data.txt"
valdatafile=$inputdir"mock_test_data.txt"
ontfile=$inputdir"drugcell_ont.txt"

mutationfile=$inputdir"cell2mutation.txt"
# drugfile=$inputdir"drug2fingerprint.txt"

cudaid=3

modeldir=ModelTest
mkdir $modeldir

source activate pytorch3drugcell

# python -u ../code/train_drugcell.py -onto $ontfile -gene2id $gene2idfile -drug2id $drug2idfile -cell2id $cell2idfile -train $traindatafile -test $valdatafile -model $modeldir -cuda $cudaid -genotype $mutationfile -fingerprint $drugfile -genotype_hiddens 6 -drug_hiddens '100,50,6' -final_hiddens 6 -epoch 100 -batchsize 5000 > train_script.log
python -u ../code/train_drugcell.py -onto $ontfile -gene2id $gene2idfile -cell2id $cell2idfile -train $traindatafile -test $valdatafile -model $modeldir -cuda $cudaid -genotype $mutationfile -genotype_hiddens 6 -final_hiddens 2 -epoch 10 -batchsize 2000 > train_script.log
