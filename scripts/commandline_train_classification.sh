#!/bin/bash
inputdir="../classification_data/"
gene2idfile=$inputdir"gene2ind.txt" 
cell2idfile=$inputdir"sample2ind.txt" #"cell2ind.txt"
# drug2idfile=$inputdir"drug2ind.txt"
traindatafile=$inputdir"sample2type_test_ILC_IDC.txt" #"mock_train_data.txt"
valdatafile=$inputdir"sample2type_train_ILC_IDC.txt" #"mock_test_data.txt"
ontfile=$inputdir"drugcell_ont.txt"

mutationfile=$inputdir"sample2binary_mutation.txt" #"cell2mutation.txt"
# drugfile=$inputdir"drug2fingerprint.txt"

cudaid=0

modeldir=ClassificationModelTest
mkdir $modeldir

source activate pytorch3drugcell

# python -u ../code/train_drugcell.py -onto $ontfile -gene2id $gene2idfile -drug2id $drug2idfile -cell2id $cell2idfile -train $traindatafile -test $valdatafile -model $modeldir -cuda $cudaid -genotype $mutationfile -fingerprint $drugfile -genotype_hiddens 6 -drug_hiddens '100,50,6' -final_hiddens 6 -epoch 100 -batchsize 5000 > train_script.log
python -u ../code/train_drugcell.py -onto $ontfile -gene2id $gene2idfile -cell2id $cell2idfile -train $traindatafile -test $valdatafile -model $modeldir -cuda $cudaid -genotype $mutationfile -genotype_hiddens 6 -final_hiddens 2 -epoch 50 -batchsize 5000 > train_script.log
