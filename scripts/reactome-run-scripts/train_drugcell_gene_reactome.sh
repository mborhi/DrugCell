#!/bin/bash
inputdir="../../classification_data/"
gene2idfile=$inputdir"gene2ind.txt"
cell2idfile=$inputdir"reactome-sample2ind.txt" 
traindatafile=$inputdir"sample2type_test_ILC_IDC.txt"
valdatafile=$inputdir"sample2type_train_ILC_IDC.txt"
ontfile=$inputdir"ReactomePathwaysStd.txt"

mutationfile=$inputdir"sample2binary_mutation.txt"

cudaid=0

modeldir=Reactome_DrugCell_Gene_Model
mkdir $modeldir

source activate pytorch3drugcell

python -u ../../code/train_drugcell.py -onto $ontfile -gene2id $gene2idfile -cell2id $cell2idfile -train $traindatafile -test $valdatafile -model $modeldir -cuda $cudaid -genotype $mutationfile -genotype_hiddens 6 -final_hiddens 2 -epoch 100 -batchsize 1000 > train_script.log