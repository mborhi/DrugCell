#!/bin/bash
inputdir="../../classification_data/"
gene2idfile=$inputdir"top_varied_gene2ind.txt"
cell2idfile=$inputdir"reactome-sample2ind.txt"
testdatafile=$inputdir"sample2type_test_ILC_IDC.txt"
model_type="f1"
model_fold="4"
mutationfile=$inputdir"top-gene-reactome-sample2binary_mutation.txt"

modelfile="./Reactome_Top_Genes_Cross_Validation_All_Metrics_"$model_fold"/model_best_"$model_type".pt"

resultdir="Reactome_Top_Genes_Classification_Cross_Val_Result_sample"
hiddendir="Reactome_Top_Genes_Classification_Cross_Val_Hidden_sample"

cudaid=1

mkdir $resultdir
mkdir $hiddendir

source activate pytorch3drugcell

python -u ../../code/predict_drugcell.py -gene2id $gene2idfile -cell2id $cell2idfile -genotype $mutationfile -hidden $hiddendir -result $resultdir -predict $testdatafile -load $modelfile -cuda $cudaid > test_cross_val_$model_type.log
