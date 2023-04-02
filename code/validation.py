import argparse
import os 
import numpy as np
import torch
import torch.utils.data as du
from torch.autograd import Variable
import torch.nn as nn
import torch.nn.functional as F
import train_drugcell
# import util
from util import *
# from drugcell_NN import *
# import argparse


def prepare_data(test_file, train_files, cell2id_mapping):
    """
    Params:
    test_file: path to the file to keep as testing
    train_files: list of paths to the files to use as training
    cell2id_mapping: mapping from the sample to id
    """
    train_feature = []
    train_label = []
    for train_file in train_files:
        # exclude the test_file
        if test_file == train_file:
            continue
        # load the train file
        with open(train_file, 'r') as fi:
            for line in fi:
                tokens = line.strip().split('\t')

                train_feature.append([cell2id_mapping[tokens[0]]])
                train_label.append([float(tokens[1])])
        
        # train_feature, train_label = load_train_data(train_file, cell2id_mapping)
        
    test_feature, test_label = load_train_data(test_file, cell2id_mapping)
    
    return (torch.Tensor(train_feature), torch.FloatTensor(train_label), torch.Tensor(test_feature), torch.FloatTensor(test_label))


def cross_validate(data_dir, onto_file):
    print(f"Starting Cross Validation...")
    filepaths = [data_dir + "/" + file for file in os.listdir(data_dir)]
    results = [] # list of dics
    iter = 0
    for filepath in filepaths:
        train_data = prepare_data(filepath, filepaths, cell2id_mapping)
        print(f"Fold {iter + 1} / 5")
        store_dir = base_dir + "/" + model_save_folder + "_" + str(iter)
        if not os.path.exists(store_dir):
            os.makedirs(store_dir)
        dG, root, term_size_map, term_direct_gene_map = load_ontology(onto_file, gene2id_mapping)
        max_corr, max_mcc, max_acc, f1, neg_f1 = train_drugcell.train_model(root, term_size_map, term_direct_gene_map, dG, train_data, gene_dim, store_dir, train_epochs, batch_size, learning_rate, num_hiddens_genotype, num_hiddens_final, cell_features)

        results.append({"f1" : f1,
                        "corr": max_corr,
                        "mcc": max_mcc,
                        "acc": max_acc,
                        "neg_f1": neg_f1
                        })
        iter += 1

    corr_avg, mcc_avg, acc_avg, f1_avg, neg_f1_avg = 0, 0, 0, 0, 0
    for result in results:
        corr_avg += result["corr"]
        mcc_avg += result["mcc"]
        acc_avg += result["acc"]
        f1_avg += result["f1"]
        neg_f1_avg += result["neg_f1"]

    
    print(f"Pearson correlation average: {corr_avg / 5}")
    print(f"MCC average: {mcc_avg / 5}")
    print(f"Accuracy average: {acc_avg / 5}")
    print(f"F1 score average: {f1_avg / 5}")
    print(f"Negative F1 score average: {neg_f1_avg / 5}")
    

    return results
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Validate dcell')

    # parser.add_argument('-datadir', help='directory containing data splits for k-fold cross validation', type=str)
    parser.add_argument('-onto', help='Ontology file used to guide the neural network', type=str)
    parser.add_argument('-splitsdir', help='Directory containing data splits', action='store', type=str)
    parser.add_argument('-train', help='Training dataset', type=str)
    parser.add_argument('-test', help='Validation dataset', type=str)
    parser.add_argument('-epoch', help='Training epochs for training', type=int, default=300)
    parser.add_argument('-lr', help='Learning rate', type=float, default=0.001)
    parser.add_argument('-batchsize', help='Batchsize', type=int, default=3000)
    parser.add_argument('-basedir', help='Base folder for trained models', type=str, default='FOLDS/')
    parser.add_argument('-modeldir', help='Folder for trained models', type=str, default='MODEL/')
    parser.add_argument('-cuda', help='Specify GPU', type=int, default=0)
    parser.add_argument('-gene2id', help='Gene to ID mapping file', type=str)
    parser.add_argument('-cell2id', help='Cell to ID mapping file', type=str)
    parser.add_argument('-genotype_hiddens', help='Mapping for the number of neurons in each term in genotype parts', type=int, default=2)
    parser.add_argument('-final_hiddens', help='The number of neurons in the top layer', type=int, default=6)
    parser.add_argument('-genotype', help='Mutation information for cell lines', type=str)

    opt = parser.parse_args()
    torch.set_printoptions(precision=5)

    # load input data
    train_data, cell2id_mapping = prepare_train_data(opt.train, opt.test, opt.cell2id)
    gene2id_mapping = load_mapping(opt.gene2id)

    # load cell features
    cell_features = np.genfromtxt(opt.genotype, delimiter=',')
    num_cells = len(cell2id_mapping)
    gene_dim = len(gene2id_mapping)

    # load ontology
    # dG, root, term_size_map, term_direct_gene_map = load_ontology(opt.onto, gene2id_mapping)

    # load the number of hiddens #######
    num_hiddens_genotype = opt.genotype_hiddens

    # num_hiddens_drug = list(map(int, opt.drug_hiddens.split(',')))

    num_hiddens_final = opt.final_hiddens
    #####################################

    CUDA_ID = opt.cuda
    model_save_folder, train_epochs, batch_size, learning_rate = opt.modeldir, opt.epoch, opt.batchsize, opt.lr
    base_dir = opt.basedir

    # data_dir = "../../classification_data/cross_validation_data/"
    data_dir = opt.splitsdir
    onto_file = opt.onto
    cross_validate(data_dir, onto_file)