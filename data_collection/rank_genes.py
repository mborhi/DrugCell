import numpy as np
import argparse
import pandas as pd

def rank_gene_variances(filename: str, percentage: float, save=True, plot=True) -> pd.DataFrame:
    df = pd.read_csv(filename, sep="\t", header=None,
                     names=["Gene", "Variance"])
    sorted_df = df.sort_values("Variance", ascending=False)
    top_variances = sorted_df.head(int(df.shape[0] * percentage))

    if save:
        df.to_csv("gene_variances.csv")
        sorted_df.to_csv("sorted_gene_variances.csv")
        top_variances.to_csv("top_gene_variances.csv")

    return sorted_df


def read_gene_file(filename: str) -> set:
    gene_set = set()
    with open(filename, 'r') as f:
        for line in f.readlines()[1:]:
            gene_set.add(line.strip().split("\t")[0])

    return gene_set
        
def find_variance(gene2patient_mut, output_dir):
    gene_var = {}
    for gene in gene2patient_mut:
        gene_var[gene] = np.var(gene2patient_mut[gene])
    
    with open(output_dir + "gene_variances.txt", 'w') as f:
        for gene in gene_var:
            line = gene + "\t" + str(gene_var[gene]) + "\n"
            f.write(line)
    
    return gene_var, "gene_variances.txt"

def get_gene_patient_matrix(gene2patient_ids, patient_count):
    gene2patient_mut = {}
    for gene in gene2patient_ids.keys():
        patient_line = [0] * patient_count
        patients_with_this_gene_mutated = gene2patient_ids[gene]
        for ind in patients_with_this_gene_mutated:
            patient_line[int(ind)] = 1
        gene2patient_mut[gene] = patient_line
    
    return gene2patient_mut

def get_gene2patient_ids(filename, patient2id):
    gene2patient_ids = {}
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines[1:]:
            tokens = line.strip().split("\t")
            # tokens[0] is sample, tokens[1].split(",")
            for gene in tokens[1].split(","):
                if gene in gene2patient_ids:
                    gene2patient_ids[gene].append(patient2id[tokens[0]])
                else :
                    gene2patient_ids[gene] = [patient2id[tokens[0]]]

    return gene2patient_ids

def get_patient2id(filename):
    patient2id = {}
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            tokens = line.strip().split("\t")
            patient2id[tokens[1]] = tokens[0]
    return patient2id

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Write supplementary file data')
    parser.add_argument('-output_dir', help='The directory to write new data files to', type=str)
    parser.add_argument('-input_dir', help='The directory to read data files from', type=str)
    parser.add_argument('-sample2mutations', help='The file containing sample, mutation pairs', type=str)
    opt = parser.parse_args()

    # Find the variance of genes
    patient2id = get_patient2id(opt.input_dir + "reactome-sample2ind.txt")
    gene2patient_ids = get_gene2patient_ids(opt.sample2mutations, patient2id)
    gene2patient_mut = get_gene_patient_matrix(gene2patient_ids, len(patient2id.keys()))
    gene_variance_dict, gene_variance_file = find_variance(gene2patient_mut, opt.output_dir)
    # Sort and write top varied genes
    sorted_df = rank_gene_variances(
        "gene_variance_file", percentage=0.20, save=False, plot=False)
    