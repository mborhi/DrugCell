import argparse
# import pandas as pd

def construct_mutation_line(mutations_list, gene2ind):
    binary = ["0"] * len(gene2ind.keys())# 3008
    for mutation in mutations_list:
        if mutation in gene2ind:
            binary[gene2ind[mutation]] = "1"
    
    return binary



def write_sample2binary_mut(filename, data_points, gene2ind):
    with open(filename, 'w') as f:
        ind = 0
        for data in data_points.keys():
            binary_list = construct_mutation_line(data_points[data], gene2ind)
            line = ",".join(binary_list) + "\n"
            f.write(line)


def read_data(filename):
    collection = {}
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            tokens = line.rstrip().split("\t")
            collection[tokens[0]] = tokens[1].split(",")
    return collection

def get_all_mutated_genes(filename, gene_filename):
    # genes_df = pd.read_csv(filename, sep="\t")
    all_genes = set()
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines[1:]:
            tokens = line.strip().split("\t")
            all_genes.add(tokens[0])
    # all_genes = genes_df["Gene"]
    with open(gene_filename, 'w') as f:
        ind = 0
        for gene in all_genes:
            line = str(ind) + "\t" + gene + "\n"
            f.write(line)
            ind += 1

def get_gene2ind(filename):
    gene2ind = {}
    with open(filename, 'r') as f:
        lines = f.readlines()
        for line in lines:
            tokens = line.rstrip().split("\t")
            gene2ind[tokens[1]] = int(tokens[0])
    
    return gene2ind

def write_sample2ind(file, samples):
    with open(file, 'w') as f:
        ind = 0
        for sample in samples:
            line = str(ind) + "\t" + sample.strip() + "\n"
            f.write(line)
            ind += 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Write supplementary file data')
    parser.add_argument('-output_dir', help='The directory to write new data files to', type=str)
    parser.add_argument('-sample2mutations', help='The file containing sample, mutation pairs', type=str)
    parser.add_argument('-gene2ind', help='The file containing sample, mutation pairs', type=str)
    opt = parser.parse_args()

    output_dir = opt.output_dir
    
    sample2mutations = read_data(opt.sample2mutations)
    gene2ind = get_gene2ind(opt.gene2ind)
    # write_sample2ind(output_dir + "reactome-sample2ind.txt", list(sample2mutations.keys()))
    # write_sample2binary_mut(output_dir + "reactome-sample2binary_mutation.txt", sample2mutations, gene2ind)
    write_sample2binary_mut(output_dir + "top-gene-reactome-sample2binary_mutation.txt", sample2mutations, gene2ind)
    # get_all_mutated_genes(output_dir + "All_Mutated_Genes.txt", output_dir + "all-cbioPortal-gene2ind.txt")
    # cbioPortal_gene2ind = get_gene2ind(output_dir + "all-cbioPortal-gene2ind.txt")
    # write_sample2binary_mut(output_dir + "all-cbioPortal-sample2binary-mut.txt", sample2mutations, cbioPortal_gene2ind)
