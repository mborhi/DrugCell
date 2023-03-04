import sys
import argparse
import networkx as nx
import networkx.algorithms.components.connected as nxacc
import networkx.algorithms.dag as nxadag

def write_hierarchy_file(filename, pathway_relations, pathway_genes):
    pathways = set(pathway_relations.keys()).union(set(pathway_genes.keys()))
    with open(filename, 'w') as f:
        """
        format: [Pathway] \\t [gene|pathway] [gene|default]
        """
        for pathway in pathways:
            # get relations
            if pathway in pathway_relations:
                relations = pathway_relations[pathway]
                for relation in relations:
                    line = pathway + "\t" + relation + "\tdefault\n"
                    f.write(line)
            if pathway in pathway_genes:
                genes = pathway_genes[pathway]
                for gene in genes:
                    f.write(pathway + "\t" + gene + "\tgene\n")

def write_all_gene2ind(filename, all_genes):
    with open(filename, 'w') as f: 
        ind = 0
        for gene in all_genes:
            f.write(str(ind) + "\t" + gene + "\n")
            ind += 1
            
def collect_pathway_data(relations_file, pathways_file):

    def clean_genes(gene_list):
        cleaned = set()
        for gene in gene_list:
            if "gene" not in gene:
                cleaned.add(gene)
        return list(cleaned)
    
    all_genes = set()

    pathway_to_genes = {}
    # Note: maybe ignore first line
    with open(pathways_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            tokens = line.strip().split("\t")
            pathway_name = tokens[1]
            genes = tokens[3:]
            if len(genes) > 0:
                genes = clean_genes(genes)
                pathway_to_genes[pathway_name] = set(genes)
                all_genes = all_genes.union(set(genes))

    # str -> set(str)
    pathway_relations = {}
    with open(relations_file, 'r') as f:
        lines = f.readlines()
        for line in lines:
            connection = line.strip().split("\t")
            pathway_one = connection[0]
            pathway_two = connection[1]
            if pathway_one in pathway_relations:
                pathway_relations[pathway_one].add(pathway_two)
            else :
                pathway_relations[pathway_one] = set([pathway_two]) # to prevent string split

    return pathway_to_genes, pathway_relations, all_genes

def get_total_genes(term, dG, subsystem_to_genes, total_genes: set):
    if term in subsystem_to_genes:
        total_genes = total_genes.union(subsystem_to_genes[term])

    for child in dG.successors(term):
        get_total_genes(child, dG, subsystem_to_genes, total_genes)

    return total_genes


def load_reactome(file_name, gene2id_mapping):

    dG = nx.DiGraph()
    term_direct_gene_map = {}
    term_size_map = {}

    file_handle = open(file_name)

    gene_set = set()

    for line in file_handle:

        line = line.rstrip().split()
        
        if line[2] == 'default':
            dG.add_edge(line[0], line[1])
        else:
            if line[1] not in gene2id_mapping:
                continue

            if line[0] not in term_direct_gene_map:
                term_direct_gene_map[ line[0] ] = set()

            term_direct_gene_map[line[0]].add(gene2id_mapping[line[1]])

            gene_set.add(line[1])

    file_handle.close()

    print('There are', len(gene_set), 'genes')

    nodes_to_remove = []
    for term in dG.nodes():
        predecessor = [p for p in dG.predecessors(term)]
        if len(predecessor) == 0:
            continue
        predecessor = predecessor[0]
        deslist = nxadag.descendants(dG, term)

        term_total_genes = get_total_genes(term, dG, term_direct_gene_map, set())
        if len(term_total_genes) == 0:
            nodes_to_remove.append(term)
            for child in deslist:
                dG.add_edge(predecessor, child)

    dG.remove_nodes_from(nodes_to_remove)

    leaves = [n for n in dG.nodes if dG.in_degree(n) == 0]
    #leaves = [n for n,d in dG.in_degree() if d==0]

    uG = dG.to_undirected()
    connected_subG_list = list(nxacc.connected_components(uG))

    root_node = 'root'
    edges = [(root_node, n) for n in leaves]
    dG.add_edges_from(edges)

    print('There are', len(leaves), 'roots:', leaves)
    print('There are', len(dG.nodes()), 'terms')
    print('There are', len(connected_subG_list), 'connected componenets')

    leaves = [n for n in dG.nodes if dG.in_degree(n) == 0]
    uG = dG.to_undirected()
    connected_subG_list = list(nxacc.connected_components(uG))
    print("Making adjustments...")
    print('There are', len(leaves), 'roots:', leaves)
    print('There are', len(dG.nodes()), 'terms')
    print('There are', len(connected_subG_list), 'connected componenets')

    for term in dG.nodes():
        
        term_gene_set = set()

        if term in term_direct_gene_map:
            term_gene_set = term_direct_gene_map[term]

        deslist = nxadag.descendants(dG, term)

        for child in deslist:
            if child in term_direct_gene_map:
                term_gene_set = term_gene_set | term_direct_gene_map[child] # this is union

        # jisoo
        if len(term_gene_set) == 0:
            # print("no terms? but there are:", get_total_genes(term, dG, term_direct_gene_map, set()))
            print('There is empty terms, please delete term:', term)
            sys.exit(1)
        else:
            term_size_map[term] = len(term_gene_set)

    if len(leaves) > 1:
        print('There are more than 1 root of ontology. Please use only one root.')
        sys.exit(1)
    if len(connected_subG_list) > 1:
        print( 'There are more than connected components. Please connect them.')
        sys.exit(1)

    return dG, leaves[0], term_size_map, term_direct_gene_map

def load_mapping(mapping_file):

    mapping = {}

    file_handle = open(mapping_file)

    for line in file_handle:
        line = line.rstrip().split("\t")
        mapping[line[1]] = int(line[0])

    file_handle.close()

    return mapping

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Write Reactome Pathway data')
    parser.add_argument('-output_dir', help='The directory to write new data files to', type=str)
    parser.add_argument('-pathway_relations', help='The file containing pathway relations', type=str)
    parser.add_argument('-pathway_genes', help='The gmt file containing pathway to gene information', type=str)
    opt = parser.parse_args()

    pathway_to_genes, pathway_relations, all_genes = collect_pathway_data(opt.pathway_relations, opt.pathway_genes)
    write_all_gene2ind(opt.output_dir + "reactome-genes2ind.txt",all_genes)
    write_hierarchy_file(opt.output_dir + "ReactomePathwaysStd.txt", pathway_relations, pathway_to_genes)
    gene2ind = load_mapping(opt.output_dir + "reactome-genes2ind.txt")

    ### Test hierarchy load
    # load_reactome(opt.output_dir + "ReactomePathwaysStd.txt", gene2ind)
