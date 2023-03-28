import sys
import torch
import networkx as nx
import networkx.algorithms.components.connected as nxacc
import networkx.algorithms.dag as nxadag
import numpy as np

def pearson_corr(x, y):
	xx = x - torch.mean(x)
	yy = y - torch.mean(y)

	return torch.sum(xx*yy) / (torch.norm(xx, 2)*torch.norm(yy,2))

def accuracy(x, y):
	preds = torch.argmax(x, dim=1).long()
	bools = torch.eq(preds, torch.squeeze(y.long()))
	acc = torch.sum(bools).item() / len(preds)
	return acc

def confusion_matrix(x, y):
	preds = torch.argmax(x, dim=1).long()
	targets = torch.squeeze(y).long()
	# TODO genearlize to N by N
	c_m = torch.zeros(2, 2)
	for i, target in enumerate(targets):
		c_m[preds[i]][target] += 1
	return c_m

def matthew_cc(c_m):
	tp = c_m[0][0]
	tn = c_m[1][1]
	fp = c_m[0][1]
	fn = c_m[1][0]
	numerator = (tp * tn) - (fp * fn)
	denom = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
	return numerator / denom

def specific_agreement(category, conf_mat):
    """Calculates Specific Agreement scores for dichotomous catgories.
    
    https://github.com/jmgirard/mReliability/wiki/Specific-agreement-coefficient

    Inputs
    ------
    category: int
        Either 0 or 1, signifying the positive and negative categories respectively. 
    conf-mat: [[int]]
        The confusion matrix.
    
    Returns
    -------
    float 
        The Specific Agreement score for the given category.
    """
    num = 2 * conf_mat[category][category]
    denom = 2 * conf_mat[category][category] + conf_mat[0][1] + conf_mat[1][0]
    return num / denom

def f1_score(c_m):
	tp = c_m[0][0]
	tn = c_m[1][1]
	fp = c_m[0][1]
	fn = c_m[1][0]
	precision = tp / (tp + fp)
	recall = tp / (tp + fn)
	return 2 * (precision * recall) / (precision + recall)

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

    min_req_genes = 0#30
    min_distinct = 0#10
    nodes_to_remove = []
    for term in dG.nodes():
        predecessor = [p for p in dG.predecessors(term)]
        if len(predecessor) == 0:
            continue
        predecessor = predecessor[0]
        deslist = nxadag.descendants(dG, term)

        children_genes = set()
        child_genes_list = []
        for child in deslist:
            child_genes = get_total_genes(
                child, dG, term_direct_gene_map, set())
            children_genes = children_genes.union(child_genes)
            child_genes_list.append(child_genes)

        # term_genes = set()
        if term in term_direct_gene_map:
            term_genes = term_direct_gene_map[term]
            term_genes = term_genes.union(children_genes)
        else:
            term_genes = children_genes
        # has at least min_req_genes in its hierarchy
        enough_genes = len(term_genes) >= min_req_genes

        # has to have at least min_distinct unique genes in its hierarchy than any child
        distinct = all(len(term_genes - child_genes) >=
                        min_distinct for child_genes in child_genes_list)

        term_total_genes = get_total_genes(term, dG, term_direct_gene_map, set())
        if len(term_total_genes) == 0 or (not (enough_genes and distinct)):
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

def load_ontology(file_name, gene2id_mapping):
	# print("file name to load ont from:", file_name)
	if "Reactome" in file_name:
		return load_reactome(file_name, gene2id_mapping)
	
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

	for term in dG.nodes():
		
		term_gene_set = set()

		if term in term_direct_gene_map:
			term_gene_set = term_direct_gene_map[term]

		deslist = nxadag.descendants(dG, term)

		for child in deslist:
			if child in term_direct_gene_map:
				term_gene_set = term_gene_set | term_direct_gene_map[child]

		# jisoo
		if len(term_gene_set) == 0:
			print('There is empty terms, please delete term:', term)
			sys.exit(1)
		else:
			term_size_map[term] = len(term_gene_set)

	leaves = [n for n in dG.nodes if dG.in_degree(n) == 0]
	#leaves = [n for n,d in dG.in_degree() if d==0]

	uG = dG.to_undirected()
	connected_subG_list = list(nxacc.connected_components(uG))

	print('There are', len(leaves), 'roots:', leaves[0])
	print('There are', len(dG.nodes()), 'terms')
	print('There are', len(connected_subG_list), 'connected componenets')

	if len(leaves) > 1:
		print('There are more than 1 root of ontology. Please use only one root.')
		sys.exit(1)
	if len(connected_subG_list) > 1:
		print( 'There are more than connected components. Please connect them.')
		sys.exit(1)

	return dG, leaves[0], term_size_map, term_direct_gene_map


# def load_train_data(file_name, cell2id, drug2id):
# 	feature = []
# 	label = []

# 	with open(file_name, 'r') as fi:
# 		for line in fi:
# 			tokens = line.strip().split('\t')

# 			feature.append([cell2id[tokens[0]], drug2id[tokens[1]]])
# 			label.append([float(tokens[2])])

# 	return feature, label

def load_train_data(file_name, cell2id):
	feature = []
	label = []

	with open(file_name, 'r') as fi:
		for line in fi:
			tokens = line.strip().split('\t')

			feature.append([cell2id[tokens[0]]])
			label.append([float(tokens[1])])
			
	return feature, label


# def prepare_predict_data(test_file, cell2id_mapping_file, drug2id_mapping_file):

# 	# load mapping files
# 	cell2id_mapping = load_mapping(cell2id_mapping_file)
# 	drug2id_mapping = load_mapping(drug2id_mapping_file)

# 	test_feature, test_label = load_train_data(test_file, cell2id_mapping, drug2id_mapping)

# 	print('Total number of cell lines = %d' % len(cell2id_mapping))
# 	print('Total number of drugs = %d' % len(drug2id_mapping))

# 	return (torch.Tensor(test_feature), torch.Tensor(test_label)), cell2id_mapping, drug2id_mapping

def prepare_predict_data(test_file, cell2id_mapping_file):

	# load mapping files
	cell2id_mapping = load_mapping(cell2id_mapping_file)

	test_feature, test_label = load_train_data(test_file, cell2id_mapping)

	print('Total number of cell lines = %d' % len(cell2id_mapping))

	return (torch.Tensor(test_feature), torch.Tensor(test_label)), cell2id_mapping


def load_mapping(mapping_file):

	mapping = {}

	file_handle = open(mapping_file)

	for line in file_handle:
		line = line.rstrip().split("\t")# added .split("\t")
		mapping[line[1]] = int(line[0])

	file_handle.close()
	
	return mapping

def load_sample2data_mapping(mapping_file):
	mapping = {}

	with open(mapping_file, "r") as f:
		lines = f.readlines()
		for line in lines:
			tokens = line.rstrip().split("\t")
			mapping[tokens[0]] = tokens[1]
	
	return mapping


# def prepare_train_data(train_file, test_file, cell2id_mapping_file, drug2id_mapping_file):

# 	# load mapping files
# 	cell2id_mapping = load_mapping(cell2id_mapping_file)
# 	drug2id_mapping = load_mapping(drug2id_mapping_file)

# 	train_feature, train_label = load_train_data(train_file, cell2id_mapping, drug2id_mapping)
# 	test_feature, test_label = load_train_data(test_file, cell2id_mapping, drug2id_mapping)

# 	print('Total number of cell lines = %d' % len(cell2id_mapping))
# 	print('Total number of drugs = %d' % len(drug2id_mapping))

# 	return (torch.Tensor(train_feature), torch.FloatTensor(train_label), torch.Tensor(test_feature), torch.FloatTensor(test_label)), cell2id_mapping, drug2id_mapping

def prepare_train_data(train_file, test_file, cell2id_mapping_file):

	# load mapping files
	cell2id_mapping = load_mapping(cell2id_mapping_file)

	train_feature, train_label = load_train_data(train_file, cell2id_mapping)
	test_feature, test_label = load_train_data(test_file, cell2id_mapping)

	print('Total number of cell lines = %d' % len(cell2id_mapping))

	return (torch.Tensor(train_feature), torch.FloatTensor(train_label), torch.Tensor(test_feature), torch.FloatTensor(test_label)), cell2id_mapping


# def build_input_vector(input_data, cell_features, drug_features):
# 	genedim = len(cell_features[0,:])
# 	drugdim = len(drug_features[0,:])
# 	feature = np.zeros((input_data.size()[0], (genedim+drugdim)))

# 	for i in range(input_data.size()[0]):
# 		feature[i] = np.concatenate((cell_features[int(input_data[i,0])], drug_features[int(input_data[i,1])]), axis=None)

# 	feature = torch.from_numpy(feature).float()
# 	return feature

def build_input_vector(input_data, cell_features):
	# print(f"build_input_vector input_data: {input_data}")
	genedim = len(cell_features[0,:])
	feature = np.zeros((input_data.size()[0], (genedim)))

	for i in range(input_data.size()[0]):
		# print(f"trying to do: {cell_features[int(input_data[i,0])]}")
		feature[i] = cell_features[int(input_data[i,0])]

	feature = torch.from_numpy(feature).float()
	return feature