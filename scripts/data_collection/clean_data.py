import numpy as np
import requests
import argparse

# classifications
CLASS_A = "Breast Invasive Ductal Carcinoma"
CLASS_B = "Breast Invasive Lobular Carcinoma"
# API base url
BASE_URL = "https://www.cbioportal.org/api"


def get_molecular_profiles(study_id):
    mpids = set()

    mpfs_req = requests.get(BASE_URL + f"/studies/{study_id}/molecular-profiles", params={
        "projection": "DETAILED"
    })
    mpfs = mpfs_req.json()
    for mpf in mpfs:
        mpids.add(mpf["molecularProfileId"])
    
    return mpids
    
def get_all_sample_gene_mutuations(mpids_list, study_sample_ids):
    study_sample_to_mutations = {}

    url = BASE_URL + "/mutations/fetch"
    mutations_data_list_req = requests.post(url,
        headers={
            "accept": "application/json",
            "Content-Type": "application/json"
        },
        params={"projection": "DETAILED"},
        json={"molecularProfileIds": mpids_list}
    )
    if mutations_data_list_req.status_code == 404:
        print(mutations_data_list_req.json())
    mutations_data_list = mutations_data_list_req.json()
    # print("mutations data list", mutations_data_list)
    for mutation_data in mutations_data_list:
        study_id = mutation_data["studyId"]
        sample_id = mutation_data["sampleId"]
        study_sample_id = study_id + ":" + sample_id
        if study_sample_id not in study_sample_ids:
            continue
        if study_sample_id in study_sample_to_mutations:
            study_sample_to_mutations[study_sample_id].add(mutation_data["gene"]["hugoGeneSymbol"])
        else :
            new_gene_set = set()
            new_gene_set.add(mutation_data["gene"]["hugoGeneSymbol"])
            study_sample_to_mutations[study_sample_id] = new_gene_set
    
    return study_sample_to_mutations
        

def collect_data(filename):
    study_sample_to_type = {}
    sample_to_gene_mutations = {}
    study_id_to_mpid = {}

    all_sample_ids = []
    all_study_ids = []
    all_patient_ids = []

    with open(filename, "r") as f:
        lines = f.readlines()
        cols = lines[0].split("\t")
        for i, line in enumerate(lines[1:]):
            split_line = line.split("\t")
            study_id = split_line[cols.index("Cancer Study")] # sometimes "Study ID"
            patient_id = split_line[cols.index("Patient ID")]
            sample_id = split_line[cols.index("Sample ID")]
            cancer_type = split_line[cols.index("Cancer Type Detailed")]
            mutation_count = split_line[cols.index("Mutation Count")]

            new_id = study_id + ":" + sample_id
            if (mutation_count.isnumeric() and int(mutation_count) > 0):
                study_sample_to_type[new_id] = 0 if cancer_type == CLASS_A else 1
            
            all_sample_ids.append(sample_id)
            all_patient_ids.append(patient_id)
            all_study_ids.append(study_id)

    # all_mpids = set()
    # for sid in all_study_ids:
    #     all_mpids = all_mpids.union(get_molecular_profiles(sid))
    # #  print("all molecular profile ids:", list(all_mpids))
    # all_mpids_list = list(all_mpids)
    all_mpids_list =  ['brca_mbcproject_2022_rna_seq_mrna_median_Zscores', 'brca_tcga_methylation_hm450', 
                       'brca_tcga_pan_can_atlas_2018_rppa_Zscores', 'brca_tcga_pub2015_mrna_median_Zscores', 
                       'breast_ink4_msk_2021_structural_variants', 'brca_mapk_hp_msk_2021_cna', 
                       'brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_all_sample_Zscores', 
                       'brca_tcga_pub2015_protein_quantification_zscores', 'brca_tcga_pub2015_rna_seq_v2_mrna', 
                       'brca_metabric_cna', 'brca_tcga_protein_quantification_zscores', 'breast_alpelisib_2020_cna', 
                       'brca_tcga_rna_seq_v2_mrna_median_Zscores', 'brca_mapk_hp_msk_2021_mutations', 
                       'brca_tcga_pub2015_mutations', 'brca_mbcproject_2022_mutations', 
                       'brca_tcga_pub2015_mrna_median_all_sample_Zscores', 'brca_tcga_pub2015_rna_seq_v2_mrna_median_all_sample_Zscores', 'brca_tcga_protein_quantification', 'breast_alpelisib_2020_structural_variants', 'brca_metabric_mrna', 'breast_msk_2018_mutations', 'brca_metabric_methylation_promoters_rrbs', 'brca_tcga_rna_seq_v2_mrna', 'brca_tcga_pub2015_rna_seq_v2_mrna_median_Zscores', 'brca_mbcproject_wagle_2017_mutations', 'brca_tcga_pan_can_atlas_2018_log2CNA', 'brca_tcga_pan_can_atlas_2018_gistic', 'brca_tcga_pan_can_atlas_2018_protein_quantification', 'breast_msk_2018_structural_variants', 'brca_tcga_linear_CNA', 'brca_tcga_pub2015_methylation_hm27', 'brca_tcga_pan_can_atlas_2018_structural_variants', 'brca_tcga_gistic', 'brca_tcga_pub2015_linear_CNA', 'brca_tcga_rna_seq_v2_mrna_median_all_sample_Zscores', 'brca_mbcproject_wagle_2017_rna_seq_v2_mrna', 'brca_hta9_htan_2022_cna', 'brca_tcga_mrna', 'brca_tcga_rppa_Zscores', 'brca_tcga_pub2015_methylation_hm450', 'brca_hta9_htan_2022_mutations', 'brca_tcga_pub2015_mrna', 'brca_tcga_phosphoprotein_quantification', 'brca_tcga_mutations', 'brca_smc_2018_mrna_seq_tpm_all_sample_Zscores', 'brca_mbcproject_2022_rna_seq_v2_mrna', 'breast_ink4_msk_2021_cna', 'brca_mbcproject_2022_structural_variants', 'brca_tcga_mrna_median_all_sample_Zscores', 'brca_smc_2018_mutations', 'breast_alpelisib_2020_mutations', 'brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna', 'brca_tcga_pan_can_atlas_2018_methylation_hm27_hm450_merge', 'brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_Zscores', 'brca_tcga_pan_can_atlas_2018_protein_quantification_zscores', 'brca_tcga_mrna_median_Zscores', 'brca_tcga_pan_can_atlas_2018_phosphoprotein_quantification', 'brca_smc_2018_mrna_seq_tpm', 'brca_tcga_pan_can_atlas_2018_mutations', 'breast_msk_2018_cna', 'brca_tcga_pan_can_atlas_2018_rna_seq_v2_mrna_median_all_sample_ref_normal_Zscores', 'breast_ink4_msk_2021_mutations', 'brca_tcga_pub2015_gistic', 'brca_metabric_mrna_median_all_sample_Zscores', 'brca_tcga_pub2015_protein_quantification', 'brca_tcga_pan_can_atlas_2018_armlevel_cna', 'brca_tcga_rppa', 'brca_mbcproject_2022_gistic', 'brca_tcga_pan_can_atlas_2018_microbiome_signature', 'brca_mbcproject_wagle_2017_gistic', 'brca_mapk_hp_msk_2021_structural_variants', 'brca_mbcproject_wagle_2017_rna_seq_mrna_median_Zscores', 'brca_tcga_pan_can_atlas_2018_rppa', 'brca_metabric_mutations']
    sample_to_mutations = get_all_sample_gene_mutuations(all_mpids_list, study_sample_to_type.keys())
    
    return study_sample_to_type, sample_to_mutations#sample_to_gene_mutations

def write_sample_to_type(filename, sample_to_type):
    with open(filename, "w") as f:
        for sample_id in sample_to_type.keys():
            line = sample_id + "\t" + str(sample_to_type[sample_id]) + "\n"
            f.write(line)

def write_sample_to_gene_mutations(filename, sample_to_gene_mutations):
    with open(filename, "w") as f:
        for sample_id in sample_to_type.keys():
            mutations = ",".join(list(sample_to_gene_mutations[sample_id]))
            line = sample_id + "\t" + mutations + "\n"
            f.write(line)

def write_test_train(filename, keys, dict):
    with open(filename, "w") as f:
        for key in keys:
            values = str(dict[key])
            line = key + "\t" + values + "\n"
            f.write(line)

def prepare_data(dict):
    keys = list(dict.keys())
    np.random.shuffle(keys)
    # split data
    train = keys[:int(len(keys) * 0.85)]
    test = keys[int(len(keys) * 0.85):]

    return train, test

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Clean data')
    parser.add_argument('-output_dir', help='The directory to write new data files to', type=str)
    parser.add_argument('-combined_file', help='The file containing all combined patient data', type=str)
    opt = parser.parse_args()

    output_dir = opt.output_dir
    sample_to_type, sample_to_gene_mutations = collect_data(opt.combined_file)
    
    sample_to_type_train, sample_to_type_test = prepare_data(sample_to_type)

    write_sample_to_type(output_dir + "sample2type_ILC_IDC.txt", sample_to_type)
    write_sample_to_gene_mutations(output_dir + "sample2mutations_ILC_IDC.txt", sample_to_gene_mutations)

    write_test_train(output_dir + "sample2type_train_ILC_IDC.txt", sample_to_type_train, sample_to_type)
    write_test_train(output_dir + "sample2type_test_ILC_IDC.txt", sample_to_type_test, sample_to_type)
