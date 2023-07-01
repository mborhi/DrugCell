# Change Log

## Remove Drug ANN, Change to Classification

The drug ANN branch is removed. Final layer only considers gene ontology VNN branch.

Model classifies each input as either being Invasive Breast Carcinoma or Breast Invasive Lobular Carcinoma. The output layer is two neurons, each representing one of the classes. Uses Cross Entropy loss. 

## Data

Data sets:

* Colorectal:
    - [cBioPortal Dashboard](https://www.cbioportal.org/study/summary?id=appendiceal_msk_2022%2Ccoad_caseccc_2015%2Ccoad_cptac_2019%2Ccoadread_dfci_2016%2Ccoadread_genentech%2Cbowel_colitis_msk_2022%2Ccoadread_tcga%2Ccoadread_tcga_pub%2Ccoadread_tcga_pan_can_atlas_2018%2Ccoadread_mskcc%2Ccoadread_mskresistance_2022%2Ccrc_apc_impact_2020%2Ccrc_dd_2022%2Ccrc_nigerian_2020%2Ccrc_msk_2017%2Crectal_msk_2022%2Crectal_msk_2019)
    - Types: Colon Adenocarcinoma, Rectal Adenocarcinoma
* Lymphoid:
    - [cBioPortal Dashboard](https://www.cbioportal.org/study?id=all_stjude_2015%2Call_stjude_2016%2Clcll_broad_2013%2Ccll_broad_2015%2Ccll_iuopa_2015%2Ccllsll_icgc_2011%2Cctcl_columbia_2015%2Cdlbcl_dfci_2018%2Cdlbc_broad_2012%2Cdlbcl_duke_2017%2Cdlbc_tcga_pan_can_atlas_2018%2Cnhl_bcgsc_2013%2Cdlbc_tcga%2Clymphoma_cellline_msk_2020%2Cmcl_idibips_2013%2Cmbn_mdacc_2013%2Cmtnn_msk_2022%2Cmm_broad%2Cnhl_bcgsc_2011%2Call_phase2_target_2018_pub%2Cpcnsl_mayo_2015)
    - Types: Diffuse Large B-Cell Lymphoma - NOS, Chronic Lymphocytic Leukemia/Small Lymphocytic Lymphoma
* Lung:
    - [cBioPortal Dashboard](https://www.cbioportal.org/study/summary?id=luad_broad%2Cluad_cptac_2020%2Clung_msk_mind_2020%2Cluad_mskimpact_2021%2Cluad_mskcc_2020%2Cluad_msk_npjpo_2021%2Cluad_mskcc_2015%2Cluad_oncosg_2020%2Cluad_tcga%2Cluad_tcga_pub%2Cluad_tcga_pan_can_atlas_2018%2Cluad_tsp%2Clung_smc_2016%2Clung_nci_2022%2Clusc_cptac_2021%2Clusc_tcga%2Clusc_tcga_pub%2Clusc_tcga_pan_can_atlas_2018%2Cnsclc_ctdx_msk_2022%2Clung_msk_2017%2Cnsclc_mskcc_2018%2Cnsclc_pd1_msk_2018%2Cnsclc_mskcc_2015%2Cnsclc_tracerx_2017%2Cnsclc_unito_2016%2Cnsclc_tcga_broad_2016%2Csclc_clcgp%2Csclc_jhu%2Csclc_ucologne_2015%2Csclc_cancercell_gardner_2017%2Clung_pdx_msk_2021%2Clung_msk_pdx)
    - Types: Lung Adenocarcinoma, Lung Squamous Cell Carcinoma
* Brain:
    - [cBioPortal Dashboard](https://www.cbioportal.org/study/summary?id=odg_msk_2017%2Clgg_tcga%2Clgg_tcga_pan_can_atlas_2018%2Cgbm_mayo_pdx_sarkaria_2019%2Cdifg_glass_2019%2Cgbm_cptac_2021%2Cgbm_columbia_2019%2Cgbm_tcga_pub2013%2Cgbm_tcga_pub%2Cgbm_tcga%2Cgbm_tcga_pan_can_atlas_2018%2Cglioma_mskcc_2019%2Cglioma_msk_2018%2Clgg_ucsf_2014%2Cmbl_broad_2012%2Cmbl_dkfz_2017%2Cmbl_icgc%2Cmbl_pcgp%2Cmbl_sickkids_2016%2Cmng_utoronto_2021%2Clgggbm_tcga_pub%2Cbrain_cptac_2020%2Cpcpg_tcga%2Cpast_dkfz_heidelberg_2013)
    - Types: Glioblastoma Multiforme, Oligodendroglioma
* Pediatric:
    [cBioPortal Dashboard](https://www.cbioportal.org/study/summary?id=all_stjude_2015%2Call_stjude_2016%2Ces_iocurie_2014%2Cmbl_pcgp%2Call_phase2_target_2018_pub%2Caml_target_2018_pub%2Ces_dfarber_broad_2014%2Cnbl_target_2018_pub%2Cpediatric_dkfz_2017%2Cmixed_pipseq_2017%2Cpptc_2019%2Crt_target_2018_pub%2Cwt_target_2018_pub)
    - Types: Neuroblastoma, Acute Myeloid Leukemia
* Stomach, Esophageal:
    - [cBioPortal Dashboard](https://www.cbioportal.org/study/summary?id=esca_broad%2Cesca_tcga_pan_can_atlas_2018%2Cegc_trap_msk_2020%2Cesca_tcga%2Cstes_tcga_pub%2Cescc_icgc%2Cescc_ucla_2014%2Cegc_mskcc_2020%2Cegc_msk_tp53_ccr_2022%2Cegc_tmucih_2015%2Cstad_oncosg_2018%2Cegc_msk_2017%2Cstad_pfizer_uhongkong%2Cstad_tcga%2Cstad_tcga_pub%2Cstad_tcga_pan_can_atlas_2018%2Cstad_utokyo%2Cstad_uhongkong)
    - Types: Stomach Adenocarcinoma, Esophageal Adenocarcinoma

Data set: [cBioPortal, Breast Invasive Ductal Carcinoma and Breast Invasive Lobular Carcinoma](https://www.cbioportal.org/study/summary?id=acbc_mskcc_2015%2Cbrca_hta9_htan_2022%2Cbrca_metabric%2Cbreast_msk_2018%2Cbrca_pareja_msk_2020%2Cbrca_mskcc_2019%2Cbreast_alpelisib_2020%2Cbrca_smc_2018%2Cbrca_bccrc_xenograft_2014%2Cbfn_duke_nus_2015%2Cbrca_bccrc%2Cbrca_broad%2Cbrca_sanger%2Cbrca_tcga_pub2015%2Cbrca_tcga%2Cbrca_tcga_pub%2Cbrca_tcga_pan_can_atlas_2018%2Cbrca_jup_msk_2020%2Cbrca_mapk_hp_msk_2021%2Cmbc_msk_2021%2Cbrca_igr_2015%2Cbreast_ink4_msk_2021%2Cbrca_cptac_2020%2Cbrca_mbcproject_wagle_2017%2Cbrca_mbcproject_2022)

Combined study TSV file: [combined_study_clinical_ILC_IDC.tsv](./custom_data/combined_study_clinical_ILC_IDC.tsv)

`sample2type_ILC_IDC.txt` - tab-delimited file.
* first column: the id of the cell line: `[study_id]:[sample_id]`
* second column: the type of cancer of the cell line, `0` if `Breast Invasive Ductal Carcinoma`, and `1` if `Breast Invasive Lobular Carcinoma`

`sample2type_train_ILC_IDC.txt` - same as `sample2type_ILC_IDC.txt` but used for training (80%)

`sample2type_test_ILC_IDC.txt` - same as `sample2type_ILC_IDC.txt` but used for testing (20%)

`sample2mutation_ILC_IDC.txt` - tab-delimited file.
* first column: the id of the cell line `[study_id]:[sample_id]`
* second column: the mutated genes in the cell line delimited by commas

### Reactome

Added all Reactome human pathway data.

To generate the necessary data to train a model using a Reactome architecture:
1. Run [`reactome.sh`](./scripts/data_collection/reactome.sh), found in the `./scripts/data_collection/` dir.
2. Run [`write_reactome_data.sh`](./scripts/data_collection/write_reactome_data.sh), found in the `./scripts/data_collection/` dir.

Model can be trained using the [`commandline_train_reactome_classification.sh`](scripts/commandline_train_reactome_classification.sh) script.


Model can be tested using the [`scripts/commandline_test_reactome_classification.sh`](scripts/commandline_test_reactome_classification.sh) script.