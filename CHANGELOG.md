# Change Log

## Remove Drug ANN, Change to Classification

The drug ANN branch is removed. Final layer only considers gene ontology VNN branch.

Model classifies each input as either being a Invasive Breast Carcinoma or Breast Invasive Lobular Carcinoma. The final output is a single neuron using sigmoid, representing the probability of the input belonging to the two classes.

## Data

Data set: cBioPortal, invasive breast carcinoma, and breast invasive lobular carcinoma.

Accessible here: [Data](https://www.cbioportal.org/study/summary?id=acbc_mskcc_2015%2Cbrca_hta9_htan_2022%2Cbrca_metabric%2Cbreast_msk_2018%2Cbrca_pareja_msk_2020%2Cbrca_mskcc_2019%2Cbreast_alpelisib_2020%2Cbrca_smc_2018%2Cbrca_bccrc_xenograft_2014%2Cbfn_duke_nus_2015%2Cbrca_bccrc%2Cbrca_broad%2Cbrca_sanger%2Cbrca_tcga_pub2015%2Cbrca_tcga%2Cbrca_tcga_pub%2Cbrca_tcga_pan_can_atlas_2018%2Cbrca_jup_msk_2020%2Cbrca_mapk_hp_msk_2021%2Cmbc_msk_2021%2Cbrca_igr_2015%2Cbreast_ink4_msk_2021%2Cbrca_cptac_2020%2Cbrca_mbcproject_wagle_2017%2Cbrca_mbcproject_2022)

Combined study and clinical data for all samples: [TSV file](https://www.cbioportal.org/study/clinicalData?id=acbc_mskcc_2015%2Cbrca_hta9_htan_2022%2Cbrca_metabric%2Cbreast_msk_2018%2Cbrca_pareja_msk_2020%2Cbrca_mskcc_2019%2Cbreast_alpelisib_2020%2Cbrca_smc_2018%2Cbrca_bccrc_xenograft_2014%2Cbfn_duke_nus_2015%2Cbrca_bccrc%2Cbrca_broad%2Cbrca_sanger%2Cbrca_tcga_pub2015%2Cbrca_tcga%2Cbrca_tcga_pub%2Cbrca_tcga_pan_can_atlas_2018%2Cbrca_jup_msk_2020%2Cbrca_mapk_hp_msk_2021%2Cmbc_msk_2021%2Cbrca_igr_2015%2Cbreast_ink4_msk_2021%2Cbrca_cptac_2020%2Cbrca_mbcproject_wagle_2017%2Cbrca_mbcproject_2022)

`sample2type.txt` - tab-delimited file.
* first column: the id of the cell line: `[study_id]:[sample_id]`
* second column: the type of cancer of the cell line, `0` if `Invasive Breast Carcinoma`, and `1` if `Breast Invasive Lobular Carcinoma`

`sample2type_train.txt` - same as `sample2type.txt` but used for training (80%)

`sample2type_est.txt` - same as `sample2type.txt` but used for testing (20%)

`sample2mutation.txt` - tab-delimited file.
* first column: the id of the cell line `[study_id]:[sample_id]`
* second column: the mutated genes in the cell line delimited by commas