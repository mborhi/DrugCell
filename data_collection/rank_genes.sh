input_dir="../classification_data/"
output_dir="gene_data"
sample2mutations=$input_dir"sample2mutations_ILC_IDC.txt"

mkdir $output_dir

source activate pytorch3drugcell

python -u ./rank_genes.py -input_dir=$input_dir -output_dir=$output_dir -sample2mutations=$sample2mutations