#########

## Go to TIGAR tool directory
cd /home/jyang/GIT/TIGAR

TIGAR_dir="/home/jyang/GIT/TIGAR"
Gene_Exp_train_file="${TIGAR_dir}/ExampleData/gene_exp.txt"
# All training sampleID need to exist in gene expression file
train_sample_ID_file="${TIGAR_dir}/ExampleData/sampleID.txt"
genofile="${TIGAR_dir}/ExampleData/example.vcf.gz"
out_dir="${TIGAR_dir}/ExampleData/output"


mafval=0.01
hweval=0.00001

### Load Anaconda3 module, set up PATHONPATH to include required libraries
module load Anaconda3/4.2.0
conda activate myenv
export PYTHONPATH=/home/rparrish/.conda/envs/myenv/lib/python3.5/site-packages/:$PYTHONPATH

### Train DPR model
${TIGAR_dir}/TIGAR_Model_Train.sh --model DPR \
--Gene_Exp ${Gene_Exp_train_file} \
--train_sampleID ${train_sample_ID_file} \
--genofile ${genofile} --chr 1 \
--genofile_type vcf --Format GT \
--maf ${mafval} \
--hwe ${hweval} \
--cvR2 1 \
--dpr 1 --ES fixed \
--TIGAR_dir ${TIGAR_dir} \
--out_dir ${out_dir} &

### Train elastic_net model
${TIGAR_dir}/TIGAR_Model_Train.sh --model elastic_net \
--Gene_Exp ${Gene_Exp_train_file} \
--train_sampleID ${train_sample_ID_file} \
--genofile ${genofile} --chr 1 \
--genofile_type vcf --Format GT \
--maf ${mafval} \
--hwe ${hweval} \
--cvR2 1 \
--alpha 0.5 \
--TIGAR_dir ${TIGAR_dir} \
--out_dir ${out_dir} 

### Predict GReX for the same group of samples used for model training
eQTL_ES_file="${TIGAR_dir}/ExampleData/eQTLweights.txt"
gene_anno_file="${TIGAR_dir}/ExampleData/gene_anno.txt"
# All test sample ID need to exist in genotype file
test_sample_ID_file="${TIGAR_dir}/ExampleData/test_sampleID.txt"
genofile="${TIGAR_dir}/ExampleData/example.vcf.gz"
out_dir="${TIGAR_dir}/ExampleData/output"

${TIGAR_dir}/TIGAR_Model_Pred.sh --chr 1 \
--weight ${eQTL_ES_file} \
--gene_anno ${gene_anno_file} \
--test_sampleID ${test_sample_ID_file} \
--genofile ${genofile} \
--genofile_type vcf --Format GT \
--TIGAR_dir ${TIGAR_dir} \
--out_dir ${out_dir}

### TWAS using individual-level data
GReX_pred_file="${TIGAR_dir}/CHR1_Pred_GReX.txt"
PED="${TIGAR_dir}/example_PED.ped"
PED_info="${TIGAR_dir}/PED_Info_SinglePheno.txt"
out_dir="${TIGAR_dir}/ExampleData/output"

${TIGAR_dir}/TIGAR_TWAS.sh --asso 1 \
--gene_exp ${GReX_pred_file} \
--PED ${PED} --PED_info ${PED_info} \
--method OLS \
--out ${out_dir}



### TWAS using summary-level data
gene_anno_file="${TIGAR_dir}/ExampleData/gene_anno.txt"
Zscore_file="${TIGAR_dir}/CHR1_GWAS_Zscore.txt.gz"
eQTL_ES_file="${TIGAR_dir}/ExampleData/eQTLweights.txt"
LD_file="${TIGAR_dir}/CHR1_LD.txt.gz"
out_dir="${TIGAR_dir}/ExampleData/output"

${TIGAR_dir}/TIGAR_TWAS.sh --asso 2 \
--gene_anno ${gene_anno_file} \
--Zscore ${Zscore_file} \
--weight ${eQTL_ES_file} \
--LD ${LD_file} \
--chr 1 \
--out ${out_dir}



