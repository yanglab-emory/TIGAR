#########

## Go to TIGAR tool directory
cd /home/jyang/GIT/TIGAR

Gene_Exp_train_file="./ExampleData/gene_exp.txt"
# All training sampleID need to exist in gene expression file
train_sample_ID_file="./ExampleData/sampleID.txt"
genofile="./ExampleData/example.vcf.gz"
out_dir="./ExampleData/output"

mafval=0.01
hweval=0.00001

### Load Anaconda3 module, set up PATHONPATH to include required libraries
module load Anaconda3/4.2.0
conda activate myenv
export PYTHONPATH=/home/rparrish/.conda/envs/myenv/lib/python3.5/site-packages/:$PYTHONPATH

### Train DPR model
./TIGAR_Model_Train.sh --model DPR \
--Gene_Exp ${Gene_Exp_train_file} \
--train_sampleID ${train_sample_ID_file} \
--genofile ${genofile} --chr 1 \
--genofile_type vcf --Format GT \
--maf ${mafval} \
--hwe ${hweval} \
--cvR2 1 \
--dpr 1 --ES fixed \
--out_dir ${out_dir} &

### Train elastic_net model
./TIGAR_Model_Train.sh --model elastic_net \
--Gene_Exp ${Gene_Exp_train_file} \
--train_sampleID ${train_sample_ID_file} \
--genofile ${genofile} --chr 1 \
--genofile_type vcf --Format GT \
--maf ${mafval} \
--hwe ${hweval} \
--cvR2 1 \
--alpha 0.5 \
--out_dir ${out_dir} 

### Predict GReX for the same group of samples used for model training
eQTL_ES_file="./ExampleData/eQTLweights.txt"
gene_anno_file="./ExampleData/gene_anno.txt"
# All test sample ID need to exist in genotype file
test_sample_ID_file="./ExampleData/test_sampleID.txt"
genofile="./ExampleData/example.vcf.gz"
out_dir="./ExampleData/output"

./TIGAR_Model_Pred.sh --chr 1 \
--weight ${eQTL_ES_file} \
--gene_anno ${gene_anno_file} \
--test_sampleID ${test_sample_ID_file} \
--genofile ${genofile} \
--genofile_type vcf --Format GT \
--out_dir ${out_dir}


