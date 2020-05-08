#########

## Go to TIGAR tool directory
cd /home/jyang/GIT/TIGAR

Gene_Exp_train_file="./example_data/Gene_Exp.txt"
train_sample_ID_file="./example_data/sampleID.txt"
genofile="./example_data/example.vcf.gz"
out_dir="./example_data/output"

mafval=0.01
hweval=0.00001

### Load Anaconda3 module, set up PATHONPATH to include required libraries
module load Anaconda3/4.2.0
conda activate myenv
export PYTHONPATH=/home/rparrish/.conda/envs/myenv/lib/python3.5/site-packages/:$PYTHONPATH

./TIGAR_Model_Train.sh --model DPR \
--Gene_Exp ${Gene_Exp_train_file} \
--train_sampleID ${train_sample_ID_file} \
--genofile ${genofile} --chr 1 \
--genofile_type vcf --Format GT \
--maf ${mafval} \
--hwe ${hweval} \
--out_dir ${out_dir}



############### Randy scripts #########
DATADIR=/mnt/YangFSS/data/rparrish/GTEx_V8

Gene_Exp_train_file=${DATADIR}/ExpressionFiles/${tissue_name}_GTEx_Exp.txt
train_sample_path=${DATADIR}/SubjectIDFiles/${tissue_name}_subjid.txt
geno_dir=${DATADIR}/GenotypeFiles/CHR${chromnum}_GTEx_WGS.vcf.gz
genoformat=GT
model=DPR
mafval=0.01
# 10e^-6
hweval=0.00001

OUTDIR=${DATADIR}/TIGAR_train/${tissue_name}
mkdir -p ${OUTDIR}

# LOAD MODULES
module load tabix/0.2.6
conda activate myenv

# 02/27/20
DPR=/home/rparrish/TIGAR/Model_Train_Pred/DPR
# DPR=/home/rparrish/DPR_Model/DPR
# export PYTHONPATH=/home/xmeng/.local/lib/python3.5/site-packages/:/home/rparrish/.local/lib/python3.5/site-packages/
export PYTHONPATH=/home/rparrish/.conda/envs/myenv/lib/python3.5/site-packages/
# /home/rparrish/.local/lib/python3.5/site-packages/
echo Job ID: "${JOB_ID}"
echo Running TIGAR train with "${model}" model for "$ind"th tissue, "${tissue_name}", on chromosome "${chromnum}" with "${NSLOTS}" cores

start_time=`date +%s`
echo Start time: "${start_time}"

TIGAR_Model_Train.sh \
--model ${model} \
--Gene_Exp ${Gene_Exp_train_file} \
--train_sample ${train_sample_path} \
--chr ${chromnum} \
--thread ${NSLOTS} \
--genofile_type vcf \
--genofile_dir ${geno_dir} \
--Format ${genoformat} \
--maf ${mafval} \
--hwe ${hweval} \
--out ${OUTDIR}
