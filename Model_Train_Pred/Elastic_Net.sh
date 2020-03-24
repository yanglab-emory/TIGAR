#!/usr/bin/bash

#####################################################################################################
# vairable needed for training
###
# model: elastic_net
# Gene_Exp: Gene annotation and Expression level file path
# train_sample: a column of sampleIDs use for training
# chr: chromosome number for corresponding training data
# thread: number of thread for multiprocessing
# geno_train: vcf or dosages
# train_dir: training data path, should be tabix(contains .gz and .tbi)
# FT: Format using for training data(GT or DS), default DS
# maf: Threshold for Minor Allele Frequency (range from 0-1),default 0.01
# hwe: Threshold of p-value for Hardy Weinberg Equilibrium exact test,default 0.001
# window: window for selecting data
# cv: cv-fold cross-validation in model selection, default 5-fold
# alpha: L1 & L2 ratio for elastic net regression, default 0.5
#        If alpha=0, lasso regression
#        If alpha=1, ridge regression
# out: output dir

# variable needed for prediction
# pred: y or n (default is n)
#       y means user wants to conduct prediction
#       n means user do not want to run prediction part
# geno_test: vcf or dosages
# test_dir: testing data path, should be tabix(contains .gz and .tbi)
# FP: Format using for testing data (GT or DS), default DS
# maf_diff: threshold of difference between training maf and testing maf, if difference is larger than this, 
#           then drop this data, default 0.2  

VARS=`getopt -o "" -a -l \
model:,Gene_Exp:,train_sample:,chr:,genofile_type:,genofile_dir:,Format:,maf:,hwe:,window:,cv:,alpha:,thread:,out: \
-- "$@"`

if [ $? != 0 ]
then
    echo "Terminating....." >&2
    exit 1
fi
 
eval set -- "$VARS"

while true
do
    case "$1" in
        --model|-model) model=$2; shift 2;;
        --Gene_Exp|-Gene_Exp) Gene_Exp=$2; shift 2;;
        --train_sample|-train_sample) train_sample=$2; shift 2;;
        --chr|-chr) chr_num=$2; shift 2;;
        --genofile_type|-genofile_type) genofile_type=$2; shift 2;;
        --genofile_dir|-genofile_dir) genofile_dir=$2; shift 2;;
        --Format|-Format) Format=$2; shift 2;;
        --maf|-maf) maf=$2; shift 2;;
        --hwe|-hwe) hwe=$2; shift 2;;
        --window|-window) window=$2; shift 2;;
        --cv|-cv) cv=$2; shift 2;;
        --alpha|-alpha) alpha=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --out|-out) out_prefix=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done

######################################################################################################
### 1. 
### Create dir & Store Result
mkdir -p ${out_prefix}/${model}_CHR${chr_num}

### Create dir to store log file from python
mkdir -p ${out_prefix}/${model}_CHR${chr_num}/log_file

### 2. 
### Extract vcf header from vcf file
zcat ${genofile_dir} | grep 'CHROM' > ${out_prefix}/CHR${chr_num}_train_names.txt

### 3.
### Training
python ./Model_Train_Pred/Elastic_Net_Train.py \
--Gene_Exp_path ${Gene_Exp} \
--train_sample ${train_sample} \
--chr_num ${chr_num} \
--thread ${thread} \
--train_dir ${genofile_dir} \
--train_names ${out_prefix}/CHR${chr_num}_train_names.txt \
--geno ${genofile_type} \
--Format ${Format} \
--maf ${maf} \
--hwe ${hwe} \
--window ${window} \
--cv ${cv} \
--alpha ${alpha} \
--out_prefix ${out_prefix}/${model}_CHR${chr_num} \
> ${out_prefix}/${model}_CHR${chr_num}/log_file/Elastic_Net_Train_log.txt

### 4.
### Remove file
rm ${out_prefix}/CHR${chr_num}_train_names.txt









