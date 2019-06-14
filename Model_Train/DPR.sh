#!/usr/bin/bash

#############################################################################################
### vairable needed for training
# model: DPR
# Gene_Exp: Gene annotation and Expression level file
# train_sample: a column of sampleIDs use for training
# chr: chromosome number for corresponding training data
# thread: number of thread for multiprocessing
# geno_train: vcf or dosages, both should be tabix
# train_dir: training data path
# FT: Format using for training data(GT or DS), default DS
# maf: Threshold for Minor Allele Frequency (range from 0-1),default 0.01
# hwe: Threshold of p-value for Hardy Weinberg Equilibrium exact test,default 0.001
# window: window for selecting data
# dpr: model when we choose DPR, default 1
# ES: effect size(fixed/additive)
# out: output dir

### variable needed for prediction
# pred: y or n (default is n)
#       y means user wants to conduct prediction
#       n means user do not want to run prediction part
# geno_test: vcf or dosages, both should be tabix
# test_dir: testing data path
# FP: Format using for testing data (GT or DS), default DS
# maf_diff: threshold of difference between training maf and testing maf, if difference is larger than this, 
#           then drop this data, default 0.2  


VARS=`getopt -o "" -a -l \
model:,Gene_Exp:,train_sample:,chr:,thread:,geno_train:,train_dir:,FT:,maf:,hwe:,window:,dpr:,ES:,pred:,test_dir:,geno_test:,FP:,maf_diff:,out: \
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
        --thread|-thread) thread=$2; shift 2;;
        --geno_train|-geno_train) geno_train=$2; shift 2;;
        --train_dir|-train_dir) train_dir=$2; shift 2;;
        --FT|-FT) FT=$2; shift 2;;
        --maf|-maf) maf=$2; shift 2;;
        --hwe|-hwe) hwe=$2; shift 2;;
        --window|-window) window=$2; shift 2;;
        --dpr|-dpr) dpr_num=$2; shift 2;;
        --ES|-ES) ES=$2; shift 2;;
        --pred|-pred) pred=$2; shift 2;;
        --test_dir) test_dir=$2; shift 2;;
        --geno_test|-geno_test) geno_test=$2; shift 2;;
        --FP|-FP) FP=$2; shift 2;;
        --maf_diff|-maf_diff) maf_diff=$2; shift 2;;
        --out|-out) out_prefix=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done

##########################################################################################
### 1.
### Extract vcf names(training dataset) from vcf file
zcat ${train_dir} | grep 'CHROM' >> ${out_prefix}/CHR${chr_num}_train_names.txt

### 2.
### Store overall result from DPR
mkdir -p ${out_prefix}/DPR_CHR${chr_num}

### Store python log file
mkdir -p ${out_prefix}/DPR_CHR${chr_num}/log_file

### 3.
### DPR input
### Store result for DPR_input
mkdir -p ${out_prefix}/DPR_CHR${chr_num}/DPR_input

### Store cross validation result
mkdir -p ${out_prefix}/DPR_CHR${chr_num}/DPR_input/CV
mkdir -p ${out_prefix}/DPR_CHR${chr_num}/DPR_input/CV/bimbam
mkdir -p ${out_prefix}/DPR_CHR${chr_num}/DPR_input/CV/pheno
mkdir -p ${out_prefix}/DPR_CHR${chr_num}/DPR_input/CV/SNP_annot

### Store bimbam file
mkdir -p ${out_prefix}/DPR_CHR${chr_num}/DPR_input/bimbam

### Store phenotype file
mkdir -p ${out_prefix}/DPR_CHR${chr_num}/DPR_input/pheno

### Store SNP annotation file
mkdir -p ${out_prefix}/DPR_CHR${chr_num}/DPR_input/SNP_annot

### Prepare for DPR input
python ./Model_Train/DPR_Train.py \
--Gene_Exp_path ${Gene_Exp} \
--train_sample ${train_sample} \
--chr_num ${chr_num} \
--thread ${thread} \
--train_dir ${train_dir} \
--train_names ${out_prefix}/CHR${chr_num}_train_names.txt \
--geno ${geno_train} \
--Format ${FT} \
--hwe ${hwe} \
--maf ${maf} \
--window ${window} \
--dpr ${dpr_num} \
--ES ${ES} \
--out_prefix ${out_prefix}/DPR_CHR${chr_num} > \
${out_prefix}/DPR_CHR${chr_num}/log_file/DPR_Train_log.txt

### remove file after finishing input preparation
rm ${out_prefix}/CHR${chr_num}_train_names.txt
rm -r ${out_prefix}/DPR_CHR${chr_num}/DPR_input/CV

### 4.
### If have testing data
### Prediction processing

if [[ "$pred"x == "y"x ]];then
  echo "Prediction Start"
  if [ -z "$test_dir" ];then
    echo "Testing Data is empty"
    echo "Please input Testing Data"
    echo "Prediction Stop"
  else
    zcat ${test_dir} | grep 'CHROM' > ${out_prefix}/CHR${chr_num}_test_names.txt

    python ./Model_Train/Prediction.py \
    --model ${model} \
    --train_result_path ${out_prefix}/DPR_CHR${chr_num}/CHR${chr_num}_DPR_training_param.txt \
    --train_info_path ${out_prefix}/DPR_CHR${chr_num}/CHR${chr_num}_DPR_training_info.txt \
    --chr_num ${chr_num} \
    --test_dir ${test_dir} \
    --test_names ${out_prefix}/CHR${chr_num}_test_names.txt \
    --thread ${thread} \
    --Format ${FP} \
    --geno ${geno_test} \
    --window ${window} \
    --maf_diff ${maf_diff} \
    --out_prefix ${out_prefix}/DPR_CHR${chr_num} > \
    ${out_prefix}/DPR_CHR${chr_num}/log_file/DPR_Prediction_log.txt

    rm ${out_prefix}/CHR${chr_num}_test_names.txt
  fi
elif [[ "$pred"x == "n"x ]];then
    echo "No Testing Data"
else
  echo "Command not found"
fi

















