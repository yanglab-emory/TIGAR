#!/usr/bin/bash

#####################################################################################################
# software requirement
# python 3
# tabix
###

# vairable needed for training
###
# model: elastic_net or DPR
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
# out: output dir

### Elastic Net only
# cv: cv-fold cross-validation in model selection, default 5-fold
# alpha: L1 & L2 ratio for elastic net regression, default 0.5
#        If alpha=0, lasso regression
#        If alpha=1, ridge regression

### DPR only
# dpr: model using for DPR
# ES: effect size (fixed/additive). Default is fixed.

# variable needed for prediction
###
# pred: y or n (default is n)
#       y means user wants to conduct prediction
#       n means user do not want to run prediction part
# geno_test: vcf or dosages
# test_dir: testing data path, should be tabix(contains .gz and .tbi)
# FP: Format using for testing data (GT or DS), default DS
# maf_diff: threshold of difference between training maf and testing maf, if difference is larger than this, 
#           then drop this data, default 0.2  

#############################################################################################################
VARS=`getopt -o "" -a -l \
model:,Gene_Exp:,train_sample:,chr:,thread:,geno_train:,train_dir:,FT:,maf:,hwe:,window:,cv:,alpha:,dpr:,ES:,pred:,test_dir:,geno_test:,FP:,maf_diff:,out: \
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
        --cv|-cv) cv=$2; shift 2;;
        --alpha|-alpha) alpha=$2; shift 2;;
        --dpr|-dpr) dpr_num=$2; shift 2;;
        --ES|-ES) ES=$2; shift 2;;
        --pred|-pred) pred=$2; shift 2;;
        --test_dir) test_dir=$2; shift 2;;
        --geno_test) geno_test=$2; shift 2;;
        --FP|-FP) FP=$2; shift 2;;
        --maf_diff|-maf_diff) maf_diff=$2; shift 2;;
        --out|-out) out_prefix=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done

#### setting default value
thread=${thread:-1}
maf=${maf:-0.01}
hwe=${hwe:-0.001}
window=${window:-$((10**6))}
cv=${cv:-5}
alpha=${alpha:-0.5}
dpr_num=${dpr_num:-1}
ES=${ES:-"fixed"}

### default for prediction part
pred=${pred:-n}
test_dir=${test_dir:-""}
maf_diff=${maf_diff:-0.2}

##########################################################################################################
# 1. Model Training & Prediction
if [[ "$model"x == "elastic_net"x ]];then
    echo "Using Elastic Net model for training."

    ./Model_Train/Elastic_Net.sh \
    --model ${model} \
    --Gene_Exp ${Gene_Exp} \
    --train_sample ${train_sample} \
    --chr ${chr_num} \
    --thread ${thread} \
    --geno_train ${geno_train} \
    --train_dir ${train_dir} \
    --FT ${FT} \
    --maf ${maf} \
    --hwe ${hwe} \
    --window ${window} \
    --cv ${cv} \
    --alpha ${alpha} \
    --pred ${pred} \
    --test_dir ${test_dir} \
    --geno_test ${geno_test} \
    --FP ${FP} \
    --maf_diff ${maf_diff} \
    --out ${out_prefix}
elif [[ "$model"x == "DPR"x ]]; then
    echo "Using DPR model for training."

    ./Model_Train/DPR.sh \
    --model ${model} \
    --Gene_Exp ${Gene_Exp} \
    --train_sample ${train_sample} \
    --chr ${chr_num} \
    --thread ${thread} \
    --geno_train ${geno_train} \
    --train_dir ${train_dir} \
    --FT ${FT} \
    --maf ${maf} \
    --hwe ${hwe} \
    --window ${window} \
    --dpr ${dpr_num} \
    --ES ${ES} \
    --pred ${pred} \
    --test_dir ${test_dir} \
    --geno_test ${geno_test} \
    --FP ${FP} \
    --maf_diff ${maf_diff} \
    --out ${out_prefix}

else
    echo "Model not found."
fi















