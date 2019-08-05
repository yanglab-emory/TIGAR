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
# genofile_tye: vcf or dosages
# genofile_dir: training data path, should be tabix(contains .gz and .tbi)
# Format: Format using for training data(GT or DS), default DS
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

#############################################################################################################
VARS=`getopt -o "" -a -l \
model:,Gene_Exp:,train_sample:,chr:,genofile_type:,genofile_dir:,Format:,maf:,hwe:,window:,cv:,alpha:,dpr:,ES:,thread:,out: \
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
        --dpr|-dpr) dpr_num=$2; shift 2;;
        --ES|-ES) ES=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
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
thread=${thread:-1}

##########################################################################################################
# 1. Model Training & Prediction
if [[ "$model"x == "elastic_net"x ]];then
    echo "Using Elastic Net model for training."

    ./Model_Train_Pred/Elastic_Net.sh \
    --model ${model} \
    --Gene_Exp ${Gene_Exp} \
    --train_sample ${train_sample} \
    --chr ${chr_num} \
    --genofile_type ${genofile_type} \
    --genofile_dir ${genofile_dir} \
    --Format ${Format} \
    --maf ${maf} \
    --hwe ${hwe} \
    --window ${window} \
    --cv ${cv} \
    --alpha ${alpha} \
    --thread ${thread} \
    --out ${out_prefix}
elif [[ "$model"x == "DPR"x ]]; then
    echo "Using DPR model for training."

    ./Model_Train_Pred/DPR.sh \
    --model ${model} \
    --Gene_Exp ${Gene_Exp} \
    --train_sample ${train_sample} \
    --chr ${chr_num} \
    --genofile_type ${genofile_type} \
    --genofile_dir ${genofile_dir} \
    --Format ${Format} \
    --maf ${maf} \
    --hwe ${hwe} \
    --window ${window} \
    --dpr ${dpr_num} \
    --ES ${ES} \
    --thread ${thread} \
    --out ${out_prefix}
else
    echo "Model not found."
fi