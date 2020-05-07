#!/usr/bin/bash

#########################################################
####### Required software
# python 3
# tabix
#########################################################

# Variable needed for training
###
# --model: Gene expression prediction model: "elastic_net" or "DPR"
# --Gene_Exp: Path for Gene annotation and Expression file
# --train_sampleID: Path for a file with sampleIDs that will be used for training
# --genofile: Path for the training genotype file (bgzipped and tabixed) 
# --chr: Chromosome number need to be specified with respect to the genotype input data
# --genofile_tye: Genotype file type: "vcf" or "dosage"
# --Format: Genotype format in VCF file that should be used: "GT" (default) for genotype data or "DS" for dosage data, only required if the input genotype file is of VCF file
# --maf: Minor Allele Frequency threshold (ranges from 0 to 1; default 0.01) to exclude rare variants
# --hwe: Hardy Weinberg Equilibrium p-value threshold (default 0.00001) to exclude variants that violated HWE
# --window: Window size around gene transcription starting sites (TSS) for selecting cis-SNPs for fitting gene expression prediction model (default 1000000 for +- 1MB region around TSS)
# --thread: Number of threads for parallel computation (default 1)
# --out_dir: Output directory (will be created if not exist)

### Variables for fitting elastic_net gene expression prediction models only
# --cv: Number of cross validation folds for tuning elastic-net penalty parameter (default 5)
# --alpha: Fixed L1 & L2 penalty ratio for elastic-net model (default 0.5)
#        If alpha=0, equivalent to lasso regression
#        If alpha=1, equivalent to ridge regression

### DPR only
# --dpr: Bayesian inference algorithm used by DPR: "1" (Variational Bayesian) or "2" (MCMC)
# --ES: Output effect size type: "fixed" (default) for fixed effects or "additive" for an addition of fixed and random effects)

######### Setting Default Values ########################
thread=${thread:-1}
maf=${maf:-0.01}
hwe=${hwe:-0.00001}
window=${window:-$((10**6))}
cv=${cv:-5}
alpha=${alpha:-0.5}
dpr_num=${dpr_num:-1} # 1 is for VB ; 2 is for MCMC
ES=${ES:-"fixed"}
thread=${thread:-1}

#################################

VARS=`getopt -o "" -a -l \
model:,Gene_Exp:,train_sampleID:,chr:,genofile_type:,genofile:,Format:,maf:,hwe:,window:,cv:,alpha:,dpr:,ES:,thread:,out_dir: \
-- "$@"`

if [ $? != 0 ]
then
    echo "Please provide input files. Terminating....." >&2
    exit 1
fi
 
eval set -- "$VARS"

while true
do
    case "$1" in
        --model|-model) model=$2; shift 2;;
        --Gene_Exp|-Gene_Exp) Gene_Exp=$2; shift 2;;
        --train_sampleID|-train_sampleID) train_sampleID=$2; shift 2;;
        --chr|-chr) chr=$2; shift 2;;
        --genofile_type|-genofile_type) genofile_type=$2; shift 2;;
        --genofile|-genofile) genofile=$2; shift 2;;
        --Format|-Format) Format=$2; shift 2;;
        --maf|-maf) maf=$2; shift 2;;
        --hwe|-hwe) hwe=$2; shift 2;;
        --window|-window) window=$2; shift 2;;
        --cv|-cv) cv=$2; shift 2;;
        --alpha|-alpha) alpha=$2; shift 2;;
        --dpr|-dpr) dpr_num=$2; shift 2;;
        --ES|-ES) ES=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --out_dir|-out_dir) out_dir=$2; shift 2;;
        --) shift;break;;
        *) echo "Wrong input arguments!";exit 1;;
        esac
done



#### Create output directory if not existed
make -p ${out_dir}

##########################################################################################################
# 1. Model Training & Prediction
if [[ "$model"x == "elastic_net"x ]];then
    echo "Train gene expression imputation models by Elastic-Net method..."

    ./Model_Train_Pred/Elastic_Net.sh \
    --model ${model} \
    --Gene_Exp ${Gene_Exp} \
    --train_sampleID ${train_sampleID} \
    --chr ${chr} \
    --genofile_type ${genofile_type} \
    --genofile ${genofile} \
    --Format ${Format} \
    --maf ${maf} \
    --hwe ${hwe} \
    --window ${window} \
    --cv ${cv} \
    --alpha ${alpha} \
    --thread ${thread} \
    --out_dir ${out_dir}
elif [[ "$model"x == "DPR"x ]]; then
    echo "Train gene expression imputation models by Nonparametric Bayesian DPR method..."

    ./Model_Train_Pred/DPR.sh \
    --model ${model} \
    --Gene_Exp ${Gene_Exp} \
    --train_sampleID ${train_sampleID} \
    --chr ${chr} \
    --genofile_type ${genofile_type} \
    --genofile ${genofile} \
    --Format ${Format} \
    --maf ${maf} \
    --hwe ${hwe} \
    --window ${window} \
    --dpr ${dpr_num} \
    --ES ${ES} \
    --thread ${thread} \
    --out_dir ${out_dir}
else
    echo "Error: Please specify --model to be either elastic_net or DPR "
fi

