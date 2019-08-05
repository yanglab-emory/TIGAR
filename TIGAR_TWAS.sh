#!/usr/bin/bash

#############################################################################################################
# software requirement
# python 3
# tabix
###

###
# asso: Method for association study:
#       1) asso==1, providing PED file (phenotype data)
#       2) asso==2, providing previous GWAS Zscore result
# Gene_Exp: Gene annotation and Expression level file path

# For asso==1
# PED : Standard genotype file
# Asso_Info : Instruction for association study
#             1) P : column names for corresponding phenotype in PED file
#             2) C : column names for covariates to regress out
# method : OLS or Logit, default is OLS

# For asso==2
# Zscore : Zscore file from previous GWAS study (tabixed)
# Weight : File contains snps effect size (Same format as prediction output file ?)
# Covar : Reference covariance matrix (Scripts is provided, see in covar_calculation.py, tabixed)
# chr : Chromosome number
# window : Window for selecting genotype data, default 10^6

# thread : Number of thread, default=1
# out : output dir

#############################################################################################################
VARS=`getopt -o "" -a -l \
asso:,Gene_Exp:,PED:,Asso_Info:,method:,block:,geno_path:,Format:,maf:,Zscore:,Weight:,Covar:,chr:,window:,thread:,out: \
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
        --asso|-asso) asso=$2; shift 2;;
        --Gene_Exp|-Gene_Exp) Gene_Exp=$2; shift 2;;
        --PED|-PED) PED=$2; shift 2;;
        --Asso_Info|-Asso_Info) Asso_Info=$2; shift 2;;
        --method|-method) method=$2; shift 2;;
        --Zscore|-Zscore) Zscore=$2; shift 2;;
        --Weight|-Weight) Weight=$2; shift 2;;
        --Covar|-Covar) Covar=$2; shift 2;;
        --chr|-chr) chr_num=$2; shift 2;;
        --window|-window) window=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --out|-out) out_prefix=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done

#### setting default value
thread=${thread:-1}

### asso==1
method=${method:-'OLS'}
PED=${PED:-""}
Asso_Info=${Asso_Info:-""}

### asso==2
window=${window:-$((10**6))}
Zscore=${Zscore:-""}
Weight=${Weight:-""}
Covar=${Covar:-""}
chr_num=${chr_num:-""}

#############################################################################################################
if [[ "$asso"x == "1"x ]];then
    echo "asso is 1"

    mkdir -p ${out_prefix}/TIGAR_TWAS_M1

    ./TWAS/Asso_Study_01.py \
    --Gene_Exp_path ${Gene_Exp} \
    --PED ${PED} \
    --Asso_Info ${Asso_Info} \
    --method ${method} \
    --thread ${thread} \
    --out_prefix ${out_prefix}/TIGAR_TWAS_M1

elif [[ "$asso"x == "2"x ]];then
    echo "asso is 2"

    mkdir -p ${out_prefix}/TIGAR_TWAS_M2

    ./TWAS/Asso_Study_02.sh \
    --Gene_Exp_path ${Gene_Exp} \
    --Zscore ${Zscore} \
    --Weight ${Weight} \
    --Covar ${Covar} \
    --chr ${chr_num} \
    --window ${window} \
    --thread ${thread} \
    --out ${out_prefix}/TIGAR_TWAS_M2
fi


