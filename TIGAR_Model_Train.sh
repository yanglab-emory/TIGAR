#!/usr/bin/bash

#########################################################
#### Required software
#       -python 3
#       -tabix
#       -DPR
#########################################################

# Variable needed for training
###
# --model: Gene expression prediction model: "elastic_net" or "DPR"
# --gene_exp: Path for Gene annotation and Expression file
# --train_sampleID: Path for a file with sampleIDs that will be used for training
# --genofile: Path for the training genotype file (bgzipped and tabixed) 
# --chr: Chromosome number need to be specified with respect to the genotype input data
# --genofile_type: Genotype file type: "vcf" or "dosage"
# --format: Genotype format in VCF file that should be used: "GT" (default) for genotype data or "DS" for dosage data, only required if the input genotype file is of VCF file
# --maf: Minor Allele Frequency threshold (ranges from 0 to 1; default 0.01) to exclude rare variants
# --hwe: Hardy Weinberg Equilibrium p-value threshold (default 0.00001) to exclude variants that violated HWE
# --window: Window size around gene region for selecting cis-SNPs for fitting gene expression prediction model (default 1000000 for +- 1MB region around gene)
# --cvR2: Take value 0 for calculating training R2 from fitted model and 1 for calculating training R2 from 5-fold cross validation
# --TIGAR_dir : Specify the directory of TIGAR source code
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

#################################
VARS=`getopt -o "" -a -l \
model:,gene_exp:,train_sampleID:,chr:,genofile_type:,genofile:,format:,maf:,hwe:,window:,cvR2:,cv:,alpha:,use_alpha:,dpr:,ES:,TIGAR_dir:,thread:,out_dir: \
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
        --gene_exp|-gene_exp) gene_exp=$2; shift 2;;
        --train_sampleID|-train_sampleID) train_sampleID=$2; shift 2;;
        --chr|-chr) chr=$2; shift 2;;
        --genofile_type|-genofile_type) genofile_type=$2; shift 2;;
        --genofile|-genofile) genofile=$2; shift 2;;
        --format|-format) format=$2; shift 2;;
        --maf|-maf) maf=$2; shift 2;;
        --hwe|-hwe) hwe=$2; shift 2;;
        --window|-window) window=$2; shift 2;;
        --cvR2|-cvR2) cvR2=$2; shift 2;;
        --cv|-cv) cv=$2; shift 2;;
        --alpha|-alpha) alpha=$2; shift 2;;
        --use_alpha|-use_alpha) use_alpha=$2; shift 2;;
        --dpr|-dpr) dpr_num=$2; shift 2;;
        --ES|-ES) ES=$2; shift 2;;
        --TIGAR_dir|-TIGAR_dir) TIGAR_dir=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --out_dir|-out_dir) out_dir=$2; shift 2;;
        --) shift;break;;
        *) echo "Wrong input arguments!";exit 1;;
        esac
done

##########################################
# Setting Default Input Argument Values 
##########################################

thread=${thread:-1}
maf=${maf:-0.01}
hwe=${hwe:-0.00001}
window=${window:-$((10**6))}
cvR2=${cvR2:-1}
cv=${cv:-5}
alpha=${alpha:-0.5}
dpr_num=${dpr_num:-1} # 1 is for VB ; 2 is for MCMC
ES=${ES:-"fixed"}
use_alpha=${use_alpha:-1}

#### Create output directory if not existed
mkdir -p ${out_dir}
mkdir -p ${out_dir}/logs

# check tabix command
if [ ! -x "$(command -v tabix)" ]; then
    echo 'Error: required tool TABIX is not available.' >&2
    exit 1
fi

# Check gene expression file
if [ ! -f "${gene_exp}" ] ; then
    echo Error: Gene expression file ${gene_exp} does not exist or is empty. >&2
    exit 1
fi

# Check training sample ID file
if [ ! -f "${train_sampleID}" ] ; then
    echo Error: Training sample ID file ${train_sampleID} does not exist or is empty. >&2
    exit 1
fi

# Check genotype file 
if [ ! -f "${genofile}" ] ; then
    echo Error: Training genotype file ${genofile} does not exist or is empty. >&2
    exit 1
fi


####################################
# 1. Model Training 
####################################

if [[ "$model"x == "elastic_net"x ]];then
    echo "Training gene expression imputation models using Elastic-Net method..."

    # Make python script executible
    if [[ ! -x ${TIGAR_dir}/Model_Train_Pred/Elastic_Net_Train.py ]] ; then
        chmod 755 ${TIGAR_dir}/Model_Train_Pred/Elastic_Net_Train.py
    fi

    mkdir -p ${out_dir}/EN_CHR${chr}

    python ${TIGAR_dir}/Model_Train_Pred/Elastic_Net_Train.py \
    --gene_exp ${gene_exp} \
    --train_sampleID ${train_sampleID} \
    --chr ${chr} \
    --genofile ${genofile} \
    --genofile_type ${genofile_type} \
    --format ${format} \
    --maf ${maf} \
    --hwe ${hwe} \
    --window ${window} \
    --cvR2 ${cvR2} \
    --cv ${cv} \
    --thread ${thread} \
    --alpha ${alpha} \
    --use_alpha ${use_alpha} \
    --TIGAR_dir ${TIGAR_dir} \
    --out_dir ${out_dir}/EN_CHR${chr} \
    > ${out_dir}/logs/CHR${chr}_EN_train_log.txt

    # set temp, weight filepaths for sorting
    temp=${out_dir}/EN_CHR${chr}/temp_CHR${chr}_EN_train_eQTLweights.txt
    weight=${out_dir}/EN_CHR${chr}/CHR${chr}_EN_train_eQTLweights.txt


elif [[ "$model"x == "DPR"x ]]; then
    echo "Training gene expression imputation models using Nonparametric Bayesian DPR method..."

    ### Store DPR Results
    mkdir -p ${out_dir}/DPR_CHR${chr}

    ### Store files for DPR under DPR_Files
    mkdir -p ${out_dir}/DPR_CHR${chr}/DPR_Files

    ### Store Cross Validation DPR input files and outputs
    if [ ${cvR2} == "1" ] ; then
        mkdir -p ${out_dir}/DPR_CHR${chr}/CV_Files
        echo "Running 5-fold cross validation to evaluate DPR model."
    else
        echo "Skipping 5-fold CV."
    fi

    # Make DPR file executible
    if [[ ! -x ${TIGAR_dir}/Model_Train_Pred/DPR ]] ; then
        chmod 755 ${TIGAR_dir}/Model_Train_Pred/DPR
    fi

    # Make python script executible
    if [[ ! -x ${TIGAR_dir}/Model_Train_Pred/DPR_Train.py ]] ; then
        chmod 755 ${TIGAR_dir}/Model_Train_Pred/DPR_Train.py
    fi

    python ${TIGAR_dir}/Model_Train_Pred/DPR_Train.py \
    --gene_exp ${gene_exp} \
    --train_sampleID ${train_sampleID} \
    --chr ${chr} \
    --genofile ${genofile} \
    --genofile_type ${genofile_type} \
    --format ${format} \
    --hwe ${hwe} \
    --maf ${maf} \
    --window ${window} \
    --cvR2 ${cvR2} \
    --dpr ${dpr_num} \
    --ES ${ES} \
    --TIGAR_dir ${TIGAR_dir} \
    --thread ${thread} \
    --out_dir ${out_dir}/DPR_CHR${chr} \
    > ${out_dir}/logs/CHR${chr}_DPR_train_log.txt

    ### 4. Remove DPR input files
    echo Removing DPR input files used for training.
    rm -fr ${out_dir}/DPR_CHR${chr}/DPR_Files

    if [ ${cvR2} == "1" ] ; then
        echo Removing DPR input files used for CV.
        rm -fr ${out_dir}/DPR_CHR${chr}/CV_Files
    fi

    # set temp, weight filepaths for sorting
    temp=${out_dir}/DPR_CHR${chr}/temp_CHR${chr}_DPR_train_eQTLweights.txt
    weight=${out_dir}/DPR_CHR${chr}/CHR${chr}_DPR_train_eQTLweights.txt


else
    echo "Error: Please specify --model as either elastic_net or DPR "
    exit 1
fi

echo "Training ${model} model job completed for CHR ${chr}."

# SORT, BGZIP, AND TABIX
echo "Sort/bgzip/tabix-ing output weight file."
head -n1 ${temp} > ${weight}

tail -n+2 ${temp} | \
sort -nk1 -nk2 >> ${weight} && \
rm ${temp}

if [ ! -f "${temp}" ] ; then
    echo "Sort successful. Bgzip/tabix-ing."
    bgzip -f ${weight} && \
    tabix -f -b 2 -e 2 -S 1 ${weight}.gz

else
    echo "Sort failed; Unable to bgzip/tabix output weights file."
    exit 1
fi



exit
