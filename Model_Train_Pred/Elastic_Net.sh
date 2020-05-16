#!/usr/bin/bash

#########################################################################
# vairable needed for training
###
# --model: elastic_net
# --Gene_Exp: Path for Gene annotation and Expression file
# --train_sampleID: Path for a file with sampleIDs that will be used for training
# --genofile: Path for the training genotype file (bgzipped and tabixed) 
# --chr: Chromosome number need to be specified with respect to the genotype input data
# --genofile_tye: Genotype file type: "vcf" or "dosage"
# --Format: Genotype format in VCF file that should be used: "GT" (default) for genotype data or "DS" for dosage data, only required if the input genotype file is of VCF file
# --maf: Minor Allele Frequency threshold (ranges from 0 to 1; default 0.01) to exclude rare variants
# --hwe: Hardy Weinberg Equilibrium p-value threshold (default 0.00001) to exclude variants that violated HWE
# --window: Window size around gene transcription starting sites (TSS) for selecting cis-SNPs for fitting gene expression prediction model (default 1000000 for +- 1MB region around TSS)
# --cvR2: Take value 0 for calculating training R2 from fitted model and 1 for calculating training R2 from 5-fold cross validation
# --TIGAR_dir : Specify the directory of TIGAR source code
# --thread: Number of threads for parallel computation (default 1)
# --out_dir: Output directory (will be created if not exist)

## Specific variables for taining elastic_net models
# --cv: Number of cross validation folds for tuning elastic-net penalty parameter (default 5)
# --alpha: Fixed L1 & L2 penalty ratio for elastic-net model (default 0.5)
#        If alpha=0, equivalent to lasso regression
#        If alpha=1, equivalent to ridge regression


VARS=`getopt -o "" -a -l \
model:,Gene_Exp:,train_sampleID:,chr:,genofile_type:,genofile:,Format:,maf:,hwe:,window:,cvR2:,cv:,alpha:,TIGAR_dir:,thread:,out_dir: \
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
        --train_sampleID|-train_sampleID) train_sampleID=$2; shift 2;;
        --chr|-chr) chr=$2; shift 2;;
        --genofile_type|-genofile_type) genofile_type=$2; shift 2;;
        --genofile|-genofile) genofile=$2; shift 2;;
        --Format|-Format) Format=$2; shift 2;;
        --maf|-maf) maf=$2; shift 2;;
        --hwe|-hwe) hwe=$2; shift 2;;
        --window|-window) window=$2; shift 2;;
        --cvR2|-cvR2) cvR2=$2; shift 2;;
        --cv|-cv) cv=$2; shift 2;;
        --alpha|-alpha) alpha=$2; shift 2;;
        --TIGAR_dir|-TIGAR_dir) TIGAR_dir=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --out_dir|-out_dir) out_dir=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done

###################################################
### 1. 
### Create dir & Store Result
mkdir -p ${out_dir}/EN_CHR${chr}

### 2. 
### Extract genotype file header
zcat ${genofile} | grep 'CHROM' > ${out_dir}/EN_CHR${chr}/geno_colnames.txt
echo ${out_dir}/EN_CHR${chr}/geno_colnames.txt

### 3.


### Training

if [[ ! -x ${TIGAR_dir}/Model_Train_Pred/Elastic_Net_Train.py ]] ; then
    chmod 755 ${TIGAR_dir}/Model_Train_Pred/Elastic_Net_Train.py
fi

python ${TIGAR_dir}/Model_Train_Pred/Elastic_Net_Train.py \
--Gene_Exp ${Gene_Exp} \
--train_sampleID ${train_sampleID} \
--chr ${chr} \
--thread ${thread} \
--genofile ${genofile} \
--geno_colnames ${out_dir}/EN_CHR${chr}/geno_colnames.txt \
--geno ${genofile_type} \
--Format ${Format} \
--maf ${maf} \
--hwe ${hwe} \
--window ${window} \
--cvR2 ${cvR2} \
--cv ${cv} \
--alpha ${alpha} \
--out_dir ${out_dir}/EN_CHR${chr} \
> ${out_dir}/EN_CHR${chr}/CHR${chr}_EN_train_Log.txt

### 4.
### Remove file
# rm -f ${out_dir}/EN_CHR${chr}/geno_colnames.txt









