#!/usr/bin/bash

#######################################################################
### Input Arguments for DPR model training
#######################################################################

# --model: DPR
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
# --thread: Number of threads for parallel computation (default 1)
# --out_dir: Output directory (will be created if not exist)

# Specific variables for training DPR models
# --dpr: Bayesian inference algorithm used by DPR: "1" (Variational Bayesian) or "2" (MCMC)
# --ES: Output effect size type: "fixed" (default) for fixed effects or "additive" for an addition of fixed and random effects)


VARS=`getopt -o "" -a -l \
model:,Gene_Exp:,train_sampleID:,chr:,genofile_type:,genofile:,Format:,maf:,hwe:,window:,cvR2:,dpr:,ES:,thread:,out_dir: \
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
        --dpr|-dpr) dpr_num=$2; shift 2;;
        --ES|-ES) ES=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --out_dir|-out_dir) out_dir=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done

#################################################
### 1.
### Store DPR Files
mkdir -p ${out_dir}/DPR_CHR${chr}

### 2.
### Store files for DPR under DPR_Files
mkdir -p ${out_dir}/DPR_CHR${chr}/DPR_Files

### 3.
### Extract column names from training genotype file
if [ -f ${genofile} ] ; then
    zcat ${genofile} | grep 'CHROM' > ${out_dir}/DPR_CHR${chr}/geno_colnames.txt
else
    echo Error: ${genofile} dose not exist or empty. >&2
    exit 1
fi

### Store Cross Validation DPR input files and outputs
if [ ${cvR2} == "1" ] ; then
    mkdir -p ${out_dir}/DPR_CHR${chr}/CV_Files
    echo "Run 5-fold cross validation to evaluate DPR model."
else
    echo "Skip 5-fold CV."
fi


### Prepare for DPR input
python ./Model_Train_Pred/DPR_Train.py \
--Gene_Exp ${Gene_Exp} \
--train_sampleID ${train_sampleID} \
--chr ${chr} \
--thread ${thread} \
--train_geno_file ${genofile} \
--geno_colnames ${out_dir}/DPR_CHR${chr}/geno_colnames.txt \
--geno ${genofile_type} \
--Format ${Format} \
--hwe ${hwe} \
--maf ${maf} \
--window ${window} \
--cvR2 ${cvR2} \
--dpr ${dpr_num} \
--ES ${ES} \
--out_dir ${out_dir}/DPR_CHR${chr} \
> ${out_dir}/DPR_CHR${chr}/CHR${chr}_DPR_train_Log.txt

### 4. Remove DPR input files
# rm -f ${out_dir}/CHR${chr}_geno_colnames.txt
# rm -fr ${out_dir}/DPR_CHR${chr}/DPR_Files


















