#!/usr/bin/bash

#######################################################################
### Input Arguments for GReX Prediction
#######################################################################

# --chr: Chromosome number need to be specified with respect to the genotype input data
# --weight: Path for SNP weight (eQTL effect size) file 
# --gene_anno : Path for gene annotation file to specify a list of gene for GReX prediction
# `--gene_anno`: Gene annotation file to specify the list of genes, which is of the same format as the first five columsn of gene expression file
# --test_sampleID: Path for a file with sampleIDs that should be contained in the genotype file
# --genofile: Path for the training genotype file (bgzipped and tabixed) 
# --genofile_tye: Genotype file type: "vcf" or "dosage"
# --format: Genotype format in VCF file that should be used: "GT" (default) for genotype data or "DS" for dosage data, only required if the input genotype file is of VCF file
# --window: Window size around gene region for selecting cis-SNPs for fitting gene expression prediction model (default 1000000 for +- 1MB region around gene)
# maf_diff: MAF difference threshold for matching SNPs from eQTL weight file and test genotype file. If SNP MAF difference is greater than maf_diff (default 0.2), , the SNP will be excluded
# --TIGAR_dir : Specify the directory of TIGAR source code
# --thread: Number of threads for parallel computation (default 1)
# --out_dir: Output directory (will be created if not exist)


#######################################################################
VARS=`getopt -o "" -a -l \
chr:,weight:,gene_anno:,genofile_type:,genofile:,test_sampleID:,format:,window:,missing_rate:,maf_diff:,TIGAR_dir:,thread:,sub_dir:,out_pred_file:,log_file:,out_dir: \
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
        --chr|-chr) chr=$2; shift 2;;
        --weight|-weight) weight_file=$2; shift 2;;
        --gene_anno|-gene_anno) gene_anno=$2; shift 2;;
        --genofile_type|-genofile_type) genofile_type=$2; shift 2;;
        --genofile|-genofile) genofile=$2; shift 2;;
        --test_sampleID|-test_sampleID) test_sampleID_file=$2; shift 2;;
        --format|-format) format=$2; shift 2;;
        --window|-window) window=$2; shift 2;;
        --missing_rate|-missing_rate) missing_rate=$2; shift 2;;
        --maf_diff|-maf_diff) maf_diff=$2; shift 2;;
        --TIGAR_dir|-TIGAR_dir) TIGAR_dir=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --sub_dir|-sub_dir) sub_dir=$2; shift 2;;
        --out_pred_file|-out_pred_file) out_pred_file=$2; shift 2;;
        --log_file|-log_file) log_file=$2; shift 2;;
        --out_dir|-out_dir) out_dir=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done

#### default value
window=${window:-$((10**6))}
missing_rate=${missing_rate:-0.2}
maf_diff=${maf_diff:-0.2}
thread=${thread:-1}
out_pred_file=${out_pred_file:-CHR${chr}_Pred_GReX}.txt
log_file=${log_file:-CHR${chr}_Pred_log.txt}

# sub_dir: whether to use subdirectory inside out_dir for output files
sub_dir=${sub_dir:-1}

#### Create output directory if not existed
mkdir -p ${out_dir}
mkdir -p ${out_dir}/logs

# sub directory in out directory
if [[ "$sub_dir"x == "1"x ]];then
    out_sub_dir=${out_dir}/Pred_CHR${chr}
else
    out_sub_dir=${out_dir}
fi

mkdir -p ${out_sub_dir}

####################################################
# check tabix command
if [ ! -x "$(command -v tabix)" ]; then
    echo 'Error: required tool TABIX is not available.' >&2
    exit 1
fi

# Check genotype file 
if [ ! -f "${genofile}" ] ; then
    echo Error: Training genotype file ${genofile} does not exist or is empty. >&2
    exit 1
fi

# Check training sample ID file
if [ ! -f "${test_sampleID_file}" ] ; then
    echo Error: Test sample ID file ${test_sampleID_file} does not exist or is empty. >&2
    exit 1
fi

# Check eQTL weight file
if [ ! -f "${weight_file}" ] ; then
    echo Error: eQTL weight file ${weight_file} does not exist or is empty. >&2
    exit 1
fi

# Check gene annotation file
if [ ! -f "${gene_anno}" ] ; then
    echo Error: gene info file ${gene_anno} does not exist or is empty. >&2
    exit 1
fi

# Make python script executible
if [[ ! -x  ${TIGAR_dir}/Model_Train_Pred/Prediction.py ]] ; then
    chmod 755  ${TIGAR_dir}/Model_Train_Pred/Prediction.py
fi

echo Predicting gene expression.

python ${TIGAR_dir}/Model_Train_Pred/Prediction.py \
--chr ${chr} \
--weight ${weight_file} \
--genofile ${genofile} \
--test_sampleID ${test_sampleID_file} \
--format ${format} \
--genofile_type ${genofile_type} \
--gene_anno ${gene_anno} \
--window ${window} \
--thread ${thread} \
--maf_diff ${maf_diff} \
--missing_rate ${missing_rate} \
--out_pred_file ${out_pred_file} \
--TIGAR_dir ${TIGAR_dir} \
--out_dir ${out_sub_dir} \
> ${out_dir}/logs/${log_file}

echo Completed prediction.





