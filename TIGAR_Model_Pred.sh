#!/usr/bin/bash

##################################################################
# software requirement
# python 3
# tabix
###

#######################################################################
### Input Arguments for GReX Prediction
#######################################################################

# --chr: Chromosome number need to be specified with respect to the genotype input data
# --weight: Path for SNP weight (eQTL effect size) file 
# --test_sampleID: Path for a file with sampleIDs that should be contained in the genotype file
# --genofile: Path for the training genotype file (bgzipped and tabixed) 
# --genofile_tye: Genotype file type: "vcf" or "dosage"
# --Format: Genotype format in VCF file that should be used: "GT" (default) for genotype data or "DS" for dosage data, only required if the input genotype file is of VCF file
# --window: Window size around gene transcription starting sites (TSS) for selecting cis-SNPs for fitting gene expression prediction model (default 1000000 for +- 1MB region around TSS)
# maf_diff: MAF difference threshold for matching SNPs from eQTL weight file and test genotype file. If SNP MAF difference is greater than maf_diff (default 0.2), , the SNP will be excluded
# --thread: Number of threads for parallel computation (default 1)
# --out_dir: Output directory (will be created if not exist)


#######################################################################
VARS=`getopt -o "" -a -l \
chr:,weight:,gene_anno:,genofile_type:,genofile:,test_sampleID:,Format:,window:,maf_diff:,thread:,out_dir: \
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
        --Format|-Format) Format=$2; shift 2;;
        --window|-window) window=$2; shift 2;;
        --maf_diff|-maf_diff) maf_diff=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --out_dir|-out_dir) out_dir=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done

#### default value
window=${window:-$((10**6))}
maf_diff=${maf_diff:-0.2}
thread=${thread:-1}

#### Create output directory if not existed
mkdir -p ${out_dir}
mkdir -p ${out_dir}/Pred_CHR${chr}

####################################################
# check tabix command
if [ ! -x "$(command -v tabix)" ]; then
    echo 'Error: required tool TABIX is not available.' >&2
    exit 1
fi

# Check genotype file 
if [ ! -f "${genofile}" ] ; then
    echo Error: Training genotype file ${genofile} dose not exist or empty. >&2
    exit 1
else
    zcat ${genofile} | grep 'CHROM' > ${out_dir}/Pred_CHR${chr}/test_geno_colnames.txt
fi

# Check training sample ID file
if [ ! -f "${test_sampleID_file}" ] ; then
    echo Error: Test sample ID file ${test_sampleID_file} dose not exist or empty. >&2
    exit 1
fi

# Check eQTL weight file
if [ ! -f "${weight_file}" ] ; then
    echo Error: eQTL weight file ${weight_file} dose not exist or empty. >&2
    exit 1
fi

# Check gene annotation file
if [ ! -f "${gene_anno}" ] ; then
    echo Error: eQTL weight file ${gene_anno} dose not exist or empty. >&2
    exit 1
fi

python ./Model_Train_Pred/Prediction.py \
--chr ${chr} \
--weight ${weight_file} \
--genofile ${genofile} \
--test_geno_colnames ${out_dir}/Pred_CHR${chr}/test_geno_colnames.txt \
--test_sampleID ${test_sampleID_file} \
--Format ${Format} \
--geno ${genofile_type} \
--gene_anno ${gene_anno} \
--window ${window} \
--thread ${thread} \
--maf_diff ${maf_diff} \
--out_dir ${out_dir}/Pred_CHR${chr} \
> ${out_dir}/Pred_CHR${chr}/CHR${chr}_Pred_Log.txt
    
# rm -f ${out_dir}/Pred_CHR${chr}/test_geno_colnames.txt







