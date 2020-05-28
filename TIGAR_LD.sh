#!/usr/bin/bash

#######################################################################
### Input Arguments for GReX Prediction
#######################################################################

###
# block : Block annotation
# --genofile: Path for the training genotype file (bgzipped and tabixed) 
# --genofile_tye: Genotype file type: "vcf" or "dosage"
# --format: Genotype format in VCF file that should be used: "GT" (default) for genotype data or "DS" for dosage data, only required if the input genotype file is of VCF file
# --chr: Chromosome number need to be specified with respect to the genotype input data# Format : GT or DS
# --maf : Filter SNPs by Minor Allele Frequency (range from 0-1),default 0.05
# --TIGAR_dir : Specify the directory of TIGAR source code
# --thread: Number of threads for parallel computation (default 1)
# --out_dir: Output directory (will be created if not exist)
###

VARS=`getopt -o "" -a -l \
genom_block:,genofile:,genofile_type:,chr:,format:,maf:,TIGAR_dir:,thread:,out_dir: \
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
        --genom_block|-genom_block) genom_block=$2; shift 2;;
        --genofile|-genofile) genofile=$2; shift 2;;
        --genofile_type|-genofile_type) genofile_type=$2; shift 2;;
        --chr|-chr) chr_num=$2; shift 2;;
        --format|-format) Format=$2; shift 2;;
        --maf|-maf) maf=$2; shift 2;;
        --TIGAR_dir|-TIGAR_dir) TIGAR_dir=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --out_dir|-out_dir) out_dir=$2; shift 2;;
        --) shift;break;;
        *) echo "Please check input arguments!";exit 1;;
        esac
done

#### setting default value
thread=${thread:-1}
maf=${maf:-0}

###############################################################################################
### 1. Run covariance calculation
python ./TWAS/covar_calculation.py \
--block ${block} \
--geno_path ${genofile_path} \
--geno ${genofile_type} \
--chr_num ${chr_num} \
--Format ${Format} \
--maf ${maf} \
--thread ${thread} \
--out ${out_prefix}

### 2. tabix output file
sort -n -k2 ${out_prefix}/CHR${chr_num}_reference_cov.txt | bgzip -c \
> ${out_prefix}/CHR${chr_num}_reference_cov.txt.gz

tabix -p vcf ${out_prefix}/CHR${chr_num}_reference_cov.txt.gz

rm ${out_prefix}/CHR${chr_num}_reference_cov.txt













