#!/usr/bin/bash

#####################################################################################################
# software requirement
# python 3
# tabix
###

###
# block : Block annotation
# geno_path : genotype information (tabixed)
# geno : vcf or dosages
# chr : Chromosome number
# Format : GT or DS
# maf : Threshold for Minor Allele Frequency (range from 0-1),default 0.05
# thread : Number of thread
# out : output dir
###

VARS=`getopt -o "" -a -l \
block:,geno_path:,geno:,chr:,Format:,maf:,thread:,out: \
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
        --block|-block) block=$2; shift 2;;
        --geno_path|-geno_path) geno_path=$2; shift 2;;
        --geno|-geno) geno=$2; shift 2;;
        --chr|-chr) chr_num=$2; shift 2;;
        --Format|-Format) Format=$2; shift 2;;
        --maf|-maf) maf=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --out|-out) out_prefix=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done

#### setting default value
thread=${thread:-1}
maf=${maf:-0.05}

###############################################################################################
### 1. Run covariance calculation
python ./TWAS/covar_calculation.py \
--block ${block} \
--geno_path ${geno_path} \
--geno ${geno} \
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













