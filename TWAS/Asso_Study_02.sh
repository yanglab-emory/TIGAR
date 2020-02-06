#!/usr/bin/bash

#####################################################################################################
# software requirement
# python 3
# tabix
###

# Gene_anno : Gene annotation file directory
# Zscore : Zscore file from previous GWAS study (tabixed)
# Weight : File contains snps effect size (Same format as prediction output file)
# Covar : Reference covariance matrix (Scripts is provided, see in covar_calculation.py, tabixed)
# chr : Chromosome number
# window : Window for selecting genotype data, default 10^6
# thread : Number of thread, default 1
# out : output directory

#####################################################################################################
VARS=`getopt -o "" -a -l \
Gene_anno:,Zscore:,Weight:,Covar:,chr:,window:,thread:,out: \
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
        --Gene_anno|-Gene_anno) Gene_anno=$2; shift 2;;
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

#####################################################################################################
### tabix training parameter file
sort -n -k2 ${Weight} | bgzip -c > ${out_prefix}/CHR${chr_num}_weight.txt.gz
tabix -p vcf ${out_prefix}/CHR${chr_num}_weight.txt.gz

### extract gene annotation file
cut -f1-5 ${Gene_anno} > ${out_prefix}/CHR${chr_num}_Gene_Annotation.txt

#zcat ${out_prefix}/${Weight}.gz | grep 'CHROM' > ${out_prefix}/CHR${chr_num}_weight_names.txt

grep 'CHROM' ${Weight} > ${out_prefix}/CHR${chr_num}_weight_names.txt

zcat ${Zscore} | grep 'CHROM' > ${out_prefix}/CHR${chr_num}_GWAS_names.txt

python ./TWAS/summary_stat.py \
--Gene ${out_prefix}/CHR${chr_num}_Gene_Annotation.txt \
--Zscore ${Zscore} \
--Zscore_names ${out_prefix}/CHR${chr_num}_GWAS_names.txt \
--Weight ${out_prefix}/CHR${chr_num}_weight.txt.gz \
--Weight_names ${out_prefix}/CHR${chr_num}_weight_names.txt \
--Covar ${Covar} \
--chr_num ${chr_num} \
--window ${window} \
--thread ${thread} \
--out_prefix ${out_prefix}

rm ${out_prefix}/CHR${chr_num}_weight.txt.gz
rm ${out_prefix}/CHR${chr_num}_weight.txt.gz.tbi
rm ${out_prefix}/CHR${chr_num}_Gene_Annotation.txt
rm ${out_prefix}/CHR${chr_num}_weight_names.txt
rm ${out_prefix}/CHR${chr_num}_GWAS_names.txt

