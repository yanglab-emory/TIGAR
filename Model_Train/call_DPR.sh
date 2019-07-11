#!/usr/bin/bash

file_dir=$1
dpr_num=$2
TargetID=$3
out_prefix=$4

BIMBAM=${file_dir}/bimbam
PHENO=${file_dir}/pheno
SNP=${file_dir}/SNP_annot

cd ${out_prefix}

echo $TargetID
DPR \
-g ${BIMBAM}/${TargetID}_bimbam.txt \
-p ${PHENO}/${TargetID}_pheno.txt \
-a ${SNP}/${TargetID}_snp_annot.txt \
-dpr ${dpr_num} \
-o DPR_${TargetID}


