#!/usr/bin/bash

file_dir=$1
dpr_num=$2
TargetID=$3
TIGAR_dir=$4

## Shell scripts to call C++ executible file DPR
# DPR=/mnt/icebreaker/data2/home/snagpal/DPR/DPR-master/DPR

if [ ! -x "$(command -v DPR)" ]; then
	echo 'Error: please add DPR executible file into your PATH.' >&2
	exit 1
else
	echo Running DPR with $TargetID ...
	cd ${file_dir}
	${TIGAR_dir}/Model_Train_Pred/DPR \
	-g ${TargetID}_bimbam.txt \
	-p ${TargetID}_pheno.txt \
	-a ${TargetID}_snp_annot.txt \
	-dpr ${dpr_num} \
	-o DPR_${TargetID}
fi


