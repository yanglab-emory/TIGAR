#!/usr/bin/bash

#############################################################################################
### vairable needed for training
# model: DPR
# Gene_Exp: Gene annotation and Expression level file
# train_sample: a column of sampleIDs use for training
# chr: chromosome number for corresponding training data
# thread: number of thread for multiprocessing
# genofile_type: vcf or dosages, both should be tabix
# genofile_dir: training data path
# Format: Format using for training data(GT or DS), default DS
# maf: Threshold for Minor Allele Frequency (range from 0-1),default 0.01
# hwe: Threshold of p-value for Hardy Weinberg Equilibrium exact test,default 0.001
# window: window for selecting data
# dpr: model when we choose DPR, default 1
# ES: effect size(fixed/additive)
# out: output dir

VARS=`getopt -o "" -a -l \
model:,Gene_Exp:,train_sample:,chr:,genofile_type:,genofile_dir:,Format:,maf:,hwe:,window:,dpr:,ES:,thread:,out: \
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
        --train_sample|-train_sample) train_sample=$2; shift 2;;
        --chr|-chr) chr_num=$2; shift 2;;
        --genofile_type|-genofile_type) genofile_type=$2; shift 2;;
        --genofile_dir|-genofile_dir) genofile_dir=$2; shift 2;;
        --Format|-Format) Format=$2; shift 2;;
        --maf|-maf) maf=$2; shift 2;;
        --hwe|-hwe) hwe=$2; shift 2;;
        --window|-window) window=$2; shift 2;;
        --dpr|-dpr) dpr_num=$2; shift 2;;
        --ES|-ES) ES=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --out|-out) out_prefix=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done

##########################################################################################
### 1.
### Extract vcf names(training dataset) from vcf file
zcat ${genofile_dir} | grep 'CHROM' > ${out_prefix}/CHR${chr_num}_train_names.txt

### 2.
### Store overall result from DPR
mkdir -p ${out_prefix}/DPR_CHR${chr_num}

### Store python log file
mkdir -p ${out_prefix}/DPR_CHR${chr_num}/log_file

### 3.
### DPR input
### Store result for DPR_input
mkdir -p ${out_prefix}/DPR_CHR${chr_num}/DPR_input

### Store cross validation result
mkdir -p ${out_prefix}/DPR_CHR${chr_num}/DPR_input/CV
mkdir -p ${out_prefix}/DPR_CHR${chr_num}/DPR_input/CV/bimbam
mkdir -p ${out_prefix}/DPR_CHR${chr_num}/DPR_input/CV/pheno
mkdir -p ${out_prefix}/DPR_CHR${chr_num}/DPR_input/CV/SNP_annot

### Store bimbam file
mkdir -p ${out_prefix}/DPR_CHR${chr_num}/DPR_input/bimbam

### Store phenotype file
mkdir -p ${out_prefix}/DPR_CHR${chr_num}/DPR_input/pheno

### Store SNP annotation file
mkdir -p ${out_prefix}/DPR_CHR${chr_num}/DPR_input/SNP_annot

### Prepare for DPR input
python ./Model_Train_Pred/DPR_Train.py \
--Gene_Exp_path ${Gene_Exp} \
--train_sample ${train_sample} \
--chr_num ${chr_num} \
--thread ${thread} \
--train_dir ${genofile_dir} \
--train_names ${out_prefix}/CHR${chr_num}_train_names.txt \
--geno ${genofile_type} \
--Format ${Format} \
--hwe ${hwe} \
--maf ${maf} \
--window ${window} \
--dpr ${dpr_num} \
--ES ${ES} \
--out_prefix ${out_prefix}/DPR_CHR${chr_num} \
> ${out_prefix}/DPR_CHR${chr_num}/log_file/DPR_Train_log.txt

### 4. remove file after finishing input preparation
rm ${out_prefix}/CHR${chr_num}_train_names.txt
rm -r ${out_prefix}/DPR_CHR${chr_num}/DPR_input/CV


















