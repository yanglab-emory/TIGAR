#!/usr/bin/bash

#####################################################################################################
# software requirement
# python 3
# tabix
###

### variable needed for prediction
# model: elastic_net or DPR
# chr: chromosome number for corresponding training data
# train_result_path: training parameter file (same format as training output)
# train_info_path: training information file (same format as training output)
# genofile_type: vcf or dosages
# genofile_dir: test genotype file
# Format: Format using for testing data (GT or DS)
# window: window for selecting data, default 10**6
# maf_diff: threshold of difference between training maf and testing maf, if difference is larger than this value, 
#           then drop data, default 0.2 

# thread: number of thread for multiprocessing

#####################################################################################################
VARS=`getopt -o "" -a -l \
model:,chr:,train_result_path:,train_info_path:,genofile_type:,genofile_dir:,test_sample:,Format:,window:,maf_diff:,thread:,out: \
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
        --chr|-chr) chr_num=$2; shift 2;;
        --train_result_path) train_result_path=$2; shift 2;;
        --train_info_path|-train_info_path) train_info_path=$2; shift 2;;
        --genofile_type|-genofile_type) genofile_type=$2; shift 2;;
        --genofile_dir|-genofile_dir) genofile_dir=$2; shift 2;;
        --test_sample|-test_sample) test_sample=$2; shift 2;;
        --Format|-Format) Format=$2; shift 2;;
        --window|-window) window=$2; shift 2;;
        --maf_diff|-maf_diff) maf_diff=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --out|-out) out_prefix=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done

#### default value
window=${window:-$((10**6))}
maf_diff=${maf_diff:-0.2}
thread=${thread:-1}

##############################################################################

zcat ${genofile_dir} | grep 'CHROM' > ${out_prefix}/CHR${chr_num}_test_names.txt

mkdir -p ${out_prefix}/${model}_CHR${chr_num}
mkdir -p ${out_prefix}/${model}_CHR${chr_num}/log_file

python ./Model_Train_Pred/Prediction.py \
--model ${model} \
--chr_num ${chr_num} \
--train_result_path ${train_result_path} \
--train_info_path ${train_info_path} \
--test_dir ${genofile_dir} \
--test_names ${out_prefix}/CHR${chr_num}_test_names.txt \
--test_sample ${test_sample} \
--Format ${Format} \
--geno ${genofile_type} \
--window ${window} \
--thread ${thread} \
--maf_diff ${maf_diff} \
--out_prefix ${out_prefix}/${model}_CHR${chr_num} \
> ${out_prefix}/${model}_CHR${chr_num}/log_file/Elastic_Net_Prediction_log.txt
    
rm ${out_prefix}/CHR${chr_num}_test_names.txt







