#!/usr/bin/bash

#######################################################################
### Input Arguments for VC_TWAS with summary stat
#######################################################################

# --chr: Chromosome number need to be specified with respect to the genotype input data
# --weight: Path for SNP weight (eQTL effect size) file 
# --gene_anno: Path for gene annotation file to specify a list of gene for GReX prediction
# --GWAS_result: Path to GWAS result file 
# --sample_size: Sample size of summary-level GWAS data
# --LD: Path to Reference LD genotype covariance file
# --window: Window size around gene transcription starting sites (TSS) for selecting cis-SNPs for fitting gene expression prediction model (default 1000000 for +- 1MB region around TSS)
# --weight_threshold: Threshold for estimated cis-eQTL effect sizes, filter SNPs with absolute cis-eQTL effect sizes smaller than threshold
# --TIGAR_dir: Specify the directory of TIGAR source code
# --thread: Number of threads for parallel computation (default 1)
# --out_dir: Output directory (will be created if not exist)


#######################################################################
VARS=`getopt -o "" -a -l \
gene_anno:,GWAS_result:,Weight:,sample_size:,weight_threshold:,LD:,chr:,window:,TIGAR_dir:,thread:,out_dir: \
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
        --gene_anno|-gene_anno) gene_anno=$2; shift 2;;
        --GWAS_result|-GWAS_result) GWAS_result=$2; shift 2;;
        --Weight|-Weight) Weight=$2; shift 2;;
        --sample_size|-sample_size) sample_size=$2; shift 2;;
        --weight_threshold|-weight_threshold) weight_threshold=$2; shift 2;;
        --LD|-LD) LD=$2; shift 2;;
        --chr|-chr) chr=$2; shift 2;;
        --window|-window) window=$2; shift 2;;
        --TIGAR_dir|-TIGAR_dir) TIGAR_dir=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --out_dir|-out_dir) out_dir=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done

#### default value
window=${window:-$((10**6))}
thread=${thread:-1}

#### Create output directory if not existed
mkdir -p ${out_dir}
mkdir -p ${out_dir}/VC_TWAS_summary_CHR${chr}

####################################################
# check tabix command
if [ ! -x "$(command -v tabix)" ]; then
    echo 'Error: required tool TABIX is not available.' >&2
    exit 1
fi

# Check Weight file 
if [ ! -f "${Weight}" ] ; then
    echo Error: Training genotype file ${Weight} dose not exist or empty. >&2
    exit 1
fi


# Check gwas result file
if [ ! -f "${GWAS_result}" ] ; then
    echo Error: eQTL weight file ${GWAS_result} dose not exist or empty. >&2
    exit 1
fi

# Check gene annotation file
if [ ! -f "${gene_anno}" ] ; then
    echo Error: eQTL weight file ${gene_anno} dose not exist or empty. >&2
    exit 1
fi

# Make python script executible
if [[ ! -x  ${TIGAR_dir}/VC_TWAS/VC_TWAS_summary.py ]] ; then
    #chmod 755  ${TIGAR_dir}/VC_TWAS/VC_TWAS.py
    chmod 755  ${TIGAR_dir}/VC_TWAS/VC_TWAS_summary.py
fi

#python ${TIGAR_dir}/VC_TWAS/VC_TWAS.py \
python ${TIGAR_dir}/VC_TWAS/VC_TWAS_summary.py \
--gene_anno ${gene_anno} \
--GWAS_result ${GWAS_result} \
--Weight ${Weight} \
--sample_size ${sample_size} \
--weight_threshold ${weight_threshold} \
--LD ${LD} \
--chr ${chr} \
--window ${window} \
--TIGAR_dir ${TIGAR_dir} \
--thread ${thread} \
--out_dir ${out_dir}/VC_TWAS_summary_CHR${chr} \
> ${out_dir}/logs/CHR${chr}_sum_VC_TWAS_log.txt

# > ${out_dir}/VC_TWAS_summary_CHR${chr}/CHR${chr}_sum_VC_TWAS_Log.txt
    
# rm -f ${out_dir}/Pred_CHR${chr}/test_geno_colnames.txt







