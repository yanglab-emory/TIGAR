#!/usr/bin/bash

########################################################################
### Input Arguments for TWAS with summary-level GWAS Z-score statistics
########################################################################

# --gene_anno : Path for gene annotation file to specify a list of gene for GReX prediction
# --chr: Chromosome number need to be specified 
# --weight: Path for SNP weight (eQTL effect size) file 
# --Zscore : Path for GWAS summary Zscore statistics
# --LD : Path for reference LD (SNP genotype covariance matrix) that should be bgzipped and tabixed. Can be generated using our `covar_calculation.py` script
# --window: Window size around gene transcription starting sites (TSS) for selecting cis-SNPs for fitting gene expression prediction model (default 1000000 for +- 1MB region around TSS)
# --TIGAR_dir : Specify the directory of TIGAR source code
# --thread: Number of threads for parallel computation (default 1)
# --out_dir: Output directory (will be created if not exist)


#############################################################
VARS=`getopt -o "" -a -l \
gene_anno:,Zscore:,weight:,LD:,chr:,window:,TIGAR_dir:,thread:,out_dir: \
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
        --Zscore|-Zscore) Zscore=$2; shift 2;;
        --weight|-weight) weight=$2; shift 2;;
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

######## Extract colum names of weight and Zscore files
zcat ${weight} | head -n1 > ${out_dir}/CHR${chr}_weight_colnames.txt
zcat ${Zscore} | head -n1 > ${out_dir}/CHR${chr}_Zscore_colnames.txt

if [[ ! -x  ${TIGAR_dir}/TWAS/Asso_Study_02.py ]] ; then
    chmod 755  ${TIGAR_dir}/TWAS/Asso_Study_02.py
fi

python ${TIGAR_dir}/TWAS/Asso_Study_02.py \
--gene_anno ${gene_anno} \
--Zscore ${Zscore} \
--Zscore_colnames ${out_dir}/CHR${chr}_Zscore_colnames.txt \
--weight ${weight} \
--weight_colnames ${out_dir}/CHR${chr}_weight_colnames.txt \
--LD ${LD} \
--chr ${chr} \
--window ${window} \
--thread ${thread} \
--out_dir ${out_dir}


# rm -f ${out_dir}/CHR${chr}_weight_colnames.txt
# rm -f ${out_dir}/CHR${chr}_Zscore_colnames.txt

