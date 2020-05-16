#!/usr/bin/bash

#######################################################################
### Input Arguments for TWAS
#######################################################################
# --asso : Specify TWAS type either by using individual-level phenotype data and predicted GReX with value `1` (default) or using GWAS summary Z-score statistics with value `2`
# --thread: Number of threads for parallel computation (default 1)
# --out_dir: Output directory (will be created if not exist)

############### Variables for asso=1 ###############
# --gene_exp: Path for predicted GReX file

# --PED : Path for PED file that contains phenotype and covariate data
# --PED_info : Specify culumn names for phenotypes and covariates that will be used in TWAS.
#             1) P : phenotype column names
#             2) C : covariate column names
# --method : `OLS` (default) for studying quantitative phenotype or multivariate phenotypes or `Logit` for studying dichotomous univariate phenotype

############### Variables for asso=2 ###############
# --gene_anno : Path for gene annotation file to specify a list of gene for GReX prediction
# --chr: Chromosome number need to be specified 
# --weight: Path for SNP weight (eQTL effect size) file 
# --Zscore : Path for GWAS summary Zscore statistics
# --LD : Path for reference LD (SNP genotype covariance matrix) that should be bgzipped and tabixed. Can be generated using our `covar_calculation.py` script
# --window: Window size around gene transcription starting sites (TSS) for selecting cis-SNPs for fitting gene expression prediction model (default 1000000 for +- 1MB region around TSS)
# --TIGAR_dir : Specify the directory of TIGAR source code


###############################################################
VARS=`getopt -o "" -a -l \
asso:,gene_exp:,gene_anno:,PED:,PED_info:,method:,Zscore:,weight:,LD:,chr:,window:,TIGAR_dir:,thread:,out_dir: \
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
        --asso|-asso) asso=$2; shift 2;;
        --gene_exp|-gene_exp) gene_exp=$2; shift 2;;
		--gene_anno|-gene_anno) gene_anno=$2; shift 2;;
        --PED|-PED) PED=$2; shift 2;;
        --PED_info|-PED_info) PED_info=$2; shift 2;;
        --method|-method) method=$2; shift 2;;
        --Zscore|-Zscore) Zscore=$2; shift 2;;
        --weight|-weight) weight=$2; shift 2;;
        --LD|-LD) LD=$2; shift 2;;
        --chr|-chr) chr_num=$2; shift 2;;
        --window|-window) window=$2; shift 2;;
        --TIGAR_dir|-TIGAR_dir) TIGAR_dir=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --out_dir|-out_dir) out_dir=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done

########## Set default value
thread=${thread:-1}
window=${window:-$((10**6))}
method=${method:-'OLS'}

############# TWAS 

if [[ "$asso"x == "1"x ]];then
    echo "Conducting TWAS using individual-level GReX and phenotype data ... "

    # Check gene expression file
    if [ ! -f "${gene_exp}" ] ; then
        echo Error: Gene expression file ${gene_exp} dose not exist or empty. >&2
        exit 1
    fi

    # Check PED file
    if [ ! -f "${PED}" ] ; then
        echo Error: PED file ${PED} dose not exist or empty. >&2
        exit 1
    fi

    # Check PED_info file
    if [ ! -f "${PED_info}" ] ; then
        echo Error: PED information file ${PED_info} dose not exist or empty. >&2
        exit 1
    fi

    ## make output directory
    mkdir -p ${out_dir}/TWAS_Ind

    #### TWAS
    if [[ ! -x  ${TIGAR_dir}/TWAS/Asso_Study_01.py ]] ; then
        chmod 755  ${TIGAR_dir}/TWAS/Asso_Study_01.py
    fi

    ${TIGAR_dir}/TWAS/Asso_Study_01.py \
    --gene_exp ${gene_exp} \
    --PED ${PED} \
    --PED_info ${PED_info} \
    --method ${method} \
    --thread ${thread} \
    --out_dir ${out_dir}/TWAS_Ind

elif [[ "$asso"x == "2"x ]];then
    echo "Conducting TWAS using summary-level GWAS Z-score statistics and reference LD covariance file ... "

    # Check gene_annotation file
    if [ ! -f "${gene_anno}" ] ; then
        echo Error: Gene expression file ${gene_exp} dose not exist or empty. >&2
        exit 1
    fi

    # Check LD file
    if [ ! -f "${LD}" ] ; then
        echo Error: Reference LD genotype covariance file ${LD} dose not exist or empty. >&2
        exit 1
    fi

    # Check Zscore file
    if [ ! -f "${Zscore}" ] ; then
        echo Error: Gene expression file ${Zscore} dose not exist or empty. >&2
        exit 1
    else
        zcat ${Zscore} | head -n1 > ${out_dir}/CHR${chr}_Zscore_colnames.txt
    fi

    # Check weight file
    if [ ! -f "${weight}" ] ; then
        echo Error: Gene expression file ${weight} dose not exist or empty. >&2
        exit 1
    else
        zcat ${weight} | head -n1 > ${out_dir}/CHR${chr}_weight_colnames.txt
    fi

    if [[ ! -x  ${TIGAR_dir}/TWAS/Asso_Study_02.sh ]] ; then
        chmod 755  ${TIGAR_dir}/TWAS/Asso_Study_02.sh
    fi

    ## make output directory
    mkdir -p ${out_dir}/TWAS_Sum

    ## TWAS
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
    --out_dir ${out_dir}/TWAS_Sum

fi


