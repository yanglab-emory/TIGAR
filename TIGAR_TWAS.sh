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
# --window: Window size around gene region for selecting cis-SNPs for fitting gene expression prediction model (default 1000000 for +- 1MB region around gene)
# --weight_threshold : for asso=2, only include SNPs with magnitude of weight greater than the weight_threshold value; default is 0
# --TIGAR_dir : Specify the directory of TIGAR source code


###############################################################
VARS=`getopt -o "" -a -l \
asso:,gene_exp:,gene_anno:,PED:,PED_info:,method:,Zscore:,weight:,LD:,chr:,window:,TIGAR_dir:,thread:,weight_threshold:,sub_dir:,out_twas_file:,out_prefix:,log_file:,out_dir:,sampleID:,test_stat:,test_sampleID: \
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
        --chr|-chr) chr=$2; shift 2;;
        --window|-window) window=$2; shift 2;;
        --TIGAR_dir|-TIGAR_dir) TIGAR_dir=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --weight_threshold|-weight_threshold) weight_threshold=$2; shift 2;;
        --test_stat|-test_stat) test_stat=$2; shift 2;;
        --sampleID|-sampleID) sampleID=$2; shift 2;;
        --test_sampleID|-test_sampleID) test_sampleID=$2; shift 2;;
        --sub_dir|-sub_dir) sub_dir=$2; shift 2;;
        --out_prefix|-out_prefix) out_prefix=$2; shift 2;;
        --out_twas_file|-out_twas_file) out_twas_file=$2; shift 2;;
        --log_file|-log_file) log_file=$2; shift 2;;
        --out_dir|-out_dir) out_dir=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done

########## Set default value
thread=${thread:-1}
window=${window:-$((10**6))}
method=${method:-'OLS'}
weight_threshold=${weight_threshold:-0}
test_stat=${test_stat:-'both'}
sub_dir=${sub_dir:-1}

# sub directory in out directory
if [[ "$sub_dir"x == "1"x ]];then
    out_sub_dir=${out_dir}/TWAS_CHR${chr}
else
    out_sub_dir=${out_dir}
fi

############# TWAS 


## make output directory
mkdir -p ${out_dir}
mkdir -p ${out_dir}/logs
mkdir -p ${out_sub_dir}

if [[ "$asso"x == "1"x ]];then
    echo "Conducting TWAS using individual-level GReX and phenotype data ... "

    out_prefix=${out_prefix:-indv_${method}_TWAS}
    out_twas_file=${out_twas_file:-${out_prefix}_assoc.txt}
    log_file=${log_file:-${out_prefix}_log.txt}

    # Check gene expression file
    if [ ! -f "${gene_exp}" ] ; then
        echo Error: Gene expression file ${gene_exp} does not exist or is empty. >&2
        exit 1
    fi

    # Check PED file
    if [ ! -f "${PED}" ] ; then
        echo Error: PED file ${PED} does not exist or is empty. >&2
        exit 1
    fi

    # Check PED_info file
    if [ ! -f "${PED_info}" ] ; then
        echo Error: PED information file ${PED_info} does not exist or is empty. >&2
        exit 1
    fi

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
    --TIGAR_dir ${TIGAR_dir} \
    --out_dir ${out_sub_dir} \
    --out_twas_file ${out_twas_file} \
    > ${out_dir}/logs/${log_file}

elif [[ "$asso"x == "2"x ]];then
    echo "Conducting TWAS using summary-level GWAS Z-score statistics and reference LD covariance file ... "

    out_prefix=${out_prefix:-CHR${chr}_sumstat_TWAS}
    out_twas_file=${out_twas_file:-${out_prefix}_assoc.txt}
    log_file=${log_file:-${out_prefix}_log.txt}

    # Check gene_annotation file
    if [ ! -f "${gene_anno}" ] ; then
        echo Error: Gene annotation file ${gene_anno} does not exist or is empty. >&2
        exit 1
    fi

    # Check LD file
    if [ ! -f "${LD}" ] ; then
        echo Error: Reference LD genotype covariance file ${LD} does not exist or is empty. >&2
        exit 1
    fi

    # Check Zscore file
    if [ ! -f "${Zscore}" ] ; then
        echo Error: Zscore file ${Zscore} does not exist or is empty. >&2
        exit 1
    fi

    # Check weight file and tabix weight file
    if [ ! -f "${weight}" ] ; then
        echo Error: Weight file ${weight} does not exist or is empty. >&2
        exit 1
    fi

    if [[ ! -x  ${TIGAR_dir}/TWAS/Asso_Study_02.py ]] ; then
        chmod 755  ${TIGAR_dir}/TWAS/Asso_Study_02.py
    fi

    ## TWAS
    python ${TIGAR_dir}/TWAS/Asso_Study_02.py \
    --gene_anno ${gene_anno} \
    --Zscore ${Zscore} \
    --weight ${weight} \
    --LD ${LD} \
    --chr ${chr} \
    --window ${window} \
    --thread ${thread} \
    --weight_threshold ${weight_threshold} \
    --test_stat ${test_stat} \
    --out_dir ${out_sub_dir} \
    --TIGAR_dir ${TIGAR_dir} \
    --out_twas_file ${out_twas_file} \
    > ${out_dir}/logs/${log_file}


fi

echo "Completed TWAS."

# ALSO TRY MOVING FINAL TWAS FILE DIRECTLY TO OUT DIRECTORY INSTEAD OF HAVING A WHOLE FOLDER FOR 1 FILE?
