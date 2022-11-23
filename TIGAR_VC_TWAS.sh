#!/usr/bin/bash

#######################################################################
### Input Arguments for VC-TWAS
#######################################################################

# --chr: Chromosome number need to be specified with respect to the genotype input data
# --weight: Path for SNP weight (eQTL effect size) file 
# --gene_anno : Path for gene annotation file to specify a list of gene for GReX prediction
# `--gene_anno`: Gene annotation file to specify the list of genes, which is of the same format as the first five columns of gene expression file
# --test_sampleID: Path for a file with sampleIDs that should be contained in the genotype file
# --genofile: Path for the training genotype file (bgzipped and tabixed) 
# --genofile_type: Genotype file type: "vcf" or "dosage"
# --format: Genotype format in VCF file that should be used: "GT" (default) for genotype data or "DS" for dosage data, only required if the input genotype file is of VCF file
# --window: Window size around gene region for selecting cis-SNPs for fitting gene expression prediction model (default 1000000 for +- 1MB region around gene)
# --phenotype_type: phenotype type, continous phenotype is "C" and binomial phenotype is "D"
# --maf: Threshold for Minor Allele Frequency (range from 0-1),default 0.01
# --hwe: Hardy Weinberg Equilibrium p-value threshold (default 0.00001) to exclude variants that violated HWE
# --weight_threshold: Threshold for estimated cis-eQTL effect sizes, filter SNPs with absolute cis-eQTL effect sizes smaller than threshold
# --TIGAR_dir : Specify the directory of TIGAR source code
# --thread: Number of threads for parallel computation (default 1)
# --PED : Path for PED file that contains phenotype and covariate data
# --PED_info : Specify culumn names for phenotypes and covariates that will be used in TWAS.
#             1) P : phenotype column names
#             2) C : covariate column names
# --out_dir: Output directory (will be created if not exist)


#######################################################################
VARS=`getopt -o "" -a -l \
chr:,weight:,test_sampleID:,gene_anno:,genofile:,genofile_type:,format:,window:,phenotype_type:,missing_rate:,maf:,hwe:,weight_threshold:,PED:,PED_info:,TIGAR_dir:,thread:,sub_dir:,log_file:,in_dir:,out_dir: \
-- "$@"`

echo "Conducting VC-TWAS using individual-level GReX and phenotype data ... "

if [ $? != 0 ]
then
    echo "Terminating....." >&2
    exit 1
fi
 
eval set -- "$VARS"

while true
do
    case "$1" in
        --chr|-chr) chr=$2; shift 2;;
        --weight|-weight) weight_file=$2; shift 2;;
        --test_sampleID|-test_sampleID) test_sampleID_file=$2; shift 2;;
        --gene_anno|-gene_anno) gene_anno=$2; shift 2;;
        --genofile|-genofile) genofile=$2; shift 2;;
        --genofile_type|-genofile_type) genofile_type=$2; shift 2;;
        --format|-format) format=$2; shift 2;;
        --window|-window) window=$2; shift 2;;
        --phenotype_type|-phenotype_type) phenotype_type=$2; shift 2;;
        --missing_rate|-missing_rate) missing_rate=$2; shift 2;;
        --maf|-maf) maf=$2; shift 2;;
        --hwe|-hwe) hwe=$2; shift 2;;
        --weight_threshold|-weight_threshold) weight_threshold=$2; shift 2;;
        --PED|-PED) PED=$2; shift 2;;
        --PED_info|-PED_info) PED_info=$2; shift 2;;
        --TIGAR_dir|-TIGAR_dir) TIGAR_dir=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --sub_dir|-sub_dir) sub_dir=$2; shift 2;;
        --log_file|-log_file) log_file=$2; shift 2;;
        --out_dir|-out_dir) out_dir=$2; shift 2;;
				--in_dir|-in_dir) in_dir=$2; shift 2;;
        --) shift;break;;
        *) echo "Internal error!";exit 1;;
        esac
done

#### default value
thread=${thread:-1}
window=${window:-$((10**6))}
missing_rate=${missing_rate:-0.2}
maf=${maf:-0.01}
hwe=${hwe:-0.00001}
log_file=${log_file:-CHR${chr}_VCTWAS_log.txt}

# sub_dir: whether to use subdirectory inside out_dir for output files
sub_dir=${sub_dir:-1}

# check if user submitted in_dir
if [[ "$in_dir"x != ""x ]];then
	# if yes, check if in_dir var ends with a backslash
  if [[ "$in_dir"x != */x ]];then
  	# if it doesn't, add backslash
    in_dir=$in_dir"/"
  fi
fi

#### Create output directory if not existed
mkdir -p ${out_dir}
mkdir -p ${out_dir}/logs

# sub directory in out directory
if [[ "$sub_dir"x == "1"x ]];then
    out_sub_dir=${out_dir}/VCTWAS_CHR${chr}
else
    out_sub_dir=${out_dir}
fi

mkdir -p ${out_sub_dir}
# mkdir -p ${out_dir}/VC_TWAS_CHR${chr}

####################################################
# check tabix command
if [ ! -x "$(command -v tabix)" ]; then
    echo 'Error: required tool TABIX is not available.' >&2
    exit 1
fi

# Check genotype file 
if [ ! -f "${in_dir}${genofile}" ]; then
    echo Error: Training genotype file ${in_dir}${genofile} does not exist or is empty. >&2
    exit 1
fi

# Check training sample ID file
if [ ! -f "${in_dir}${test_sampleID_file}" ]; then
    echo Error: Test sample ID file ${in_dir}${test_sampleID_file} does not exist or is empty. >&2
    exit 1
fi

# Check eQTL weight file
if [ ! -f "${in_dir}${weight_file}" ]; then
    echo Error: eQTL weight file ${in_dir}${weight_file} does not exist or is empty. >&2
    exit 1
fi

# Check gene annotation file
if [ ! -f "${in_dir}${gene_anno}" ]; then
    echo Error: eQTL weight file ${in_dir}${gene_anno} does not exist or is empty. >&2
    exit 1
fi

# Make python script executible
# if [[ ! -x  ${TIGAR_dir}/VC_TWAS/VC_TWAS_ver1.py ]] ; then
#     chmod 755  ${TIGAR_dir}/VC_TWAS/VC_TWAS.py
# fi

# Make python script executible
if [[ ! -x  ${TIGAR_dir}/VC_TWAS/VC_TWAS.py ]]; then
    chmod 755 ${TIGAR_dir}/VC_TWAS/VC_TWAS.py
fi

python ${TIGAR_dir}/VC_TWAS/VC_TWAS.py \
--chr ${chr} \
--weight ${in_dir}${weight_file} \
--test_sampleID ${in_dir}${test_sampleID_file} \
--gene_anno ${in_dir}${gene_anno} \
--genofile ${in_dir}${genofile} \
--genofile_type ${genofile_type} \
--format ${format} \
--window ${window} \
--phenotype_type ${phenotype_type} \
--missing_rate ${missing_rate} \
--maf ${maf} \
--hwe ${hwe} \
--weight_threshold ${weight_threshold} \
--thread ${thread} \
--PED ${in_dir}${PED} \
--PED_info ${in_dir}${PED_info} \
--TIGAR_dir ${TIGAR_dir} \
--out_dir ${out_sub_dir} \
> ${out_dir}/logs/${log_file}

echo "Completed TWAS."


