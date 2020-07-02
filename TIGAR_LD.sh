#!/usr/bin/bash

#######################################################################
### Input Arguments for GReX Prediction
#######################################################################

###
# genom_block : Genome segmentation block annotation based on LD
# --genofile: Path for the reference genotype file (bgzipped and tabixed) 
# --genofile_type: Genotype file type: "vcf" or "dosage"
# --format: Genotype format in VCF file that should be used: "GT" (default) for genotype data or "DS" for dosage data, only required if the input genotype file is of VCF file
# --chr: Chromosome number need to be specified with respect to the genotype input data
# --maf: Minor Allele Frequency threshold (ranges from 0 to 1; default 0.01) to exclude rare variants
# --TIGAR_dir : Specify the directory of TIGAR source code
# --thread: Number of threads for parallel computation (default 1)
# --out_dir: Output directory (will be created if not exist)
###

VARS=`getopt -o "" -a -l \
genome_block:,genofile:,genofile_type:,chr:,format:,maf:,TIGAR_dir:,thread:,out_dir:,sampleID: \
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
        --genome_block|-genome_block) genome_block=$2; shift 2;;
        --genofile|-genofile) genofile=$2; shift 2;;
        --genofile_type|-genofile_type) genofile_type=$2; shift 2;;
        --chr|-chr) chr=$2; shift 2;;
        --format|-format) format=$2; shift 2;;
        --maf|-maf) maf=$2; shift 2;;
        --TIGAR_dir|-TIGAR_dir) TIGAR_dir=$2; shift 2;;
        --thread|-thread) thread=$2; shift 2;;
        --out_dir|-out_dir) out_dir=$2; shift 2;;
        --sampleID|-sampleID) sampleID=$2; shift 2;;
        --) shift;break;;
        *) echo "Please check input arguments!";exit 1;;
        esac
done

#### setting default value
thread=${thread:-1}
maf=${maf:-0}

###
mkdir -p ${out_dir}
# mkdir -p ${out_dir}/RefLD

# check tabix command
if [ ! -x "$(command -v tabix)" ]; then
    echo 'Error: required tool TABIX is not available.' >&2
    exit 1
fi

# Check genotype file 
if [ ! -f "${genofile}" ] ; then
    echo Error: Reference genotype file ${genofile} dose not exist or empty. >&2
    exit 1
fi


################################################
### 1. Calculate covariance matrix (LD) of genetic variants

# Make python script executible
if [[ ! -x ${TIGAR_dir}/TWAS/Get_LD.py ]] ; then
    chmod 755 ${TIGAR_dir}/TWAS/Get_LD.py
fi

python ${TIGAR_dir}/TWAS/Get_LD.py \
--genome_block ${genome_block} \
--genofile ${genofile} \
--genofile_type ${genofile_type} \
--chr ${chr} \
--format ${format} \
--maf ${maf} \
--thread ${thread} \
--out_dir ${out_dir} \
--sampleID ${sampleID}

#/RefLD

### 2. TABIX output file
# Check genotype file 
if [ ! -f ${out_dir}/CHR${chr}_reference_cov.txt ] ; then
    echo Error: Reference LD covariance file ${out_dir}/CHR${chr}_reference_cov.txt was not generated successfully. >&2
    exit 1
else
    sort -n -k2 ${out_dir}/CHR${chr}_reference_cov.txt | bgzip -c > ${out_dir}/CHR${chr}_reference_cov.txt.gz
    tabix -f -p vcf ${out_dir}/CHR${chr}_reference_cov.txt.gz
    rm -f ${out_dir}/CHR${chr}_reference_cov.txt
fi


# Check genotype file 
if [ ! -f ${out_dir}/CHR${chr}_reference_cov.txt.gz.tbi ] ; then
    echo Error: Tabix reference LD covariance file ${out_dir}/CHR${chr}_reference_cov.txt.gz failed. >&2
    exit 1
fi

echo Generate reference LD covariance file successfully ... 










