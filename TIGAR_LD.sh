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
genome_block:,genofile:,genofile_type:,chr:,format:,maf:,TIGAR_dir:,thread:,sub_dir:,log_file:,out_dir:,out_ld_file:,sampleID: \
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
        --sub_dir|-sub_dir) sub_dir=$2; shift 2;;
        --log_file|-log_file) log_file=$2; shift 2;;
        --out_dir|-out_dir) out_dir=$2; shift 2;;
        --out_ld_file|-out_ld_file) out_ld_file=$2; shift 2;;
        --sampleID|-sampleID) sampleID=$2; shift 2;;
        --) shift;break;;
        *) echo "Please check input arguments!";exit 1;;
        esac
done

#### setting default value
thread=${thread:-1}
maf=${maf:-0}
log_file=${log_file:-CHR${chr}_LD_log.txt}
out_ld_file=${out_ld_file:-CHR${chr}_reference_cov}

# sub_dir: whether to use subdirectory inside out_dir for output files
sub_dir=${sub_dir:-1}

#### Create output directory if not existed
mkdir -p ${out_dir}
mkdir -p ${out_dir}/logs

# sub directory in out directory
if [[ "$sub_dir"x == "1"x ]];then
    out_sub_dir=${out_dir}/RefLD
else
    out_sub_dir=${out_dir}
fi

mkdir -p ${out_sub_dir}

################################################

# check tabix command
if [ ! -x "$(command -v tabix)" ]; then
    echo 'Error: required tool TABIX is not available.' >&2
    exit 1
fi

# Check genotype file 
if [ ! -f "${genofile}" ]; then
    echo Error: Training genotype file ${genofile} does not exist or is empty. >&2
    exit 1
  # check chromosome column to see if it starts with "chr"; will cause "No tabix data" error for all genes
  startswithchr=`zgrep -E '^#' -v ${genofile} | head -n10 | grep -E '^chr' -c`
  # warn but don't stop job?
  if ((startswithchr != 0)); then 
    echo Error: A check of the first 10 non-header lines in your genotype file shows that CHROM column values seem to start with \"chr\"! >&2
    exit 1
  fi    
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
--out_dir ${out_sub_dir} \
--out_ld_file ${out_ld_file} \
--sampleID ${sampleID} \
--TIGAR_dir ${TIGAR_dir} \
> ${out_dir}/logs/${log_file}


### 2. TABIX output file
# Check ld file 
if [ ! -f ${out_sub_dir}/${out_ld_file}.txt ] ; then
    echo Error: Reference LD covariance file ${out_sub_dir}/${out_ld_file}.txt was not generated successfully. >&2
    exit 1
# elif [ ! -f ${out_sub_dir}/${out_ld_file}_block_0.txt.gz ] ; then
elif [ ! -f ${out_sub_dir}/${out_ld_file}_block_0.txt ] ; then
    # if 0th block file does not exist something went wrong
    echo Error: Reference LD covariance block files not generated successfully. >&2
    exit 1
else
    # # bgzip file with header in order to append concatenated, numberd, and bgzipped block files lines from stdout
    # bgzip -f ${out_sub_dir}/${out_ld_file}.txt

    # file should contain header and the output block index starts at 0, but just in case setting nblocks equal to number of lines in the file
    nblocks=$(grep -c $ "${genome_block}")

    # get list of block paths
    block_paths=()
    for block_i in $(seq 0 $nblocks); do
        # block_path=${out_sub_dir}/${out_ld_file}_block_${block_i}.txt.gz
        block_path=${out_sub_dir}/${out_ld_file}_block_${block_i}.txt
        # if block_path exists
        if [ -f ${block_path} ] ; then
            block_paths+=( ${block_path} )
        fi
    done

    echo Concatenating and bgzipping output LD block files.

    start_row=1
    for block_path in ${block_paths[@]}; do
        # # number the lines then add bgzipped lines to final output
        # gunzip -c ${block_path} | \
        #     nl -nln -s$'\t' -v${start_row} | \
        #     bgzip -c >> ${out_sub_dir}/${out_ld_file}.txt.gz
        # # get starting row number for next block
        # start_row=$(( $(zgrep -c $ "${block_path}") + start_row ))

        # number the lines then add bgzipped lines to final output
        nl -nln -s$'\t' -v${start_row} ${block_path} \
            >> ${out_sub_dir}/${out_ld_file}.txt
        # get starting row number for next block
        start_row=$(( $(grep -c $ "${block_path}") + start_row ))

        # remove block input file
        rm -f ${block_path}
    done
    echo 'Bgzipping final output LD file.'
    # bgzip -d ${out_sub_dir}/${out_ld_file}.txt.gz
    bgzip -f ${out_sub_dir}/${out_ld_file}.txt

    echo 'Tabix-ing output LD file.'
    # tabix
    tabix -f -s3 -b4 -e4 -S1 ${out_sub_dir}/${out_ld_file}.txt.gz
fi    


# Check tabix file
if [ ! -f ${out_sub_dir}/${out_ld_file}.txt.gz.tbi ] ; then
    echo Error: Tabix reference LD covariance file ${out_sub_dir}/${out_ld_file}.txt.gz failed. >&2
    exit 1
else 
    echo Successfully generated reference LD covariance file: ${out_sub_dir}/${out_ld_file}.txt.gz
fi



