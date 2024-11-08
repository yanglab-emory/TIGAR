
TIGAR_dir=/Users/randyparr/github/TIGAR
gene_anno_file="${TIGAR_dir}/ExampleData/gene_anno.txt"
eQTL_ES_file="${TIGAR_dir}/ExampleData/eQTLweights.txt.gz"
Zscore_file="${TIGAR_dir}/ExampleData/CHR1_GWAS_Zscore.txt.gz"
LD_file="${TIGAR_dir}/ExampleData/CHR1_reference_cov.txt.gz"

out_dir=/Users/randyparr/tigartesting

chmod 755 ${TIGAR_dir}/TIGAR_TWAS.sh

${TIGAR_dir}/TIGAR_TWAS.sh \
--asso 2 \
--gene_anno ${gene_anno_file} \
--chr 1 \
--weight ${eQTL_ES_file} \
--Zscore ${Zscore_file} \
--LD ${LD_file} \
--thread 1 \
--out_dir ${out_dir} \
--TIGAR_dir ${TIGAR_dir}

/usr/bin/bash: bad interpreter: No such file or directory

gene_anno=${gene_anno_file}
weight=${eQTL_ES_file}
Zscore=${Zscore_file}
LD=${LD_file}
asso=2
chr=1

######



out_dir=/Users/randyparr/tigartesting

# chmod 755 ${TIGAR_dir}/TIGAR_TWAS.sh
# chmod 755 ${TIGAR_dir}/TWAS/*.py

TIGAR_dir=/Users/randyparr/github/TIGAR
gene_exp="${TIGAR_dir}/ExampleData/CHR1_Pred_GReX.txt"
# PED="${TIGAR_dir}/ExampleData/example_PED_binary.ped"
PED="${TIGAR_dir}/ExampleData/example_PED.ped"
PED_info="${TIGAR_dir}/ExampleData/PED_Info_SinglePheno.txt"
method="Logit"
asso="1"
thread=1
out_sub_dir="/Users/randyparr/tigartesting"


out_prefix=indv_${method}_TWAS
out_twas_file=${out_prefix}_assoc.txt
log_file=${out_prefix}_log.txt

${TIGAR_dir}/TIGAR_TWAS.sh \
--asso ${asso} \
--gene_exp ${gene_exp} \
--PED ${PED} \
--PED_info ${PED_info} \
--method ${method} \
--thread ${thread} \
--out_dir ${out_dir} \
--sub_dir 0 \
--TIGAR_dir ${TIGAR_dir}



out_dir=/Users/randyparr/tigartesting
TIGAR_dir=/Users/randyparr/github/TIGAR
genofile="${TIGAR_dir}/ExampleData/example.vcf.gz"
gene_anno="${TIGAR_dir}/ExampleData/gene_anno.txt"
test_sampleID="${TIGAR_dir}/ExampleData/test_sampleID.txt"
weight="${TIGAR_dir}/ExampleData/eQTLweights.txt.gz"
chr=1
weight_file=${weight}
test_sampleID_file=${test_sampleID}
format=GT
genofile_type=vcf

${TIGAR_dir}/TIGAR_GReX_Pred.sh \
--gene_anno ${gene_anno} \
--test_sampleID ${test_sampleID} \
--chr 1 \
--weight ${weight} \
--genofile ${genofile} \
--genofile_type vcf \
--format GT \
--thread 2 \
--out_dir ${out_dir} \
--TIGAR_dir ${TIGAR_dir}


