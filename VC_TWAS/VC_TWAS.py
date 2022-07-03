#! /bin/env python

############################################################
# import packages needed
import argparse
import operator
import multiprocessing
import subprocess
import sys

from time import time

import pandas as pd
import numpy as np

##########################################################
# time calculation
start_time = time()

#################################################
# Input Arguments for VC_TWAS 
# --chr: Chromosome number need to be specified with respect to the genotype input data
# --weight: Path for SNP weight (eQTL effect size) file 
# --test_sampleID: Path for a file with sampleIDs that should be contained in the genotype file
# --genofile: Path for the training genotype file (bgzipped and tabixed) 
# --genofile_type: Genotype file type: 'vcf' or 'dosage'
# --test_geno_colnames: File with column heads of genotype file
# --format: Genotype format in VCF file that should be used: 'GT' (default) for genotype data or 'DS' for dosage data, only required if the input genotype file is of VCF file
# --window: Window size around gene transcription starting sites (TSS) for selecting cis-SNPs for fitting gene expression prediction model (default 1000000 for +- 1MB region around TSS)
# --phenotype_type: phenotype type, continous phenotype is 'C' and binomial phenotype is 'D'
# -- maf: threshold for seleting genotype data
# --weight_threshold: threshold for weights estimated from DPR and Elastic Net(absolute value), default threshold is 10^-4
# --thread: Number of threads for parallel computation (default 1)
# --PED: PED file, phenotype 
# --PED_info: association information file for PED
# --out_dir: Output directory (will be created if not exist)

#######################################
# parse input arguments
parser = argparse.ArgumentParser(description='VC TWAS')

# Specify tool directory
parser.add_argument('--TIGAR_dir', type=str)

# Specified chromosome number
parser.add_argument('--chr', type=str, dest='chrm')

# eQTL weight file path
parser.add_argument('--weight', type=str, dest='w_path')

# Test sampleID path
parser.add_argument('--test_sampleID', type=str, dest='sampleid_path')

# Gene annotation file path
parser.add_argument('--gene_anno', type=str, dest='annot_path')

# Test genotype file path
parser.add_argument('--genofile', type=str, dest='geno_path')

# Specified input file type(vcf or dosages)
parser.add_argument('--genofile_type', type=str)

# 'DS' or 'GT' for VCF genotype file
parser.add_argument('--format', type=str, dest='data_format')

# window
parser.add_argument('--window', type=int)

# phenotype_type
parser.add_argument('--phenotype_type', type=str)

# missing rate: threshold for excluding SNPs with too many missing values
parser.add_argument('--missing_rate', type=float)

# maf threshold for seleting genotype data
parser.add_argument('--maf', type=float)

# hwe threshold for seleting genotype data
parser.add_argument('--hwe', type=float)

# weight_threshold
parser.add_argument('--weight_threshold', type=float)

# number of thread
parser.add_argument('--thread', type=int)

# PED file path
parser.add_argument('--PED', type=str, dest='ped_path')

# Association Information file
parser.add_argument('--PED_info', type=str, dest='pedinfo_path')

# output dir
parser.add_argument('--out_dir', type=str)

args = parser.parse_args()

sys.path.append(args.TIGAR_dir)
sys.path.append(args.TIGAR_dir + '/VC_TWAS')

#############################################
# Import TIGAR functions
import TIGARutils as tg
import SKAT

#############################################
# check input arguments
if args.genofile_type == 'vcf':
	if (args.data_format != 'GT') and (args.data_format != 'DS'):
		raise SystemExit('Please specify the genotype data format used by the vcf file (--format ) as either "GT" or "DS".\n')

elif args.genofile_type == 'dosage':
	args.data_format = 'DS'

else:
	raise SystemExit('Please specify the type input genotype file type (--genofile_type) as either "vcf" or "dosage".\n')

if (args.phenotype_type != 'C') and (args.phenotype_type != 'D'):
	raise SystemExit('Please specify phenotype type (--phenotype_type) as either "C" for continous or "D" for dichotomous.\n')

out_twas_path = args.out_dir + '/CHR' + args.chrm + '_indv_VC_TWAS.txt'

###############################################################
# Print input arguments
print(
'''********************************
Input Arguments
Gene annotation file specifying genes for TWAS: {annot_path}
PED phenotype/covariate data file: {ped_path}
PED information file: {pedinfo_path}
Test sampleID file: {sampleid_path}
Chromosome: {chrm}
cis-eQTL weight file: {w_path}
Test genotype file: {geno_path}
Genotype file used for testing is type: {genofile_type}
Genotype data format: {data_format}
Gene testing region SNP inclusion window: +-{window}
MAF threshold for SNP inclusion: {maf}
HWE p-value threshold for SNP inclusion: {hwe}
SNP weight inclusion threshold: {weight_threshold}
{pheno_type_str} phenotype used for SKAT.
Number of threads: {thread}
Output directory: {out_dir}
Output TWAS results file: {out_path}
********************************'''.format(
	**args.__dict__,
	pheno_type_str = {'C':'Continuous', 'D':'Dichotomous'}[args.phenotype_type],
	out_path = out_twas_path))

# tg.print_args(args)

#############################################
## STARTUP 

# sampleID startup
sampleID, sample_size, geno_info, ped_cols, n_pheno, pheno, cov = tg.sampleid_startup(**args.__dict__)

# read in PED file
PED = pd.read_csv(
	args.ped_path,
	sep='\t',
	usecols=['IND_ID', *ped_cols])
PED = PED[PED.IND_ID.isin(sampleID)]
PED = tg.optimize_cols(PED)

# get phenotype, covariate for SKAT test
phenotype = PED[pheno]
covariate = PED[cov]

# read in gene annotation file
print('Reading gene annotation file.')
Gene, TargetID, n_targets = tg.read_gene_annot_exp(**args.__dict__)

# get weight file info
weight_info = tg.weight_file_info(add_cols=['MAF','b','beta'], drop_cols=['ES'], **args.__dict__)

# prep output
print('Creating file: ' + out_twas_path + '\n')

out_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName','n_snps']

if (args.phenotype_type == 'D') and (sample_size < 2000):
	out_cols += ['p_value_resampling','p_value_noadj']
else:
	out_cols += ['p_value']

pd.DataFrame(columns=out_cols).to_csv(
	out_twas_path,
	sep='\t',
	header=True,
	index=None,
	mode='w')

print('********************************\n')

##########################################################
@tg.error_handler
def thread_process(num):
	target = TargetID[num]
	print('num=' + str(num) + '\nTargetID=' + target)
	Gene_Info = Gene.iloc[[num]].reset_index(drop=True)

	start = str(max(int(Gene_Info.GeneStart)-args.window, 0))
	end = str(int(Gene_Info.GeneEnd)+args.window)

	# check that both files have data for target
	tabix_query = tg.tabix_query_files(start, end, **args.__dict__)
	if not tabix_query:
		print('No cis-eQTL weights and/or genotype data for TargetID: ' + target + '\n')
		return None

	# read in weight file
	Weight = tg.read_tabix(start, end, target=target, **weight_info)[['snpID','ES','MAF']]

	# read in genotype file
	print('Reading genotype data.')
	Geno = tg.read_tabix(start, end, sampleID, **geno_info)

	# filter out variants that exceed missing rate threshold
	Geno = tg.handle_missing(Geno, sampleID, args.missing_rate)

	# calculate MAF
	Geno = tg.calc_maf(Geno, sampleID, 0, op=operator.ge)

	# get, filter p_HWE
	Geno = tg.calc_p_hwe(Geno, sampleID, args.hwe)

	# get flipped snpIDs
	Geno['snpIDflip'] = tg.get_snpIDs(Geno, flip=True)

	# Handle overlaps
	snp_overlap = np.intersect1d(Weight.snpID, Geno[['snpID','snpIDflip']])
	n_snps = snp_overlap.size

	if not n_snps:
		print('No overlapping test SNPs between weight and genotype file for TargetID: ' + target + '\n')
		return None

	# filter out non-matching snpID rows
	Weight = Weight[Weight.snpID.isin(snp_overlap)]
	Geno = Geno[np.any(Geno[['snpID','snpIDflip']].isin(snp_overlap), axis=1)].reset_index(drop=True)

	# if not in Weight.snpIDs, assumed flipped; if flipped, 1 - MAF
	flip = np.where(Geno.snpID.isin(Weight.snpID.values), True, False)

	if not np.all(flip):
		# set correct snpID, MAF
		Geno['snpID'] = np.where(flip, Geno.snpID, Geno.snpIDflip)
		Geno['MAF'] = np.where(flip, Geno.MAF, 1 - Geno.MAF)

		# reshape flip for setting sampleIDs, set correct sampleID values
		flip = flip.reshape((Geno.shape[0], 1))
		Geno[sampleID] = np.where(flip, Geno[sampleID], 2 - Geno[sampleID])

	Geno = Geno.drop(columns=['CHROM','POS','REF','ALT','snpIDflip'])

	# filter by MAF
	Geno = Geno[Geno['MAF'] > args.maf]

	if Geno.empty:
		print('No overlapped SNPs exceed the specified MAF threshold of ' + str(args.maf) + ' for TargetID: ' + target + '\n')
		return None

	# MERGE DATAFRAMES
	Geno_Weight = Geno.merge(
		Weight,
		left_on='snpID',
		right_on='snpID',
		how='outer')

	# GENOTYPE MATRIX
	genotype_mat = Geno_Weight[sampleID].T.values

	# WEIGHT VECTOR
	weight_mat = Geno_Weight['ES'].values

	# INITIALIZE OUTPUT DATAFRAME
	Result = Gene_Info.copy()
	Result['n_snps'] = n_snps

	# SKAT TEST
	print('Performing SKAT test.')
	p_value = SKAT.SKAT(
		genotype_mat,
		phenotype,
		covariate,
		weight_mat,
		args.phenotype_type)

	# Pvalue output
	if (args.phenotype_type == 'D') and (len(p_value) == 2):
		Result['p_value_resampling'] = p_value[0]
		Result['p_value_noadj'] = p_value[1]
	else: 
		Result['p_value'] = p_value

	Result['TargetID'] = target

	Result.to_csv(
		out_twas_path,
		sep='\t',
		index=None,
		header=None,
		mode='a')

	print('Target TWAS completed.\n')

##############################################################
# thread process
if __name__ == '__main__':
	print('Starting VC-TWAS for ' + str(n_targets) + ' target genes.\n')
	pool = multiprocessing.Pool(args.thread)
	pool.map(thread_process,[num for num in range(n_targets)])
	pool.close()
	pool.join()
	print('Done.')

############################################################
# time calculation
elapsed_sec = time() - start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)
