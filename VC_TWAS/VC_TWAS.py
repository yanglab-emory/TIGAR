#!/usr/bin/env python

############################################################
# import packages needed
import argparse
import multiprocessing
import operator
import os
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

parser.add_argument('--chr', type=str, dest='chrm', 
	choices=[str(i + 1) for i in range(22)],
	required=True, 
	help='chromosome number')
parser.add_argument('--format', type=str, dest='data_format', choices=['GT', 'DS'], default='GT', 
	help='data format of VCF genotype data (DS, GT [default])')
parser.add_argument('--gene_anno', type=str, dest='annot_path', required=True)
parser.add_argument('--genofile', type=str, dest='geno_path', required=True)
parser.add_argument('--genofile_type', type=str, choices=['vcf', 'dosage'], 
	help='filetype of genofile (vcf, dosages)')
parser.add_argument('--hwe', type=float, default=0.00001, 
	help='threshold p-value for Hardy Weinberg Equilibrium exact test (default: 0.00001)')
parser.add_argument('--log_file', type=str, default='')
parser.add_argument('--maf', type=float, default=0.01, 
	help='folded Minor Allele Frequency threshold; range from 0-0.5 (default: 0.01)')
parser.add_argument('--missing_rate', type=float, default=0.2, 
	help='missing rate threshold for excluding SNPs with too many missing values (default: 0.2)')
parser.add_argument('--out_dir', type=str, default=os.getcwd())
parser.add_argument('--out_prefix', type=str, default='')
parser.add_argument('--out_twas_file', type=str, default='')
parser.add_argument('--PED', type=str, dest='ped_path', required=True)
parser.add_argument('--PED_info', type=str, dest='pedinfo_path', required=True)
parser.add_argument('--phenotype_type', type=str, choices=['C','D'], 
	help='phenotype type (C: continuous, D: dichotomous/binomial' )
parser.add_argument('--sub_dir', type=int, choices=[0, 1], default=1)
parser.add_argument('--test_sampleID', type=str, dest='sampleid_path', required=True)
parser.add_argument('--thread', type=int, default=1)
parser.add_argument('--TIGAR_dir', type=str, help='tool directory', 
	default=os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
parser.add_argument('--weight', type=str, dest='w_path', required=True)
parser.add_argument('--weight_threshold', type=float, default=0.0001,
	help='weight magnitude threshold for SNP inclusion; include only SNPs with magnitude of weight greater than this value when conducting TWAS(default: 0.0001 [SNPs with |weight|>0.0001 included])')
parser.add_argument('--window', type=int, default=1000000, 
	help='size around gene region for selecting cis-SNPs for fitting gene expression prediction model (default: 1000000 [ie, +-1MB region around gene])')

args = parser.parse_args()
sys.path.append(args.TIGAR_dir)
sys.path.append(args.TIGAR_dir + '/VC_TWAS')

#############################################
# Import TIGAR functions
import TIGARutils as tg
import SKAT

#############################################################
# set output file names
if not args.out_prefix:
	args.out_prefix = 'CHR' + args.chrm + '_VCTWAS'

if not args.log_file:
	args.log_file = args.out_prefix + '_log.txt'

if not args.out_twas_file:
	args.out_twas_file = args.out_prefix + '_indv_assoc.txt'

# sub-directory in out directory
if args.sub_dir:
	out_sub_dir = os.path.join(args.out_dir, 'VCTWAS_CHR' + args.chrm)
else:
	out_sub_dir = args.out_dir
out_sub_dir = tg.get_abs_path(out_sub_dir)

# Check tabix command
tg.check_tabix()

# Check input files
tg.check_input_files(args)

# Make output, log directories
os.makedirs(out_sub_dir, exist_ok=True)
os.makedirs(os.path.join(args.out_dir, 'logs'), exist_ok=True)

# set stdout to log
sys.stdout = open(os.path.join(args.out_dir, 'logs', args.log_file), 'w')

#############################################
tmp_twas_path = out_sub_dir + '/temp_' + args.out_twas_file
out_twas_path = out_sub_dir + '/' + args.out_twas_file

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
print('Creating file: ' + tmp_twas_path + '\n')

out_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName','n_snps']

if (args.phenotype_type == 'D') and (sample_size < 2000):
	out_cols += ['p_value_resampling','p_value_noadj']
else:
	out_cols += ['p_value']

pd.DataFrame(columns=out_cols).to_csv(
	tmp_twas_path,
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
		tmp_twas_path,
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

	tg.sort_tabix_output(tmp_twas_path, out_twas_path, do_tabix=0)

############################################################
# time calculation
elapsed_sec = time() - start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)

sys.stdout.close()
