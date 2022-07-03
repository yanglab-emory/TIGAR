# %%
#! /bin/env python

###################################################################
# Import packages needed
import argparse
import multiprocessing
import subprocess
import sys
from time import time
import numpy as np
import pandas as pd

###############################################################
# time calculation
start_time = time()

###############################################################
# parse input arguments
parser = argparse.ArgumentParser(description='VC_TWAS summary statistics')
### Gene annotation file
parser.add_argument("--gene_anno", type=str, dest='annot_path')

### GWAS result file
parser.add_argument('--GWAS_result', type=str, dest='gwas_path')

### Weight
parser.add_argument('--weight', type=str, dest='w_path')

### sample size 
parser.add_argument('--sample_size', type=int)

### weight threshold
parser.add_argument('--weight_threshold', type=float)

### Reference covariance file
parser.add_argument('--LD', type=str, dest='ld_path')

### chromosome number
parser.add_argument('--chr', type=str, dest='chrm')

### window
parser.add_argument('--window', type=int)

### Number of thread
parser.add_argument('--thread', type=int)

### Output dir
parser.add_argument('--out_dir', type=str)

parser.add_argument('--TIGAR_dir', type=str)

args = parser.parse_args()
sys.path.append(args.TIGAR_dir)
sys.path.append(args.TIGAR_dir + '/VC_TWAS')

import SKAT
import TIGARutils as tg

#############################################################
# Print input arguments to log
out_sum_VCTWAS_path = args.out_dir + '/CHR' + args.chrm + '_sum_VC_TWAS.txt'

print(
'''********************************
Input Arguments
Gene annotation file specifying genes for VC-TWAS: {annot_path}
GWAS summary statistics file: {gwas_path}
cis-eQTL weight file: {w_path}
sample—size:{sample_size}
SNP weight inclusion threshold:{weight_threshold}
Reference LD genotype covariance file: {ld_path}
Chromosome: {chrm}
Gene training region SNP inclusion window: +-{window}
Number of threads: {thread}
Output directory: {out_dir}
Output TWAS results file: {out_path}
********************************'''.format(
	**args.__dict__,
	out_path = out_sum_VCTWAS_path))

# tg.print_args(args)

###############################################################
# read in gene annotation file
print('Reading gene annotation file.')
Gene, TargetID, n_targets = tg.read_gene_annot_exp(**args.__dict__)

# read in headers for Weight, GWAS files
print('Reading file headers.\n')
weight_info = tg.weight_file_info(add_cols=['MAF','b','beta'], drop_cols=['ES'], **args.__dict__)
gwas_info = tg.gwas_file_info(**args.__dict__)

# PREP OUTPUT - print output headers to files
print('Creating file: ' + out_sum_VCTWAS_path + '\n')
out_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName','n_snps','Pvalue']
pd.DataFrame(columns=out_cols).to_csv(
	out_sum_VCTWAS_path,
	sep='\t',
	index=None,
	header=True,
	mode='w')

###############################################################
# thread function
@tg.error_handler
def thread_process(num):
	target = TargetID[num]
	print('num=' + str(num) + '\nTargetID=' + target)
	Gene_Info = Gene.iloc[[num]].reset_index(drop=True)

	# get start and end positions to tabix
	start = str(max(int(Gene_Info.GeneStart)-args.window,0))
	end = str(int(Gene_Info.GeneEnd)+args.window)

	# check that both files have data for target
	tabix_query = tg.tabix_query_files(start, end, **args.__dict__)
	if not tabix_query:
		print('No cis-eQTL weights and/or genotype data for TargetID: ' + target + '\n')
		return None

	# tabix Weight file
	Weight = tg.read_tabix(start, end, target=target, **weight_info)

	# tabix gwas result file
	GWAS = tg.read_tabix(start, end, **gwas_info)

	# get flipped snpIDs
	GWAS['snpIDflip'] = tg.get_snpIDs(GWAS, flip=True)

	snp_overlap = np.intersect1d(Weight.snpID, GWAS[['snpID','snpIDflip']])

	if not snp_overlap.size:
		print('No overlapping test SNPs that have magnitude of cis-eQTL weights greater than threshold value and with GWAS results for TargetID: ' + target + '.\n')
		return None

	# filter out non-matching snpID rows
	Weight = Weight[Weight.snpID.isin(snp_overlap)]
	GWAS = GWAS[np.any(GWAS[['snpID','snpIDflip']].isin(snp_overlap), axis=1)].reset_index(drop=True)	

	flip = np.where(GWAS.snpID.isin(Weight.snpID), 1, -1)

	# if not in Weight.snpIDs, assumed flipped; if flipped, flip Beta sign
	if not np.all(flip == 1):
		GWAS['snpID'] = np.where(flip == 1, GWAS.snpID, GWAS.snpIDflip)
		GWAS['BETA'] = flip * GWAS['BETA']

	# drop unneeded columns
	GWAS = GWAS.drop(columns=['CHROM','POS','REF','ALT','snpIDflip'])

	GW = Weight.merge(GWAS, 
		left_on='snpID', 
		right_on='snpID')

	snp_search_ids = GW.snpID

	# Read in reference covariance matrix file by snpID
	MCOV = tg.get_ld_data(args.ld_path, snp_search_ids)

	if MCOV.empty:
	  print('No reference covariance information for target SNPs for TargetID: ' + target + '\n')
	  return None

	snp_sd, V, D = tg.get_ld_matrix(MCOV, return_diag=True)

	# filter GW to include only snpIDs also in MCOV
	GW = GW[GW.snpID.isin(MCOV.snpID)]
	n_snps = str(GW.snpID.size)

	#prepare input for Summary statistics
	beta_var = np.power(GW['SE'].values, 2)
	beta_estimate = GW['BETA'].values
	weight_temp = GW['ES'].values

	###VC-TWAS SUMMARY
	p_val = SKAT.SKAT_summary(beta_var, beta_estimate, weight_temp, args.sample_size, V, D)

	### create output dataframe
	Result = Gene_Info.copy()
	Result['n_snps'] = n_snps
	Result['Pvalue'] = p_val
	
	# write to file
	Result.to_csv(
		out_sum_VCTWAS_path,
		sep='\t',
		index=None,
		header=None,
		mode='a')

	print('Target VC-TWAS summary completed.\n')

###############################################################
# thread process
if __name__ == '__main__':
	print('Starting VC-TWAS summary statistics with GWAS result for ' + str(n_targets) + ' target genes.\n')
	pool = multiprocessing.Pool(args.thread)
	pool.imap(thread_process,[num for num in range(n_targets)])
	pool.close()
	pool.join()
	print('Done.')

###############################################################
# time calculation
elapsed_sec = time()-start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)
