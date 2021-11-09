# %%
#!/usr/bin/env python

###################################################################
# Import packages needed
import argparse
import multiprocessing
import os
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

parser.add_argument('--gene_anno', type=str, dest='annot_path', required=True)
parser.add_argument('--chr', type=str, dest='chrm', required=True, 
	help='chromosome number')
parser.add_argument('--GWAS_result', type=str, dest='gwas_path', required=True)
parser.add_argument('--LD', type=str, dest='ld_path', required=True)
parser.add_argument('--log_file', type=str, default='')
parser.add_argument('--out_dir', type=str, default=os.getcwd())
parser.add_argument('--out_prefix', type=str, default='')
parser.add_argument('--out_twas_file', type=str, default='')
parser.add_argument('--sample_size', type=int, required=True)
parser.add_argument('--sub_dir', type=int, choices=[0, 1], default=1)
parser.add_argument('--thread', type=int, default=1)
parser.add_argument('--TIGAR_dir', type=str, help='tool directory', 
	default=os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
parser.add_argument('--weight', type=str, dest='w_path', required=True, 
	help='path to SNP weight (eQTL effect size) file; output file _eQTLweights.txt from model training can be used here')
parser.add_argument('--weight_threshold', type=float, default=0, 
	help='weight magnitude threshold for SNP inclusion; include only SNPs with magnitude of weight greater than this value when conducting TWAS(default: 0 [all SNPs included])')
parser.add_argument('--window', type=int, default=1000000, 
	help='size around gene region for selecting cis-SNPs for fitting gene expression prediction model (default: 1000000 [ie, +-1MB region around gene])')

args = parser.parse_args()
sys.path.append(args.TIGAR_dir)
sys.path.append(args.TIGAR_dir + '/VC_TWAS')

import SKAT
import TIGARutils as tg

#############################################################
# set output file names
if not args.out_prefix:
	args.out_prefix = 'CHR' + args.chrm + '_VCTWAS_sum'

if not args.log_file:
	args.log_file = args.out_prefix + '_log.txt'

if not args.out_twas_file:
	args.out_twas_file = args.out_prefix + '_assoc.txt'

# sub-directory in out directory
if args.sub_dir:
	out_sub_dir = os.path.join(args.out_dir, 'VCTWAS_sum_CHR' + args.chrm)
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


#############################################################
# Print input arguments to log

# out_twas_path = args.out_dir + '/CHR' + args.chrm + '_sum_VC_TWAS.txt'
tmp_twas_path = out_sub_dir + '/temp_' + args.out_twas_file
out_twas_path = out_sub_dir + '/' + args.out_twas_file


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
	out_path = out_twas_path))

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
print('Creating file: ' + tmp_twas_path + '\n')
out_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName','n_snps','Pvalue']
pd.DataFrame(columns=out_cols).to_csv(
	tmp_twas_path,
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
		tmp_twas_path,
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

	# sort output
	tg.sort_tabix_output(tmp_twas_path, out_twas_path, do_tabix=0)

###############################################################
# time calculation
elapsed_sec = time()-start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)

sys.stdout.close()
