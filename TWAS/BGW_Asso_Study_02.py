#!/usr/bin/env python

###################################################################
# Import packages needed
import argparse
import multiprocessing
import os
import subprocess
import sys

from io import StringIO
from time import time

import numpy as np
import pandas as pd

from scipy.stats import chi2

###############################################################
# time calculation
start_time = time()

###############################################################
# parse input arguments
parser = argparse.ArgumentParser(description='Asso Study 02')

# Specify tool directory
parser.add_argument('--TIGAR_dir', type=str)

# Gene annotation file path
parser.add_argument('--gene_anno', type=str, dest='annot_path')

# # chromosome number
# parser.add_argument('--chr', type=str, dest='chrm')

# Weight file path
parser.add_argument('--weight_prefix', type=str, dest='w_path_pre')
parser.add_argument('--weight_suffix', type=str, dest='w_path_suf')

# GWAS Z score file path
# parser.add_argument('--Zscore', type=str, dest='z_path')

parser.add_argument('--Zscore_prefix', type=str, dest='z_path_pre')
parser.add_argument('--Zscore_suffix', type=str, dest='z_path_suf')


# Reference covariance file path
parser.add_argument('--LD_prefix', type=str, dest='ld_path_pre')
parser.add_argument('--LD_suffix', type=str, dest='ld_path_suf')

# parser.add_argument('--LD', type=str, dest='ld_path')

# window
# parser.add_argument('--window',type=float)

# Weight threshold to include SNP in TWAS
parser.add_argument('--weight_threshold',type=float)

# specify 'FUSION', 'SPrediXcan', or 'both': Zscore test statistic to use
parser.add_argument('--test_stat', type=str)

# Number of threads
parser.add_argument('--thread',type=int)

# Output dir
parser.add_argument('--out_dir', type=str)

# output file
parser.add_argument('--out_twas_file', type=str)

args = parser.parse_args()

sys.path.append(args.TIGAR_dir)

###############################################################
## Import TIGAR functions, define other functions
import TIGARutils as tg

def get_pval(z): return np.format_float_scientific(1-chi2.cdf(z**2, 1), precision=15, exp_digits=0)

def get_V_cor(V_cov):
	V_cov = V_cov.copy()
	v = np.sqrt(np.diag(V_cov))
	outer_v = np.outer(v, v)
	V_cor = V_cov / outer_v
	V_cor[V_cov == 0] = 0
	return V_cor

def get_z_denom(V, w):
	return np.sqrt(np.linalg.multi_dot([w, V, w]))

def get_spred_zscore(V_cov, w, Z_gwas, snp_sd):
	Z_twas = snp_sd.dot(w * Z_gwas) / get_z_denom(V_cov, w)
	return Z_twas, get_pval(Z_twas)
	
def get_fusion_zscore(V_cov, w, Z_gwas, snp_sd=None):
	V_cor = get_V_cor(V_cov)
	Z_twas = np.vdot(Z_gwas, w) / get_z_denom(V_cor, w)
	return Z_twas, get_pval(Z_twas)

def get_burden_zscore(test_stat, get_zscore_args):
	if test_stat =='FUSION':
		return get_fusion_zscore(*get_zscore_args)
	if test_stat == 'SPrediXcan':
		return get_spred_zscore(*get_zscore_args)


def bgw_weight_file_info(w_path, weight_threshold=0, add_cols=[], drop_cols=['ID'], **kwargs):
	info_dict = tg.get_cols_dtype(
		tg.get_header(w_path, zipped=False), 
		cols=['CHROM','POS','REF','ALT','Trans','PCP','beta'], 
		add_cols=add_cols, 
		drop_cols=drop_cols, 
		get_id=True,  
		ind_namekey=True)
	return {'path': w_path, 
		'sampleID': [], 
		'target_ind': info_dict['file_cols'].index('Trans'), 
		'data_format': 'bgw_weight', 
		'genofile_type': 'bgw_weight', 
		'weight_threshold': weight_threshold,
		**info_dict}


def read_bgw_weight(w_path, weight_threshold=0, **kwargs):
	# ensure file at w_path exists
	if not os.path.isfile(w_path):
		print('No valid weight file for target.\n')
		raise tg.NoTargetDataError

	# get weight info
	w_info = bgw_weight_file_info(w_path)

	# initialize string
	proc_out = ''

	# read lines into proc_out
	with open(w_path, 'r') as w_file:
		for line in w_file:
			row = line.split('\t')
			line = '\t'.join([row[x] for x in w_info['col_inds']])
			line += '' if line.endswith('\n') else '\n'
			proc_out += line

	# read data into dataframe
	Weight = pd.read_csv(
		StringIO(proc_out),
		sep='\t',
		low_memory=False,
		header=None,
		comment='#',
		names=w_info['cols'],
		dtype=w_info['dtype'])

	Weight = tg.optimize_cols(Weight)

	# get snpIDs
	Weight['snpID'] = tg.get_snpIDs(Weight)
	Weight = Weight.drop_duplicates(['snpID'], keep='first').reset_index(drop=True)

	# get ES column
	Weight['ES'] = Weight['PCP'] * Weight['beta']
	Weight = Weight.drop(columns=['PCP','beta'])

	if weight_threshold:
		# filter out weights below threshold
		Weight = Weight[operator.gt(np.abs(Weight['ES']), weight_threshold)].reset_index(drop=True)

	if Weight.empty:
		print('No valid weight file for target.\n')
		raise tg.NoTargetDataError

	Weight = Weight.sort_values(by=['CHROM','POS']).reset_index(drop=True)

	return Weight


def read_bgw_zscore(num, w_df_info, snpIDs, **kwargs):
	chrm = w_df_info['CHROM'][num]
	start = w_df_info['min'][num]
	end = w_df_info['max'][num]

	z_path = args.z_path_pre + chrm + args.z_path_suf
	zscore_info = tg.zscore_file_info(z_path, chrm)

	Zscore = tg.read_tabix(start, end, **zscore_info)
	Zscore['snpIDflip'] = tg.get_snpIDs(Zscore, flip=True)

	snp_overlap = np.intersect1d(snpIDs, Zscore[['snpID','snpIDflip']])

	if not snp_overlap.size:
		return None

	# filter out non-matching snpID rows
	Zscore = Zscore[np.any(Zscore[['snpID','snpIDflip']].isin(snp_overlap), axis=1)].reset_index(drop=True)

	# if not in snpIDs, assumed flipped; if flipped, flip Zscore sign
	flip = np.where(Zscore.snpID.isin(snpIDs.values), 1, -1)

	if not np.all(flip == 1):
		Zscore['snpID'] = np.where(flip == 1, Zscore.snpID, Zscore.snpIDflip)
		Zscore['Zscore'] = flip * Zscore['Zscore']

	# drop unneeded columns
	Zscore = Zscore.drop(columns=['CHROM','POS','REF','ALT','snpIDflip'])

	return Zscore

def get_chrm_ld_data(chrm, snp_ids, ld_path_pre, ld_path_suf, **kwargs):
	ld_path = ld_path_pre + chrm + ld_path_suf
	MCOV = tg.get_ld_data(ld_path, snp_ids).reset_index()
	MCOV['CHROM'] = chrm
	return MCOV



#############################################################
# Print input arguments to log
# out_twas_path = args.out_dir + '/CHR' + args.chrm + '_sumstat_assoc.txt'

out_twas_path = args.out_dir + '/' + args.out_twas_file

print(
'''********************************
Input Arguments
Gene annotation file specifying genes for TWAS: {annot_path}
cis-eQTL weight file: {w_path_pre}[target]{w_path_suf}
GWAS summary statistics Z-score file: {z_path}
Reference LD genotype covariance file: {ld_path_pre}[CHRM]{ld_path_suf}
Gene training region SNP inclusion window: +-{window}
SNP weight inclusion threshold: {weight_threshold}
Test statistic to use: {test_stat_str}
Number of threads: {thread}
Output directory: {out_dir}
Output TWAS results file: {out_path}
********************************'''.format(
	**args.__dict__,
	test_stat_str = 'FUSION and SPrediXcan' if args.test_stat=='both' else args.test_stat,
	out_path = out_twas_path))

# tg.print_args(args)

###############################################################
### Read in gene annotation 
print('Reading gene annotation file.')
# Gene, TargetID, n_targets = tg.read_gene_annot_exp(args.annot_path, args.chrm)
Gene, TargetID, n_targets = tg.read_gene_annot_exp(**args.__dict__)

# read in headers for Weight and Zscore files; get the indices and dtypes for reading files into pandas
print('Reading file headers.\n')
# weight_info = tg.weight_file_info(**args.__dict__)
# zscore_info = tg.zscore_file_info(**args.__dict__)

# PREP OUTPUT - print output headers to files
print('Creating file: ' + out_twas_path + '\n')
out_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName','n_snps']
if args.test_stat == 'both':
	out_cols += ['FUSION_Z','FUSION_PVAL','SPred_Z','SPred_PVAL']
else:
	out_cols += ['Zscore','PVALUE']

pd.DataFrame(columns=out_cols).to_csv(
	out_twas_path,
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

	# read weight file
	w_path = args.w_path_pre + target + args.w_path_suf
	Weight = read_bgw_weight(w_path, **args.__dict__)

	# get start and end positions to tabix, per chromosome
	w_df_info = Weight.groupby('CHROM')['POS'].agg(['min','max']).reset_index().astype(str)

	# read in Zscore files
	Zscore = pd.concat([read_bgw_zscore(num, w_df_info, Weight.snpID) for num in range(w_df_info.shape[0])]).reset_index(drop=True)

	# merge Zscore and Weight dataframes on snpIDs
	ZW = Weight.merge(Zscore[['snpID','Zscore']], 
		left_on='snpID', 
		right_on='snpID', 
		how='inner')

	# snp_search_ids = ZW.snpID.values

	# Read in reference covariance matrix file by snpID
	MCOV = pd.concat([get_chrm_ld_data(str(chrm), ZW[ZW.CHROM == chrm].snpID.values, **args.__dict__) for chrm in np.unique(ZW.CHROM)]).reset_index(drop=True)


	

	# Read in reference covariance matrix file by snpID
	MCOV = tg.get_ld_data(args.ld_path, snp_search_ids)

	if MCOV.empty:
	  print('No reference covariance information for target SNPs for TargetID: ' + target + '\n')
	  return None

	# get the snp variance and covariance matrix
	snp_sd, V_cov = tg.get_ld_matrix(MCOV)

	ZW = ZW[ZW.snpID.isin(MCOV.snpID)]
	n_snps = str(ZW.snpID.size)

	print('Running TWAS.\nN SNPs=' + n_snps)

	### create output dataframe
	Result = Gene_Info.copy()
	Result['n_snps'] = n_snps

	### calculate zscore(s), pvalue(s)
	# get_zscore_args = [V_cov, ZW, snp_sd]
	get_zscore_args = [V_cov, ZW.ES.values, ZW.Zscore.values, snp_sd]

	if args.test_stat == 'both':
		Result['FUSION_Z'], Result['FUSION_PVAL'] = get_fusion_zscore(*get_zscore_args)

		Result['SPred_Z'], Result['SPred_PVAL'] = get_spred_zscore(*get_zscore_args)

	else:
		Result['TWAS_Zscore'], Result['PVALUE'] = get_burden_zscore(args.test_stat, get_zscore_args)

	# write to file
	Result.to_csv(
		out_twas_path,
		sep='\t',
		index=None,
		header=None,
		mode='a')

	print('Target TWAS completed.\n')

###############################################################
# thread process
if __name__ == '__main__':
	print('Starting TWAS for ' + str(n_targets) + ' target genes.\n')
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







