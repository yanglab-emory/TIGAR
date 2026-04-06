#!/usr/bin/env python

###################################################################
# Import packages needed
import argparse
import functools
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
parser = argparse.ArgumentParser(description='BGW-TWAS Association Study')

# Specify tool directory
parser.add_argument('--TIGAR_dir', type=str)

# Gene annotation file path
parser.add_argument('--gene_list', type=str, dest='annot_path')

# # chromosome number
# parser.add_argument('--chr', type=str, dest='chrm')

# use target dir (ie, weight files will be in [weight_prefix]/[target]/ instead of [weight_prefix]])
parser.add_argument('--target_dir', type=int, dest='target_dir', default=0)

# Weight file path
parser.add_argument('--weight_prefix', type=str, dest='w_path_pre')
parser.add_argument('--weight_suffix', type=str, dest='w_path_suf')

# GWAS Z score file path
parser.add_argument('--Zscore', type=str, dest='z_path', default='')

parser.add_argument('--Zscore_prefix', type=str, dest='z_path_pre', default='')
parser.add_argument('--Zscore_suffix', type=str, dest='z_path_suf', default='')


# Reference covariance file path
parser.add_argument('--LD_prefix', type=str, dest='ld_path_pre', default='')
parser.add_argument('--LD_suffix', type=str, dest='ld_path_suf', default='')
# , deprecated=True ## py 3.9 and up
parser.add_argument('--TIGAR_LD_prefix', type=str, dest='tigar_ld_path_pre', default='')
parser.add_argument('--TIGAR_LD_suffix', type=str, dest='tigar_ld_path_suf', default='')


# plink
parser.add_argument('--plink_LD', type=str, dest='plink_pre_chrm', default='')
# alias=['plink_LD_pre'] ## alias added later
parser.add_argument('--plink_LD_suf', type=str, dest='plink_suf', default='')
# parser.add_argument('--LD_suffix', type=str, dest='tigar_ld_path_suf')


# are plink files per-chromosome (need to merge?)
parser.add_argument('--plink_per_chr', type=bool, dest='do_plink_per_chrm', default=1)

# parser.add_argument('--LD', type=str, dest='ld_path')
# log file name
parser.add_argument('--log_file', type=str, default='')

# window
# parser.add_argument('--window',type=float)

# Filter to use only cis- or trans- snps
parser.add_argument('--snp_type', type=str, choices=['both','cis','trans'], default='both')

# Weight threshold to include SNP in TWAS
parser.add_argument('--weight_threshold', type=float, default=0)

# specify 'FUSION', 'SPrediXcan', or 'both': Zscore test statistic to use
parser.add_argument('--test_stat', type=str, default='')

# Number of threads
parser.add_argument('--thread', type=int, default=1)

# Output dir
parser.add_argument('--out_dir', type=str)

# output prefix
parser.add_argument('--out_prefix', type=str, default='')

# output file
parser.add_argument('--out_twas_file', type=str, default='')

# 
parser.add_argument('--job_suf', type=str)


parser.add_argument('--clean_output', type=int, choices=[0, 1], default=1, 
	help='clean_output (0: keep temp files, 1: remove temp files [default])')

# 'CHR' + args.chrm + '_BGWTWAS_assoc.txt'

args = parser.parse_args()

sys.path.append(args.TIGAR_dir)
sys.path.append(args.TIGAR_dir + '/TWAS')

###############################################################
# set output file names
if not args.out_prefix:
	args.out_prefix = 'BGWTWAS'

if not args.log_file:
	args.log_file = args.out_prefix + '_log.txt'

if not args.out_twas_file:
	args.out_twas_file = args.out_prefix + '_asso.txt'

# random temp directory names for plink files
if args.do_plink and (not args.job_suf):
	args.job_suf = str(int.from_bytes(os.urandom(2), byteorder='big'))

# directories
abs_out_dir = tg.get_abs_path(args.out_dir)
args.temp_out_dir = abs_out_dir + '/TWAS_temp_' + args.job_suf + '/'
out_twas_path = abs_out_dir + '/' + args.out_twas_file

# Make output, log directories
os.makedirs(args.out_dir, exist_ok=True)
os.makedirs(os.path.join(args.out_dir, 'logs'), exist_ok=True)
os.makedirs(args.temp_out_dir, exist_ok=True)

# set stdout to log
sys.stdout = open(os.path.join(args.out_dir, 'logs', args.log_file), 'w')

# Check tabix command
try:
	subprocess.check_call(['which','tabix'], stdout=subprocess.DEVNULL)
except:
	raise SystemExit('Error: Required tool TABIX is not available.\n')

# Check gene list file
try:
	os.path.isfile(args.annot_path)
except:
	SystemExit('Error: Gene list file does not exist.')

# handle using deprecated --LD_prefix, --LD_suffix
if (args.ld_path_pre != '') and (args.tigar_ld_path_pre == ''):
	args.tigar_ld_path_pre = args.ld_path_pre

if (args.ld_path_suf != '') and (args.tigar_ld_path_suf == ''):
	args.tigar_ld_path_suf = args.ld_path_suf

# make sure user actually input something for LD
# okay for args.plink_suf to be empty
if all([arg == '' for arg in [args.tigar_ld_path_pre, args.tigar_ld_path_suf, args.plink_pre_chrm]]):
	raise SystemExit('Error: No user input for LD. Must specify --plink_LD or  {--TIGAR_LD_prefix and --TIGAR_LD_suffix}. --plink_LD_suf is optional and may be empty if the file names end with the chromosome number (Ex: "PLINK_chr[chrm #].bim", etc.). \n')

# do plink?
args.do_plink = (args.tigar_ld_path_pre == '')

# if not do_plink_per_chrm, only one plink_pre
if args.do_plink and (not args.do_plink_per_chrm):
	args.plink_pre = args.plink_pre_chr

# set custom defaults for test_stat based on plink
if args.test_stat == '':
	args.test_stat = 'FUSION' if args.do_plink else 'both'

# if user input something for the test_stat
if (args.do_plink) and (args.test_stat != 'FUSION'):
		raise SystemExit('Error: User entered --test_stat is' + args.test_stat + '; "both" and "SPrediXcan" are not valid options when using PLINK for LD; the test_stat must be "FUSION". this is the default setting when PLINK LD is specified.\n')



###############################################################
## Import TIGAR functions, define other functions
import TIGARutils as tg

# from Asso_Study_02 import pval, get_V_cor, zscore_denom, spred_zscore, fusion_zscore, burden_zscore

def pval(z): return np.format_float_scientific(1-chi2.cdf(z**2, 1), precision=15, exp_digits=0)

def get_V_cor(V_cov):
	V_cov = V_cov.copy()
	v = np.sqrt(np.diag(V_cov))
	outer_v = np.outer(v, v)
	V_cor = V_cov / outer_v
	V_cor[V_cov == 0] = 0
	return V_cor

def zscore_denom(V, w):
	return np.sqrt(np.linalg.multi_dot([w, V, w]))

def spred_zscore(V_cov, w, Z_gwas, snp_sd):
	Z_twas = snp_sd.dot(w * Z_gwas) / zscore_denom(V_cov, w)
	return Z_twas, pval(Z_twas)
	
def fusion_zscore(V_cov, w, Z_gwas, **kwargs):
	V_cor = get_V_cor(V_cov)
	Z_twas = np.vdot(Z_gwas, w) / zscore_denom(V_cor, w)
	return Z_twas, pval(Z_twas)

def fusion_zscore_vcor(V_cor, w, Z_gwas, **kwargs):
	Z_twas = np.vdot(Z_gwas, w) / zscore_denom(V_cor, w)
	return Z_twas, pval(Z_twas)	

def burden_zscore(test_stat, get_zscore_args):
	if test_stat =='FUSION':
		return fusion_zscore(*get_zscore_args)
	if test_stat == 'SPrediXcan':
		return spred_zscore(*get_zscore_args)


def bgw_weight_file_info(w_path, weight_threshold=0, add_cols=[], drop_cols=['ID'], **kwargs):
	# cols=['CHROM','POS','REF','ALT','Trans','PCP','beta'],
	# got_header=tg.get_header(w_path, zipped=False, rename={'#CHR':'CHROM','CHR':'CHROM','CPP':'PCP'})
	# print(got_header)
	info_dict = tg.get_cols_dtype(
		tg.get_header(w_path, zipped=False, rename={'#CHR':'CHROM','CHR':'CHROM','CPP':'PCP','Beta':'BETA', 'beta':'BETA'}), 
		cols=['CHROM','POS','REF','ALT','Trans','PCP','BETA'], 
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


def read_bgw_weight(w_path, snp_type='both', weight_threshold=0, **kwargs):
	# ensure file at w_path exists
	if not os.path.isfile(w_path):
		print('No valid weight file for target at: ' + w_path + '\n')
		# print('No valid weight file for target.\n')
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
	Weight['ES'] = Weight['PCP'] * Weight['BETA']
	Weight = Weight.drop(columns=['PCP','BETA'])

	# filter cis- or trans- only
	if not snp_type == 'both':
		Weight = Weight[Weight.Trans == {'cis':0, 'trans':1}[snp_type]]

	# filter out weights below threshold
	if weight_threshold:
		Weight = Weight[operator.gt(np.abs(Weight['ES']), weight_threshold)].reset_index(drop=True)

	if Weight.empty:
		print('No valid weight file for target.\n')
		raise tg.NoTargetDataError

	Weight = Weight.sort_values(by=['CHROM','POS']).reset_index(drop=True)

	return Weight

def get_z_path(chrm):
	return args.z_path_pre + str(chrm) + args.z_path_suf

def read_bgw_zscore(z_path, snp_ids, **kwargs):

	# z_path = args.z_path_pre + chrm + args.z_path_suf
	zscore_info = tg.zscore_file_info(z_path, '')

	# Zscore = tg.read_tabix(start, end, **zscore_info)

	# get region strings
	regs_lst = list(tg.get_regions_list(snp_ids))
	N = len(regs_lst)

	# get function for filtering line
	filter_line = functools.partial(tg.filter_other_line, col_inds=zscore_info['col_inds'])

	# read in specific snps
	try:
		regs_str = ' '.join(regs_lst)
		proc_out = tg.call_tabix_regions(z_path, regs_str, filter_line = filter_line)

	except OSError:
		# argument may be too long for OS; if so try subset instead of getting all regions at once
		print('Subseting regions to tabix.')
		n = 2500
		while n:
			try: 
				regs_str_lst = [' '.join(regs_lst[i:i+n]) for i in range(0, N, n)]
				proc_out = b''.join([tg.call_tabix_regions(z_path, regs_str, filter_line = filter_line) for regs_str in regs_str_lst])
			except OSError:
				n -= 500
				pass
			else:
				n = 0

	# read data into dataframe
	Zscore = pd.read_csv(
		StringIO(proc_out.decode('utf-8')),
		sep='\t',
		low_memory=False,
		header=None,
		names=zscore_info['cols'],
		dtype=zscore_info['dtype'])

	Zscore = tg.optimize_cols(Zscore)

	# get snpIDs
	Zscore['snpID'] = tg.get_snpIDs(Zscore)
	Zscore = Zscore.drop_duplicates(['snpID'], keep='first').reset_index(drop=True)
	Zscore['snpIDflip'] = tg.get_snpIDs(Zscore, flip=True)

	return Zscore


def get_chrm_ld_data(chrm, snp_ids, tigar_ld_path_pre, tigar_ld_path_suf, **kwargs):
	ld_path = tigar_ld_path_pre + chrm + tigar_ld_path_suf
	MCOV = tg.get_ld_data(ld_path, snp_ids).reset_index()
	MCOV['CHROM'] = chrm
	return MCOV


def get_multi_chrm_ld_matrix(MCOV, return_diag=False):
	MCOV = MCOV.copy()
	
	MCOV['COV'] =  MCOV['COV'].apply(lambda x:np.fromstring(x, dtype=np.float32, sep=','))

	inds = MCOV.index
	n_inds = inds.size
	V_upper = np.zeros((n_inds, n_inds))

	for ii, idx_i in enumerate(inds):
		cov_i = MCOV.COV.at[idx_i]
		N = cov_i.size

		# index of last row with matching chrom
		last_chrm_ind = np.where(MCOV.CHROM == MCOV.CHROM.at[idx_i])[0][-1]

		# n_chrm_inds = chrm_jj.size
		possible_jj = slice(ii, last_chrm_ind + 1)
		possible_idx_j = inds[possible_jj]  # These are the snp_inds at j
		idxs_j_minus_i = possible_idx_j - idx_i # This is inds[j] - inds[i], as below
		# Get the *valid* indices, those less than N (0-indexed)
		valid_jj = np.where(idxs_j_minus_i < N)[0]

		V_upper[ii, valid_jj + ii] = cov_i[idxs_j_minus_i[valid_jj]]

	snp_Var = V_upper.diagonal()
	V = V_upper + V_upper.T - np.diag(snp_Var)
	snp_sd = np.sqrt(snp_Var)

	if return_diag:
		return snp_sd, V, snp_Var

	else:
		return snp_sd, V


def call_PLINK_merge_2chrms(target, mrg_chrms, temp_out_dir, clean_output=1, **kwargs):
	file_1_str = temp_out_dir + target + '_CHR' + str(mrg_chrms[0])
	file_2_str = ' '.join(temp_out_dir + target + '_CHR' + str(mrg_chrms[1]) + x for x in ['.bed', '.bim', '.fam'])

	cmd = ['plink --bfile ' + file_1_str + ' --bmerge ' + file_2_str + ' --make-bed --out ' + args.temp_out_dir + target]

	try:
		proc = subprocess.check_call(cmd, stdout=subprocess.PIPE, shell=True)

	except subprocess.CalledProcessError:
		print('Unable to merge LD files for TargetID: ' + target + ', CHRMS ' + str(mrg_chrms[0]) + ', ' + str(mrg_chrms[1]) + '\n')
		if clean_output:
			clean_plink_output(target, temp_out_dir)
		raise NoTargetDataError

	return None


#############################################################
# Print input arguments to log
# out_twas_path = args.out_dir + '/CHR' + args.chrm + '_sumstat_assoc.txt'


# cis-eQTL weight file: {w_path_pre}[target]{w_path_suf}

print(
'''********************************
Input Arguments
Gene annotation file specifying genes for TWAS: {annot_path}
cis-eQTL weight file: {in_w_dir}/[target]{w_path_suf}
GWAS summary statistics Z-score file: {in_z_path}
{ld_str}
SNP weight inclusion threshold: {weight_threshold}
Using {snp_type_str} SNPs for analysis
Test statistic to use: {test_stat_str}
Number of threads: {thread}
Output directory: {out_dir}
Output TWAS results file: {out_path}
Clean output: {clean_str}
********************************'''.format(
	**args.__dict__,
	snp_type_str = 'cis- and trans-' if args.test_stat=='both' else args.snp_type + '-',
	test_stat_str = 'FUSION and SPrediXcan' if args.test_stat=='both' else args.test_stat,
	ld_str='PLINK files for reference LD:' + args.plink_pre_chrm + '[CHRM]' + args.plink_suf if args.do_plink else 'Reference TIGAR LD genotype covariance files:' + args.tigar_ld_path_pre + '[chr]' + args.tigar_ld_path_suf,
	in_w_dir = args.target_dir + '/[target]' if args.target_dir else args.w_path_pre,
	in_z_path = args.z_path_pre + '[CHRM]' + args.z_path_suf if args.z_path == '' else args.z_path,
	clean_str = bool(args.clean_output),
	out_path = out_twas_path))


	# in_w_dir = args.w_path_pre if not args.target_dir else args.w_path_pre + '/[target]',
tg.print_args(args)

###############################################################
### Read in gene annotation 
print('Reading gene annotation file.')
Gene = pd.read_csv(
	args.annot_path,
	sep='\t',
	header=None)[0].drop_duplicates()
n_targets = Gene.size

# read in headers for Weight and Zscore files; get the indices and dtypes for reading files into pandas

# PREP OUTPUT - print output headers to files
print('Creating file: ' + out_twas_path + '\n')
out_cols = ['GeneName','n_snps','n_trans_snps']
if args.test_stat == 'both':
	out_cols += ['FUSION_Z','FUSION_PVAL','SPred_Z','SPred_PVAL']
else:
	# out_cols += ['Zscore','PVALUE']
	# out_cols += [{'FUSION':'FUSION', 'SPrediXcan':'SPred'}[args.test_stat] + x for x in ['_Z', '_PVAL']]
	test_stat_str = 'SPred' if args.test_stat == 'SPrediXcan' else 'FUSION'
	out_cols += [test_stat_str + x for x in ['_Z', '_PVAL']]

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
	target = Gene[num]
	print('num=' + str(num) + '\nTargetID=' + target)

	# read weight file
	if not args.target_dir:
		w_path = args.w_path_pre + '/' + target + args.w_path_suf
	else:
		w_path = args.w_path_pre + '/' + target + '/' + target + args.w_path_suf

	Weight = read_bgw_weight(w_path, **args.__dict__)

	# get start and end positions to tabix, per chromosome
	# w_df_info = Weight.groupby('CHROM')['POS'].agg(['min','max']).reset_index().astype(str)

	# read in Zscore files
	if args.z_path is not '':
		Zscore = read_bgw_zscore(args.z_path, Weight.snpID.values)
		# print(Weight.snpID.values)
		# print(Zscore)
	else:
		Zscore = pd.concat([read_bgw_zscore(get_z_path(chrm), Weight[Weight.CHROM == chrm].snpID.values) for chrm in np.unique(Weight.CHROM)]).reset_index(drop=True)
		# print(Zscore)
		# print(get_z_path(chrm))

	# merge Weight, Zscore file on SNPs with Weight as reference
	ZW, snp_overlap = tg.merge_ref_z_on_snps(Weight, Zscore, target)
	ZW = ZW.drop_duplicates(['snpID'], keep='first').reset_index(drop=True)

	# snp_overlap = np.intersect1d(Weight.snpID, Zscore[['snpID','snpIDflip']])
	# # print(Weight)
	# # print(Zscore)
	# if not snp_overlap.size:
	# 	print('No overlapping test SNPs that have magnitude of cis-eQTL weights greater than threshold value and with GWAS Zscore for TargetID: ' + target + '\n')
	# 	return None

	# # filter out non-matching snpID rows
	# Weight = Weight[Weight.snpID.isin(snp_overlap)]
	# Zscore = Zscore[np.any(Zscore[['snpID','snpIDflip']].isin(snp_overlap), axis=1)].reset_index(drop=True)

	# # if not in Weight.snpIDs, assumed flipped; if flipped, flip Zscore sign
	# flip = np.where(Zscore.snpID.isin(Weight.snpID), 1, -1)

	# if not np.all(flip == 1):
	# 	Zscore['snpID'] = np.where(flip == 1, Zscore.snpID, Zscore.snpIDflip)
	# 	Zscore['Zscore'] = flip * Zscore['Zscore']

	# # drop unneeded columns
	# Zscore = Zscore.drop(columns=['CHROM','POS','REF','ALT','snpIDflip'])

	# # merge Zscore and Weight dataframes on snpIDs
	# ZW = Weight.merge(Zscore[['snpID','Zscore']], 
	# 	left_on='snpID', 
	# 	right_on='snpID', 
	# 	how='inner')

	if args.do_plink:
		try: 
			# if separate files for each chrm
			if args.do_plink_per_chrm:

				mrg_chrms = []

				# loop through unique chromosomes
				print('call_PLINK_extract')
				for chrm in np.unique(ZW.CHROM):
					try:
						# get start, end positions and plink file path for chrm
						start_pos = ZW[ZW.CHROM == chrm].POS.iloc[0]
						end_pos = ZW[ZW.CHROM == chrm].POS.iloc[-1]
						plink_pre = args.plink_pre_chrm + str(chrm) + args.plink_suf

						# create file
						plink_extract_target = tg.call_PLINK_extract(start_pos, end_pos, target, chrm, plink_pre, set_var_ids=True, **args.__dict__)
						
						mrg_chrms += [chrm]

					# don't stop if just one chromosome has no data
					except tg.NoTargetDataError:
						pass

				n_merge = len(mrg_chrms)

				# no ref data for any chrm
				if (not n_merge):
					print('No reference covariance information for target SNPs for TargetID: ' + target + '.\n')
					return None

				# only one chromosome, no need for merge
				elif (n_merge == 1):
					old_pre = args.temp_out_dir + target + '_CHR' + str(mrg_chrms[0])
					new_pre = args.temp_out_dir + target
					os.rename(old_pre + '.bed', new_pre + '.bed')
					os.rename(old_pre + '.bim', new_pre + '.bim')
					os.rename(old_pre + '.fam', new_pre + '.fam')

				# two chromosome, different plink command
				elif (n_merge == 2):
					plink_merge_targets = call_PLINK_merge_2chrms(target, mrg_chrms, **args.__dict__)

				# merge multiple 
				else:
					plink_merge_targets = call_PLINK_merge_mchrms(target, mrg_chrms, **args.__dict__)

			# only one plink LD file for all chrms
			else:
				target_ranges = []
				for chrm in np.unique(ZW.CHROM):
					start_pos = ZW[ZW.CHROM == chrm].POS.iloc[0]
					end_pos = ZW[ZW.CHROM == chrm].POS.iloc[-1]
					target_ranges += [' '.join([str(x) for x in [chrm, start_pos, end_pos]])]
				print('call_PLINK_extract_ranges')
				plink_extract_target = tg.call_PLINK_extract_ranges(target_ranges, target, **args.__dict__)

			## clean up the PLINK output
			print('read_format_ref_bim')
			target_ref = tg.read_format_ref_bim(target, **args.__dict__)

			# filter again, with target_ref this time
			print('merge on snps')
			ZW, snp_overlap_ld = tg.merge_ref_z_on_snps(target_ref, ZW, target)

			# get reference covariance matrix
			print('VCOR')
			V_cor = tg.get_plink_LD_matrix(target, snp_overlap_ld, **args.__dict__)

			# get_zscore_args = [V_cov, ZW.ES.values, ZW.Zscore.values, snp_sd]
			zscore_args = {'V_cor': V_cor, 'w': ZW.ES.values, 'Z_gwas': ZW.Zscore.values}

		except Exception as e:
			raise e

		finally:
			tg.clean_plink_output(target, **args.__dict__)

	# TIGAR LD files
	else:

		# Read in reference covariance matrix file by snpID
		MCOV = pd.concat([get_chrm_ld_data(str(chrm), ZW[ZW.CHROM == chrm].snpID.values, **args.__dict__) for chrm in np.unique(ZW.CHROM)]).reset_index(drop=True)

		if MCOV.empty:
			print('No reference covariance information for target SNPs for TargetID: ' + target + '\n')
			return None

		# set index
		MCOV.index = MCOV.row
		MCOV = MCOV.drop(columns='row')

		# get the snp variance and covariance matrix
		snp_sd, V_cov = get_multi_chrm_ld_matrix(MCOV)

		ZW = ZW[ZW.snpID.isin(MCOV.snpID)]
		ZW = ZW.drop_duplicates(['snpID'], keep='first').reset_index(drop=True)

		### calculate zscore(s), pvalue(s)
		ztest_args = {'V_cov': V_cov, 'w': ZW.ES.values, 'Z_gwas': ZW.Zscore.values, 'snp_sd': snp_sd}


	# number of snps
	n_snps = str(ZW.snpID.size)
	n_trans_snps = str(np.sum(ZW.Trans))

	print('Running TWAS.\nN SNPs=' + n_snps)

	### create output dataframe
	Result = pd.DataFrame.from_records({'GeneName': target}, index=[0])
	Result['n_snps'] = n_snps
	Result['n_trans_snps'] = n_trans_snps

	if args.do_plink:
		# PLINK LD: FUSION Zscore from V_cor
		Result['FUSION_Z'], Result['FUSION_PVAL'] = fusion_zscore_vcor(**zscore_args)
	else:
		# TIGAR LD
		if args.test_stat == 'both':
			# FUSION Zscore from V_cov
			Result['FUSION_Z'], Result['FUSION_PVAL'] = fusion_zscore(**ztest_args)
			Result['SPred_Z'], Result['SPred_PVAL'] = spred_zscore(**ztest_args)
		else:
			# FUSION Zscore from V_cov
			Result['TWAS_Zscore'], Result['PVALUE'] = burden_zscore(args.test_stat, ztest_args)

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







