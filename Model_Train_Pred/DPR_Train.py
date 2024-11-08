#!/usr/bin/env python

#############################################################
# Import packages needed
import argparse
import multiprocessing
import operator
import os
import subprocess
import sys

from time import time

import numpy as np
import pandas as pd

import scipy.stats as stats
from sklearn.model_selection import KFold
import statsmodels.api as sm

#############################################################
# time calculation
start_time = time()

#############################################################
# parse input arguments
parser = argparse.ArgumentParser(description='DPR Training')

parser.add_argument('--chr', type=str, dest='chrm', 
	choices=[str(i + 1) for i in range(22)],
	required=True, 
	help='chromosome number')
parser.add_argument('--cvR2', type=int, choices=[0, 1], default=1, 
	help='cvR2 (0: no CV, 1: CV)')
parser.add_argument('--cvR2_threshold', type=float, default=0.005, 
	help='threshold cvR2 value for training (default: 0.005)')
parser.add_argument('--dpr', type=str, choices=['1', '2'], default=1, 
	help='Bayesian inference algorithm used by DPR (1: Variational Bayesian [default], 2: MCMC)')
parser.add_argument('--ES', type=str, choices=['additive', 'fixed'], default='fixed', 
	help='output effect size (additive: fixed + random, fixed: fixed effects [default])')
parser.add_argument('--format', type=str, dest='data_format', choices=['DS','GT'], default='GT')
parser.add_argument('--gene_exp', type=str, dest='geneexp_path', required=True)
parser.add_argument('--genofile', type=str, dest='geno_path', required=True)
parser.add_argument('--genofile_type', type=str, choices=['vcf', 'dosage'], 
	help='filetype of genofile (vcf, dosages)')
parser.add_argument('--hwe', type=float, default=0.00001, 
	help='threshold p-value for Hardy Weinberg Equilibrium exact test (default: 0.00001)')
parser.add_argument('--job_suf', type=str, default='', help=
	'suffix added to temp directories for DPR intermmediate files (default: out_prefix value)')
parser.add_argument('--log_file', type=str, default='')
parser.add_argument('--maf', type=float, default=0.01, 
	help='folded Minor Allele Frequency threshold; range from 0-0.5 (default: 0.01)')
parser.add_argument('--missing_rate', type=float, default=0.2, 
	help='missing rate threshold for excluding SNPs with too many missing values (default: 0.2)')
parser.add_argument('--out_dir', type=str, default=os.getcwd())
parser.add_argument('--out_info_file', type=str, default='')
parser.add_argument('--out_prefix', type=str, default='')
parser.add_argument('--out_weight_file', type=str, default='')
parser.add_argument('--sub_dir', type=int, choices=[0, 1], default=1)
parser.add_argument('--thread', type=int, default=1)
parser.add_argument('--TIGAR_dir', type=str, help='tool directory', 
	default=os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
parser.add_argument('--train_sampleID', type=str, dest='sampleid_path', required=True)
parser.add_argument('--window', type=int, default=1000000, 
	help='size around gene region for selecting cis-SNPs for fitting gene expression prediction model (default: 1000000 [ie, +-1MB region around gene])')

args = parser.parse_args()
sys.path.append(args.TIGAR_dir)
import TIGARutils as tg

#############################################################
# set output file names
if not args.out_prefix:
	args.out_prefix = 'CHR' + args.chrm + '_DPR_train'

if not args.log_file:
	args.log_file = args.out_prefix + '_log.txt'

if not args.out_info_file:
	args.out_info_file = args.out_prefix + '_GeneInfo.txt'

if not args.out_weight_file:
	args.out_weight_file = args.out_prefix + '_eQTLweights.txt'

if not args.job_suf:
	args.job_suf = args.out_prefix

# sub-directory in out directory
if args.sub_dir:
	out_sub_dir = os.path.join(args.out_dir, 'DPR_CHR' + args.chrm)
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
## Store files for DPR training
os.makedirs(os.path.join(out_sub_dir, 'DPR_Files' + args.job_suf), exist_ok=True)
## Store Cross Validation DPR files
if args.cvR2:
	os.makedirs(os.path.join(out_sub_dir, 'CV_Files' + args.job_suf), exist_ok=True)

# set stdout to log
sys.stdout = open(os.path.join(args.out_dir, 'logs', args.log_file), 'w')

#############################################################
# DEFINE, IMPORT FUNCTIONS
# import TIGARutils as tg

# directory for temporary output files; need to be defined here for some functions
dpr_file_dir = out_sub_dir + '/DPR_Files' + args.job_suf + '/'
dpr_file_dir_cv = out_sub_dir + '/CV_Files' + args.job_suf + '/'

# preps dpr input files, runs DPR, reads in dpr output
def prep_call_dpr(Bimbam_df, Pheno_df, dpr_file_dir, targetid):
	## PATHS FOR DPR INPUT
	bimbam_pth = dpr_file_dir + targetid + '_bimbam.txt'
	pheno_pth = dpr_file_dir + targetid + '_pheno.txt'

	#  OUTPUT FILES FOR DPR INPUT
	out_args = {'header':False, 'index':None, 'sep':'\t', 'mode':'w', 'float_format':'%f'}
	Bimbam_df.to_csv(bimbam_pth, **out_args)
	Pheno_df.to_csv(pheno_pth, **out_args)

	# CALL DPR
	try:
		DPR_call_args = [DPR_path, 
			'-g', bimbam_pth, 
			'-p', pheno_pth, 
			'-dpr', args.dpr, 
			'-notsnp',
			'-o', 'DPR_' + targetid]

		subprocess.check_call(
			DPR_call_args,
			cwd=dpr_file_dir,
			stdout=subprocess.DEVNULL,
			stderr=subprocess.DEVNULL)

	except subprocess.CalledProcessError as err:
		raise err

	finally:
		os.remove(bimbam_pth)
		os.remove(pheno_pth)

	# READ IN AND PROCESS DPR OUTPUT
	dpr_out_pth = dpr_file_dir + 'output/DPR_' + targetid + '.param.txt'

	DPR_Out = pd.read_csv(
		dpr_out_pth,
		sep='\t',
		header=0,
		names=['CHROM','snpID','POS','n_miss','b','beta','gamma'],
		usecols=['snpID','b','beta'],
		dtype={'snpID':object, 'b':np.float64, 'beta':np.float64})

	os.remove(dpr_out_pth)

	DPR_Out = tg.optimize_cols(DPR_Out)

	# GET EFFECT SIZE
	if args.ES == 'fixed':
		DPR_Out['ES'] = DPR_Out['beta']

	elif args.ES == 'additive':
		DPR_Out['ES'] = DPR_Out['beta'] + DPR_Out['b']

	return DPR_Out

# calculated r2 of prediction based on ES weights, centered geno, centered pheno
def calc_r2(Weight_df, Geno_df, Pheno_df, cv=False):

	if cv:
		Weight_df = Weight_df.fillna(0)

	### calculate predicted value for test set
	Pred = np.dot(Geno_df, Weight_df)
	lm = sm.OLS(Pheno_df.values, sm.add_constant(Pred)).fit()

	if cv:
		return lm.rsquared

	# else, return Pvalue, R2 for final training
	return lm.f_pvalue, lm.rsquared


# function to do the ith cross validation step
# def do_cv(i, target, Bimbam_df, Pheno_df, cv_trainID, cv_testID):
def do_cv(target, Bimbam_df, Pheno_df, BimbamC_df, PhenoC_df, cv_train_test_id, k_folds=5):

	k_fold_R2 = []

	for i in range(k_folds):
		target_cv = target + '_CV' + str(i+1)

		trainID, testID = cv_train_test_id[i]

		# PREP INPUT, CALL DPR
		try:
			DPR_Out = prep_call_dpr(
				Bimbam_df[['snpID','REF','ALT', *trainID]], 
				Pheno_df.loc[trainID], 
				dpr_file_dir_cv, 
				target_cv)[['snpID','ES']]

			### for R2 calculation
			cv_rsquared = calc_r2(
				DPR_Out.ES,
				BimbamC_df[BimbamC_df.snpID.isin(DPR_Out.snpID)][testID].T,
				PhenoC_df.loc[testID], 
				cv=True)

			k_fold_R2.append(cv_rsquared)

		except subprocess.CalledProcessError as err:
			print('DPR failed in CV' + str(i+1))
			k_fold_R2.append(0)

	avg_r2_cv = sum(k_fold_R2) / k_folds
		
	# RETURN R2 RESULT
	return avg_r2_cv

def remove_dir(path):
	try:
		call_args = ['rm', '-rf', tg.get_abs_path(path)]
		subprocess.check_call(
			call_args,
			cwd=tg.get_abs_path(os.path.join(path, '..')),
			stdout=subprocess.DEVNULL)
	except subprocess.CalledProcessError as err:
		raise err


#############################################################
# set absolute paths
DPR_path = tg.get_abs_path(args.TIGAR_dir) + '/Model_Train_Pred/DPR'
# DPR_path = TIGAR_dir + '/Model_Train_Pred/DPR'

tmp_weight_path = out_sub_dir + '/temp_' + args.out_weight_file
out_weight_path = out_sub_dir + '/' + args.out_weight_file
out_info_path = out_sub_dir + '/' +  args.out_info_file

#############################################################
# Print input arguments to log
print(
'''********************************
Input Arguments
Gene Annotation and Expression file: {geneexp_path}
Training sampleID file: {sampleid_path}
Chromosome: {chrm}
Training genotype file: {geno_path}
Genotype file used for training is type: {genofile_type}
Genotype data format: {data_format}
Gene training region SNP inclusion window: +-{window}
Excluding SNPs if missing rate exceeds: {missing_rate}
MAF threshold for SNP inclusion: {maf}
HWE p-value threshold for SNP inclusion: {hwe}
{cvR2_str1} DPR model by 5-fold cross validation{cvR2_str2}.
DPR model: {dpr} - {dpr_type}
Output Effect-size type: {ES}
Number of threads: {thread}
Output directory: {out_dir}
Output training weights file: {out_weight}
Output training info file: {out_info}
********************************'''.format(
	**args.__dict__,
	dpr_type = {'1':'Variational Bayesian', '2':'MCMC'}[args.dpr],
	cvR2_str1 = {0:'Skipping evaluation of', 1:'Evaluating'}[args.cvR2],
	cvR2_str2 = {0:'', 1:' with inclusion threshold Avg.CVR2 >' + str(args.cvR2_threshold)}[args.cvR2],
	out_weight = out_weight_path,
	out_info = out_info_path))

# tg.print_args(args)

#############################################################
# Prepare DPR input

### Read in Gene Expression/Annotation file
### First five columns should be fixed:
### 1.CHROM
### 2.GeneStart
### 3.GeneEnd
### 4.TargetID [i.e.GeneID, treated as unique annotation for each gene]
### 5.GeneName

# Startup for training jobs: get column header info, sampleIDs
print('Reading genotype, expression file headers, sample IDs.\n')
sampleID, sample_size, geno_info, exp_info = tg.sampleid_startup(**args.__dict__)

# Read in gene expression info
print('Reading gene expression data.\n')
GeneExp, TargetID, n_targets = tg.read_gene_annot_exp(**exp_info)

# PREP CROSS VALIDATION SAMPLES - Split sampleIDs for cross validation
if args.cvR2:
	print('Splitting sample IDs randomly for 5-fold cross validation by average R2...\n')
	CV_train_test_ID = [(sampleID[x], sampleID[y]) for x,y in KFold(n_splits=5).split(sampleID)]
else:
	print('Skipping splitting samples for 5-fold cross validation...\n')

# PREP OUTPUT - print output headers to files
print('Creating file: ' + tmp_weight_path + '\n')
weight_cols = ['CHROM','POS', 'snpID', 'REF','ALT','TargetID','MAF','p_HWE','ES','b','beta']
pd.DataFrame(columns=weight_cols).to_csv(
	tmp_weight_path,
	sep='\t',
	header=True,
	index=None,
	mode='w')

print('Creating file: ' + out_info_path + '\n')
info_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName','sample_size','n_snp', 'n_effect_snp','CVR2','TrainPVALUE','TrainR2','CVR2_threshold']
pd.DataFrame(columns=info_cols).to_csv(
	out_info_path,
	sep='\t',
	header=True,
	index=None,
	mode='w')

print('********************************\n')

##############################################################
# thread function
@tg.error_handler
def thread_process(num):
	target = TargetID[num]
	print('num=' + str(num) + '\nTargetID=' + target)
	Expr = GeneExp.iloc[[num]]

	start = str(max(int(Expr.GeneStart) - args.window, 0))
	end = str(int(Expr.GeneEnd) + args.window)

	# READ IN AND PROCESS GENOTYPE DATA 
	# file must be bgzipped and tabix
	Geno = tg.read_tabix(start, end, sampleID, **geno_info)

	# filter out variants that exceed missing rate threshold
	Geno = tg.handle_missing(Geno, sampleID, args.missing_rate)
	
	# get, filter maf
	Geno = tg.calc_maf(Geno, sampleID, args.maf)

	# get, filter p_HWE
	Geno = tg.calc_p_hwe(Geno, sampleID, args.hwe)

	if Geno.empty:
		print('No valid data for target.')
		return None

	n_snp = Geno.snpID.size

	# Geno, Expr not centered since 
	## 1) DPR script centers input 
	## 2) DPR script sometimes segfaults when reading in genotype files when data was centered

	# PREP INPUT FILES, CALL DPR, READ IN DPR OUTPUT
	Bimbam = Geno[['snpID','REF','ALT', *sampleID]]
	Pheno = Expr[sampleID].T

	# center data for R2 calculation
	BimbamC = tg.center(Bimbam, sampleID).drop(columns=['REF', 'ALT'])
	PhenoC = tg.center(Expr[sampleID]).T

	# 5-FOLD CROSS-VALIDATION
	if args.cvR2:
		print('Running 5-fold CV.')

		avg_r2_cv = do_cv(target, Bimbam, Pheno, BimbamC, PhenoC, CV_train_test_ID)
		print('Average R2 for 5-fold CV: {:.4f}'.format(avg_r2_cv))

		if avg_r2_cv < args.cvR2_threshold:
			print('Average R2 < ' + str(args.cvR2_threshold) + '; Skipping DPR training for TargetID: ' + target + '\n')
			return None

	else:
		print('Skipping evaluation by 5-fold CV average R2...')
		avg_r2_cv = 0

	# FINAL MODEL TRAINING
	print('Running DPR training.')

	try:
		DPR_Out = prep_call_dpr(Bimbam, Pheno, dpr_file_dir, target)

	except subprocess.CalledProcessError as err:
		print('DPR failed for target.\n')
		return None

	# FILTER FOR SNPS WITH ES!=0
	DPR_Out = DPR_Out[DPR_Out.ES != 0]
	n_effect_snp = DPR_Out.ES.size

	# R2 CALCULATION
	BimbamC = BimbamC[BimbamC.snpID.isin(DPR_Out.snpID)][sampleID].T

	train_pvalue, train_rsquared = calc_r2(DPR_Out.ES, BimbamC, PhenoC)

	# OUTPUT TARGET WEIGHTS TO FILE
	# initialize df with MAF, pHWE, other info
	Weight = Geno[['CHROM','POS','REF','ALT','snpID','p_HWE','MAF']].copy()
	Weight['TargetID'] = target

	# merge with dpr output weights, reorder columns using existing col list
	Weight = Weight.merge(
		DPR_Out, 
		left_on='snpID',
		right_on='snpID',
		how='inner')[weight_cols]

	Weight.to_csv(
		tmp_weight_path,
		sep='\t',
		header=None,
		index=None,
		mode='a')

	# OUTPUT TARGET TRAINING INFO TO FILE
	# initialize dataframe for storing training info
	Info = Expr[
		['CHROM','GeneStart','GeneEnd','TargetID','GeneName']].copy()
	Info['sample_size'] = sample_size
	Info['n_snp'] = n_snp
	Info['n_effect_snp'] = n_effect_snp
	Info['CVR2'] = avg_r2_cv
	Info['TrainPVALUE'] = train_pvalue
	Info['TrainR2'] = train_rsquared
	Info['CVR2_threshold'] = args.cvR2_threshold if args.cvR2 else 0

	# output training info
	Info.to_csv(
		out_info_path,
		sep='\t',
		header=None,
		index=None,
		mode='a')

	print('Target DPR training completed.\n')

##############################################################
# start thread  process

# if (args.thread < int(len(EXP)/100) | args.thread > len(EXP)):
	# args.thread = (int(len(EXP)/100)+1)*100

if __name__ == '__main__':
	try:
		print('Starting DPR training for ' + str(n_targets) + ' target genes.\n')
		pool = multiprocessing.Pool(args.thread)
		pool.imap(thread_process,[num for num in range(n_targets)])
		pool.close()
		pool.join()
		print('Done.')

		# sort tabix output
		tg.sort_tabix_output(tmp_weight_path, out_weight_path)

	finally:
		# remove DPR file directories
		remove_dir(os.path.join(out_sub_dir, 'DPR_Files' + args.job_suf))

		if args.cvR2:
			remove_dir(os.path.join(out_sub_dir, 'CV_Files' + args.job_suf))		


############################################################
# time calculation
elapsed_sec = time() - start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)

sys.stdout.close()
