#!/usr/bin/env python

#############################################################
# Import packages needed
import argparse
import multiprocessing
import operator
import os
import subprocess
import sys

from io import StringIO
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

# Specify tool directory
parser.add_argument('--TIGAR_dir', type=str)

# for Gene annotation and Expression-level file path
parser.add_argument('--gene_exp', type=str, dest='geneexp_path')

# training genotype file sampleIDs path
parser.add_argument('--train_sampleID', type=str, dest='sampleid_path')

# specified chromosome number
parser.add_argument('--chr', type=str, dest='chrm')

# Genotype file path
parser.add_argument('--genofile', type=str, dest='geno_path')

# specified input file type (vcf or dosages)
parser.add_argument('--genofile_type', type=str)

# format of genotype data 'DS' or 'GT'
parser.add_argument('--format', type=str, dest='data_format')

# window
parser.add_argument('--window', type=int)

# missing rate: threshold for excluding SNPs with too many missing values
parser.add_argument('--missing_rate', type=float)

# maf
parser.add_argument('--maf', type=float)

# p-value for HW test
parser.add_argument('--hwe', type=float)

# cvR2
## 0 do not run cvR2
## 1 run cvR2 [default]
parser.add_argument('--cvR2', type=int)

# threshold cvR2 value for training (default: 0.005)
parser.add_argument('--cvR2_threshold', type=float)

# Bayesian inference algorithm used by DPR: 
## '1' (Variational Bayesian)
## '2' (MCMC)
parser.add_argument('--dpr', type=str)

# output effect-size
## 'fixed' (fixed effects) [default]
## 'additive' (fixed + random)
parser.add_argument('--ES', type=str)

# file paths
parser.add_argument('--out_weight_file', type=str)
parser.add_argument('--out_info_file', type=str)

# suffix to directories for DPR intermmediate files
parser.add_argument('--job_suf', type=str)

# threads to use
parser.add_argument('--thread', type=int)

# output dir
parser.add_argument('--out_dir', type=str)

args = parser.parse_args()

sys.path.append(args.TIGAR_dir)

#############################################################
# DEFINE, IMPORT FUNCTIONS
import TIGARutils as tg

# directory for temporary output files; need to be defined here for some functions
abs_out_dir = tg.get_abs_path(args.out_dir)
dpr_file_dir = abs_out_dir + '/DPR_Files' + args.job_suf + '/'
dpr_file_dir_cv = abs_out_dir + '/CV_Files' + args.job_suf + '/'

# preps dpr input files, runs DPR, reads in dpr output
def prep_call_dpr(bimbam_df, pheno_df, snpannot_df, dpr_file_dir, targetid):
	## PATHS FOR DPR INPUT
	bimbam_pth = dpr_file_dir + targetid + '_bimbam.txt'
	pheno_pth = dpr_file_dir + targetid + '_pheno.txt'
	snpannot_pth = dpr_file_dir + targetid + '_snp_annot.txt'

	#  OUTPUT FILES FOR DPR INPUT
	out_args = {'header':False, 'index':None, 'sep':'\t', 'mode':'w', 'float_format':'%f'}
	bimbam_df.to_csv(bimbam_pth, **out_args)
	pheno_df.to_csv(pheno_pth, **out_args)
	snpannot_df.to_csv(snpannot_pth, **out_args)

	# CALL DPR
	try:
		DPR_call_args = [DPR_path, 
			'-g', bimbam_pth, 
			'-p', pheno_pth, 
			'-a', snpannot_pth, 
			'-dpr', args.dpr, 
			'-notsnp',
			'-o', 'DPR_' + targetid]

		subprocess.check_call(
			DPR_call_args,
			cwd=dpr_file_dir,
			stdout=subprocess.DEVNULL)

	except subprocess.CalledProcessError as err:
		raise err

	finally:
		os.remove(bimbam_pth)
		os.remove(pheno_pth)
		os.remove(snpannot_pth)

	# READ IN AND PROCESS DPR OUTPUT
	dpr_out_pth = dpr_file_dir + 'output/DPR_' + targetid + '.param.txt'

	dpr_out = pd.read_csv(
		dpr_out_pth,
		sep='\t',
		header=0,
		names=['CHROM','snpID','POS','n_miss','b','beta','gamma'],
		usecols=['snpID','b','beta'],
		dtype={'snpID':object, 'b':np.float64, 'beta':np.float64})

	os.remove(dpr_out_pth)

	dpr_out = tg.optimize_cols(dpr_out)

	# GET EFFECT SIZE
	if args.ES == 'fixed':
		dpr_out['ES'] = dpr_out['beta']

	elif args.ES == 'additive':
		dpr_out['ES'] = dpr_out['beta'] + dpr_out['b']

	return dpr_out


# calculated r2 of prediction based on out_weights_df, genotype data in bimbam_test_df , actual values pheno_test_df
def calc_r2(out_weights_df, bimbam_test_df, pheno_test_df, cv=False):

	# filter by snp overlap
	snp_overlap = np.intersect1d(out_weights_df.snpID, bimbam_test_df.snpID)
	out_weights_df = out_weights_df[out_weights_df.snpID.isin(snp_overlap)]
	bimbam_test_df = bimbam_test_df[bimbam_test_df.snpID.isin(snp_overlap)]

	# genotype, weight data for prediction
	test_geno_weights = out_weights_df.merge(
		bimbam_test_df,
		left_on='snpID',
		right_on='snpID',
		how='outer').set_index(['snpID'])
	test_geno_weights = tg.optimize_cols(test_geno_weights)

	snp_weights = test_geno_weights['ES']

	if cv:
		snp_weights = snp_weights.fillna(0)

	test_geno = test_geno_weights.drop(columns=['ES']).T

	### calculate predicted value for test set
	target_pred = np.dot(test_geno, snp_weights)

	lm = sm.OLS(pheno_test_df.values, sm.add_constant(target_pred)).fit()

	if cv:
		return lm.rsquared

	# else, return Pvalue, R2 for final training
	return lm.f_pvalue, lm.rsquared


# function to do the ith cross validation step
def do_cv(i, target, Geno_df, Expr_df, snp_annot_df, cv_trainID, cv_testID):
	target_cv = target + '_CV' + str(i+1)

	trainID = cv_trainID[i]
	testID = cv_testID[i]

	bimbam_train = Geno_df[['snpID','REF','ALT', *trainID]]

	pheno_train = Expr_df[trainID].T

	# PREP INPUT, CALL DPR
	try:
		dpr_out_cv = prep_call_dpr(
			bimbam_train, 
			pheno_train, 
			snp_annot_df,
			dpr_file_dir_cv, 
			target_cv)

	except subprocess.CalledProcessError as err:
		print('DPR failed in CV' + str(i+1) + ' for TargetID: ' + target)
		return 0
	
	### for R2 calculation
	out_weights_cv = dpr_out_cv[['snpID','ES']]
	bimbam_test = Geno_df[['snpID', *testID]]
	pheno_test = Expr_df[testID].T

	cv_rsquared = calc_r2(
		out_weights_cv, 
		bimbam_test, 
		pheno_test, 
		cv=True)

	# RETURN R2 RESULT
	return(cv_rsquared)

#############################################################
# set absolute paths
DPR_path = tg.get_abs_path(args.TIGAR_dir) + '/Model_Train_Pred/DPR'

# check input arguments
if args.genofile_type == 'vcf':
	if (args.data_format != 'GT') and (args.data_format != 'DS'):
		raise SystemExit('Please specify the genotype data format used by the vcf file (--format ) as either "GT" or "DS".\n')

elif args.genofile_type == 'dosage':
	args.data_format = 'DS'

else:
	raise SystemExit('Please specify the type input genotype file type (--genofile_type) as either "vcf" or "dosage".\n')

out_weight_path = args.out_dir + '/temp_' + args.out_weight_file
out_info_path = args.out_dir + '/' +  args.out_info_file

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
	kf = KFold(n_splits=5)
	kf_splits = [(sampleID[x], sampleID[y]) for x,y in kf.split(sampleID)]
	CV_trainID, CV_testID = zip(*kf_splits)
else:
	print('Skipping splitting samples for 5-fold cross validation...\n')

# PREP OUTPUT - print output headers to files
print('Creating file: ' + out_weight_path + '\n')
weight_cols = ['CHROM','POS', 'snpID', 'REF','ALT','TargetID','MAF','p_HWE','ES','b','beta']
pd.DataFrame(columns=weight_cols).to_csv(
	out_weight_path,
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

	# Geno, Exprr not centered since 
	## 1) DPR script centers input 
	## 2) DPR script sometimes segfaults when reading in genotype files when data was centered

	snp_annot = Geno[['snpID','POS','CHROM']]

	# 5-FOLD CROSS-VALIDATION
	if args.cvR2:
		print('Running 5-fold CV.')
		do_cv_args = [target, Geno, Expr, snp_annot, CV_trainID, CV_testID]

		k_fold_R2 = [do_cv(i, *do_cv_args) for i in range(5)]

		avg_r2_cv = sum(k_fold_R2) / 5

		print('Average R2 for 5-fold CV: {:.4f}'.format(avg_r2_cv))

		if avg_r2_cv < args.cvR2_threshold:
			print('Average R2 < ' + str(args.cvR2_threshold) + '; Skipping DPR training for TargetID: ' + target + '\n')
			return None

	else:
		avg_r2_cv = 0
		print('Skipping evaluation by 5-fold CV average R2...')


	# FINAL MODEL TRAINING
	print('Running DPR training.')

	# PREP INPUT FILES, CALL DPR, READ IN DPR OUTPUT
	bimbam = Geno[['snpID','REF','ALT', *sampleID]]
	pheno = Expr[sampleID].T

	try:
		dpr_out = prep_call_dpr(bimbam, pheno, snp_annot, dpr_file_dir, target)

	except subprocess.CalledProcessError as err:
		print('DPR failed for TargetID: ' + target + '\n')
		return None

	# FILTER FOR SNPS WITH ES!=0
	n_snp = dpr_out.ES.size

	dpr_out = dpr_out[dpr_out.ES != 0]

	n_effect_snp = dpr_out.ES.size

	# R2 CALCULATION
	out_weights = dpr_out[['snpID','ES']]

	bimbam = bimbam.drop(columns=['REF','ALT'])

	train_pvalue, train_rsquared = calc_r2(out_weights, bimbam, pheno)

	# OUTPUT TARGET WEIGHTS TO FILE
	# initialize df with MAF, pHWE, other info
	Weight = Geno[['CHROM','POS','REF','ALT','snpID','p_HWE','MAF']].copy()
	Weight = Weight[Weight.snpID.isin(dpr_out.snpID)]
	Weight['TargetID'] = target

	# merge with dpr output weights, reorder columns using existing col list
	Weight = Weight.merge(
		dpr_out, 
		left_on='snpID',
		right_on='snpID',
		how='outer')[weight_cols]

	Weight.to_csv(
		out_weight_path,
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
	print('Starting DPR training for ' + str(n_targets) + ' target genes.\n')
	pool = multiprocessing.Pool(args.thread)
	pool.imap(thread_process,[num for num in range(n_targets)])
	pool.close()
	pool.join()
	print('Done.')


############################################################
# time calculation
elapsed_sec = time()-start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)


