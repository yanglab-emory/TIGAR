#!/usr/bin/env python

#########################################################
# Import packages needed
import argparse
import multiprocessing
import operator
import os
import subprocess
import sys
import warnings

from time import time

import numpy as np
import pandas as pd

### import grid search for model selection
from sklearn.model_selection import KFold

### For Elastic Net Regression
from sklearn.linear_model import ElasticNetCV

### For OLS regression in cross validation
from scipy import stats
import statsmodels.api as sm

warnings.filterwarnings('ignore')

###############################################################
# time calculation
start_time = time()

###############################################################
# parse input arguments
parser = argparse.ArgumentParser(description='Elastic Net Training')

parser.add_argument('--alpha', type=float, default=0.5, 
	help='ratio of L1 and L2 for EN model training (default: 0.5)')
parser.add_argument('--chr', type=str, dest='chrm', 
	choices=[str(i + 1) for i in range(22)],
	required=True, 
	help='chromosome number')
parser.add_argument('--cv', type=int, default=5, 
	help='k for k-fold cross validation for EN model training (default: 5)')
parser.add_argument('--cvR2', type=int, choices=[0, 1], default=1, 
	help='cvR2 (0: no CV, 1: CV [default])')
parser.add_argument('--cvR2_threshold', type=float, default=0.005, 
	help='threshold cvR2 value for training (default: 0.005)')
parser.add_argument('--format', type=str, dest='data_format', choices=['GT', 'DS'], default='GT', 
	help='data format of VCF genotype data (DS, GT [default])')
parser.add_argument('--gene_exp', type=str, dest='geneexp_path', required=True)
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
parser.add_argument('--out_info_file', type=str, default='')
parser.add_argument('--out_prefix', type=str, default='')
parser.add_argument('--out_weight_file', type=str, default='')
parser.add_argument('--sub_dir', type=int, choices=[0, 1], default=1)
parser.add_argument('--thread', type=int, default=1)
parser.add_argument('--TIGAR_dir', type=str, help='tool directory', 
	default=os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
parser.add_argument('--train_sampleID', type=str, dest='sampleid_path', required=True)
parser.add_argument('--use_alpha', type=int, default=1, 
	help='use specified alpha value? (0: [0.1, 0.5, 0.9, 1], 1: --alpha value [default])')
parser.add_argument('--window', type=int, default=1000000, 
	help='size around gene region for selecting cis-SNPs for fitting gene expression prediction model (default: 1000000 [ie, +-1MB region around gene])')

args = parser.parse_args()
sys.path.append(args.TIGAR_dir)
import TIGARutils as tg

#############################################################
args.alpha = [0.1, 0.5, 0.9, 1] if not args.use_alpha else args.alpha

# set output file names
if not args.out_prefix:
	args.out_prefix = 'CHR' + args.chrm + '_EN_train'

if not args.log_file:
	args.log_file = args.out_prefix + '_log.txt'

if not args.out_info_file:
	args.out_info_file = args.out_prefix + '_GeneInfo.txt'

if not args.out_weight_file:
	args.out_weight_file = args.out_prefix + '_eQTLweights.txt'

# sub-directory in out directory
if args.sub_dir:
	out_sub_dir = os.path.join(args.out_dir, 'EN_CHR' + args.chrm)
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

###############################################################
# DEFINE, IMPORT FUNCTIONS

### Elastic Net
### Input: 
### 1.train: training data, independent variables values for training(SNPs Genotype) with last column of response variable values(Gene Expression Level)
### 2.test (optional): testing data, same format as training.
### 3.Alpha (optional): ratio for L1 and L2 penalty in elastic net regression,default=args.alpha (default args.alpha=0.5)
###          Alpha=0: Lasso Regression
###          Alpha=1: Ridge Regression
###          0 < Alpha < 1: Elastic Net Regression
### 4.k (optional): k-fold cross validation,default=args.cv (default args.cv=5)
### Return:
### if test specified: 
### 1. returns 5-folds cross validation R2 calculation
### if test not specified:
### 1.Regression coefficent for training data
### 2.Training R-square
### 3.Parameter selected by cross validation
### 4.Corresponding mean cross validation score

### Using grid search and cross-validation to find the best lambda(penalty)
def elastic_net(train, test=None, k=args.cv, Alpha=args.alpha):
	train = train.copy()
	trainX = train.iloc[:,0:-1]
	trainY = train.iloc[:,-1]

	if test is not None:
		test = test.copy()
		testX = test.iloc[:,0:-1]
		testY = test.iloc[:,-1]

	else:
		testX = trainX
		testY = trainY

	reg = ElasticNetCV(
		l1_ratio=Alpha,
		fit_intercept=False,
		alphas=np.arange(0,1.01,0.01),
		selection='random',
		cv=k).fit(trainX,trainY)

	Lambda = reg.alpha_
	cvm = np.min(reg.mse_path_)
	beta = reg.coef_

	predY = reg.predict(testX)

	lm = sm.OLS(testY, sm.add_constant(predY)).fit()

	Rsquared = lm.rsquared

	if test is not None:
		return Rsquared

	Pvalue = lm.f_pvalue

	return beta, Rsquared, Pvalue, Alpha, Lambda, cvm


# function to do the ith cross validation step
def do_cv(i, Geno_Exp_df, cv_trainID, cv_testID):
	Geno_Exp_df = Geno_Exp_df.copy()
	train_geno_exp = Geno_Exp_df.loc[cv_trainID[i]].dropna()
	test_geno_exp = Geno_Exp_df.loc[cv_testID[i]].dropna()

	cv_rsquared = elastic_net(train_geno_exp, test_geno_exp)
	return cv_rsquared

###############################################################	
tmp_weight_path = out_sub_dir + '/temp_' + args.out_weight_file
out_weight_path = out_sub_dir + '/' + args.out_weight_file
out_info_path = out_sub_dir + '/' +  args.out_info_file

###############################################################
# Print input arguments
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
{cvR2_str1} Elastic-Net model by 5-fold cross validation{cvR2_str2}.
Number of cross-validation folds used to tune Elastic-Net penalty parameter (lambda): {cv}
Ratio for L1 & L2 penalty used by Elastic-Net regression (alpha): {alpha}
Number of threads: {thread}
Output directory: {out_dir}
Output training weights file: {out_weight}
Output training info file: {out_info}
********************************'''.format(
	**args.__dict__,
	cvR2_str1 = {0:'Skipping evaluation of', 1:'Evaluating'}[args.cvR2],
	cvR2_str2 = {0:'', 1:' with inclusion threshold Avg.CVR2 >' + str(args.cvR2_threshold)}[args.cvR2],
	out_weight = out_weight_path,
	out_info = out_info_path))

# tg.print_args(args)

###############################################################
# Training Processing
### Read in Gene annotation and Expression level file (text file)
### First five columns should be fixed:
### 1.Chrmosome Number
### 2.GeneStart Posistion
### 3.GeneEnd Position
### 4.TargetID (i.e.GeneID, treated as unique annotation for each gene)
### 5.Gene Name

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
print('Creating file: ' + tmp_weight_path + '\n')
weight_cols = ['CHROM','POS','snpID','REF','ALT','TargetID','MAF','p_HWE','ES']
pd.DataFrame(columns=weight_cols).to_csv(
	tmp_weight_path,
	header=True,
	index=None,
	sep='\t',
	mode='w')

print('Creating file: ' + out_info_path + '\n')
info_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName','sample_size','n_snp','n_effect_snp','CVR2','TrainPVALUE','TrainR2','k-fold','alpha','Lambda','cvm','CVR2_threshold']
pd.DataFrame(columns=info_cols).to_csv(
	out_info_path,
	header=True,
	index=None,
	sep='\t',
	mode='w')

print('********************************\n')

###############################################################
# thread function
@tg.error_handler
def thread_process(num):
	target = TargetID[num]
	print('num=' + str(num) + '\nTargetID=' + target)
	Expr = GeneExp.iloc[[num]]

	start = str(max(int(Expr.GeneStart) - args.window,0))
	end = str(int(Expr.GeneEnd) + args.window)

	# READ IN AND PROCESS GENOTYPE DATA 
	# file must be bgzipped and tabix
	Geno = tg.read_tabix(start, end, sampleID, **geno_info)

	# filter out variants that exceed missing rate threshold
	Geno = tg.handle_missing(Geno, sampleID, args.missing_rate)

	# get, filter maf
	Geno = tg.calc_maf(Geno, sampleID, args.maf, op=operator.ge)

	# get, filter p_HWE
	Geno = tg.calc_p_hwe(Geno, sampleID, args.hwe, op=operator.ge)

	n_snp = Geno['snpID'].size

	# center data
	Geno = tg.center(Geno, sampleID)
	Expr = tg.center(Expr, sampleID)

	# merge geno, expression files, transpose
	Geno_Exp = pd.concat([
		Geno.set_index(['snpID'])[sampleID],
		Expr.set_index(['TargetID'])[sampleID]
		]).T

	# 5-FOLD CROSS-VALIDATION
	# Evaluate whether Elastic Net model valid
	if args.cvR2:
		print('Running 5-fold CV.')
		# print('Starting with Elastic-Net penalty parameter')
		do_cv_args = [Geno_Exp, CV_trainID, CV_testID]
		k_fold_R2 = [do_cv(i, *do_cv_args) for i in range(5)]

		avg_r2_cv = sum(k_fold_R2) / 5

		print('Average R2 for 5-fold CV: {:.4f}'.format(avg_r2_cv))

		if avg_r2_cv < args.cvR2_threshold:
			print('Average R2 < ' + str(args.cvR2_threshold) + '; Skipping Elastic-Net training for TargetID: ' + target + '\n')
			return None

	else:
		avg_r2_cv = 0
		print('Skipping evaluation by 5-fold CV average R2...')

	# FINAL MODEL TRAINING
	print('Running Elastic-Net training.')
	# initialize Weight dataframe
	Weight = Geno[['CHROM','POS','snpID','REF','ALT','p_HWE','MAF']].copy()

	Weight['TargetID'] = target

	# do elastic net training
	Weight['ES'], R2, Pvalue, Alpha, Lambda, cvm = elastic_net(Geno_Exp)

	# filter
	Weight = Weight[Weight['ES'] != 0]
	n_effect_snp = Weight.ES.size

	# reorder columns for output
	Weight = Weight[weight_cols]

	Weight.to_csv(
		tmp_weight_path,
		header=False,
		index=None,
		sep='\t',
		mode='a')

	# output training information, result from elastic net
	Info = Expr[['CHROM','GeneStart','GeneEnd','TargetID','GeneName']].copy()

	Info['sample_size'] = sample_size
	Info['n_snp'] = n_snp
	Info['n_effect_snp'] = n_effect_snp
	Info['CVR2'] = avg_r2_cv
	Info['TrainPVALUE'] = Pvalue if not np.isnan(Pvalue) else 'NaN'
	Info['TrainR2'] = R2 if n_effect_snp else 0
	Info['k_fold'] = args.cv
	Info['alpha'] = Alpha
	Info['lambda'] = Lambda
	Info['cvm'] = cvm
	Info['CVR2_threshold'] = args.cvR2_threshold if args.cvR2 else 0

	Info.to_csv(
		out_info_path,
		header=None,
		index=None,
		sep='\t',
		mode='a')

	print('Target Elastic-Net training completed.\n')


###############################################################
# start thread process
if __name__ == '__main__':
	print('Starting Elastic-Net training for ' + str(n_targets) + ' target genes.\n')
	pool = multiprocessing.Pool(args.thread)
	pool.imap(thread_process,[num for num in range(n_targets)])
	pool.close()
	pool.join()
	print('Done.')

	# sort tabix output
	tg.sort_tabix_output(tmp_weight_path, out_weight_path)



###############################################################
# time calculation
elapsed_sec = time()-start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)

sys.stdout.close()
