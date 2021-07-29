#!/usr/bin/env python

#########################################################
# Import packages needed
import argparse
import multiprocessing
import operator
import subprocess
import sys
import warnings

from io import StringIO
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

# Specify tool directory
parser.add_argument('--TIGAR_dir', type=str)

# Gene Annotation and Expression level file path
parser.add_argument('--gene_exp', type=str, dest='geneexp_path')

# Training sampleID path
parser.add_argument('--train_sampleID', type=str, dest='sampleid_path')

# Specified chromosome number
parser.add_argument('--chr', type=str, dest='chrm')

# Training genotype file path
parser.add_argument('--genofile', type=str, dest='geno_path')

# Specified input file type(vcf or dosages)
parser.add_argument('--genofile_type', type=str)

# 'DS' or 'GT' for VCF genotype file
parser.add_argument('--format', type=str, dest='data_format')

# window
parser.add_argument('--window', type=int)

# missing rate: threshold for excluding SNPs with too many missing values
parser.add_argument('--missing_rate', type=float)

# Folded Minor Allele Frequency (range from 0-0.5) threshold
parser.add_argument('--maf', type=float)

# p-value for Hardy Weinberg Equilibrium exact test threshold
parser.add_argument('--hwe', type=float)

# cvR2 (0: no CV, 1: CV)
parser.add_argument('--cvR2', type=int)

# threshold cvR2 value for training (default: 0.005)
parser.add_argument('--cvR2_threshold', type=float)

# k-fold cross validation (for EN model training)
parser.add_argument('--cv', type=int)

# Ratio of L1 and L2 (for EN model training)
parser.add_argument('--alpha', type=float)

# Use specified args.alpha? (0: [0.1, 0.5, 0.9, 1], 1: args.alpha)
parser.add_argument('--use_alpha', type=int)

# file paths
parser.add_argument('--out_weight_file', type=str)
parser.add_argument('--out_info_file', type=str)

# Number of thread
parser.add_argument('--thread', type=int)

# output dir
parser.add_argument('--out_dir', type=str)

args = parser.parse_args()

args.alpha = [0.1, 0.5, 0.9, 1] if not args.use_alpha else args.alpha

sys.path.append(args.TIGAR_dir)

###############################################################
# DEFINE, IMPORT FUNCTIONS
import TIGARutils as tg

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
def do_cv(i, target_geno_exp_df, cv_trainID, cv_testID):
	target_geno_exp_df = target_geno_exp_df.copy()
	train_geno_exp = target_geno_exp_df.loc[cv_trainID[i]].dropna()
	test_geno_exp = target_geno_exp_df.loc[cv_testID[i]].dropna()

	cv_rsquared = elastic_net(train_geno_exp, test_geno_exp)
	return cv_rsquared

###############################################################
# check input arguments
if args.genofile_type == 'vcf':
	if (args.data_format != 'GT') and (args.data_format != 'DS'):
		raise SystemExit('Please specify the genotype data format used by the vcf file (--format ) as either "GT" or "DS".\n')
		
elif args.genofile_type == 'dosage':
	args.data_format = 'DS'

else:
	raise SystemExit('Please specify the type input genotype file type (--genofile_type) as either "vcf" or "dosage".\n')
	
out_train_weight_path = args.out_dir + '/temp_' + args.out_weight_file
out_train_info_path = args.out_dir + '/' +  args.out_info_file

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
	out_weight = out_train_weight_path,
	out_info = out_train_info_path))

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
# sampleID, sample_size, exp_info, geno_info = tg.train_startup(**args.__dict__)
sampleID, sample_size, geno_info, exp_info = tg.sampleid_startup(**args.__dict__)

# Read in gene expression info
print('Reading gene expression data.\n')
GeneExp, TargetID, n_targets = tg.read_gene_annot_exp(**args.__dict__, **exp_info)

# PREP CROSS VALIDATION SAMPLES - Split sampleIDs for cross validation
if args.cvR2:
	print('Splitting sample IDs randomly for 5-fold cross validation by average R2...\n')
	kf = KFold(n_splits=5)
	kf_splits = [(sampleID[x], sampleID[y]) for x,y in kf.split(sampleID)]
	CV_trainID, CV_testID = zip(*kf_splits)

else:
	print('Skipping splitting samples for 5-fold cross validation...\n')

# PREP OUTPUT - print output headers to files
print('Creating file: ' + out_train_weight_path + '\n')
weight_out_cols = ['CHROM','POS','ID','REF','ALT','TargetID','MAF','p_HWE','ES']
pd.DataFrame(columns=weight_out_cols).to_csv(
	out_train_weight_path,
	header=True,
	index=None,
	sep='\t',
	mode='w')

print('Creating file: ' + out_train_info_path + '\n')
info_out_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName','sample_size','n_snp','n_effect_snp','CVR2','TrainPVALUE','TrainR2','k-fold','alpha','Lambda','cvm','CVR2_threshold']
pd.DataFrame(columns=info_out_cols).to_csv(
	out_train_info_path,
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
	target_exp = GeneExp.iloc[[num]]

	start = str(max(int(target_exp.GeneStart) - args.window,0))
	end = str(int(target_exp.GeneEnd) + args.window)

	# READ IN AND PROCESS GENOTYPE DATA 
	# file must be bgzipped and tabix
	# target_geno = tg.read_genotype(start, end, sampleID, **geno_info, **args.__dict__)
	target_geno = tg.read_tabix(start, end, sampleID, **geno_info)

	# filter out variants that exceed missing rate threshold
	target_geno = tg.handle_missing(target_geno, sampleID, args.missing_rate)

	# get, filter maf
	target_geno = tg.calc_maf(target_geno, sampleID, args.maf, op=operator.ge)

	# get, filter p_HWE
	target_geno = tg.calc_p_hwe(target_geno, sampleID, args.hwe, op=operator.ge)

	n_snp = target_geno['snpID'].size

	# center data
	target_geno = tg.center(target_geno, sampleID)
	target_exp = tg.center(target_exp, sampleID)

	# merge geno, expression files, transpose
	target_geno_exp = pd.concat([
		target_geno.set_index(['snpID'])[sampleID],
		target_exp.set_index(['TargetID'])[sampleID]
		]).T

	# 5-FOLD CROSS-VALIDATION
	# Evaluate whether Elastic Net model valid
	if args.cvR2:
		print('Running 5-fold CV.')
		# print('Starting with Elastic-Net penalty parameter')
		do_cv_args = [target_geno_exp, CV_trainID, CV_testID]
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
	# initialize target_weights dataframe
	target_weights = target_geno[['CHROM','POS','snpID','REF','ALT','p_HWE','MAF']].copy()

	target_weights['TargetID'] = target

	# do elastic net training
	target_weights['ES'], R2, Pvalue, Alpha, Lambda, cvm = elastic_net(target_geno_exp)

	# filter
	target_weights = target_weights[target_weights['ES']!=0]

	# reorder columns for output
	target_weights = target_weights[['CHROM','POS','snpID','REF','ALT','TargetID','MAF','p_HWE','ES']]

	target_weights.to_csv(
		out_train_weight_path,
		header=False,
		index=None,
		sep='\t',
		mode='a')

	# output training information, result from elastic net
	train_info = target_exp[['CHROM','GeneStart','GeneEnd','TargetID','GeneName']].copy()

	train_info['sample_size'] = sample_size
	train_info['n_snp'] = n_snp
	train_info['n_effect_snp'] = n_effect_snp = target_weights.ES.size
	train_info['CVR2'] = avg_r2_cv
	train_info['TrainPVALUE'] = Pvalue if not np.isnan(Pvalue) else 'NaN'
	train_info['TrainR2'] = R2 if n_effect_snp else 0
	train_info['k_fold'] = args.cv
	train_info['alpha'] = Alpha
	train_info['lambda'] = Lambda
	train_info['cvm'] = cvm
	train_info['CVR2_threshold'] = args.cvR2_threshold if args.cvR2 else 0

	train_info.to_csv(
		out_train_info_path,
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

###############################################################
# time calculation
elapsed_sec = time()-start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)

