#!/usr/bin/env python

######################################################
# Import packages needed
import argparse
import multiprocessing
import os
import sys

from time import time

import numpy as np
import pandas as pd

# For OLS and Logistics regression
import statsmodels.api as sm
from scipy import stats
stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)

############################################################
# time calculation
start_time = time()

##########################################################
# parse input arguments
parser = argparse.ArgumentParser(description='Asso Study 01')

parser.add_argument('--gene_exp', type=str, dest='geneexp_path', required=True)
parser.add_argument('--log_file', type=str, default='')
parser.add_argument('--method', type=str, choices=['OLS','Logit'], default='OLS', 
	help='method (OLS: quantitative phenotype or multivariate phenotypes [default], Logit: dichotomous univariate phenotype')
parser.add_argument('--out_dir', type=str, default=os.getcwd())
parser.add_argument('--out_prefix', type=str, default='')
parser.add_argument('--out_twas_file', type=str, default='')
parser.add_argument('--PED', type=str, dest='ped_path', required=True)
parser.add_argument('--PED_info', type=str, dest='pedinfo_path', required=True)
parser.add_argument('--sub_dir', type=int, choices=[0, 1], default=1)
parser.add_argument('--thread', type=int, default=1)
parser.add_argument('--TIGAR_dir', type=str, help='tool directory', 
	default=os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

args = parser.parse_args()
sys.path.append(args.TIGAR_dir)
import TIGARutils as tg

##################################################
# set output file names
if not args.out_prefix:
	args.out_prefix = 'indv_' + args.method + '_TWAS'

if not args.log_file:
	args.log_file = args.out_prefix + '_log.txt'

if not args.out_twas_file:
	args.out_twas_file = args.out_prefix + '_assoc.txt'

# sub-directory in out directory
if args.sub_dir:
	out_sub_dir = os.path.join(args.out_dir, 'TWAS_indv')
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

##################################################
# Import TIGAR functions, define other functions

# For single phenotype
def regression_single(method,X,Y,Annot_df: pd.DataFrame,target):
	Result = Annot_df.copy()

	# add intercept column for design matrix
	newX = sm.add_constant(X)
	
	# regression
	if method=='OLS':
		lm = sm.OLS(Y,newX).fit()
		Result['R2'] = lm.rsquared
		
	elif method=='Logit':
		lm = sm.Logit(Y-1,newX).fit()
		Result['R2'] = lm.prsquared

	Result['BETA'] = lm.params.get(target)
	Result['BETA_SE'] = lm.bse.get(target)
	Result['T_STAT'] = lm.tvalues.get(target)
	Result['PVALUE'] = lm.pvalues.get(target)
	Result['N'] = len(X)
	
	return Result

# For multiple phenotype
def regression_multi(X,Y,Annot_df: pd.DataFrame):
	Result = Annot_df.copy()

	lm = sm.OLS(Y,X).fit()

	Result['R2'] = lm.rsquared
	Result['F_STAT'] = lm.fvalue
	Result['PVALUE'] = lm.f_pvalue
	Result['N'] = len(X)
	
	return Result

###########################################################
# Print input arguments
# out_twas_path = args.out_dir + '/indv_' + args.method + '_assoc.txt'
out_twas_path = out_sub_dir + '/' + args.out_twas_file
tmp_twas_path = out_sub_dir + '/tmp_' + args.out_twas_file


print(
'''********************************
Input Arguments
Predicted GReX data file: {geneexp_path}
PED phenotype/covariate data file: {ped_path}
PED information file: {pedinfo_path}
Regression model used for association test: {method}
Number of threads: {thread}
Output directory: {out_dir}
Output TWAS results file: {out_path}
********************************'''.format(
	**args.__dict__, 
	out_path = out_twas_path))

# tg.print_args(args)

############################################################

# get sampleIDs, ped column info, expression file info
sampleID, sample_size, exp_info, ped_cols, n_pheno, pheno, cov = tg.sampleid_startup(**args.__dict__)
print('Phenotypes to be studied: ' + ', '.join(pheno) + '\n')
print('Covariates to be used: ' + ', '.join(cov) + '\n')

print('Reading gene expression/annotation data.\n')
AnnotExp, TargetID, n_targets = tg.read_gene_annot_exp(**exp_info)

# Read in PED file
PED = pd.read_csv(
	args.ped_path,
	sep='\t',
	usecols=['IND_ID', *ped_cols])
PED = PED[PED.IND_ID.isin(sampleID)]
PED = tg.optimize_cols(PED)

# get separate Annot, Exp dataframes
Exp = (AnnotExp[sampleID]).T
Exp.columns = TargetID
Exp['IND_ID'] = Exp.index
Exp = Exp.reset_index(drop=True)

Annot = AnnotExp[AnnotExp.columns[0:5]]

###################################################
# Thread Process

# Single Phenotype
@tg.error_handler
def thread_single(num):
	target = TargetID[num]

	target_data = PEDExp[[*pheno,*cov,target]].dropna(axis=0, how='any')
	target_annot = Annot.iloc[[num]]

	X = target_data[[*cov, target]]
	Y = target_data[pheno]

	Result = regression_single('OLS',X,Y,target_annot,target)

	Result.to_csv(
		tmp_twas_path,
		sep='\t',
		header=None,
		index=None,
		mode='a')


# Multiple Phenotype
@tg.error_handler
def thread_multi(num):
	target = TargetID[num]
	target_data = Resid_Exp[[target, *pheno]].dropna(axis=0, how='any')
	target_annot = Annot.iloc[[num]]

	X = target_data[pheno]
	Y = target_data[target]

	Result = regression_multi(X,Y,target_annot)

	Result.to_csv(
		tmp_twas_path,
		sep='\t',
		header=None,
		index=None,
		mode='a')

###################################################
# Association Study
if __name__ == '__main__':
	if n_pheno == 1:
		PEDExp = PED.merge(
			Exp,
			left_on='IND_ID',
			right_on='IND_ID',
			how='outer').drop(columns=['IND_ID'])

		# output columns to dataframe
		out_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName',
			'R2','BETA','BETA_SE','T_STAT','PVALUE','N'] 
		pd.DataFrame(columns=out_cols).to_csv(
			tmp_twas_path,
			sep='\t',
			header=True,
			index=None,
			mode='w')

		pool = multiprocessing.Pool(args.thread)
		pool.imap(thread_single,[num for num in range(n_targets)])
		pool.close()
		pool.join()

	elif n_pheno > 1:
		Resid = PED[['IND_ID']].copy()
		
		for i in range(n_pheno):
			Resid[pheno[i]] = sm.OLS(PED[pheno[i]],
				sm.add_constant(PED[cov])).fit().resid.values

		Resid_Exp = Resid.merge(
			Exp,
			left_on='IND_ID',
			right_on='IND_ID',
			how='outer').drop(columns=['IND_ID'])

		# output columns to dataframe
		out_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName',
			'R2','F_STAT','PVALUE','N']   
		pd.DataFrame(columns=out_cols).to_csv(
			tmp_twas_path,
			sep='\t',
			header=True,
			index=None,
			mode='w')
		
		pool = multiprocessing.Pool(args.thread)
		pool.imap(thread_multi,[num for num in range(n_targets)])
		pool.close()
		pool.join()
	print('Done.')

	# sort output
	tg.sort_tabix_output(tmp_twas_path, out_twas_path, do_tabix=0)

############################################################
# time calculation
elapsed_sec = time() - start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)

sys.stdout.close()
