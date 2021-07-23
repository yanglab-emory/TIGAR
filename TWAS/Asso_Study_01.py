#!/usr/bin/env python

######################################################
# Import packages needed
import argparse
import multiprocessing
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

# Specify tool directory
parser.add_argument('--TIGAR_dir' ,type=str)

# Gene annotation and Expression level file
parser.add_argument('--gene_exp' ,type=str ,dest='geneexp_path')

# PED file path 
parser.add_argument('--PED' ,type=str ,dest='ped_path')

# Association Information file path
parser.add_argument('--PED_info' ,type=str ,dest='pedinfo_path')

# Method to use for regression
parser.add_argument('--method' ,type=str)

# number of thread
parser.add_argument('--thread' ,type=int)

# output dir
parser.add_argument('--out_dir' ,type=str)

args = parser.parse_args()

sys.path.append(args.TIGAR_dir)

##################################################
# Import TIGAR functions, define other functions
import TIGARutils as tg

# For single phenotype
def regression_single(method,X,Y,annotdf: pd.DataFrame,target):
	result = annotdf.copy()

	# add intercept column for design matrix
	newX = sm.add_constant(X)
	
	# regression
	if method=='OLS':
		lm = sm.OLS(Y,newX).fit()
		result['R2'] = lm.rsquared
		
	elif method=='Logit':
		lm = sm.Logit(Y-1,newX).fit()
		result['R2'] = lm.prsquared

	result['BETA'] = lm.params.get(target)
	result['BETA_SE'] = lm.bse.get(target)
	result['T_STAT'] = lm.tvalues.get(target)
	result['PVALUE'] = lm.pvalues.get(target)
	result['N'] = len(X)
	
	return result

# For multiple phenotype
def regression_multi(X,Y,annotdf: pd.DataFrame):
	result = annotdf.copy()

	lm = sm.OLS(Y,X).fit()

	result['R2'] = lm.rsquared
	result['F_STAT'] = lm.fvalue
	result['PVALUE'] = lm.f_pvalue
	result['N'] = len(X)
	
	return result

###########################################################
# Print input arguments
out_twas_path = args.out_dir + '/indv_' + args.method + '_assoc.txt'

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

############################################################
# Read in PED file
PED = pd.read_csv(
	args.ped_path,
	sep='\t').rename(columns={'#FAM_ID':'FAM_ID'})
PED = tg.optimize_cols(PED)

# Gene annotation and Expression level file
GeneAnnotExp = pd.read_csv(
	args.geneexp_path,
	sep='\t',
	low_memory=False)
GeneAnnotExp = tg.optimize_cols(GeneAnnotExp)

# Read in Association information
# P:phenotype
# C:covariate
Asso_Info = pd.read_csv(
	args.pedinfo_path,
	sep='\t',
	header=None,
	names=['Ind','Var'])

# phenotype
pheno = Asso_Info[Asso_Info.Ind=='P'].Var.values
n_pheno = pheno.size

if not n_pheno:
	raise SystemExit('No phenotype column name is provided by --PED_info.')

print('Phenotypes to be studied: ' + ', '.join([x for x in pheno.tolist()]) + '\n')

# covariates
cov = Asso_Info[Asso_Info.Ind=='C'].Var.values
if not cov.size:
	raise SystemExit('No covariates provided.')

print('Covariates to be used: ' + ', '.join([x for x in cov.tolist()]) + '\n')

TargetID = GeneAnnotExp.TargetID
n_targets = TargetID.size

if not n_targets:
	raise SystemExit('There is no GREx data in gene expression file provided by --gene_exp ')

sampleID = np.intersect1d(PED.IND_ID, GeneAnnotExp.columns[5:])

if not sampleID.size:
	raise SystemExit('There is no overlapped sample IDs between gene expression file and PED file.')

# Organizing PED and Gene-expression file
ped_cols = np.concatenate((['IND_ID'], pheno, cov))
PED = PED[PED.IND_ID.isin(sampleID)][ped_cols]

gene_cols = ['CHROM', 'GeneStart', 'GeneEnd', 'TargetID', 'GeneName',*sampleID]
GeneAnnotExp = GeneAnnotExp[gene_cols]


GeneExp = (GeneAnnotExp[sampleID]).T
GeneExp.columns = TargetID
GeneExp['IND_ID'] = GeneExp.index
GeneExp = GeneExp.reset_index(drop=True)

GeneAnnot = GeneAnnotExp[GeneAnnotExp.columns[0:5]]

###################################################
# Thread Process

# Single Phenotype
@tg.error_handler
def thread_single(num):
	target = TargetID[num]
	target_cols = np.concatenate((pheno, cov, [target]))
	target_data = PEDExp[target_cols].dropna(axis=0, how='any')
	target_annot = GeneAnnot.iloc[[num]]

	X_cols = np.append(cov, target)
	X = target_data[X_cols]
	Y = target_data[pheno]

	out = regression_single('OLS',X,Y,target_annot,target)

	out.to_csv(
		out_twas_path,
		sep='\t',
		header=None,
		index=None,
		mode='a')


# Multiple Phenotype
@tg.error_handler
def thread_multi(num):
	target = TargetID[num]
	target_cols = np.insert(pheno,0,target)
	target_data = resid_exp[target_cols].dropna(axis=0, how='any')
	target_annot = GeneAnnot.iloc[[num]]

	X = target_data[pheno]
	Y = target_data[target]

	out = regression_multi(X,Y,target_annot)

	out.to_csv(
		out_twas_path,
		sep='\t',
		header=None,
		index=None,
		mode='a')

###################################################
# Association Study
if __name__ == '__main__':
	if n_pheno == 1:
		PEDExp = PED.merge(
			GeneExp,
			left_on='IND_ID',
			right_on='IND_ID',
			how='outer').drop(columns=['IND_ID'])

		# output columns to dataframe
		out_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName',
			'R2','BETA','BETA_SE','T_STAT','PVALUE','N'] 
		pd.DataFrame(columns=out_cols).to_csv(
				out_twas_path,
				sep='\t',
				header=True,
				index=None,
				mode='w')

		pool = multiprocessing.Pool(args.thread)
		pool.imap(thread_single,[num for num in range(n_targets)])
		pool.close()
		pool.join()

	elif n_pheno > 1:
		resid = PED[['IND_ID']].copy()
		
		for i in range(n_pheno):
			resid[pheno[i]] = sm.OLS(PED[pheno[i]],
				sm.add_constant(PED[cov])).fit().resid.values

		resid_exp = resid.merge(
			GeneExp,
			left_on='IND_ID',
			right_on='IND_ID',
			how='outer').drop(columns=['IND_ID'])

		# output columns to dataframe
		out_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName',
			'R2','F_STAT','PVALUE','N']   
		pd.DataFrame(columns=out_cols).to_csv(
				out_twas_path,
				sep='\t',
				header=True,
				index=None,
				mode='w')
		
		pool = multiprocessing.Pool(args.thread)
		pool.imap(thread_multi,[num for num in range(n_targets)])
		pool.close()
		pool.join()
	print('Done.')

############################################################
# time calculation
elapsed_sec = time()-start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)



