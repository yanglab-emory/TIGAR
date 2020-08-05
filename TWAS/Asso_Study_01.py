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
### time calculation
start_time = time()

##########################################################
### variables need
parser = argparse.ArgumentParser(description='Asso Study 01')

### Specify tool directory
parser.add_argument('--TIGAR_dir',type=str)

### for Gene annotation and Expression level file
parser.add_argument('--gene_exp',type=str,dest='geneexp_path')

### for PED file
parser.add_argument('--PED',type=str,dest='ped_path')

### Association Information file
parser.add_argument('--PED_info',type=str,dest='pedinfo_path')

### Method use for regression
parser.add_argument('--method',type=str)

### number of thread
parser.add_argument('--thread',type=int)

### output dir
parser.add_argument('--out_dir',type=str)

### sampleID path
parser.add_argument('--sampleID',type=str,dest='sampleid_path')

args = parser.parse_args()

sys.path.append(args.TIGAR_dir)

##################################################
## Import TIGAR functions, define other functions
import TIGARutils as tg

# For single phenotype
def regression_single(method,X,Y,annotdf: pd.DataFrame,target):
    result = annotdf.copy()

    ### add intercept column for design matrix
    newX = sm.add_constant(X)
    
    ### regression
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
# Print out variables or path using
print("Predicted GReX data file: "+args.geneexp_path + "\n")
print("PED phenotype and covariate data file using: "+args.ped_path+ "\n")
print("PED phenotype and covariates information file: "+args.pedinfo_path+ "\n")
print("Regression model used for association test: "+ args.method + "\n")
print("Number of threads: "+str(args.thread)+ "\n")
print("Output directory: "+args.out_dir+ "\n")
if args.sampleid_path:
    print("SampleID path: " + args.sampleid_path + "\n")
############################################################
# Read in PED file
PED = pd.read_csv(
    args.ped_path,
    sep='\t').set_index(['IND_ID']).rename(columns={'#FAM_ID':'FAM_ID'})
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
Asso_Info=pd.read_csv(
    args.pedinfo_path,
    sep='\t',
    header=None,
    names=['Ind','Var'])

### phenotype
pheno = Asso_Info[Asso_Info.Ind=='P'].Var.values
n_pheno = pheno.size

if not n_pheno:
    raise SystemExit("No phenotype column name is provided by --PED_info.")

print("Phenotypes to be studied: ", pheno)

### covariates
cov = Asso_Info[Asso_Info.Ind=='C'].Var.values
if not cov.size:
    raise SystemExit("No covariates provided.")

print("Covariates to be used: ",cov)

TargetID = GeneAnnotExp.TargetID.values
n_targets = TargetID.size

if not n_targets:
    raise SystemExit("There is no GREx data in gene expression file provided by --gene_exp ")


pedannot_sampleids = np.intersect1d(PED.IND_ID.values, GeneAnnotExp.columns[5:].values)

# if user specified path with sampleids, and at least one sampid from that file, get intersection
if args.sampleid_path:
    spec_sampleids = pd.read_csv(args.sampleid_path, sep='\t', header=None)
    spec_sampleids = spec_sampleids[0].drop_duplicates()

    sampleID = np.intersect1d(pedannot_sampleids, spec_sampleids)

else:
    sampleID = pedannot_sampleids 

if not sampleID.size:
    raise SystemExit("There is no overlapped sample IDs between gene expression file and PED file.")

# Organizing PED and Gene-expression file
ped_cols = np.concatenate((pheno, cov))
PED = PED[PED.index.isin(sampleID)][ped_cols]

gene_cols = np.concatenate((['CHROM', 'GeneStart', 'GeneEnd', 'TargetID', 'GeneName'], sampleID))
GeneAnnotExp = GeneAnnotExp[gene_cols]

GeneExp = (GeneAnnotExp.set_index(['TargetID'])[sampleID]).T
GeneAnnot = GeneAnnotExp[GeneAnnotExp.columns[0:5]]

###################################################
# Thread Process

# Single Phenotype
def thread_single(num):
    try:
        target = TargetID[num]
        target_cols = np.concatenate((pheno, cov, [target]))
        target_data = PEDExp[target_cols].dropna(axis=0, how='any')
        target_annot = GeneAnnot.iloc[[num]]

        X_cols = np.append(cov, target)
        X = target_data[X_cols]
        Y = target_data[pheno]

        out = regression_single('OLS',X,Y,target_annot,target)

        out.to_csv(
            args.out_dir+"/indv_" + args.method+"_assoc.txt",
            sep='\t',
            header=None,
            index=None,
            mode='a')

    except Exception as e:
        e_type, e_obj, e_tracebk = sys.exc_info()
        e_line_num = e_tracebk.tb_lineno

        e, e_type, e_line_num = [str(x) for x in [e, e_type, e_line_num]]

        print('Caught a type '+ e_type +' exception for TargetID='+target+', num=' + str(num) + ' on line'+e_line_num+':\n' + e )

    finally:
        # print info to log do not wait for buffer to fill up
        sys.stdout.flush()

# Multiple Phenotype
def thread_multi(num):
    try:
        target = TargetID[num]
        target_cols = np.append(pheno,target)
        target_data = resid_exp[target_cols].dropna(axis=0, how='any')
        target_annot = GeneAnnot.iloc[[num]]

        X = target_data[pheno]
        Y = target_data[target]

        out = regression_multi(X,Y,target_annot)

        out.to_csv(
            args.out_dir+"/indv_" + args.method+"_assoc.txt",
            sep='\t',
            header=None,
            index=None,
            mode='a')

    except Exception as e:
        e_type, e_obj, e_tracebk = sys.exc_info()
        e_line_num = e_tracebk.tb_lineno
        
        e, e_type, e_line_num = [str(x) for x in [e, e_type, e_line_num]]

        print('Caught a type '+ e_type +' exception for TargetID='+target+', num=' + str(num) + ' on line '+e_line_num+':\n' + e )

    finally:
        # print info to log do not wait for buffer to fill up
        sys.stdout.flush()

###################################################
# Association Study
if __name__ == '__main__':
    if n_pheno == 1:
        PEDExp = PED.merge(
            GeneExp,
            left_index=True, 
            right_index=True,
            how='outer')

        # output columns to dataframe
        out_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName',
            'R2','BETA','BETA_SE','T_STAT','PVALUE','N'] 
        pd.DataFrame(columns=out_cols).to_csv(
                args.out_dir+"/indv_"+args.method+"_assoc.txt",
                sep='\t',
                header=True,
                index=None,
                mode='w')

        pool = multiprocessing.Pool(args.thread)
        pool.map(thread_single,[num for num in range(n_targets)])
        pool.close()
        pool.join()

    elif n_pheno > 1:
        resid = pd.DataFrame(index=PED.index.copy())
        
        for i in range(n_pheno):
            resid[pheno[i]] = sm.OLS(PED[pheno[i]],
                sm.add_constant(PED[cov])).fit().resid.values

        resid_exp = resid.merge(
            GeneExp,
            left_index=True,
            right_index=True,
            how='outer')

        # output columns to dataframe
        out_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName',
            'R2','F_STAT','PVALUE','N']   
        pd.DataFrame(columns=out_cols).to_csv(
                args.out_dir+"/indv_"+args.method+"_assoc.txt",
                sep='\t',
                header=True,
                index=None,
                mode='w')
        
        pool = multiprocessing.Pool(args.thread)
        pool.map(thread_multi,[num for num in range(n_targets)])
        pool.close()
        pool.join()

############################################################
### time calculation
elapsed_sec = time()-start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print("Computation time (DD:HH:MM:SS): " + elapsed_time)



