#!/usr/bin/env python

##############################################################################################
# Import packages needed
import argparse
import warnings
import time
import pandas as pd
import numpy as np
from numpy import *
from dfply import *
from scipy import stats
import multiprocessing

# For OLS and Logistics regression
import statsmodels.api as sm
from scipy import stats
stats.chisqprob = lambda chisq, df: stats.chi2.sf(chisq, df)

##############################################################################################
# For single phenotype
def regression_single(method,X,Y,TargetID):
    ### add intercept column for design matrix
    newX = sm.add_constant(X)
    
    result = pd.DataFrame()
    result['TargetID'] = np.array(TargetID).ravel()
    
    ### regression
    if method=='OLS':
        lm = sm.OLS(Y,newX).fit()
        result['R2'] = np.array(lm.rsquared)
        
    elif method=='Logit':
        lm = sm.Logit(Y-1,newX).fit()
        result['R2'] = np.array(lm.prsquared)
    
    
    result['BETA'] = pd.DataFrame(lm.params).loc[TargetID]
    result['BETA_SE'] = pd.DataFrame(lm.bse).loc[TargetID]
    result['T_STAT'] = pd.DataFrame(lm.tvalues).loc[TargetID]
    result['PVALUE'] = pd.DataFrame(lm.pvalues).loc[TargetID]
    result['N'] = np.array(len(X))
    
    return result


###############################################################################################
# For multiple phenotype
def regression_multi(X,Y,TargetID):
    lm = sm.OLS(Y,X).fit()
    
    result = pd.DataFrame()
    result['TargetID'] = np.array(TargetID).ravel()
    result['R2'] = np.array(lm.rsquared).ravel()
    result['F_STAT'] = np.array(lm.fvalue).ravel()
    result['F_PVALUE'] = np.array(lm.f_pvalue).ravel()
    result['N'] = np.array(len(X))
    
    return result

######################################################################################################
### variables need
parser = argparse.ArgumentParser(description='manual to this script')

### for Gene annotation and Expression level file
parser.add_argument('--Gene_Exp_path',type=str,default = None)

### for PED file
parser.add_argument('--PED',type=str,default = None)

### Association Information file
parser.add_argument('--Asso_Info',type=str,default = None)

### Method use for regression
parser.add_argument('--method',type=str,default = None)

### number of thread
parser.add_argument('--thread',type=int,default = None)

### output dir
parser.add_argument('--out_prefix',type=str,default=None)

args = parser.parse_args()

#####################################################################################################
# Print out variables or path using
print("Gene-expression file:"+args.Gene_Exp_path)
print("PED file using:"+args.PED)
print("Assosiation Information file:"+args.Asso_Info)
print("Regression method using:"+args.method)
print("Number of thread:"+str(args.thread))
print("Output dir:"+args.out_prefix)
####################################################################################################
# Read in PED file
PED = pd.read_csv(args.PED,sep='\t').rename(columns={'#FAM_ID':'FAM_ID'})

# Gene annotation and Expression level file
Genecode = pd.read_csv(args.Gene_Exp_path,sep='\t')

# Read in Association information
# P:phenotype
# C:covariate
Asso_Info=pd.read_csv(args.Asso_Info,sep='\t',header=None)
Asso_Info.columns=['Ind','Var']

### phenotype
pheno = Asso_Info >> mask(Asso_Info.Ind=='P') >> select(Asso_Info.Var) 
pheno = np.array(pheno).ravel()
if len(pheno)==0:
    raise SystemExit("No phenotype provided.")
else:
    print("Phenotype Used:",pheno)
### covariance
cov = Asso_Info >> mask(Asso_Info.Ind=='C') >> select(Asso_Info.Var)
cov = np.array(cov).ravel()
if len(cov)==0:
    raise SystemExit("No covariates provided.")
else:
    print("Covariance used:",cov)

TargetID = np.array(Genecode.TargetID)
if len(TargetID)==0:
    raise SystemExit("No Gene in gene-expression file")

sampleID = np.intersect1d(np.array(PED.IND_ID),np.array(Genecode.columns[5:]))
if len(sampleID)==0:
    raise SystemExit("No overlaping IND_ID between gene-expression and PED file.")

# Organizing PED and Gene-expression file
PED = PED >> mask(PED.IND_ID.isin(sampleID)) >> select(PED.IND_ID,PED[pheno],PED[cov])
PED[PED=='X']=np.nan

Genecode=Genecode >> select(Genecode[Genecode.columns[0:5]],Genecode[sampleID])

Gene_temp = (Genecode[Genecode.columns[4:]]).T
Gene_temp.columns = Gene_temp.loc['TargetID']
Gene_temp = Gene_temp.drop(['TargetID'])
Gene_temp['IND_ID'] = Gene_temp.index
Gene_temp = Gene_temp.reset_index(drop=True)

##################################################################################################
# Thread Process

# Single Phenotype
def thread_single(num):
    Target_temp = (Target >> select(Target[pheno],Target[cov],Target[TargetID[num]])).dropna(axis=0,how='any')
    X = (Target >> select(Target[cov],Target[TargetID[num]]))
    Y = Target[pheno]
    
    lm = regression_single(args.method,X,Y,TargetID[num])
        
    Gene_annot = Genecode >> mask(Genecode.TargetID==TargetID[num]) >> select(Genecode.columns[0:5])
    
    out = Gene_annot.merge(lm,left_on='TargetID',right_on='TargetID',how='outer')
    
    out.to_csv(args.out_prefix+"/association_study_"+args.method+".txt",sep='\t',header=None,index=None,mode='a')

# Multiple Phenotype
def thread_multi(num):
    Target_temp = Target >> select(Target[TargetID[num]],Target[pheno])
    Target_temp = pd.DataFrame(Target_temp.dropna(axis=0,how='any'),dtype='float')
    
    X = Target_temp[pheno]
    Y = Target_temp[TargetID[num]]
    
    lm = regression_multi(X,Y,TargetID[num])
    
    Gene_annot = Genecode >> mask(Genecode.TargetID==TargetID[num]) >> select(Genecode.columns[0:5])
    out = Gene_annot.merge(lm,left_on='TargetID',right_on='TargetID',how='outer')
    
    out.to_csv(args.out_prefix+"/association_study_"+args.method+".txt",sep='\t',header=None,index=None,mode='a')

###################################################################################################
# Association Study
if len(pheno)==1:
    Target = PED.merge(Gene_temp,left_on='IND_ID',right_on='IND_ID',how='outer')
    Target = pd.DataFrame((Target >> drop(Target.IND_ID)),dtype='float')
    
    out_temp = pd.DataFrame(columns=['CHROM','GeneStart','GeneEnd','GeneName','TargetID',
                                     'R2','BETA','BETA_SE','T_STAT','PVALUE','N'])
    
    out_temp.to_csv(args.out_prefix+"/association_study_"+args.method+".txt",sep='\t',header=True,index=None,mode='a')

    pool = multiprocessing.Pool(args.thread)
    pool.map(thread_single,[num for num in range(len(TargetID))])
    pool.close()
    pool.join()
        
elif len(pheno)>1:
    res = pd.DataFrame()
    res['IND_ID'] = np.array(PED.IND_ID).ravel()
    
    for i in range(len(pheno)):
        res[pheno[i]] = np.array(sm.OLS(PED[pheno[i]],sm.add_constant(PED[cov])).fit().resid).ravel()
    
    out_temp = pd.DataFrame(columns=['CHROM','GeneStart','GeneEnd','GeneName','TargetID',
                                     'R2','F_STAT','F_PVALUE','N'])
    
    out_temp.to_csv(args.out_prefix+"/association_study_"+args.method+".txt",sep='\t',header=True,index=None,mode='a')
    
    Target = res.merge(Gene_temp,left_on='IND_ID',right_on='IND_ID',how='outer')
    Target = pd.DataFrame((Target >> drop(Target.IND_ID)),dtype='float')
    
    pool = multiprocessing.Pool(args.thread)
    pool.map(thread_multi,[num for num in range(len(TargetID))])
    pool.close()
    pool.join()












