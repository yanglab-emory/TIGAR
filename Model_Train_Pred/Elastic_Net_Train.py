#!/usr/bin/env python

#########################################################
# Import packages needed
import argparse
import io
from io import StringIO
import multiprocessing
import operator
import subprocess
import sys
from time import time
import warnings

import numpy as np
import pandas as pd

### import grid search for model selection
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold

### For Elastic Net Regression
from sklearn.linear_model import ElasticNet
from sklearn.metrics import r2_score

### For OLS regression in cross validation
from scipy import stats
import statsmodels.api as sm

warnings.filterwarnings("ignore")

######################################################
### time calculation
start_time = time()

####################################################
# Parse input variables
parser = argparse.ArgumentParser(description='Elastic Net Training')

### Specify tool directory
parser.add_argument('--TIGAR_dir',type=str)

### Gene Annotation and Expression level file
parser.add_argument('--gene_exp',type=str,dest='geneexp_path')

### Training sampleID
parser.add_argument('--train_sampleID',type=str,dest='sampleid_path')

### Specified chromosome number
parser.add_argument('--chr',type=str)

### Number of thread
parser.add_argument('--thread',type=int)

### Training genotype files
parser.add_argument('--genofile',type=str,dest='geno_path')

### Specified input file type(vcf or dosages)
parser.add_argument('--genofile_type',type=str)

### 'DS' or 'GT' for VCF genotype file
parser.add_argument('--format',type=str)

### for data selection
### Folded Minor Allele Frequency (range from 0-0.5)
parser.add_argument('--maf',type=float)

### p-value for Hardy Weinberg Equilibrium exact test
parser.add_argument('--hwe',type=float)

### window
parser.add_argument('--window',type=int)

### cvR2
parser.add_argument('--cvR2',type=int)

#########
# Specific for EN model training
### k-fold cross validation
parser.add_argument('--cv',type=int)

### Ratio of L1 and L2
parser.add_argument('--alpha',type=float)
#########

### output dir
parser.add_argument('--out_dir',type=str)

args = parser.parse_args()

sys.path.append(args.TIGAR_dir)

#######################################################
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
    trainX = train.iloc[:,0:-1]
    trainY = train.iloc[:,-1]

    if test is not None:
        testX = test.iloc[:,0:-1]
        testY = test.iloc[:,-1]

    else:
        testX = trainX
        testY = trainY

    clf = GridSearchCV(
        ElasticNet(
            l1_ratio=Alpha,
            fit_intercept=False),
        [{'alpha':np.arange(0,1.01,0.01)}],
        cv=k).fit(trainX,trainY)

    Lambda = clf.best_params_['alpha']

    reg = ElasticNet(
            l1_ratio=Alpha,
            alpha=Lambda).fit(trainX,trainY)

    predY = reg.predict(testX)

    lm = sm.OLS(testY, sm.add_constant(predY)).fit()

    if test is not None:
        return lm.rsquared

    beta = reg.coef_
    Rsquared = r2_score(trainY,predY)
    Pvalue = lm.f_pvalue
    cvm = clf.best_score_

    return beta, Rsquared, Pvalue, Lambda, cvm


# function to do the ith cross validation step
def do_cv(i, target_geno_exp_df, cv_trainID, cv_testID):
    train_geno_exp = target_geno_exp_df.loc[cv_trainID[i]].dropna()
    test_geno_exp = target_geno_exp_df.loc[cv_testID[i]].dropna()

    cv_rsquared = elastic_net(train_geno_exp, test_geno_exp)
    return cv_rsquared


###############################################################
# variable checking
### Print all variables for model training
print("********************************\n   Input Arguments\n********************************\n")

print("Gene Annotation and Expression data file: "+args.geneexp_path + "\n")
print("Training sampleID file: "+args.sampleid_path+ "\n")
print("Chrmosome: "+args.chr+ "\n")
print("Training genotype file: "+args.geno_path+ "\n")

if args.genofile_type=='vcf':
    print("VCF genotype file is used for training with genotype format: " + args.format + "\n")
    gcol_sampleids_strt_ind = 9
elif args.genofile_type=='dosage':
    print("Dosage genotype file is used for Training."+ "\n")
    args.format = 'DS'
    gcol_sampleids_strt_ind = 5    
else:
    raise SystemExit("Please specify input genotype file as either 'vcf' or 'dosage'."+ "\n")

print("Gene region size: window = "+str(args.window)+ "\n")

if args.cvR2:
    print("Evaluating Elastic-Net model by 5-fold cross validation.\n")

print("Threshold for MAF: "+str(args.maf)+ "\n")
print("Threshold for HWE p-value: "+str(args.hwe)+ "\n")

print("Using "+str(args.cv)+"-fold for cross-validation to tune penalty parameter: lambda"+ "\n")
print("The ratio for L1 & L2 penalty used by Elastic-Net regression: alpha = "+str(args.alpha)+ "\n")

print("Number of threads: "+str(args.thread)+ "\n")
print("Output dir: "+args.out_dir+ "\n")


out_train_weight_path = args.out_dir+'/CHR'+args.chr+'_EN_train_eQTLweights.txt'
print("Training weights output file: " + out_train_weight_path +"\n")

out_train_info_path = args.out_dir+'/CHR'+args.chr+'_EN_train_GeneInfo.txt'
print("Training info file: " + out_train_info_path +"\n")

print("********************************\n\n")
###############################################################

### Training Processing

### Read in Gene annotation and Expression level file (text file)
### First five columns should be fixed:
### 1.Chrmosome Number
### 2.GeneStart Posistion
### 3.GeneEnd Position
### 4.TargetID (i.e.GeneID, treated as unique annotation for each gene)
### 5.Gene Name

# Gene Expression header, sampleIDs
exp_cols = tg.get_header(args.geneexp_path)
exp_sampleids = exp_cols[5:]

# genofile header, sampleIDs
g_cols = tg.call_tabix_header(args.geno_path)
gcol_sampleids = g_cols[gcol_sampleids_strt_ind:]

# geno, exp overlapping sampleIDs
gcol_exp_sampleids = np.intersect1d(gcol_sampleids, exp_sampleids)

if not gcol_exp_sampleids.size:
    raise SystemExit("The gene expression file and genotype file have no sampleIDs in common.")

# if user specified path with sampleids, and at least one sampid from that file, get intersection; otherwise use the overlapping sampleIDs in the genotype and expression file
if args.sampleid_path:
    spec_sampleids = pd.read_csv(
        args.sampleid_path,
        sep='\t',
        header=None)

    spec_sampleids = spec_sampleids[0].drop_duplicates()

    sampleID = np.intersect1d(spec_sampleids, gcol_exp_sampleids)

else:
    sampleID = gcol_exp_sampleids

sample_size = sampleID.size

if not sample_size:
    raise SystemExit("There is no overlapped sample IDs between gene expression file, genotype file, and sampleID file.")

# get columns to read in
g_cols_ind, g_dtype = tg.genofile_cols_dtype(g_cols, args.genofile_type, sampleID)
e_cols_ind, e_dtype = tg.exp_cols_dtype(exp_cols, sampleID)

### EXPRESSION FILE
# ### Extract expression level by chromosome
GeneExp_chunks = pd.read_csv(
    args.geneexp_path, 
    sep='\t', 
    iterator=True, 
    chunksize=10000,
    usecols=e_cols_ind,
    dtype=e_dtype)

GeneExp = pd.concat([x[x['CHROM']==args.chr] for x in GeneExp_chunks]).reset_index(drop=True)

if GeneExp.empty:
    raise SystemExit('There are no valid gene expression training data for chromosome '+args.chr+".")

GeneExp = tg.optimize_cols(GeneExp)

TargetID = GeneExp.TargetID.values
n_targets = TargetID.size

### PREP CROSS VALIDATION SAMPLES
### Split sampleIDs for cross validation
if args.cvR2:
    print("Evaluate Elastic-Net model by average R2 of 5-fold cross validation ... "+ "\nSplit sample IDs randomly into 5 folds ..."+ "\n")

    kf = KFold(n_splits=5)
    kf_splits = [(sampleID[x], sampleID[y]) for x,y in kf.split(sampleID)]
    CV_trainID, CV_testID = zip(*kf_splits)

else:
    print("Skipping sample split for 5-fold cross validation ...")


### PREP OUTPUT
weight_out_cols = ['CHROM','POS','ID','REF','ALT','TargetID','MAF','p_HWE','ES']
pd.DataFrame(columns=weight_out_cols).to_csv(
    out_train_weight_path,
    header=True,
    index=None,
    sep='\t',
    mode='w')

info_out_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName','sample_size','n_snp','n_effect_snp','CVR2','TrainPVALUE','TrainR2','k-fold','alpha','Lambda','cvm']
pd.DataFrame(columns=info_out_cols).to_csv(
    out_train_info_path,
    header=True,
    index=None,
    sep='\t',
    mode='w')

###############################################################
### thread function
def thread_process(num):
    try:
        target = TargetID[num]
        print("\nnum="+str(num)+"\nTargetID="+target)
        target_exp = GeneExp.iloc[[num]]

        start = str(max(int(target_exp.GeneStart)-args.window,0))
        end = str(int(target_exp.GeneEnd)+args.window)

        ### select corresponding vcf file by tabix
        print("Loading genotype data for Gene: " + target +"\n")
        g_proc_out = tg.call_tabix(args.geno_path, args.chr, start, end)

        if not g_proc_out:
            print("There is no genotype data for this Gene:"+target+"\n") 
            return None

        ### Recode subprocess output in 'utf-8'
        target_geno = pd.read_csv(StringIO(g_proc_out.decode('utf-8')),
                sep='\t',
                low_memory=False,
                header=None,
                usecols=g_cols_ind,
                dtype=g_dtype)
        target_geno.columns = [g_cols[i] for i in target_geno.columns]
        target_geno = tg.optimize_cols(target_geno)

        # get snpIDs
        target_geno['snpID'] = tg.get_snpIDs(target_geno)
        target_geno = target_geno.drop_duplicates(['snpID'],keep='first').reset_index(drop=True)

        n_snp = target_geno['snpID'].size

        # prep vcf file
        if args.genofile_type=='vcf':
            target_geno = tg.check_prep_vcf(target_geno, args.format, sampleID)

        # reformat sample values
        target_geno[sampleID] = target_geno[sampleID].apply(lambda x:tg.reformat_sample_vals(x,args.format), axis=0)

        # get, filter maf
        target_geno = tg.calc_maf(target_geno, sampleID, args.maf, op=operator.ge)

        # get, filter p_HWE
        target_geno = tg.calc_p_hwe(target_geno, sampleID, args.hwe, op=operator.ge)

        # merge geno, expression files, transpose
        target_geno_exp = pd.concat([
            target_geno.set_index(['snpID'])[sampleID],
            target_exp.set_index(['TargetID'])[sampleID]
            ]).T

        # 5-FOLD CROSS-VALIDATION
        ### Evaluate whether Elastic Net model valid
        if args.cvR2:
            print("Running 5-fold CV to tune Elastic-Net penalty parameter for Gene: "+target+"\n")
            do_cv_args = [target_geno_exp, CV_trainID, CV_testID]
            k_fold_R2 = [do_cv(i, *do_cv_args) for i in range(5)]

            avg_r2_cv = sum(k_fold_R2)/5

            if avg_r2_cv < 0.005:
                print("Average R2 by 5-fold CV = " + str(avg_r2_cv) + ", less than 0.005 for "+target+"\nSkip running Elastic-Net regression model for gene " + target +"\n")
                return None

        else:
            avg_r2_cv = 0
            print("Skipping evaluation by 5-fold CV R2 ..." + "\n")


        # FINAL MODEL TRAINING
        print("Train Elastic-Net imputation model for Gene: " + target + "\n")
        # initialize target_weights dataframe
        target_weights = target_geno[['CHROM','POS','snpID','REF','ALT','p_HWE','MAF']].copy()

        target_weights['TargetID'] = target

        # do elastic net training
        target_weights['ES'],R2,Pvalue,Lambda,cvm = elastic_net(target_geno_exp)

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

        ### output training information, result from elastic net
        train_info = target_exp[['CHROM','GeneStart','GeneEnd','TargetID','GeneName']].copy()

        train_info['sample_size'] = sample_size
        train_info['n_snp'] = n_snp
        train_info['n_effect_snp'] = n_effect_snp = target_weights.ES.size
        train_info['CVR2'] = avg_r2_cv
        train_info['TrainPVALUE'] = Pvalue if not np.isnan(Pvalue) else 'NaN'
        train_info['TrainR2'] = R2 if n_effect_snp else 0
        train_info['k_fold'] = args.cv
        train_info['alpha'] = args.alpha
        train_info['lambda'] = Lambda
        train_info['cvm'] = cvm

        train_info.to_csv(
            out_train_info_path,
            header=None,
            index=None,
            sep='\t',
            mode='a')

    except Exception as e:
        e_type, e_obj, e_tracebk = sys.exc_info()
        e_line_num = e_tracebk.tb_lineno

        e, e_type, e_line_num = [str(x) for x in [e, e_type, e_line_num]]

        print('Caught a type '+ e_type +' exception for TargetID='+target+', num=' + str(num) + ' on line '+e_line_num+':\n' + e )

    finally:
        # print info to log do not wait for buffer to fill up
        sys.stdout.flush()

##################################################################
### start thread process
if __name__ == '__main__':
    print("Starting Elastic-Net training for "+str(n_targets)+" target genes.")
    pool = multiprocessing.Pool(args.thread)
    pool.map(thread_process,[num for num in range(n_targets)])
    pool.close()
    pool.join()
    print('Done.\n')

#########################################################
### time calculation
elapsed_sec = time()-start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print("Computation time (DD:HH:MM:SS): " + elapsed_time)

