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
parser.add_argument('--TIGAR_dir',type=str)

# for Gene annotation and Expression-level file path
parser.add_argument('--gene_exp',type=str,dest='geneexp_path')

# training genotype file sampleIDs path
parser.add_argument('--train_sampleID',type=str,dest='sampleid_path')

# specified chromosome number
parser.add_argument('--chr',type=str)

# Genotype file path
parser.add_argument('--genofile',type=str,dest='geno_path')

# specified input file type (vcf or dosages)
parser.add_argument('--genofile_type',type=str)

# format of genotype data 'DS' or 'GT'
parser.add_argument('--format',type=str)

# window
parser.add_argument('--window',type=int)

# missing rate: threshold for excluding SNPs with too many missing values
parser.add_argument('--missing_rate',type=float)

# maf
parser.add_argument('--maf',type=float)

# p-value for HW test
parser.add_argument('--hwe',type=float)

# cvR2
## 0 do not run cvR2
## 1 run cvR2 [default]
parser.add_argument('--cvR2',type=int)

# Bayesian inference algorithm used by DPR: 
## '1' (Variational Bayesian)
## '2' (MCMC)
parser.add_argument('--dpr',type=str)

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
parser.add_argument('--thread',type=int)

# output dir
parser.add_argument('--out_dir',type=str)

args = parser.parse_args()

sys.path.append(args.TIGAR_dir)

#############################################################
# DEFINE, IMPORT FUNCTIONS
import TIGARutils as tg

# preps dpr input files, runs DPR, reads in dpr output
def prep_call_dpr(bimbam_df, pheno_df, snpannot_df, dpr_file_dir, targetid):
    ## PATHS FOR DPR INPUT
    bimbam_pth = dpr_file_dir + targetid + '_bimbam.txt'
    pheno_pth = dpr_file_dir + targetid + '_pheno.txt'
    snpannot_pth = dpr_file_dir + targetid + '_snp_annot.txt'

    #  OUTPUT FILES FOR DPR INPUT
    bimbam_df.to_csv(
        bimbam_pth,
        header=False,
        index=None,
        sep='\t',
        mode='w',
        float_format='%f')

    pheno_df.to_csv(
        pheno_pth,
        header=False,
        index=None,
        sep='\t',
        mode='w',
        float_format='%f')

    snpannot_df.to_csv(
        snpannot_pth,
        header=False,
        index=None,
        sep='\t',
        mode='w',
        float_format='%f')

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
        dtype={'snpID': object, 'b': np.float64, 'beta': np.float64})

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
def do_cv(i, target, target_geno_df, target_exp_df, snp_annot_df, cv_trainID, cv_testID, ):
    dpr_file_dir_cv = abs_out_dir + '/CV_Files' + args.job_suf + '/'
    target_cv = target + '_CV' + str(i+1)

    trainID = cv_trainID[i]
    testID = cv_testID[i]

    bimbam_train = target_geno_df[np.concatenate((['snpID','REF','ALT'],trainID))]

    pheno_train = target_exp_df[trainID].T

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
    bimbam_test = target_geno_df[np.append(['snpID'],testID)]
    pheno_test = target_exp_df[testID].T

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
abs_out_dir = tg.get_abs_path(args.out_dir)

# check input arguments
if args.genofile_type == 'vcf':
    gcol_sampleids_strt_ind = 9

    if (args.format != 'GT') and (args.format != 'DS'):
        raise SystemExit('Please specify the genotype data format used by the vcf file (--format ) as either "GT" or "DS".\n')
        
elif args.genofile_type == 'dosage':
    gcol_sampleids_strt_ind = 5
    args.format = 'DS'

else:
    raise SystemExit('Please specify the type input genotype file type (--genofile_type) as either "vcf" or "dosage".\n')

out_train_weight_path = args.out_dir + '/temp_' + args.out_weight_file
out_train_info_path = args.out_dir + '/' +  args.out_info_file
#############################################################
# Print input arguments to log
print(
'''********************************
Input Arguments

Gene Annotation and Expression file: {geneexp_path}

Training sampleID file: {sampleid_path}

Chromosome: {chr}

Training genotype file: {geno_path}

Genotype file used for training is type: {genofile_type}

Genotype data format: {format}

Gene training region SNP inclusion window: +-{window}

Excluding SNPs if missing rate exceeds: {missing_rate}

MAF threshold for SNP inclusion: {maf}

HWE p-value threshold for SNP inclusion: {hwe}

{cvR2_str} DPR model by 5-fold cross validation.

DPR model: {dpr} - {dpr_type}

Output Effect-size type: {ES}

Number of threads: {thread}

Output directory: {out_dir}

Output training weights file: {out_weight}

Output training info file: {out_info}
********************************'''.format(
    **args.__dict__,
    dpr_type = {'1':'Variational Bayesian', '2':'MCMC'}[args.dpr],
    cvR2_str = {0:'Skipping evaluation of', 1:'Evaluating'}[args.cvR2],
    out_weight = out_train_weight_path,
    out_info = out_train_info_path))

#############################################################
# Prepare DPR input

### Read in Gene annotation and Expression level file (text file)
### First five columns should be fixed:
### 1.Chrmosome Number
### 2.GeneStart Posistion
### 3.GeneEnd Position
### 4.TargetID (i.e.GeneID, treated as unique annotation for each gene)
### 5.Gene Name

# gene Expression header, sampleIDs
print('Reading genotype, expression file headers.\n')
exp_cols = tg.get_header(args.geneexp_path)
exp_sampleids = exp_cols[5:]

# genofile header, sampleIDs
try:
    g_cols = tg.call_tabix_header(args.geno_path)
except: 
    g_cols = tg.get_header(args.geno_path, zipped=True)


gcol_sampleids = g_cols[gcol_sampleids_strt_ind:]

# geno, exp overlapping sampleIDs
gcol_exp_sampleids = np.intersect1d(gcol_sampleids, exp_sampleids)

if not gcol_exp_sampleids.size:
    raise SystemExit('The gene expression file and genotype file have no sampleIDs in common.')

# get sampleIDs
print('Reading sampleID file.\n')
spec_sampleids = pd.read_csv(
    args.sampleid_path,
    sep='\t',
    header=None)[0].drop_duplicates()

print('Matching sampleIDs.\n')
sampleID = np.intersect1d(spec_sampleids, gcol_exp_sampleids)

sample_size = sampleID.size

if not sample_size:
    raise SystemExit('There are no overlapped sample IDs between the gene expression file, genotype file, and sampleID file.')

# get columns to read in
g_cols_ind, g_dtype = tg.genofile_cols_dtype(g_cols, args.genofile_type, sampleID)
e_cols_ind, e_dtype = tg.exp_cols_dtype(exp_cols, sampleID)

# extract expression level for chromosome
print('Reading gene expression data.\n')
try:
    GeneExp_chunks = pd.read_csv(
        args.geneexp_path, 
        sep='\t', 
        iterator=True, 
        chunksize=10000,
        usecols=e_cols_ind,
        dtype=e_dtype)

    GeneExp = pd.concat([x[x['CHROM']==args.chr] for x in GeneExp_chunks]).reset_index(drop=True)

except:
    GeneExp_chunks = pd.read_csv(
        args.geneexp_path, 
        sep='\t', 
        iterator=True, 
        header=None,
        chunksize=10000,
        usecols=e_cols_ind)

    GeneExp = pd.concat([x[x[0]==args.chr] for x in GeneExp_chunks]).reset_index(drop=True).astype(e_dtype)

    GeneExp.columns = [exp_cols[i] for i in GeneExp.columns]

if GeneExp.empty:
    raise SystemExit('There are no valid gene expression training data for chromosome ' + args.chr + '\n')

GeneExp = tg.optimize_cols(GeneExp)

TargetID = GeneExp.TargetID
n_targets = TargetID.size

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
weight_out_cols = ['CHROM','POS', 'snpID', 'REF','ALT','TargetID','MAF','p_HWE','ES','b','beta']
pd.DataFrame(columns=weight_out_cols).to_csv(
    out_train_weight_path,
    sep='\t',
    header=True,
    index=None,
    mode='w')

print('Creating file: ' + out_train_info_path + '\n')
info_out_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName','sample_size','n_snp', 'n_effect_snp','CVR2','TrainPVALUE','TrainR2']
pd.DataFrame(columns=info_out_cols).to_csv(
    out_train_info_path,
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
    target_exp = GeneExp.iloc[[num]]

    start = str(max(int(target_exp.GeneStart) - args.window, 0))
    end = str(int(target_exp.GeneEnd) + args.window)

    # READ IN AND PROCESS GENOTYPE DATA 
    target_geno = read_genotype(args.geno_path, args.chr, start, end, g_cols, g_cols_ind, g_dtype, args.genofile_type, args.format, sampleID)
    # # Requirement for input vcf file: Must be bgzip and tabix
    # ### select corresponding vcf file by tabix
    # print('Reading genotype data.')
    # g_proc_out = tg.call_tabix(args.geno_path, args.chr, start, end)

    # if not g_proc_out:
    #     print('There is no genotype data for TargetID: ' + target + '\n')
    #     return None

    # print('Preparing DPR input.')
    # target_geno = pd.read_csv(StringIO(g_proc_out.decode('utf-8')),
    #         sep='\t',
    #         low_memory=False,
    #         header=None,
    #         usecols=g_cols_ind,
    #         dtype=g_dtype,
    #         na_values=['.', '.|.', './.'])
    # target_geno.columns = [g_cols[i] for i in target_geno.columns]
    # target_geno = tg.optimize_cols(target_geno)

    # # get snpIDs
    # target_geno['snpID'] = tg.get_snpIDs(target_geno)
    # target_geno = target_geno.drop_duplicates(['snpID'],keep='first').reset_index(drop=True)

    # # prep vcf file
    # if args.genofile_type == 'vcf':
    #     target_geno = tg.check_prep_vcf(target_geno, args.format, sampleID)

    # # reformat sample values
    # target_geno = tg.reformat_sample_vals(target_geno, args.format, sampleID)
    
    # filter out variants that exceed missing rate threshold
    target_geno = tg.handle_missing(target_geno, sampleID, args.missing_rate)
    
    # get, filter maf
    target_geno = tg.calc_maf(target_geno, sampleID, args.maf)

    # get, filter p_HWE
    target_geno = tg.calc_p_hwe(target_geno, sampleID, args.hwe)

    # target_geno, target_expr not centered since 
    ## 1) DPR script centers input 
    ## 2) DPR script sometimes segfaults when reading in genotype files when data was centered

    snp_annot = target_geno[['snpID','POS','CHROM']]

    # 5-FOLD CROSS-VALIDATION
    if args.cvR2:
        print('Running 5-fold CV.')
        do_cv_args = [target, target_geno, target_exp, snp_annot, CV_trainID, CV_testID]

        k_fold_R2 = [do_cv(i, *do_cv_args) for i in range(5)]

        avg_r2_cv = sum(k_fold_R2) / 5

        print('Average R2 for 5-fold CV: {:.4f}'.format(avg_r2_cv))

        if avg_r2_cv < 0.005:
            print('Average R2 < 0.005; Skipping DPR training for TargetID: ' + target + '\n')
            return None

    else:
        avg_r2_cv = 0
        print('Skipping evaluation by 5-fold CV average R2...')


    # FINAL MODEL TRAINING
    print('Running DPR training.')
    dpr_file_dir = abs_out_dir + '/DPR_Files' + args.job_suf + '/'
    
    bimbam = target_geno[np.concatenate((['snpID','REF','ALT'],sampleID))]

    pheno = target_exp[sampleID].T

    # PREP INPUT FILES, CALL DPR, READ IN DPR OUTPUT
    try:
        dpr_out = prep_call_dpr(bimbam, pheno, snp_annot, dpr_file_dir, target)

    except subprocess.CalledProcessError as err:
        print('DPR failed for TargetID: ' + target + '\n')
        return None

    # FILTER FOR SNPS WITH ES!=0
    n_snp = dpr_out.ES.size

    dpr_out = dpr_out[dpr_out.ES!=0]

    n_effect_snp = dpr_out.ES.size

    # R2 CALCULATION
    out_weights = dpr_out[['snpID','ES']]

    bimbam = bimbam.drop(columns=['REF','ALT'])

    train_pvalue, train_rsquared = calc_r2(out_weights, bimbam, pheno)

    # OUTPUT TARGET WEIGHTS TO FILE
    # initialize df with MAF, pHWE, other info
    target_weights = target_geno[['CHROM','POS','REF','ALT','snpID','p_HWE','MAF']].copy()
    target_weights = target_weights[target_weights.snpID.isin(dpr_out.snpID)]
    target_weights['TargetID'] = target

    # merge with dpr output weights, reorder columns using existing col list
    target_weights = target_weights.merge(
        dpr_out, 
        left_on='snpID',
        right_on='snpID',
        how='outer')[weight_out_cols]

    target_weights.to_csv(
        out_train_weight_path,
        sep='\t',
        header=None,
        index=None,
        mode='a')        

    # OUTPUT TARGET TRAINING INFO TO FILE
    # initialize dataframe for storing training info
    train_info = target_exp[
        ['CHROM','GeneStart','GeneEnd','TargetID','GeneName']].copy()
    train_info['sample_size'] = sample_size
    train_info['n_snp'] = n_snp
    train_info['n_effect_snp'] = n_effect_snp
    train_info['CVR2'] = avg_r2_cv
    train_info['TrainPVALUE'] = train_pvalue
    train_info['TrainR2'] = train_rsquared

    # output training info
    train_info.to_csv(
        out_train_info_path,
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


