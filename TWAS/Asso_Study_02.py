#!/usr/bin/env python

###################################################################
# Import packages needed
import argparse
import multiprocessing
import subprocess
import sys

from io import StringIO
from time import time

import numpy as np
import pandas as pd

from scipy.stats import chi2

###############################################################
# time calculation
start_time = time()

###############################################################
# parse input arguments
parser = argparse.ArgumentParser(description='Asso Study 02')

# Specify tool directory
parser.add_argument('--TIGAR_dir',type=str)

# Gene annotation file path
parser.add_argument('--gene_anno',type=str,dest='annot_path')

# chromosome number
parser.add_argument('--chr',type=str)

# Weight file path
parser.add_argument('--weight',type=str,dest='w_path')

# GWAS Z score file path
parser.add_argument('--Zscore',type=str,dest='z_path')

# Reference covariance file path
parser.add_argument('--LD',type=str,dest='ld_path')

# window
parser.add_argument('--window',type=float)

# Weight threshold to include SNP in TWAS
parser.add_argument('--weight_threshold',type=float)

# specify 'FUSION', 'SPrediXcan', or 'both': Zscore test statistic to use
parser.add_argument('--test_stat',type=str)

# Number of threads
parser.add_argument('--thread',type=int)

# Output dir
parser.add_argument('--out_dir',type=str)

args = parser.parse_args()

sys.path.append(args.TIGAR_dir)

###############################################################
## Import TIGAR functions, define other functions
import TIGARutils as tg

# Determine index of required columns,
# assigns correct dtype to correct index
# 'ID' column does not exist in TIGAR training output from previous versions
# check to see if it is in the file, if not will need to generate snpIDs later

# return correct snpID and Zscore value
# change sign of Zscore value if matching snpID is flipped wrt Weight snpID
def handle_flip(df: pd.DataFrame, origID, flipID, origValCol, orig_overlap, flip_overlap):
    orig = df[origID].values
    flip = df[flipID].values
    origval = df[origValCol].values

    ids = np.empty_like(orig)
    val = np.empty_like(origval)

    for i in range(len(df)):
        if orig[i] in orig_overlap:
            ids[i], val[i] = orig[i], origval[i]
        elif flip[i] in flip_overlap:
            ids[i], val[i] = flip[i], -origval[i]

    return ids, val

def get_pval(z): return np.format_float_scientific(1-chi2.cdf(z**2, 1), precision=15, exp_digits=0)

def get_V_cor(V_cov):
    V_cov = V_cov.copy()
    v = np.sqrt(np.diag(V_cov))
    outer_v = np.outer(v, v)
    V_cor = V_cov / outer_v
    V_cor[V_cov == 0] = 0
    return V_cor

# def get_z_denom(V, ZW):
#     return np.sqrt(np.linalg.multi_dot([ZW.ES.values, V_cov, ZW.ES.values]))

# def get_spred_zscore(V_cov, ZW, snp_sd):
#     ZW = ZW.copy()
#     zscore = snp_sd.dot(ZW.ES.values * ZW.Zscore.values) / get_z_denom(V_cov, ZW)
#     return zscore, get_pval(zscore)
    
# def get_fusion_zscore(V_cov, ZW, snp_sd=None):
#     ZW = ZW.copy()
#     V_cor = get_V_cor(V_cov)
#     zscore = np.vdot(ZW.Zscore.values, ZW.ES.values) / get_z_denom(V_cor, ZW)
#     return zscore, get_pval(zscore)

def get_z_denom(V, w):
    return np.sqrt(np.linalg.multi_dot([w, V, w]))

def get_spred_zscore(V_cov, w, Z_gwas, snp_sd):
    Z_twas = snp_sd.dot(w * Z_gwas) / get_z_denom(V_cov, w)
    return Z_twas, get_pval(Z_twas)
    
def get_fusion_zscore(V_cov, w, Z_gwas, snp_sd=None):
    V_cor = get_V_cor(V_cov)
    Z_twas = np.vdot(Z_gwas, w) / get_z_denom(V_cor, w)
    return Z_twas, get_pval(Z_twas)

def get_burden_zscore(test_stat, get_zscore_args):
    if test_stat =='FUSION':
        return get_fusion_zscore(*get_zscore_args)
    if test_stat == 'SPrediXcan':
        return get_spred_zscore(*get_zscore_args)

#############################################################
# Print input arguments to log
out_twas_path = args.out_dir + '/CHR' + args.chr + '_sumstat_assoc.txt'

print(
'''********************************
Input Arguments

Gene annotation file specifying genes for TWAS: {annot_path}

Chromosome: {chr}

cis-eQTL weight file: {w_path}

GWAS summary statistics Z-score file: {z_path}

Reference LD genotype covariance file: {ld_path}

Gene training region SNP inclusion window: +-{window}

SNP weight inclusion threshold: {weight_threshold}

Test statistic to use: {test_stat_str}

Number of threads: {thread}

Output directory: {out_dir}

Output TWAS results file: {out_path}
********************************'''.format(
    **args.__dict__,
    test_stat_str = 'FUSION and SPrediXcan' if args.test_stat=='both' else args.test_stat,
    out_path = out_twas_path))

###############################################################
### Read in gene annotation 
print('Reading gene annotation file.')
Gene_chunks = pd.read_csv(
    args.annot_path, 
    sep='\t', 
    iterator=True, 
    chunksize=10000,
    dtype={'CHROM':object,'GeneStart':np.int64,'GeneEnd':np.int64,'TargetID':object,'GeneName':object}, 
    usecols=['CHROM','GeneStart','GeneEnd','TargetID','GeneName'])

Gene = pd.concat([x[x['CHROM']==args.chr] for x in Gene_chunks] ).reset_index(drop=True)

Gene = tg.optimize_cols(Gene)

TargetID = np.array(Gene.TargetID)
n_targets = TargetID.size

# read in headers for Weight and Zscore files
print('Reading file headers.\n')
w_cols = tg.get_header(args.w_path, zipped=True)
z_cols = tg.get_header(args.z_path, zipped=True)

# get the indices and dtypes for reading files into pandas
w_cols_ind, w_dtype = tg.weight_cols_dtype(w_cols)
z_cols_ind, z_dtype = tg.zscore_cols_dtype(z_cols)

# PREP OUTPUT - print output headers to files
print('Creating file: ' + out_twas_path + '\n')
out_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName','n_snps']
if args.test_stat == 'both':
    out_cols = out_cols + ['FUSION_Z','FUSION_PVAL','SPred_Z','SPred_PVAL']
else:
    out_cols = out_cols + ['Zscore','PVALUE']

pd.DataFrame(columns=out_cols).to_csv(
    out_twas_path,
    sep='\t',
    index=None,
    header=True,
    mode='w')

###############################################################
# thread function
@tg.error_handler
def thread_process(num):
    target = TargetID[num]
    print('num=' + str(num) + '\nTargetID=' + target)
    Gene_info = Gene.iloc[[num]].reset_index(drop=True)

    # get start and end positions to tabix
    start = str(max(int(Gene_info.GeneStart)-args.window,0))
    end = str(int(Gene_info.GeneEnd)+args.window)

    # tabix Weight file
    # print('Reading weight data.')
    w_proc_out = tg.call_tabix(args.w_path, args.chr, start, end)

    if not w_proc_out:
        print('No test SNPs with non-zero cis-eQTL weights for TargetID: ' + target + '\n')
        return None

    # tabix Zscore file
    # print('Reading Zscore data.')
    z_proc_out = tg.call_tabix(args.z_path, args.chr, start, end)

    if not z_proc_out:
        print('No test SNPs with GWAS Zscore for TargetID: ' + target + '\n')
        return None

    # parse tabix output for Weight, filtered by target, threshold
    Weight_chunks = pd.read_csv(
        StringIO(w_proc_out.decode('utf-8')),
        sep='\t',
        header=None,
        low_memory=False,
        iterator=True, 
        chunksize=10000,
        usecols=w_cols_ind,
        dtype=w_dtype)

    Weight = pd.concat([x[ (x[w_cols_ind[4]]==target) & (abs(x[w_cols_ind[5]]) > args.weight_threshold )] for x in Weight_chunks]).reset_index(drop=True)

    if Weight.empty:
        print('No test SNPs with cis-eQTL weights with magnitude that exceeds specified weight threshold for TargetID: ' + target + '\n')
        return None

    Weight.columns = [w_cols[i] for i in tuple(Weight.columns)]
    Weight = tg.optimize_cols(Weight)

    # 'ID' snpID column does not exist in TIGAR training output from previous versions
    # check to see if it is in the file, if not will need to generate snpIDs later
    if 'ID' in Weight.columns:
        Weight = Weight.rename(columns={'ID':'snpID'})

    if not 'snpID' in Weight.columns:
        Weight['snpID'] = tg.get_snpIDs(Weight)

    # parse tabix output for Zscore
    Zscore = pd.read_csv(
        StringIO(z_proc_out.decode('utf-8')),
        sep='\t',
        header=None,
        low_memory=False,
        usecols=z_cols_ind,
        dtype=z_dtype)

    Zscore.columns = [z_cols[i] for i in tuple(Zscore.columns)]
    Zscore = tg.optimize_cols(Zscore)

    Zscore['IDorig'] = tg.get_snpIDs(Zscore)
    Zscore['IDflip'] = tg.get_snpIDs(Zscore, flip=True)

    # check for overlapping SNPs in Weight, Zscore data
    snp_overlap_orig = np.intersect1d(Weight.snpID, Zscore.IDorig)
    snp_overlap_flip = np.intersect1d(Weight.snpID, Zscore.IDflip)
    snp_overlap = np.concatenate((snp_overlap_orig, snp_overlap_flip))

    if not snp_overlap.size:
        print('No overlapping test SNPs that have magnitude of cis-eQTL weights greater than threshold value and with GWAS Zscore for TargetID: ' + target + '\n')
        return None

    # filter dataframes by overlapping SNPs
    Weight = Weight[Weight.snpID.isin(snp_overlap)]
    Zscore = Zscore[Zscore.IDorig.isin(snp_overlap_orig) | Zscore.IDflip.isin(snp_overlap_flip)]
    Zscore = Zscore.drop(columns=['CHROM','POS','REF','ALT'])

    # if handle any flipped matches in Zscore using Weight snpIDs as reference
    if (snp_overlap_orig.size > 0) and (snp_overlap_flip.size > 0):
        Zscore['snpID'], Zscore['Zscore'] = handle_flip(Zscore,'IDorig','IDflip','Zscore',snp_overlap_orig, snp_overlap_flip)
    elif snp_overlap_orig.size == snp_overlap.size:   
        Zscore['snpID'], Zscore['Zscore'] = Zscore['IDorig'], Zscore['Zscore']
    else:
        Zscore['snpID'], Zscore['Zscore'] = Zscore['IDflip'], -Zscore['Zscore']

    # merge Zscore and Weight dataframes on snpIDs
    ZW = Weight.merge(Zscore[['snpID','Zscore']], 
        left_on='snpID', 
        right_on='snpID').drop_duplicates(['snpID'], keep='first').reset_index(drop=True)

    snp_search_ids = ZW.snpID.values

    # Read in reference covariance matrix file by snpID
    MCOV = tg.get_ld_data(args.ld_path, snp_search_ids)

    if MCOV.empty:
      print('No reference covariance information for target SNPs for TargetID: ' + target + '\n')
      return None

    # get the snp variance and covariance matrix
    snp_sd, V_cov = tg.get_ld_matrix(MCOV)

    ZW = ZW[ZW.snpID.isin(MCOV.snpID)]
    n_snps = str(ZW.snpID.size)

    print('Running TWAS.\nN SNPs=' + n_snps)

    ### create output dataframe
    result = Gene_info.copy()
    result['n_snps'] = n_snps

    ### calculate zscore(s), pvalue(s)
    # get_zscore_args = [V_cov, ZW, snp_sd]
    get_zscore_args = [V_cov, ZW.ES.values, ZW.Zscore.values, snp_sd]

    if args.test_stat == 'both':
        result['FUSION_Z'], result['FUSION_PVAL'] = get_fusion_zscore(*get_zscore_args)

        result['SPred_Z'], result['SPred_PVAL'] = get_spred_zscore(*get_zscore_args)

    else:
        result['TWAS_Zscore'], result['PVALUE'] = get_burden_zscore(args.test_stat, get_zscore_args)

    # write to file
    result.to_csv(
        out_twas_path,
        sep='\t',
        index=None,
        header=None,
        mode='a')

    print('Target TWAS completed.\n')

###############################################################
# thread process
if __name__ == '__main__':
    print('Starting TWAS for ' + str(n_targets) + ' target genes.\n')
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







