#!/usr/bin/env python

###############################################################
# Import packages needed
import argparse
import multiprocessing
import operator
import subprocess
import sys

from io import StringIO
from time import time

import pandas as pd
import numpy as np

###############################################################
# time calculation
start_time = time()

###############################################################
# parse input arguments
parser = argparse.ArgumentParser(description='Get LD')

# Specify tool directory
parser.add_argument('--TIGAR_dir',type=str)

# chromosome block information path
parser.add_argument('--genome_block',type=str,dest='block_path')

# sampleID path
parser.add_argument('--sampleID',type=str,dest='sampleid_path')

# chromosome number
parser.add_argument('--chr',type=str)

# genotype file path
parser.add_argument('--genofile',type=str,dest='geno_path')

# specified input file type (vcf or doasges)
parser.add_argument('--genofile_type',type=str)

# 'DS' or 'GT'
parser.add_argument('--format',type=str)

# maf threshold for seleting genotype data to calculate covariance matrix
parser.add_argument('--maf',type=float)

# number of threads
parser.add_argument('--thread',type=int)

# output dir
parser.add_argument('--out_dir',type=str)

args = parser.parse_args()

sys.path.append(args.TIGAR_dir)

###############################################################
# DEFINE, IMPORT FUNCTIONS
import TIGARutils as tg

# limit to 4 decimal places max, strip trailing 0s
def cov_print_frmt(x): return ('%.4f' % x).rstrip('0').rstrip('.')

###############################################################
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
    
out_ref_cov_path = args.out_dir+'/CHR'+args.chr+'_reference_cov.txt'

###############################################################
# Print input arguments
print(
'''********************************
Input Arguments

Genome block annotation (based on LD structure) file: {block_path}

SampleID file: {sampleid_path}

Chromosome: {chr}

Reference genotype file: {geno_path}

Genotype file used is type: {genofile_type}

Genotype data format: {format}

MAF threshold for SNP inclusion: {maf}

Number of threads: {thread}

Output directory: {out_dir}

Output reference covariance results file: {out_rc}
********************************'''.format(
    **args.__dict__,
    out_rc = out_ref_cov_path))

###############################################################
# Read in block information
print('Reading block annotation file.')
chr_blocks = pd.read_csv(
    args.block_path,
    sep='\t',
    usecols=['CHROM','Start','End'],
    dtype={'CHROM':object,'Start':object,'End':object})
chr_blocks = chr_blocks[chr_blocks['CHROM']==args.chr].reset_index(drop=True)
chr_blocks = tg.optimize_cols(chr_blocks)

n_blocks = len(chr_blocks)

print('Reading genotype file header.\n')
g_cols = tg.call_tabix_header(args.geno_path)
gcol_sampleids = g_cols[gcol_sampleids_strt_ind:]

# get sampleIDs
print('Reading sampleID file.\n')
spec_sampleids = pd.read_csv(
    args.sampleid_path,
    sep='\t',
    header=None)[0].drop_duplicates()

print('Matching sampleIDs.\n')
sampleID = np.intersect1d(gcol_sampleids, spec_sampleids)

if not sampleID.size:
    raise SystemExit('There are no sampleID in both the genotype data and the specified sampleID file.')

# get the indices and dtypes for reading files into pandas
g_cols_ind, g_dtype = tg.genofile_cols_dtype(g_cols, args.genofile_type, sampleID)

# write columns out to file
print('Creating file: ' + out_ref_cov_path + '\n')
out_cols = ['#snpID', 'CHROM', 'POS', 'COV']
pd.DataFrame(columns=out_cols).to_csv(
    out_ref_cov_path,
    sep='\t',
    index=None,
    header=True,
    mode='w')

print('********************************\n')

###############################################################
@tg.error_handler
def thread_process(num):
    block = chr_blocks.loc[num]
    print('num=' + str(num))
    
    print('Reading genotype data.')
    g_proc_out = tg.call_tabix(args.geno_path, args.chr, block.Start, block.End)

    if not g_proc_out:
        print('There is no genotype data in this block.\n')
        return None

    block_geno = pd.read_csv(StringIO(g_proc_out.decode('utf-8')),
        sep='\t',
        low_memory=False,
        header=None,
        usecols=g_cols_ind,
        dtype=g_dtype)

    block_geno.columns = [g_cols[i] for i in block_geno.columns]
    block_geno = tg.optimize_cols(block_geno)

    # get snpIDs
    block_geno['snpID'] = tg.get_snpIDs(block_geno)
    block_geno = block_geno.drop_duplicates(['snpID'],keep='first').reset_index(drop=True)

    # prep vcf file
    if args.genofile_type=='vcf':
        block_geno = tg.check_prep_vcf(block_geno, args.format, sampleID)

    # reformat sample values
    block_geno[sampleID]=block_geno[sampleID].apply(lambda x:tg.reformat_sample_vals(x,args.format), axis=0)

    # calculate, filter maf
    block_geno = tg.calc_maf(block_geno, sampleID, args.maf)

    # get covariance matrix
    mcovar = np.cov(block_geno[sampleID])

    for i in range(len(block_geno)):
        snpinfo = block_geno.iloc[i][['snpID', 'CHROM', 'POS']]
        snpcov = ','.join([cov_print_frmt(x) for x in mcovar[i,i:]])
        covar_info = pd.DataFrame(np.append(snpinfo, snpcov)).T
        covar_info.to_csv(
            out_ref_cov_path,
            sep='\t',
            index=None,
            header=None,
            mode='a')

    print('Block LD calculation completed for block.\n')


##################################################################
# thread process

if __name__ == '__main__':
    print('Starting LD calculation for ' + str(n_blocks) + ' blocks.\n')
    pool = multiprocessing.Pool(args.thread)
    pool.imap(thread_process,[num for num in range(n_blocks)])
    pool.close()
    pool.join()
    print('Done.')

################################################################
# time calculation
elapsed_sec = time()-start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)

