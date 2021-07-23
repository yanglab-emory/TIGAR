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
parser.add_argument('--TIGAR_dir', type=str)

# chromosome block information path
parser.add_argument('--genome_block', type=str, dest='block_path')

# sampleID path
parser.add_argument('--sampleID', type=str, dest='sampleid_path')

# chromosome number
parser.add_argument('--chr', type=str, dest='chrm')

# genotype file path
parser.add_argument('--genofile', type=str, dest='geno_path')

# specified input file type (vcf or doasges)
parser.add_argument('--genofile_type', type=str)

# 'DS' or 'GT'
parser.add_argument('--format', type=str, dest='data_format')

# maf threshold for seleting genotype data to calculate covariance matrix
parser.add_argument('--maf', type=float)

# number of threads
parser.add_argument('--thread', type=int)

# output dir
parser.add_argument('--out_dir', type=str)

args = parser.parse_args()

sys.path.append(args.TIGAR_dir)

###############################################################
# DEFINE, IMPORT FUNCTIONS
import TIGARutils as tg

# limit to 4 decimal places max, strip trailing 0s
# def cov_print_frmt(x): return ('%.4f' % x).rstrip('0').rstrip('.')
def cov_fmt(x): return ('%.4f' % x).rstrip('0').rstrip('.')

# trim array by positionin matrix (length should be rownumber:total for each row);
# format each element in each row, join all together separated by comma
def cov_str(cov_lst): return [','.join([cov_fmt(x) for x in row]) for row in [cov_lst[i][i:len(cov_lst)] for i in range(len(cov_lst))]]

###############################################################
# check input arguments
if args.genofile_type == 'vcf':
	if (args.data_format != 'GT') and (args.data_format != 'DS'):
		raise SystemExit('Please specify the genotype data format used by the vcf file (--format ) as either "GT" or "DS".\n')
		
elif args.genofile_type == 'dosage':
	args.data_format = 'DS'

else:
	raise SystemExit('Please specify the type input genotype file type (--genofile_type) as either "vcf" or "dosage".\n')
	
out_ref_cov_path = args.out_dir + '/CHR' + args.chrm + '_reference_cov.txt'

###############################################################
# Print input arguments
print(
'''********************************
Input Arguments
Genome block annotation (based on LD structure) file: {block_path}
SampleID file: {sampleid_path}
Chromosome: {chrm}
Reference genotype file: {geno_path}
Genotype file used is type: {genofile_type}
Genotype data format: {data_format}
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
# read in block file
chr_blocks = pd.read_csv(
	args.block_path,
	sep='\t',
	usecols=['CHROM', 'Start', 'End'],
	dtype={'CHROM':object, 'Start':object, 'End':object})
chr_blocks = chr_blocks[chr_blocks['CHROM'] == args.chrm].reset_index(drop=True)
chr_blocks = tg.optimize_cols(chr_blocks)

n_blocks = len(chr_blocks)

# Startup for get LD job: get column header info, sampleIDs
sampleID, geno_cols_info = tg.genosampid_startup(**args.__dict__)

# write columns out to file
print('Creating file: ' + out_ref_cov_path + '\n')
out_cols = ['#snpID', 'CHROM', 'POS', 'block', 'blockpos',  'COV']
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
	
	# read in and process genotype data
	# file must be bgzipped and tabix
	block_geno = tg.read_genotype(block.Start, block.End, sampleID, **geno_cols_info, **args.__dict__)

	# calculate, filter maf
	block_geno = tg.calc_maf(block_geno, sampleID, args.maf)

	# get upper covariance matrix
	mcovar = np.triu(np.cov(block_geno[sampleID].values)).tolist()

	# output values
	block_geno = block_geno[['snpID', 'CHROM', 'POS']]
	block_geno['COV'] = cov_str(mcovar)
	block_geno.to_csv(
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
