#!/usr/bin/env python

###############################################################
# Import packages needed
import argparse
import multiprocessing
import operator
import os
import subprocess
import sys

from time import time

import pandas as pd
import numpy as np

###############################################################
# time calculation
start_time = time()

###############################################################
# parse input arguments
parser = argparse.ArgumentParser(description='Get LD')

parser.add_argument('--chr', type=str, dest='chrm', required=True, 
	help='chromosome number')
parser.add_argument('--format', type=str, dest='data_format', choices=['GT', 'DS'], default='GT', 
	help='data format of VCF genotype data (DS, GT [default])')
parser.add_argument('--genofile', type=str, dest='geno_path', required=True)
parser.add_argument('--genofile_type', type=str, choices=['vcf', 'dosage'], 
	help='filetype of genofile (vcf, dosages)')
parser.add_argument('--genome_block', type=str, dest='block_path', required=True)
parser.add_argument('--maf', type=float, default=0, 
	help='folded Minor Allele Frequency threshold; range from 0-0.5 (default: 0)')
parser.add_argument('--out_dir', type=str, default=os.getcwd())
parser.add_argument('--out_ld_file', type=str, default='')
parser.add_argument('--sampleID', type=str, dest='sampleid_path', required=True)
parser.add_argument('--thread', type=int, default=1)
parser.add_argument('--TIGAR_dir', type=str, help='tool directory', 
	default=os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

## new
parser.add_argument('--log_file', type=str, default='')
parser.add_argument('--out_prefix', type=str, default='')
parser.add_argument('--sub_dir', type=int, choices=[0, 1], default=1)

args = parser.parse_args()
sys.path.append(args.TIGAR_dir)
import TIGARutils as tg

#############################################################
# set output file names
if not args.out_prefix:
	args.out_prefix = 'CHR' + args.chrm

if not args.log_file:
	args.log_file = args.out_prefix + '_LD_log.txt'

if not args.out_ld_file:
	args.out_ld_file = args.out_prefix + '_reference_cov'

# sub-directory in out directory
if args.sub_dir:
	out_sub_dir = os.path.join(args.out_dir, 'RefLD')
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

###############################################################
# DEFINE, IMPORT FUNCTIONS
import TIGARutils as tg

# limit to 4 decimal places max, strip trailing 0s
def cov_fmt(x): return ('%.4f' % x).rstrip('0').rstrip('.')

# trim array by positionin matrix (length should be rownumber:total for each row);
# format each element in each row, join all together separated by comma
def cov_str(cov_lst): return [','.join([cov_fmt(x) for x in row]) for row in [cov_lst[i][i:len(cov_lst)] for i in range(len(cov_lst))]]

def out_block_path(num):
	return(out_sub_dir + '/' + args.out_ld_file + '_block_' + str(num) + '.txt')

# def out_block_path(num):
# 	return(out_sub_dir + '/' + args.out_ld_file + '_block_' + str(num) + '.txt.gz')


# def merge_tabix_output(out_file, n_blocks):
# 	try:
# 		out_dir = tg.get_abs_path(os.path.dirname(out_file))

# 		block_paths = '"(' + ' '.join([out_block_path(i) for i in range(n_blocks)]) + ')"'
# 		# block_paths = ' '.join([out_block_path(i) for i in range(n_blocks)])

# 		call_args = [
# 			TIGAR_dir + '/TWAS/merge_tabix_ld_output.sh', 
# 			out_file, block_paths]

# 		subprocess.check_call(
# 			call_args,
# 			cwd=out_dir,
# 			stdout=subprocess.DEVNULL)

# 	except subprocess.CalledProcessError as err:
# 		raise err

def merge_tabix_output(out_file, n_blocks):
	try:
		out_dir = tg.get_abs_path(os.path.dirname(out_file))

		block_path_pre = out_sub_dir + '/' + args.out_ld_file + '_block'

		call_args = [
			TIGAR_dir + '/TWAS/merge_tabix_ld_output.sh', 
			out_file, 
			block_path_pre,
			str(n_blocks - 1)]

		subprocess.check_call(
			call_args,
			cwd=out_dir,
			stdout=subprocess.DEVNULL)

	except subprocess.CalledProcessError as err:
		raise err


###############################################################
# check input arguments
if args.genofile_type == 'vcf':
	if (args.data_format != 'GT') and (args.data_format != 'DS'):
		raise SystemExit('Please specify the genotype data format used by the vcf file (--format ) as either "GT" or "DS".\n')
		
elif args.genofile_type == 'dosage':
	args.data_format = 'DS'

else:
	raise SystemExit('Please specify the type input genotype file type (--genofile_type) as either "vcf" or "dosage".\n')

out_refcovld_path = out_sub_dir + '/' + args.out_ld_file + '.txt'

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
	out_rc = out_refcovld_path))

# tg.print_args(args)

###############################################################
# Read in block information
print('Reading block annotation file.')

# read in block file
Blocks = pd.read_csv(
	args.block_path,
	sep='\t',
	usecols=['CHROM','Start','End'],
	dtype={'CHROM':object, 'Start':object, 'End':object})
Blocks = Blocks[Blocks['CHROM'] == args.chrm].reset_index(drop=True)
Blocks = tg.optimize_cols(Blocks)
n_blocks = len(Blocks)

# Startup for get LD job: get column header info, sampleIDs
sampleID, sample_size, geno_info = tg.sampleid_startup(**args.__dict__)

# write columns out to file
print('Creating file: ' + out_refcovld_path + '\n')
out_cols = ['#0','snpID','CHROM','POS','COV']
pd.DataFrame(columns=out_cols).to_csv(
	out_refcovld_path,
	sep='\t',
	index=None,
	header=True,
	mode='w')

print('********************************\n')

###############################################################
@tg.error_handler
def thread_process(num):
	Block = Blocks.loc[num]
	print('num=' + str(num))
	
	# read in and process genotype data; file must be bgzipped/tabix
	Geno = tg.read_tabix(Block.Start, Block.End, sampleID, **geno_info)

	# calculate, filter maf
	Geno = tg.calc_maf(Geno, sampleID, args.maf)

	# get upper covariance matrix
	mcovar = np.triu(np.cov(Geno[sampleID].values)).tolist()

	# output values
	Geno = Geno[['snpID', 'CHROM', 'POS']]
	Geno['COV'] = cov_str(mcovar)
	Geno.to_csv(
		out_block_path(num),
		sep='\t',
		index=None,
		header=None,
		# compression='gzip',
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

	merge_tabix_output(out_refcovld_path, n_blocks)


################################################################
# time calculation
elapsed_sec = time()-start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)

sys.stdout.close()
