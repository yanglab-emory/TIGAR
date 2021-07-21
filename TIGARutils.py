#!/usr/bin/env python

#########################################################
import functools
import operator
import os
import re
import subprocess
import sys
import traceback

from io import StringIO
from itertools import groupby

import pandas as pd
import numpy as np
#########################################################
## FUNCTIONS:

# error_handler

# calc_maf
# call_tabix
# call_tabix_header
# format_elapsed_time
# get_header

# exp_cols_dtype
# genofile_cols_dtype
# weight_cols_dtype
# zscore_cols_dtype
# gwas_cols_dtype
# MCOV_cols_dtype

## Handling covariance files:
# ld_cols
# get_ld_regions_list
# call_tabix_regions
# get_ld_regions_data
# get_ld_data
# get_ld_matrix

# get_snpIDs
# optimize_cols
# reformat_sample_vals
# reformat_vcf
# check_prep_vcf
# substr_in_strarray

#########################################################

# wrapper for thread_process functions; adds error catching/logging for failed targets
def error_handler(func):
	@functools.wraps(func)
	def wrapper(num, *args, **kwargs):     
		try:
			return func(num, *args, **kwargs)

		except Exception as e:
			e_info = sys.exc_info()
			e_type = e_info[0].__name__
			
			# don't print traceback info for wrapper
			e_tracebk = ''.join(traceback.format_tb(e_info[2])[1:])

			print('Caught an exception for num={}:\n  {}: {}\nTraceback:\n{}'.format(num, e_type, e, e_tracebk))

		finally:
			sys.stdout.flush()

	return wrapper


# returns absolute path
def get_abs_path(x): return os.path.abspath(os.path.expanduser(os.path.expandvars(x)))


# wrapper for genotype functions; adds error handling for when an empty dataframe is read in during concatenation
def empty_df_handler(func):
	@functools.wraps(func)
	def wrapper(num, *args, **kwargs):     
		try:
			return func(num, *args, **kwargs)

		except (pd.errors.EmptyDataError, pd.errors.ParserError, AttributeError) as e:
			e_info = sys.exc_info()
			e_type = e_info[0].__name__
			
			# don't print traceback info for wrapper
			e_tracebk = ''.join(traceback.format_tb(e_info[2])[1:])
			return None

	return wrapper




# train startup
def train_startup(geno_path, genofile_type, geneexp_path, sampleid_path, **kwargs):

	print('Reading file headers.\n')

	# expression file header
	exp_cols = get_header(geneexp_path)
	exp_sampleids = exp_cols[5:]

	# genotype file header
	try:
		geno_cols = call_tabix_header(geno_path)
	except: 
		geno_cols = get_header(geno_path, zipped=True)

	# get sampleids in the genotype file
	geno_sampleids_strt_ind = {'dosage': 5, 'vcf': 9}[genofile_type]
	geno_sampleids = geno_cols[geno_sampleids_strt_ind:]

	## read in sampleids file
	print('Reading sampleID file.\n')
	spec_sampleids = pd.read_csv(
		sampleid_path,
		sep='\t',
		header=None)[0].drop_duplicates()

	## match sampleids between files
	print('Matching sampleIDs.\n')
	sampleID = functools.reduce(np.intersect1d, (exp_sampleids, geno_sampleids, spec_sampleids))

	sample_size = sampleID.size

	if not sample_size:
		raise SystemExit('There are no overlapped sample IDs between the gene expression file, genotype file, and sampleID file.')

	## columns to read-in
	exp_cols_info = exp_cols_dtype(exp_cols, sampleID)
	geno_cols_info = genofile_cols_dtype(geno_cols, genofile_type, sampleID)

	return sampleID, sample_size, exp_cols_info, geno_cols_info


def read_genexp(geneexp_path, cols, col_inds, dtype, chrm, **kwargs):
	print('Reading gene expression data.\n')

	try:
		GeneExp_chunks = pd.read_csv(
			geneexp_path, 
			sep='\t', 
			iterator=True, 
			chunksize=10000,
			usecols=col_inds,
			dtype=dtype)

		GeneExp = pd.concat([x[x['CHROM']==chrm] for x in GeneExp_chunks]).reset_index(drop=True)

	except:
		GeneExp_chunks = pd.read_csv(
			geneexp_path, 
			sep='\t', 
			iterator=True, 
			header=None,
			chunksize=10000,
			usecols=col_inds)

		GeneExp = pd.concat([x[x[0]==chrm] for x in GeneExp_chunks]).reset_index(drop=True).astype(dtype)

		GeneExp.columns = [exp_cols[i] for i in GeneExp.columns]

	if GeneExp.empty:
		raise SystemExit('There are no valid gene expression training data for chromosome ' + chrm + '\n')

	GeneExp = optimize_cols(GeneExp)

	TargetID = GeneExp.TargetID
	n_targets = TargetID.size

	return GeneExp, TargetID, n_targets


# divide range of positions into intervals of size n or less, convert to tabixable string
def get_gt_regions_list(chrm, start, end, n):
	start = int(start)
	end = int(end)
	n_regions = ((end - start) // n) + 1
	for i in range(n_regions):
		yield chrm + ':' + str(start + (n * i) + i) + '-' + str(min(start + (n * (i + 1)) + i, end))


# get proc_out from function and parse data for regions
@empty_df_handler
def get_gt_regions_data(regs_str, path, g_cols, g_cols_ind, g_dtype):

	proc_out = call_tabix_regions(path, regs_str)

	try:
		regs_data = pd.read_csv(
			StringIO(proc_out.decode('utf-8')),
			sep='\t',
			low_memory=False,
			header=None,
			usecols=g_cols_ind,
			dtype=g_dtype)
	except pd.errors.ParserError as e:
		regs_chunks = pd.read_csv(
			StringIO(proc_out.decode('utf-8')),
			sep='\t',
			low_memory=False,
			header=None,
			iterator=True,
			chunksize=5000,
			usecols=g_cols_ind,
			dtype=g_dtype)
		regs_data = pd.concat([chunk for chunk in regs_chunks]).reset_index(drop=True)

	regs_data.columns = [g_cols[i] for i in regs_data.columns]
	regs_data = optimize_cols(regs_data)

	# get snpIDs
	regs_data['snpID'] = get_snpIDs(regs_data)
	regs_data = regs_data.drop_duplicates(['snpID'],keep='first').reset_index(drop=True)

	return regs_data

# get genotype data in format needed
@empty_df_handler
def prep_gt_regions_data(regs_data, genofile_type, data_format, sampleID): 
	# prep vcf file
	if genofile_type == 'vcf':
		regs_data = check_prep_vcf(regs_data, data_format, sampleID)

	# reformat sample values
	regs_data = reformat_sample_vals(regs_data, data_format, sampleID)
	
	return regs_data

# get proc_out from function and parse data for regions, then prep
@empty_df_handler
def get_prep_gt_regions_data(regs_str, path, g_cols, g_cols_ind, g_dtype, genofile_type, data_format, sampleID):

	regs_data = get_gt_regions_data(regs_str, path, g_cols, g_cols_ind, g_dtype)

	regs_data = prep_gt_regions_data(regs_data, genofile_type, data_format, sampleID)

	return regs_data

# # read in and prep genotype data
# def read_genotype(path, chrm, start, end, g_cols, g_cols_ind, g_dtype, genofile_type, data_format, sampleID):

# 	regs_args = [path, g_cols, g_cols_ind, g_dtype]
# 	prep_args = [genofile_type, data_format, sampleID]

# 	try:
# 		regs_str = chrm + ':' + start + '-' + end
# 		gt_data = get_prep_gt_regions_data(regs_str, *regs_args, *prep_args)
# 	except (MemoryError, pd.errors.ParserError):
# 		# data may be too large; if so try subset instead of getting all SNPs at once
# 		n = 25000
# 		# *** DECREASE THIS
# 		while n:
# 			try: 
# 				regs_str_lst = get_gt_regions_list(chrm, start, end, n)
# 				gt_data = pd.concat([get_prep_gt_regions_data(regs_str, *regs_args, *prep_args) for regs_str in regs_str_lst])
# 			except (MemoryError, pd.errors.ParserError):
# 				if (n > 10000):
# 					n -= 10000
# 				elif (n > 3000):
# 					n -= 1000
# 				else:
# 					n -= 500
# 				pass
# 			else:
# 				n = 0
# 	return gt_data


def filter_vcf_line(line: bytes, data_format, col_inds, split_multi = True):
	# split line into list
	row = line.split(b'\t')

	# convert data_format to byte
	bformat = str.encode(data_format)

	# get index of data format
	data_ind = row[8].split(b':').index(bformat)
	row[8] = bformat

	# may be multiallelic; only first alt allele will be used unless split_multi; later data with bad values will be filtered out
	alt_alleles = row[4].split(b',')
	row[4] = alt_alleles[0]

	# filter sample columns to include only data in desired format; sampleIDs start at file column index 9; in new row sampleIDs now start at 4, ALT is column 3
	row = [row[x] if x <= 8 else row[x].split(b':')[data_ind] for x in col_inds]

	# turn multi-allelic lines into multiple biallelic lines
	if split_multi & (len(alt_alleles) > 1):
		sample_str = b'\t'.join(row[4:])
		line = bytearray()
		for j in range(1, len(alt_alleles) + 1):
			str_j = sample_str
			for k in range(1, len(alt_alleles) + 1):
				# substitute alt alleles besides the jth one with missing
				str_j = re.sub(str(k).encode(), b'.', str_j) if (k != j) else str_j
			# set jth alt allele to 1
			str_j = re.sub(str(j).encode(), b'1', str_j)
			# join row info information
			line_j = b'\t'.join([*row[0:3], alt_alleles[j-1], str_j])
			# append linebreak if needed
			line_j = line_j if line_j.endswith(b'\n') else line_j + b'\n'
			line += line_j

	else:
		# row to bytestring
		line = b'\t'.join(row)

		# append linebreak if needed
		line = line if line.endswith(b'\n') else line + b'\n'

	return line

def read_genotype(start, end, sampleID, geno_cols_info, chrm, geno_path, genofile_type, data_format, **kwargs):

	# subprocess command
	command_str = ' '.join(['tabix', path, chrm + ':' + start + '-' + end])

	proc = subprocess.Popen(
		[command_str],
		shell=True,
		stdout=subprocess.PIPE,
		bufsize=1)

	# initialize bytearray
	proc_out = bytearray()

	# while subprocesses running, read lines into byte array
	while proc.poll() is None:
		line = proc.stdout.readline()
		if len(line) == 0:
			break
		if genofile_type == 'vcf':
			line = filter_vcf_line(line, data_format, col_inds)
		proc_out += line

	# read in lines still remaining after subprocess completes
	for line in proc.stdout:
		if genofile_type == 'vcf':
			line = filter_vcf_line(line, data_format, col_inds)
		proc_out += line

	# read data into dataframe
	geno_data = pd.read_csv(
		StringIO(proc_out.decode('utf-8')),
		sep='\t',
		low_memory=False,
		header=None,
		names=geno_cols_info['cols'],
		dtype=geno_cols_info['dtype'])

	geno_data = optimize_cols(geno_data)

	# get snpID, filter out duplicates
	geno_data['snpID'] = get_snpIDs(geno_data)
	geno_data = geno_data.drop_duplicates(['snpID'], keep='first').reset_index(drop=True)

	# remove non-valid GT values; ie those from multiallelic rows
	if data_format == 'GT':
		valid_GT = ['.|.', '0|0', '0|1', '1|0', '1|1', 
		'./.', '0/0', '0/1', '1/0', '1/1']
		geno_data = geno_data[np.all(geno_data[sampleID].isin(valid_GT), axis=1)].reset_index(drop=True)

	return geno_data


# Call tabix, read in lines into byte array
def call_tabix(path, chrm, start, end, add_command_str = ''):

	regs_str = chrm + ':' + start + '-' + end

	command_str = ' '.join(['tabix', path, regs_str, add_command_str])
	# ["tabix " + path + " " + chrm + ":" + start + "-" + end]

	proc = subprocess.Popen(
		[command_str],
		shell=True,
		stdout=subprocess.PIPE,
		bufsize=1)

	proc_out = bytearray()

	# while subprocesses running, read lines into byte array
	while proc.poll() is None:
		line =  proc.stdout.readline()
		if len(line) == 0:
			break
		proc_out += line

	# read in lines still remaining after subprocess completes
	for line in proc.stdout:
		proc_out += line

	return proc_out


# Call tabix to get header, read in lines into byte array;
# method for reading in headers for vcf files
def call_tabix_header(path, out='tuple', rename={}):
	# create dictionary for renaming columns
	# the first column in the vcf file is expected to be parsed as #CHROM;  automatically add that to the rename list
	rename = {**{'#CHROM':'CHROM'}, **rename}

	proc = subprocess.Popen(
		["tabix -H "+path],
		shell=True,
		stdout=subprocess.PIPE,
		bufsize=1)

	proc_out = bytearray()

	# while subprocesses running, read lines into byte array
	while proc.poll() is None:
		line =  proc.stdout.readline()
		if len(line) == 0:
			break
		# filter out lines starting with ##, which denotes comment in vcf file
		if not line.startswith(b"##"):
			proc_out += line

	# read in lines still remaining after subprocess completes
	for line in proc.stdout:
		if not line.startswith(b"##"):
			proc_out += line

	# decode bytes, use pandas to read in, rename columns from dictionary
	header = pd.read_csv(
		StringIO(proc_out.decode('utf-8')),
		sep='\t',
		error_bad_lines=False).rename(columns=rename)

	if out=='tuple':
		return tuple(header)

	elif out=='list':
		return list(header)

	return header


# return human readable elapsed time string 
def format_elapsed_time(time_secs):
	val = abs(int(time_secs))

	day = val // (3600*24)
	hour = val % (3600*24) // 3600
	mins = val % 3600 // 60
	secs = val % 60

	res = '%02d:%02d:%02d:%02d' % (day, hour, mins, secs)

	if int(time_secs) < 0:
		res = "-%s" % res

	return res


# get header from non-vcf file
def get_header(path, out='tuple', zipped=False, rename={}):

	# zipped files assumed to be bgzipped, tabixed files
	compress_type = 'gzip' if zipped else None

	# the first column in the file is expected to be parsed as #CHROM;  automatically add that to the rename list
	rename = {**{'#CHROM':'CHROM'}, **rename}

	# decode bytes, use pandas to read in, rename columns from dictionary
	header = pd.read_csv(
		path,
		sep='\t',
		header=0,
		compression=compress_type,
		low_memory=False,
		nrows=0).rename(columns=rename)

	if out=='tuple':
		return tuple(header)

	elif out=='list':
		return list(header)

	return header


# for testing on systems without tabix; not currently used by any scripts
def get_vcf_header(path, out='tuple'):
	
	proc = subprocess.Popen(
		["zgrep -m1 -E 'CHROM' "+path],
		shell=True,
		stdout=subprocess.PIPE,
		bufsize=1)

	proc_out = bytearray()
	while proc.poll() is None:
		line =  proc.stdout.readline()
		if len(line) == 0:
			break
		proc_out += line

	for line in proc.stdout:
		proc_out += line

	header = pd.read_csv(
		StringIO(proc_out.decode('utf-8')),
		sep='\t',
		error_bad_lines=False).rename(columns={'#CHROM':'CHROM'})   

	if out=='tuple':
		return tuple(header)

	elif out=='list':
		return list(header)

	return header

# determine indices of file cols to read in, dtype of each col
def get_cols_dtype(file_cols, cols, sampleid=None, genofile_type=None, add_cols=[], drop_cols=[], get_id=False, ind_namekey=False, ret_dict=False, **kwargs):

	# sampleid handling
	if sampleid is not None:
		cols = cols + sampleid.tolist()
		sampleid_dtype = object if genofile_type == 'vcf' else np.float64
		sampleid_dict = {x:sampleid_dtype for x in sampleid}
	else:
		sampleid_dict = {}

	dtype_dict = {
		'ALT': object,
		'b': np.float64,
		'beta': np.float64,
		'BETA': np.float64,
		'CHROM': object,
		'COV': object,
		'ES': np.float64,
		'FILTER': object,
		'FORMAT': object,
		'INFO': object,
		'GeneEnd': np.int64,
		'GeneName': object,
		'GeneStart': np.int64,
		'ID': object,
		'MAF': np.float64,
		'POS': np.int64,
		'QUAL': object,
		'REF': object,
		'SE': np.float64,
		'snpID': object,
		'TargetID': object,
		'Zscore': np.float64,
		**sampleid_dict}
	
	# cols set up
	cols = cols + add_cols

	# if type == 'vcf':
	#     cols.insert(4, 'snpID')

	if get_id:
		if ('snpID' in file_cols):
			cols.append('snpID')

		elif ('ID' in file_cols):
			cols.append('ID')

	cols = [x for x in cols if (x not in drop_cols)]

	# create output
	col_inds = tuple(sorted([file_cols.index(x) for x in cols]))

	if ind_namekey:
		# return dtype dict with keys as column names
		out_dtype_dict = dtype_dict
	else:
		# return dtype dict with keys as column index
		ind_dtype_dict = {file_cols.index(x):dtype_dict[x] for x in cols}
		out_dtype_dict = ind_dtype_dict

	if ret_dict:
		return {
			'file_cols': file_cols,
			'cols': [file_cols[i] for i in col_inds], 
			'col_inds': col_inds, 
			'dtype': out_dtype_dict}
	
	return file_cols_ind, out_dtype_dict

def exp_cols_dtype(file_cols, sampleid, **kwargs):
	return get_cols_dtype(file_cols, 
		cols=['CHROM', 'GeneStart', 'GeneEnd', 'TargetID', 'GeneName'], 
		sampleid=sampleid, ret_dict=True)

def genofile_cols_dtype(file_cols, genofile_type, sampleid, ind_namekey=True, **kwargs):
	return get_cols_dtype(file_cols, 
		cols=['CHROM','POS','REF','ALT'],
		sampleid=sampleid, genofile_type=genofile_type, 
		ind_namekey=ind_namekey, ret_dict=True)

def weight_cols_dtype(file_cols, add_cols=[], drop_cols=[], get_id=True, **kwargs):
	return get_cols_dtype(file_cols, 
		cols=['CHROM','POS','REF','ALT','TargetID','ES'], 
		add_cols=add_cols, drop_cols=drop_cols, 
		get_id=get_id)

def zscore_cols_dtype(file_cols, **kwargs):
	return get_cols_dtype(file_cols, 
		cols=['CHROM','POS','REF','ALT','Zscore'])

def gwas_cols_dtype(file_cols, **kwargs):
	return get_cols_dtype(file_cols, 
		cols=['CHROM','POS','REF','ALT','BETA','SE'])

def MCOV_cols_dtype(file_cols, add_cols=[], drop_cols=[], get_id=True, **kwargs):
	return get_cols_dtype(file_cols, 
		cols=['CHROM','POS','REF','ALT','COV'], 
		add_cols=add_cols, drop_cols=drop_cols, 
		get_id=get_id)


# get header of ld file, get indices of columns to read in
def get_ld_cols(path):
	# get header
	file_cols = tuple(pd.read_csv(
		path,
		sep='\t',
		header=0,
		compression='gzip',
		low_memory=False,
		nrows=0).rename(columns={'#snpID':'snpID', '#ID':'snpID', 'ID':'snpID', '#0':'row', '#row':'row', '0':'row'}))

	# ld files have gone through a lot of format revisions hence the possible need to rename
	cols = ['row', 'snpID', 'COV']

	file_cols_ind = tuple([file_cols.index(x) for x in cols])

	return file_cols, file_cols_ind


# yields formatted tabix regions strings
def get_ld_regions_list(snp_ids):

	# 'chrm:' prefix for region string 
	chrm = snp_ids[0].split(':')[0] + ':'

	# snp pos values as integers
	pos_vals = [int(snp.split(':')[1]) for snp in snp_ids]

	# get intervals of start,end positions; convert to tabix string
	for x, y in groupby(enumerate(pos_vals), lambda p: p[1]-p[0]):
		y = list(y)

		# chrm:start-end
		yield chrm + str(y[0][1]) + '-' + str(y[-1][1])


# call tabix using regions string
def call_tabix_regions(path, regs_str):

	proc = subprocess.Popen(
		['tabix '+path+' '+regs_str],
		shell=True,
		stdout=subprocess.PIPE,
		bufsize=1)
	proc_out = bytearray()

	# process while subprocesses running
	while proc.poll() is None:
		line =  proc.stdout.readline()
		if len(line) == 0:
			break
		proc_out += line

	# leftover lines
	for line in proc.stdout:
		proc_out += line

	return proc_out


# get proc_out from function and parse data for regions
def get_ld_regions_data(regs_str, path, snp_ids, ld_cols, ld_cols_ind):

	proc_out = call_tabix_regions(path, regs_str)

	regs_data = pd.read_csv(
		StringIO(proc_out.decode('utf-8')),
		sep='\t',
		low_memory=False,
		header=None,
		names=ld_cols,
		usecols=ld_cols_ind, 
		dtype={
			'snpID': object,
			'row': np.int32,
			'COV': object}
		).drop_duplicates(['snpID'], keep='first')

	regs_data = regs_data[regs_data.snpID.isin(snp_ids)]

	return regs_data


# read in covariance data for snps
def get_ld_data(path, snp_ids):

	# get columns names, indices for ld file
	ld_cols, ld_cols_ind = get_ld_cols(path)

	# format tabix regions from snp_ids; 'chrm:start-end'
	regs_lst = list(get_ld_regions_list(snp_ids))
	N = len(regs_lst)

	# arguments to pass
	regs_args = [path, snp_ids, ld_cols, ld_cols_ind]
	try:
		regs_str = ' '.join(regs_lst)
		cov_data = get_ld_regions_data(regs_str, *regs_args)

	except OSError:
		# argument may be too long for OS; if so try subset instead of getting all regions at once
		# print('Subseting regions to tabix.')
		n = 2500
		while n:
			try: 
				regs_str_lst = [' '.join(regs_lst[i:i+n]) for i in range(0, N, n)]
				cov_data = pd.concat([get_ld_regions_data(regs_str, *regs_args) for regs_str in regs_str_lst])
			except OSError:
				n -= 500
				pass
			else:
				n = 0

	return cov_data.set_index('row')


def get_ld_matrix(MCOV):
	MCOV = MCOV.copy()
	
	MCOV['COV'] =  MCOV['COV'].apply(lambda x:np.fromstring(x, dtype=np.float32, sep=','))

	inds = MCOV.index
	n_inds = inds.size
	V_upper = np.zeros((n_inds, n_inds))
	
	for i in range(n_inds):
		cov_i = MCOV.COV.loc[inds[i]]
		N = cov_i.size
		
		for j in range(i,n_inds):
			if inds[j] - inds[i] < N:
				V_upper[i,j] = cov_i[inds[j]-inds[i]]
			else:
				V_upper[i,j] = 0

	snp_Var = V_upper.diagonal()              
	V = V_upper + V_upper.T - np.diag(snp_Var)
	snp_sd = np.sqrt(snp_Var)

	return snp_sd, V


# return snp ids; join CHROM, POS, REF, ALT columns into : separated string
def get_snpIDs(df: pd.DataFrame, flip=False):
	chrom = df['CHROM'].astype('str').values
	pos = df['POS'].astype('str').values
	ref = df['REF'].values
	alt = df['ALT'].values
	if flip:
		return [':'.join(i) for i in zip(chrom,pos,alt,ref)]
	else:
		return [':'.join(i) for i in zip(chrom,pos,ref,alt)]

# Decrease memory by downcasting 'CHROM' column to integer, integer and float columns to minimum size that will not lose info
def optimize_cols(df: pd.DataFrame):

	# if 'CHROM' not convertable to integer, assume it's 'X', 'Y', 'MT', etc.
	if 'CHROM' in df.columns:
		try: 
			df['CHROM'] = df['CHROM'].astype(str).astype(np.int8)
		except ValueError:
			pass

	ints = df.select_dtypes(include=['int64']).columns.tolist()
	df[ints] = df[ints].apply(pd.to_numeric, downcast='integer')

	floats = df.select_dtypes(include=['float64']).columns.tolist()
	df[floats] = df[floats].apply(pd.to_numeric, downcast='float')

	return df




# Reform vcf file
### input each sample genotype
### For GT data_format:
###  code '0|0' or '0/0' as 0
###  code ('0|1' or '1|0')  or ('0/1' or '1/0') as 1
###  code '1|1' or '1/1' as 2
###  code '.|.' or './.' as nan(missing)

### For DS data_format:
### code '.' as nan(missing)
def reformat_sample_vals(df: pd.DataFrame, data_format, sampleID):
	df = df.reset_index(drop=True).copy()
	vals = df[sampleID].values
	if data_format=='GT':
		vals[(vals=='0|0')|(vals=='0/0')] = 0
		vals[(vals=='1|0')|(vals=='1/0')|(vals=='0|1')|(vals=='0/1')] = 1
		vals[(vals=='1|1')|(vals=='1/1')] = 2
		vals[(vals=='.|.')|(vals=='./.')] = np.nan
	elif data_format=='DS':
		vals[(vals=='.')] = np.nan
	vals = vals.astype(np.float32)
	df = pd.concat([df.drop(columns=sampleID), pd.DataFrame(vals, columns=sampleID)], axis=1)
	return df


# # reformats a vcf dataframe
# def reformat_vcf(df: pd.DataFrame, data_format, sampleID, uniqfrmts, singleformat=True):
# 	# df = df.copy()
# 	if singleformat:
# 		frmt_ind = uniqfrmts[0].split(':').index(data_format)
# 		df[sampleID]=df[sampleID].applymap(lambda x: x.split(':')[frmt_ind])
# 	else:
# 		# reformats sample values in row to include only specified format
# 		def vals_by_format(row):
# 			frmt_ind = row.FORMAT.split(':').index(data_format)
# 			return row[sampleID].apply(lambda y: y.split(':')[frmt_ind])

# 		# apply to each row        
# 		df[sampleID] = df.apply(lambda x: vals_by_format(x), axis=1)

# 	return df

# reformats sample values in row to include only specified format

@empty_df_handler
def check_prep_vcf(df: pd.DataFrame, data_format, sampleID):
	# df = df.copy()

	# check that all rows include data in the args.format format
	rowfrmts = np.unique(df.FORMAT.values)
	frmt_in_all = np.all(substr_in_strarray(data_format,rowfrmts))

	if not frmt_in_all:
		raise Exception("Exception in check_prep_vcf(): Specified genotype format, format=" + data_format + ", does not exist in all rows of the FORMAT column for this section of the input VCF file.")

	if rowfrmts.size > 1:
		def vals_by_format(row):
			frmt_ind = row.FORMAT.split(':').index(data_format)
			return row[sampleID].apply(lambda y: y.split(':')[frmt_ind])

		df[sampleID] = df.apply(lambda x: vals_by_format(x), axis=1)
		
	# else assume rowfrmts.size == 1 
	# if contains ':', needs to be reformatted
	elif (':' in rowfrmts[0]): 
		frmt_ind = rowfrmts[0].split(':').index(data_format)
		df[sampleID] = df[sampleID].applymap(lambda x: x.split(':')[frmt_ind])

	# if doesnt contain ':' but isn't equivalent to data_format then something's very wrong
	elif (rowfrmts[0] != data_format):
		raise Exception('Exception in check_prep_vcf(): There is only one format in the FORMAT column for this section of the input VCF file, format_in_vcf='+str(rowfrmts[0]) +', which contains the specified genotype format, format=' + data_format + ', but it cannot be parsed. ')

	df = df.drop(columns=['FORMAT'])

	# handle multi-allelic
	if data_format == 'GT':
		df = handle_multi_allele(df, sampleID)

	return df

# TIGAR currently only works with bi-allelic GT data
@empty_df_handler
def handle_multi_allele(df: pd.DataFrame, sampleID):
	df = df.copy()

	valid_GT = ['.|.', '0|0', '0|1', '1|0', '1|1', 
		'./.', '0/0', '0/1', '1/0', '1/1']

	# remove rows where any value contains an alt alleles; remaining values should be biallelic or missing
	# n_snps = df.shape[0]

	drop_index = df.loc[~np.all(df[sampleID].isin(valid_GT), axis=1)].index
	df = df.drop(index=drop_index)
	df = df.reset_index(drop=True)

	# df = df[np.all(df[sampleID].isin(valid_GT), axis=1)].reset_index(drop=True)

	# n_snps_removed = n_snps - df.shape[0]

	# if n_snps_removed:
	#     print('Multi-allelic GT calls: Removed {} variants with alternate call >1.'.format(n_snps_removed))

	# reformat ALT and snpIDs for remaining rows to contain only first alt allele
	df[['snpID','ALT']] = df[['snpID','ALT']].applymap(lambda x: x.split(',')[0])

	return df


# returns a boolean array; whether substring is in a np object array
def substr_in_strarray(substr, strarray):
   return np.frompyfunc(lambda x: substr in x, 1,1)(strarray)


# count the number of non-nan values
def count_notnan(x):
	return x.size - np.count_nonzero(np.isnan(x))


# drop variants with missing rate that exceeds threshold
def handle_missing(df: pd.DataFrame, sampleID, missing_rate, filter=True, op=operator.le):
	df = df.copy()

	# if all sample data for a row is NaN, drop the row
	drop_index = df.loc[df[sampleID].count(axis=1) == 0].index
	df = df.drop(index=drop_index)

	# calculate missing rate for each snp
	df['missing_rate'] = df[sampleID].apply(lambda x: np.count_nonzero(np.isnan(x))/x.size, axis=1)

	# downcast floats
	samp_miss_cols = np.append(sampleID,'missing_rate')
	df[samp_miss_cols] = df[samp_miss_cols].apply(pd.to_numeric, downcast='float')

	if filter:
		df = df[op(df.missing_rate, missing_rate)].reset_index(drop=True)

	return df


# calculate maf
def calc_maf(df: pd.DataFrame, sampleID, maf, filter=True, op=operator.gt):
	df = df.copy()

	# calculate MAF
	df['MAF'] = df[sampleID].apply(lambda x:np.nansum(x)/(2 * count_notnan(x)), axis=1)

	### Dealing with NaN - impute missing with mean
	samp_maf_cols = np.append(sampleID,'MAF')
	df[samp_maf_cols] = df[samp_maf_cols].apply(lambda x: x.fillna(2*x.MAF), axis=1)

	# downcast floats
	df[samp_maf_cols] = df[samp_maf_cols].apply(pd.to_numeric, downcast='float')

	if filter:
		df = df[op(df.MAF, maf)].reset_index(drop=True)

	return df


# Calculate, filter p val of HWE
def calc_p_hwe(df: pd.DataFrame, sampleID, pval, filter=True, op=operator.gt):
	df = df.copy()

	df['p_HWE'] = df[sampleID].apply(lambda x:p_HWE(x.dropna()), axis=1)
	df['p_HWE'] = pd.to_numeric(df['p_HWE'], downcast='float')

	if filter:
		df = df[op(df.p_HWE, pval)].reset_index(drop=True)

	return df


### Prepare for HWE input
def p_HWE(sample_row):
	vals = sample_row.values
	if not vals.size:
		p_hwe = np.nan
	else:
		N_hets = vals[(vals >= 0.5) & (vals < 1.5)].size
		N_aa = vals[(vals >= 0) & (vals < 0.5)].size
		N_AA = vals[(vals >= 1.5) & (vals <= 2)].size

		p_hwe = calc_HWE(N_hets, N_AA, N_aa)

	return p_hwe


###### Calculating p-value for Hardy Weinberg Equilibrium exact test

### gij denote number of minor alleles for ith SNP and jth sample
### 0 <= gij< 0.5 denote as 0
### 0.5 <= gij < 1.5 denote as 1
### 1.5 <= gij <2 denote as 2

### Input value:
### 1.obs_hets: Observed heterozygosity = Number of 1 in each SNPs(i.e. 0.5 <= gij < 1.5)
### 2.obs_hom1: Observed AA homozygosity = Number of 0 in each SNPs(i.e. 0 <= gij< 0.5)
### 3.obs_hom2: Observed aa homozygosity = Number of 2 in each SNPs(i.e. 1.5 <= gij <= 2)

### Output: p-value for Hardy Weinberg Equilibrium exact test
def calc_HWE(obs_hets, obs_hom1, obs_hom2):
	if obs_hom1 < 0 or obs_hom2 < 0 or obs_hets < 0:
		raise Exception("FATAL ERROR - SNP-HWE: Current genotype configuration (%s  %s %s) includes negative count" % (obs_hets, obs_hom1, obs_hom2))

	obs_homc = obs_hom2 if obs_hom1 < obs_hom2 else obs_hom1
	obs_homr = obs_hom1 if obs_hom1 < obs_hom2 else obs_hom2

	rare_copies = 2*obs_homr + obs_hets
	genotypes   = obs_hets + obs_homc + obs_homr

	het_probs = [0.0]*(rare_copies + 1)

	#start at midpoint
	mid = int(rare_copies*(2*genotypes - rare_copies)/(2*genotypes))

	#check to ensure that midpoint and rare alleles have same parity
	if (rare_copies & 1)^(mid & 1):
		mid += 1

	curr_hets = mid
	curr_homr = (rare_copies - mid) / 2
	curr_homc = genotypes - curr_hets - curr_homr

	het_probs[mid] = 1.0
	sum_het_probs = float(het_probs[mid])

	for curr_hets in range(mid,1,-2):
		het_probs[curr_hets - 2] = het_probs[curr_hets]*curr_hets*(curr_hets - 1.0)/(4.0*(curr_homr + 1.0)*(curr_homc + 1.0))

		sum_het_probs += het_probs[curr_hets - 2];

		# 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
		curr_homr += 1
		curr_homc += 1

	curr_hets = mid
	curr_homr = (rare_copies - mid)/2
	curr_homc = genotypes - curr_hets - curr_homr

	for curr_hets in range(mid,rare_copies-1,2):
		het_probs[curr_hets + 2] = het_probs[curr_hets]*4.0*curr_homr*curr_homc/((curr_hets + 2.0)*(curr_hets + 1.0))

		sum_het_probs += het_probs[curr_hets + 2]

		#add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
		curr_homr -= 1
		curr_homc -= 1

	for i in range(0,rare_copies + 1):
		het_probs[i] /= sum_het_probs

	#alternate p-value calculation for p_hi/p_lo
	p_hi = float(het_probs[obs_hets])
	for i in range(obs_hets,rare_copies+1):
		p_hi += het_probs[i]

	p_lo = float(het_probs[obs_hets])
	for i in range(obs_hets-1,-1,-1):
		p_lo += het_probs[i]

	p_hi_lo = 2.0 * p_hi if p_hi < p_lo else 2.0 * p_lo

	p_hwe = 0.0
	#  p-value calculation for p_hwe
	for i in range(0,rare_copies + 1):
		if het_probs[i] > het_probs[obs_hets]:
			continue
		p_hwe += het_probs[i]

	p_hwe = 1.0 if p_hwe > 1.0 else p_hwe

	return p_hwe

# center
def center(df: pd.DataFrame, sampleID):
	df = df.copy()
	df[sampleID] = df[sampleID].apply(lambda x: x - np.mean(x), axis=1)
	return df

# print args (for testing)
def print_args(args):
	for key, value in args.__dict__.items():
		if isinstance(value, str):
			print('args.', key, ' = \'', value, '\'', sep='')
		else:
			print('args.', key, ' = ', value, sep='')
	print('\n')
