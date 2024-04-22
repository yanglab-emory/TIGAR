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

# used to catch exceptions that dont require a traceback
class NoTargetDataError(Exception):
	pass

# wrapper for thread_process functions; adds error catching/logging for failed targets
def error_handler(func):
	@functools.wraps(func)
	def wrapper(num, *args, **kwargs):     
		try:
			return func(num, *args, **kwargs)

		except NoTargetDataError:
			return None

		except Exception as e:
			e_info = sys.exc_info()
			e_type = e_info[0].__name__
			
			# don't print traceback info for wrapper
			e_tracebk = ''.join(traceback.format_tb(e_info[2])[1:])

			print('Caught an exception for num={}:\n  {}: {}\nTraceback:\n{}'.format(num, e_type, e, e_tracebk))

		finally:
			sys.stdout.flush()

	return wrapper

def fatal_error_handler(func):
	@functools.wraps(func)
	def wrapper(num, *args, **kwargs):     
		try:
			return func(num, *args, **kwargs)

		except NoTargetDataError:
			return None

		except Exception as e:
			e_info = sys.exc_info()
			e_type = e_info[0].__name__
			
			# don't print traceback info for wrapper
			e_tracebk = ''.join(traceback.format_tb(e_info[2])[1:])

			print('Caught an exception for num={}:\n  {}: {}\nTraceback:\n{}'.format(num, e_type, e, e_tracebk))

			sys.exit()

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
	

def ped_startup(ped_path, pedinfo_path):
	Asso_Info = pd.read_csv(
		pedinfo_path,
		sep='\t',
		header=None,
		names=['Ind','Var'])
	pheno = Asso_Info[Asso_Info.Ind=='P'].Var.values
	n_pheno = pheno.size
	cov = Asso_Info[Asso_Info.Ind=='C'].Var.values

	if not n_pheno:
		raise SystemExit('No phenotypes provided by --PED_info.')
	if not cov.size:
		raise SystemExit('No covariates is provided by --PED_info.')

	pheno = pheno.tolist()
	cov = cov.tolist()

	# get header of PED file
	ped_header = get_header(ped_path)
	ped_cols = np.intersect1d(ped_header, [*pheno, *cov])

	# check that all phenos, covs in ped_header
	missing_cols = [x for x in [*pheno, *cov] if not x in ped_cols]
	if len(missing_cols):
		print('WARNING: The following phenotypes and/or covariates given in PED_info are missing in the PED file: ' + ', '.join(missing_cols) + '\n')
		pheno = [x for x in pheno if not x in missing_cols]
		cov = [x for x in cov if not x in missing_cols]

	ped_sampleid = pd.read_csv(
		ped_path,
		sep='\t',
		usecols=['IND_ID'])['IND_ID'].drop_duplicates()

	# ped_info = {'path':ped_path, 'file_cols':ped_header, 'cols':ped_cols, 'dtype': {x:object for x in ped_cols}, 'col_inds':tuple(sorted([ped_header.index(x) for x in ped_cols]))}

	return ped_sampleid, ped_cols, n_pheno, pheno, cov


# def sampleid_startup(geno_path, sampleid_path, chrm, genofile_type, data_format, geneexp_path=None, **kwargs):
def sampleid_startup(chrm=None, genofile_type=None, data_format=None, sampleid_path=0, geno_path=0, geneexp_path=0, ped_path=0, pedinfo_path=0, **kwargs):

	sampleid_lst = []

	## PED files
	if (ped_path and pedinfo_path):
		ped_sampleid, ped_cols, n_pheno, pheno, cov = ped_startup(ped_path, pedinfo_path)
		sampleid_lst.append(ped_sampleid)

	## read in sampleids file
	if sampleid_path:
		print('Reading sampleID file.\n')
		spec_sampleids = pd.read_csv(
			sampleid_path,
			sep='\t',
			header=None)[0].drop_duplicates()
		sampleid_lst.append(spec_sampleids)

	## read in file headers
	print('Reading file headers.\n')

	## genotype file
	if geno_path:
		try:
			geno_cols = call_tabix_header(geno_path)
		except: 
			geno_cols = get_header(geno_path, zipped=True)
		# get sampleids in the genotype file
		geno_sampleids_strt_ind = {'dosage': 5, 'vcf': 9}[genofile_type]
		geno_sampleids = geno_cols[geno_sampleids_strt_ind:]
		sampleid_lst.append(geno_sampleids)

	## expression file
	if geneexp_path:
		exp_cols = get_header(geneexp_path)
		exp_sampleids = exp_cols[5:]
		sampleid_lst.append(exp_sampleids)

	## match sampleids between files
	print('Matching sampleIDs.\n')
	sampleID = np.array(functools.reduce(np.intersect1d, tuple(sampleid_lst)))

	if not sampleID.size:
		raise SystemExit('There are no overlapped sample IDs in the input files.')

	print('Running job for ' + str(sampleID.size) + ' matched sampleIDs.')

	## return values
	return_lst = [sampleID, sampleID.size]

	if geno_path:
		geno_info = {'path': geno_path, 'chrm': chrm, 'genofile_type': genofile_type, 'data_format': data_format, 
			**get_cols_dtype(geno_cols, sampleid=sampleID,
				cols=['CHROM','POS','REF','ALT'],
				genofile_type=genofile_type, 
				ind_namekey=True)}
		return_lst.append(geno_info)

	if geneexp_path:
		exp_info = {'geneexp_path': geneexp_path, 'chrm': chrm, 
			**get_cols_dtype(exp_cols, sampleid=sampleID, 
				cols=['CHROM','GeneStart','GeneEnd','TargetID','GeneName'])}
		return_lst.append(exp_info)

	if (ped_path and pedinfo_path):
		return_lst = [*return_lst, ped_cols, n_pheno, pheno, cov]

	return tuple(return_lst)


def gwas_file_info(gwas_path, chrm, **kwargs):
	info_dict = get_cols_dtype(
		get_header(gwas_path, zipped=True), 
		cols=['CHROM','POS','REF','ALT','BETA','SE'], 
		ind_namekey=True)
	return {'path': gwas_path, 
		'chrm': chrm, 
		'sampleID': [], 
		'data_format': 'gwas', 
		'genofile_type': 'gwas', 
		**info_dict}

def bgw_weight_file_info(w_path, chrm, weight_threshold=0, add_cols=[], drop_cols=['ID'], **kwargs):
	info_dict = get_cols_dtype(
		get_header(w_path, zipped=True), 
		cols=['CHROM','POS','REF','ALT','Trans','PCP','beta'], 
		add_cols=add_cols, 
		drop_cols=drop_cols, 
		get_id=True,  
		ind_namekey=True)
	return {'path': w_path, 
		'chrm': chrm, 
		'sampleID': [], 
		'target_ind': info_dict['file_cols'].index('Trans'), 
		'data_format': 'bgw_weight', 
		'genofile_type': 'bgw_weight', 
		'weight_threshold': weight_threshold,
		**info_dict}

def weight_file_info(w_path, chrm, weight_threshold=0, add_cols=[], drop_cols=['ID'], **kwargs):
	info_dict = get_cols_dtype(
		get_header(w_path, zipped=True), 
		cols=['CHROM','POS','REF','ALT','TargetID','ES'], 
		add_cols=add_cols, 
		drop_cols=drop_cols, 
		get_id=True,  
		ind_namekey=True)
	return {'path': w_path, 
		'chrm': chrm, 
		'sampleID': [], 
		'target_ind': info_dict['file_cols'].index('TargetID'), 
		'data_format': 'weight', 
		'genofile_type': 'weight', 
		'weight_threshold': weight_threshold,
		**info_dict}

def zscore_file_info(z_path, chrm, **kwargs):
	info_dict = get_cols_dtype(
		get_header(z_path, zipped=True),
		cols=['CHROM','POS','REF','ALT','Zscore'],
		ind_namekey=True)
	return {'path':z_path, 
		'chrm': chrm, 
		'sampleID': [], 
		'data_format': 'zscore', 
		'genofile_type': 'zscore', 
		**info_dict}


# def bgw_weight_cols_dtype(file_cols, add_cols=[], drop_cols=[], get_id=True, ret_dict=True, ind_namekey=True, **kwargs):
# 	return get_cols_dtype(file_cols, 
# 		cols=['CHROM','POS','REF','ALT','Trans','PCP','beta'], 
# 		add_cols=add_cols, drop_cols=drop_cols, 
# 		get_id=get_id, ret_dict=ret_dict, ind_namekey=ind_namekey)

# read in annotation file or gene expression file
def read_gene_annot_exp(chrm=None, geneexp_path=None, annot_path=None, cols=['CHROM','GeneStart','GeneEnd','TargetID','GeneName'], col_inds=[0,1,2,3,4], dtype={'CHROM':object,'GeneStart':np.int64,'GeneEnd':np.int64,'TargetID':object,'GeneName':object}, **kwargs):

	# set the path
	path = geneexp_path if geneexp_path is not None else annot_path

	if (chrm is not None):
		try:
			Gene_chunks = pd.read_csv(
				path, 
				sep='\t', 
				iterator=True, 
				chunksize=10000,
				usecols=cols,
				dtype=dtype)
			Gene = pd.concat([x[x['CHROM']==chrm] for x in Gene_chunks]).reset_index(drop=True)
		except:
			Gene_chunks = pd.read_csv(
				path, 
				sep='\t', 
				iterator=True, 
				chunksize=10000,
				usecols=cols)
			Gene = pd.concat([x[x['CHROM']==chrm] for x in Gene_chunks]).reset_index(drop=True).astype(dtype)
			# Gene.columns = [cols[i] for i in Gene.columns]

		if Gene.empty:
			raise SystemExit('There are no valid gene expression/annotation data for chromosome ' + chrm + '.\n')
	else: 
		try:
			Gene = pd.read_csv(
				path, 
				sep='\t', 
				usecols=cols,
				dtype=dtype)
		except:
			Gene = pd.read_csv(
				path, 
				sep='\t', 
				header=None,
				usecols=cols).reset_index(drop=True).astype(dtype)
			# Gene.columns = [cols[i] for i in Gene.columns]

	Gene = optimize_cols(Gene)

	TargetID = Gene.TargetID
	n_targets = TargetID.size

	return Gene, TargetID, n_targets


## line filter functions for read_tabix

def filter_vcf_line(line: bytes, bformat, col_inds, split_multi_GT):
	try:
		# split line into list
		row = line.split(b'\t')
		# get index of data format
		data_fmts = row[8].split(b':')
		# may be multiallelic; only first alt allele will be used unless split_multi_GT; later data with bad values will be filtered out
		alt_alleles = row[4].split(b',')
		row[4] = alt_alleles[0]
		# filter sample columns to include only data in desired format; sampleIDs start at file column index 9; in new row sampleIDs now start at 4, ALT is column 3
		if (len(data_fmts) > 1):
			data_ind = data_fmts.index(bformat)
			row[8] = bformat
			row = [row[x] if x <= 8 else row[x].split(b':')[data_ind] for x in col_inds]
		else:
			row = [row[x] for x in col_inds]
		# turn multi-allelic lines into multiple biallelic lines
		if split_multi_GT & (len(alt_alleles) > 1):
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
			line += b'' if line.endswith(b'\n') else b'\n'
		return line
	except:
		return b''

def filter_weight_line(line: bytes, btarget: bytes, target_ind, col_inds):
	# split line into list
	row = line.split(b'\t')
	# check if row is for correct target
	if (row[target_ind].startswith(btarget)):
		line = b'\t'.join([row[x] for x in col_inds])
		line += b'' if line.endswith(b'\n') else b'\n'
		return line
	else: 
		return b''


def filter_other_line(line: bytes, col_inds):
	# split line into list
	row = line.split(b'\t')
	# filter out unneeded columns
	line = b'\t'.join([row[x] for x in col_inds])
	line += b'' if line.endswith(b'\n') else b'\n'
	return line


def read_tabix(start, end, sampleID, chrm, path, file_cols, col_inds, cols, dtype, genofile_type=None, data_format=None, target_ind=5, target=None, weight_threshold=0, raise_error=True, **kwargs):

	# subprocess command
	command_str = ' '.join(['tabix', path, chrm + ':' + start + '-' + end])

	proc = subprocess.Popen(
		[command_str],
		shell=True,
		stdout=subprocess.PIPE,
		bufsize=1)

	# initialize bytearray
	proc_out = bytearray()

	# set correct filter function by file type
	if genofile_type == 'vcf':
		bformat = str.encode(data_format)
		filter_line = functools.partial(filter_vcf_line, bformat=bformat, col_inds=col_inds, split_multi_GT=data_format == 'GT')
	elif genofile_type == 'weight':
		btarget = str.encode(target)
		filter_line = functools.partial(filter_weight_line, btarget=btarget, target_ind=target_ind, col_inds=col_inds)
	elif genofile_type == 'bgw_weight':
		btarget = b'0'
		filter_line = functools.partial(filter_weight_line, btarget=btarget, target_ind=target_ind, col_inds=col_inds)
	else:
		filter_line = functools.partial(filter_other_line, col_inds=col_inds)

	# while subprocesses running, read lines into byte array
	while proc.poll() is None:
		line = proc.stdout.readline()
		if len(line) == 0:
			break
		proc_out += filter_line(line)
	# read in lines still remaining after subprocess completes
	for line in proc.stdout:
		proc_out += filter_line(line)

	if not proc_out and raise_error:
		print('No tabix data for target.\n')
		raise NoTargetDataError

	# read data into dataframe
	df = pd.read_csv(
		StringIO(proc_out.decode('utf-8')),
		sep='\t',
		low_memory=False,
		header=None,
		names=cols,
		dtype=dtype)
	
	# kill process
	proc.kill()

	# filter out rows where all sampleID values are nan
	if len(sampleID):
		df = df[df[sampleID].count(axis=1) != 0].reset_index(drop=True)

	df = optimize_cols(df)

	# get snpID
	if (genofile_type != 'weight') or (not 'snpID' in cols):
		df['snpID'] = get_snpIDs(df)

	# filter out duplicated snpIDs
	df = df.drop_duplicates(['snpID'], keep='first').reset_index(drop=True)

	# weight file handling
	if (genofile_type == 'weight'):
		## figure out a way to handle target ids with decimal places when both files dont necessarily have those
		# if not np.all(df['TargetID'] == target):
		# partial_targetids = np.unique(df['TargetID'])

		# VC_TWAS uses ES = b + beta
		if (not 'ES' in cols):
			if (('b' in cols) and ('beta' in cols)):
				df['ES'] = df['b'] + df['beta']
			if (('PCP' in cols) and ('beta' in cols)):
				df['ES'] = np.prod(df['PCP'], df['beta'])

		if weight_threshold:
			# filter out weights below threshold
			df = df[operator.gt(np.abs(df['ES']), weight_threshold)].reset_index(drop=True)

			if df.empty and raise_error:
				print('No test SNPs with cis-eQTL weights with magnitude that exceeds specified weight threshold for TargetID: ' + target + '.\n')
				raise NoTargetDataError

	# remove rows with non-valid GT values; ie those from multiallelic rows
	if (data_format == 'GT'):
		valid_GT = ['.|.', '0|0', '0|1', '1|0', '1|1', 
		'./.', '0/0', '0/1', '1/0', '1/1']
		df = df[np.all(df[sampleID].isin(valid_GT), axis=1)].reset_index(drop=True)

	# process dosage, vcf file sample df
	if (data_format == 'GT') or (data_format == 'DS'):
		df = reformat_sample_vals(df, data_format, sampleID)

	if df.empty and raise_error:
		print('No valid tabix data for target.\n')
		raise NoTargetDataError

	return df


def tabix_query_file(path, reg_str):
	proc = subprocess.Popen(['tabix '+path+reg_str+' | head -n1 '],
		shell=True, stdout=subprocess.PIPE, bufsize=1)
	try:
		outs, err = proc.communicate(timeout=15)
		ret_val = len(outs) > 1
	except TimeoutExpired:
		proc.kill()
		outs, errs = proc.communicate()
		ret_val = len(outs) > 1
	except:
		proc.kill()
		ret_val = False
	finally:
		# kill process
		proc.kill()

	return ret_val


def tabix_query_files(start, end, chrm, geno_path=None, gwas_path=None, w_path=None, z_path=None, **kwargs):
	paths = [path for path in [geno_path, gwas_path, w_path, z_path] if path is not None]
	reg_str = ' ' + chrm + ':' + start + '-' + end

	# for each path test if length of first line != 0
	return np.all([tabix_query_file(path, reg_str)  for path in paths])


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
	
	# kill process
	proc.kill()

	return proc_out


# Call tabix to get header, read in lines into byte array;
# method for reading in headers for vcf files
def call_tabix_header(path, out='tuple', rename={}):
	# create dictionary for renaming columns
	# the first column in the vcf file is expected to be parsed as #CHROM;  automatically add that to the rename list
	rename = {**{'#CHROM':'CHROM'}, **rename}

	proc = subprocess.Popen(
		['tabix -H ' + path],
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
		if not line.startswith(b'##'):
			proc_out += line
	# read in lines still remaining after subprocess completes
	for line in proc.stdout:
		if not line.startswith(b'##'):
			proc_out += line

	# decode bytes, use pandas to read in, rename columns from dictionary
	header = pd.read_csv(
		StringIO(proc_out.decode('utf-8')),
		sep='\t',
		error_bad_lines=False).rename(columns=rename)
	
	# kill process
	proc.kill()

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
	
	# kill process
	proc.kill()   

	if out=='tuple':
		return tuple(header)

	elif out=='list':
		return list(header)

	return header

# determine indices of file cols to read in, dtype of each col
def get_cols_dtype(file_cols, cols, sampleid=None, genofile_type=None, add_cols=[], drop_cols=[], get_id=False, ind_namekey=False, **kwargs):

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
		'PCP': np.float64,
		'POS': np.int64,
		'QUAL': object,
		'REF': object,
		'SE': np.float64,
		'snpID': object,
		'TargetID': object,
		'Trans': np.int64,
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
		out_dtype_dict = {x:dtype_dict[x] for x in cols}
	else:
		# return dtype dict with keys as column index
		ind_dtype_dict = {file_cols.index(x):dtype_dict[x] for x in cols}
		out_dtype_dict = ind_dtype_dict

	# if ret_dict:
	return {
		'file_cols': file_cols,
		'cols': [file_cols[i] for i in col_inds], 
		'col_inds': col_inds, 
		'dtype': out_dtype_dict}
	
	# return col_inds, out_dtype_dict

# def exp_cols_dtype(file_cols, sampleid, ind_namekey=True, ret_dict=True, **kwargs):
# 	return get_cols_dtype(file_cols, 
# 		cols=['CHROM', 'GeneStart', 'GeneEnd', 'TargetID', 'GeneName'], 
# 		sampleid=sampleid, ind_namekey=ind_namekey, ret_dict=ret_dict)

# def genofile_cols_dtype(file_cols, genofile_type, sampleid, ind_namekey=True, ret_dict=True, **kwargs):
# 	return get_cols_dtype(file_cols, 
# 		cols=['CHROM','POS','REF','ALT'],
# 		sampleid=sampleid, genofile_type=genofile_type, 
# 		ind_namekey=ind_namekey, ret_dict=ret_dict)

# def weight_cols_dtype(file_cols, add_cols=[], drop_cols=[], get_id=True, ret_dict=True, ind_namekey=True, **kwargs):
# 	return get_cols_dtype(file_cols, 
# 		cols=['CHROM','POS','REF','ALT','TargetID','ES'], 
# 		add_cols=add_cols, drop_cols=drop_cols, 
# 		get_id=get_id, ret_dict=ret_dict, ind_namekey=ind_namekey)

# def zscore_cols_dtype(file_cols, ret_dict=True, ind_namekey=True, **kwargs):
# 	return get_cols_dtype(file_cols, 
# 		cols=['CHROM','POS','REF','ALT','Zscore'], ret_dict=ret_dict, ind_namekey=ind_namekey)

def gwas_cols_dtype(file_cols, ind_namekey=True, **kwargs):
	return get_cols_dtype(file_cols, 
		cols=['CHROM','POS','REF','ALT','BETA','SE'], ind_namekey=ind_namekey)

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

# yields formatted tabix regions strings
def get_regions_list(snp_ids):

	# split into groups by chromosome
	for chrm, grp in groupby(snp_ids, lambda x: x.split(':')[0]):
		# snp pos values as integers
		pos_vals = [int(snp.split(':')[1]) for snp in grp]

		# get intervals of start,end positions; convert to tabix string
		for x, y in groupby(enumerate(pos_vals), lambda p: p[1]-p[0]):
			y = list(y)

			# chrm:start-end
			yield chrm + ':' + str(y[0][1]) + '-' + str(y[-1][1])


# call tabix using regions string
def call_tabix_regions(path, regs_str, filter_line = lambda x:x ):

	proc = subprocess.Popen(
		['tabix '+path+' '+regs_str],
		shell=True,
		stdout=subprocess.PIPE,
		bufsize=1)
	proc_out = bytearray()

	# process while subprocesses running
	while proc.poll() is None:
		line = proc.stdout.readline()
		if len(line) == 0:
			break
		proc_out += filter_line(line)

	# leftover lines
	for line in proc.stdout:
		proc_out += filter_line(line)

	# kill process
	proc.kill()

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


def get_ld_matrix(MCOV, return_diag=False):
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
				V_upper[i,j] = cov_i[inds[j] - inds[i]]
			else:
				V_upper[i,j] = 0

	snp_Var = V_upper.diagonal()
	V = V_upper + V_upper.T - np.diag(snp_Var)
	snp_sd = np.sqrt(snp_Var)

	if return_diag:
		return snp_sd, V, snp_Var

	else:
		return snp_sd, V


# return snp ids; join CHROM, POS, REF, ALT columns into : separated string
def get_snpIDs(df: pd.DataFrame, flip=False):
	chrms = df['CHROM'].astype('str').values
	pos = df['POS'].astype('str').values
	ref = df['REF'].values
	alt = df['ALT'].values
	if flip:
		return [':'.join(i) for i in zip(chrms,pos,alt,ref)]
	else:
		return [':'.join(i) for i in zip(chrms,pos,ref,alt)]


# gives flipped snpIDs for an array/list of snpIDs
def flip_snpIDs(snpIDs):
	return np.array([':'.join([y[0],y[1],y[3],y[2]]) for y in [x.split(':') for x in snpIDs]])


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


# returns a boolean array; whether substring is in a np object array
def substr_in_strarray(substr, strarray):
   return np.frompyfunc(lambda x: substr in x, 1,1)(strarray)


# drop variants with missing rate that exceeds threshold
def handle_missing(df: pd.DataFrame, sampleID, missing_rate, filter=True, op=operator.le):
	df = df.copy()

	# exclude rows where all sample data is NaN
	df = df[df[sampleID].count(axis=1) > 0].reset_index(drop=True)
	vals = df[sampleID].values

	# calculate missing rate for each snp, read into series and downcast
	MissingRate = pd.to_numeric(
		pd.Series(
			np.apply_along_axis(lambda x: np.count_nonzero(np.isnan(x))/x.size, 1, vals), 
			name='missing_rate'),
		downcast='float')

	# re combine dataframe
	df = pd.concat([df.drop(columns=sampleID), pd.DataFrame(vals, columns=sampleID), MissingRate], axis=1)

	if filter:
		df = df[op(df.missing_rate, missing_rate)].reset_index(drop=True)

	return df

# function to perform MAF calculation/imputation on each row
def row_maf_impute(x):
	# calculate MAF
	MAF_val = np.array(
		np.nansum(x) / (2 * (x.size - np.count_nonzero(np.isnan(x)))), 
		dtype=x.dtype)
	# impute missing values
	x[np.where(np.isnan(x))] = 2 * MAF_val
	return np.append(x, MAF_val)

def calc_maf(df: pd.DataFrame, sampleID, maf, filter=True, op=operator.gt, filter_bid=False):
	df = df.reset_index(drop=True).copy()
	vals = df[sampleID].values

	# calculate MAF, impute missing
	sample_MAF = np.apply_along_axis(row_maf_impute, 1, vals)

	# re combine dataframe (faster than editing in place)
	df = pd.concat([df.drop(columns=sampleID), pd.DataFrame(sample_MAF, columns=[*sampleID, 'MAF'])], axis=1)	

	if filter:
		df = df[op(df.MAF, maf)].reset_index(drop=True)
		if filter_bid:
			df = df[op(1 - df.MAF, maf)].reset_index(drop=True)

	return df

# Calculate, filter p val of HWE
def calc_p_hwe(df: pd.DataFrame, sampleID, pval, filter=True, op=operator.gt):
	df = df.reset_index(drop=True).copy()
	vals = df[sampleID].values

	HWE = pd.to_numeric(
		pd.Series(
			np.apply_along_axis(lambda x:prep_p_HWE(x[np.where(~np.isnan(x))]), 1, vals), 
			name='p_HWE'),
		downcast='float')

	# re combine dataframe
	df = pd.concat([df.drop(columns=sampleID), pd.DataFrame(vals, columns=sampleID), HWE], axis=1)	

	if filter:
		df = df[op(df.p_HWE, pval)].reset_index(drop=True)

	return df

### Prepare for HWE input
def prep_p_HWE(vals):
	if not vals.size:
		p_hwe = np.nan
	else:
		N_hets = vals[(vals >= 0.5) & (vals < 1.5)].size
		N_aa = vals[(vals >= 0) & (vals < 0.5)].size
		N_AA = vals[(vals >= 1.5) & (vals <= 2)].size

		p_hwe = HWE(N_hets, N_AA, N_aa)

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
def HWE(obs_hets, obs_hom1, obs_hom2):
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

	# #alternate p-value calculation for p_hi/p_lo
	# p_hi = float(het_probs[obs_hets])
	# for i in range(obs_hets,rare_copies+1):
	# 	p_hi += het_probs[i]

	# p_lo = float(het_probs[obs_hets])
	# for i in range(obs_hets-1,-1,-1):
	# 	p_lo += het_probs[i]

	# p_hi_lo = 2.0 * p_hi if p_hi < p_lo else 2.0 * p_lo

	p_hwe = 0.0
	#  p-value calculation for p_hwe
	for i in range(0,rare_copies + 1):
		if het_probs[i] > het_probs[obs_hets]:
			continue
		p_hwe += het_probs[i]

	p_hwe = 1.0 if p_hwe > 1.0 else p_hwe

	return p_hwe


# center
def center(df: pd.DataFrame, sampleID=None):
	df = df.reset_index(drop=True).copy()
	if sampleID is None:
		vals = df.values
		vals = np.apply_along_axis(lambda x: x - np.mean(x), 1, vals)
		df = pd.DataFrame(vals, columns=df.columns)
	else: 
		vals = df[sampleID].values
		vals = np.apply_along_axis(lambda x: x - np.mean(x), 1, vals)
		df = pd.concat([df.drop(columns=sampleID), 
			pd.DataFrame(vals, columns=sampleID)], axis=1)	
	return df


# print args (for testing)
def print_args(args):
	for key, value in args.__dict__.items():
		if isinstance(value, str):
			print('args.', key, ' = \'', value, '\'', sep='')
		else:
			print('args.', key, ' = ', value, sep='')
	print('\n')
