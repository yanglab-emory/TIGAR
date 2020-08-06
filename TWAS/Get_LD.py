#!/usr/bin/env python

###############################################################
# Import packages needed
import argparse
import io
from io import StringIO
import multiprocessing
import operator
import subprocess
import sys
from time import time

import pandas as pd
import numpy as np

###############################################################
### time calculation
start_time = time()

###############################################################
### variable needed
parser = argparse.ArgumentParser(description='Get LD')

### Specify tool directory
parser.add_argument('--TIGAR_dir',type=str)

### chromosome block information
parser.add_argument('--genome_block',type=str,dest='block_path')

### genotype file dir
parser.add_argument('--genofile',type=str,dest='geno_path')

### specified input file type(vcf or doasges)
parser.add_argument('--genofile_type',type=str)

### chromosome number
parser.add_argument('--chr',type=str)

### 'DS' or 'GT'
parser.add_argument('--format',type=str)

### maf threshold for seleting genotype data to calculate covariance matrix
parser.add_argument('--maf',type=float)

### number of thread
parser.add_argument('--thread',type=int)

### output dir
parser.add_argument('--out_dir',type=str)

### sampleID path
parser.add_argument('--sampleID',type=str,dest='sampleid_path')

args = parser.parse_args()

sys.path.append(args.TIGAR_dir)

###############################################################
# DEFINE, IMPORT FUNCTIONS
import TIGARutils as tg


# limit to 4 decimal places max, strip trailing 0s
def cov_print_frmt(x): return ("%.4f" % x).rstrip('0').rstrip('.')

###############################################################
### variable checking
print("\n \nGenome block annotation based on LD structure : " + args.block_path + "\n")
print("Reference genotype file: " + args.geno_path + "\n")

if args.genofile_type=='vcf':
    print("Using genotype data of format " + args.format + " in the VCF file.\n")
    # cols in file:
    # CHROM     POS     ID  REF     ALT     QUAL    FILTER  INFO    FORMAT  Sample1 ... Samplen
    # cols read in:
    # CHROM     POS     REF     ALT     FORMAT  Sample1 ... Samplen
    gcol_sampleids_strt_ind = 9

elif args.genofile_type=='dosage':
    print("Using genotype data from the dosage type file."+ "\n")
    args.format = 'DS'
    # cols in file:
    # CHROM     POS     ID  REF     ALT Sample1 ... Samplen
    # cols read in (excludes ID):
    # CHROM     POS     REF     ALT Sample1 ... Samplen
    gcol_sampleids_strt_ind = 5

else:
    raise SystemExit("Genotype file type (--genofile_type) should be specified as either vcf or dosage."+ "\n")

print("Chromosome number: " + args.chr + "\n")
print("Number of threads: " + str(args.thread)+ "\n")
print("Only using variants with MAF > " + str(args.maf) + " to calculate reference LD covariance file."+ "\n")
print("Output directory: " + args.out_dir + "\n")

out_ref_cov_path = args.out_dir+'/CHR'+args.chr+'_reference_cov.txt'
print("Reference covariance results file: " + out_ref_cov_path+"\n")

if args.sampleid_path:
    print("SampleID path: " + args.sampleid_path + "\n")

###############################################################
### Read in block information
chr_blocks = pd.read_csv(
    args.block_path,
    sep='\t',
    usecols=['CHROM','Start','End'],
    dtype={'CHROM':object,'Start':object,'End':object})
chr_blocks = chr_blocks[chr_blocks['CHROM']==args.chr].reset_index(drop=True)
chr_blocks = optimize_cols(chr_blocks)

n_blocks = len(chr_blocks)

g_cols = tg.call_tabix_header(args.geno_path)
gcol_sampleids = g_cols[gcol_sampleids_strt_ind:]

# if user specified path with sampleids, and at least one sampid from that file, get intersection
if args.sampleid_path:
    spec_sampleids = pd.read_csv(
        args.sampleid_path,
        sep='\t',
        header=None)
    spec_sampleids = spec_sampleids[0].drop_duplicates()

    sampleID = np.intersect1d(gcol_sampleids, spec_sampleids)

else:
    sampleID = gcol_sampleids

if not sampleID.size:
    raise SystemExit('There are no sampleID in both the genotype data and the specified sampleID file.')

# get the indices and dtypes for reading files into pandas
g_cols_ind, g_dtype = tg.genofile_cols_dtype(g_cols, args.genofile_type, sampleID)

# write columns out to file
out_cols = ['#CHROM', 'POS', 'REF', 'ALT', 'snpID','COV']
pd.DataFrame(columns=out_cols).to_csv(
    out_ref_cov_path,
    sep='\t',
    index=None,
    header=True,
    mode='w')

###############################################################
def thread_process(num):
    try:
        block = chr_blocks.loc[num]
        
        g_proc_out = tg.call_tabix(args.geno_path, args.chr, block.Start, block.End)

        if not g_proc_out:
            print("There is no genotype data in this block.")
            return None

        block_geno = pd.read_csv(StringIO(g_proc_out.decode('utf-8')),
            sep='\t',
            low_memory=False,
            header=None,
            usecols=g_cols_ind,
            dtype=g_dtype)

        block_geno.columns = [g_cols[i] for i in block_geno.columns]
        block_geno = optimize_cols(block_geno)

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
            snpinfo = block_geno.iloc[i][['CHROM','POS','REF','ALT','snpID']]
            snpcov = ','.join([cov_print_frmt(x) for x in mcovar[i,i:]])
            covar_info = pd.DataFrame(np.append(snpinfo, snpcov)).T
            covar_info.to_csv(
                out_ref_cov_path,
                sep='\t',
                index=None,
                header=None,
                mode='a')

    except Exception as e:
        e_type, e_obj, e_tracebk = sys.exc_info()
        e_line_num = e_tracebk.tb_lineno

        e, e_type, e_line_num = [str(x) for x in [e, e_type, e_line_num]]

        print('Caught a type '+ e_type +' exception for block num='+str(num)+' on line '+e_line_num+':\n' + e )

    finally:
        # print info to log do not wait for buffer to fill up
        sys.stdout.flush()

##################################################################
### thread process

if __name__ == '__main__':
    print("Starting LD calculation for "+str(n_blocks)+" blocks.")
    pool = multiprocessing.Pool(args.thread)
    pool.map(thread_process,[num for num in range(n_blocks)])
    pool.close()
    pool.join()

################################################################
### time calculation
elapsed_sec = time()-start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print("Computation time (DD:HH:MM:SS): " + elapsed_time)

