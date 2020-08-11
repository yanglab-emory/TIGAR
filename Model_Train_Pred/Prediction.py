#!/usr/bin/env python

###############################################################
# import packages needed
import argparse
import operator
import multiprocessing
import subprocess
import sys

from io import StringIO
from time import time

import numpy as np
import pandas as pd

###############################################################
# time calculation
start_time=time()

###############################################################
# parse input arguments
parser = argparse.ArgumentParser(description='Prediction')

# Specify tool directory
parser.add_argument('--TIGAR_dir',type=str)

# eQTL weight file path
parser.add_argument('--weight',type=str,dest='w_path')

# Test sampleID path
parser.add_argument('--test_sampleID',type=str,dest='sampleid_path')

# Specified chromosome number
parser.add_argument('--chr',type=str)

# Test genotype file path
parser.add_argument('--genofile',type=str,dest='geno_path')

# Specified input file type (vcf or dosages)
parser.add_argument('--genofile_type',type=str)

# 'DS' or 'GT' for VCF genotype file
parser.add_argument('--format',type=str)

# window
parser.add_argument('--window',type=int)

# Gene annotation file path
parser.add_argument('--gene_anno',type=str,dest='annot_path')

# number of thread
parser.add_argument('--thread',type=int)

# Threshold of difference of maf between training data and testing data
parser.add_argument('--maf_diff',type=float)

# output dir
parser.add_argument('--out_dir',type=str)

args = parser.parse_args()

sys.path.append(args.TIGAR_dir)
###############################################################
import TIGARutils as tg

# return correct snpID and GT value
# 2-GT value if matching snpID is flipped wrt Weight snpID
def handle_flip_pred(df: pd.DataFrame, sampleID, orig_overlap, flip_overlap):
    df = df.copy()
    orig = df['IDorig'].values
    flip = df['IDflip'].values
    origmaf = df['MAF'].values
    df = df[sampleID]

    ids = np.empty_like(orig)
    maf = np.empty_like(origmaf)
    outdf = pd.DataFrame(columns = sampleID)

    for i in range(len(df)):
        if orig[i] in orig_overlap:
            ids[i], maf[i] = orig[i], origmaf[i].astype('float')
            outdf = outdf.append(df.iloc[i]).astype('float')
        elif flip[i] in flip_overlap:
            ids[i], maf[i] = flip[i], 1-origmaf[i].astype('float')
            outdf = outdf.append(df.iloc[i].apply(lambda x: 2-x)).astype('float')

    return ids, maf, outdf

#######################################################################
# Input Arguments for GReX Prediction

# --chr: Chromosome number need to be specified with respect to the genotype input data
# --weight: Path for SNP weight (eQTL effect size) file 
# --test_sampleID: Path for a file with sampleIDs that should be contained in the genotype file
# --genofile: Path for the testing genotype file (bgzipped and tabixed) 
# --genofile_type: Genotype file type: 'vcf' or 'dosage'
# --genofile_colnames: File with column heads of genotype file
# --format: Genotype format in VCF file that should be used: 'GT' (default) for genotype data or 'DS' for dosage data, only required if the input genotype file is of VCF file
# --window: Window size around gene transcription starting sites (TSS) for selecting cis-SNPs for fitting gene expression prediction model (default 1000000 for +- 1MB region around TSS)
# --maf_diff: MAF difference threshold for matching SNPs from eQTL weight file and test genotype file. If SNP MAF difference is greater than maf_diff (default 0.2), , the SNP will be excluded
# --thread: Number of threads for parallel computation (default 1)
# --out_dir: Output directory (will be created if not exist)

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

out_pred_path = args.out_dir + '/CHR' + args.chr + '_Pred_GReX.txt'

###############################################################
# Print input arguments
print(
'''********************************
Input Arguments

Gene annotation file specifying genes for prediction: {annot_path}

Prediction sampleID file: {sampleid_path}

Chromosome: {chr}

cis-eQTL weight file: {w_path}

Prediction genotype file: {geno_path}

Genotype file used for prediction is type: {genofile_type}

Genotype data format: {format}

Gene prediction region SNP inclusion window: +-{window}

Excluding SNPs matched between eQTL weight file and prediction genotype file if MAF difference exceeds: |{maf_diff}|

Number of threads: {thread}

Output directory: {out_dir}

Output prediction results file: {out_path}
********************************'''.format(
    **args.__dict__,
    out_path = out_pred_path))

###############################################################
# Load eQTL weights (ES)
w_cols = tg.get_header(args.w_path)
w_cols_ind, w_dtype = tg.weight_cols_dtype(w_cols, ['MAF'])

print('Reading eQTL weights file.')
Weight_chunks=pd.read_csv(
    args.w_path, 
    sep='\t', 
    iterator=True, 
    chunksize=10000,
    usecols=w_cols_ind,
    dtype=w_dtype)

Weight = pd.concat([x[x['CHROM']==args.chr] for x in Weight_chunks]).reset_index(drop=True)

if Weight.empty:
    raise SystemExit('There are no valid eQTL weights.')

# Weight.columns = [w_cols[i] for i in tuple(Weight.columns)]
Weight = tg.optimize_cols(Weight)

if 'ID' in Weight.columns:
    Weight.rename(columns={'ID':'snpID'})

if not 'snpID' in Weight.columns:
    Weight['snpID'] = tg.get_snpIDs(Weight)

Weight = Weight.drop(columns=['CHROM','POS','REF','ALT'])

# Load annotation file
print('Reading gene annotation file.')
Gene_chunks = pd.read_csv(
    args.annot_path, 
    sep='\t', 
    iterator=True, 
    chunksize=10000,
    dtype={'CHROM':object,'GeneStart':np.int64,'GeneEnd':np.int64,'TargetID':object,'GeneName':object}, 
    usecols=['CHROM','GeneStart','GeneEnd','TargetID','GeneName'])
Gene = pd.concat([x[x['CHROM'] == args.chr] for x in Gene_chunks]).reset_index(drop=True)

if Gene.empty:
    raise SystemExit('There are no valid annotations.')

Gene = tg.optimize_cols(Gene)

# Load genotype column names of test genotype file
print('Reading genotype file header.\n')
g_cols = tg.call_tabix_header(args.geno_path)
gcol_sampleids = g_cols[gcol_sampleids_strt_ind:]

# get intersection
# Load test sample IDs
print('Reading sampleID file.\n')
spec_sampleids = pd.read_csv(
    args.sampleid_path,
    sep='\t',
    header=None)[0].drop_duplicates()

print('Matching sampleIDs.\n')
sampleID = np.intersect1d(gcol_sampleids, spec_sampleids)

if not sampleID.size:
    raise SystemExit('There are no sampleID in both the genotype data and the specified sampleID file.')

# get the column indices and dtypes for reading genofile into pandas
g_cols_ind, g_dtype = tg.genofile_cols_dtype(g_cols, args.genofile_type, sampleID)

# Define TargetID
TargetID = Gene.TargetID.values
n_targets = TargetID.size

print('Creating output file: ' + out_pred_path + '\n')
out_cols = np.concatenate((
    ['CHROM','GeneStart','GeneEnd','TargetID','GeneName'], sampleID))
pd.DataFrame(columns=out_cols).to_csv(
    out_pred_path,
    sep='\t', 
    index=None, 
    header=True, 
    mode='w')

print('********************************\n')

###############################################################
# thread function
def thread_process(num):
    try:
        target = TargetID[num]
        print('num=' + str(num) + '\nTargetID=' + target)
        Gene_info = Gene.iloc[[num]].reset_index(drop=True)

        start = str(max(int(Gene_info.GeneStart)-args.window, 0))
        end = str(int(Gene_info.GeneEnd)+args.window)

        print('Reading genotype data.')
        g_proc_out = tg.call_tabix(args.geno_path, args.chr, start, end)
        
        if not g_proc_out:
            print('There is no genotype data for TargetID: ' + target + '\n')
            return None
        
        target_geno = pd.read_csv(StringIO(g_proc_out.decode('utf-8')),
                sep='\t',
                low_memory=False,
                header=None,
                usecols=g_cols_ind,
                dtype=g_dtype)
        target_geno.columns = [g_cols[i] for i in target_geno.columns]
        target_geno = tg.optimize_cols(target_geno)

        # Intersect SNPs from eQTL weight file and test genotype file
        # initial filter to reduce amount of dataframe processing
        print('Getting weight data for target.')
        target_weight = Weight[Weight.TargetID==target][['snpID', 'ES', 'MAF']]

        if target_weight.empty:
            print('No cis-eQTL weights for TargetID: ' + target + '\n')
            return None       

        # vcf files may have data in multiple formats, check if this is the case and remove unnecessary formats. requires that all rows have data in the user-specified format
        if args.genofile_type=='vcf':
            target_geno = tg.check_prep_vcf(target_geno, args.format, sampleID)

        # reformat values in target_geno data frame
        target_geno[sampleID]=target_geno[sampleID].apply(lambda x:tg.reformat_sample_vals(x,args.format), axis=0)

        # calculate MAF
        target_geno = tg.calc_maf(target_geno, sampleID, 0, op=operator.ge)

        # Get original and flipped snpIDs, filter out duplicates
        target_geno['IDorig'] = tg.get_snpIDs(target_geno)
        target_geno['IDflip'] = tg.get_snpIDs(target_geno, flip=True)
        target_geno = target_geno.drop(columns=['CHROM','POS','REF','ALT'])
        target_geno = target_geno.drop_duplicates(['IDorig'], keep='first')

        # get overlapping snps
        snp_overlap_orig = np.intersect1d(target_weight.snpID, target_geno.IDorig)
        snp_overlap_flip = np.intersect1d(target_weight.snpID, target_geno.IDflip)
        snp_overlap = np.concatenate((snp_overlap_orig, snp_overlap_flip))

        if not snp_overlap.size:
            print('No overlapping test SNPs between weight and genotype file for TargetID: ' + target + '\n')
            return None

        # print('Number of SNPs overlapped between eQTL weight file and test genotype file: ' + str(snp_overlap.size))

        target_weight = target_weight[target_weight.snpID.isin(snp_overlap)]
        target_geno = target_geno[target_geno.IDorig.isin(snp_overlap_orig) | target_geno.IDflip.isin(snp_overlap_flip)]

        # HANDLE FLIPPED SNPS 
        if (snp_overlap_orig.size > 0) and (snp_overlap_flip.size > 0):
            # assume mix of flipped, non-flipped
            target_geno['snpID'], target_geno['MAF_test'], target_geno[sampleID] = handle_flip_pred(target_geno, sampleID,snp_overlap_orig, snp_overlap_flip)

        elif snp_overlap_orig.size == snp_overlap.size:
            # assume all non-flipped 
            target_geno[['snpID','MAF_test']] = target_geno[['IDorig','MAF']]

        else:
            # assume all flipped
            target_geno['snpID'], target_geno['MAF_test'] = target_geno['IDflip'], 1-target_geno['MAF'].astype('float')
            target_geno[sampleID] = target_geno[sampleID].apply(lambda x: 2-x).astype('float')

        target_geno = target_geno.drop(columns=['IDorig','IDflip','MAF'])

        # merge target_geno, target_weight
        Pred = target_geno.merge(
            target_weight, 
            left_on='snpID', 
            right_on='snpID', 
            how='outer')

        Pred['diff'] = np.abs(Pred['MAF'].astype('float') - Pred['MAF_test'].astype('float'))
        
        Pred = Pred[Pred['diff'] <= args.maf_diff].drop(columns=['MAF','MAF_test','diff'])

        if Pred.empty:
            print('All SNP MAFs for training data and testing data differ by a magnitude greater than '+str(args.maf_diff) + ' for TargetID: ' + target + '\n')
            return None

        # print('Number of SNPs used for prediction after filtering by maf_diff: ' + str(Pred.snpID.size))

        print('Predicting GReX.\nN SNPs=' + str(Pred.snpID.size))

        # output results
        pred_grex = pd.DataFrame(
            data=np.dot(Pred[sampleID].T, Pred['ES'].values),
            index=sampleID).T

        result = pd.concat([Gene_info, pred_grex], axis=1)

        result.to_csv(
            out_pred_path,
            sep='\t',
            index=None,
            header=None,
            mode='a')

        print('Target prediction completed.\n')

    except Exception as e:
        e_type, e_obj, e_tracebk = sys.exc_info()
        e_line_num = e_tracebk.tb_lineno

        print('Caught a type {} exception for TargetID={}, num={} on line {}:\n{}'.format(e_type, target, num, e_line_num, e))

    finally:
        # print info to log do not wait for buffer to fill up
        sys.stdout.flush()       

###############################################################
# thread process
if __name__ == '__main__':
    print('Starting prediction for ' + str(n_targets) + ' target genes.\n')
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


