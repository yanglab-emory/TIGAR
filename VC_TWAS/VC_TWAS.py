#!/usr/bin/env python

############################################################
# import packages needed
import argparse
import operator
import multiprocessing
import subprocess
import sys
import traceback

from io import StringIO
from time import time

import pandas as pd
import numpy as np

##########################################################
# time calculation
start_time = time()

#################################################
# Input Arguments for VC_TWAS 
# --chr: Chromosome number need to be specified with respect to the genotype input data
# --weight: Path for SNP weight (eQTL effect size) file 
# --test_sampleID: Path for a file with sampleIDs that should be contained in the genotype file
# --genofile: Path for the training genotype file (bgzipped and tabixed) 
# --genofile_type: Genotype file type: 'vcf' or 'dosage'
# --test_geno_colnames: File with column heads of genotype file
# --format: Genotype format in VCF file that should be used: 'GT' (default) for genotype data or 'DS' for dosage data, only required if the input genotype file is of VCF file
# --window: Window size around gene transcription starting sites (TSS) for selecting cis-SNPs for fitting gene expression prediction model (default 1000000 for +- 1MB region around TSS)
# --phenotype_type: phenotype type, continous phenotype is 'C' and binomial phenotype is 'D'
# -- maf: threshold for seleting genotype data
# --weight_threshold: threshold for weights estimated from DPR and Elastic Net(absolute value), default threshold is 10^-4
# --thread: Number of threads for parallel computation (default 1)
# --PED: PED file, phenotype 
# --PED_info: association information file for PED
# --out_dir: Output directory (will be created if not exist)

#######################################
# parse input arguments
parser = argparse.ArgumentParser(description='VC TWAS')

# Specify tool directory
parser.add_argument('--TIGAR_dir',type=str)

# Specified chromosome number
parser.add_argument('--chr',type=str)

# eQTL weight file path
parser.add_argument('--weight',type=str,dest='w_path')

# Test sampleID path
parser.add_argument('--test_sampleID',type=str,dest='sampleid_path')

# Gene annotation file path
parser.add_argument('--gene_anno',type=str,dest='annot_path')

# Test genotype file path
parser.add_argument('--genofile',type=str,dest='geno_path')

# Specified input file type(vcf or dosages)
parser.add_argument('--genofile_type',type=str)

# 'DS' or 'GT' for VCF genotype file
parser.add_argument('--format',type=str)

# window
parser.add_argument('--window',type=int)

# phenotype_type
parser.add_argument('--phenotype_type',type=str)

# maf threshold for seleting genotype data
parser.add_argument('--maf',type=float)

# hwe threshold for seleting genotype data
parser.add_argument('--hwe',type=float)

# weight_threshold
parser.add_argument('--weight_threshold',type=float)

# number of thread
parser.add_argument('--thread',type=int)

# PED file path
parser.add_argument('--PED',type=str,dest='ped_path')

# Association Information file
parser.add_argument('--PED_info',type=str,dest='pedinfo_path')

# output dir
parser.add_argument('--out_dir',type=str)

args = parser.parse_args()

sys.path.append(args.TIGAR_dir)
sys.path.append(args.TIGAR_dir + '/VC_TWAS')

#############################################
# Import TIGAR functions
import TIGARutils as tg
import SKAT

#############################################
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

if (args.phenotype_type != 'C') and (args.phenotype_type != 'D'):
    raise SystemExit('Please specify phenotype type (--phenotype_type) as either "C" for continous or "D" for dichotomous.\n')

out_twas_path = args.out_dir + '/CHR' + args.chr + '_indv_VC_TWAS.txt'

###############################################################
# Print input arguments
print(
'''********************************
Input Arguments

Gene annotation file specifying genes for TWAS: {annot_path}

PED phenotype/covariate data file: {ped_path}

PED information file: {pedinfo_path}

Test sampleID file: {sampleid_path}

Chromosome: {chr}

cis-eQTL weight file: {w_path}

Test genotype file: {geno_path}

Genotype file used for testing is type: {genofile_type}

Genotype data format: {format}

Gene testing region SNP inclusion window: +-{window}

MAF threshold for SNP inclusion: {maf}

HWE p-value threshold for SNP inclusion: {hwe}

SNP weight inclusion threshold: {weight_threshold}

{pheno_type_str} phenotype used for SKAT.

Number of threads: {thread}

Output directory: {out_dir}

Output TWAS results file: {out_path}
********************************'''.format(
    **args.__dict__,
    pheno_type_str = {'C':'Continuous', 'D':'Dichotomous'}[args.phenotype_type],
    out_path = out_twas_path))

#############################################
# read in eQTL weights (ES) file

print('Reading weight file.')
# read in headers for Weight file
w_cols = tg.get_header(args.w_path)

# get the indices and dtypes for reading file into pandas
w_cols_ind, w_dtype = tg.weight_cols_dtype(w_cols, ['MAF','b','beta'], ['ES'])

# read in file in chunks
print('Reading eQTL weights file.')
Weight_chunks  = pd.read_csv(
    args.w_path, 
    sep='\t', 
    iterator=True, 
    chunksize=10000,
    usecols=w_cols_ind,
    dtype=w_dtype)

Weight = pd.concat([x[x['CHROM']==args.chr] for x in Weight_chunks]).reset_index(drop=True)

if Weight.empty:
    raise SystemExit('There are no valid eQTL weights.')

Weight = tg.optimize_cols(Weight)

if 'ID' in Weight.columns:
    Weight.rename(columns={'ID':'snpID'})

if not 'snpID' in Weight.columns:
    Weight['snpID'] = tg.get_snpIDs(Weight)

# sum additive and fixed effect
Weight['ES_sum'] = Weight['b'] + Weight['beta']

# filter SNPs with small abs effect size
Weight = Weight[abs(Weight.ES_sum) > args.weight_threshold]

# drop columns no longer needed
Weight = Weight.drop(columns=['CHROM','POS','REF','ALT'])

# read in gene annotation file
print('Reading gene annotation file.')
Gene_chunks = pd.read_csv(
    args.annot_path, 
    sep='\t', 
    iterator=True, 
    chunksize=10000,
    dtype={'CHROM':object,'GeneStart':np.int64,'GeneEnd':np.int64,'TargetID':object,'GeneName':object}, 
    usecols=['CHROM','GeneStart','GeneEnd','TargetID','GeneName'])

Gene = pd.concat([x[x['CHROM']==args.chr] for x in Gene_chunks]).reset_index(drop=True)

Gene = tg.optimize_cols(Gene)

# define target ids
TargetID = Gene.TargetID.values
n_targets = TargetID.size

# read in ped file
print('Reading PED, PED info files.')
PED = pd.read_csv(
    args.ped_path,
    sep='\t').rename(columns={'#FAM_ID':'FAM_ID'})
PED = tg.optimize_cols(PED)

# read in ped info
# P:phenotype
# C:covariate
Asso_Info = pd.read_csv(
    args.pedinfo_path,
    sep='\t',
    header=None,
    names=['Ind','Var'])

# get phenotype, covariate for SKAT test
phenotype = PED[Asso_Info[Asso_Info.Ind=='P'].Var]
covariate = PED[Asso_Info[Asso_Info.Ind=='C'].Var]

# read genotype file header
print('Reading genotype file header.\n')
g_cols = tg.call_tabix_header(args.geno_path)
gcol_sampleids = g_cols[gcol_sampleids_strt_ind:]

# get sampleids to use
# intersect samples in phenotype and genotype file
gcol_pheno_sampleids = np.intersect1d(PED.IND_ID, gcol_sampleids)

if not gcol_pheno_sampleids.size:
    raise SystemExit('The phenotype file and genotype file have no sampleIDs in common.')

# load sampleids
print('Reading sampleID file.\n')
spec_sampleids = pd.read_csv(
    args.sampleid_path,
    sep='\t',
    header=None)[0].drop_duplicates()

# intersect samples in phenotype, user-specified sampleids, and genotype file
print('Matching sampleIDs.\n')
sampleID = np.intersect1d(spec_sampleids, gcol_pheno_sampleids) 

n_samples = sampleID.size

if not n_samples:
    raise SystemExit('The phenotype file, genotype file, and sampleID file have no sampleIDs in common.')

# get genotype  columns, dtype to read in
g_cols_ind, g_dtype = tg.genofile_cols_dtype(g_cols, args.genofile_type, sampleID)

# prep output
print('Creating file: ' + out_twas_path + '\n')

out_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName','n_snps']

if args.phenotype_type == 'D' and n_samples < 2000:
    out_cols = out_cols + ['p_value_resampling', 'p_value_noadj']

else:
    out_cols = out_cols + ['p_value']

pd.DataFrame(columns=out_cols).to_csv(
    out_twas_path,
    sep='\t',
    header=True,
    index=None,
    mode='w')

print('********************************\n')

##########################################################

def thread_process(num):
    try:
        target = TargetID[num]
        print('num=' + str(num) + '\nTargetID=' + target)
        Gene_info = Gene.iloc[[num]].reset_index(drop=True)

        # make sure the weight file has SNPs for this target before processing genotype data
        target_weight = Weight[Weight.TargetID==target][['snpID', 'ES_sum', 'MAF']]
        target_weight = target_weight.drop_duplicates(['snpID'],keep='first') 

        if target_weight.empty:
            print('No test SNPs with non-zero cis-eQTL weights for TargetID: ' + target + '\n')
            return None        

        # READING IN GENOTYPE FILE WITHIN A THREAD_PROCESS
        start = str(max(int(Gene_info.GeneStart)-args.window,0))
        end = str(int(Gene_info.GeneEnd) + args.window)    

        g_proc_out = tg.call_tabix(args.geno_path, args.chr, start, end)

        if not g_proc_out:
            print('No test SNPs with GWAS Zscore for TargetID: ' + target + '\n')
            return None  

        target_geno = pd.read_csv(StringIO(g_proc_out.decode('utf-8')),
                sep='\t',
                low_memory=False,
                header=None,
                usecols=g_cols_ind,
                dtype=g_dtype)

        target_geno.columns = [g_cols[i] for i in target_geno.columns]

        # OPTIMIZE THE DATAFRAME
        target_geno = tg.optimize_cols(target_geno)

        # GET SNP IDS FOR A DATAFRAME
        target_geno['snpID'] = tg.get_snpIDs(target_geno)
        target_geno = target_geno.drop_duplicates(['snpID'],keep='first').reset_index(drop=True)
        
        # REFORMAT VCF FILES
        if args.genofile_type=='vcf':
            target_geno = tg.check_prep_vcf(target_geno, args.format, sampleID)

        # REFORMAT DATAFRAME VALUES - immediately after the reformat vcf step, used for both dosage and vcf flies
        target_geno[sampleID] = target_geno[sampleID].apply(lambda x:tg.reformat_sample_vals(x,args.format), axis=0)

        # CALCULATE AND FILTER MAF - only include MAF >= 0, snpIDs may be flipped so don't do final MAF filter
        target_geno = tg.calc_maf(target_geno, sampleID, 0, op=operator.gt)

        # CALCULATE AND FILTER pHWE
        target_geno = tg.calc_p_hwe(target_geno, sampleID, args.hwe)
            
        # MATCH SNPS - intersect SNPs from eQTL weight file and test genotype file
        snp_overlap = np.intersect1d(target_weight.snpID,target_geno.snpID)
            
        if not snp_overlap.size:
            print('No SNPs overlapped between eQTL weight file and test genotype file. Flipping REF and ALT for genotype file.')
            # GET FLIPPED SNP IDS:
            target_geno = target_geno.drop(columns=['snpID'])
            target_geno['snpID'] = tg.get_snpIDs(target_geno, flip=True)
            snp_overlap_flip = np.intersect1d(target_weight.snpID,target_geno.snpID)
                
            if not snp_overlap_flip.size:
                print('No SNPs overlapped between eQTL weight file and test genotype file after flipping for TargetID:' + target + '\n')
                return None

            snp_overlap = snp_overlap_flip
            target_weight['ES_sum'] = -target_weight['ES_sum'].astype('float')
            target_geno['MAF'] = 1-target_geno['MAF'].astype('float')
        
        # print('Number of SNPs overlapped between eQTL weight file and test genotype file before MAF filter: ' + str(snp_overlap.size) + '\n')

        target_geno = target_geno[np.append(sampleID,['snpID','MAF'])]

        # DO FINAL MAF FILTER - must be after possible flipping
        target_geno  = target_geno[target_geno['MAF'] > args.maf]

        # SNP OVERLAP AFTER FILTER
        snp_overlap = np.intersect1d(target_weight.snpID,target_geno.snpID)
        n_snps = snp_overlap.size

        if not n_snps:
            print('No overlapped SNPs exceed the specified MAF threshold of ' + str(args.maf) + ' for TargetID: ' + target + '\n')
            return None           

        print('Running TWAS.\nN SNPs=' + str(n_snps))

        # FILTER GENOTYPE AND WEIGHT DF BY SNPOVERLAP
        target_geno = target_geno[target_geno.snpID.isin(snp_overlap)]
        target_weight = target_weight[target_weight.snpID.isin(snp_overlap)]

        # MERGE DATAFRAMES
        geno_weight = target_geno.merge(
            target_weight,
            left_on='snpID',
            right_on='snpID',
            how='outer')

        # GENOTYPE MATRIX
        genotype_mat = geno_weight[sampleID].T.values

        # WEIGHT VECTOR
        weight_mat = geno_weight['ES_sum'].values

        # INITIALIZE OUTPUT DATAFRAME
        result = Gene_info.copy()
        result['n_snps'] = n_snps

        # SKAT TEST
        p_value = SKAT.SKAT(
            genotype_mat,
            phenotype,
            covariate,
            weight_mat,
            args.phenotype_type)

        if args.phenotype_type == 'C': 
            result['p_value'] = p_value
            result['TargetID'] = target

        elif args.phenotype_type == 'D':
            if len(p_value) == 2:
                result['p_value_resampling'] = p_value[0]
                result['p_value_noadj'] = p_value[1]
                result['TargetID'] = target

            else:
                result['p_value'] = p_value
                result['TargetID'] = target
                    
        result.to_csv(
            out_twas_path,
            sep='\t',
            index=None,
            header=None,
            mode='a')

        print('Target TWAS completed.\n')

    except Exception as e:
        e_info = sys.exc_info()

        e_type = e_info[0].__name__
        e_line = e_info[2].tb_lineno
        e_tracebk = ''.join(traceback.format_tb(e_info[2]))

        print('Caught a "{}" type exception for TargetID={}, num={} on line {}:\n  {}\nTraceback:\n{}'.format(
            e_type, target, num, e_line, e, e_tracebk))

    finally:
        # print info to log do not wait for buffer to fill up
        sys.stdout.flush()

##############################################################
# thread process
if __name__ == '__main__':
    print('Starting VC-TWAS for ' + str(n_targets) + ' target genes.\n')
    pool = multiprocessing.Pool(args.thread)
    pool.map(thread_process,[num for num in range(n_targets)])
    pool.close()
    pool.join()
    print('Done.')

############################################################
# time calculation
elapsed_sec = time() - start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print('Computation time (DD:HH:MM:SS): ' + elapsed_time)
