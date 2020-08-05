#!/usr/bin/env python

###################################################################
# Import packages needed
import argparse
import io
from io import StringIO
import multiprocessing
import subprocess
import sys
from time import time

import numpy as np
import pandas as pd

from scipy.stats import chi2

#########################################################################
### timecalculation
start_time = time()

############################################################
### variable needed
parser = argparse.ArgumentParser(description='Asso Study 02')

### Specify tool directory
parser.add_argument('--TIGAR_dir',type=str)

### Gene annotation file
parser.add_argument("--gene_anno",type=str,dest='annot_path')

### GWAS Z score file
parser.add_argument('--Zscore',type=str,dest='z_path')

### Header of GWAS Z score file
parser.add_argument('--Zscore_colnames',type=str,dest='zcol_path')

### Weight
parser.add_argument('--weight',type=str,dest='w_path')

### Reference covariance file
parser.add_argument('--LD',type=str,dest='ld_path')

### chromosome number
parser.add_argument('--chr',type=str)

### window
parser.add_argument('--window',type=float)

### Number of thread
parser.add_argument('--thread',type=int)

### Weight threshold to include SNP in TWAS
parser.add_argument('--weight_threshold',type=float)

# specify "FUSION" or "S_PrediXcan" Zscore test statistic to use
parser.add_argument('--test_stat',type=str,default='FUSION')

### Output dir
parser.add_argument('--out_dir',type=str)

args = parser.parse_args()

sys.path.append(args.TIGAR_dir)

##################################################
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

################################################################################################
### variable checking
print("Gene annotation file to specify the list of genes for TWAS : " + args.annot_path + "\n")
print("GWAS summary statistics Z-score file : " + args.z_path+ "\n")
print("cis-eQTL weight file : " + args.w_path + "\n")
print("Reference LD genotype covariance file: " + args.ld_path + "\n")
print("Chromosome number : " + args.chr + "\n")
print("Test gene region including SNPs within +- window = " + str(args.window) + " base pair of GeneStart/GeneEnd positions \n")
print("Number of threads : " + str(args.thread) + "\n")
print("SNP weight inclusion threshold : " + str(args.weight_threshold) + "\n")
print("Test statistic to use: " + args.test_stat + "\n")
print("Output directory : " + args.out_dir + "\n")


##################################################
### Read in gene annotation 
print("Reading gene annotation file.")
Gene_chunks = pd.read_csv(
    args.annot_path, 
    sep='\t', 
    iterator=True, 
    chunksize=10000,
    dtype={'CHROM':object,'GeneStart':np.int64,'GeneEnd':np.int64,'TargetID':object,'GeneName':object}, 
    usecols=['CHROM','GeneStart','GeneEnd','TargetID','GeneName'])

Gene = pd.concat([x[x['CHROM'] == args.chr] for x in Gene_chunks] ).reset_index(drop=True)

Gene = tg.optimize_cols(Gene)

TargetID = np.array(Gene.TargetID)
n_targets = TargetID.size

# read in headers for Weight and Zscore files
w_cols = tg.get_header(args.w_path, zipped=True)
z_cols = tg.get_header(args.z_path, zipped=True)

# get the indices and dtypes for reading files into pandas
w_cols_ind, w_dtype = tg.weight_cols_dtype(w_cols)
z_cols_ind, z_dtype = tg.zscore_cols_dtype(z_cols)


print("Creating file: "+'CHR'+args.chr+'_sumstat_assoc.txt')
out_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName','n_snps','Zscore','PVALUE']
pd.DataFrame(columns=out_cols).to_csv(
    args.out_dir+'/CHR'+args.chr+'_sumstat_assoc.txt',
    sep='\t',
    index=None,
    header=True,
    mode='w')


##################################################
def thread_process(num):
    try: 
        print("\nnum="+str(num)+"\nTargetID="+TargetID[num])
        Gene_info = Gene.iloc[[num]]

        # get start and end positions to tabix
        start = str(max(int(Gene_info.GeneStart)-args.window,0))
        end = str(int(Gene_info.GeneEnd)+args.window)

        # tabix Weight file
        w_proc_out = tg.call_tabix(args.w_path, args.chr, start, end)

        if not w_proc_out:
            print("No test SNPs with non-zero cis-eQTL weights in window of gene="+TargetID[num])
            return None

        # tabix Zscore file
        z_proc_out = tg.call_tabix(args.z_path, args.chr, start, end)

        if not z_proc_out:
            print("No test SNPs with GWAS Zscore for gene="+TargetID[num])
            return None

        # parse tabix output for Weight, filtered by TargetID[num], threshold
        Weight_chunks = pd.read_csv(
            StringIO(w_proc_out.decode('utf-8')),
            sep='\t',
            header=None,
            low_memory=False,
            iterator=True, 
            chunksize=10000,
            usecols=w_cols_ind,
            dtype=w_dtype)

        Weight = pd.concat([x[ (x[w_cols_ind[4]]==TargetID[num]) & (abs(x[w_cols_ind[5]]) > args.weight_threshold )] for x in Weight_chunks]).reset_index(drop=True)

        if Weight.empty:
            print("No test SNPs with cis-eQTL weights with magnitude that exceeds specified threshold for gene="+TargetID[num])
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
            print("No overlapping test SNPs that have magnitude of cis-eQTL weights greater than threshold value and with GWAS Zscore for gene="+TargetID[num])
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
            right_on='snpID').reset_index(drop=True)

        snp_search_ids = ZW.snpID.values

        ### Read in reference covariance matrix file by snpID
        MCOV_chunks = pd.read_csv(
            args.ld_path, 
            sep='\t', 
            compression='gzip', 
            iterator=True, 
            chunksize=10000,
            usecols=['snpID','COV'], 
            dtype={'snpID': object, 'COV': object})

        MCOV = pd.concat([x[x.snpID.isin(snp_search_ids)] for x in MCOV_chunks]).drop_duplicates(['snpID'], keep='first')

        if MCOV.empty:
          print("No reference covariance information for target SNPs for gene="+TargetID[num])
          return None
        
        MCOV['COV'] = MCOV['COV'].apply(lambda x:np.array(x.split(',')).astype('float'))

        # construct covariance matrix
        inds = MCOV.index
        n_inds = inds.size
        V_upper = np.zeros((n_inds,n_inds))
        
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
        
        # filter ZW to include only snpIDs also in MCOV
        ZW = ZW[ZW.snpID.isin(MCOV.snpID)]
        n_snps = str(ZW.snpID.size)

        print("TWAS for gene="+TargetID[num])
        print("N SNPs="+n_snps)

        ### Calculate burden Z score
        z_denom = np.sqrt(np.linalg.multi_dot([ZW.ES.values, V, ZW.ES.values]))

        if args.test_stat == 'FUSION':
            z_numer = np.vdot(ZW.Zscore.values, ZW.ES.values)

        elif args.test_stat == 'S_PrediXcan':
            snp_SD = np.sqrt(snp_Var)
            z_numer = snp_SD.dot(ZW.ES.values * ZW.Zscore.values)

        burden_Z = z_numer/z_denom
        
        if np.isnan(burden_Z):
        	print("Could not calculate burden_Z: NaN value.")
        	return None

        ### p-value for chi-square test
        pval = 1-chi2.cdf(burden_Z**2,1)

        ### create output dataframe
        result = Gene_info.copy()
        result['n_snps'] = n_snps
        result['TWAS_Zscore'] = burden_Z
        result['PVALUE'] = pval

        ### write to file
        result.to_csv(
            args.out_dir+'/CHR'+args.chr+'_sumstat_assoc.txt',
            sep='\t',
            index=None,
            header=None,
            mode='a',
            float_format='%.4f')

        print('TWAS complete.')

    except Exception as e:
        e_type, e_obj, e_tracebk = sys.exc_info()
        e_line_num = e_tracebk.tb_lineno

        e, e_type, e_line_num = [str(x) for x in [e, e_type, e_line_num]]
        
        print('Caught a type '+ e_type +' exception for TargetID='+TargetID[num]+', num=' + str(num) + ' on line '+e_line_num+':\n' + e )

    finally:
        # print info to log do not wait for buffer to fill up
        sys.stdout.flush()

###############################################################
### thread process

if __name__ == '__main__':
    print("Starting TWAS for "+str(n_targets)+" target genes.")
    pool = multiprocessing.Pool(args.thread)
    pool.map(thread_process,[num for num in range(n_targets)])
    pool.close()
    pool.join()
    print("Done.")


############################################################
### time calculation
elapsed_sec = time()-start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print("Computation time (DD:HH:MM:SS): " + elapsed_time)







