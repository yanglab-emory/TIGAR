#!/usr/bin/env python

###################################################################
# Import packages needed
import argparse
import time
import subprocess
import pandas as pd
import numpy as np
import io
from io import StringIO
from scipy.stats import chi2
import sys
import multiprocessing

#########################################################################
### timecalculation
start_time=time.clock()

############################################################
### variable needed
parser = argparse.ArgumentParser(description='Help: ')

### Gene annotation file
parser.add_argument("--gene_anno",type=str,default=None,dest='annot_path')

### GWAS Z score file
parser.add_argument('--Zscore',type=str,default=None,dest='z_path')

### Header of GWAS Z score file
parser.add_argument('--Zscore_colnames',type=str,default=None,dest='zcol_path')

### Weight
parser.add_argument('--weight',type=str,default=None,dest='w_path')

### Header of Weight file
parser.add_argument('--weight_colnames',type=str,default=None,dest='wcol_path')

### Reference covariance file
parser.add_argument('--LD',type=str,default=None,dest='ld_path')

### chromosome number
parser.add_argument('--chr',type=str,default=None)

### window
parser.add_argument('--window',type=float,default=None)

### Number of thread
parser.add_argument('--thread',type=int,default=None)

### Weight threshold to include SNP in TWAS
parser.add_argument('--threshold',type=float,default=None)

### Output dir
parser.add_argument('--out_dir',type=str,default=None)

args = parser.parse_args()

################################################################################################
### variable checking
print("Gene annotation file to specify the list of genes for TWAS : "+args.annot_path + "\n")
print("GWAS summary statistics Z-score file : " + args.z_path+ "\n")
print("cis-eQTL weight file : " + args.w_path + "\n")
print("Reference LD genotype covariance file:"+args.ld_path + "\n")
print("Chromosome number : "+args.chr+ "\n")
print("Test gene region including SNPs within +- window = "+str(args.window) + " base pair of GeneStart/GeneEnd positions \n")
print("Number of threads : "+str(args.thread) + "\n")
print("SNP weight inclusion threshold : "+str(args.threshold) + "\n")
print("Output directory : " + args.out_dir + "\n")

##################################################
# Call tabix, read in lines into byt array
def call_tabix(path, chr, start, end):
    proc = subprocess.Popen(
        ["tabix "+path+" "+chr+":"+start+"-"+end],
        shell=True,
        stdout=subprocess.PIPE,
        bufsize=1)
    proc_out = bytearray()
    while proc.poll() is None:
        line =  proc.stdout.readline()
        if len(line) == 0:
            break
        proc_out += line
    return proc_out

# Determine index of required columns,
# assigns correct dtype to correct index
# 'ID' column does not exist in TIGAR training output from previous versions
# check to see if it is in the file, if not will need to generate snpIDs later
def default_cols_dtype(file_cols, df_name):
  df_dict = {
  'Weight': {'cols': ['CHROM','POS','REF','ALT','TargetID','ES'], 
             'dtype': [object,np.int64,object,object,object,np.float64]},
  'Zscore': {'cols': ['CHROM','POS','REF','ALT','Zscore'], 
             'dtype': [object,np.int64,object,object,np.float64]}
  }

  if (df_name =='Weight') & ('ID' in file_cols):
    df_dict['Weight']['cols'].append('ID')
    df_dict['Weight']['dtype'].append('object')

  file_cols_ind = tuple(map(lambda col: file_cols.index(col), df_dict[df_name]['cols']))
  file_dtype = {c:d for c,d in zip(file_cols_ind, df_dict[df_name]['dtype'])}

  return file_cols_ind, file_dtype

# Decrease memory by downcasting 'CHROM' column to integer, integer and float columns to minimum size that will not lose info
def optimize_cols(df: pd.DataFrame):
  if 'CHROM' in df.columns:
    df['CHROM'] = df['CHROM'].astype(str).astype(int)

  ints = df.select_dtypes(include=['int64']).columns.tolist()
  df[ints] = df[ints].apply(pd.to_numeric, downcast='integer')

  floats = df.select_dtypes(include=['float64']).columns.tolist()
  df[floats] = df[floats].apply(pd.to_numeric, downcast='float')

  return df

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

def get_snpIDs(df: pd.DataFrame, flip=False):
    chroms = df['CHROM'].astype('str').values
    pos = df['POS'].astype('str').values
    ref = df['REF'].values
    alt = df['ALT'].values
    if flip:
        return [':'.join(i) for i in zip(chroms,pos,alt,ref)]
    else:
        return [':'.join(i) for i in zip(chroms,pos,ref,alt)]  
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

Gene = optimize_cols(Gene)

TargetID = np.array(Gene.TargetID)
n_targets = len(TargetID)

# read in headers for Weight and Zscore files
w_cols = tuple(pd.read_csv(args.wcol_path,sep='\t'))
z_cols = tuple(pd.read_csv(args.zcol_path,sep='\t'))

# get the indices and dtypes for reading files into pandas
w_cols_ind, w_dtype = default_cols_dtype(w_cols,'Weight')
z_cols_ind, z_dtype = default_cols_dtype(z_cols,'Zscore')

print("Creating data frame:"+'CHR'+args.chr+'_sumstat_assoc.txt')
pd.DataFrame(columns=['CHROM','GeneStart','GeneEnd','TargetID','GeneName','Zscore','PVALUE']).to_csv(args.out_dir+'/CHR'+args.chr+'_sumstat_assoc.txt',
                     sep='\t',index=None,header=True,mode='w')

def thread_process(num):
    try: 
        print("\nnum="+str(num)+"\nTargetID="+TargetID[num])
        Gene_info = Gene.iloc[[num]]

        # get start and end positions to tabix
        start = str(int(Gene_info.GeneStart)-args.window)
        end = str(int(Gene_info.GeneEnd)+args.window)

        # tabix Weight file
        w_proc_out = call_tabix(args.w_path, args.chr, start, end)

        if not w_proc_out:
            print("No test SNPs with non-zero cis-eQTL weights in window of gene="+TargetID[num])
            return None

        # tabix Zscore file
        z_proc_out = call_tabix(args.z_path, args.chr, start, end)

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

        Weight = pd.concat([x[ (x[w_cols_ind[4]]==TargetID[num]) & (abs(x[w_cols_ind[5]]) > args.threshold )  ] for x in Weight_chunks]).reset_index(drop=True)

        if Weight.empty:
            # print("No test SNPs with non-zero cis-eQTL weights for gene="+TargetID[num])
            print("No test SNPs with cis-eQTL weights with magnitude that exceeds specified threshold for gene="+TargetID[num])
            return None

        Weight.columns = [w_cols[i] for i in tuple(Weight.columns)]
        Weight = optimize_cols(Weight)

        # Weight['snpID'] = (Weight['CHROM'].astype('str')
        #     +':'+Weight['POS'].astype('str')
        #     +':'+Weight.REF
        #     +':'+Weight.ALT)


        # 'ID' snpID column does not exist in TIGAR training output from previous versions
        # check to see if it is in the file, if not will need to generate snpIDs later
        if 'ID' in Weight.columns:
            Weight.rename(columns={'ID':'snpID'})
        else:
            Weight['snpID'] = get_snpIDs(Weight)

        # parse tabix output for Zscore
        Zscore = pd.read_csv(
            StringIO(z_proc_out.decode('utf-8')),
            sep='\t',
            header=None,
            low_memory=False,
            usecols=z_cols_ind,
            dtype=z_dtype)

        Zscore.columns = [z_cols[i] for i in tuple(Zscore.columns)]
        Zscore = optimize_cols(Zscore)

        # Zscore['IDorig'] = (Zscore['CHROM'].astype('str')
        #     +':'+Zscore['POS'].astype('str')
        #     +':'+Zscore.REF
        #     +':'+Zscore.ALT)
        Zscore['IDorig'] = get_snpIDs(Zscore)

        # Zscore['IDflip'] = (Zscore['CHROM'].astype('str')
        #     +':'+Zscore['POS'].astype('str')
        #     +':'+Zscore.ALT
        #     +':'+Zscore.REF)
        Zscore['IDflip'] = get_snpIDs(Zscore, flip=True)

        # check for overlapping SNPs in Weight, Zscore data
        snp_overlap_orig = np.intersect1d(Weight.snpID, Zscore.IDorig)
        snp_overlap_flip = np.intersect1d(Weight.snpID, Zscore.IDflip)
        snp_overlap = np.concatenate((snp_overlap_orig, snp_overlap_flip))

        if snp_overlap.size == 0:
            print("No overlapping test SNPs with both non-zero cis-eQTL weights and with GWAS Zscore for gene="+TargetID[num])
            return None

        # filter dataframes by overlapping SNPs
        Weight = Weight[Weight.snpID.isin(snp_overlap)]
        Zscore = Zscore[Zscore.IDorig.isin(snp_overlap_orig) | Zscore.IDflip.isin(snp_overlap_flip)]
        Zscore = Zscore.drop(['CHROM','POS','REF','ALT'], axis=1)

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
        n_inds = len(inds)
        V_upper = np.zeros((n_inds,n_inds))
        
        for i in range(n_inds):
            cov_i = MCOV.COV.loc[inds[i]]
            N = len(cov_i)
            
            for j in range(i,n_inds):
                if inds[j] - inds[i] < N:
                    V_upper[i,j] = cov_i[inds[j]-inds[i]]
                else:
                    V_upper[i,j] = 0
                     
        V=V_upper+V_upper.T-np.diag(V_upper.diagonal())
        
        # filter ZW to include only snpIDs also in MCOV
        ZW = ZW[ZW.snpID.isin(MCOV.snpID)]
        n_snps = str(len(ZW.snpID))

        print("TWAS for gene="+TargetID[num])
        print("N SNPs="+n_snps)

        ### Calculate burden Z score
        burden_Z = np.asscalar(np.mat(ZW.Zscore)*np.mat(ZW.ES).T/np.sqrt(np.mat(ZW.ES)*V*np.mat(ZW.ES).T))
        
        if np.isnan(burden_Z):
        	print("Could not calculate burden_Z: NaN value.")
        	return None

        ### p-value for chi-square test
        pval = 1-chi2.cdf(burden_Z**2,1)

        ### create output row
        result = Gene_info.copy()
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
        print('Caught exception for TargetID='+TargetID[num]+', num='+str(num)+':' )
        print(e)

    finally:
        # print info to log do not wait for buffer to fill up
        sys.stdout.flush()


###############################################################
### thread process

if __name__ == '__main__':
    print("Starting TWAS for "+str(n_targets)+" target genes.")
    pool = multiprocessing.Pool(args.thread)
    pool.imap(thread_process,[num for num in range(n_targets)])
    pool.close()
    pool.join()
    print("Done.")


############################################################
### time calculation
time=round((time.clock()-start_time)/60,2)

print("Time:"+str(time)+' minutes')







