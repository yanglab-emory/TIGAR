# %%
#!/usr/bin/env python

###################################################################
# Import packages needed
import argparse
from io import StringIO
import multiprocessing
import subprocess
import sys
from time import time
import numpy as np
import pandas as pd
import SKAT
import TIGARutils as tg


# %%
###############################################################
# time calculation
start_time = time()

###############################################################
# parse input arguments
parser = argparse.ArgumentParser(description='VC_TWAS summary statistics')
### Gene annotation file
parser.add_argument("--gene_anno",type=str,dest='annot_path')

### GWAS result file
parser.add_argument('--GWAS_result',type=str,dest='gwas_path')

### Weight
parser.add_argument('--Weight',type=str,dest='w_path')

### sample size 
parser.add_argument('--sample_size',type=int)

### weight threshold
parser.add_argument('--weight_threshold',type=float)

### Reference covariance file
parser.add_argument('--LD',type=str,dest='ld_path')

### chromosome number
parser.add_argument('--chr',type=str)

### window
parser.add_argument('--window',type=int)

### Number of thread
parser.add_argument('--thread',type=int)

### Output dir
parser.add_argument('--out_dir',type=str)

args = parser.parse_args()
#sys.path.append(args.TIGAR_dir)
#sys.path.append(args.TIGAR_dir + '/VC_TWAS')

#############################################################
# Print input arguments to log
out_sum_VCTWAS_path = args.out_dir + '/CHR' + args.chr + '_sum_VC_TWAS.txt'

print(
'''********************************
Input Arguments

Gene annotation file specifying genes for VC-TWAS: {annot_path}

GWAS summary statistics file: {gwas_path}

cis-eQTL weight file: {w_path}

sample—size:{sample_size}

SNP weight inclusion threshold:{weight_threshold}

Reference LD genotype covariance file: {ld_path}

Chromosome: {chr}

Gene training region SNP inclusion window: +-{window}

Number of threads: {thread}

Output directory: {out_dir}

Output TWAS results file: {out_path}
********************************'''.format(
    **args.__dict__,
    out_path = out_sum_VCTWAS_path))

###############################################################
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

###############################################################
### Read in gene annotation 
print('Reading gene annotation file.')
Gene_chunks = pd.read_csv(
    args.annot_path, 
    sep='\t', 
    iterator=True, 
    chunksize=10000,
    dtype={'CHROM':object,'GeneStart':np.int64,'GeneEnd':np.int64,'TargetID':object,'GeneName':object}, 
    usecols=['CHROM','GeneStart','GeneEnd','TargetID','GeneName'])
Gene = pd.concat([x[x['CHROM']==args.chr] for x in Gene_chunks] ).reset_index(drop=True)

TargetID = np.array(Gene.TargetID)
n_targets = TargetID.size

# read in headers for Weight, GWAS, MCOV result files
w_cols = tg.get_header(args.w_path, zipped=True)
gwas_cols = tg.get_header(args.gwas_path, zipped=True)
mcov_cols = tg.get_header(args.ld_path, zipped=True)

# get the indices and dtypes for reading files into pandas
w_cols_ind, w_dtype = tg.weight_cols_dtype(w_cols, ['MAF','b','beta'], ['ES'])
gwas_cols_ind, gwas_dtype = tg.gwas_cols_dtype(gwas_cols)
MCOV_cols_ind, MCOV_dtype = tg.MCOV_cols_dtype(mcov_cols)


# PREP OUTPUT - print output headers to files
print('Creating file: ' + out_sum_VCTWAS_path + '\n')
out_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName','Pvalue']
pd.DataFrame(columns=out_cols).to_csv(
    out_sum_VCTWAS_path,
    sep='\t',
    index=None,
    header=True,
    mode='w')

# %%
###############################################################
# thread function
def thread_process(num):
    try: 
        target = TargetID[num]
        print('num=' + str(num) + '\nTargetID=' + target)
        Gene_info = Gene.iloc[[num]].reset_index(drop=True)

        # get start and end positions to tabix
        start = str(max(int(Gene_info.GeneStart)-args.window,0))
        end = str(int(Gene_info.GeneEnd)+args.window)
        # tabix Weight file
        # print('Reading weight data.')
        w_proc_out = tg.call_tabix(args.w_path, args.chr, start, end)

        if not w_proc_out:
            print('No test SNPs with non-zero cis-eQTL weights for TargetID: ' + target + '\n')
            return None

        # parse tabix output for Weight, filtered by target
        Weight_chunks = pd.read_csv(
            StringIO(w_proc_out.decode('utf-8')),
            sep='\t',
            header=None,
            low_memory=False,
            iterator=True, 
            chunksize=10000,
            usecols=w_cols_ind,
            dtype=w_dtype)
        Weight = pd.concat([x[ (x[w_cols_ind[4]]==target)] for x in Weight_chunks]).reset_index(drop=True)

        Weight.columns = [w_cols[i] for i in tuple(Weight.columns)]
        # add random effect and fixed effect together
        Weight["ES_sum"]=Weight["beta"]+Weight["b"]
        # parse tabix output for Weight, filtered by threshold
        Weight = Weight[abs(Weight.ES_sum) > args.weight_threshold]

        if Weight.empty:
            print('No test SNPs with cis-eQTL weights with magnitude that exceeds specified weight threshold for TargetID: ' + target + '\n')
            return None
    
        # 'ID' snpID column does not exist in TIGAR training output from previous versions
        # check to see if it is in the file, if not will need to generate snpIDs later
        if 'ID' in Weight.columns:
            Weight = Weight.rename(columns={'ID':'snpID'})

        if not 'snpID' in Weight.columns:
            Weight['snpID'] = tg.get_snpIDs(Weight)
            Weight =Weight.drop_duplicates(['snpID'],keep='first').reset_index(drop=True)
        
        
        # tabix gwas result file
        # print('Reading gwas data.')
        gwas_proc_out = tg.call_tabix(args.gwas_path, args.chr, start, end)

        if not gwas_proc_out:
            print('No test SNPs with GWAS Result for TargetID: ' + target + '\n')
            return None
        # parse tabix output for gwas
        GWAS_result = pd.read_csv(
            StringIO(gwas_proc_out.decode('utf-8')),
            sep='\t',
            header=None,
            low_memory=False,
            usecols=gwas_cols_ind,
            dtype=gwas_dtype)

        GWAS_result.columns = [gwas_cols[i] for i in tuple(GWAS_result.columns)]
        #check ID 
        if 'ID' in GWAS_result.columns:
            GWAS_result = GWAS_result.rename(columns={'ID':'snpID'})
        
        if not 'snpID' in GWAS_result.columns:
            GWAS_result['snpID'] = tg.get_snpIDs(GWAS_result)
            GWAS_result =GWAS_result.drop_duplicates(['snpID'],keep='first').reset_index(drop=True)
        
        #get flip snpID
        GWAS_result['snpID_flip'] = tg.get_snpIDs(GWAS_result, flip=True)  
    
        # check for overlapping SNPs in Weight, gwas data
        snp_overlap_orig = np.intersect1d(np.array(Weight.snpID), np.array(GWAS_result.snpID))
        snp_overlap_flip = np.intersect1d(Weight.snpID, np.array(GWAS_result.snpID_flip))
        snp_overlap = np.concatenate((snp_overlap_orig, snp_overlap_flip))   

        # filter dataframes by overlapping SNPs
        Weight = Weight[Weight.snpID.isin(snp_overlap)]
        GWAS_result = GWAS_result[GWAS_result.snpID.isin(snp_overlap_orig) | GWAS_result.snpID_flip.isin(snp_overlap_flip)]

     
        #handle flip
        if (snp_overlap_orig.size > 0) and (snp_overlap_flip.size > 0):
            GWAS_result['snpID'], GWAS_result['BETA'] = handle_flip(GWAS_result,'snpID','snpID_flip','BETA',snp_overlap_orig, snp_overlap_flip)
         elif snp_overlap_orig.size == snp_overlap.size:   
            GWAS_result['snpID'], GWAS_result['BETA'] = GWAS_result['snpID'], GWAS_result['BETA']
        else:
            GWAS_result['snpID'], GWAS_result['BETA'] = GWAS_result['snpID_flip'], -GWAS_result['BETA']
        

        GW = Weight.merge(GWAS_result[['snpID','BETA','SE']], 
            left_on='snpID', 
            right_on='snpID').drop_duplicates(['snpID'], keep='first').reset_index(drop=True)
        GW = GW.sort_values(by='POS')
            snp_search_ids = GW.snpID


        # Read in reference covariance matrix file 
        # print('Reading reference covariance data.')
        MCOV_proc_out = tg.call_tabix(args.ld_path, args.chr, start, end)
        MCOV = pd.read_csv(
                    StringIO(MCOV_proc_out.decode('utf-8')),
                    sep='\t',
                    header=None,
                    low_memory=False,
                    usecols=MCOV_cols_ind,
                    dtype=MCOV_dtype)
        MCOV.columns = [mcov_cols[i] for i in tuple(MCOV.columns)]

        ###check SNPID
        if 'ID' in MCOV.columns:
            MCOV = MCOV.rename(columns={'ID':'snpID'})

        if not 'snpID' in MCOV.columns:
            MCOV['snpID'] = tg.get_snpIDs(MCOV)


        if MCOV.empty:
          print('No reference covariance information for target SNPs for TargetID: ' + target + '\n')
          return None
        
        MCOV = MCOV.drop_duplicates(['snpID'],keep='first').reset_index(drop=True)
        MCOV['COV'] = MCOV['COV'].apply(lambda x:np.array(x.split(',')).astype('float'))
        MCOV = MCOV.sort_values(by='POS')

        # overlap mcov weight gwas, final snp ID
        snp_overlap_final= np.intersect1d(np.array(snp_search_ids),np.array(MCOV.snpID))

        if len(snp_overlap_final) == 0:
            print("No overlap snp in weight file, cov matrix & gwas result:"+TargetID[num])
            return None
        inds = MCOV[MCOV.snpID.isin(snp_overlap_final)].index
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
        D= np.diag(V)
        # filter GW to include only snpIDs also in MCOV
        GW = GW[GW.snpID.isin(snp_overlap_final)]
        #pepare input for Summary statistics
        beta_var = pow(GW["SE"].values,2)
        beta_estimate = GW["BETA"].values
        weight_temp = GW["ES_sum"].values
        ###VC-TWAS SUMMARY
        p_val =SKAT.SKAT_summary(beta_var, beta_estimate, weight_temp, args.sample_size, V, D)

        ### create output dataframe
        result = Gene_info.copy()    
        result['Pvalue'] = p_val      
        
        # write to file
        result.to_csv(
            out_sum_VCTWAS_path,
            sep='\t',
            index=None,
            header=None,
            mode='a')

        print('Target VC-TWAS summary completed.\n')

    except Exception as e:
        e_type, e_obj, e_tracebk = sys.exc_info()
        e_line_num = e_tracebk.tb_lineno

        print('Caught a type {} exception for TargetID={}, num={} on line {}:\n{}\n'.format(e_type, target, num, e_line_num, e))

    finally:
        # print info to log do not wait for buffer to fill up
        sys.stdout.flush()

###############################################################

# %%
###############################################################
# thread process
if __name__ == '__main__':
    print('Starting VC-TWAS summary statistics with GWAS result for ' + str(n_targets) + ' target genes.\n')
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

