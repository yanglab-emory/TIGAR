#/usr/bin/env python

#############################################################
# Import packages needed
import argparse
import io
from io import StringIO
import multiprocessing
import operator
import subprocess
import sys
from time import time

import numpy as np
import pandas as pd

import scipy.stats as stats
from sklearn.model_selection import KFold
import statsmodels.api as sm

#############################################################
### time calculation
start_time = time()

#############################################################
### variables need
parser = argparse.ArgumentParser(description='DPR Training')

### Specify tool directory
parser.add_argument('--TIGAR_dir',type=str)

### for Gene annotation and Expression level file
parser.add_argument('--gene_exp',type=str,dest='geneexp_path')

### for traininread
parser.add_argument('--thread',type=int)

### training vcf fileg sampleID
parser.add_argument('--train_sampleID',type=str,dest='sampleid_path')

### specified chromosome number
parser.add_argument('--chr',type=str)

### Number of Th
parser.add_argument('--genofile',type=str,dest='geno_path')

### specified input file type(vcf or dosages)
parser.add_argument('--genofile_type',type=str)

### 'DS' or 'GT'
parser.add_argument('--format',type=str)

### p-value for HW test
parser.add_argument('--hwe',type=float)

### maf
parser.add_argument('--maf',type=float)

### window
parser.add_argument('--window',type=int)

### cvR2
parser.add_argument('--cvR2',type=int)

### model to run DPR
# Bayesian inference algorithm used by DPR: "1" (Variational Bayesian) or "2" (MCMC)
parser.add_argument('--dpr',type=str)

### define effect-size
# Output effect size type: "fixed" (default) for fixed effects or "additive" for an addition of fixed and random effects)
parser.add_argument('--ES', type=str)

### output dir
parser.add_argument('--out_dir',type=str)

args = parser.parse_args()

sys.path.append(args.TIGAR_dir)

DPR_path = args.TIGAR_dir + '/Model_Train_Pred/DPR'
#############################################################
# DEFINE, IMPORT FUNCTIONS
import TIGARutils as tg

# preps dpr input files, runs DPR, reads in dpr output
def prep_call_dpr(bimbam_df, pheno_df, snpannot_df, dpr_file_dir, targetid):
    ## PATHS FOR DPR INPUT
    bimbam_pth = dpr_file_dir + targetid + '_bimbam.txt'
    pheno_pth = dpr_file_dir + targetid + '_pheno.txt'
    snpannot_pth = dpr_file_dir + targetid + '_snp_annot.txt'

    ## OUTPUT FILES FOR DPR INPUT
    bimbam_df.to_csv(
        bimbam_pth,
        header=False,
        index=None,
        sep='\t',
        mode='w')

    pheno_df.to_csv(
        pheno_pth,
        header=False,
        index=None,
        sep='\t',
        mode='w')

    snpannot_df.to_csv(
            snpannot_pth,
            header=False,
            index=None,
            sep='\t',
            mode='w')

    ## CALL DPR
    try:
        DPR_call_args = [DPR_path, 
            '-g', bimbam_pth, 
            '-p', pheno_pth, 
            '-a', snpannot_pth, 
            '-dpr', args.dpr, 
            '-o', 'DPR_'+targetid]

        subprocess.check_call(DPR_call_args, cwd=dpr_file_dir)

    except subprocess.CalledProcessError as err:
        raise err

    ## READ IN AND PROCESS DPR OUTPUT
    dpr_out = pd.read_csv(
            dpr_file_dir + 'output/DPR_'+targetid+'.param.txt',
            sep='\t',
            header=0,
            names=['CHROM','snpID','POS','n_miss','b','beta','gamma'],
            usecols=['snpID','b','beta'],
            dtype={'snpID': object, 'b': np.float64, 'beta': np.float64})

    dpr_out = tg.optimize_cols(dpr_out)

    ### GET EFFECT SIZE
    if args.ES == 'fixed':
        dpr_out['ES'] = dpr_out['beta']

    elif args.ES == 'additive':
        dpr_out['ES'] = dpr_out['beta'] + dpr_out['b']

    return dpr_out


# calculated r2 of prediction based on out_weights_df, genotype data in bimbam_test_df , actual values pheno_test_df
def calc_r2(out_weights_df, bimbam_test_df, pheno_test_df, cv=False):

    # filter by snp overlap
    snp_overlap = np.intersect1d(out_weights_df.snpID, bimbam_test_df.snpID)
    out_weights_df = out_weights_df[out_weights_df.snpID.isin(snp_overlap)]
    bimbam_test_df = bimbam_test_df[bimbam_test_df.snpID.isin(snp_overlap)]

    # genotype, weight data for prediction
    test_geno_weights = out_weights_df.merge(
        bimbam_test_df,
        left_on='snpID',
        right_on='snpID',
        how='outer').set_index(['snpID'])
    test_geno_weights = tg.optimize_cols(test_geno_weights)

    snp_weights = test_geno_weights['ES']

    if cv:
        snp_weights = snp_weights.fillna(0)

    test_geno = test_geno_weights.drop(columns=['ES']).T

    ### calculate predicted value for test set
    target_pred = np.dot(test_geno, snp_weights)

    lm = sm.OLS(pheno_test_df.values, sm.add_constant(target_pred)).fit()

    if cv:
        return lm.rsquared

    # else, return Pvalue, R2 for final training
    return lm.f_pvalue, lm.rsquared


# function to do the ith cross validation step
def do_cv(i, target, target_geno_df, target_exp_df, snp_annot_df, cv_trainID, cv_testID, ):
    dpr_file_dir_cv = args.out_dir + '/CV_Files/'
    target_cv = target +'_CV'+str(i+1)

    trainID = cv_trainID[i]
    testID = cv_testID[i]

    bimbam_train = target_geno_df[np.concatenate((['snpID','REF','ALT'],trainID))]

    pheno_train = target_exp_df[trainID].T

    ### PREP INPUT, CALL DPR
    try:
        dpr_out_cv = prep_call_dpr(
            bimbam_train, 
            pheno_train, 
            snp_annot_df,
            dpr_file_dir_cv, 
            target_cv)

    except subprocess.CalledProcessError as err:
        print("DPR failed in CV"+str(i+1)+" for TargetID: "+target)
        return 0
    
    ### for R2 calculation
    out_weights_cv = dpr_out_cv[['snpID','ES']]
    bimbam_test = target_geno_df[np.append(['snpID'],testID)]
    pheno_test = target_exp_df[testID].T

    cv_rsquared = calc_r2(
        out_weights_cv, 
        bimbam_test, 
        pheno_test, 
        cv=True)

    # RETURN R2 RESULT
    return(cv_rsquared)


#############################################################
### check input arguments
print("********************************\n   Imput Arguments\n********************************\n")

print("Gene Annotation and Expression file: "+args.geneexp_path + "\n")

if args.sampleid_path:
    print("Training sampleID file: "+args.sampleid_path+ "\n")

print("Chromosome number: "+args.chr+ "\n")
print("Number of threads: "+str(args.thread)+ "\n")
print("Training genotype file: "+args.geno_path+ "\n")
# print("Column names of genotype file:"+args.gcol_path+ "\n")

if args.genofile_type=='vcf':
    print("Genotype file used for training is VCF type with format: " + args.format + "\n")
    gcol_sampleids_strt_ind = 9

elif args.genofile_type=='dosage':
    print("Genotype file used for training is dosage type."+ "\n")
    args.format = 'DS'
    gcol_sampleids_strt_ind = 5

else:
    raise SystemExit("Please specify input genotype file type (--genofile_type) as either 'vcf' or 'dosage'."+ "\n")

print("Gene region size: window ="+str(args.window)+ "\n")
print("Evaluate DPR model by 5-fold cross validation: cvR2 = "+str(args.cvR2) + "\n")
print("Threshold for MAF: "+str(args.maf)+ "\n")
print("Threshold for HWE p-value: "+str(args.hwe)+ "\n")
print("Runing DPR Model: dpr="+args.dpr+ "\n")
print("Output Effect-size type: "+args.ES+ "\n")
print("Output directory: "+args.out_dir+ "\n")

out_train_weight_path = args.out_dir+'/CHR'+args.chr+'_DPR_train_eQTLweights.txt'
print("Training weights output file: " + out_train_weight_path +"\n")

out_train_info_path = args.out_dir+'/CHR'+args.chr+'_DPR_train_GeneInfo.txt'
print("Training info file: " + out_train_info_path +"\n")

print("********************************\n\n")
#############################################################
# Prepare DPR input

### Read in Gene annotation and Expression level file (text file)
### First five columns should be fixed:
### 1.Chrmosome Number
### 2.GeneStart Posistion
### 3.GeneEnd Position
### 4.TargetID (i.e.GeneID, treated as unique annotation for each gene)
### 5.Gene Name
# Gene Expression header, sampleIDs
exp_cols = tg.get_header(args.geneexp_path)
exp_sampleids = exp_cols[5:]

# genofile header, sampleIDs
g_cols = tg.call_tabix_header(args.geno_path)
gcol_sampleids = g_cols[gcol_sampleids_strt_ind:]

# geno, exp overlapping sampleIDs
gcol_exp_sampleids = np.intersect1d(gcol_sampleids, exp_sampleids)

if not gcol_exp_sampleids.size:
    raise SystemExit("The gene expression file and genotype file have no sampleIDs in common.")

# if user specified path with sampleids, and at least one sampid from that file, get intersection; otherwise use the overlapping sampleIDs in the genotype and expression file
if args.sampleid_path:
    spec_sampleids = pd.read_csv(
        args.sampleid_path,
        sep='\t',
        header=None)

    spec_sampleids = spec_sampleids[0].drop_duplicates()

    sampleID = np.intersect1d(spec_sampleids, gcol_exp_sampleids)

else:
    sampleID = gcol_exp_sampleids

sample_size = sampleID.size

if not sample_size:
    raise SystemExit("There is no overlapped sample IDs between gene expression file, genotype file, and sampleID file.")

# get columns to read in
g_cols_ind, g_dtype = tg.genofile_cols_dtype(g_cols, args.genofile_type, sampleID)
e_cols_ind, e_dtype = tg.exp_cols_dtype(exp_cols, sampleID)

### EXPRESSION FILE
# ### Extract expression level by chromosome
GeneExp_chunks = pd.read_csv(
    args.geneexp_path, 
    sep='\t', 
    iterator=True, 
    chunksize=10000,
    usecols=e_cols_ind,
    dtype=e_dtype)

GeneExp = pd.concat([x[x['CHROM']==args.chr] for x in GeneExp_chunks]).reset_index(drop=True)

if GeneExp.empty:
    raise SystemExit('There are no valid gene expression training data for chromosome '+args.chr+".")

GeneExp = tg.optimize_cols(GeneExp)

TargetID = GeneExp.TargetID.values
n_targets = TargetID.size

### PREP CROSS VALIDATION SAMPLES
### Split sampleIDs for cross validation
if args.cvR2:
    print("Evaluate DPR model by average R2 of 5-fold cross validation ... "+ "\nSplit sample IDs randomly into 5 folds ..."+ "\n")

    kf = KFold(n_splits=5)
    kf_splits = [(sampleID[x], sampleID[y]) for x,y in kf.split(sampleID)]
    CV_trainID, CV_testID = zip(*kf_splits)

else:
    print("Skip 5-fold cross validation ...")

### PREP OUTPUT
## print output headers to files
weight_out_cols = ['CHROM','POS', 'snpID', 'REF','ALT','TargetID','MAF','p_HWE','ES','b','beta']
pd.DataFrame(columns=weight_out_cols).to_csv(
    out_train_weight_path,
    sep='\t',
    header=True,
    index=None,
    mode='w')

info_out_cols = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName','sample_size','n_snp', 'n_effect_snp','CVR2','TrainPVALUE','TrainR2']
pd.DataFrame(columns=info_out_cols).to_csv(
    out_train_info_path,
    sep='\t',
    header=True,
    index=None,
    mode='w')

###############################################################
### thread function
def thread_process(num):
    try:
        target = TargetID[num]
        print("\nnum="+str(num)+"\nTargetID="+target)
        target_exp = GeneExp.iloc[[num]]

        start=str(max(int(target_exp.GeneStart)-args.window,0))
        end=str(int(target_exp.GeneEnd)+args.window)

        # READ IN AND PROCESS GENOTYPE DATA 
        # Requirement for input vcf file: Must be bgzip and tabix
        ### select corresponding vcf file by tabix
        print("Loading genotype data for Gene: " + target +"\n")
        g_proc_out = tg.call_tabix(args.geno_path, args.chr, start, end)

        if not g_proc_out:
            print("There is no genotype data for gene TargetID="+target+".")
            return None

        print("Preparing DPR input for Gene: " + target + "\n")
        target_geno = pd.read_csv(StringIO(g_proc_out.decode('utf-8')),
                sep='\t',
                low_memory=False,
                header=None,
                usecols=g_cols_ind,
                dtype=g_dtype)
        target_geno.columns = [g_cols[i] for i in target_geno.columns]
        target_geno = tg.optimize_cols(target_geno)

        # get snpIDs
        target_geno['snpID'] = tg.get_snpIDs(target_geno)
        target_geno = target_geno.drop_duplicates(['snpID'],keep='first').reset_index(drop=True)

        # prep vcf file
        if args.genofile_type=='vcf':
            target_geno = tg.check_prep_vcf(target_geno, args.format, sampleID)

        # reformat sample values
        target_geno[sampleID]=target_geno[sampleID].apply(lambda x: tg.reformat_sample_vals(x,args.format), axis=0)

        # get, filter maf
        target_geno = tg.calc_maf(target_geno, sampleID, args.maf)

        # get, filter p_HWE
        target_geno = tg.calc_p_hwe(target_geno, sampleID, args.hwe)

        snp_annot = target_geno[['snpID','POS','CHROM']]

        # 5-FOLD CROSS-VALIDATION
        if args.cvR2:
            print("Running 5-fold CV for Gene: "+target+"\n")
            do_cv_args = [target, target_geno, target_exp, snp_annot, CV_trainID, CV_testID]

            k_fold_R2 = [do_cv(i, *do_cv_args) for i in range(5)]

            avg_r2_cv = sum(k_fold_R2)/5

            if avg_r2_cv < 0.005:
                print('Average R2 by 5-fold CV =' + str(avg_r2_cv) + ' < 0.005 for ' + target +'\nSkip running DPR for gene '+ target +'\n')
                return None

        else:
            avg_r2_cv = 0
            print("Skip evaluation by 5-fold CV R2 ..." + "\n")


        # FINAL MODEL TRAINING
        print("Running DPR training for Gene:"+target +"\n")
        dpr_file_dir = args.out_dir + '/DPR_Files/'
        
        bimbam = target_geno[np.concatenate((['snpID','REF','ALT'],sampleID))]

        pheno = target_exp[sampleID].T

        # PREP INPUT FILES, CALL DPR, READ IN DPR OUTPUT
        try:
            dpr_out = prep_call_dpr(bimbam, pheno, snp_annot, dpr_file_dir, target)

        except subprocess.CalledProcessError as err:
            print("DPR failed for TargetID: " + target)
            return None

        # FILTER FOR SNPS WITH ES!=0
        n_snp = dpr_out.ES.size

        dpr_out = dpr_out[dpr_out.ES!=0]

        n_effect_snp = dpr_out.ES.size

        # R2 CALCULATION
        out_weights = dpr_out[['snpID','ES']]

        bimbam = bimbam.drop(columns=['REF','ALT'])

        train_pvalue, train_rsquared = calc_r2(out_weights, bimbam, pheno)

        # OUTPUT TARGET WEIGHTS TO FILE
        # initialize df with MAF, pHWE, other info
        target_weights = target_geno[['CHROM','POS','REF','ALT','snpID','p_HWE','MAF']].copy()
        target_weights = target_weights[target_weights.snpID.isin(dpr_out.snpID)]
        target_weights['TargetID'] = target

        # merge with dpr output weights, reorder columns using existing col list
        target_weights = target_weights.merge(
            dpr_out, 
            left_on='snpID',
            right_on='snpID',
            how='outer')[weight_out_cols]

        target_weights.to_csv(
            out_train_weight_path,
            sep='\t',
            header=None,
            index=None,
            mode='a')        

        # OUTPUT TARGET TRAINING INFO TO FILE
        # initialize dataframe for storing training info
        train_info = target_exp[
            ['CHROM','GeneStart','GeneEnd','TargetID','GeneName']].copy()
        train_info['sample_size'] = sample_size
        train_info['n_snp'] = n_snp
        train_info['n_effect_snp'] = n_effect_snp
        train_info['CVR2'] = avg_r2_cv
        train_info['TrainPVALUE'] = train_pvalue
        train_info['TrainR2'] = train_rsquared

        # output training info
        train_info.to_csv(
            out_train_info_path,
            sep='\t',
            header=None,
            index=None,
            mode='a')

    except Exception as e:
        e_type, e_obj, e_tracebk = sys.exc_info()
        e_line_num = e_tracebk.tb_lineno

        e, e_type, e_line_num = [str(x) for x in [e, e_type, e_line_num]]

        print('Caught a type '+ e_type +' exception for TargetID='+target+', num=' + str(num) + ' on line '+e_line_num+':\n' + e )

    finally:
        # print info to log do not wait for buffer to fill up
        sys.stdout.flush() 

##############################################################
### start thread  process

# if (args.thread < int(len(EXP)/100) | args.thread > len(EXP)):
    # args.thread = (int(len(EXP)/100)+1)*100

if __name__ == '__main__':
    print("Starting DPR training for "+str(n_targets)+" target genes.")
    pool = multiprocessing.Pool(args.thread)
    pool.map(thread_process,[num for num in range(n_targets)])
    pool.close()
    pool.join()
    print('Done.\n')


############################################################
### time calculation
elapsed_sec = time()-start_time
elapsed_time = tg.format_elapsed_time(elapsed_sec)
print("Computation time (DD:HH:MM:SS): " + elapsed_time)


