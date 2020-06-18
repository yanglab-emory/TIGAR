#!/usr/bin/env python

###################################################################
# Import packages needed
import argparse
import time
import subprocess
from subprocess import *
import pandas as pd
import numpy as np
from numpy import *
from dfply import *
import io
from io import StringIO
from io import *
from scipy.stats import chi2

import multiprocessing

#########################################################################
### timecalculation
start_time=time.clock()

############################################################
### variable needed
parser = argparse.ArgumentParser(description='Help: ')

### Gene annotation file
parser.add_argument("--gene_anno",type=str,default = None)

### GWAS Z score file
parser.add_argument('--Zscore',type=str,default = None)

### Header of GWAS Z score file
parser.add_argument('--Zscore_colnames',type=str,default=None)

### Weight
parser.add_argument('--weight',type=str,default = None)

### Header of Weight file
parser.add_argument('--weight_colnames',type=str,default=None)

### Reference covariance file
parser.add_argument('--LD',type=str,default = None)

### chromosome number
parser.add_argument('--chr',type=int,default = None)

### window
parser.add_argument('--window',type=int,default = None)

### Number of thread
parser.add_argument('--thread',type=int,default = None)

### Output dir
parser.add_argument('--out_dir',type=str,default = None)

args = parser.parse_args()

################################################################################################
### variable checking
print("Gene annotation file to specify the list of genes for TWAS : "+args.gene_anno + "\n")
print("GWAS summary statistics Z-score file : " + args.Zscore+ "\n")
print("cis-eQTL weight file : " + args.weight + "\n")
print("Reference LD genotype covariance file:"+args.LD + "\n")
print("Chromosome number : "+str(args.chr)+ "\n")
print("Test gene region including SNPs within +- window = "+str(args.window) + " base pair of GeneStart/GeneEnd positions \n")
print("Number of threads : "+str(args.thread) + "\n")
print("Output directory : " + args.out_dir + "\n")

##################################################
### Read in gene annotation 
Gene = pd.read_csv(args.gene_anno,sep='\t')

Gene = (Gene >> mask(Gene['CHROM'].astype('str')==str(args.chr))).reset_index(drop=True)

TargetID = np.array(Gene.TargetID)

pd.DataFrame(columns=['CHROM','GeneStart','GeneEnd','TargetID','GeneName','Zscore','PVALUE']).to_csv(args.out_dir+'/CHR'+str(args.chr)+'_sumstat_assoc.txt',
                     sep='\t',index=None,header=True,mode='w')

Weight_names = pd.read_csv(args.weight_colnames,sep='\t')
Zscore_names = pd.read_csv(args.Zscore_colnames,sep='\t')

def thread_process(num):
    Gene_temp = Gene >> mask(Gene.TargetID == TargetID[num])

    start=max(int(Gene_temp.GeneStart)-args.window,0)
    end=max(int(Gene_temp.GeneEnd)+args.window,0)
    
    Zscore_process = subprocess.Popen(["tabix"+" "+args.Zscore+" "+str(args.chr)+":"+str(start)+"-"+str(end)],
                                      shell=True,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
    Zscore_out = Zscore_process.communicate()[0]

    Weight_process = subprocess.Popen(["tabix"+" "+args.weight+" "+str(args.chr)+":"+str(start)+"-"+str(end)],
                                      shell=True,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
    Weight_out = Weight_process.communicate()[0]

    if (len(Zscore_out) == 0) | (len(Weight_out) == 0):
        print("There is no test SNP with non-zero cis-eQTL weights for gene:"+TargetID[num])
        return None

    print("TWAS for gene:"+TargetID[num])
    Zscore = pd.read_csv(StringIO(Zscore_out.decode('utf-8')),sep='\t',header=None,low_memory=False)
    Zscore.columns = np.array(tuple(Zscore_names))

    Weight = pd.read_csv(StringIO(Weight_out.decode('utf-8')),sep='\t',header=None,low_memory=False)
    Weight.columns = np.array(tuple(Weight_names))

    
    Zscore['snpID'] = (Zscore['CHROM'].astype('str')+':'+Zscore['POS'].astype('str')
                       +':'+Zscore.REF+':'+Zscore.ALT)

    Weight['snpID'] = (Weight['CHROM'].astype('str')+':'+Weight['POS'].astype('str')
                       +':'+Weight.REF+':'+Weight.ALT)

    snp_overlap = np.intersect1d(np.array(Weight.snpID),np.array(Zscore.snpID))

    if snp_overlap.size == 0:
        print("There are no test SNPs with both non-zero cis-eQTL weights and GWAS summary statistics for gene:"+TargetID[num])
        return None

    # merge Zscore and Weight File
    ZW = (Weight >> mask(Weight.snpID.isin(snp_overlap))).merge((Zscore
                                                                 >>mask(Zscore.snpID.isin(snp_overlap))
                                                                 >>select(Zscore.snpID,Zscore.Zscore)),
                                                                left_on='snpID',
                                                                right_on='snpID')
    ZW = ZW.drop_duplicates(['snpID'],keep='first')
    ZW = ZW.sort_values(by='POS')
    ZW = ZW.reset_index(drop=True)
    
    ### Read in reference covariance matrix file
    covar_process = subprocess.Popen(["tabix"+" "+args.LD+" "+str(args.chr)+":"+str(min(ZW.POS))+"-"+str(max(ZW.POS))],
                                     shell=True,
                                     stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
    
    out=covar_process.communicate()[0]
    if len(out)==0:
        print("No reference covariance information for Gene:"+TargetID[num])
        return None
    
    MCOV = pd.read_csv(StringIO(out.decode('utf-8')),sep='\t',header=None,low_memory=False)
    MCOV.columns = ['CHROM','POS','ID','REF','ALT','COV']
    
    MCOV['snpID'] = (MCOV['CHROM'].astype('str')+':'+MCOV['POS'].astype('str')
                     +':'+MCOV.REF+':'+MCOV.ALT)
    
    snp_overlap_cov = np.intersect1d(snp_overlap,np.array(MCOV.snpID))
    
    if snp_overlap_cov.size == 0:
         print("No reference covariance information for target SNPs for gene:"+TargetID[num])
         return None

    MCOV = MCOV.drop_duplicates(['snpID'],keep='first')

    MCOV['COV'] = MCOV['COV'].apply(lambda x:np.array(x.split(',')).astype('float'))
    MCOV = MCOV.sort_values(by='POS')
    MCOV = MCOV.reset_index(drop=True)
    
    indicates = (MCOV >> mask(MCOV.snpID.isin(ZW.snpID))).index
    
    V_temp = np.zeros((len(indicates),len(indicates)))
    
    for i in range(len(indicates)):
        cov_temp = MCOV.COV.loc[indicates[i]]
        N = len(cov_temp)
        
        for j in range(i,len(indicates)):
            if indicates[j] - indicates[i] < N:
                V_temp[i,j] = cov_temp[indicates[j]-indicates[i]]
            else:
                V_temp[i,j] = 0
                
    V=V_temp+V_temp.T-np.diag(V_temp.diagonal())
    
    ZW = ZW >> mask(ZW.snpID.isin(MCOV.snpID))
    ### Calculate burden Z score
    burden_Z = mat(ZW.Zscore)*mat(ZW.ES).T/sqrt(mat(ZW.ES)*V*mat(ZW.ES).T)
    
    if np.isnan(burden_Z):
    	print("Could not calculate burden_Z: NaN value.")
    	return None
    
    result = pd.DataFrame()
    result['CHROM'] = np.array(args.chr).ravel()
    result['GeneStart'] = np.array(Gene_temp.GeneStart).ravel()
    result['GeneEnd'] = np.array(Gene_temp.GeneEnd).ravel()
    result['TargetID'] = np.array(TargetID[num]).ravel()
    result['GeneName'] = np.array(Gene_temp.GeneName).ravel()
    result['TWAS_Zscore'] = burden_Z
    ### p-value for chi-square test
    result['PVALUE'] = 1-chi2.cdf(burden_Z**2,1)

    result.to_csv(args.out_dir+'/CHR'+str(args.chr)+'_sumstat_assoc.txt',
                  sep='\t',index=None,header=None,mode='a')


###############################################################
### thread process
pool = multiprocessing.Pool(args.thread)

pool.map(thread_process,[num for num in range(len(TargetID))])

pool.close()
pool.join()

############################################################
### time calculation
time=round((time.clock()-start_time)/60,2)

# print(str(time)+' minutes')







