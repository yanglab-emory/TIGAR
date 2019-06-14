#!/usr/bin/env python

##############################################################################################
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

##############################################################################################
### timecalculation
start_time=time.clock()

##############################################################################################
### variable needed
parser = argparse.ArgumentParser(description='manual to this script')

### Gene annotation file
parser.add_argument("--Gene",type=str,default = None)

### GWAS Z score file
parser.add_argument('--Zscore',type=str,default = None)

### Header of GWAS Z score file
parser.add_argument('--Zscore_names',type=str,default=None)

### Weight
parser.add_argument('--Weight',type=str,default = None)

### Header of Weight file
parser.add_argument('--Weight_names',type=str,default=None)

### Reference covariance file
parser.add_argument('--Covar',type=str,default = None)

### chromosome number
parser.add_argument('--chr_num',type=int,default = None)

### window
parser.add_argument('--window',type=int,default = None)

### Number of thread
parser.add_argument('--thread',type=int,default = None)

### Output dir
parser.add_argument('--out_prefix',type=str,default = None)

args = parser.parse_args()

################################################################################################
### variable checking
print("Gene annotation file:"+args.Gene)
print("GWAS test result:"+args.Zscore)
print("Snps weight:"+args.Weight)
print("Reference covariance file:"+args.Covar)
print("Chromosome number:"+str(args.chr_num))
print("window="+str(args.window))
print("Number of thread:"+str(args.thread))
print("Output dir:"+args.out_prefix)

################################################################################################
### Read in gene annotation 
Gene = pd.read_csv(args.Gene,sep='\t')

Gene = (Gene >> mask(Gene['CHROM'].astype('str')==str(args.chr_num))).reset_index(drop=True)

TargetID = np.array(Gene.TargetID)

pd.DataFrame(columns=['CHROM','GeneStart','GeneEnd','TargetID','Zscore','Pvalue']).to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_association_study.txt',
                                                                                          sep='\t',index=None,header=True,mode='a')

Weight_names = pd.read_csv(args.Weight_names,sep='\t')
Zscore_names = pd.read_csv(args.Zscore_names,sep='\t')

def thread_process(num):
    Gene_temp = Gene >> mask(Gene.TargetID == TargetID[num])

    start=max(int(Gene_temp.GeneStart)-args.window,0)
    end=max(int(Gene_temp.GeneEnd)+args.window,0)
    
    Zscore_process = subprocess.Popen(["tabix"+" "+args.Zscore+" "+str(args.chr_num)+":"+str(start)+"-"+str(end)],
                                      shell=True,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
    Zscore_out = Zscore_process.communicate()[0]

    Weight_process = subprocess.Popen(["tabix"+" "+args.Weight+" "+str(args.chr_num)+":"+str(start)+"-"+str(end)],
                                      shell=True,
                                      stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE)
    Weight_out = Weight_process.communicate()[0]
    
    if len(Zscore_out)==0 | len(Weight_out) == 0:
        print("No GWAS test result or weight for this gene:"+TargetID[num])
    else:
        print("Running association study for gene:"+TargetID[num])
        Zscore = pd.read_csv(StringIO(Zscore_out.decode('utf-8')),sep='\t',header=None,low_memory=False)
        Zscore.columns = np.array(tuple(Zscore_names))

        Weight = pd.read_csv(StringIO(Weight_out.decode('utf-8')),sep='\t',header=None,low_memory=False)
        Weight.columns = np.array(tuple(Weight_names))

        
        Zscore['snpID'] = (Zscore['CHROM'].astype('str')+':'+Zscore['POS'].astype('str')
                           +':'+Zscore.REF+':'+Zscore.ALT)

        Weight['snpID'] = (Weight['CHROM'].astype('str')+':'+Weight['POS'].astype('str')
                           +':'+Weight.REF+':'+Weight.ALT)

        snp_overlap = np.intersect1d(np.array(Weight.snpID),np.array(Zscore.snpID))
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
        covar_process = subprocess.Popen(["tabix"+" "+args.Covar+" "+str(args.chr_num)+":"+str(min(ZW.POS))+"-"+str(max(ZW.POS))],
                                         shell=True,
                                         stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE)
        
        out=covar_process.communicate()[0]
        if len(out)==0:
            print("No reference covariance information for Gene:"+TargetID[num])
        
        else:
            MCOV = pd.read_csv(StringIO(out.decode('utf-8')),sep='\t',header=None,low_memory=False)
            MCOV.columns = ['CHROM','POS','ID','REF','ALT','COV']
            
            MCOV['snpID'] = (MCOV['CHROM'].astype('str')+':'+MCOV['POS'].astype('str')
                             +':'+MCOV.REF+':'+MCOV.ALT)
            
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

            result = pd.DataFrame()
            result['CHROM'] = np.array(args.chr_num).ravel()
            result['GeneStart'] = np.array(Gene_temp.GeneStart).ravel()
            result['GeneEnd'] = np.array(Gene_temp.GeneEnd).ravel()
            result['TargetID'] = np.array(TargetID[num]).ravel()
            result['Zscore'] = burden_Z
            ### p-value for chi-square test
            result['Pvalue'] = 1-chi2.cdf(burden_Z**2,1)

            result.to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_association_study.txt',
                          sep='\t',index=None,header=None,mode='a')

########################################################################################################
### thread process
pool = multiprocessing.Pool(args.thread)

pool.map(thread_process,[num for num in range(len(TargetID))])

pool.close()
pool.join()

#########################################################################################################
### time calculation
time=round((time.clock()-start_time)/60,2)

print(str(time)+' minutes')







