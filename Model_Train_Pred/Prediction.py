#!/usr/bin/env python

############################################################
# import packages needed
import argparse
import time
import os
import subprocess
from subprocess import *
import io
from io import StringIO
import pandas as pd
import numpy as np
from numpy import *
from dfply import *
import multiprocessing

##########################################################
### time calculation
start_time=time.clock()

###########################################################
# Construct Dataframe for Analysis

### input each sample genotype
### For GT Format:
###  code '0|0' or '0/0' as 0
###  code ('0|1' or '1|0')  or ('0/1' or '1/0') as 1
###  code '1|1' or '1/1' as 2
###  code '.|.' or './.' as nan(missing)

### For DS Format:
### code '.' as nan(missing)
def geno_reform(data,Format):
    if Format=='GT':
        data[(data=='0|0')|(data=='0/0')]=0
        data[(data=='1|0')|(data=='1/0')|(data=='0|1')|(data=='0/1')]=1
        data[(data=='1|1')|(data=='1/1')]=2
        data[(data=='.|.')|(data=='./.')]=nan
    elif Format=='DS':
        data[(data==".")]=nan
    return data

### For vcf input
### Split input dataframe by Format. ex, '0|0:0.128'
### Input:
### 1. data:The first nine columns fixed
### 2. Format: GT or DS

### Output:
### 1. The First six columns of output dataframe should be:
###    1) CHROM
###    2) POS
###    3) ID (i.e. rsID)
###    4) REF
###    5) ALT
###    6) snpID (CHROM:POS:REF:ALT)
###    7) p_HWE:p-value for Hardy Weinberg Equilibrium exact test
###    8) MAF: Minor Allele Frequency (range from 0~1)
###    9) Samples gene variance splited by Format (GT or DS)
### 2. Output dataframe will be selected by p_HWE and MAF:
###    Default threshold is p_HWE > 10**(-3) and MAF > 0.01

def CHR_Reform_vcf(data,Format,sampleID):
    data['snpID']=(data['CHROM'].astype('str')+':'+data['POS'].astype('str')+':'
                  +data.REF+':'+data.ALT)
    
    CHR=data >> select(data[['CHROM','POS','ID','REF','ALT','snpID']],data[sampleID])
    
    CHR=CHR.drop_duplicates(['snpID'],keep='first')
        
    indicate=data.FORMAT[0].split(":").index(Format)
    CHR[sampleID]=CHR[sampleID].applymap(lambda x:x.split(":")[indicate])
    
    CHR[sampleID]=CHR[sampleID].apply(lambda x:geno_reform(x,Format),axis=0)
        
    ### Calculating folded MAF and p_HWE by SNPs
    temp=pd.DataFrame((CHR >> select(CHR[sampleID])),dtype=np.float)
    
    ### Calculate MAF(range from 0-1)
    CHR['MAF']=temp.apply(lambda x:sum(x)/(2*len(x.dropna())),axis=1)

    
    ### Dealing with NaN
    CHR[np.hstack(([sampleID,'MAF']))]=CHR[np.hstack(([sampleID,'MAF']))].apply(lambda x:x.fillna(2*x.MAF),axis=1)

    CHR=CHR >> mask(CHR.MAF>=0)
    
    return CHR.rename(columns={'MAF':'MAF_test'})

### For dosages input
### Input:
### 1. data:The first five columns fixed
### 2. Format: GT or DS
def CHR_Reform_DS(data,sampleID):
    data['snpID']=(data['CHROM'].astype('str')+':'+data['POS'].astype('str')+':'
                  +data.REF+':'+data.ALT)

    data=data.drop_duplicates(['snpID'],keep='first')
    
    data[data[sampleID].astype('str')=='.']=nan
    
    data[sampleID]=data[sampleID].astype('float')

    data['MAF']=data[sampleID].apply(lambda x:sum(x)/(2*len(x.dropna())),axis=1)
    
    data[np.hstack(([sampleID,'MAF']))]=data[np.hstack(([sampleID,'MAF']))].apply(lambda x:x.fillna(2*x.MAF),axis=1)

    data=data >> mask(data.MAF>=0)
    
    return data.rename(columns={'MAF':'MAF_test'})

#######################################################################
### Input Arguments for GReX Prediction
#######################################################################

# --chr: Chromosome number need to be specified with respect to the genotype input data
# --weight: Path for SNP weight (eQTL effect size) file 
# --test_sampleID: Path for a file with sampleIDs that should be contained in the genotype file
# --genofile: Path for the training genotype file (bgzipped and tabixed) 
# --genofile_tye: Genotype file type: "vcf" or "dosage"
# --test_geno_colnames : File with column heads of genotype file
# --Format: Genotype format in VCF file that should be used: "GT" (default) for genotype data or "DS" for dosage data, only required if the input genotype file is of VCF file
# --window: Window size around gene transcription starting sites (TSS) for selecting cis-SNPs for fitting gene expression prediction model (default 1000000 for +- 1MB region around TSS)
# --maf_diff: MAF difference threshold for matching SNPs from eQTL weight file and test genotype file. If SNP MAF difference is greater than maf_diff (default 0.2), , the SNP will be excluded
# --thread: Number of threads for parallel computation (default 1)
# --out_dir: Output directory (will be created if not exist)

#######################################
### variables need
parser = argparse.ArgumentParser(description='manual to this script')

### eQTL weight file
parser.add_argument('--weight',type=str,default = None)

### Test sampleID
parser.add_argument('--test_sampleID',type=str,default = None)

### Specified chromosome number
parser.add_argument('--chr',type=int,default = None)

### Test genotype files
parser.add_argument('--genofile',type=str,default = None)
parser.add_argument('--test_geno_colnames', type=str, default = None)

### Specified input file type(vcf or dosages)
parser.add_argument('--geno',type=str,default = None)

### 'DS' or 'GT' for VCF genotype file
parser.add_argument('--Format',type=str,default=None)

### window
parser.add_argument('--window',type=int,default=None)

### Gene annotation file
parser.add_argument('--gene_anno',type=str,default = None)

### number of thread
parser.add_argument('--thread',type=int,default = None)

### Threshold of difference of maf between training data and testing data
parser.add_argument('--maf_diff',type=float,default=None)

### output dir
parser.add_argument('--out_dir',type=str,default=None)

args = parser.parse_args()

### check input command
print("********************************\n   Imput Arguments\n********************************\n")
print("Chrmosome : "+str(args.chr)+ "\n")
print("eQTL weight file : "+args.weight+ "\n")
print("Test gene annotation file : "+args.gene_anno+ "\n")
print("Test sampleID file : "+args.test_sampleID+ "\n")
print("Test genotype file : "+args.genofile+ "\n")
# print("Column names of test genotype file:"+args.geno_colnames+ "\n")

if args.geno=='vcf':
    print("VCF genotype file is used for prediction with genotype format : " + args.Format + "\n")
elif args.geno=='dosage':
    print("Dosage genotype file is used for prediction."+ "\n")
else:
    raise SystemExit("Please specify input test genotype file as either 'vcf' or 'dosage'."+ "\n")

print("Gene region size : window = "+str(args.window)+ "\n")

print("MAF difference threshold for matching SNPs from eQTL weight file and test genotype file : "+str(args.maf_diff)+ "\n")

print("Number of threads : "+str(args.thread)+ "\n")
print("Output dir : "+args.out_dir+ "\n")

############# Load eQTL weights (ES) #################################
Result=pd.read_csv(args.weight,sep='\t')
Result['CHROM']=Result['CHROM'].astype('int')
Result['POS']=Result['POS'].astype('int')
Result['snpID']=(Result['CHROM'].astype('str')+':'+Result['POS'].astype('str')
                 +':'+Result.REF+':'+Result.ALT)
if len(Result.ES)==0:
    raise SystemExit('There is no valid eQTL weights.')

#############  Training Information
Gene_Info=pd.read_csv(args.gene_anno,sep='\t')

# Load genotype column names of test genotype file
test_names=pd.read_csv(args.test_geno_colnames,sep='\t').rename(columns={'#CHROM':'CHROM'})
test_names=np.array(tuple(test_names))

# Load test sample IDs
sampleID=pd.read_csv(args.test_sampleID,sep='\t',header=None)
sampleID=np.array(sampleID).ravel()

### Define TargetID
TargetID=unique(np.array(Result.TargetID))

out_header=pd.DataFrame(columns=np.hstack((['CHROM','GeneStart','GeneEnd','TargetID','GeneName'],
                                           sampleID)))

out_header.to_csv(args.out_dir+'/CHR'+str(args.chr) + '_Pred_GReX.txt',
                  sep='\t', index=None, header=True, mode='w')

### thread function
def thread_process(num):
    Info_temp=Gene_Info >> mask(Gene_Info.TargetID==TargetID[num])
    Info_temp=Info_temp.reset_index(drop=True)
    
    start = max(int(Info_temp.GeneStart)-args.window,0)
    end = max(int(Info_temp.GeneEnd)+args.window,0)
    
    test_process=subprocess.Popen(["tabix"+" "+args.genofile+" "+ str(args.chr)+":"+str(start)+"-"+str(end)],
                                  shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    
    out=test_process.communicate()[0]
   # print(str(len(out)))
    
    if len(out)==0:
        print("There is no test genotype data for Gene:" + TargetID[num] + "\n")
    else:
        print("Predict GReX for Gene : "+TargetID[num])
        Test_temp=pd.read_csv(StringIO(out.decode('utf-8')),sep='\t',header=None,low_memory=False)
        Test_temp.columns = np.array(tuple(test_names))
        Test_temp=Test_temp.reset_index(drop=True)
    
        if args.geno=='vcf':
            if args.Format not in unique(Test_temp.FORMAT)[0].split(":"):
                raise SystemExit('Given genotype Format (e.g., GT or ES) need to be specified in the FORMAT column of test VCF genotype file.')
            else:
                Chrom = CHR_Reform_vcf(Test_temp,args.Format,sampleID)
        elif args.geno=='dosages':
            Chrom = CHR_Reform_DS(Test_temp,sampleID)
        
        ### Intersect SNPs from eQTL weight file and test genotype file
        beta_temp = Result >> mask(Result.TargetID==TargetID[num]) >> select(Result.snpID,Result.ES,Result.MAF)
        beta_temp = beta_temp.drop_duplicates(['snpID'],keep='first')
        
        overlapID = np.intersect1d(np.array(beta_temp.snpID),np.array(Test_temp.snpID))
        print("Number of SNPs overlapped between eQTL weight file and test genotype file : "+str(len(overlapID)) + "\n")
        
        if len(overlapID)==0:
            Chrom=Chrom >> drop(Chrom.snpID)
            Chrom['snpID']=(Chrom['CHROM'].astype('str')+':'+Chrom['POS'].astype('str')
                            +':'+Chrom.ALT+':'+Chrom.REF)

            overlapID_flip= np.intersect1d(np.array(beta_temp.snpID),np.array(Chrom.snpID))
            print("overlapID_flip:"+str(len(overlapID_flip)))
            
            if len(overlapID_flip)==0:
                print("No corresponding training parameter for Gene:"+TargetID[num])
            else:
                overlapID=overlapID_flip
                beta_temp['ES']=-beta_temp['ES'].astype('float')
                Chrom['MAF_test']=1-Chrom['MAF_test'].astype('float')

        if len(overlapID)!=0:
            Chrom = Chrom >> select(Chrom.snpID,Chrom.MAF_test,Chrom[sampleID])
            ### Store prediction result
            pred=pd.DataFrame()
            pred['sampleID']=np.array(sampleID).ravel()
        
            Pred = (Chrom
                    >> mask(Chrom.snpID.isin(overlapID))).merge((beta_temp
                                                                 >> mask(beta_temp.snpID.isin(overlapID))),
                                                                left_on='snpID',
                                                                right_on='snpID',
                                                                how='outer')
            
            Pred['diff'] = abs(Pred['MAF'].astype('float')-Pred['MAF_test'].astype('float'))
            Pred = Pred >> mask(Pred['diff']<=args.maf_diff) >> drop(Pred[['MAF','MAF_test','diff']])
            print("Number of SNPs used for prediction after filtered by maf_diff : "+str(len(Pred.snpID)))
            
            if len(Pred.snpID)==0:
                print("SNP MAF between training and testing samples differs greater than "+str(args.maf_diff) + "\n")
            else:
                testX_temp=Pred.T
                testX_temp.columns=testX_temp.loc['snpID']
                testX=pd.DataFrame(testX_temp.drop(['snpID','ES']),dtype='float')
                pred[TargetID[num]]=mat(testX)*(mat(Pred['ES'].astype('float')).T)
                
                pred=pred.T.reset_index(drop=True)
                pred.columns=pred.loc[0]
                pred=pred.drop([0])
                pred['TargetID']=TargetID[num]
                
                out=(Info_temp[['CHROM','GeneStart','GeneEnd',
                                'TargetID','GeneName']].merge(pred,left_on='TargetID',
                                                              right_on='TargetID',how='outer'))
                out.to_csv(args.out_dir+'/CHR'+str(args.chr)+'_Pred_GReX.txt',
                           sep='\t',index=None,header=None,mode='a')
        else:
            print("There is no matched SNPs from eQTL weight file and test genotype file. \n")

########################################################################################################
# thread begin
if (args.thread < int(len(TargetID)/100) | args.thread > len(TargetID)):
    args.thread = (int(len(EXP)/100)+1)*100

pool = multiprocessing.Pool(args.thread)

pool.map(thread_process,[num for num in range(len(TargetID))])

pool.close()

pool.join()

########################################################################################################
### time calculation
time=round((time.clock()-start_time)/60,2)

# print(str(time)+' minutes')













