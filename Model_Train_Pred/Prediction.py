#!/usr/bin/env python

##########################################################################################
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

#######################################################################################
### time calculation
start_time=time.clock()

#################################################################################################
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

#########################################################################################
### variables need
parser = argparse.ArgumentParser(description='manual to this script')

### specified training method
parser.add_argument('--model',type=str,default = None)

### for Training result
parser.add_argument('--train_result_path',type=str,default = None)
### for DPR Training information
parser.add_argument('--train_info_path',type=str,default = None)

### specified chromosome number
parser.add_argument('--chr_num',type=int,default = None)

### training vcf fir
parser.add_argument('--test_dir',type=str,default = None)
parser.add_argument('--test_names',type=str,default = None)

### SampleIDs
parser.add_argument('--test_sample',type=str,default=None)

### number of thread
parser.add_argument('--thread',type=int,default = None)

### specified input file type(vcf or dosages)
parser.add_argument('--geno',type=str,default = None)

### 'DS' or 'GT'
parser.add_argument('--Format',type=str,default=None)

### window
parser.add_argument('--window',type=int,default=None)

### Threshold of difference of maf between training data and testing data
parser.add_argument('--maf_diff',type=float,default=None)

### output dir
parser.add_argument('--out_prefix',type=str,default=None)

args = parser.parse_args()

### check input command
print(args.train_result_path)
print(args.train_info_path)

if args.model!="elastic_net" and args.model!="DPR":
	raise SystemExit("Model "+args.model+" not found.")
else:
	print("Training Model:"+args.model)

print("chromosome number:"+str(args.chr_num))
print("Testing data:"+args.test_dir)
print(args.test_names)
print("Test sampleID:"+args.test_sample)
print("Number of thread:"+str(args.thread))

if args.geno=='vcf':
    print("Using "+args.Format+" Format for Prediction.")
elif args.geno=='dosages':
    print("Using DS Format for Prediction.")
else:
	raise SystemExit("Geno file can not identify.")

print("Window:"+str(args.window))
print("MAF threshold for dropping testing snps:"+str(args.maf_diff))
print("Output dir:"+args.out_prefix)

############################################################################################
### Training Parameter
Result=pd.read_csv(args.train_result_path,sep='\t')
Result['CHROM']=Result['CHROM'].astype('int')
Result['POS']=Result['POS'].astype('int')

Result['snpID']=(Result['CHROM'].astype('str')+':'+Result['POS'].astype('str')
                 +':'+Result.REF+':'+Result.ALT)

if len(Result.ES)==0:
    raise SystemExit('No training parameters, please check input files.')

### Training Information
Train_Info=pd.read_csv(args.train_info_path,sep='\t')

test_names=pd.read_csv(args.test_names,sep='\t').rename(columns={'#CHROM':'CHROM'})
test_names=np.array(tuple(test_names))

sampleID=pd.read_csv(args.test_sample,sep='\t',header=None)
sampleID=np.array(sampleID).ravel()

### Define TargetID
TargetID=unique(np.array(Result.TargetID))

out_header=pd.DataFrame(columns=np.hstack((['CHROM','GeneStart','GeneEnd','GeneName','TargetID'],
                                           sampleID)))

out_header.to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_'+args.model+'_prediction.txt',
                  sep='\t',index=None,header=True,mode='a')


### thread function
def thread_process(num):
    Info_temp=Train_Info >> mask(Train_Info.TargetID==TargetID[num])
    Info_temp=Info_temp.reset_index(drop=True)
    
    start = max(int(Info_temp.GeneStart)-args.window,0)
    end = max(int(Info_temp.GeneEnd)+args.window,0)
    
    test_process=subprocess.Popen(["tabix"+" "+args.test_dir+" "+str(args.chr_num)+":"+str(start)+"-"+str(end)],
                                  shell=True,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    
    out=test_process.communicate()[0]
    print(str(len(out)))
    
    if len(out)==0:
        print("No testing data for Gene:"+TargetID[num])
    else:
        print("Running prediction for Gene:"+TargetID[num])
        Test_temp=pd.read_csv(StringIO(out.decode('utf-8')),sep='\t',header=None,low_memory=False)
        Test_temp.columns = np.array(tuple(test_names))
        Test_temp=Test_temp.reset_index(drop=True)
    
        if args.geno=='vcf':
            if args.Format not in unique(Test_temp.FORMAT)[0].split(":"):
                print("Format needed for training is not provided in input vcf file.")
                raise SystemExit('Geno file can not identify.')
            else:
                Chrom = CHR_Reform_vcf(Test_temp,args.Format,sampleID)
        elif args.geno=='dosages':
            Chrom = CHR_Reform_DS(Test_temp,sampleID)
        
        ### Deal with training result
        beta_temp = Result >> mask(Result.TargetID==TargetID[num]) >> select(Result.snpID,Result.ES,Result.MAF)
        beta_temp = beta_temp.drop_duplicates(['snpID'],keep='first')
        
        overlapID = np.intersect1d(np.array(beta_temp.snpID),np.array(Test_temp.snpID))
        print("overlapID:"+str(len(overlapID)))
        
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
            print("Overall ID used:"+str(len(Pred.snpID)))
            
            if len(Pred.snpID)==0:
                print("Differences of maf between training and testing are greater than "+str(args.maf_diff))
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
                                'GeneName','TargetID']].merge(pred,left_on='TargetID',
                                                              right_on='TargetID',how='outer'))
                out.to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_'+args.model+'_prediction.txt',
                           sep='\t',index=None,header=None,mode='a')
        else:
            print("Prediction Stop")

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

print(str(time)+' minutes')













