#/usr/bin/env python

###############################################################################################
# Import packages needed
import argparse
import time
import subprocess
from subprocess import *
import shlex
import io
from io import StringIO
import pandas as pd
import numpy as np
from numpy import *
from dfply import *

from sklearn.model_selection import KFold
import statsmodels.api as sm
import scipy.stats as stats

import multiprocessing

#################################################################################################
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

### Calculating p-value for Hardy Weinberg Equilibrium exact test

### gij denote value for gene variance for jth samples in ith SNPs
### 0 <= gij< 0.5 denote as 0
### 0.5 <= gij < 1.5 denote as 1
### 1.5 <= gij <2 denote as 2

### Input value:
### 1.obs_hets: Observed heterozygosity = Number of 1 in each SNPs(i.e. 0.5 <= gij < 1.5)
### 2.obs_hom1: Observed AA homozygosity = Number of 0 in each SNPs(i.e. 0 <= gij< 0.5)
### 3.obs_hom2: Observed aa homozygosity = Number of 2 in each SNPs(i.e. 1.5 <= gij <= 2)

### Output: p-value for Hardy Weinberg Equilibrium exact test

def HWE(obs_hets,obs_hom1,obs_hom2):
    if obs_hom1 < 0 or obs_hom2 < 0 or obs_hets < 0:
        raise Exception("FATAL ERROR - SNP-HWE: Current genotype configuration (%s  %s %s) includes negative count" % (obs_hets, obs_hom1, obs_hom2))

    obs_homc = obs_hom2 if obs_hom1 < obs_hom2 else obs_hom1
    obs_homr = obs_hom1 if obs_hom1 < obs_hom2 else obs_hom2

    rare_copies = 2*obs_homr + obs_hets
    genotypes   = obs_hets + obs_homc + obs_homr

    het_probs = [0.0]*(rare_copies + 1)

    #start at midpoint
    mid = int(rare_copies*(2*genotypes - rare_copies)/(2*genotypes))

    #check to ensure that midpoint and rare alleles have same parity
    if (rare_copies & 1)^(mid & 1):
        mid += 1

    curr_hets = mid
    curr_homr = (rare_copies - mid) / 2
    curr_homc = genotypes - curr_hets - curr_homr

    het_probs[mid] = 1.0
    sum = float(het_probs[mid])

    for curr_hets in range(mid,1,-2):
        het_probs[curr_hets - 2] = het_probs[curr_hets]*curr_hets*(curr_hets - 1.0)/(4.0*(curr_homr + 1.0)*(curr_homc + 1.0))

        sum += het_probs[curr_hets - 2];

        # 2 fewer heterozygotes for next iteration -> add one rare, one common homozygote
        curr_homr += 1
        curr_homc += 1

    curr_hets = mid
    curr_homr = (rare_copies - mid)/2
    curr_homc = genotypes - curr_hets - curr_homr

    for curr_hets in range(mid,rare_copies-1,2):
        het_probs[curr_hets + 2] = het_probs[curr_hets]*4.0*curr_homr*curr_homc/((curr_hets + 2.0)*(curr_hets + 1.0))

        sum += het_probs[curr_hets + 2]

        #add 2 heterozygotes for next iteration -> subtract one rare, one common homozygote
        curr_homr -= 1
        curr_homc -= 1

    for i in range(0,rare_copies + 1):
        het_probs[i] /= sum

    #alternate p-value calculation for p_hi/p_lo
    p_hi = float(het_probs[obs_hets])
    for i in range(obs_hets,rare_copies+1):
        p_hi += het_probs[i]

    p_lo = float(het_probs[obs_hets])
    for i in range(obs_hets-1,-1,-1):
        p_lo += het_probs[i]

    p_hi_lo = 2.0 * p_hi if p_hi < p_lo else 2.0 * p_lo

    p_hwe = 0.0
    #  p-value calculation for p_hwe
    for i in range(0,rare_copies + 1):
        if het_probs[i] > het_probs[obs_hets]:
            continue
        p_hwe += het_probs[i]

    p_hwe = 1.0 if p_hwe > 1.0 else p_hwe

    return p_hwe


### Prepare for HWE input
def p_HWE(data):
    if len(data)==0:
        p_hwe=nan
    else:
        N_hets=len(data[(data>=0.5)&(data<1.5)])
        N_aa=len(data[(data>=0)&(data<0.5)])
        N_AA=len(data[(data>=1.5)&(data<=2)])

        p_hwe=HWE(N_hets,N_AA,N_aa)

    return p_hwe


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

def CHR_Reform_vcf(data,Format,p_hwe,maf):
    names = data.columns
    sampleID = names[9:]

    data['snpID']=(data['CHROM'].astype('str')+":"+data['POS'].astype('str')
                   +":"+data.REF+":"+data.ALT)

    CHR=data >> select(data[['CHROM','POS','ID','REF','ALT','snpID']],data[sampleID])
    CHR=CHR.drop_duplicates(['snpID'],keep='first')

    indicate=data.FORMAT[0].split(":").index(Format)
    CHR[sampleID]=CHR[sampleID].applymap(lambda x:x.split(":")[indicate])

    CHR[sampleID]=CHR[sampleID].apply(lambda x:geno_reform(x,Format),axis=0)

    ### Calculating folded MAF and p_HWE by SNPs
    temp=pd.DataFrame((CHR >> select(CHR[sampleID])),dtype=np.float)

    ### Calculating p_HWE
    CHR['p_HWE']=temp.apply(lambda x:p_HWE(x.dropna()),axis=1)

    ### Calculate MAF(range from 0-1)
    CHR['MAF']=temp.apply(lambda x:sum(x)/(2*len(x.dropna())),axis=1)
    ### Dealing with NaN
    CHR[np.hstack(([sampleID,'MAF']))] = CHR[np.hstack(([sampleID,'MAF']))].apply(lambda x:x.fillna(2*x.MAF),axis=1)

    return CHR >> mask((CHR.p_HWE > p_hwe)&(CHR.MAF > maf))

### For dosages input
### Input:
### 1. data:The first five columns fixed
### 2. Format: DS
def CHR_Reform_DS(data,p_hwe,maf):
    names=data.columns
    sampleID=names[5:]

    data['snpID']=(data['CHROM'].astype('str')+':'+data['POS'].astype('str')
                   +':'+data.REF+':'+data.ALT)
    data=data.drop_duplicates(['snpID'],keep='first')

    data[data[sampleID].astype('str')=='.']=nan

    data[sampleID]=data[sampleID].astype('float')

    data['p_HWE']=data[sampleID].apply(lambda x:p_HWE(x.dropna()),axis=1)

    data['MAF']=data[sampleID].apply(lambda x:sum(x)/(2*len(x.dropna())),axis=1)

    data[np.hstack(([sampleID,'MAF']))] = data[np.hstack(([sampleID,'MAF']))].apply(lambda x:x.fillna(2*x.MAF),axis=1)

    return data >> mask((data.p_HWE > p_hwe)&(data.MAF > maf))

#######################################################################################################################

### Merge SNPs information with Gene Expression Level by sampleID

def Info_Gene(Chr,Exp):
    sampleID_temp = np.intersect1d(np.array(Chr.columns),np.array(Exp.columns))

    sampleID = np.delete(sampleID_temp,np.where(sampleID_temp=='CHROM'))

    CHR_target = Chr >> select(Chr[['CHROM','POS','REF','ALT','snpID']],
                               Chr[sampleID])

    Exp_target = Exp >> select(Exp.TargetID,Exp[sampleID])

    info = pd.concat([CHR_target,Exp_target],sort=False)

    return info,sampleID

#######################################################################################################################
### variables need
parser = argparse.ArgumentParser(description='manual to this script')

### for Gene annotation and Expression level file
parser.add_argument('--Gene_Exp_path',type=str,default = None)

### for training sampleID
parser.add_argument('--train_sample',type=str,default = None)

### specified chromosome number
parser.add_argument('--chr_num',type=int,default = None)

### Number of Thread
parser.add_argument('--thread',type=int,default = None)

### training vcf file
parser.add_argument('--train_dir',type=str,default = None)
parser.add_argument('--train_names',type=str,default = None)

### specified input file type(vcf or dosages)
parser.add_argument('--geno',type=str,default = None)

### 'DS' or 'GT'
parser.add_argument('--Format',type=str,default=None)

### p-value for HW test
parser.add_argument('--hwe',type=float,default=None)
### maf
parser.add_argument('--maf',type=float,default=None)

### window
parser.add_argument('--window',type=int,default=None)

### model to run DPR
parser.add_argument('--dpr',type=int,default=None)

### define effect-size
parser.add_argument('--ES', type=str,default=None)

### output dir
parser.add_argument('--out_prefix',type=str,default=None)

args = parser.parse_args()

### check input command
print("Gene annotation and Expressiop level file:"+args.Gene_Exp_path)
print("Training sampleID path:"+args.train_sample)
print("Chromosome number:"+str(args.chr_num))
print("Number of thread:"+str(args.thread))
print("Training file:"+args.train_dir)
print(args.train_names)
if args.geno=='vcf':
    print("Using "+args.Format+" Format for Training")
elif args.geno=='dosages':
    print("Using DS Format for Training.")
print("window="+str(args.window))
print("Threshold for MAF:"+str(args.maf))
print("Threshold for p-value of HW test:"+str(args.hwe))
print("Runing DPR Model with dpr="+str(args.dpr))
print("Using effect-size:"+args.ES)
print("Output dir:"+args.out_prefix)

#####################################################################################################################
# Prepare DPR input

### Read in Gene annotation and Expression level file (text file)
### First five columns should be fixed:
### 1.Chrmosome Number
### 2.GeneStart Posistion
### 3.GeneEnd Position
### 4.TargetID (i.e.GeneID, treated as unique annotation for each gene)
### 5.Gene Name
Gene_Exp = pd.read_csv(args.Gene_Exp_path,sep='\t',low_memory=False)
Gene_header = np.array(Gene_Exp.columns)
Gene_header[0:5] = ['CHROM','GeneStart','GeneEnd','TargetID','GeneName']
Gene_Exp.columns = Gene_header

train_sampleID = pd.read_csv(args.train_sample,sep='\t',header=None)
train_sampleID = np.array(train_sampleID).ravel()

### seperate sampleIDs for cross validation
CV_trainID = []
CV_testID = []

kf=KFold(n_splits=5)
for train_index,test_index in kf.split(train_sampleID):
    CV_trainID.append(np.array(','.join(train_sampleID[train_index])))
    CV_testID.append(np.array(','.join(train_sampleID[test_index])))

CV_trainID = pd.DataFrame(CV_trainID)
CV_testID = pd.DataFrame(CV_testID)

CV_trainID = CV_trainID.apply(lambda x:x.str.split(","))
CV_testID = CV_testID.apply(lambda x:x.str.split(","))

### Extract expression level by chromosome
EXP = Gene_Exp >> mask(Gene_Exp['CHROM'].astype('str')==str(args.chr_num))
EXP = EXP >> select(EXP[Gene_Exp.columns[0:5]],EXP[train_sampleID])
EXP = EXP.reset_index(drop=True)
if len(EXP.CHROM)==0:
    raise SystemExit("No training data.")

train_names=pd.read_csv(args.train_names,sep='\t').rename(columns={'#CHROM':'CHROM'})

### Initialize output dataframe
param_out=pd.DataFrame(columns=['CHROM','POS','REF','ALT','TargetID','n_miss',
                                'b','beta','ES','gamma','p_HWE','MAF'])

param_out.to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_DPR_training_param.txt',
                 sep='\t',header=True,index=None,mode='a')

Info_out=pd.DataFrame(columns=['CHROM','GeneStart','GeneEnd','TargetID','GeneName',
                               'snp_size','effect_snp_size','sample_size','5-fold-CV-R2','TrainPVALUE','Train-R2'])

Info_out.to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_DPR_training_info.txt',
                sep='\t',header=True,index=None,mode='a')

def thread_process(num):
    Exp_temp = pd.DataFrame(EXP.loc[num]).T
    TargetID = np.array(Exp_temp.TargetID)[0]

    start=max(int(Exp_temp.GeneStart)-args.window,0)
    end=max(int(Exp_temp.GeneEnd)+args.window,0)

    # Requirement for input vcf file:Must be bgzip and tabix
    ### select corresponding vcf file by tabix
    train_process=subprocess.Popen(["tabix"+" "+args.train_dir+" "+str(args.chr_num)+":"+str(start)+"-"+str(end)],
                                   shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
    out=train_process.communicate()[0]

    if len(out)==0:
        print("No corresponding vcf data for this Gene:"+TargetID)

    else:
        print("Preparing DPR input for Gene:"+TargetID)

        Chr_temp=pd.read_csv(StringIO(out.decode('utf-8')),sep='\t',header=None,low_memory=False)
        Chr_temp.columns=np.array(tuple(train_names))
        Chr_temp=Chr_temp.reset_index(drop=True)

        if args.geno=='vcf':
            if args.Format not in unique(Chr_temp.FORMAT)[0].split(":"):
                print("Format needed for training is not provided in input vcf file.")
            else:
                Chrom = CHR_Reform_vcf(Chr_temp,args.Format,args.hwe,args.maf)
        elif args.geno=='dosages':
            Chrom = CHR_Reform_DS(Chr_temp,args.hwe,args.maf)
        else:
            print("geno file can not identify")


        Info,sampleID = Info_Gene(Chrom,Exp_temp)

        ### Evaluate vhether DPR model vaild
        ### create bimbam file for 5-folds cross validation

        CV_file_dir=args.out_prefix+"/DPR_input/CV"
        CV_output=args.out_prefix+"/DPR_input"

        k_fold_R2=[]
        k_fold_R2.append(0)
        for i in range(5):
            trainID = np.intersect1d(np.array(CV_trainID.loc[i][0]),sampleID)
            testID = np.intersect1d(np.array(CV_testID.loc[i][0]),sampleID)
            
            Info_trainCV = Info >> select(Info.snpID,Info.POS,Info.CHROM,Info.REF,Info.ALT,Info[trainID],Info.TargetID)
            ### bimbam file
            bimbam_train = (Info_trainCV 
                            >> select(Info_trainCV.snpID,Info_trainCV.REF,Info_trainCV.ALT,Info_trainCV[trainID])).dropna(axis=0,how='any')
            bimbam_train.to_csv(CV_file_dir+'/bimbam/'+TargetID+'_CV'+str(i+1)+'_bimbam.txt',
                                header=False,index=None,sep='\t',mode='a')
            ### phenotype file
            pheno_train = (Info_trainCV 
                           >> mask(Info_trainCV.TargetID==TargetID)>> drop(Info_trainCV.TargetID)).dropna(axis=1,how='any')
            pheno_train.T.to_csv(CV_file_dir+'/pheno/'+TargetID+'_CV'+str(i+1)+'_pheno.txt',
                                 header=False,index=None,sep='\t',mode='a')
            ### SNP annotation file
            SNP_annot_train = (Info_trainCV >> select(Info_trainCV.snpID,Info_trainCV.POS,Info_trainCV.CHROM)).dropna(axis=0,how='any')
            SNP_annot_train.to_csv(CV_file_dir+'/SNP_annot/'+TargetID+'_CV'+str(i+1)+'_snp_annot.txt',
                                   header=False,index=None,sep='\t',mode='a')
            ### call DPR
            TargetID_CV = TargetID+'_CV'+str(i+1)
            stop_CV=0
            try:
                subprocess.check_call(shlex.split('./Model_Train_Pred/call_DPR.sh'+' '+CV_file_dir+' '+str(args.dpr)+' '+TargetID_CV+' '+CV_output))
            except subprocess.CalledProcessError as err:
                stop_CV=1
                print("DPR failed in CV"+str(i+1)+" for TargetID:"+TargetID)
                
            if stop_CV==1:
                continue
            else:
                ### Read in cross validation training result
                result_CV = pd.read_csv(CV_output+'/output/DPR_'+TargetID_CV+'.param.txt',sep='\t').rename(columns={'rs':'snpID'})

                ### overall effect size
                if args.ES=='fixed':
                    result_CV['ES'] = result_CV['beta']
                elif args.ES=='additive':
                    result_CV['ES'] = result_CV['b']+result_CV['beta']

                result_CV = result_CV >> select(result_CV.snpID,result_CV.ES)
            
                ### calculate predicted value for test set
                Info_testCV = Info >> select(Info.snpID,Info.POS,Info.CHROM,Info.REF,Info.ALT,Info[testID],Info.TargetID)
                
                bimbam_test = (Info_testCV 
                               >> select(Info_testCV.snpID,Info_testCV[testID])).dropna(axis=0,how='any')
                pheno_test = (Info_testCV 
                              >> mask(Info_testCV.TargetID==TargetID)>> drop(Info_testCV.TargetID)).dropna(axis=1,how='any')
                pheno_test = pheno_test.reset_index(drop=True)
                pheno_test = pd.DataFrame(pheno_test,dtype=np.float) 

                overall = bimbam_test.merge(result_CV,left_on='snpID',right_on='snpID',how='outer')
                overall['ES'].fillna(0,inplace=True)

                Ypred=np.array(mat(pd.DataFrame(overall[testID],dtype=np.float)).T*mat(overall.ES).reshape((len(overall.snpID),1))).ravel()

                lm = sm.OLS(np.array(pheno_test.loc[0]),sm.add_constant(Ypred)).fit()
                k_fold_R2.append(lm.rsquared)
        
        if sum(k_fold_R2)/5 < 0.01:
            print("DPR model is not valid for Gene:"+TargetID)
            print(str(sum(k_fold_R2)/5))
        
        else:
            print("Running DPR training for Gene:"+TargetID)
            print(str(sum(k_fold_R2)/5))
            file_dir=args.out_prefix+'/DPR_input'
            print("Running model training for Gene:"+TargetID)
            bimbam = (Info >> select(Info.snpID,Info.REF,Info.ALT,Info[sampleID])).dropna(axis=0,how='any')
            
            bimbam.to_csv(file_dir+'/bimbam/'+TargetID+'_bimbam.txt',
                          header=False,index=None,sep='\t',mode='a')
            pheno = (Info >> mask(Info.TargetID==TargetID)>> drop(Info.TargetID)).dropna(axis=1,how='any')
            pheno.T.to_csv(file_dir+'/pheno/'+TargetID+'_pheno.txt',
                           header=False,index=None,sep='\t',mode='a')
            
            SNP_annot = (Info >> select(Info.snpID,Info.POS,Info.CHROM)).dropna(axis=0,how='any')
            SNP_annot.to_csv(file_dir+'/SNP_annot/'+TargetID+'_snp_annot.txt',
                             header=False,index=None,sep='\t',mode='a')
            
            stop_DPR=0
            try:
                subprocess.check_call(shlex.split('./Model_Train_Pred/call_DPR.sh'+' '+file_dir+' '+str(args.dpr)+' '+TargetID+' '+args.out_prefix))
            except subprocess.CalledProcessError as err:
                stop_DPR=1
                print("DPR failed for TargetID:"+TargetID)
            
            if stop_DPR==0:
                Info_Train=pd.DataFrame()
                Info_Train['TargetID']=np.array(TargetID).ravel()
                
                result=pd.read_csv(args.out_prefix+'/output/DPR_'+TargetID+'.param.txt',sep='\t')
                result['TargetID']=TargetID

                if args.ES=='fixed':
                    result['ES']=result['beta']
                elif args.ES=='additive':
                    result['ES']=result['beta']+result['b']

                Info_Train['snp_size']=np.array(len(result.ES)).ravel()
                ### only keep snps with ES!=0
                result = result >> mask(result.ES!=0)
                Info_Train['effect_snp_size']=np.array(len(result.ES)).ravel()

                ### output DPR training result to a single file
                result=result.rename(columns={'chr':'CHROM','rs':'snpID','ps':'POS'})

                result=result >> select(result.CHROM,result.snpID,result.POS,result.TargetID,
                                        result.n_miss,result.b,result.beta,result.ES,result.gamma)

                ### Store Filter information
                Filter = Chrom >> select(Chrom.snpID,Chrom.p_HWE,Chrom.MAF)
        
                ### output training parameter file
                param = result.merge((Filter
                                      >> mask(Filter.snpID.isin(np.array(result.snpID)))),
                                     left_on='snpID',right_on='snpID',how='outer')
                param['REF'] = param['snpID'].apply(lambda x:x.split(":")[2])
                param['ALT'] = param['snpID'].apply(lambda x:x.split(":")[3])

                param = param >> select(param.CHROM,param.POS,param.REF,param.ALT,param.TargetID,
                                        param.n_miss,param.b,param.beta,param.ES,param.gamma,param.p_HWE,param.MAF)
                param['CHROM'] = param['CHROM'].astype('int')
                param['POS'] = param['POS'].astype('int')

                param.to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_DPR_training_param.txt',
                             sep='\t',header=None,index=None,mode='a')

                result = result >> select(result.snpID,result.ES)

                ### for R2 calculation
                bimbam = bimbam >> drop(bimbam.REF,bimbam.ALT)

                ### read in phenotype file
                pheno = pd.DataFrame(pheno,dtype=np.float).reset_index(drop=True)
                Info_Train['sample_size']=np.array(len(pheno)).ravel()

                ID = np.intersect1d(np.array(result.snpID),np.array(bimbam.snpID))

                pred_temp = (result
                             >> mask(result.snpID.isin(ID))).merge((bimbam
                                                                    >> mask(bimbam.snpID.isin(ID))),
                                                                   left_on='snpID',right_on='snpID',
                                                                   how='outer')
                pred = pred_temp.T
                pred.columns = pred.loc['snpID']
                pred = (pred.drop(['snpID','ES'])).reset_index(drop=True)
                pred = pd.DataFrame(pred,dtype='float')
                pheno_pred = np.array(mat(pred)*mat(pred_temp.ES).T).ravel()

                lm_final = sm.OLS(np.array(pheno_pred),sm.add_constant(np.array(pheno.loc[0]))).fit()
                Info_Train['5-fold-CV-R2'] = np.array(sum(k_fold_R2)/5).ravel()
                Info_Train['TrainPVALUE'] = np.array(lm_final.f_pvalue).ravel()
                Info_Train['Train-R2'] = np.array(lm_final.rsquared).ravel()

                Info = (Exp_temp
                        >> select(Exp_temp[['CHROM','GeneStart','GeneEnd',
                                            'GeneName','TargetID']])).merge(Info_Train,
                                                                            left_on='TargetID',
                                                                            right_on='TargetID',
                                                                            how='outer')
                Info.to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_DPR_training_info.txt',
                            sep='\t',header=None,index=None,mode='a')

################################################################################################################
### Start thread
if (args.thread < int(len(EXP)/100) | args.thread > len(EXP)):
    args.thread = (int(len(EXP)/100)+1)*100

pool = multiprocessing.Pool(args.thread)

pool.map(thread_process,[num for num in range(len(EXP))])

pool.close()
pool.join()

#################################################################################################################
### time calculation
time=round((time.clock()-start_time)/60,2)

print(str(time)+' minutes')


