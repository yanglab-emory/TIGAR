#!/usr/bin/env python

##############################################################################################
# Import packages needed
import argparse
import warnings
import time
import subprocess
from subprocess import *
import io
from io import StringIO
from io import *
import pandas as pd
import numpy as np
from numpy import *
from dfply import *
import multiprocessing

### import grid search for model selection
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold


### For Elastic Net Regression
from sklearn.linear_model import ElasticNet
from sklearn.metrics import r2_score

### For OLS regression in cross validation
import statsmodels.api as sm
from scipy import stats

###############################################################################################
### time calculation
start_time=time.clock()

###############################################################################################
warnings.filterwarnings("ignore")
###############################################################################################
# Construct Dataframe for Analysis
###############################################################################################

# Reform vcf file
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

################################################################################################
### Calculating p-value for Hardy Weinberg Equilibrium exact test

### Input value:
### 1.obs_hets: Observed heterozygosity = Number of 1 in each SNPs(i.e. 0.5 <= gij < 1.5)
### 2.obs_hom1: Observed AA homozygosity = Number of 0 in each SNPs(i.e. 0 <= gij< 0.5)
### 3.obs_hom2: Observed aa homozygosity = Number of 2 in each SNPs(i.e. 1.5 <= gij <= 2)

### gij denote value for gene variance for jth samples in ith SNPs
### 0 <= gij< 0.5 denote as 0
### 0.5 <= gij < 1.5 denote as 1
### 1.5 <= gij <2 denote as 2

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

#######################################################################################################
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
    sampleID = data.columns[9:]

    data['snpID']=(data['CHROM'].astype('str')+":"+data['POS'].astype('str')
                   +":"+data.REF+":"+data.ALT)
        
    CHR = data >> select(data[['CHROM','POS','ID','REF','ALT','snpID']],data[sampleID])

    CHR=CHR.drop_duplicates(['snpID'],keep='first')
        
    indicate=data.FORMAT[0].split(":").index(Format)
    CHR[sampleID]=CHR[sampleID].applymap(lambda x:x.split(":")[indicate])
    
    CHR[sampleID]=CHR[sampleID].apply(lambda x:geno_reform(x,Format),axis=0)
        
    ### Calculating MAF and p_HWE by SNPs
    temp=pd.DataFrame((CHR >> select(CHR[sampleID])),dtype=np.float)
    
    ### Calculating p_HWE
    CHR['p_HWE']=temp.apply(lambda x:p_HWE(x.dropna()),axis=1)
    
    ### Calculate MAF(range from 0-1)
    CHR['MAF']=temp.apply(lambda x:sum(x)/(2*len(x.dropna())),axis=1)
    ### Dealing with NaN
    CHR[np.hstack(([sampleID,'MAF']))] = CHR[np.hstack(([sampleID,'MAF']))].apply(lambda x:x.fillna(2*x.MAF),axis=1)
    
    return CHR >> mask((CHR.p_HWE>=p_hwe)&(CHR.MAF>=maf))

### For dosages input
### Input:
### 1. data:The first five columns fixed
###    1) CHROM
###    2) POS
###    3) ID (i.e. rsID)
###    4) REF
###    5) ALT
### 2. Format: DS
### 3. Output dataframe will be selected by p_HWE and MAF:
###    Default threshold is p_HWE > 10**(-3) and MAF > 0.01
def CHR_Reform_DS(data,p_hwe,maf):
    sampleID=data.columns[5:]
    
    data['snpID']=(data['CHROM'].astype('str')+':'+data['POS'].astype('str')
                   +':'+data.REF+':'+data.ALT)
    data=data.drop_duplicates(['snpID'],keep='first')

    data[data[sampleID].astype('str')=='.']=nan
    
    data[sampleID]=data[sampleID].astype('float')
    
    data['p_HWE']=data[sampleID].apply(lambda x:p_HWE(x.dropna()),axis=1)
    
    data['MAF']=data[sampleID].apply(lambda x:sum(x)/(2*len(x.dropna())),axis=1)
    
    data[np.hstack(([sampleID,'MAF']))] = data[np.hstack(([sampleID,'MAF']))].apply(lambda x:x.fillna(2*x.MAF),axis=1)
    
    return data >> mask((data.p_HWE>=p_hwe)&(data.MAF>=maf))

##############################################################################################

### Merge SNPs information with Gene Expression Level by sampleID

### Output dataframe:
### sample ID as index
### SNps & TargetID as columns
def Info_Gene(Chr,Exp):
    sampleID_temp = np.intersect1d(np.array(Chr.columns),np.array(Exp.columns))
    
    sampleID=np.delete(sampleID_temp,np.where(sampleID_temp=='CHROM'))
    
    CHR_target = Chr >> select(Chr.snpID,Chr[sampleID])
    Exp_target = Exp >> select(Exp.TargetID,Exp[sampleID])
    
    info_temp = pd.concat([CHR_target,Exp_target],sort=False).T
    
    snpID = np.array(info_temp.loc['snpID'].dropna())
    
    TargetID = Exp.TargetID
    
    info_temp.columns = np.hstack((snpID,TargetID))
    
    info = pd.DataFrame(info_temp.drop(['snpID','TargetID']),
                        dtype=np.float)
    
    return info,snpID,len(sampleID)

########################################################################################
# Model Selection

### 5-folds cross validation R2 calculation
def CV_R2(trainX,trainY,testX,testY,k,Alpha):
    clf = GridSearchCV(ElasticNet(l1_ratio=Alpha,fit_intercept=False),
                       [{'alpha':np.arange(0,1.01,0.01)}],cv=k).fit(trainX,trainY)
    
    reg = ElasticNet(l1_ratio=Alpha,alpha=clf.best_params_['alpha']).fit(trainX,trainY)
    
    lm = sm.OLS(testY,sm.add_constant(reg.predict(testX))).fit()
    
    return lm.rsquared

### Elastic Net
### Input:
### 1.trainX: independent variables values for training(SNPs Genotype)
### 2.trainY: Response variable values(Gene Expression Level)
### 3.Alpha: ratio for L1 and L2 penalty in elastic net regression,default=0.5
###          Alpha=0: Lasso Regression
###          Alpha=1: Ridge Regression
###          0 < Alpha < 1: Elastic Net Regression
### 4.k: k-fold cross validation,default=5

### Return:
### 1.Regression coefficent for training data
### 2.Training R-square
### 3.Parameter selected by cross validation
### 4.Corresponding mean cross validation score

### Using grid search and cross-validation to find the best lambda(penalty)

def elastic_net(trainX,trainY,k,Alpha):
    clf = GridSearchCV(ElasticNet(l1_ratio=Alpha,fit_intercept=False),
                       [{'alpha':np.arange(0,1.01,0.01)}],cv=k).fit(trainX,trainY)
    
    reg = ElasticNet(l1_ratio=Alpha,alpha=clf.best_params_['alpha']).fit(trainX,trainY)
    Ypred = reg.predict(trainX)
    
    lm = sm.OLS(trainY,sm.add_constant(Ypred)).fit()

    return reg.coef_,r2_score(trainY,Ypred),lm.f_pvalue,clf.best_params_['alpha'],clf.best_score_

#########################################################################################
# Model Training

### variables need
parser = argparse.ArgumentParser(description='manual to this script')

### for Gene annotation and Expression level file
parser.add_argument('--Gene_Exp_path',type=str,default = None)

### for training sampleID
parser.add_argument('--train_sample',type=str,default = None)

### specified chromosome number
parser.add_argument('--chr_num',type=int,default = None)

### number of thread
parser.add_argument('--thread',type=int,default = None)

### training input
parser.add_argument('--train_dir',type=str,default = None)
parser.add_argument('--train_names',type=str,default = None)

### specified input file type(vcf or dosages)
parser.add_argument('--geno',type=str,default = None)

### 'DS' or 'GT'
parser.add_argument('--Format',type=str,default=None)

### for data selection
### Folded Minor Allele Frequency (range from 0-0.5)
parser.add_argument('--maf',type=float,default=None)
### p-value for Hardy Weinberg Equilibrium exact test
parser.add_argument('--hwe',type=float,default=None)
### window
parser.add_argument('--window',type=int,default=None)

### for model trainin
### k-fold cross validation
parser.add_argument('--cv',type=int,default=None)
### Ratio of L1 and L2
parser.add_argument('--alpha',type=float,default=None)

### output dir
parser.add_argument('--out_prefix',type=str,default=None)

args = parser.parse_args()
##############################################################################################
# variable checking

### Print all variables for model training
print("Gene annotation and Expressiop level file:"+args.Gene_Exp_path)
print("Training sampleID path:"+args.train_sample)
print("Chrmosome number:"+str(args.chr_num))
print("Number of thread:"+str(args.thread))
print("Training file:"+args.train_dir)
print(args.train_names)

if args.geno=='vcf':
    print("Using "+args.Format+" Format for Training")
elif args.geno=='dosages':
    print("Using DS Format for Training.")
else:
    raise SystemExit("Geno file can not identify.")


print("output dir:"+args.out_prefix)
print("Threshold for MAF:"+str(args.maf))
print("Threshold for p-value of HW test:"+str(args.hwe))
print("window="+str(args.window))

print("Using "+str(args.cv)+"-fold for cross-validation")
print("L1 & L2 ratio in elastic net regression:"+str(args.alpha))

########################################################################################

### Training Processing

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

### Initialized output
output_param=pd.DataFrame(columns=['CHROM','POS','ID','REF','ALT','TargetID','MAF','p_HWE','ES'])

output_param.to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_elastic_net_training_param.txt',
                    header=True,index=None,sep='\t',mode='a')

output_Info=pd.DataFrame(columns=['CHROM','GeneStart','GeneEnd','TargetID','GeneName',
                                  'sample_size','snp_size','effect_snp_size','5-fold-CV-R2','TrainPVALUE','Train-R2',
                                  'k-fold','alpha','Lambda','cvm'])

output_Info.to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_elastic_net_training_info.txt',
                   header=True,index=None,sep='\t',mode='a')

### Read in header for vcf file
train_names=pd.read_csv(args.train_names,sep='\t').rename(columns={'#CHROM':'CHROM'})

def thread_process(num):
    Exp_temp=pd.DataFrame(EXP.loc[num]).T
    TargetID = np.array(Exp_temp.TargetID)[0]

    start=max(int(Exp_temp.GeneStart)-args.window,0)
    end=max(int(Exp_temp.GeneEnd)+args.window,0)

    ### select corresponding vcf file by tabix
    train_process=subprocess.Popen(["tabix"+" "+args.train_dir+" "+str(args.chr_num)+":"+str(start)+"-"+str(end)],
                                   shell=True,
                                   stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE)
    out=train_process.communicate()[0]
    
    if len(out)==0:
        print("No Training Data for this Gene:"+TargetID)

    else:
        ### Recode subprocess output in 'utf-8'
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
        
        print("Running training for Gene:"+TargetID)
        Information,SNPs,sample_size = Info_Gene(Chrom,Exp_temp)
        
        ### Evaluate whether Elastic Net model valid
        k_fold_R2=[]
        for i in range(5):
            Info_CV_train=Information.loc[np.array(CV_trainID.loc[i][0])].dropna()
            Info_CV_test=Information.loc[np.array(CV_testID.loc[i][0])].dropna()

            k_fold_R2.append(CV_R2(Info_CV_train[SNPs],Info_CV_train[TargetID],
                             Info_CV_test[SNPs],Info_CV_test[TargetID],args.cv,args.alpha))

        if sum(k_fold_R2)/5 < 0.01:
            print("Elastic Net model is not valid for Gene:"+TargetID)
            print(str(sum(k_fold_R2)/5))

        else:
            print(str(sum(k_fold_R2)/5))
            Beta_temp=pd.DataFrame()
            
            Beta_temp['SNPs']=SNPs
            Beta_temp['TargetID']=TargetID
            Beta_temp['beta'],R2,Pvalue,Lambda,cvm=elastic_net(Information[SNPs],
                                                               Information[TargetID],
                                                               args.cv,args.alpha)

            Beta=(Chrom[['CHROM','POS','ID','REF','ALT','snpID',
                         'p_HWE','MAF']].merge(Beta_temp,left_on='snpID',
                                               right_on='SNPs',how='outer'))
            
            Beta=Beta >> drop(Beta[['SNPs','snpID']])

            Beta = (Beta
                    >> mask(Beta.beta != 0) 
                    >> select(Beta.CHROM,Beta.POS,Beta.ID,Beta.REF,
                              Beta.ALT,Beta.TargetID,Beta.MAF,Beta.p_HWE,Beta.beta))
            Beta=Beta.rename(columns={'beta':'ES'})
        
            Beta.to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_elastic_net_training_param.txt',
                        header=False,index=None,sep='\t',mode='a')
        
            ### Store training information
            ### Store result from elastic net
            Info_Train=pd.DataFrame()
        
            Info_Train['TargetID']=np.array(TargetID).ravel()
            Info_Train['sample_size']=sample_size
            Info_Train['snp_size']=len(SNPs)
            Info_Train['effect_snp_size']=len(Beta.ES)
            Info_Train['5-fold-CV-R2'] = np.array(sum(k_fold_R2)/5).ravel()
            Info_Train['TrainPVALUE'] = np.array(Pvalue).ravel()
            
            if len(Beta.ES)==0:
                Info_Train['Train-R2']=0
            else:
                Info_Train['Train-R2']=R2
            
            Info_Train['k_fold']=args.cv
            Info_Train['alpha']=args.alpha
            Info_Train['Lambda']=Lambda
            Info_Train['cvm']=cvm
            
            Info = (Exp_temp >> select(Exp_temp[['CHROM','GeneStart','GeneEnd',
                                                 'GeneName','TargetID']])).merge(Info_Train,
                                                                                 left_on='TargetID',
                                                                                 right_on='TargetID',
                                                                                 how='outer')
            Info.to_csv(args.out_prefix+'/CHR'+str(args.chr_num)+'_elastic_net_training_info.txt',
                        header=None,index=None,sep='\t',mode='a')

###############################################################################################################
### Start thread
if (args.thread < int(len(EXP)/100) | args.thread > len(EXP)):
    args.thread = (int(len(EXP)/100)+1)*100

pool = multiprocessing.Pool(args.thread)

pool.map(thread_process,[num for num in range(len(EXP))])

pool.close()
pool.join()

####################################################################################################################
### time calculation
time=round((time.clock()-start_time)/60,2)

print(str(time)+' minutes')
