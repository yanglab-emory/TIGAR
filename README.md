## TIGAR
**TIGAR** stands for **Transcriptome-Integrated Genetic Association Resource**, which is developed using *Python* and *BASH* scripts for conducting **Transcriptome-Wide Association Studies (TWAS)** by training gene expression imputation models by **nonparametric Bayesian Dirichlet Process Regression (DPR)** method with reference dataset. 

1. **TIGAR** can train gene expression imputation models by both **nonparametric Bayesian DPR** and **Elastic-Net (PrediXcan)** methods with reference dataset that contain transcriptomic and genetic data of the same samples.
2. Impute **Genetically Regulated gene eXpression (GReX)** from *Individual-level* genetic data
3. Conduct **TWAS** using both *Individual-level* and *Summary-level* GWAS data for studying *Univariate* and *Multivariate* phenotypes.

### Software Setup

#### 1. Setup Executable Files
- Change all `*.sh` and `*.py` files into executable files
```
cd [TIGAR_directory]
chmod 755 *.sh ./Model_Train_Pred/*.sh ./Model_Train_Pred/*.py ./Model_Train_Pred/DPR ./TWAS/*.sh ./TWAS/*.py
```
- Add the executable file `./Model_Train_Pred/DPR` to your `${PATH}` directory. 
	Assuming `~/bin/` is a directory added to your `${PATH}`, you can accommodate the following example command
```
cp ./Model_Train_Pred/DPR ~/bin/
```

#### 2. Install BGZIP, TABIX, Python 3.5 and the following python libraries
- [BGZIP](http://www.htslib.org/doc/bgzip.html) 
- [TABIX](http://www.htslib.org/doc/tabix.html) 
- Python 3.5 libraries
   - dfply
   - io
   - subprocess
   - multiprocess

### Input Files
Example input files provided under `./example_data/` are generated artificially. All input files are *Tab Delimited Text Files*.


#### 1. Gene Expression File (`./example_data/Gene_Exp.txt`)
- First 5 columns are *Chromosome number, Gene start position, Gene end position, Target gene ID, Gene name* (optional, could be the same as Target gene ID).
- Gene expression data start from the 6th column. Each column denotes the corresponding gene expression value per sample. 

| CHROM | GeneStart | GeneEnd |   TargetID      | GeneName | Sample1 | Sample...|
|:-----:|:---------:|:-------:|:---------------:|:--------:|:-------:|:--------:|
|   1   |    100    |   200   |     ENSG0000    |     X    |   0.2   |     ...  |


#### 2. Genotype File
- [VCF](http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/) file (`./example_data/example.vcf.gz`)
	- Sorted by chromosome and base pair position, bgzipped by `bgzip`, and tabixed by `tabix`.
	- Example tabix command for VCF file: `tabix -f -p vcf *.vcf.gz`.
	- First 9 columns are Chromosome number, Base pair position, Variant ID, Reference allele, Alternative allele, Quality score, Filter status, Variant information, Genotype data format
	- Genotype data start from the 10th column. Each column denotes the genotype data per sample.

	| CHROM | POS |  ID | REF | ALT | QUAL | FILTER | INFO | FORMAT |  Sample1 | Sample...|
	|:-----:|:---:|:---:|:---:|:---:|:----:|:------:|:----:|:------:|:--------:|:--------:|
	|   1   | 100 | rs1 |  C  |  T  |   .  |  PASS  |   .  |  GT:DS | 0/0:0.01 |    ...   |

- Dosage file
	- The first 5 columns are of the same format as VCF file.
	- Dosage genotype data start from the 6th column that are values ranging from 0 to 2, denoting the number of minor alleles. Each column denotes the genotype data per sample.

	| CHROM | POS |  ID | REF | ALT | Sample1 | Sample...|
	|:-----:|:---:|:---:|:---:|:---:|:-------:|:--------:|
	|   1   | 100 | rs** |  C  |  T  |   0.01  |    ...   |

#### 3. Phenotype [PED](http://zzz.bwh.harvard.edu/plink/data.shtml#ped) File (`./example_data/example_PED.ped`)
- First five columns are *Family ID, Individual ID, Paternal ID, Maternal ID, Sex* (1=male; 2=female; other=unknown). The combination of family and individual ID should uniquely identify a person.
- Phenotype is the *Sixth* column. A PED file must have 1 and only 1 phenotype in the *Sixth* column.  
- Other covariates start from the *Seventh* column

| FAM_ID | IND_ID | FAT_ID | MOT_ID | SEX | PHENO | COV1 | COV...|
|:------:|:------:|:------:|:------:|:---:|:-----:|:---:|:---:|
|   11A  |   11A  |    X   |    X   |  1  |  0.2  | 0.3 |...|

#### 4. Association Study Information File (`./example_data/Asso_Info_*.txt`)
- Specify column headers of phenotype and covariates for TWAS 
- Two columns with the first column specifying the Phenotype (P) and Covariate variables (C) from the PED file, and the second column specifying the corresponding column headers in the PED file. 

|P|PHENO|
|:-----:|:---:|
|C|COV1|
|C|COV2|
|C|SEX|

#### 5. Zscore File (`./example_data/CHR1_GWAS_Zscore.txt.gz`)
- Contain GWAS summary statistics (*Zscore*) for TWAS
- First 4 columns are of the same format as VCF file.
- Zscore statistic value must be provided under the column header of `Zscore`.
- Sorted by Chromosome and Base pair position, bgzipped by `bgzip`, and tabixed by `tabix`. Example tabix command, `tabix -f -p vcf *_Zscore.txt.gz`.


| CHROM | POS | REF | ALT | Zscore |
|:-----:|:---:|:---:|:---:|:------:|
|   1   | 100 |  C  |  T  |  0.01  |

#### 6. Variant Weight File (`./example_data/weight.txt`)
- Contain cis-eQTL effect sizes (i.e., SNP weights) estimated from reference data for TWAS
- First 5 columns have to be of the following format, specifying *Chromosome number, Base pair position, Reference allele, Alternative allele, and Target gene ID*. 
- Column `ES` denotes variant weights (estimated eQTL effect sizes from reference dataset) with respect to the Target gene ID.
- Each row denotes the variant weight per variant for testing the association of the Target gene ID.
- Variants will be matched with those from Zscore file by their unique `CHROM:POS:REF:ALT`.
- Headers of the first 5 columns and variant weights must be the same as shown below:

| CHROM | POS | REF | ALT |     TargetID    |  ES  |
|:-----:|:---:|:---:|:---:|:---------------:|:----:|
|   1   | 100 |  C  |  T  |     ENSG0000    |  0.2 |


#### 7. Gene Annotation File (`./example_data/Gene_annotation.txt`)
- Provide a list of genes for TWAS 
- Same format as the first five columns of the Gene Expression File.

| CHROM | GeneStart | GeneEnd |     TargetID    | GeneName | 
|:-----:|:---------:|:-------:|:---------------:|:--------:|
|   1   |    100    |   200   |     ENSG0000    |     X    |

#### 8. Genome Block Annotation File (`./example_data/block_annotation_EUR.txt`)
- Genome block annotation file need to be specified for generating LD coefficients of your reference data that can be used to conduct TWAS with summary-level GWAS statistics
- A tab delimited text file with 4 columns `CHROM Start End File`, denoting the Chromosome number, Starting position, Ending position, and corresponding reference VCF file name under specified `--geno_path`. 
- Reference VCF files shall be of one file per Chromosome, or one file for all genome-wide variants. Example genome block annotation file for European samples is provided `./TIGAR/example_data/block_annotation_EUR.txt`. 

| CHROM |   Start   |    End  |        File     |
|:-----:|:---------:|:-------:|:---------------:|
|   1   |    100    | 20000   |  CHR1.vcf.gz    |

- Block annotation files of other ethnicities can be adopted from the genome segmentation generated by `LDetect`, https://bitbucket.org/nygcresearch/ldetect-data/src/master/.

#### 9. LD File
- Provide LD coefficients generated from reference data for TWAS with GWAS summary statistics



### Example Usage 
#### 1. Train gene expression imputation models per Chromosome

- Variables to specify
	- `--model`: Gene expression imputation model: `elastic_net` or `DPR`
	- `--Gene_Exp`: Path for Gene annotation and Expression file
	- `--train_sampleID`: Path for a file with sampleIDs that will be used for training
	- `--genofile`: Path for the training genotype file (bgzipped and tabixed) 
	- `--chr`: Chromosome number need to be specified with respect to the genotype input data (default: `1`)
	- `--genofile_tye`: Genotype file type: `vcf` or `dosage`
	- `--Format`: Genotype format in VCF file that should be used: `GT` (default) for genotype data or `DS` for dosage data, only required if the input genotype file is of VCF file
	- `--maf`: Minor Allele Frequency threshold (ranges from 0 to 1; default `0.01`) to exclude rare variants
	- `--hwe`: Hardy Weinberg Equilibrium p-value threshold (default `0.00001`) to exclude variants that violated HWE
	- `--window`: Window size around gene transcription starting sites (TSS) for selecting cis-SNPs for fitting gene expression imputation model (default `1000000` for +- 1MB region around TSS)
	- `--thread`: Number of threads for parallel computation (default `1`)
	- `--out_dir`: Output directory (will be created if not exist)


- Train *nonparametric Bayesian DPR* imputation model
	- Variables to specify for training nonparametric Bayesian DPR imputation model
		- `--dpr`: Bayesian inference algorithm used by DPR: `1` (Variational Bayesian, faster but may less accurate) or `2` (MCMC, slower but accurate)
		- `--ES`: Output effect size type: `fixed` (default) for fixed effects or `additive` for additive fixed and random effects
	- Example bash command
```
# Setup input file paths
Gene_Exp_train_file="./example_data/Gene_Exp.txt"
train_sample_ID_file="./example_data/sampleID.txt"
genofile="./example_data/example.vcf.gz"
out_dir="./example_data/output"

# Call TIGAR model training shell script
./TIGAR_Model_Train.sh --model DPR \
--Gene_Exp ${Gene_Exp_train_file} \
--train_sampleID ${train_sample_ID_file} \
--genofile ${genofile} --chr 1 \
--genofile_type vcf --Format GT \
--maf 0.01 \
--hwe 0.0001 \
--out_dir ${out_dir}
```

- Train *Elastic-Net* imputation model
	- Variables to specify for training Elastic-Net imputation model
		- `--cv`: Number of cross validation folds for tuning elastic-net penalty parameter (default `5`)
		- `--alpha`: Fixed L1 & L2 penalty ratio for elastic-net model (default `0.5`)
			- If alpha=0, equivalent to lasso regression
			- If alpha=1, equivalent to ridge regression

	- Example bash command
```
./TIGAR_Model_Train.sh --model elastic_net \
--Gene_Exp ${Gene_Exp_train_file} --train_sample ${train_sample_path} \
--chr 1 --train_dir ${train_dir} \
--geno_train vcf --FT DS \
--out ${out_prefix}
```


#### 2. Predict GReX
```
./TIGAR_Model_Pred.sh --chr 1 \
--train_result_path ${train_result_path} \
--train_info_path ${train_info_path} \
--genofile_dir ${genofile_dir} \
--genofile_type vcf --Format GT \
--out ${out_prefix}
```

#### 3. TWAS

Using individual-level GWAS data. Take the output `*_GReX_prediction.txt` from gene expression prediction as the input for `--Gene_EXP` here. 
```
./TIGAR_TWAS.sh --asso 1 \
--Gene_EXP ${Gene_Exp_prediction_file} --PED ${PED} --Asso_Info ${asso_Info} \
--out ${out_prefix}
```

Using summary-level GWAS data. Take the output `*_training_param.txt` from imputation model training as the input Weight file here. The first five columns of the gene expression file will be taken as gene annotation file here for `--Gene_anno`. The same gene expression file can be used as input for `--Gene_anno`. 

```
./TIGAR_TWAS.sh --asso 2 \
--Gene_anno ${Gene_anno_file} --Zscore ${Zscore} --Weight ${Weight} \
--Covar ${Ref_Covariance_file} --chr 22 \
--out ${out_prefix}
```

#### 4. Generate reference covariance files
```
.TIGAR_Covar.sh --block ${block_annotation} \
--geno_path ${geno_path} --geno vcf \
--chr 22 --Format GT \
--out ${out_prefix}
```

### Reference
- PrediXcan : https://github.com/hakyimlab/PrediXcan  
- DPR: https://github.com/biostatpzeng/DPR
