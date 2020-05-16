## TIGAR
**TIGAR** stands for **Transcriptome-Integrated Genetic Association Resource**, which is developed using *Python* and *BASH* scripts for conducting **Transcriptome-Wide Association Studies (TWAS)** by training gene expression imputation models by **nonparametric Bayesian Dirichlet Process Regression (DPR)** method with reference dataset. 

1. **TIGAR** can train gene expression imputation models by both **nonparametric Bayesian DPR** and **Elastic-Net (PrediXcan)** methods with reference dataset that contain transcriptomic and genetic data of the same samples.
2. Impute **Genetically Regulated gene eXpression (GReX)** from *Individual-level* genetic data
3. Conduct **TWAS** using both *Individual-level* and *Summary-level* GWAS data for studying *Univariate* and *Multivariate* phenotypes.



### Software Setup

#### 1. Install BGZIP, TABIX, Python 3.5 and the following python libraries
- [BGZIP](http://www.htslib.org/doc/bgzip.html) 
- [TABIX](http://www.htslib.org/doc/tabix.html) 
- Python 3.5 libraries
   - panda, numpy, dfply, io, argparse
   - time, shlex, warnings
   - subprocess, multiprocess
   - sklearn, statsmodels, scipy

#### 2. Setup Executable Files
- Change `*.sh` files under **TIGAR** directory into executable files
```
cd [TIGAR_directory]
chmod 755 *.sh 
```


### Input Files
Example input files provided under `./ExampleData/` are generated artificially. All input files are *Tab Delimited Text Files*.


#### 1. Gene Expression File (`./ExampleData/gene_exp.txt`)
- First 5 columns are *Chromosome number, Gene start position, Gene end position, Target gene ID, Gene name* (optional, could be the same as Target gene ID).
- Gene expression data start from the 6th column. Each column denotes the corresponding gene expression value per sample. 

| CHROM | GeneStart | GeneEnd |   TargetID      | GeneName | Sample1 | Sample...|
|:-----:|:---------:|:-------:|:---------------:|:--------:|:-------:|:--------:|
|   1   |    100    |   200   |     ENSG0000    |     X    |   0.2   |     ...  |


#### 2. Genotype File
- [VCF](http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/) file (`./ExampleData/example.vcf.gz`)
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

#### 3. Phenotype [PED](http://zzz.bwh.harvard.edu/plink/data.shtml#ped) File (`./ExampleData/example_PED.ped`)
- First five columns are *Family ID, Individual ID, Paternal ID, Maternal ID, Sex* (1=male; 2=female; other=unknown). The combination of family and individual ID should uniquely identify a person.
- Phenotype is the *Sixth* column. A PED file must have 1 and only 1 phenotype in the *Sixth* column.  
- Other covariates start from the *Seventh* column

| FAM_ID | IND_ID | FAT_ID | MOT_ID | SEX | PHENO | COV1 | COV...|
|:------:|:------:|:------:|:------:|:---:|:-----:|:---:|:---:|
|   11A  |   11A  |    X   |    X   |  1  |  0.2  | 0.3 |...|

#### 4. Phenotype and Covariate Information File (`./ExampleData/PED_Info_*.txt`)
- Specify column headers of phenotype and covariates for TWAS 
- Two columns with the first column specifying the Phenotype (P) and Covariate variables (C) from the PED file, and the second column specifying the corresponding column headers in the PED file. 

|P|PHENO|
|:-----:|:---:|
|C|COV1|
|C|COV2|
|C|SEX|

#### 5. Zscore File (`./ExampleData/CHR1_GWAS_Zscore.txt.gz`)
- Contain GWAS summary statistics (*Zscore*) for TWAS
- First 4 columns are of the same format as VCF file.
- Zscore statistic value must be provided under the column header of `Zscore`.
- Sorted by Chromosome and Base pair position, bgzipped by `bgzip`, and tabixed by `tabix`. Example tabix command, `tabix -f -p vcf *_Zscore.txt.gz`.


| CHROM | POS | REF | ALT | Zscore |
|:-----:|:---:|:---:|:---:|:------:|
|   1   | 100 |  C  |  T  |  0.01  |

#### 6. Variant Weight (eQTL effect size) File (`./ExampleData/eQTLweights.txt`)
- Contain cis-eQTL effect sizes (i.e., SNP weights) estimated from reference data for TWAS. The output `*_eQTLweights.txt` file by `./TIGAR_Model_Train.sh` can be used here as the variant weight file.
- First 5 columns must be of the following format, specifying *Chromosome number, Base pair position, Reference allele, Alternative allele, and Target gene ID*. 
- Must contain a column named `ES` to denote variant weights (estimated eQTL effect sizes from reference dataset) with respect to the Target gene ID.
- Each row denotes the variant weight (eQTL effect size) for testing the association of the Target gene ID.
- Variants will be matched with those from Zscore file by their unique `CHROM:POS:REF:ALT`.

| CHROM | POS | REF | ALT |     TargetID    |  ES  |
|:-----:|:---:|:---:|:---:|:---------------:|:----:|
|   1   | 100 |  C  |  T  |     ENSG0000    |  0.2 |


#### 7. Gene Annotation File (`./ExampleData/gene_anno.txt`)
- Provide a list of genes for TWAS 
- Same format as the first five columns of the Gene Expression File.

| CHROM | GeneStart | GeneEnd |     TargetID    | GeneName | 
|:-----:|:---------:|:-------:|:---------------:|:--------:|
|   1   |    100    |   200   |     ENSG0000    |     X    |

#### 8. LD File (`./ExampleData/`)
- Provide LD coefficients generated from reference data for TWAS with GWAS summary statistics

#### 9. Genome Block Annotation File (`./ExampleData/block_annotation_EUR.txt`)
- Genome block annotation file need to be specified for generating LD coefficients of your reference data that can be used to conduct TWAS with summary-level GWAS statistics
- A tab delimited text file with 4 columns `CHROM Start End File`, denoting the Chromosome number, Starting position, Ending position, and corresponding reference VCF file name under specified `--geno_path`. 
- Reference VCF files shall be of one file per Chromosome, or one file for all genome-wide variants. Example genome block annotation file for European samples is provided `./TIGAR/ExampleData/block_annotation_EUR.txt`. 

| CHROM |   Start   |    End  |        File     |
|:-----:|:---------:|:-------:|:---------------:|
|   1   |    100    | 20000   |  CHR1.vcf.gz    |

- Block annotation files of other ethnicities can be adopted from the genome segmentation generated by `LDetect`, https://bitbucket.org/nygcresearch/ldetect-data/src/master/.


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
	- `--cvR2`: Take value `1` for evaluating training model by R2 from 5-fold cross validation (only output trained model if R2_CV < 0.005, default) or value `0` to skip this 5-fold cross validation evaluation (output trained model anyway)
	- `--TIGAR_dir` : Specify the directory of TIGAR source code
	- `--thread`: Number of threads for parallel computation (default `1`)
	- `--out_dir`: Output directory (will be created if not exist)


- Train *nonparametric Bayesian DPR* imputation model
	- Variables to specify for training nonparametric Bayesian DPR imputation model
		- `--dpr`: Bayesian inference algorithm used by DPR: `1` (Variational Bayesian, faster but may less accurate) or `2` (MCMC, slower but accurate)
		- `--ES`: Output effect size type: `fixed` (default) for fixed effects or `additive` for additive fixed and random effects
	- Example bash command
```
# Setup input file paths
TIGAR_dir="/home/jyang/GIT/TIGAR"
Gene_Exp_train_file="${TIGAR_dir}/ExampleData/gene_exp.txt"
train_sample_ID_file="${TIGAR_dir}/ExampleData/sampleID.txt"
genofile="${TIGAR_dir}/ExampleData/example.vcf.gz"
out_dir="${TIGAR_dir}/ExampleData/output"

# Call TIGAR model training shell script
${TIGAR_dir}/TIGAR_Model_Train.sh --model DPR \
--Gene_Exp ${Gene_Exp_train_file} \
--train_sampleID ${train_sample_ID_file} \
--genofile ${genofile} --chr 1 \
--genofile_type vcf --Format GT \
--maf 0.01 \
--hwe 0.0001 \
--cvR2 1 \
--dpr 1 --ES fixed \
--TIGAR_dir ${TIGAR_dir} \
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
${TIGAR_dir}/TIGAR_Model_Train.sh --model elastic_net \
--Gene_Exp ${Gene_Exp_train_file} \
--train_sampleID ${train_sample_ID_file} \
--genofile ${genofile} --chr 1 \
--genofile_type vcf --Format GT \
--maf 0.01 \
--hwe 0.0001 \
--cvR2 1 \
--alpha 0.5 \
--TIGAR_dir ${TIGAR_dir} \
--out_dir ${out_dir}
```

- Model training output files
	- `*_eQTLweights.txt` is the file storing all eQTL effect sizes (`ES`) estimated from the gene expression imputation model per gene (`TargetID`)
	- `*_GeneInfo.txt` is the file storing information about the fitted gene expression imputation model per gene (per row), including gene annotation (`CHROM GeneStart GeneEnd TargetID GeneName`), sample size, number of SNP used in the model training (`n_snp`), number of SNPs with non-zero eQTL effect sizes (`n_effect_snp`), imputation R2 by 5-fold cross validation (`CVR2`), imputation R2 using all given training samples (`TrainR2`)
	- `*_Log.txt` is the file storing all log messages for model training.

#### 2. Predict GReX
- Predict GReX value with given variant weights (eQTL effect sizes) from trained gene expression imputation models and individual-level genotype data of test samples
- Variables to specify
	- `--chr`: Chromosome number need to be specified with respect to the genotype input data
	- `--weight`: Path for SNP weight (eQTL effect size) file 
	- `--test_sampleID`: Path for a file with sampleIDs that should be contained in the genotype file
	- `--genofile`: Path for the training genotype file (bgzipped and tabixed) 
	- `--genofile_tye`: Genotype file type: `vcf` or `dosage`
	- `--Format`: Genotype format in VCF file that should be used: `GT` (default) for genotype data or `DS` for dosage data, only required if the input genotype file is of VCF file
	- `--window`: Window size around gene transcription starting sites (TSS) for selecting cis-SNPs for fitting gene expression prediction model (default `1000000` for `+- 1MB` region around TSS)
	- `--maf_diff`: MAF difference threshold for matching SNPs from eQTL weight file and test genotype file. If SNP MAF difference is greater than `maf_diff` (default `0.2`), the SNP will be excluded
	- `--TIGAR_dir` : Specify the directory of TIGAR source code
	- `--thread`: Number of threads for parallel computation (default `1`)
	- `--out_dir`: Output directory (will be created if not exist)
- Example bash command
```
eQTL_ES_file="${TIGAR_dir}/ExampleData/eQTLweights.txt"
gene_anno_file="${TIGAR_dir}/ExampleData/gene_anno.txt"
test_sample_ID_file="${TIGAR_dir}/ExampleData/test_sampleID.txt"

${TIGAR_dir}/TIGAR_GReX_Pred.sh --chr 1 \
--weight ${eQTL_ES_file} \
--gene_anno ${gene_anno_file} \
--test_sampleID ${test_sample_ID_file} \
--genofile ${genofile} \
--genofile_type vcf --Format GT \
--TIGAR_dir ${TIGAR_dir} \
--out_dir ${out_dir}
```
- GReX prediction output
	- `*_Pred_GReX.txt` : File containing predicted GReX values in the same format as gene expression file.
	- `*_Pred_Log.txt` : File containing log messages

#### 3. TWAS
- Variables to specify for TWAS
	- `--asso` : Specify TWAS type either by using individual-level phenotype data and predicted GReX with value `1` (default) or using GWAS summary Z-score statistics with value `2`
	- `--thread` : Number of threads for parallel computation (default `1`)
	- `--out_dir`: Output directory (will be created if not exist)

- Variables to specify for using *individual-level* GWAS data. 
	- `--gene_exp` : Path for predicted GReX file
	- `--PED` : Path for PED file that contains phenotype and covariate data
	- `--PED_info` : Specify culumn names for phenotypes and covariates that will be used in TWAS.
		- Rows started with value `P` : followed with phenotype column names
		- Rows started with value `C` : followed with covariate column names
	- `--method` : `OLS` (default) for studying quantitative phenotype or multivariate phenotypes or `Logit` for studying dichotomous univariate phenotype

```
GReX_pred_file="${TIGAR_dir}/CHR1_Pred_GReX.txt"
PED="${TIGAR_dir}/example_PED.ped"
PED_info="${TIGAR_dir}/PED_Info_SinglePheno.txt"

${TIGAR_dir}/TIGAR_TWAS.sh --asso 1 \
--gene_exp ${GReX_pred_file} \
--PED ${PED} --PED_info ${PED_info} \
--method OLS \
--out ${out_dir}
```

- Variables to specify for using *summary-level* GWAS data.
	- `--gene_anno` : Path for gene annotation file to specify a list of gene for TWAS. The first six columns of the gene expression file can be used here.
	- `--chr`: Chromosome number need to be specified to conduct TWAS per chromosome 
	- `--weight`: Path for SNP weight (eQTL effect size) file. Output file by model training can be used here. 
	- `--Zscore` : Path for GWAS summary Zscore statistics.
	- `--LD` : Path for reference LD (SNP genotype covariance matrix) that should be bgzipped and tabixed. Can be generated by using our `${TIGAR_dir}/TWAS/Get_LD.py` script with individual-level genotype data of reference samples.
	- `--window`: Window size around gene transcription starting sites (TSS) for selecting cis-SNPs for TWAS (default `1000000` for `+- 1MB` region around TSS)
	- `--TIGAR_dir` : Specify the directory of **TIGAR** source code

```
gene_anno_file="${TIGAR_dir}/ExampleData/gene_anno.txt"
Zscore_file="${TIGAR_dir}/CHR1_GWAS_Zscore.txt.gz"
eQTL_ES_file="${TIGAR_dir}/ExampleData/eQTLweights.txt"
LD_file="${TIGAR_dir}/CHR1_LD.txt.gz"
out_dir="${TIGAR_dir}/ExampleData/output"

${TIGAR_dir}/TIGAR_TWAS.sh --asso 2 \
--gene_anno ${gene_anno_file} \
--Zscore ${Zscore_file} \
--weight ${eQTL_ES_file} \
--LD ${LD_file} \
--chr 1 \
--out ${out_dir}
```

#### 4. Generate reference covariance files
```
${TIGAR_dir}/TIGAR_LD.sh --block ${block_annotation} \
--geno_path ${geno_path} --geno vcf \
--chr 1 --Format GT \
--out ${out_prefix}
```

### Reference
- PrediXcan : https://github.com/hakyimlab/PrediXcan  
- DPR: https://github.com/biostatpzeng/DPR
