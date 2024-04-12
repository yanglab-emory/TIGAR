# TIGAR
**TIGAR** stands for **Transcriptome-Integrated Genetic Association Resource**, which is developed using *Python* and *BASH* scripts for conducting **Transcriptome-Wide Association Studies (TWAS)** by training gene expression imputation models by **nonparametric Bayesian Dirichlet Process Regression (DPR)** and **Elastic-Net (PrediXcan)** methods with reference transcriptomic panels. 

1. **TIGAR** can train gene expression imputation models by both **nonparametric Bayesian DPR** and **Elastic-Net (PrediXcan)** methods with reference dataset that contain transcriptomic and genetic data of the same samples.
2. Impute **Genetically Regulated gene eXpression (GReX)** from *Individual-level* genetic data.
3. Conduct **TWAS** using both *Individual-level* and *Summary-level* GWAS data for studying *Univariate* and *Multivariate* phenotypes. Both *FUSION* and *S-PrediXcan* Z-score test statistics will calculated in the TWAS results.
4. Conduct **VC-TWAS** using *Individual-level* and *Summary-level* GWAS data using cis-eQTL effect sizes estimated from reference panels for *Univariate* phenotypes.
5. **Note:** please use our most recently updated scripts on Github to conduct TWAS.
6. *Cis-eQTL Effect Sizes* (i.e., weights) estimated from GTEx V8 reference data of 48 tissue types by *Nonparametric Bayesian DPR* method, and *Reference LD Files* from GTEx V8 EUR samples, are available from [Synapse:syn16804296](https://www.synapse.org/#!Synapse:syn16804296/wiki/611467).
7. Please cite our TIGAR papers:
	- [*TIGAR: An Improved Bayesian Tool for Transcriptomic Data Imputation Enhances Gene Mapping of Complex Traits.* AJHG, 2019. Volume 105, ISSUE 2, P267-282.](https://www.cell.com/ajhg/fulltext/S0002-9297(19)30205-8)
	- [*TIGAR-V2: Efficient TWAS tool with nonparametric Bayesian eQTL weights of 49 tissue types from GTEx V8.* HGG Advances, 2022. Volume 3, Issue 1.](https://doi.org/10.1016/j.xhgg.2021.100068)

---

![TIGAR framework flowchart including TWAS steps of training gene expression prediction models from reference data, predicting GReX with individual-level GWAS data, and testing gene-based association with both individual-level and summary-level GWAS data.](TIGAR_flow.png?raw=true)

---

- [Software Setup](#software-setup)
- [Input Files](#input-files)
	- [1. Genotype File](#1-genotype-file)
		- [VCF](#vcf-variant-call-format)
		- [Dosage](#dosage-file)
	- [2. SampleID File](#2-sampleid-file)
	- [3. Gene Expression File](#3-gene-expression-file)
	- [4. Gene Annotation File](#4-gene-annotation-file)
	- [5. Weight (eQTL effect size) File](#5-weight-eqtl-effect-size-file)	
	- [6. PED](#6-ped)
	- [7. PED info](#7-ped-info)
	- [8. Genome Block Annotation File](#8-genome-block-annotation-file)
	- [9. LD File](#9-ld-file)
	- [10. Zscore File](#10-zscore-file)
<!-- - [Commands](#commands) -->
- [Example Usage](#example-usage)
	- [1. Train Gene Expression Imputation Models](#1-train-gene-expression-imputation-models-per-chromosome)
		- [DPR Model](#train-nonparametric-bayesian-dpr-imputation-model)
		- [Elastic-Net Model](#train-elastic-net-imputation-model)
	- [2. Predict GReX](#2-predict-grex)
	- [3. TWAS](#3-twas)
		- [TWAS with individual-level GWAS data](#twas-using-individual-level-gwas-data)
		- [TWAS with summary-level GWAS data](#twas-using-summary-level-gwas-data)
	- [4. Generate Reference LD Genotype Covariance Files](#4-generate-reference-ld-genotype-covariance-files)
	- [5. VC-TWAS](#5-vc-twas)
		- [VC-TWAS with individual-level GWAS data](#vc-twas-using-individual-level-gwas-data)
		- [VC-TWAS with summary-level GWAS data](#vc-twas-using-summary-level-gwas-data)
	- [6. BGWTWAS](#6-bgw-twas)
- [Updates](#updates)
- [Reference](#reference)

---

## Software Setup

### 1. Install BGZIP, TABIX, Python 3.5, and the following Python libraries
- [BGZIP](http://www.htslib.org/doc/bgzip.html) 
- [TABIX](http://www.htslib.org/doc/tabix.html) 
- Python 3.5 modules/libraries:
	- [pandas](https://pandas.pydata.org/)
	- [numpy](https://numpy.org/)
	- [scipy](https://www.scipy.org/)
	- [sklearn](https://scikit-learn.org/stable/)
	- [statsmodels](https://www.statsmodels.org/stable/index.html)
   <!-- - Python Standard Library: argparse, io, math, multiprocessing, operator, subprocess, sys, time, warnings -->

#### Python Environment Setup (using conda or pip)

##### Using conda
- Installation:
```bash
# create the environment tigarenv
conda create --name tigarenv python=3.5 pandas numpy scipy scikit-learn statsmodels
# deactivate the conda environment
conda deactivate
```
- Activate Environment:
```bash
# activate the environment
conda activate tigarenv
# set the PYTHONPATH
export PYTHONPATH=${CONDA_PREFIX}/lib/python3.5/site-packages/:$PYTHONPATH
```
After running TIGAR use `conda deactivate` to deactivate the conda environment.


##### Using pip
- Installation:
```bash
# install pip
# install virtualenv
python3 -m pip install --user virtualenv
# cd to preferred install_directory
cd ${install_dir}
# create the virtual environment tigarenv in the current directory
python3 -m virtualenv tigarenv --python=python3.5
# activate the environment
source ${install_dir}/tigarenv/bin/activate
# install the packages
python3 -m pip install numpy==1.15.2 pandas==0.23.4 scikit-learn==0.20.0 scipy==1.1.0 statsmodels==0.9.0
# deactivate the environment
deactivate
```

- Activate Environment:
```bash
# activate the environment
source ${install_dir}/tigarenv/bin/activate
# set the PYTHONPATH
PYTHONPATH=${install_dir}/tigarenv/lib/python3.5/site-packages/:$PYTHONPATH
```
After running TIGAR use `deactivate` to deactivate the environment.


### 2. Make files executable
- Make `*.sh` files in the **TIGAR** directory executable
```bash
TIGAR_dir="/home/jyang/GIT/TIGAR"
chmod 755 ${TIGAR_dir}/*.sh 
```
- Make `${TIGAR_dir}/Model_Train_Pred/DPR` file executable
```bash
chmod 755 ${TIGAR_dir}/Model_Train_Pred/DPR
```

## Input Files
Example input files provided under `./ExampleData/` are generated artificially. All input files are *Tab Delimited Text Files*.


### 1. Genotype File

#### VCF ([Variant Call Format](http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/))

- `--genofile_type vcf`
	- `--format GT`
	- `--format DS`

| CHROM | POS |  ID | REF | ALT | QUAL | FILTER | INFO | FORMAT |  Sample1 | Sample...|
|:-----:|:---:|:---:|:---:|:---:|:----:|:------:|:----:|:------:|:--------:|:--------:|
|   1   | 100 | rs1 |  C  |  T  |   .  |  PASS  |   .  |  GT:DS | 0/0:0.01 |    ...   |

- Sorted by chromosome and base pair position, bgzipped by `bgzip`, and tabixed by `tabix`
- **The values in the CHROM column should not start with "chr"**
- Example tabix command for a VCF file: `tabix -f -p vcf *.vcf.gz`.
- First 9 columns are Chromosome number, Base pair position, Variant ID, Reference allele, Alternative allele, Quality score, Filter status, Variant information, Genotype data format
- Genotype data starts from the 10th column. Each column denotes the genotype data per sample.
- Example: `./ExampleData/example.vcf.gz`



#### Dosage file

- `--genofile_type dosage`

| CHROM | POS |  ID | REF | ALT | Sample1 | Sample...|
|:-----:|:---:|:---:|:---:|:---:|:-------:|:--------:|
|   1   | 100 | rs** |  C  |  T  |   0.01  |    ...   |

- First 5 columns have the same format as the first 5 columns of a VCF file.
- Dosage genotype data start from the 6th column that are values ranging from 0 to 2, denoting the number of minor alleles. Each column denotes the genotype data per sample.


#### Usage
```
--genofile
--genofile_type
--format
```

- Train Gene Expression Imputation Models
- Predict GReX
- Generate Reference LD Genotype Covariance Files
- VC-TWAS


### 2. SampleID File
- Headerless, single-column file containing sampleIDs to use.
- Examples:
	- `./ExampleData/sampleID.txt`
	- `./ExampleData/test_sampleID.txt`

#### Usage
`--train_sampleID`

- Train Gene Expression Imputation Models

`--test_sampleID`

- Predict GReX
- VC-TWAS

`--sampleID`

- Generate Reference LD Genotype Covariance Files


### 3. Gene Expression File

| CHROM | GeneStart | GeneEnd |   TargetID      | GeneName | Sample1 | Sample...|
|:-----:|:---------:|:-------:|:---------------:|:--------:|:-------:|:--------:|
|   1   |    100    |   200   |     ENSG0000    |     X    |   0.2   |     ...  |

- First 5 columns are *Chromosome number, Gene start position, Gene end position, Target gene ID, Gene name* (optional, could be the same as Target gene ID).
- Gene expression data start from the 6th column. Each column denotes the corresponding gene expression value per sample. 
- Example: `./ExampleData/gene_exp.txt`

#### Usage
`--gene_exp`

- Train Gene Expression Imputation Models
- TWAS (with individual-level GWAS data)


### 4. Gene Annotation File

| CHROM | GeneStart | GeneEnd |     TargetID    | GeneName | 
|:-----:|:---------:|:-------:|:---------------:|:--------:|
|   1   |    100    |   200   |     ENSG0000    |     X    |

- Provides a list of genes for TWAS 
- Same format as the first five columns of the Gene Expression File.
- Example: `./ExampleData/gene_anno.txt`

#### Usage
`--gene_anno`

- Predict GReX
- TWAS (with summary-level GWAS data)
- VC-TWAS


### 5. Weight (eQTL effect size) File

| CHROM | POS | REF | ALT |     TargetID    |  ES  |
|:-----:|:---:|:---:|:---:|:---------------:|:----:|
|   1   | 100 |  C  |  T  |     ENSG0000    |  0.2 |

- Contain cis-eQTL effect sizes (i.e., SNP weights) estimated from reference data for TWAS. The output `*_eQTLweights.txt` file by `./TIGAR_Model_Train.sh` can be used here as the variant weight file.
- First 5 columns must be be: *Chromosome number, Base pair position, Reference allele, Alternative allele, and Target gene ID*. 
- Must contain a column named `ES` to denote variant weights (estimated eQTL effect sizes from a reference dataset) with for the gene (TargetID).
- Each row denotes the variant weight (eQTL effect size) for testing the association of the gene (TargetID) with a phenotype
- Variants will be matched with those from the Zscore file by their unique `CHROM:POS:REF:ALT` snpID
- Example: `./ExampleData/eQTLweights.txt`

#### Usage
`--weight`

- Predict GReX
- TWAS (with summary-level GWAS data)
- VC-TWAS


### 6. PED

| FAM_ID | IND_ID | FAT_ID | MOT_ID | SEX | PHENO | COV1 | COV...|
|:------:|:------:|:------:|:------:|:---:|:-----:|:---:|:---:|
|   11A  |   11A  |    X   |    X   |  1  |  0.2  | 0.3 |...|

- Phenotype [PED](http://zzz.bwh.harvard.edu/plink/data.shtml#ped) File
- First five columns are *Family ID, Individual ID, Paternal ID, Maternal ID, Sex* (1=male; 2=female; other=unknown). The combination of family and individual ID should uniquely identify a person.
- Phenotype is the *6th* column. A PED file must have 1 and only 1 phenotype in the *6th* column.  
- Other covariates start from the *7th* column
- Example: `./ExampleData/example_PED.ped`

#### Usage
`--PED`

- TWAS (with individual-level GWAS data)
- VC-TWAS


### 7. PED Info
|     |     |
|:---:|:---:|
| P | PHENO |
| C | COV1  |
| C | COV2  |
| C | SEX   |

- Phenotype and Covariate Information File
- Used to specify phenotype columns and covariate columns in the `PED` file to use for the TWAS.
- Headerless, two-column, tab-delimited file. Specify one column per row. The first column of each row is the type (`P` for phenotype, `C` for covariate) and the second column in each row is the name of the column in the `PED` file.
- Example: `./ExampleData/PED_Info_*.txt`

#### Usage
`--PED_info`

- TWAS (with individual-level GWAS data)
- VC-TWAS


### 8. Genome Block Annotation File

| CHROM |   Start   |    End  |
|:-----:|:---------:|:-------:|
|   1   |    100    | 20000   |

- Genome block annotation file used for generating the reference LD covariance files (the LD file is required for TWAS with summary-level GWAS statistics)
- A tab-delimited text file with 3 columns `CHROM Start End`, denoting the chromosome number, block start position, and block ending position
- Block annotation files can be created from genome segmentation files generated by [LDetect](https://bitbucket.org/nygcresearch/ldetect-data/src/master/)
- Example: `./ExampleData/example_genome_block_CHR1.txt`

#### Usage
`--genom_block`

- Generate Reference LD Genotype Covariance Files


### 9. LD File
- Contains LD covariance coefficients generated from reference data.
- Used for TWAS with GWAS summary statistics
- **LD files for TWAS should be made with the same reference data used for model training.**

#### Usage
`--LD`

- TWAS (with summary-level GWAS data)


### 10. Zscore File

| CHROM | POS | REF | ALT | Zscore |
|:-----:|:---:|:---:|:---:|:------:|
|   1   | 100 |  C  |  T  |  0.01  |

- Contains GWAS summary statistics (*Zscore*) for TWAS
- First 4 columns have the same format as the first 4 columns of a VCF file.
- The name of the column containing the Zscore statistic values to use must be `Zscore`.
- The file must be sorted by chromosome and base pair position, bgzipped by `bgzip`, and tabixed by `tabix`. Example tabix command: `tabix -f -p vcf *_Zscore.txt.gz`.
- Example: `./ExampleData/CHR1_GWAS_Zscore.txt.gz`

#### Usage
`--Zscore`

- TWAS (with summary-level GWAS data)

### 11. GWAS result File

| CHROM | POS | REF | ALT | BETA | SE | 
|:-----:|:---:|:---:|:---:|:------:|
|   1   | 100 |  C  |  T  |  0.01  |  0.52  |

- Contains GWAS summary statistics (*BETA*) and (*SE*) for TWAS
- First 4 columns have the same format as the first 4 columns of a VCF file.
- The name of the column containing the BETA and SE to use must be `BETA` and `SE`.
- The file must be sorted by chromosome and base pair position, bgzipped by `bgzip`, and tabixed by `tabix`. Example tabix command: `tabix -f -p vcf *_GWAS.txt.gz`.
- Example: `./ExampleData/sample_GWAS_Result.txt.gz`

#### Usage
`--GWAS_result`

- VC-TWAS (with summary-level GWAS data)





<!-- - Example: -->

<!-- 
## Commands
### 1. Training 

 -->

## Example Usage 
### 1. Train gene expression imputation models per Chromosome

#### Arguments for both imputation models
- `--model`: Gene expression imputation model: `elastic_net` or `DPR`
- `--gene_exp`: Path to Gene annotation and Expression file
- `--train_sampleID`: Path to a file with sampleIDs that will be used for training
- `--chr`: Chromosome number need to be specified with respect to the genotype input data
- `--genofile`: Path to the training genotype file (bgzipped and tabixed)
- `--genofile_type`: Genotype file type: `vcf` or `dosage`
- `--format`: (Required if `genofile_type` is `vcf`) Genotype data format that should be used: `GT` or `DS`
	- `GT`: genotype data
	- `DS`: dosage data
- `--missing rate`: Missing rate threshold. If the rate of missing values for a SNP exceeds this value the SNP will be excluded. Otherwise, missing values for non-excluded SNPs will be imputed as the mean. (default: `0.2`)
- `--maf`: Minor Allele Frequency threshold (ranges from 0 to 1) to exclude rare variants (default: `0.01`
- `--hwe`: Hardy Weinberg Equilibrium p-value threshold to exclude variants that violated HWE (default: `0.00001`)
- `--window`: Window size (in base pairs) around gene region from which to include SNPs (default: `1000000` [`+- 1MB` region around gene region])
- `--cvR2`: Whether to perform 5-fold cross validation by average R2: `0` or `1` (default: `1`)
	- `0`: Skip 5-fold cross validation evaluation
	- `1`: Do 5-fold cross validation and only do final training if the average R2 of the 5-fold validations is at least >= 0.005
- `--thread`: Number of simultaneous *processes* to use for parallel computation (default: `1`)
- `--out_dir`: Output directory (will be created if it does not exist)
- `--TIGAR_dir`: Specify the directory of **TIGAR** source code

#### Train *nonparametric Bayesian DPR* imputation model

##### Additional arguments
- `--dpr`: Bayesian inference algorithm used by DPR: `1` or `2` (default: `1`)
	- `1`: Variational Bayesian, faster but may be less accurate
	- `2`: MCMC, slower but more accurate
- `--ES`: Output effect size type: `fixed` or `additive` (default: `fixed`)
	- `fixed`: use fixed effects only
	- `additive`: use the sum of fixed and random effects

##### Example Command
```bash
# Setup input file paths
Gene_Exp_train_file="${TIGAR_dir}/ExampleData/gene_exp.txt"
train_sample_ID_file="${TIGAR_dir}/ExampleData/sampleID.txt"
genofile="${TIGAR_dir}/ExampleData/example.vcf.gz"
out_dir="${TIGAR_dir}/ExampleData/output"

# Call TIGAR model training shell script
${TIGAR_dir}/TIGAR_Model_Train.sh \
--model DPR \
--gene_exp ${Gene_Exp_train_file} \
--train_sampleID ${train_sample_ID_file} \
--chr 1 \
--genofile ${genofile} \
--genofile_type vcf \
--format GT \
--maf 0.01 \
--hwe 0.0001 \
--cvR2 1 \
--dpr 1 \
--ES fixed \
--thread 2 \
--out_dir ${out_dir} \
--TIGAR_dir ${TIGAR_dir}
```

#### Train *Elastic-Net* imputation model

##### Additional arguments
- `--cv`: Number of cross validation folds for tuning elastic-net penalty parameter (default: `5`)
- `--alpha`: Fixed L1 & L2 penalty ratio for elastic-net model (default: `0.5`)
	- alpha=0: equivalent to lasso regression
	- alpha=1: equivalent to ridge regression

##### Example Command
```bash
Gene_Exp_train_file="${TIGAR_dir}/ExampleData/gene_exp.txt"
train_sample_ID_file="${TIGAR_dir}/ExampleData/sampleID.txt"
genofile="${TIGAR_dir}/ExampleData/example.vcf.gz"
out_dir="${TIGAR_dir}/ExampleData/output"

${TIGAR_dir}/TIGAR_Model_Train.sh \
--model elastic_net \
--gene_exp ${Gene_Exp_train_file} \
--train_sampleID ${train_sample_ID_file} \
--chr 1 \
--genofile ${genofile} \
--genofile_type vcf \
--format GT \
--maf 0.01 \
--hwe 0.0001 \
--cvR2 1 \
--alpha 0.5 \
--thread 2 \
--out_dir ${out_dir} \
--TIGAR_dir ${TIGAR_dir}
```

#### Output
- `CHR${chr}_${model}_eQTLweights.txt` is the file storing all eQTL effect sizes (`ES`) estimated from the gene expression imputation model per gene (`TargetID`)
- `CHR${chr}_${model}_GeneInfo.txt` is the file storing information about the fitted gene expression imputation model per gene (per row), including gene annotation (`CHROM GeneStart GeneEnd TargetID GeneName`), sample size, number of SNPs used in the model training (`n_snp`), number of SNPs with non-zero eQTL effect sizes (`n_effect_snp`), imputation R2 by 5-fold cross validation (`CVR2`), imputation R2 using all given training samples (`TrainR2`)
- `CHR${chr}_${model}_train_log.txt` is the file storing all log messages for model training.


### 2. Predict GReX
Predict GReX value with given variant weights (eQTL effect sizes) from trained gene expression imputation models and individual-level genotype data of test samples

#### Arguments
- `--gene_anno`: Gene annotation file to specify the list of genes, which is of the same format as the first five columsn of gene expression file 
- `--test_sampleID`:  Path to a file with sampleIDs that should be contained in the genotype file which should be used for the prediction
- `--chr`: Chromosome number need to be specified with respect to the genotype input data
- `--weight`: Path to SNP weight (eQTL effect size) file
- `--genofile`: Path to the training genotype file (bgzipped and tabixed) 
- `--genofile_type`: Genotype file type: `vcf` or `dosage`
- `--format`: (Required if `genofile_type` is `vcf`) Genotype data format that should be used: `GT` or `DS`
	- `GT`: genotype data
	- `DS`: dosage data
- `--window`: Window size (in base pairs) around gene region from which to include SNPs (default: `1000000` [`+- 1MB` region around gene region])
- `--missing rate`: Missing rate threshold. If the rate of missing values for a SNP exceeds this value the SNP will be excluded. Otherwise, missing values for non-excluded SNPs will be imputed as the mean. (default: `0.2`)
- `--maf_diff`: MAF difference threshold. If the difference in MAF between a matching SNP in the eQTL weight file and test genotype file is greater than `maf_diff`, that SNP will be excluded. (default: `0.2`)
- `--thread`: Number of simultaneous *processes* to use for parallel computation (default: `1`)
- `--out_dir`: Output directory (will be created if it does not exist)
- `--TIGAR_dir`: Specify the directory of **TIGAR** source code

#### Output
- `CHR${chr}_Pred_GReX.txt`: File containing predicted GReX values in the same format as the input gene expression file
- `CHR${chr}_Pred_log.txt`: File containing log messages

#### Example Command
```bash
gene_anno_file="${TIGAR_dir}/ExampleData/gene_anno.txt"
test_sample_ID_file="${TIGAR_dir}/ExampleData/test_sampleID.txt"
eQTL_ES_file="${TIGAR_dir}/ExampleData/eQTLweights.txt.gz"

${TIGAR_dir}/TIGAR_GReX_Pred.sh \
--gene_anno ${gene_anno_file} \
--test_sampleID ${test_sample_ID_file} \
--chr 1 \
--weight ${eQTL_ES_file} \
--genofile ${genofile} \
--genofile_type vcf \
--format GT \
--thread 2 \
--out_dir ${out_dir} \
--TIGAR_dir ${TIGAR_dir}
```


### 3. TWAS

#### Arguments for both TWAS types
- `--asso`: Specify TWAS type: `1` or `2`
	- `1`: *individual-level* TWAS using individual-level phenotype data and predicted GReX (generated by `TIGAR_GReX_Pred.sh`)
	- `2`: *summary-level* TWAS using GWAS summary Z-score statistics and reference LD (generated by `TIGAR_LD.sh`)
- `--thread`: Number of simultaneous *processes* to use for parallel computation (default: `1`)
- `--out_dir`: Output directory (will be created if it does not exist)
- `--TIGAR_dir`: Specify the directory of **TIGAR** source code

#### TWAS using *individual-level* GWAS data

##### Additional arguments
- `--gene_exp`: Path to predicted GReX file
- `--PED`: Path to PED file that contains phenotype and covariate data
- `--PED_info`: A two-column, tab-delimited file specifying phenotype columns and covariate columns in the `PED` file to use for the TWAS. Specify one column per row. Each row should start with the type of column (`P` for phenotype, `C` for covariate) and the column name.
- `--method`: Method to use based on type of phenotype(s): `OLS` or `Logit` (default: `OLS`)
	- `OLS`: quantitative phenotype or multivariate phenotypes
	- `Logit`: dichotomous univariate phenotype

##### Output
- `indv_${method}_assoc.txt`
- `indv_OLS_TWAS_log.txt`: File containing log messages

##### Example Command
```bash
GReX_pred_file="${TIGAR_dir}/ExampleData/CHR1_Pred_GReX.txt"
PED="${TIGAR_dir}/ExampleData/example_PED.ped"
PED_info="${TIGAR_dir}/ExampleData/PED_Info_SinglePheno.txt"

${TIGAR_dir}/TIGAR_TWAS.sh \
--asso 1 \
--gene_exp ${GReX_pred_file} \
--PED ${PED} \
--PED_info ${PED_info} \
--method OLS \
--thread 2 \
--out_dir ${out_dir} \
--TIGAR_dir ${TIGAR_dir}
```

#### TWAS using *summary-level* GWAS data

##### Additional arguments
- `--gene_anno`: Path to gene annotation file to specify a list of gene for TWAS. The first six columns of the gene expression file can be used here.
- `--chr`: Chromosome number need to be specified to conduct TWAS per chromosome 
- `--weight`: Path to SNP weight (eQTL effect size) file. Output file `_eQTLweights.txt` from model training can be used here. 
- `--Zscore`: Path to GWAS summary Zscore statistics.
- `--LD`: Path to reference LD (SNP genotype covariance matrix) that should be bgzipped and tabixed. Can be generated by using our `${TIGAR_dir}/TWAS/Get_LD.py` script with individual-level genotype data from a reference sample.
- **LD files for TWAS should be made with the same reference data used for model training.**
- **The tool assumes that the reference/alternate allele orientation of SNPs in the weight file match the LD file. Zscore SNPs are "flipped" to match the weight file.**
- `--window`: Window size (in base pairs) around gene region from which to include SNPs (default: `1000000` [`+- 1MB` region around gene region])
-  `--weight_threshold`: include only SNPs with magnitude of weight greater than this value when conducting TWAS (default: `0`)
- `--test_stat`: burden Z test statistic to calculate: `FUSION`, `SPrediXcan`, or `both` (default `both`)
	- `FUSION`: calculate output Zscore and Pvalue results based on method used by FUSION
	- `SPrediXcan`: calculate output Zscore and Pvalue results based on method used by S-PrediXcan
	- `both`: calculate and output Zscore and Pvalue results for both methods

##### Output
- `CHR${chr}_sumstat_assoc.txt`
- `CHR${chr}_TWAS_log.txt`

##### Example Command
```bash
gene_anno_file="${TIGAR_dir}/ExampleData/gene_anno.txt"
eQTL_ES_file="${TIGAR_dir}/ExampleData/eQTLweights.txt.gz"
Zscore_file="${TIGAR_dir}/ExampleData/CHR1_GWAS_Zscore.txt.gz"
LD_file="${TIGAR_dir}/ExampleData/CHR1_reference_cov.txt.gz"

${TIGAR_dir}/TIGAR_TWAS.sh \
--asso 2 \
--gene_anno ${gene_anno_file} \
--chr 1 \
--weight ${eQTL_ES_file} \
--Zscore ${Zscore_file} \
--LD ${LD_file} \
--thread 2 \
--out_dir ${out_dir} \
--TIGAR_dir ${TIGAR_dir}
```


### 4. Generate reference LD genotype covariance files
Reference LD genotype covaraince file can be generated using reference genotype data, which should be generated per chromosome with corresponding genome segmentation block anntoation file.

**LD files for TWAS should be made with the same reference data used for model training.**

#### Arguments
- `--genom_block`: Path to genome segmentation block annotation based on LD
- `--sampleID`: Path to a file with sampleIDs to be used for generating reference LD files. 
- `--chr`: Specify chromosome number of the given reference genotype file
- `--genofile`: Path to the reference genotype file (bgzipped and tabixed). It is recommended that genotype files with data for multiple chromosomes be split into per-chromosome files.
- `--genofile_type`: Genotype file type: `vcf` or `dosage`
- `--format`: (Required if `genofile_type` is `vcf`) Genotype data format that should be used: `GT` or `DS`
	- `GT`: genotype data
	- `DS`: dosage data
- `--maf`: Minor Allele Frequency threshold (ranges from 0 to 1) to exclude rare variants (default: `0.01`)
- `--thread`: Number of simultaneous *processes* to use for parallel computation (default: `1`)
- `--out_dir`: Output directory (will be created if it does not exist)
- `--TIGAR_dir`: Specify the directory of **TIGAR** source code

#### Output
- `CHR${chr}_reference_cov.txt.gz`
- `CHR${chr}_reference_cov.txt.gz.tbi`
- `CHR${chr}_LD_log.txt`

#### Example Command
```bash
block_annotation="${TIGAR_dir}/ExampleData/example_genome_block_CHR1.txt"
sample_id="${TIGAR_dir}/sampleID.txt"

${TIGAR_dir}/TIGAR_LD.sh \
--genome_block ${block_annotation} \
--sampleID ${sample_id}
--chr 1 \
--genofile ${genofile} \
--genofile_type vcf \
--format GT \
--maf 0.01 \
--thread 2 \
--out_dir ${out_dir} \
--TIGAR_dir ${TIGAR_dir}
```

### 5. VC-TWAS 
VC-TWAS with cis-eQTL effect sizes estimated from nonparametric Bayesian DPR method for Univariate phenotypes.

#### Arguments for both VC-TWAS types
- `--gene_anno`: Gene annotation file to specify the list of genes, which is of the same format as the first five columns of gene expression file 
- `--weight`: Path to SNP weight (eQTL effect size) file. The weight file must be generated using the DPR model.
- `--weight_threshold`: Threshold for estimated cis-eQTL effect sizes, filter SNPs with absolute cis-eQTL effect sizes smaller than threshold
- `--chr`: Chromosome number need to be specified with respect to the genotype input data
- `--window`: Window size around gene region from which to include SNPs (default `1000000` for `+- 1MB` region around gene region)
- `--thread`: Number of simultaneous *processes* to use for parallel computation (default: `1`)
- `--out_dir`: Output directory (will be created if it does not exist)
- `--TIGAR_dir`: Specify the directory of **TIGAR** source code

#### VC-TWAS using *individual-level* GWAS data

##### Additional arguments
- `--PED`: Path to PED file that contains phenotype and covariate data
- `--PED_info`: A two-column, tab-delimited file specifying phenotype columns and covariate columns in the `PED` file to use for the TWAS. Specify one column per row. Each row should start with the type of column (`P` for phenotype, `C` for covariate) and the column name.
- `--test_sampleID`: Path to a file with sampleIDs that should be contained in the genotype file
- `--genofile`: Path to the training genotype file (bgzipped and tabixed) 	
- `--genofile_type`: Genotype file type: `vcf` or `dosage`
- `--format`: (Required if `genofile_type` is `vcf`) Genotype data format that should be used: `GT` or `DS`
	- `GT`: genotype data
	- `DS`: dosage data
- `--maf`: Minor Allele Frequency threshold (ranges from 0 to 1) to exclude rare variants (default: `0.01`))
- `--phenotype_type`: phenotype type: `C` or `D`
	- `C`: continous phenotype
	- `D`: binomial (dichotomous) phenotype

#### Output
- `CHR${chr}_indv_VC_TWAS.txt`: File containing VC_TWAS results with TargetID gene information
- `CHR${chr}_VC_TWAS_log.txt`: File containing log messages	

#### Example Command
```bash
gene_anno_file="${TIGAR_dir}/ExampleData/gene_anno.txt"
PED_info="${TIGAR_dir}/ExampleData/PED_Info_SinglePheno.txt"
test_sample_ID_file="${TIGAR_dir}/ExampleData/test_sampleID.txt"
eQTL_ES_file="${TIGAR_dir}/ExampleData/eQTLweights.txt"
genofile="${TIGAR_dir}/ExampleData/example.vcf.gz"

# continuous phenotype
PED="${TIGAR_dir}/ExampleData/example_PED.ped"
${TIGAR_dir}/VC_TWAS.sh \
--gene_anno ${gene_anno_file} \
--PED ${PED} \
--PED_info ${PED_info} \
--test_sampleID ${test_sample_ID_file} \
--chr 1 \
--weight ${eQTL_ES_file} \
--genofile ${genofile} \
--genofile_type vcf \
--format GT \
--weight_threshold 0.0001 \
--phenotype_type C \
--thread 2 \
--out_dir ${out_dir} \
--TIGAR_dir ${TIGAR_dir}

# dichotomous phenotype
PED="${TIGAR_dir}/ExampleData/example_PED_binary.ped"
${TIGAR_dir}/VC_TWAS.sh \
--gene_anno ${gene_anno_file} \
--PED ${PED} \
--PED_info ${PED_info} \
--test_sampleID ${test_sample_ID_file} \
--chr 1 \
--weight ${eQTL_ES_file} \
--genofile ${genofile} \
--genofile_type vcf \
--format GT \
--weight_threshold 0.0001 \
--phenotype_type D \
--thread 2 \
--out_dir ${out_dir} \
--TIGAR_dir ${TIGAR_dir}
```

#### VC-TWAS using *summary-level* GWAS data

##### Additional arguments
- `--GWAS_result`: Path to GWAS result file 
- `--sample_size`: Sample size of summary-level GWAS data
- `--LD`: Path to Reference LD genotype covariance file

#### Output
- `CHR${chr}_sum_VC_TWAS.txt`: File containing VC_TWAS results with TargetID gene information
- `CHR${chr}_sum_VC_TWAS_log.txt`: File containing log messages	

#### Example Command
```bash
gene_anno_file="${TIGAR_dir}/ExampleData/gene_anno.txt"
eQTL_ES_file="${TIGAR_dir}/ExampleData/eQTLweights.txt.gz"
ld_file="${TIGAR_dir}/ExampleData/eQTLweights.txt.gz"
gwas_file="${TIGAR_dir}/ExampleData/sample_GWAS_Result.txt.gz"
${TIGAR_dir}/VC_TWAS_summary.sh \
--gene_anno ${gene_anno_file} \
--GWAS_result ${gwas_file} \
--Weight ${eQTL_ES_file} \
--sample_size ${sample_size} \
--weight_threshold 0.0001 \
--LD ${ld_file} \
--chr 1 \
--window 1000000 \
--thread 1 \
--out_dir ${out_dir}
```


### 6. BGW-TWAS 
Association study using cis- and trans-eQTL effect sizes estimate by [BGW-TWAS](https://github.com/yanglab-emory/BGW-TWAS).

#### Arguments for both VC-TWAS types
- `--gene_list`: Gene list file to specify the names of genes
- `--weight_prefix`: Start of path to SNP weight (eQTL effect size) file generated by BGTWAS. This string is the path (including directory) before the gene name.
- `--weight_suffix`: End of the path to SNP weight (eQTL effect size) file generated by BGTWAS. This string is the path (including extension) after the gene name.
- `--LD_prefix`: Start of path to LD file generated by TIGAR. This string is the path (including directory) before the chromosome number.
- `--LD_suffix`: End of the path to LD file generated by TIGAR. This string is the path (including extension) after the chromosome number.
- `--Zscore_prefix`: Start of path to Zscore summary statistics file. This string is the path (including directory) before the chromosome number.
- `--Zscore_suffix`: End of the path to Zscore summary statistics file. This string is the path (including extension) after the chromosome number.
- `--weight_threshold`: Threshold for estimated cis-eQTL effect sizes, filter SNPs with absolute cis-eQTL effect sizes smaller than threshold. (default: `0`; Note that BGWTWAS already filters weights.)
- `--thread`: Number of simultaneous *processes* to use for parallel computation (default: `1`)
- `--out_dir`: Output directory (will be created if it does not exist)
- `--out_prefix`: Start of log and output file names. (default: `BGWTWAS`)
- `--out_twas_file`: Output file name. (default: `${out_prefix}_asso.txt`, `BGWTWAS_asso.txt`)
- `--out_log_file`: Output log file. (default: `${out_prefix}_log.txt`, `BGWTWAS_log.txt`)
- `--TIGAR_dir`: Specify the directory of **TIGAR** source code

#### Example Command
```bash
out_dir=${TIGAR_dir}/ExampleData/

gene_list=${TIGAR_dir}/ExampleData/gene_list.txt

weight_prefix='${TIGAR_dir}/ExampleData/BGWTWASweights/'
weight_suffix='_BGW_eQTL_weights.txt'

LD_prefix='${TIGAR_dir}/ExampleData/CHR'
LD_suffix='_reference_cov.txt.gz'

Zscore_prefix='${TIGAR_dir}/ExampleData/CHR'
Zscore_suffix='_Zscore.txt.gz'

python ${TIGAR_dir}/TWAS/BGWTWAS_Asso_Study.py \
	--gene_list ${gene_list} \
	--weight_prefix ${weight_prefix} \
	--weight_suffix ${weight_suffix} \
	--Zscore_prefix ${Zscore_prefix} \
	--Zscore_suffix ${Zscore_suffix} \
	--LD_prefix ${LD_prefix} \
	--LD_suffix ${LD_suffix} \
	--thread 1 \
	--TIGAR_dir ${TIGAR_dir} \
	--out_dir ${out_dir}
```


## Updates
- added `--missing_rate` option (for excluding SNPs with many missing values) to model training and prediction scripts(default: `0.2`)
- weight files now automatically sorted and bgzipped/tabixed during training step
- Removed 'dfply' Python package dependency
- Added VC-TWAS method
- Added `SPrediXcan` test statistic calculation to TWAS with summary-level GWAS data (in addition to existing `FUSION` test statistic); added option (`--test_stat`) to select which test statistic to use (default: `both`)
- Added `--weight-threshold` option (for excluding SNPs with small effect-sizes) to summary-level TWAS (default: `0`)
- Added `--sampleID` argument to `TIGAR_LD.sh`
- Reduced size of output LD files
- Improved memory usage of Python scripts (especially for model-training and TWAS with summary-level GWAS data)
- Improved speed (especially for model-training)
- Improved error handling for parallelized functions
- Added/improved log output in all Python scripts
- Added time elapsed calculation for logging
- Removed intermediate `call_DPR.sh` script
- Added script to do summary-level association test with BGWTWAS weights.

## Reference
- [PrediXcan](https://github.com/hakyimlab/PrediXcan)  
- [DPR](https://github.com/biostatpzeng/DPR)
