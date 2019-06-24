## TIGAR
"TIGAR" standing for Transcriptome-Intergrated Genetic Association Resource, which is developed using Python and BASH scripts. TIGAR can fit both Elastic-Net and nonparametric Beyesian model (Dirichlet Process Regression, i.e. DPR) for gene expression imputation, impute genetically regulated gene expression (GReX) from genotype data, and conduct transcriptome-wide association studies (TWAS) using both individual-level and summary-level GWAS data for univariate and multivariate phenotypes.

### Software

TABIX, Python 3.5 and the following python libraries are required for running TIGAR
1. TABIX: http://www.htslib.org/doc/tabix.html 
2. python 3.5 
   - dfply
   - io
   - subprocess
   - multiprocess

### Example Usage 
- More details are available in the TIGAR_Manual.pdf
- Train gene expression imputation model
```
./TIGAR_Model_Train.sh --model DPR \
--Gene_Exp ${Gene_Exp_path} --train_sample ${train_sample_path} \
--chr 1 --train_dir ${train_dir} \
--geno_train vcf --FT DS \
--out ${out_prefix}
```

- Predict GReX
```

```

- TWAS
```
./TIGAR_TWAS.sh --asso 1 \
--Gene_EXP ${Gene_Exp_path} --PED ${PED} --Asso_Info ${asso_Info} \
--out ${out_prefix}
```

### Input file format
The example data provided here are generated artificially.


#### 1. Gene Expression File (Gene_Exp_combination.txt)
| CHROM | GeneStart | GeneEnd | TargetID/GeneID | GeneName | sample1 |
|:-----:|:---------:|:-------:|:---------------:|:--------:|:-------:|
|   1   |    100    |   200   |     ENSG0000    |     X    |   0.2   |


#### 2. Genotype File
1) vcf file
- http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/

| CHROM | POS |  ID | REF | ALT | QUAL | FILTER | INFO | FORMAT |  sample1 |
|:-----:|:---:|:---:|:---:|:---:|:----:|:------:|:----:|:------:|:--------:|
|   1   | 100 | rs1 |  C  |  T  |   .  |  PASS  |   .  |  GT:DS | 0/0:0.01 |

2) dosages file

| CHROM | POS |  ID | REF | ALT | sample1 |
|:-----:|:---:|:---:|:---:|:---:|:-------:|
|   1   | 100 | rs1 |  C  |  T  |   0.01  |

#### 3. PED File
- http://zzz.bwh.harvard.edu/plink/data.shtml#ped

| FAM_ID | IND_ID | FAT_ID | MOT_ID | SEX | PHENO | COV |
|:------:|:------:|:------:|:------:|:---:|:-----:|:---:|
|   11A  |   11A  |    X   |    X   |  1  |  0.2  | 0.3 |

#### 4.Zscore File

| CHROM | POS | REF | ALT | Zscore |
|:-----:|:---:|:---:|:---:|:------:|
|   1   | 100 |  C  |  T  |  0.01  |

#### 5. Genome block annotation file
| CHROM | Start | End | File |
|:-----:|:---------:|:-------:|:---------------:|
|   1   |    100    | 20000   |  CHR1.vcf.gz    |

### Reference
- Elastic Net: https://github.com/hakyimlab/PrediXcan  
- DPR: https://github.com/biostatpzeng/DPR
