#### 1. Gene Expression File (Gene_Exp_combination.txt)
| CHROM | GeneStart | GeneEnd | TargetID/GeneID | GeneName | sample1 |
|:-----:|:---------:|:-------:|:---------------:|:--------:|:-------:|
|   1   |    100    |   200   |     ENSG0000    |     X    |   0.2   |


#### 2. Genotype File
1) vcf 
- http://www.internationalgenome.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40/

| CHROM | POS |  ID | REF | ALT | QUAL | FILTER | INFO | FORMAT |  sample1 |
|:-----:|:---:|:---:|:---:|:---:|:----:|:------:|:----:|:------:|:--------:|
|   1   | 100 | rs1 |  C  |  T  |   .  |  PASS  |   .  |  GT:DS | 0/0:0.01 |

2) dosages file

| CHROM | POS |  ID | REF | ALT | sample1 |
|:-----:|:---:|:---:|:---:|:---:|:-------:|
|   1   | 100 | rs1 |  C  |  T  |   0.01  |

#### 3. PED 
- http://zzz.bwh.harvard.edu/plink/data.shtml#ped

| FAM_ID | IND_ID | FAT_ID | MOT_ID | SEX  | PHENO | COV |
|:------:|:------:|:------:|:------:|:----:|:----:|:-----:|
|   11A  |   11A  |    X   |    X   |   1  |  0.2 |  0.3  |











