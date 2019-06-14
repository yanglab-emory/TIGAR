#### 1. Gene Expression File
| CHROM | GeneStart | GeneEnd | TargetID/GeneID | GeneName | sample1 |
|:-----:|:---------:|:-------:|:---------------:|:--------:|:-------:|
|   1   |    100    |   200   |     ENSG0000    |     X    |   0.2   |

#### 2. Genotype File
- vcf  

| CHROM | POS |  ID | REF | ALT | QUAL | FILTER | INFO | FORMAT |  sample1 |
|:-----:|:---:|:---:|:---:|:---:|:----:|:------:|:----:|:------:|:--------:|
|   1   | 100 | rs1 |  C  |  T  |   .  |  PASS  |   .  |  GT:DS | 0/0:0.01 |

TIGAR accept GT format like 0/0 or 0|0.

- dosages

| CHROM | POS |  ID | REF | ALT | sample1 |
|:-----:|:---:|:---:|:---:|:---:|:-------:|
|   1   | 100 | rs1 |  C  |  T  |   0.01  |

#### 3. PED

| FAM_ID | IND_ID | FAT_ID | MOT_ID | COV1 | COV2 | PHENO |
|:------:|:------:|:------:|:------:|:----:|:----:|:-----:|
|   11A  |   11A  |    X   |    X   |   1  |  0.2 |  0.3  |
