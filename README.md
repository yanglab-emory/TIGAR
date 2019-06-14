# TIGAR
"TIGAR" standing for Transcriptome-Intergrated Genetic Association Resource, which is developed using Python and BASH scripts.TIGAR can fit both Elastic-Net and nonparametric Beyesian model(Dirichlet Process Regression, i.e. DPR), impute transcriptomic data, and conduct genetic association studies using both individual-level and summary-level GWAS data for univariate and multivariate phenotypes.

## Software
- python 3.5 
- TABIX

## Usage Example
- Model Train and Prediction
```
./TIGAR_Model_Train.sh --model DPR \
--Gene_Exp ${Gene_Exp_path} --train_sample ${train_sample_path} \
--chr 1 --train_dir ${train_dir} \
--geno_train vcf --FT DS \
--pred y \
--geno_test dosages --FP DS \
--out ${out_prefix}
```
- TWAS
```
./TIGAR_TWAS.sh --asso 1 \
--Gene_EXP ${Gene_Exp_path} --PED ${PED} --Asso_Info ${asso_Info} \
--out ${out_prefix}
```

## Reference
- Elastic Net: https://github.com/hakyimlab/PrediXcan  
- DPR: https://github.com/biostatpzeng/DPR
