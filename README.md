## TIGAR
"TIGAR" standing for Transcriptome-Intergrated Genetic Association Resource, which is developed using Python and BASH scripts. TIGAR can fit both Elastic-Net and nonparametric Beyesian model (Dirichlet Process Regression, i.e. DPR) for gene expression imputation, impute genetically regulated gene expression (GReX) from genotype data, and conduct transcriptome-wide association studies (TWAS) using both individual-level and summary-level GWAS data for univariate and multivariate phenotypes.

### Software

- Please add the executable file `./TIGAR/Model_Train/DPR` to your linux `${PATH}` directory. Assuming `~/bin/` is a directory added to your `${PATH}` environmental variable, you can accomodate the following example command

```
cp ./TIGAR/Model_Train/DPR ~/bin/
```

- BGZIP, TABIX, Python 3.5 and the following python libraries are required for running TIGAR
1. BGZIP: http://www.htslib.org/doc/bgzip.html 
2. TABIX: http://www.htslib.org/doc/tabix.html 
3. python 3.5 
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

### Reference
- Elastic Net: https://github.com/hakyimlab/PrediXcan  
- DPR: https://github.com/biostatpzeng/DPR
