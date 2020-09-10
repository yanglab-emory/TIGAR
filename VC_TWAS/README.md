### Python scripts used for VC-TWAS

1. VC-TWAS with cis-eQTL effect sizes estimated from nonparametric Bayesian DPR method for Univariate phenotypes `VC-TWAS.py`
- additional scripts
	- `SKAT.py`
	- `qfc_checked.py`

2.VC_TWAS individual GWAS data
- SKAT.SKAT(genotype, phenotype, covariate, weight, phenotype_type)
	- Format: Phenotype, covariate should be pandas DataFrame. genotype and weights should be numpy.ndarray.
	- phenotype_type: "C" for continuous phenotype, "D" for binary phenotype.

3.VC_TWAS summary GWAS data
- without the knowledge of the sum of square phenotype, will estimate the sum of square phenotype based on the SNPs within GENE region.
	- SKAT.SKAT_summary(beta_var, beta_estimate, weight, sample_size, COV, D)
		- beta_estimate: Estimated coefficient from GWAS result for certain SNPs.
		- beta_var: Variance of estimated coefficient for certain SNPs.
		- sample_size: The sample size of GWAS result data.
		- COV: The covariance matrix
		- D: Diag matrix of COV

- withthe knowledge of the sum of square phenotype
	- SKAT.SKAT_summary_withy(y_square,beta_estimate, weight, sample_size, COV, D)
		- y_square:the sum of square phenotype

- estimate the sum of square phenotype based on selected region.
	- SKAT.SKAT_summary_gety(beta_var,beta_estimate, sample_size, COV, D)
		- NOTE: beta_var,beta_estimate,COV,D  are based on the selected region.
