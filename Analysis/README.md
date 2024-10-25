# TIGAR-Analysis
This directory includes helpful scripts for analyzing TIGAR output.

These scripts require additional packages not used by the main TIGAR tool and are primarily written in R. Requirements are plot-specific with the exception of the `ggplot2` package, which is used in all scripts.

<!-- primarily makes use of R visualization packages, particularly ggplot2. -->

---

- [Software Setup](#software-setup)
- [Manhattan Plots](#manhattan-plots)



---


## Software Setup

### 1. Installing R packages

#### within R
```R
install.packages("ggplot2")
```

#### from bash
```bash
Rscript -e 'install.packages("ggplot2")'
```

---

## Manhattan Plots

Additional **required** R packages: `ggrepel`

```R
install.packages('ggrepel')
```

The `ggrepel` package is used to add labels to specified significant genes. If you do not wish to add those labels, the entire`geom_label_repel(...) +` command can be commented out. 


The dataframe passed to the function `manhattan_plot()` requires the following columns: `CHROM`, `POS`, `Pvalue`, and `label_text`.

<!-- The dataframe should also contain any other columns you plan to use for additional column-dependent plot functions (ex. If you later want to facet by the factor column `cohort`, the passed dataframe must have that column).


The `label_text` specifies the string to use in a label for the gene in that row. Only include a string in the `label_text` column if you want that label to show in the plot; use empty string (i.e. '', "") for all other rows. Labeled genes are assumed to be significant. 

## modify dataframe to be in correct format
	## required columns: 'CHROM', 'POS', 'Pvalue', 'label_text'
	## if you plan to facet, the dataframe passed to manhattan_plot() must have that column
	## you can include other columns
	## only include a string in the 'label_text' column if you want that label to show in the plot; use empty string (i.e. '', "") for all other rows -->


TIGAR TWAS output does not have a column called `POS`; users may choose to rename the `GeneStart` or `GeneEnd` columns or may create avother POS column (e.g. center point between `GeneStart` and `GeneEnd`) for the position argument

TIGAR TWAS output does not have a column called `Pvalue` (case-sensitive); users may rename the `PVALUE` column from individual-level TWAS output or rename one of the `FUSION_PVAL`, `SPred_PVAL`, or `PVALUE` columns (depending on output specification) from summary-level TWAS output

TIGAR TWAS output does not have a column called `label_text`; this column should only contain labels for genes you want labeled; all other rows should have the `label_text` column equal to a blank string (ie '' or ""). The color of the points for labeled genes will be red. Currently the code assumes labeled genes are significant; non-significant genes with a non-empty string in the `label_text` column will NOT be labeled.

Additional columns may be included for later use with other ggplot features. For example, a user might include a column called "TWAStype" with levels 'DPR' and 'EN'; later to create separate facets of the manhattan plot for the two different TWAS types:
```R
plot1 <- manhattan_plot(plot_data)
plot1 + facet_grid(TWAStype ~ .)
```


---


## QQ Plots

Functions to create QQ-plots of P-values.

The `qq_plot()` function does not require additional packages. The `qq_plot_pruned()` function requires the `fastqq` package.


The `qq_plot(data, obs_col)` function requires a dataframe with a numeric column of interest specified by `obs_col`. The column specified by `obs_col` will be used as the observed values. The default value of `obs_col` is `'Pvalue'`. 


If plotting is slow due to the number of points being plotted, the `qq_plot_pruned()` function can be used. It uses the the `drop_dense()` function from the `fastqq` package to prune non-improtant values. It has the same `data` and `obs_col` arguments as the `qq_plot()` function and an additional `N_hard` argument for the desired upper bound to the number of points. The default value of `N_hard` is `10000`, the same default used by the source package.


<!-- Recommended (not required) R package: `fastqq` (`drop_dense()` function, may not be required) -->
<!-- 
library('fastqq')

`fastqq::drop_dense()` -->

<!-- If plotting is slow due to the number of points being plotted, the `drop_dense()` function from the `fastqq` package can be used to prune non-improtant values. The `qq_plot_pruned()` function  -->
<!-- 
```R
# N_hard=10000 is the default desired upper bound to the number of points
pruned_data <- fastqq::drop_dense(observ, expect, N_hard=10000)
```
 -->








