#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)

TIGAR_dir <- './TIGAR/'
out_dir <- paste0(TIGAR_dir, 'Analysis/')

source(paste0(TIGAR_dir, 'Analysis/manhattan_plot.R'))


# data is output from all chromosomes from two different datasets 
plot_data_path <- paste0(TIGAR_dir, 'ExampleData/manplot_TWAS_data.txt')
plot_data <- read.table(plot_data_path, header=TRUE, sep='\t')

## modify dataframe to be in correct format
	## required columns: 'CHROM', 'POS', 'Pvalue', 'label_text'
	## if you plan to facet, the dataframe passed to manhattan_plot() must have that column
	## you can include other columns
	## only include a string in the 'label_text' column if you want that label to show in the plot; use empty string (i.e. '', "") for all other rows

plot_data$POS <- plot_data$GeneStart

# remove excess whitespace from gene names (helpful for making labels/matching genes by name)
plot_data$GeneName <- trimws(as.character(plot_data$GeneName))

# recommend calculating Pvalue from the Zscore (if original output may have rounding issues, be read-in as a character, etc.)
plot_data$Pvalue <- exp(pchisq(plot_data[['SPred_Z']]^2, 1, lower.tail=FALSE, log.p=TRUE))
# # can also use the SPred pvalue for plotting; rename 'SPred_PVAL'
# colnames(plot_data)[which(colnames(plot_data) == 'SPred_PVAL')] <- 'Pvalue'

# if your TargetID may be an ensemble ID with version number but this may not be desirable (e.g. for matching between datasets where one does not include the version, etc.); this is a convenient way to remove the version number
plot_data$TargetID_w_version <- plot_data$TargetID # ex: ENSG00000223972.4
plot_data$TargetID <- do.call(rbind, strsplit(plot_data$TargetID, '.', fixed=TRUE))[, 1] # ex: ENSG00000223972

# NOTE: labels may not show up in plot viewing window-check the output file

# label these genes are labeled
id_genes <- c('MMRN1','MAPK8IP1P2','GPNMB','MAPT','PLEKHN1','PRR26')

plot_data$label_text <- plot_data$GeneName
plot_data[with(plot_data, !(GeneName %in% id_genes)), 'label_text'] <- ''

# only label significant genes
p1 <- manhattan_plot(plot_data); p1
p1 + facet_grid(dataset ~ .)

ggsave(paste0(out_dir, 'manplot_1-1.pdf'), p1) # ggsave will automatically format the output based on file extension
ggsave(paste0(out_dir, 'manplot_1-2.png'), p1 + facet_grid(dataset ~ .)) # example: facet by 'dataset' column
ggsave(paste0(out_dir, 'manplot_1-3.jpg'), p1, width=4.5, height=6.5) # example specify height

# label all nonsig and significant genes
p2 <- manhattan_plot(plot_data, label_nonsig=TRUE); p2
p2 + facet_grid(dataset ~ .)

ggsave(paste0(out_dir, 'manplot_2-1.pdf'), p2) # ggsave will automatically format the output based on file extension
ggsave(paste0(out_dir, 'manplot_2-2.png'), p2 + facet_grid(dataset ~ .)) # example: facet by 'dataset' column
ggsave(paste0(out_dir, 'manplot_2-3.jpg'), p2, width=4.5, height=6.5) # example specify height


## example: drop most significant gene (if it is much, much more significant than other genes)
# plot_data <- plot_data[-with(plot_data, which(Pvalue == min(Pvalue))), ]

