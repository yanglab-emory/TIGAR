#!/usr/bin/env Rscript
options(stringsAsFactors=FALSE)

## load packages
library(ggplot2)
library(ggrepel)

## modify dataframe to be in correct format
	## required columns: 'CHROM', 'POS', 'Pvalue', 'label_text'
	## if you plan to facet, the dataframe passed to manhattan_plot() must have that column
	## you can include other columns
	## only include a string in the 'label_text' column if you want that label to show in the plot; use empty string (i.e. '', "") for all other rows


# TIGAR TWAS output does not have a column called "POS"; users may choose to rename the GeneStart or GeneEnd columns or may create avother POS column (e.g. center point between GeneStart and GeneEnd) for the position argument
# TIGAR TWAS output does not have a column called "Pvalue" (case-sensitive); users may rename the PVALUE column from individual-level TWAS output or rename one of the FUSION_PVAL, SPred_PVAL, or PVALUE columns (depending on output specification) from summary-level TWAS output
# TIGAR TWAS output does not have a column called "label_text"; this column should only contain labels for genes you want labeled; all other rows should have the label_text column equal to a blank string (ie '' or ""). The color of the points for labeled genes will be red. The code assumes labeled genes are significant. 

# Additional columns may be included for later use with other ggplot features. For example, a user might include a 'TWAStype' column with levels 'DPR' and 'EN'; later to create separate facets of the manhattan plot for the two different TWAS types.
#		p1 <- manhattan_plot(mplot_data)
#		p1 + facet_grid(dataset ~ .)



## function to do the manhattan plot
manhattan_plot <- function(
	## data (columns: 'CHROM','POS','Pvalue','label_text')
		data, 
	## genes:
		nonsig_colors=c('#44B5AD','#3A948E','#36807A','#2f615d'), # color for non-sig genes
		nonsig_label_color= '#C5E063', # color for non-sig. genes with labels
		label_nonsig=FALSE, # by default never color/label non-significant genes
		sig_color='#FD923F', # color for sig. genes without labels
		sig_label_color='#D92B26', # color for sig. genes with labels
		point_alpha=0.9, # transparency for genes
	## significance level:
		sig_level=2.5e-6, # significance level value
		sig_level_line_col='black', # line color
		sig_linetype='dashed', # linetype
		sig_line_size=1, # line size
	## plot theme, other plot variables:
		chr_vec=1:22, # chromosomes to plot
		chr_gap=500, # gap between chromosomes in x-axis
		theme=theme_bw(), # ggplot2 theme (can use custom theme)
		plot_bg_col=NULL, # background color if different from theme
		# panel_border=NULL,#element_blank(), # ggplot panel border (default: blank)
		panel_border=element_rect(fill=NA, colour='#333333', size=0.175),
		strip_background=element_rect(colour='black', size=0.175),
		text_size=10, # text size
	## point labelling:
		geom_label_size=2, # label text size
		label_fill='white', # label background color
		label_col='black', # label border color
		label_seg_col='black', # color of line from label to point
		min_segment_length=0.01, # minimum length of line from label to point
		segment_size=0.2, # line from label to point
		label_force=2, # force of repulsion between overlapping text labels
		point_padding=1e-06, # padding around genes
		seed=NA, # optional seed for generating label positions
		max_iter=15000 # number of iterations to use to generate label positions
	){

	# setup dataframe for plotting; get plot positions from chromosome positions
	mplot_data <- NULL # dataframe
	endPos <- 0  # place on x-axis where last chromosome ended
	x_axis_chr_breaks <- NULL # chromosome label positions
	x_axis_chr_labels <- NULL # list of chromosomes to label
	for (chr in chr_vec) {
		# get data for chr
		temp <- data[data$CHROM==chr, ]
		if (nrow(temp) > 0) {
			# append chromosome to list of chromosomes to label
			x_axis_chr_labels <- c(x_axis_chr_labels, chr)
			# get unique positions for this chr
			uniq_pos <- sort(unique(temp$POS))
			uniq_pos <- setNames(order(uniq_pos), uniq_pos)
			# set POS to order value of the unique positions
			temp$POS <- uniq_pos[as.character(temp$POS)]
			# get plot positions for genes on this chr
			temp$plotPos <- (temp$POS - min(temp$POS, na.rm=TRUE) ) + endPos + 1
			# get new end position based on max position
			endPos <- max(temp$plotPos, na.rm=TRUE) + chr_gap
			# append label position
			x_axis_chr_breaks <- c(x_axis_chr_breaks, mean(temp$plotPos, na.rm=TRUE) )
			# add rows to mplot_data
			mplot_data <- rbind(mplot_data, temp)
		}
	}
	# set min, max values for axes
	min_x <- min(mplot_data$plotPos)
	max_x <- max(mplot_data$plotPos)
	max_y <- max(-log10(mplot_data$Pvalue))
	# all empty label_text column if it's not in the dataframe
	if ( !('label_text' %in% colnames(mplot_data)) ) {
		mplot_data$label_text <- ''
	}

	# plot
	# different initial plot based on whether label_nonsig is true
		if (label_nonsig) {
			p <- ggplot(data=mplot_data, 
						aes(x=plotPos, y=-log10(Pvalue), label=label_text)) + 
				# non-sig. genes without labels:
				geom_point(data=subset(mplot_data, Pvalue >= sig_level & nchar(label_text) == 0),
					aes(x=plotPos, y=-log10(Pvalue), color=factor(CHROM)),
					size=1, alpha=point_alpha) + 
				scale_color_manual(values=rep(nonsig_colors, 22)) + 
				# non-sig. genes with labels:
				geom_point(data=subset(mplot_data, Pvalue >= sig_level & nchar(label_text) > 0),
					aes(x=plotPos, y=-log10(Pvalue)), 
					color=nonsig_label_color,
					size=1.5, alpha=1) +
				# add labels
				geom_label_repel(data=subset(mplot_data, Pvalue >= sig_level & label_text!=''), 
					min.segment.length=min_segment_length,
					segment.size=segment_size,
					segment.color=label_seg_col,
					box.padding=1.1,
					size=geom_label_size, 
					alpha=1,
					# reduce labels covering up other non-sig genes/significance line
					ylim=c(0.35, 0.85) * -log10(sig_level), 
					xlim=c(min_x, max_x),
					force=label_force,
					point.padding=point_padding,
					max.iter=max_iter,
					colour=label_col,
					fill=label_fill,
					seed=seed)
		} else {
			p <- ggplot(data=mplot_data, 
						aes(x=plotPos, y=-log10(Pvalue), label=label_text)) + 
				# non-sig. genes:
				geom_point(data=subset(mplot_data, Pvalue >= sig_level),
					aes(x=plotPos, y=-log10(Pvalue), color=factor(CHROM)),
					size=1, alpha=point_alpha) + 
				scale_color_manual(values=rep(nonsig_colors, 22))
		}

	p <- p +
		# sig. genes
		geom_point(data=subset(mplot_data, Pvalue < sig_level),
			aes(x=plotPos, y=-log10(Pvalue), fill=factor(CHROM)),
			size=ifelse(subset(mplot_data, Pvalue < sig_level)$label_text=='', 1.25, 1.5),
			color=ifelse(subset(mplot_data, Pvalue < sig_level)$label_text=='', sig_color, sig_label_color),
			alpha=ifelse(subset(mplot_data, Pvalue < sig_level)$label_text=='', point_alpha, 1)) +
		# add labels
		geom_label_repel(data=subset(mplot_data, Pvalue < sig_level), 
			min.segment.length=min_segment_length,
			segment.size=segment_size,
			segment.color=label_seg_col,
			box.padding=1.1,
			size=geom_label_size, 
			alpha=1,
			ylim=c(-log10(sig_level), max_y),
			xlim=c(min_x, max_x),
			force=label_force,
			point.padding=point_padding,
			max.iter=max_iter,
			colour=label_col,
			fill=label_fill,
			seed=seed) +
		# significance level line
		geom_hline(yintercept=-log10(sig_level), 
			linetype=sig_linetype, 
			size=sig_line_size, 
			color=sig_level_line_col) +
		# remove legend
		guides(color='none', fill='none') + 
		# set axes titles
		labs(x='Chromosome', y=bquote(-"log"[10]("p-value"))) + 
		# x-axis labels, breaks
		scale_x_continuous(breaks=x_axis_chr_breaks, 
			labels=x_axis_chr_labels, 
			expand=c(0.01, 0)) + 
		# don't clip to extent of plot panel
		coord_cartesian(clip='off') +
		# pad y-axis
		scale_y_continuous(expand=c(0.05, 0)) +
		theme +
		# convenient theme options for a manhattan plot
		theme(
			text=element_text(size=text_size),
			# text=element_text(size=text_size, face='bold'), 
			# axis.title=element_text(face='bold', size=text_size),
			axis.text.x=element_text(size=text_size-1, 
				face='plain', angle=-90, 
				vjust=0.5, hjust=0),
			axis.text.y=element_text(#face='bold', 
					size=text_size-1),
			panel.grid.major.x=element_blank(),
			panel.grid.minor.x=element_blank(),
			panel.border=panel_border,
			plot.background=element_rect(fill=plot_bg_col),
			plot.tag=element_text(face='bold', size=text_size+5, family='Helvetica'),
			strip.background=strip_background
			)
}
