


## 


qq_plot <- function(data, obs_col='Pvalue') {

	CI <- 0.95

	observ <- -log10(sort(data[, obs_col]))
	N_obs <- length(observ)
	expect <- -log10(ppoints(N_obs))
	max_expect <- max(expect)

	# dataframe for plotting
	plot_data <- data.frame(
		observ=observ,
		expect=expect,
		eq_max_expect=0
	)
	plot_data[plot_data$expect==max_expect, 'eq_max_expect'] <- 1

	p <- ggplot(plot_data, 
			aes(x=expect, y=observ)) + 
		geom_point(size=1, alpha=0.6) +	
		geom_segment(
			data=plot_dat[plot_dat$eq_max_expect==1, ],
			aes(x=0, xend=expect, y=0, yend=expect),
			alpha=0.5,
			color='grey30',
			lineend='round') +
		theme_bw() + 
		theme(aspect.ratio=1,
			legend.position='bottom') + 
		labs(x=bquote("Expected "-"log"[10]("p-value")), 
			y=bquote("Observed "-"log"[10]("p-value")))

}


qq_plot_pruned <- function(data, obs_col='Pvalue', N_hard=10000) {

	CI <- 0.95

	observ <- -log10(sort(data[, obs_col]))
	N_obs <- length(observ)
	expect <- -log10(ppoints(N_obs))

	# prune dense points
	pruned_data <- fastqq::drop_dense(observ, expect, N_hard)
	colnames(pruned_data) <- c('observ', 'expect')

	# new number of observations
	N_obs <- nrow(pruned_data)
	max_expect <- max(pruned_data$expect)

	# dataframe for plotting
	plot_data <- data.frame(
		observ=pruned_data$observ,
		expect=pruned_data$expect,
		eq_max_expect=0
	)
	plot_data[plot_data$expect==max_expect, 'eq_max_expect'] <- 1

	p <- ggplot(plot_data, 
			aes(x=expect, y=observ)) + 
		geom_point(size=1, alpha=0.6) +	
		geom_segment(
			data=plot_dat[plot_dat$eq_max_expect==1, ],
			aes(x=0, xend=expect, y=0, yend=expect),
			alpha=0.5,
			color='grey30',
			lineend='round') +
		theme_bw() + 
		theme(aspect.ratio=1,
			legend.position='bottom') + 
		labs(x=bquote("Expected "-"log"[10]("p-value")), 
			y=bquote("Observed "-"log"[10]("p-value")))

}



CI <- 0.95

observ <- -log10(sort(data[, 'pvalue']))
N_obs <- length(observ)
expect <- -log10(ppoints(N_obs))

# prune dense points
pruned_data <- fastqq::drop_dense(observ, expect)
colnames(pruned_data) <- c('observ', 'expect')

# new number of observations
N_obs <- nrow(pruned_data)
max_expect <- max(pruned_data$expect)

# dataframe for plotting
plot_data <- data.frame(
	observ=pruned_data$observ,
	expect=pruned_data$expect,
	clower=-log10(qbeta(p=(1 - CI) / 2,
		shape1=seq(N_obs),
		shape2=rev(seq(N_obs)))),
	cupper=-log10(qbeta(
		p=(1 + CI) / 2,
		shape1=seq(N_obs),
		shape2=rev(seq(N_obs)))),
	eq_max_expect=0
)
plot_data[plot_data$expect==max_expect, 'eq_max_expect'] <- 1




p1 <- ggplot(plot_data, 
		aes(x=expect, y=observ)) + 
	geom_point(size=1, alpha=0.6) +	
	geom_segment(
		data=plot_dat[plot_dat$eq_max_expect==1, ],
		aes(x=0, xend=expect, y=0, yend=expect),
		alpha=0.5,
		color='grey30',
		lineend='round') +
	theme_bw() + 
	theme(aspect.ratio=1,
		legend.position='bottom') + 
	labs(x=bquote("Expected "-"log"[10]("p-value")), 
		y=bquote("Observed "-"log"[10]("p-value")))



