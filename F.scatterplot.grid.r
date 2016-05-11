##########################
### Scatterplot.grid.R ###
##########################

# Either copy and paste entire script into R, or run script from R:
# $ R
# > source("scatterplot.grid.R")

# Load pacakages
>>>>>>> External Changes
library(reshape)
library(ggplot2)
require(gridExtra)

# Specify output path

pdf(file="/Users/Nick/Desktop/grid_FPKM_scatter.pdf", width=10, height=6)

merged<-read.table("/Users/Nick/Documents/Bioinformatics/Data/Merged_reps/Plots/merged.txt",header=FALSE, row.names=1)

sse<-read.table("/Users/Nick/Documents/Bioinformatics/Data/Six1+Six1-Eya1/Plots/merged.txt",header=FALSE, row.names=1)

ese<-read.table("/Users/Nick/Documents/Bioinformatics/Data/Eya1+Six1-Eya1/Plots/merged.txt",header=FALSE, row.names=1)

# Six1-Eya1 data
con<-log2((merged[,9])+1)
exp<-log2((merged[,13])+1)
sig<-merged[,17]
dir<-merged[,18]

# Six1 data
sse_con<-log2((sse[,9])+1)
sse_exp<-log2((sse[,13])+1)
sse_sig<-sse[,17]
sse_dir<-sse[,18]

# Eya1 data
ese_con<-log2((ese[,9])+1)
ese_exp<-log2((ese[,13])+1)
ese_sig<-ese[,17]
ese_dir<-ese[,18]

# Set plot variables for each condition 
ese_plot<-ggplot(ese, aes(ese_con, ese_exp , colour=ese_dir, size=ese_sig) ) +
	xlim(0, 15) + 
	ylim(0, 15) +
    xlab("Control log(2) (1+FPKM)") +
    ylab("Eya1+Six1-Eya1 overexpression log(2) (1+FPKM)") +
	# Shades for density of points 
  	geom_point(alpha = 0.6) +
	# point size
  	scale_size_manual(values=c(I(1), I(1.8))) +
  	scale_colour_manual(values=c('no'="grey12", 'down'="red4", 'up'="green4"))

sse_plot<-ggplot(sse, aes(sse_con, sse_exp , colour=sse_dir, size=sse_sig) ) +
	xlim(0, 15) + 
	ylim(0, 15) +
    xlab("Control log(2) (1+FPKM)") +
    ylab("Six1+Six1-Eya1 overexpression log(2) (1+FPKM)") +
	# Shades for density of points 
  	geom_point(alpha = 0.6) +
	# point size
  	scale_size_manual(values=c(I(1), I(1.8))) +
  	scale_colour_manual(values=c('no'="grey12", 'down'="red4", 'up'="green4"))
		 
merged_plot<-ggplot(merged, aes(con, exp , colour=dir, size=sig) ) +
	xlim(0, 15) + 
	ylim(0, 15) +
    xlab("Control log(2) (1+FPKM)") +
    ylab("Six1+Eya1+Six1-Eya1 overexpression log(2) (1+FPKM)") +
	# Shades for density of points 
  	geom_point(alpha = 0.6) +
	# point size
  	scale_size_manual(values=c(I(1), I(1.8))) +
  	scale_colour_manual(values=c('no'="grey12", 'down'="red4", 'up'="green4"))

# Define legends 
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
   
mylegend<-g_legend(merged_plot)

# Arrange scatterplots in a grid
grid.arrange(arrangeGrob(merged_plot + theme(legend.position="none"),
                         sse_plot + theme(legend.position="none"),
						 ese_plot + theme(legend.position="none"),
                         nrow=1))
				
dev.off()
						   
			