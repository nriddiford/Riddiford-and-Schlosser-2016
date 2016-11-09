###################
### qPCR.plot.R ###
###################

library(ggplot2)
library(gridExtra)

# This script will take in qPCR dat arranged as follows: 
# 	Gene	FC	expt	se
# e.g.	Sox2	1.02	RNA-Seq	0

# and produce two barplots one for six1 and one for eya1
# Copy and past into an active R session, or run as: 
# R
# > source("qPCR.plot.R")

dodge <- position_dodge(width = 0.6)

# Six1

six1 <- read.table("Six1_qPCR.txt", header=TRUE)

# Manual re-order
six1$Gene <- factor(six1$Gene, c("Sox2","Hes5","Hes9","Atoh1","Ngn1","Pou4f1","Gfi1a","Tlx1","Isl2"))

s<-ggplot(data=six1, aes(x=Gene, y=FC, fill=expt, group=expt)) +

    xlab("\nGene") +
    ylab("log2 Fold Change\n") +
    
	geom_bar(colour="black", stat="identity", position=dodge, width=0.6) + 
 	
	theme(
		    legend.position="none", # remove legend
	        panel.grid.major.x = element_blank(), # remove x grids
			panel.grid.minor.x = element_blank(),
	        panel.grid.major.y = element_line( size = .25, color = "#003366"), 
			panel.grid.minor.y = element_line(size = .1, colour = "#003366"),
			# panel.background = element_rect(fill = "transparent",colour = NA), # transparent backround and plot
			plot.background = element_rect(fill = "transparent",colour = NA)
			) +
	
	geom_errorbar(aes(ymax = FC + se, ymin = FC, group=expt), position = dodge, width = 0.25) +
	scale_fill_manual(values=c("lightcoral", "deepskyblue3"))
 		

# Eya1 

eya1 <- read.table("Eya1_qPCR.txt", header=TRUE)

# Manual re-order
eya1$Gene <- factor(eya1$Gene ,c("Sox2","Hes5","Hes9","Atoh1","Ngn1","Pou4f1","Gfi1a","Tlx1"))

e<-ggplot(data=eya1, aes(x=Gene, y=FC, fill=expt, group=expt)) +

    xlab("\nGene") +
    ylab("log2 Fold Change\n") +
    
	geom_bar(colour="black", stat="identity", position=dodge, width=0.6) + 
 	
	theme(
		    legend.position="none", # remove legend
	        panel.grid.major.x = element_blank(), # remove x grids
			panel.grid.minor.x = element_blank(),
	        panel.grid.major.y = element_line( size = .25, color = "#003366"), 
			panel.grid.minor.y = element_line(size = .1, colour = "#003366"),
			panel.background = element_rect(fill = "transparent",colour = NA), # transparent backround and plot
			plot.background = element_rect(fill = "transparent",colour = NA)
			) +
	
	geom_errorbar(aes(ymax = FC + se, ymin = FC, group=expt), position = dodge, width = 0.25) +
	scale_fill_manual(values=c("lightcoral", "seagreen3"))
		
# Merge

grid.arrange(s, e, ncol=2)



