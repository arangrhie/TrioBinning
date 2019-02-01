library(ggplot2)
library(ggsci)
library(plyr)
library(scales) 

# Draw haplotype-specific kmer plot

draw_blobplot <- function(dataset, hapA, hapB, total, x_lab=names(hapA), y_lab=names(hapB), out, w=6.5, h=6, large_w=13, large_h=12) {
  max_total=max(max(hapA), max(hapB))
  #max(hapA)
  #max(hapB)
  #min(total)
  #median(total)
  
  max_total = max_total * 1.01
  
  ggplot(dataset, aes(x=hapA, y=hapB, color=Assembly, size=total)) +
    geom_point(shape=16) + theme_bw() + scale_color_brewer(palette = "Set2") +
    scale_fill_brewer(palette = "Set1") +     # Set1 -> Red / Blue. Set2 -> Green / Orange.
    scale_x_continuous(labels=comma, limits = c(0, max_total)) +
    scale_y_continuous(labels=comma, limits = c(0, max_total)) +
    scale_size_continuous(labels=comma, range = c(1, 10), name = "Total k-mers") +
    theme(legend.text = element_text(size=11),
          legend.position = c(0.80,0.50),  # Modify this if the legend is covering your favorite circle
          legend.background = element_rect(size=0.5, linetype="solid", colour ="black"),
          axis.title=element_text(size=14,face="bold"),
          axis.text=element_text(size=12, face="bold")) +
    guides( size = guide_legend(order = 1),
            colour = guide_legend(override.aes = list(size=5), order = 2)) +
    xlab(x_lab) + ylab(y_lab)
  
  ggsave(file = paste(out, '.png', sep=""), height = h, width = w)
  ggsave(file = paste(out, '.pdf', sep=""), height = large_h, width = large_w)
  
}

# Zebrafinch example

setwd("/Users/rhiea/FALCON-Phase/zebrafinch") # Change this to the proper path

# triocanu
kmercounts=read.table("triocanu.counts", header=TRUE)
draw_blobplot(out="triocanu", kmercounts, kmercounts$R.JL230, kmercounts$G.JL39, kmercounts$Total, x_lab = "R-JL230", y_lab = "G-JL39")

# falcon-unzip
kmercounts=read.table("falconunzip.counts", header=TRUE)
draw_blobplot(out="falcon-unzip", kmercounts, kmercounts$R.JL230, kmercounts$G.JL39, kmercounts$Total, x_lab = "R-JL230", y_lab = "G-JL39")

# falcon-phase
kmercounts=read.table("falcon-phase.counts", header=TRUE)
draw_blobplot(out="falcon-phase", kmercounts, kmercounts$R.JL230, kmercounts$G.JL39, kmercounts$Total, x_lab = "R-JL230", y_lab = "G-JL39")

# falcon-phase-scaff
kmercounts=read.table("falcon-phase-scaff.counts", header=TRUE)
draw_blobplot(out="falcon-phase-scaff", kmercounts, kmercounts$R.JL230, kmercounts$G.JL39, kmercounts$Total, x_lab = "R-JL230", y_lab = "G-JL39")
