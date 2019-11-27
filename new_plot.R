#!/usr/bin/env Rscript

library("argparser")
library("ggplot2")
library("plyr")
library("scales")

################################################################################
draw_blobplot <- function(counts.file,
                          mer.size,
                          out.prefix,
                          ctg.sz.cutoff,
                          w,
                          h,
                          custom.x.lab,
                          custom.y.lab) {

    counts.matrix = read.delim(counts.file, stringsAsFactors = F, check.names = F)
    counts.matrix$"length" = counts.matrix$"Total" + mer.size
    if (! is.null(ctg.sz.cutoff) & !is.na(ctg.sz.cutoff)) {
        counts.matrix = counts.matrix[counts.matrix$"length"<= ctg.sz.cutoff, ]
    }

    hapA = names(counts.matrix)[3]
    hapB = names(counts.matrix)[4]
    custom.x.lab = if(is.null(custom.x.lab) | is.na(custom.x.lab)) hapA else custom.x.lab
    custom.x.lab = sprintf("%s (%d-mer)", custom.x.lab, mer.size)
    custom.y.lab = if(is.null(custom.y.lab) | is.na(custom.y.lab)) hapB else custom.y.lab
    custom.y.lab = sprintf("%s (%d-mer)", custom.y.lab, mer.size)

    max.display = min(max(counts.matrix[[hapA]]), max(counts.matrix[[hapB]]))
    max.display = ceiling(max.display / 1000) * 1000 # round to nearest 1000
    if (max(counts.matrix[[hapA]]) < max(counts.matrix[[hapB]])) {
        warning(sprintf("Not displaying part of %s due to sparse data beyond %d,
                         %.02f%% data discarded",
                        hapB, max.display,
                        100 * length(which(counts.matrix[[hapB]] > max.display)) / nrow(counts.matrix)))
    } else if (max(counts.matrix[[hapA]]) > max(counts.matrix[[hapB]])) {
        warning(sprintf("Not displaying part of %s due to sparse data beyond %d,
                         %.02f%% data discarded",
                        hapA, max.display,
                        100 * length(which(counts.matrix[[hapA]] > max.display)) / nrow(counts.matrix)))
    }

    MAGNIFY.FACTOR = 2

    ggplot(counts.matrix, aes(x=get(hapA), y=get(hapB), color=Assembly, size=length)) +
        geom_point(shape=16) +
        theme_bw() +
        scale_fill_brewer(palette = "Set1") +     # Set1 -> Red / Blue. Set2 -> Green / Orange.
        scale_color_brewer(name="Assembly", palette = "Set2") +
        scale_x_continuous(labels=comma, limits = c(0, max.display)) +
        scale_y_continuous(labels=comma, limits = c(0, max.display)) +
        scale_radius(labels=comma, range = c(1,20), name = "Contig size") +
        theme(legend.text = element_text(size=11),
              legend.position = c(0.80,0.50),  # Modify this if the legend is covering your favorite circle
              legend.background = element_rect(size=0.5, linetype="solid", colour ="black"),
              axis.title=element_text(size=14,face="bold"),
              axis.text=element_text(size=12, face="bold")) +
        guides( size = guide_legend(order = 1),
                colour = guide_legend(override.aes = list(size=5), order = 2)) +
        xlab(custom.x.lab) +
        ylab(custom.y.lab)

    # ggsave(file = paste0(out.prefix, '.png'), height = h, width = w)
    ggsave(file = paste0(out.prefix, '.pdf'),
           height = MAGNIFY.FACTOR * h, width = MAGNIFY.FACTOR * w)
}
################################################################################
parser <- arg_parser("Make bubble plots")
parser <- add_argument(parser, "--file", type="character", help="hapmer count file", default=NULL)
parser <- add_argument(parser, "--mersize", type="integer", help="mer size used for generating counts file", default=21)
parser <- add_argument(parser, "--ctgszcutoff", type="integer", help="contig size cutoff (contig longer than this won't be used", default=NULL)
parser <- add_argument(parser, "--output", type="character", help="prefix to file name of plot", default=NULL)
parser <- add_argument(parser, "--xlab", type="character", help="xlab", default=NULL)
parser <- add_argument(parser, "--ylab", type="character", help="ylab", default=NULL)
parser <- add_argument(parser, "--width", type="float", help="width of plot",  default=6.5)
parser <- add_argument(parser, "--height", type="float", help="height of plot", default=6)
args   <- parse_args(parser)

draw_blobplot(counts.file = args$"file",
              mer.size = args$"mersize",
              out.prefix = args$"output",
              ctg.sz.cutoff = args$"ctgszcutoff",
              custom.x.lab = args$"xlab",
              custom.y.lab = args$"ylab",
              w=args$"width",
              h=args$"height")
