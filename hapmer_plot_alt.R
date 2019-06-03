#!/usr/bin/env Rscript

require("ggplot2")
require("argparser")

parser <- arg_parser("Make bubble plots")
parser <- add_argument(parser, "--file", type="character", help="hapmer count file", default="NULL")
parser <- add_argument(parser, "--output", type="character", help="file name of plot", default="NULL")
parser <- add_argument(parser, "--xlab", type="character", help="xlab", default="NULL")
parser <- add_argument(parser, "--ylab", type="character", help="ylab", default="NULL")
parser <- add_argument(parser, "--ydim", type="float", help="width of plot",  default=7)
parser <- add_argument(parser, "--xdim", type="float", help="height of plot", default=8)
args   <- parse_args(parser)
 
dat<-read.table(args$file, header=TRUE)
plot<-ggplot(data=dat, aes(x=dat[,3], y=dat[,4], size=dat[,5], colour=dat[,1]))+geom_point()
plot<-plot + xlab(names(dat)[3])
plot<-plot + ylab(names(dat)[4])
plot<-plot + scale_size_continuous(name="K-mers")
plot<-plot + scale_colour_brewer(name="ctgs", palette="Dark2")


if(args$xlab != "NULL"){
	plot<-plot + xlab(args$xlab)
}
if(args$ylab != "NULL"){
	plot<-plot + ylab(args$ylab)
}


plot<-plot + theme_bw(12)
ggsave(args$output, width=args$xdim, height=args$ydim)
