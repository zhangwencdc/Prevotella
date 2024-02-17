#################################################################
# Function: Bloom Figure
# Call: Rscript bloom.R -i abund_file -o outfile
# R packages used: optparse, reshape2, ggplot2,RColorBrewer, grDevices
# Last update: 2022-06-14, Zhang Wen

#################################################################
# install necessary libraries
p <- c("optparse","reshape2","ggplot2","RColorBrewer","grid","scales","vegan","agricolae","gridExtra","dplyr","ggrepel","gggenes","ggsignif","pheatmap","reshape2","picante","corrplot","ape","multcomp","patchwork","factoextra")
usePackage <- function(p) {
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep=TRUE, repos="http://mirrors.opencas.cn/cran/")
  suppressWarnings(suppressMessages(invisible(require(p, character.only=TRUE))))
}
invisible(lapply(p, usePackage))
## clean R environment
rm(list = ls())
setwd('./')
## parsing arguments
args <- commandArgs(trailingOnly=TRUE)
# make option list and parse command line
option_list <- list(
  make_option(c("-i", "--abund_file"), type="character", help="Input feature table with relative abundance (*.Abd) [Required]"),

  make_option(c("-o", "--out_dir"), type="character", default='MetaFigure', help="Output directory [default %default]"),
  make_option(c("-p", "--prefix"), type="character", default='Out', help="Output file prefix [default %default]"),
  make_option(c("-t", "--threshold"), type="double", default=0.01, help="Average value threshold [Optional, default %default]"),
  make_option(c("--cutoff_positivesample", "-c"), type="numeric", default=0, help="cutoff for positive sample number  [Optional, default %default]"),
  make_option(c("--cutoff_readpercentage", "-r"), type="numeric", default=0, help=" cutoff for read number/percentage  [Optional, default %default]"),
  make_option(c("--overturn", "-f"), type="character", default=F, help=" Each row for a sample (Default); -turn T Each column for a sample  [Optional, default %default]")
)
opts <- parse_args(OptionParser(option_list=option_list), args=args)
# paramenter checking
if(is.null(opts$abund_file)) stop('Please input a feature table (*.Abd)')
# load data
matrixfile <- opts$abund_file
mapfile <- opts$meta_data
ave_t <- opts$threshold
outpath <- opts$out_dir
cut_i<-opts$cutoff_positivesample
cut_j<-opts$cutoff_readpercentage
turn<-opts$overturn


			library("ggplot2") # load related packages
			library("grid")
			library("scales")
			library("vegan")
			library("agricolae")
			library("gridExtra")
			library("dplyr")
			library("ggrepel")
			library("gggenes")
			library("ggsignif")
			library("pheatmap")
			library(reshape2)
			library(picante)
			library(corrplot)
			library(ape)
			library(multcomp)
			library(patchwork)
			library(factoextra)
otu <- read.table(matrixfile,header = T, row.names = 1,sep="\t")
if(is.null(opts$overturn)==F) {otu<-otu}else{otu<-t(otu)}

# 按日期从早到晚对数据进行排序
otu <- otu[order(otu$Time), ]

p<-ggplot(data=otu,aes(x=factor(Time,levels=unique(Time)),y=Percentage,group=People))+geom_line(aes(color=People),linewidth=1.5)+geom_point(aes(color=People),size=4)
ggsave(paste(outpath,"_line.pdf", sep=""), p, width = 36, height = 10)

p1<-ggplot(data=otu,aes(x=log10(Percentage),fill=People))+geom_density(aes(color=People),alpha=0.3)
ggsave(paste(outpath,"_abundance.pdf", sep=""), p1, width = 15, height = 10)