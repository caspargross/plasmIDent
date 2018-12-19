#!/usr/bin/env Rscript

######################
# Create CIRCOS data #
#  tracks from rgi   #
######################

require(data.table)
require(seqinr)
require(ggplot2)

args <- commandArgs(trailingOnly=TRUE)
fasta <- args[1]
rgi <- args[2]
cov <- args[3]
gc <- args[4]
padding <- as.integer(args[5])

seq <- read.fasta(fasta)
dt_rgi <- fread(rgi)

dt_cov <- fread(cov)
setnames(dt_cov, c("contig_name", "start", "stop", "cov"))
dt_cov[,contig_name := as.character(contig_name)]

dt_gc <- fread(gc)
setnames(dt_gc, c("contig_name", "start", "stop", "gc"))
dt_gc[,contig_name := as.character(contig_name)]

val_cov <- dt_cov[,median(cov), by = contig_name]$V1
val_gc <- dt_gc[,mean(gc), by = contig_name]$V1

getID <- function(s, x=1) {
  split <- transpose(strsplit(s, "_"))
  n <- length(split)
  return(split[[n-x]])
}

dt_rgi[, Contig := getID(Contig),]

dt_ar <- dt_rgi[,.(ar_genes=toString(Best_Hit_ARO)), by=Contig]

dt <- data.table(id = getName(seq), 
                 length = sapply(getSequence(seq), function (x) length(x) + padding),
                 cov = val_cov,
                 gc = val_gc)

dt[,contig := getID(id, 0),]

dt <- merge(dt, dt_ar, by.x = 'contig', by.y = 'Contig', all.x = T)
fwrite(dt, sep="\t", file ="contig_summary.txt")
