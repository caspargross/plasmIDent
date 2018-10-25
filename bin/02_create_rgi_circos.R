#!usr/bin/env Rscript

######################
# Create CIRCOS data #
#  tracks from rgi   #
######################

library(seqinr)
source("00_functions.R")

args = commandArgs(trailingOnly=TRUE)
file = args[1]

rgi <- readGff(file)

# Export RGI
fwrite(rgi[,.(seqname, start, end, Name)], sep="\t", file ="rgi.txt", col.names = F)
fwrite(rgi[,.(seqname, start, end)], sep="\t", file ="rgi_span.txt", col.names = F)
