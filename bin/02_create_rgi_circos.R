#!/usr/bin/env Rscript

######################
# Create CIRCOS data #
#  tracks from rgi   #
######################

require(data.table)

args <- commandArgs(trailingOnly=TRUE)
file <- args[1]
rgi <- fread(file)


# Export RGI
fwrite(rgi[,.(Contig, Start, Stop, Best_Hit_ARO)], sep="\t", file ="rgi.txt", col.names = F)
fwrite(rgi[,.(Contig, Start, Stop)], sep="\t", file ="rgi_span.txt", col.names = F)
