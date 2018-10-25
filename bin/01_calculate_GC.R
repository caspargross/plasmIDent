#!/usr/bin/env Rscript

######################
# Create CIRCOS data #
# Calculate GC cont. #
######################

require(seqinr)
source("00_functions.R")

# Get command line parameters
args = commandArgs(trailingOnly=TRUE)

# Load Plasmid sequences
contig <- read.fasta(args[1])

# Calculate GC content for Window sizes 50 and 2000
GC50 <- GCwindow(getSequence(plC)[[1]], 50)
GC1000 <- GCwindow(getSequence(plC)[[1]], 1000)

# Write values to file
vals1000 <- aggregateCircos(GC1000, 1000)
vals50 <- aggregateCircos(GC50, 50)

lapply(vals50, write, "gc50.txt", append=TRUE)
lapply(vals1000, write, "gc1000.txt", append=TRUE)
