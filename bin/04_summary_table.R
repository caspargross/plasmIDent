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

#fasta <- "/home/caspar/Documents/Uni/MasterThesis/analysis/VRE/assembly_lr/01_VRE/04_assembled_genomes/01_VRE_unicycler_final_assembly.fasta"
#rgi <- "/home/caspar/Documents/Uni/MasterThesis/analysis/VRE/plasmid_finder/work/85/b5f40d5b64e074f45f06d3f22dc6a3/VRE01_rgi.txt"
#cov <- "/home/caspar/Documents/Uni/MasterThesis/analysis/VRE/plasmid_finder/work/6c/887dc4e6824f1c1814d91254ea2bfa/cov.txt"
#gc <- "/home/caspar/Documents/Uni/MasterThesis/analysis/VRE/plasmid_finder/work/c9/a747c7c63fde5e23d65817c247f63a/gc1000.txt"

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

#setkey(dt_ar, Contig)

dt <- data.table(id = getName(seq), 
                 length = sapply(getSequence(seq), length),
                 cov = val_cov,
                 gc = val_gc)

#setkey(dt, id)
dt[,contig := getID(id, 0),]

dt <- merge(dt, dt_ar, by.x = 'contig', by.y = 'Contig', all.x = T)
# Export RGI
fwrite(dt, sep="\t", file ="contig_summary.txt")
