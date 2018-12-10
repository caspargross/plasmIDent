#!/usr/bin/env Rscript
######################
# Create CIRCOS data #
# Calculate GC cont. #
######################

require(seqinr)

# Calculates the GC content in a sliding window for a character vector
GCwindow <- function(seq, width, circular=T, bases=c("c", "g")) {
  n <- length(seq)
  gc <- double(n)
  w <- floor(width/2)
  j <- 1
  sumgc <- 0
  for (i in 1:length(seq)){
    
    if (seq[i] %in% c("c", "g")) sumgc <- sumgc + 1
    
    if ((i %% width) == 0 ) {
      gc[j] <- sumgc/width
      sumgc <- 0
      j <- j + 1
    }  
  }
  return(gc)
}



# Write a vector with numerical values as Circos data output
# Width parameter aggregates the data (mean value) to reduce the number of datapoints
aggregateCircos <- function(v, width, chr = "C") {
  
  n <- floor(length(v) / width)
  out <- character(n)
  
  for (i in 1:n){
    end <- min(i * width, length(v))
    start <- max(0, end - width)
    val <- v[i]
    out[i] <- paste(chr, format(start, scientific=F), format(end, scientific=F), val, sep="\t") 
  }
  return(out)
}

# Get command line parameters
args = commandArgs(trailingOnly=TRUE)

# Load Plasmid sequences
fasta <- read.fasta(args[1])
#fasta <- read.fasta("data/testAssembly.fasta")

for (i in 1:length(fasta)) {
  print(paste("Reading contig", getName(fasta[[i]])))
  
  # Calculate GC content for Window sizes 50 and 1000
  GC50 <- GCwindow(getSequence(fasta)[[i]], 50)
  GC1000 <- GCwindow(getSequence(fasta)[[i]], 1000)
  
  # Write values to file
  vals1000 <- aggregateCircos(GC1000, 1000, getName(fasta[[i]]))
  vals50 <- aggregateCircos(GC50, 50, getName(fasta[[i]]))
  
  lapply(vals50, write, "gc50.txt", append=TRUE)
  lapply(vals1000, write, "gc1000.txt", append=TRUE)
  
}


