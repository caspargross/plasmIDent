#!/usr/bin/env Rscript
######################
# Create CIRCOS data #
# Calculate GC cont. #
######################

require(seqinr)

# Calculates the GC content in a sliding window for a character vector
GCwindow <- function(seq, width, circular=T, bases=c("c", "g")) {
  n <- floor(length(seq)/width)
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
# Also removes sequence padding
aggregateCircos <- function(v, width, chr = "C", len) {
  
  n <- length(v) - floor(padding/width)
  out <- character(n)
  
  offset = padding %% width
  j <- 1
  
  
  for (i in 1:n){
    end <- min(i * width, length(v))
    start <- max(0, end - width)
    val <- v[i]
    if ((start + width) > padding) {
      if (((end - width )< (len - (2*padding)))) {
        prstart <- format(max(0, start-padding), scientific = F)
        prend <- format(min(len, end-padding), scientific = F)
        out[j] <- paste(chr, prstart, prend, val, sep="\t") 
        j <- j+1
      }
    }
    
  }
  return(out)
}

# Get command line parameters
args = commandArgs(trailingOnly=TRUE)

# Load Plasmid sequences
fasta <- read.fasta(args[1])
padding <- as.integer(args[2])
#fasta <- read.fasta("/home/caspar/Documents/Uni/MasterThesis/pipelines/plasmidIdentifier/out_test/testData/alignment/testData_padded.fasta")
#padding <- 1990

for (i in 1:length(fasta)) {
  print(paste("Reading contig", getName(fasta[[i]])))
  
  # Calculate GC content for Window sizes 50 and 1000
  GC50 <- GCwindow(getSequence(fasta)[[i]], 50)
  GC1000 <- GCwindow(getSequence(fasta)[[i]], 500)
  
  len <- length(getSequence(fasta[[i]]))
  
  # Write values to file
  vals1000 <- aggregateCircos(GC1000, 500, getName(fasta[[i]]), len)
  vals50 <- aggregateCircos(GC50, 50, getName(fasta[[i]]), len)
  
  # Shift values for padded fasta
  print(vals1000)
  print(vals50)
   
  lapply(vals50, write, "gc50.txt", append=TRUE)
  lapply(vals1000, write, "gc1000.txt", append=TRUE)
  
  #testplot(GC50, GC1000)
  
}

testplot <- function (val50, val1000, w1 = 50, w2 = 500, pad=padding) {
  
  plot((1:length(val50))*w1, val50, type="l", col="black")
  plot((1:length(val1000))*w2, val1000, type="l", add=TRUE, col="red")

}
    



