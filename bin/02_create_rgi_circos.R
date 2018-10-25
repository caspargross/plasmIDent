#!/usr/bin/env Rscript

######################
# Create CIRCOS data #
#  tracks from rgi   #
######################

require(data.table)

readGff <- function(fileName, sep='=') {
  
  dt <- fread(fileName, header=F)  
  dt <- na.omit(dt)

  setnames(dt, c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))
  att <- dt[,tstrsplit(attribute, ';')]
  
  
  # Extract column names 
  att_names <- unname(unlist(att[1, .(sapply(.SD, function(x) {tstrsplit(trimws(x), sep)[[1]][1]})),]))
  setnames(att, att_names)
  att_names <- att_names[!is.na(att_names)]
  
  # Combine data.tables
  att[,(att_names) := sapply(.SD, function(x) {tstrsplit(trimws(x), sep, fixed = T, keep=2)}), .SDcols = att_names]

  dt <- cbind(dt[,!"attribute"], att[, ..att_names])
  
  return(dt)
}


args = commandArgs(trailingOnly=TRUE)
file = args[1]

rgi <- readGff(file)

# Export RGI
fwrite(rgi[,.(seqname, start, end, Name)], sep="\t", file ="rgi.txt", col.names = F)
fwrite(rgi[,.(seqname, start, end)], sep="\t", file ="rgi_span.txt", col.names = F)
