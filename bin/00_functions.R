####################
# Helper functions #
####################

# Calculates the GC content in a sliding window for a character vector
GCwindow <- function(seq, width, circular=T, bases=c("c", "g")) {
  n <- length(seq)
  gc <- double(n)
  w <- floor(width/2)
  
  for (i in 1:length(seq)){

    # Extract subsequence in the current window
    # Start at the start of sequence if circular = T
    if (i < w) {
      window <- c(seq[(n-(w-i)):n], seq[1:i+w]) 
    } else if ((i+w) > n) {
      window <- c(seq[(-w + i):n], seq[1:(w-(n-i))])
    } else {
      window <- seq[(i-w):(i+w)]
    }
    
    # Calculate base occurences
    t <- table(window)
    b1 <- t[names(t)==bases[1]]
    b2 <- t[names(t)==bases[2]]
    
    gc[i] <- sum(b1, b2)/width
  }
  return(gc)
}

# Create a table with midpoint and GC value for given window
aggregateGC <- function(v, width) {
  n <- floor(length(v) / width)
  w <- floor(width/2)
  dt <- data.table(position=integer(n), gc=numeric(n))
  
  for (i in 1:n){
    end <- (i*width)
    dt[i, position:= end-w]
    start = max(0, end-(2*w))
    dt[i, gc:= mean(v[start:end])]
  }
  return(dt)
}


# Write a vector with numerical values as Circos data output
# Width parameter aggregates the data (mean value) to reduce the number of datapoints
aggregateCircos <- function(v, width, chr = "C") {
  
  n <- floor(length(v) / width)
  out <- character(n)
  
  for (i in 1:n){
   end <- min(i * width, length(v))
   start <- max(0, end - width)
   val <- mean(v[start:end])
   out[i] <- paste(chr, format(start, scientific=F), format(end, scientific=F), val) 
  }
  return(out)
}

readGff <- function(fileName, sep='=') {
  require(data.table)
  
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
