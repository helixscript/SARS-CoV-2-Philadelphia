library(ShortRead)
library(dplyr)
library(parallel)

fragLength <- 250
stepSize <- 5

s <- readFasta('VOCs.ff')@sread
t <- s[1]


qual <- paste0(rep('A', fragLength), collapse = '')

n <- 1
if(file.exists('R1.fastq')) file.remove('R1.fastq')
invisible(lapply(seq(1, width(t), by = stepSize), function(x){
  if(x+fragLength-1 <= width(t)){
    seq <- as.character(subseq(t, x, x+fragLength-1))
    write(paste0('@', paste0('s', n), '\n', seq, '\n+\n', qual), file = 'R1.fastq', append = TRUE)
   n <<- n+1
  }
}))

R1 <- readFastq('R1.fastq')
system('gzip R1.fastq')
writeFastq(reverseComplement(R1), file = 'R2.fastq.gz', compress = TRUE)

