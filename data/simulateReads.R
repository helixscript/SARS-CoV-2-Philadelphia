library(ShortRead)
library(dplyr)
library(parallel)

fragLength <- 500
reads <- 1000000

s <- readFasta('VOCs.ff')@sread
t <- s[2]

set.seed(1)
o <- data.frame(start = sample(1:width(t), reads, replace = TRUE))
o$end <- o$start + fragLength
o <- dplyr::filter(o, end <= width(t))
o$id <- paste0('s', 1:nrow(o))

qual <- paste0(rep('A', fragLength+1), collapse = '')
cluster <- makeCluster(25)

tmpFile <- function(){ paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') }

clusterExport(cluster, c('t', 'qual', 'tmpFile'))

invisible(parLapply(cluster, split(o, ntile(1:nrow(o), 25)), function(x){
  library(Biostrings)
  f <- tmpFile()
  invisible(lapply(split(x, 1:nrow(x)), function(x2){
    seq <- as.character(subseq(t, x2$start, x2$end))
    write(paste0('@', x2$id, '\n', seq, '\n+\n', qual), file = f, append = TRUE)
   }))
}))

system(paste('cat', paste(list.files(pattern = '.tmp$'), collapse = ' '), ' > R1.fastq'))
o <- readFastq('R1.fastq')
file.remove('R1.fastq')
i <- as.character(o@id)
o <- o[order(i)]
writeFastq(o, file = 'R1.fastq.gz', compress = TRUE)
file.remove(list.files(pattern = '*.tmp$'))

o2 <- reverseComplement(o)
writeFastq(o2, file = 'R2.fastq.gz', compress = TRUE)
