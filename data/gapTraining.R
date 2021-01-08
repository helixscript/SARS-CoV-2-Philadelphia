library(ShortRead)
library(dplyr)

s <- readFasta('VOCs.ff')@sread
t <- s[2]

o <- Reduce('append', lapply(c(11292, 21767, 21991), function(n){
       Reduce('append', lapply(c((n-100):(n+100)), function(x){
         subseq(t, x, x+100)
       }))
     }))

names(o) <- paste0('s', 1:length(o))
o2 <- reverseComplement(o)

qual <- paste0(rep('A', 101), collapse = '')

file.remove(c('R1.fastq', 'R2.fastq'))
invisible(lapply(1:length(o), function(x){
  write(paste0('@', names(o[x]), '\n', as.character(o[x]), '\n+\n', qual), file = 'R1.fastq', append = TRUE)
  write(paste0('@', names(o2[x]), '\n', as.character(o2[x]), '\n+\n', qual), file = 'R2.fastq', append = TRUE)
}))

system('gzip R1.fastq R2.fastq')
system('mv R1.fastq.gz sequencing/simulated/UK_gaps_R1.fastq.gz')
system('mv R2.fastq.gz sequencing/simulated/UK_gaps_R2.fastq.gz')
