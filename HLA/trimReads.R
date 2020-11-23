library(dplyr)
library(ShortRead)

f <- tibble(file = list.files('seqData', full.names = FALSE),
            sample = unlist(lapply(strsplit(file, '_'), '[', 1)))

invisible(lapply(split(f, f$sample), function(x){
  x <- arrange(x, file)
  R1 <- trimTailw(readFastq(paste0('seqData/', x[1,]$file)), 2, 'A', 5)
  R1 <- R1[width(R1) >= 75]
  R1@id <- BStringSet(sub('\\s+.+', '', as.character(R1@id)))
  R2 <- trimTailw(readFastq(paste0('seqData/', x[2,]$file)), 2, 'A', 5)
  R2 <- R2[width(R2) >= 75]
  R2@id <- BStringSet(sub('\\s+.+', '', as.character(R2@id)))
  
  i <- base::intersect(as.character(R1@id), as.character(R2@id))
  R1 <- R1[as.character(R1@id) %in% i]
  R2 <- R2[as.character(R2@id) %in% i]
  
  writeFastq(R1, paste0('trimmedSeqData/', x[1,]$file), compress = TRUE)
  writeFastq(R2, paste0('trimmedSeqData/', x[2,]$file), compress = TRUE)
}))
