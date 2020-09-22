library(ShortRead)
library(dplyr)
library(lubridate)
options(stringsAsFactors = FALSE)

gisaid <- readFasta('GISAID/sequences_2020-08-07_17-50.fasta')

g <- readFasta('../../summaries/concensusSeqs95.fasta')

m <- Reduce('append', list(gisaid, g))
writeFasta(m, file = 'm.fasta')