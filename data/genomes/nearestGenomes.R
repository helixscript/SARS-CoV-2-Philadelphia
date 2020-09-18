library(ShortRead)
library(dplyr)
library(lubridate)
options(stringsAsFactors = FALSE)

gisaid <- readFasta('GISAID/sequences_2020-08-07_17-50.fasta')