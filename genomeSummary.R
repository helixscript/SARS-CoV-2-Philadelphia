library(dplyr)
library(readr)

r <- bind_rows(lapply(list.files(path = 'summaries/sampleSummaries', 
                                 pattern = '*.tsv', 
                                 recursive = TRUE, 
                                 full.names = TRUE), function(x){
                                   o <- read_delim(x, delim = '\t', 
                                                   col_types = cols(Subject = col_character(), 
                                                                    sampleDate = col_date(),
                                                                    genomes = col_character(),
                                                                    percentRefReadCoverage = col_number(),
                                                                    percentRefReadCoverage5 = col_number()))
                                   t <- unlist(strsplit(x, '/'))
                                   o$trial <- t[length(t)-1]
                                   o$sample <- sub('\\-\\d+[a-z]?$', '', o$exp, perl = TRUE)
                                   select(o, trial, Subject, exp, sample, type, sampleDate, largestContig, lineage, percentRefReadCoverage, percentRefReadCoverage5)
     })) %>%
     group_by(sample) %>%
     top_n(1, percentRefReadCoverage5) %>%
     dplyr::slice(1) %>%
     ungroup() %>%
     select(trial, Subject, sample, lineage, percentRefReadCoverage, percentRefReadCoverage5)



nrow(subset(r, percentRefReadCoverage  >= 90))
nrow(subset(r, percentRefReadCoverage5 >= 90))
nrow(subset(r, percentRefReadCoverage  >= 95))
nrow(subset(r, percentRefReadCoverage5 >= 95))
