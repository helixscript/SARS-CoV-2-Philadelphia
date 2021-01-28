library(dplyr)
library(lubridate)
library(stringr)

seqDataDir <- 'data/sequencing'

samples <- read.table('data/samples.tsv', sep= '\t', header = TRUE, quote = '', stringsAsFactors = FALSE)

f <- list.files(seqDataDir, pattern = '^VSP', recursive = TRUE, full.names = TRUE)
d <- tibble(path = f[grepl('000000000', f)])

d$exp <- sapply(d$path, function(x){
           x <- unlist(strsplit(x, '/'))
           x <- x[length(x)]
           unlist(strsplit(x, '_'))[1]
         })

d$date <-  sapply(d$path, function(x){
             x <- unlist(strsplit(x, '/'))
             x <- x[length(x)-1]
             as.character(ymd(unlist(strsplit(x, '_'))[1]))
          })

d$days <- sapply(d$date, function(x) as.integer(ymd(x)))

d$VSP <- str_extract(d$exp, 'VSP\\d+')
d$subject <- samples[match(d$VSP, samples$VSP),]$patient_id
d$trial <- samples[match(d$VSP, samples$VSP),]$trial_id
d$sampleDate <- samples[match(d$VSP, samples$VSP),]$sample_date
d$run <- unlist(lapply(str_match_all(d$path, '000000000\\-([^\\/]+)'), '[', 2))

d <- arrange(d, desc(days)) %>% select(date, run, exp, trial, subject, sampleDate) %>% distinct()

s <- bind_rows(lapply(list.files('summaries/sampleSummaries', pattern = '*.tsv$', recursive = TRUE, full.names = TRUE), function(x){
       o <- read.table(x, sep = '\t', header = TRUE)
       o$Subject <- as.character(o$Subject)
       o
      }))

d$percentRefReadCoverage <- s[match(d$exp, s$exp),]$percentRefReadCoverage
d$lineage <- s[match(d$exp, s$exp),]$lineage

openxlsx::write.xlsx(d, file = 'summaries/seqRunSummary.xlsx')
