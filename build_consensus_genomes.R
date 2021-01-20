library(dplyr)
library(Biostrings)


# Collate sampe summary tables.
summary <- bind_rows(lapply(list.files('summaries/sampleSummaries', full.names = TRUE), function(x){
  t <- read.table(x, sep = '\t', header = TRUE)
  t$Subject <- as.character(t$Subject)
  t$largestContig <- as.numeric(t$largestContig)
  t$percentRefReadCoverage <- as.numeric(sub('%', '', t$percentRefReadCoverage))
  t$percentRefReadCoverage5 <- as.numeric(sub('%', '', t$percentRefReadCoverage5))
  t$sampleDate2 <- gsub('-', '', t$sampleDate)
  t
})) 

# %>% filter(percentRefReadCoverage5 >= 90, ! grepl('Control|PCR|Mock|simulated', sampleType, ignore.case = TRUE))


# Select a sample to represent each subject, sample type, time point combination.
# Select composite samples when available otherwise select the samples with the greated coverage.
# Break ties with read coverage percentages.

representativeSampleSummary <- function(summary, minPercentRefReadCoverage5){
  bind_rows(lapply(split(summary, paste(summary$Subject, summary$sampleType, summary$sampleDate2)), function(x){
    r <- subset(x, type == 'composite')
    if(nrow(r) > 0){
      r <- top_n(r, 1, wt = percentRefReadCoverage5) %>% dplyr::slice(1)
    } else{
      r <- top_n(x, 1, wt = percentRefReadCoverage5) %>% dplyr::slice(1)
    }
    r
  })) %>% dplyr::filter(percentRefReadCoverage5 >= minPercentRefReadCoverage5)
}

representativeSampleSummary_90 <- representativeSampleSummary(summary, 90)


# Retrieve consensus sequences from representative samples.
retrieveConcensusSeqs <- function(summary){
  Reduce('append', lapply(1:nrow(summary), function(x){
    x <- summary[x,]
    f <- paste0('summaries/VSPdata/', x$exp, '.', ifelse(x$type == 'composite', 'composite', 'experiment'), '.RData')
    if(! file.exists(f)) stop('Can not locate VSP data file -- ', f)
    load(f)
    d <- DNAStringSet(opt$concensusSeq)
    names(d) <- paste0(x$Subject, '|', x$sampleType, '|', x$sampleDate2)
    d
  }))
}

concensusSeqs90_5 <- retrieveConcensusSeqs(representativeSampleSummary_90)

writeXStringSet(concensusSeqs90_5, file = 'summaries/allGenomes_90_5.fasta')

# Run in pangolin CONDA environment.
# conda activate pangolin
# pangolin -t 25 -o summaries/allGenomes_90_5.pangolin summaries/allGenomes_90_5.fasta

p <- read.table('summaries/allGenomes_90_5.pangolin/lineage_report.csv', sep = ',', header = TRUE)
unlink('summaries/allGenomes_90_5.pangolin', recursive = TRUE)
openxlsx::write.xlsx(p[,1:3], file = 'summaries/allGenomes_90_5.pangolin.xlsx')

