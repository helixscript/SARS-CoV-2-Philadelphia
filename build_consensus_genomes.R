library(dplyr)
library(Biostrings)
source('lib/lib.R')

# Collate sampe summary table and omit problematic samples and subjects.
summary <- removeProblematicSamples(bind_rows(lapply(list.files('summaries/sampleSummaries', full.names = TRUE), function(x){
  t <- read.table(x, sep = '\t', header = TRUE)
  t$Subject <- as.character(t$Subject)
  t$largestContig <- as.numeric(t$largestContig)
  t$percentRefReadCoverage  <- as.numeric(sub('%', '', t$percentRefReadCoverage))
  t$percentRefReadCoverage5 <- as.numeric(sub('%', '', t$percentRefReadCoverage5))
  t$sampleDate2 <- gsub('-', '', t$sampleDate)
  t
})))


# Identify the best data object for each subject -- typically the composite object.
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

system('./pangolin.sh')

p <- read.table('summaries/allGenomes_90_5.pangolin/lineage_report.csv', sep = ',', header = TRUE)
unlink('summaries/allGenomes_90_5.pangolin', recursive = TRUE)
openxlsx::write.xlsx(p[,1:3], file = 'summaries/allGenomes_90_5.pangolin.xlsx')

