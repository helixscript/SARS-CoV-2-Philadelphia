library(tidyverse)
library(Biostrings)
source('lib/lib.R')

samples <- removeProblematicSamples(read.table('data/samples.tsv', sep= '\t', header = TRUE, quote = '', stringsAsFactors = FALSE))

invisible(lapply(list.dirs('summaries/trials'), function(x){
  x <- unlist(strsplit(x, '/'))
  
  if(length(x) == 2) return()

  d <- bind_rows(lapply(list.files(file.path('summaries', 'sampleSummaries', x[3]), full.names = TRUE, pattern = '.tsv$', recursive = TRUE), function(y){
    t <- read.table(y, sep = '\t', header = TRUE)
    t$Subject <- as.character(t$Subject)
    t$trial_id <- x[3]
    t$percentRefReadCoverage5 <- as.numeric(sub('%', '', t$percentRefReadCoverage5))
    t$sampleDate2 <- gsub('-', '', t$sampleDate)
    t
  })) 
  
  # Collate the major variant tables for genomes with >= 50% coverage. Genomes with less coverage 
  # are prone to odd variant calls due to alignment edge effects.
  
  if(nrow(d) == 0) return()
  
  o <- representativeSampleSummary(d, minPercentRefReadCoverage5 = 50) %>% pull(exp)
  
  if(length(o) > 0){
    r <- bind_rows(lapply(o, function(e){
      f <- file.path('summaries', 'VSPdata', paste0(e, ifelse(grepl('-', e), '.experiment.RData', '.composite.RData')))
      if(! file.exists(f)) stop('Error - could not locate ', f)
      load(f)
      if(nrow(opt$variantTableMajor) == 0) return(tibble())
      tibble(experiment = e, 
             position = opt$variantTableMajor$POS,
             variant = paste0(opt$variantTableMajor$REF, opt$variantTableMajor$POS, opt$variantTableMajor$ALT),
             gene = opt$variantTableMajor$genes,
             type = opt$variantTableMajor$type)
    }))
  
    if(nrow(r) > 0){
      bind_rows(lapply(split(r, r$position), function(e){
                  tibble(position = e$position[1],
                        variants = paste0(unique(e$variant), collapse = ', '),
                        gene = paste0(unique(e$gene), collapse = ', '),
                        types = paste0(unique(e$type), collapse = ', '),
                        experiments = n_distinct(e$experiment))
                 })) %>% arrange(desc(experiments)) %>%
      openxlsx::write.xlsx(file.path(x[1], x[2], x[3], 'variantSummary.xlsx'))
    }
  }
  
  # Identify the best data object for each subject -- typically the composite object.
  representativeSampleSummary_90 <- representativeSampleSummary(d, 90)
  
  if(nrow(representativeSampleSummary_90) > 0){
    # Create FASTA files for genomes with >= 90% coverage.
    concensusSeqs90_5 <- retrieveConcensusSeqs(representativeSampleSummary_90)
    writeXStringSet(concensusSeqs90_5, file = file.path(x[1], x[2], x[3], 'highQualGenomes.fasta'))
  
    lineages <- bind_rows(lapply(names(concensusSeqs90_5), function(x){
      x <- unlist(strsplit(x, '\\|'))
      tibble(sample = paste(x[1:4], collapse = '|'), lineage = x[5])
    }))
  
    openxlsx::write.xlsx(lineages, file = file.path(x[1], x[2], x[3], 'lineages.xlsx'))
  }
}))
