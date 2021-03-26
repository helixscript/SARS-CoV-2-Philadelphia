removeProblematicSamples <- function(x){
  if('exp' %in% names(x))        x <- dplyr::filter(x, ! grepl('VSP0069', exp))
  if('Subject' %in% names(x))    x <- dplyr::filter(x, Subject != '237')
  if('patient_id' %in% names(x)) x <- dplyr::filter(x, patient_id != '237')
  if('patient_id' %in% names(x)) x <- dplyr::filter(x, ! grepl('MisC', patient_id, ignore.case = TRUE))
  #if('trial_id' %in% names(x))   x <- dplyr::filter(x, ! grepl('ODoherty', trial_id, ignore.case = TRUE))
  x
}

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


# Retrieve consensus sequences from representative samples.
retrieveConcensusSeqs <- function(summary){
  Reduce('append', lapply(1:nrow(summary), function(x){
    x <- summary[x,]
    f <- paste0('summaries/VSPdata/', x$exp, '.', ifelse(x$type == 'composite', 'composite', 'experiment'), '.RData')
    if(! file.exists(f)) stop('Can not locate VSP data file -- ', f)
    load(f)
    d <- DNAStringSet(opt$concensusSeq)
    
    names(d) <- gsub('\\s+', '_', paste0(x$trial_id, '|', x$Subject, '|', x$sampleType, '|', x$sampleDate2, '|', x$lineage))
    d
  }))
}