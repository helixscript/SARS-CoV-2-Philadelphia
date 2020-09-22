library(tidyverse)
library(ShortRead)
library(lubridate)
library(parallel)

blastn             <- '/home/everett/ext/blast+/bin/blastn'
makeblastdb        <- '/home/everett/ext/blast+/bin/makeblastdb'
refGenomes         <- '/home/everett/projects/SARS-CoV-2-Philadelphia/data/genomes/GISAID/sequences_2020-08-07_17-50.fasta'
refGenomesMetaData <- '/home/everett/projects/SARS-CoV-2-Philadelphia/data/genomes/GISAID/metadata_2020-08-07_19-18.tsv'
query              <- '/home/everett/projects/SARS-CoV-2-Philadelphia/summaries/concensusSeqs95.fasta'
patientMetaData    <- '/home/everett/projects/SARS-CoV-2-Philadelphia/data/subjectMetaData.tsv'

minGenomeAlignmentPercentage <- 95

tmpFile <- function(){ paste0('tmp.', paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = '')) }

refGenomesMetaData <- read.table(refGenomesMetaData, sep = '\t', comment.char = '', quote = '', fill = TRUE, header = TRUE)
refGenomesMetaData$date2 <- ymd(refGenomesMetaData$date)
refGenomesMetaData$region2 <- paste0(refGenomesMetaData$country_exposure, '|', refGenomesMetaData$division_exposure)

patientMetaData <- read.table(patientMetaData, sep = '\t', comment.char = '', quote = '', fill = TRUE, header = TRUE, check.names = FALSE)
patientMetaData$date2 <- mdy(patientMetaData$`Symptom Onset Date`)

refGenomes <- readFasta(refGenomes)

q <- readFasta(query)

cluster <- makeCluster(10)
clusterExport(cluster, c('blastn', 'makeblastdb', 'refGenomes', 'refGenomesMetaData', 'query', 
                         'patientMetaData', 'minGenomeAlignmentPercentage', 'tmpFile'))

r <- bind_rows(parLapply(cluster, split(q, 1:length(q)), function(x){
  library(ShortRead)
  library(dplyr)
  library(stringr)
  
  patientID <- unlist(str_split(as.character(x@id), '\\|'))[1]
  patientID <- sub('\\-TCE$', '', patientID)
  
  patientsymptomOnset <- patientMetaData[match(patientID, patientMetaData$`Subject ID`),]$date2
  m <- refGenomesMetaData[refGenomesMetaData$date2 < patientsymptomOnset,]
  
  t <- tmpFile()
  
  refStrains <- refGenomes[as.character(refGenomes@id) %in% m$strain] 

  # Make sure that the reference strains are withing +/- 5% of the length of the genome we are comapring.
  # refStrains <- refStrains[width(refStrains) >= (width(x) - width(x)*0.05) & width(refStrains) <= (width(x) + width(x)*0.05)]
  
  if(length(refStrains) <= 1) return(tibble()) 
  
  writeFasta(refStrains, file = paste0(t, '.ref_fasta'))
  system(paste0(makeblastdb, ' -in ', paste0(t, '.ref_fasta'), ' -dbtype nucl'))
  
  writeFasta(x, file = paste0(t, '.query_fasta'))
  system(paste0(blastn, ' -query ', paste0(t, '.query_fasta'),  ' -db ', paste0(t, '.ref_fasta'), ' -out ', 
                paste0(t, '.blast'), ' -outfmt 6 -max_target_seqs 1000' ))
  
  b <- read.table(paste0(t, '.blast'), sep = '\t', header = FALSE)
  names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
  b$alignmentLength <- b$qend - b$qstart + 1
  b$alignmentPercent <- (b$alignmentLength / width(x))*100
  b <- subset(b, alignmentPercent >= minGenomeAlignmentPercentage)
  
  minMisMatch <- min(b$mismatch)
  invisible(file.remove(list.files(pattern = paste0('^',t))))
  b$patientsymptomOnset <- patientsymptomOnset
  b$refGenomeSearched <- length(refStrains) 
  b <- left_join(b, dplyr::select(refGenomesMetaData, strain, region2), by = c("sseqid" = "strain"))
  subset(b, mismatch == minMisMatch)
}))

stopCluster(cluster)

nearestGenomes <- bind_rows(lapply(split(r, r$qname), function(x){
  a <- sort((table(x$region2)/n_distinct(x$sseqid))*100, decreasing = TRUE)
  
  if(length(a) > 3){
    b <- a[1:3]
  } else if (length(a) == 2){
    b <- a[1:2]
  } else {
    b <- a
  }
  
  p <- paste(sprintf("%.1f%%", b), names(b), collapse = ';  ')
  
  if(length(a) > 3) p <- paste0(p, '; Other ', sprintf("%.1f%%", 100 - sum(a[1:3])))
  
  tibble(subject = x$qname[1], 
         patientsymptomOnset = x$patientsymptomOnset[1],
         genomesSearched = x$refGenomeSearched[1],
         nGenomeHits = n_distinct(x$sseqid),
         nUniqueGenomeHits = n_distinct(as.character(refGenomes[as.character(refGenomes@id) %in% x$sseqid]@sread)),
         nMismatches = x$mismatch[1], 
         minPercentCoverage = sprintf("%.1f%%", min(x$alignmentPercent)), 
         genomeHits = p) 
}))

openxlsx::write.xlsx(nearestGenomes, file = 'summaries/nearestGenomes.xlsx')
