library(tidyverse)
library(parallel)
library(Biostrings)
library(lubridate)
source('lib/lib.R')
options(stringsAsFactors = FALSE)
overWriteSubjectReports <- FALSE
CPUs <- 10

# Read in sample data table.
samples <- removeProblematicSamples(read.table('data/samples.tsv', sep= '\t', header = TRUE, quote = '', stringsAsFactors = FALSE))

# Read in input data values.
sampleInputs <- read.table('data/sampleInputs.tsv', sep= '\t', header = TRUE, stringsAsFactors = FALSE)


# Trim leading and trailing white space from metadata.
trimLeadingTrailingWhtSpace <- function(x) gsub('^\\s+|\\s+$', '', x)
samples$sample_id   <- sapply(samples$sample_id, trimLeadingTrailingWhtSpace)
samples$trial_id    <- sapply(samples$trial_id, trimLeadingTrailingWhtSpace)
samples$patient_id  <- sapply(samples$patient_id, trimLeadingTrailingWhtSpace)
samples$sample_date <- sapply(samples$sample_date, trimLeadingTrailingWhtSpace)
samples$sample_type <- sapply(samples$sample_type, trimLeadingTrailingWhtSpace)
samples$VSP         <- sapply(samples$VSP, trimLeadingTrailingWhtSpace)


# Standardize date formatting
samples$sample_date <- as.character(mdy(samples$sample_date))


# Create a vector of VSP ids with sequencning data.
availableVSPs <- unique(str_extract(list.files('summaries/VSPdata'), 'VSP\\d+'))


# Identify potential sample duplications.
sampleUniquenessCheck <- 
  group_by(filter(samples, VSP %in% availableVSPs), patient_id, sample_date, sample_type) %>%
  summarise(n = n(), VSPids = paste(VSP, collapse = ', ')) %>%
  ungroup() %>%
  filter(n > 1)


# Create subject summary reports for all subjects with experimental data (VSP data objects in VSP_data/).
cluster <- makeCluster(CPUs)
clusterExport(cluster, c('samples', 'sampleInputs', 'overWriteSubjectReports'))
if(! dir.exists('summaries/patientReportsData')) dir.create('summaries/patientReportsData')


# Start subject report log.
write(date(), file = 'logs/buid_subject_reports.log', append = FALSE)


# Create a table of patient ids with a chunking vector.
d <- dplyr::select(samples, patient_id, trial_id) %>% dplyr::distinct() %>% dplyr::mutate(s = ntile(1:n(), CPUs))


#d <- d[grepl('DOH', d$patient_id),]

# Build patient reports.
invisible(parLapply(cluster, split(d, d$s), function(p){
#invisible(lapply(split(d, d$s), function(p){  
  library(tidyverse)
  
  invisible(lapply(split(p, paste(p$patient_id, p$trial_id)), function(r){
    
    x <- subset(samples, patient_id == r$patient_id & trial_id == r$trial_id)
  
    dir1 <- file.path('summaries/patientReports', r$trial_id)
    dir2 <- file.path('summaries/patientReportsData', r$trial_id)
    
    if(! dir.exists(dir1)) dir.create(dir1)
    if(! dir.exists(dir2)) dir.create(dir2)
    
    if(overWriteSubjectReports == FALSE & file.exists(file.path(dir1, paste0(r$patient_id, '.pdf')))) return()
    
    files <- list.files('summaries/VSPdata', pattern = paste0(x$VSP, collapse = '|'), full.names = TRUE)

    if(length(files) == 0){
      write(paste0('[.] No VSP data files for subject ', x$patient_id[1]), file = 'logs/buid_subject_reports.log', append = TRUE)
      return()
    }
    
    dat <- lapply(files, function(f){
           # Here we load VSP data files where there are two types of files, analyses of individual sequencing 
           # experiments and composite analyses where multiple experiments where combined.
             load(f)
             opt$vsp <- str_extract(f, 'VSP\\d+')
             opt$seq_sample <-  str_extract(f, 'VSP\\d+\\-?\\d+?[abm]?')
           
             d <- subset(samples, VSP == opt$vsp)
             opt$sample_id <- d$sample_id[1]
             opt$trial_id <- d$trial_id[1]
             opt$subject <- d$patient_id[1]
             opt$date <- d$sample_date[1]
             opt$sampleType = d$sample_type[1]
             opt$genomesPerMicroLiter = as.numeric(d$N1_copies_per_ul[1])
             if(! grepl('composite', f)){
               opt$inputGenomes <- opt$genomesPerMicroLiter * subset(sampleInputs, Sample_Name == opt$seq_sample)$uL_Inputs
             } else {
               opt$inputGenomes <- NA
             }
             opt$concensusSeq = as.character(opt$concensusSeq)
             return(opt)
           })
         
    names(dat) <- unlist(lapply(dat, '[[', 'seq_sample'))
         
    save(dat, file = file.path(dir2, paste0(r$patient_id, '.RData')))
    
    result = tryCatch({
              rmarkdown::render('report.Rmd',
                                output_file = file.path(dir1, paste0(r$patient_id, '.pdf')),
                                params = list('date'  = format(Sys.time(), "%Y-%m-%d"),
                                              'title' = paste0('COVID-19 subject ', r$patient_id)))
              }, error = function(e) {
                 write(paste0('[!] Failed to create subject report for ', dat[[1]]$subject), file = 'logs/buid_subject_reports.log', append = TRUE)
              })
    }))
}))




# Create plots of sequence coverage vs. sample inputs using different scale cutoffs to better visualize the data.

samples$N1_copies_per_ul <- as.numeric(samples$N1_copies_per_ul)
d <- samples[!is.na(samples$N1_copies_per_ul), c('VSP', 'N1_copies_per_ul')]
d <- subset(d, VSP %in% stringr::str_extract(list.files('summaries/VSPdata'), 'VSP\\d+'))
d <- bind_rows(lapply(split(d, d$VSP), function(x){
      files <- list.files('summaries/VSPdata')
      files <- files[grep(paste0('^', x$VSP), files)]
      p <- unlist(lapply(files, function(f){
        load(file.path('summaries/VSPdata', f))
        opt$refGenomePercentCovered_5reads
      }))
      x$coverage5 <- max(p)
      x$files <- length(p)
      x
    }))
  
invisible(lapply(c('Inf', '1e6', '1e4'), function(x){
            p1 <- ggplot(subset(d, N1_copies_per_ul <= as.numeric(x)), aes(N1_copies_per_ul, coverage5)) + 
                  theme_bw() +
                  scale_y_continuous(labels = scales::percent) +
                  geom_point() +
                  labs(x = 'Copies per ul', y = 'SARS-CoV2 coverage (>= 5 reads per position)') +
                  ggtitle(paste0('Concentration limit: ', x))
    
            ggsave(p1, file = paste0('summaries/inputGenomesVsCoverage_', x, '_cutoff.pdf'), units = 'in', height = 5)
          }))
