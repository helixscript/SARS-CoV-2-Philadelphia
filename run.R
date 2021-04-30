library(parallel)
library(dplyr)
library(ggplot2)
library(parallel)
library(Biostrings)
library(lubridate)
library(stringr)

CPUs          <- 30
overwrite     <- FALSE
updateGenomes <- FALSE
Rscript       <- '/home/opt/R-3.4.0/bin/Rscript'
mafftPath     <- '/home/everett/ext/mafft/bin/mafft'
softwareDir   <- '/home/everett/projects/SARS-CoV-2-Philadelphia'
workDir       <- file.path(softwareDir, 'scratch')
VOIs          <- c("B.1.1.7", "B.1.351", "B.1.427", "B.1.429", "B.1.525", "B.1.526", "P.1", "P.2")

remoteOutputDir <- '/data/SARS-CoV-2'
sequenceDataDir <- '/data/SARS-CoV-2/sequencing'
sampleData      <- '/data/SARS-CoV-2/samples.tsv'
sampleInputData <- '/data/SARS-CoV-2/sampleInputs.tsv'
VSPdataDir      <- file.path(softwareDir, 'summaries/VSPdata')


# Check remote system for working file.
if('working' %in% list.files(remoteOutputDir)) q()


# Create remote working file to let users know the pipeline is running.
system('ssh microb120 touch /media/lorax/data/SARS-CoV-2/working')


# Internal pipepline settings.
options(stringsAsFactors = FALSE)
overWriteSubjectReports <- TRUE
source(paste0(softwareDir, '/lib/lib.R'))


# Clear previous log files.
system(paste0('rm ', softwareDir, '/summaries/logs/*'))
logFile <- paste0(softwareDir, '/summaries/logs/log')
write(date(), file = logFile)


# Read in sample data table.
samples <- removeProblematicSamples(read.table(sampleData, sep= '\t', header = TRUE, quote = '', stringsAsFactors = FALSE))


# Create a list of black listed runs.
d <- list.files(sequenceDataDir, recursive = TRUE, pattern = '^VSP', full.names = TRUE)
r <- sapply(d, function(x){
  o <- unlist(strsplit(x, '/'))
  runID <- o[length(o)-1]
  if(toupper(substr(runID, 1, 1)) == 'X'){
    return(runID)
  } else {
    return(NA)
  }
})

r <- unique(r[!is.na(r)])


# Identify VSP objects created from black listed runs and remove them from the current VSP repository
# Remove reports if reports come from an experiment in a black listed run and the data object currently exists. 
if(length(r) > 0){
  for(i in r){
    o <- unique(stringr::str_extract(d[grep(i, d)], 'VSP\\d+\\-\\d+[a-z]?'))
    o2 <- unique(sub('\\-\\d+[a-z]?', '', o))
    
    f <- file.path(softwareDir, 'summaries', 'VSPdata',  paste0(o2, '.composite.RData'))
    f <- f[file.exists(f)]
    if(length(f) > 0){
      message('Removing report: ', f)
      file.remove(f)
    }
    
    f <- file.path(softwareDir, 'summaries', 'VSPdata', paste0(o, '.experiment.RData'))
    f <- f[file.exists(f)]
    if(length(f) > 0){
      updateGenomes <- TRUE
      file.remove(f)
      o <- select(subset(samples, VSP %in% stringr::str_extract(f, 'VSP\\d+')), trial_id, patient_id)
      message('Removing reports: ', as.character(o))
      file.remove(paste0(softwareDir, '/summaries/trials/', o$trial_id, '/', o$patient_id, '.pdf'))
      file.remove(paste0(softwareDir, '/summaries/sampleSummaries/', o$trial_id, '/', o$patient_id, '.sampleSummary.tsv'))
    }
  }
}


# Create a list of all available sequencing data.
f1 <- unname(sapply(list.files(sequenceDataDir, recursive = TRUE, pattern = '^VSP', full.names = TRUE), function(x){
  o <- unlist(strsplit(x, '/'))
  runID <- o[length(o)-1]
  if(! toupper(substr(runID, 1, 1)) == 'X'){
    return(unlist(strsplit(o[length(o)], '_'))[1])
  } else {
    return(NA)
  }
}))
f1 <- unique(f1[!is.na(f1)])


# Create a list of all VSP data objects already created.
f2 <- list.files(VSPdataDir, recursive = TRUE, pattern = '^VSP', full.names = FALSE)
f2 <- unique(unlist(lapply(strsplit(f2, '\\.'), '[', 1)))

# Identify new incoming experiments. 
f3 <- f1[! f1 %in% f2]


if(length(f3) > 0){
  updateGenomes <- TRUE
  VSPs <- unique(sub('\\-\\d+[a-z]?', '', f3))

  # Remove patient reports and sample summaries so that they are recreated.
  invisible(sapply(VSPs, function(x){
    x <- subset(samples, VSP == x)
    if(nrow(x) == 0) return() 
    if(nrow(x) == 1 & file.exists(file.path(softwareDir, 'summaries', 'trials', x$trial_id, paste0(x$patient_id, '.pdf')))){
       message('Removing reports for: ', x$trial_id, ' / ', x$patient_id)
       file.remove(file.path(softwareDir, 'summaries', 'trials', x$trial_id, paste0(x$patient_id, '.pdf')))
       file.remove(file.path(softwareDir, 'summaries', 'sampleSummaries', x$trial_id, paste0(x$patient_id, '.sampleSummary.tsv')))
     }
  }))

  d <- data.frame(vsp = VSPs)
  d$s <- ntile(1:nrow(d), CPUs)
  
  cluster <- makeCluster(CPUs)
  clusterExport(cluster, c('Rscript', 'overwrite', 'sequenceDataDir', 'softwareDir', 'sequenceDataDir', 
                           'VSPdataDir', 'sampleData', 'sampleInputData', 'workDir', 'tmpFile'))
  
  invisible(parLapply(cluster, split(d, d$s), function(x){
  #invisible(lapply(split(d, d$s), function(x){
  library(dplyr)
  library(stringr)
    
  logFile <- paste0(softwareDir, '/summaries/logs/build_VSP_data.', x$s[1], '.log')
  write(date(), file = logFile, append = FALSE)
  
  invisible(sapply(x$vsp, function(vsp){
    
    # All the files to draw R1 and R2 read files from for VSP.
    sampleFiles <- list.files(sequenceDataDir, recursive = TRUE, pattern = paste0('^', vsp, '-'), full.names = FALSE)
  
    # Remove black listed runs.
    sampleFiles <- sampleFiles[! toupper(substr(sampleFiles, 1, 1)) == 'X']
    if(length(sampleFiles) == 0) return()
    if(! length(sampleFiles) %% 2 == 0) return()
    
    # Split experiment files into a list where each key contains the R1 and R2 sequencing files.
    s <- unlist(lapply(strsplit(unlist(lapply(strsplit(sampleFiles, '/'), '[', 2)), '_'), '[', 1))
    
    sampleRunFiles <- split(sampleFiles, s)

    write(paste0('Starting ', vsp, ' with ', length(sampleRunFiles), ' seq experiments.'), file = logFile, append = TRUE)
    
    if(length(sampleRunFiles) > 1){
      outputFile <- paste0(softwareDir, '/summaries/VSPdata/', vsp, '.composite.RData')
      
      ### if(overwrite == FALSE & file.exists(outputFile)){
      ###     write(paste0(outputFile , ' exists -- skipping.'), file = logFile, append = TRUE)
      ### } else {
        if(file.exists(outputFile)) file.remove(outputFile)
        
        comm <- paste0(Rscript, ' ', softwareDir, '/assemble.R --outputFile ', outputFile, ' --workDir ', file.path(workDir, tmpFile()), 
                       ' --softwareDir ', softwareDir, ' --R1 ', 
                       paste0(sequenceDataDir, '/', sampleFiles[grepl('_R1_', sampleFiles)], collapse = ','), ' --R2 ',
                       paste0(sequenceDataDir, '/', sampleFiles[grepl('_R2_', sampleFiles)], collapse = ','))
                  
        write(paste0(date(), ' (', outputFile, ') ', comm, '\n'), file = logFile, append = TRUE)
        
        tryCatch({
           system(comm)
        }, error = function(e) {
           write(paste0('\n[!] (', outputFile, ') data object creation failed, system error caught: "', e, '".\n'), file = logFile, append = TRUE)
        })
        
        if(! file.exists(outputFile)) write(paste0('\n[!] (', outputFile, ') data object creation failed.'), file = logFile, append = TRUE)
      ### }
    }
    
    invisible(mapply(function(n, x){
      outputFile <- paste0(softwareDir, '/summaries/VSPdata/', n, '.experiment.RData')
         
       if(overwrite == FALSE & file.exists(outputFile)){
         write(paste0(outputFile , ' exists -- skipping.'), file = logFile, append = TRUE)
       } else {
         comm <-  paste0(Rscript, ' ', softwareDir, '/assemble.R --outputFile ', outputFile, ' --workDir ', file.path(workDir, tmpFile()), ' --R1 ', 
                         paste0(sequenceDataDir, '/', x[grepl('_R1_', x)], collapse = ','), ' --R2 ',
                         paste0(sequenceDataDir, '/', x[grepl('_R2_', x)], collapse = ','))
         
         write(paste0(date(), ' (', outputFile, ') ', comm, '\n'), file = logFile, append = TRUE)
         
         tryCatch({
                     system(comm)
                  }, error = function(e) {
                     write(paste0('\n[!] (', outputFile, ') data object creation failed, system error caught: "', e, '".\n'), file = logFile, append = TRUE)
                  })
         
           if(! file.exists(outputFile)) write(paste0('\n[!] (', outputFile, ') data object creation failed.\n'), file = logFile, append = TRUE)
         }
      }, names(sampleRunFiles), sampleRunFiles))
    
      write(paste0('Completed ', vsp, '\n'), file = logFile, append = TRUE)
    }))
  }))

  stopCluster(cluster)
  system(paste0('find ', softwareDir, '/summaries/logs -name build_VSP_data* | xargs grep failed > ', softwareDir, '/summaries/logs/VSP_data.failed'))
}

if(updateGenomes == FALSE)
{
  system('ssh microb120 rm /media/lorax/data/SARS-CoV-2/working')
  write(c('done'), file = logFile, append = TRUE)
  q()
}


# Prepare data for subject report creation
# -------------------------------------------------------------------------------------------------

# Read in input data values.
sampleInputs <- read.table(sampleInputData, sep= '\t', header = TRUE, stringsAsFactors = FALSE)


# Trim leading and trailing white space from metadata.
trimLeadingTrailingWhtSpace <- function(x) gsub('^\\s+|\\s+$', '', x)
samples$sample_id   <- sapply(samples$sample_id, trimLeadingTrailingWhtSpace)
samples$trial_id    <- sapply(samples$trial_id, trimLeadingTrailingWhtSpace)
samples$patient_id  <- sapply(samples$patient_id, trimLeadingTrailingWhtSpace)
samples$sampleCollection_date <- sapply(samples$sampleCollection_date, trimLeadingTrailingWhtSpace)
samples$sample_type <- sapply(samples$sample_type, trimLeadingTrailingWhtSpace)
samples$VSP         <- sapply(samples$VSP, trimLeadingTrailingWhtSpace)

# Create a vector of VSP ids with sequencning data.
availableVSPs <- unique(str_extract(list.files(paste0(softwareDir, '/summaries/VSPdata')), 'VSP\\d+'))

# Find missing reports.
d <- unique(bind_rows(lapply(availableVSPs, function(x){
  x <- subset(samples, VSP == x)
  if(nrow(x) != 1) return(tibble())
  if(! file.exists(file.path(softwareDir, 'summaries', 'trials', x$trial_id, paste0(x$patient_id, '.pdf'))) | 
     ! file.exists(file.path(softwareDir, 'summaries', 'sampleSummaries', x$trial_id, paste0(x$patient_id, '.sampleSummary.tsv')))){
    return(select(x, patient_id, trial_id))
  } else {
    return(tibble())
  }
}))) 



# Build missing subject reports.
# -------------------------------------------------------------------------------------------------

if(nrow(d) > 0) {
  d <- dplyr::mutate(d, s = ntile(1:n(), CPUs))
 
  CPUs <- 10 
  cluster <- makeCluster(CPUs)
  clusterExport(cluster, c('samples', 'sampleInputs', 'overWriteSubjectReports', 'softwareDir',
                           'sequenceDataDir', 'VSPdataDir', 'sampleData', 'sampleInputData'))
  
  if (!dir.exists(paste0(softwareDir, '/summaries/patientReportsData'))) 
    dir.create(paste0(softwareDir, '/summaries/patientReportsData'))
  
  invisible(parLapply(cluster, split(d, d$s), function(p) {
  #invisible(lapply(split(d, d$s), function(p){
    library(tidyverse)
    library(Biostrings)
    library(lubridate)
    
    invisible(lapply(split(p, paste(p$patient_id, p$trial_id)), function(r) {
      x <- subset(samples, patient_id == r$patient_id & trial_id == r$trial_id)
      
      dir1 <- file.path(softwareDir, 'summaries/trials', r$trial_id)
      dir2 <- file.path(softwareDir, 'summaries/patientReportsData', r$trial_id)
      dir3 <- file.path(softwareDir, 'summaries/sampleSummaries', r$trial_id)
      
      if (!dir.exists(dir1)) dir.create(dir1)
      
      ### message(file.path(dir1, paste0(r$patient_id, '.pdf')))
      
      if (overWriteSubjectReports == FALSE & file.exists(file.path(dir1, paste0(r$patient_id, '.pdf')))) return()
      
      files <- list.files(paste0(softwareDir, '/summaries/VSPdata'), pattern = paste0(x$VSP, collapse = '|'), full.names = TRUE)
      
      if (length(files) == 0) {
        ### write(paste0('[.] No VSP data files for subject ', x$patient_id[1]), file = 'logs/buid_subject_reports.log', append = TRUE)
        return()
      }
      
      if (!dir.exists(dir2)) dir.create(dir2)
      if (!dir.exists(dir3)) dir.create(dir3)
      
      dat <- lapply(files, function(f) {
        # Here we load VSP data files where there are two types of files, analyses of individual sequencing
        # experiments and composite analyses where multiple experiments where combined.
        
        load(f)
        opt$vsp <- stringr::str_extract(f, 'VSP\\d+')
        opt$seq_sample <-  stringr::str_extract(f, 'VSP\\d+\\-?\\d+?[abm]?')
        
        d <- subset(samples, VSP == opt$vsp)
        opt$sample_id <- d$sample_id[1]
        opt$trial_id <- d$trial_id[1]
        opt$subject <- d$patient_id[1]
        opt$date <- d$sampleCollection_date[1]
        opt$sampleType = d$sample_type[1]
        opt$genomesPerMicroLiter = as.numeric(d$N1_copies_per_ul[1])
        
        if (!grepl('composite', f)) {
          opt$inputGenomes <-
            opt$genomesPerMicroLiter * subset(sampleInputs, Sample_Name == opt$seq_sample)$uL_Inputs
        } else {
          opt$inputGenomes <- NA
        }
        opt$concensusSeq = as.character(opt$concensusSeq)
        return(opt)
      })
      
      names(dat) <- unlist(lapply(dat, '[[', 'seq_sample'))
      
      save(dat, file = file.path(dir2, paste0(r$patient_id, '.RData')))
      
      result = tryCatch({
        rmarkdown::render(file.path(softwareDir, 'report.Rmd'),
          output_file = file.path(dir1, paste0(r$patient_id, '.pdf')),
          params = list('date'  = format(Sys.time(), "%Y-%m-%d"), 'title' = paste0('COVID-19 subject ', r$patient_id)))
      }, error = function(e) {
        write(paste0('[!] Failed to create subject report for ', r$patient_id, ' from the ', r$trial_id, ' trial.'),
          file = file.path(softwareDir, 'logs/buid_subject_reports.log'), append = TRUE)
      })
    }))
  }))
  
  stopCluster(cluster)
}



# Calcualte run statistics by collating tables created during the creation of paient reports.
#---------------------------------------------------------------------------------------------

f <- list.files(sequenceDataDir, pattern = '^VSP', recursive = TRUE, full.names = TRUE)
d <- tibble(path = f)

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
d$sampleDate <- samples[match(d$VSP, samples$VSP),]$sampleCollection_date
d$run <- unlist(lapply(strsplit(sub(paste0(sequenceDataDir, '/'), '', d$path), '/'), '[', 1))

d <- arrange(d, desc(days)) %>% select(date, run, exp, trial, subject, sampleDate) %>% distinct()

s <- bind_rows(lapply(list.files(paste0(softwareDir, '/summaries/sampleSummaries'), pattern = '*.tsv$', recursive = TRUE, full.names = TRUE), function(x){
  o <- read.table(x, sep = '\t', header = TRUE)
  o$Subject <- as.character(o$Subject)
  o
}))

# Add read coverage statistics and lineages from sample summaries created with subject reports.
d$percent_1xCoverage <- s[match(d$exp, s$exp),]$percentRefReadCoverage
d$percent_5xCoverage <- s[match(d$exp, s$exp),]$percentRefReadCoverage5
d$lineage <- s[match(d$exp, s$exp),]$lineage

# Only report lineages for samples with >= 95% ref genome coverage w/ >= 5 reads per position.
n <- as.numeric(stringr::str_extract(d$percent_5xCoverage, '[\\d\\.]+'))
n[is.na(n)] <- 0
d[n < 95,]$lineage <- NA

# Write out sequencing run report.
openxlsx::write.xlsx(d, file = file.path(softwareDir, 'summaries', 'seqRunSummary.xlsx'))


# Tally genomes for reporting. Exclude sample ids starting with VSP9 because samples 
# in this upper range are used for controls and simulated data sets.
d$VSP <- stringr::str_extract(d$exp, 'VSP\\d+')
d <- d[! grepl('^VSP9', d$VSP),]
totalSequencedSamples <- n_distinct(d$VSP)

i <- as.numeric(sub('%', '', d$percent_5xCoverage)) >= 95
totalSequencedSamples_highQual <- n_distinct(d[i,]$VSP)



# Collate sample summary table and omit problematic samples and subjects.
# -------------------------------------------------------------------------------------------------

summary <- removeProblematicSamples(bind_rows(lapply(list.files(file.path(softwareDir, 'summaries', 'sampleSummaries'), full.names = TRUE, pattern = '.tsv$', recursive = TRUE), function(x){
  t <- read.table(x, sep = '\t', header = TRUE)
  t$trial_id <- samples[match(sub('\\-\\d+[a-z]?$', '', t$exp, perl = TRUE), samples$VSP), ]$trial_id
  t$Subject <- as.character(t$Subject)
  t$largestContig <- as.numeric(t$largestContig)
  t$percentRefReadCoverage  <- as.numeric(sub('%', '', t$percentRefReadCoverage))
  t$percentRefReadCoverage5 <- as.numeric(sub('%', '', t$percentRefReadCoverage5))
  t$sampleDate2 <- gsub('-', '', t$sampleDate)
  t
}))) %>% filter(! grepl('none|simulated', trial_id, ignore.case = TRUE))


# Identify the best data object for each subject -- typically the composite object.
representativeSampleSummary_95 <- representativeSampleSummary(summary, 95)
concensusSeqs95_5 <- retrieveConcensusSeqs(representativeSampleSummary_95)


genomeMetaData <- bind_rows(lapply(names(concensusSeqs95_5), function(x){
  x <- unlist(strsplit(x, '\\|'))
  tibble(trial         = x[1],
         patient_id    = x[2],
         sample_type   = x[3],
         sample_date   = as.character(ymd(x[4])),
         lab_id        = x[5],
         mean_coverage = as.integer(stringr::str_extract(x[6], '\\d+')),
         lineage       = x[7],
         genome_id     = paste0('h-CoV-19/USA/', x[5], '/', year(ymd(x[4]))))
})) %>% filter(! grepl('^VSP9', lab_id))

genomeMetaData$rationale <- samples[match(genomeMetaData$lab_id, samples$VSP),]$Rationale
openxlsx::write.xlsx(genomeMetaData, file = file.path(softwareDir, 'summaries/highQualGenomes/genomeMetaData.xlsx'))


# Extract lineages from genome ids and write them out to an Excel file.
lineages <- bind_rows(lapply(names(concensusSeqs95_5), function(x){
  x <- unlist(strsplit(x, '\\|'))
  tibble(sample = paste(x[1:3], collapse = '|'), VSP = x[5],  date = x[4], lineage = x[7])
}))

openxlsx::write.xlsx(lineages, file = file.path(softwareDir, 'summaries/highQualGenomes/lineages.xlsx'))




# Create a longitudinal bar plot of identified lineages based on sample arrival date.
# -------------------------------------------------------------------------------------------------

d <- bind_rows(lapply(1:nrow(lineages), function(x){
       x <- lineages[x,]
       o <- ymd(x$date)
       x$days <- as.integer(o)
       x$dateLabel <- paste0(month(o), '/', year(o))
       x
    })) %>%
    dplyr::filter(! is.na(days)) %>%
    dplyr::arrange(days) %>%
    dplyr::mutate(dateLabel = factor(dateLabel, levels = unique(dateLabel)))
  

# d <- d[!is.na(d$days),]
# d <- d[order(d$days),]
# d$dateLabel <- factor(d$dateLabel, levels = unique(d$dateLabel))
# ### d$lineage <- factor(d$lineage, levels = unique(d$lineage))

n <- 19
o <- sort(table(subset(d, ! lineage %in% VOIs)$lineage), decreasing = TRUE)
a <- c(VOIs, names(o)[1:(n-length(VOIs))])
d$lineage <- as.character(d$lineage)
d <- d[! d$lineage == 'NA',]
d$lineage <- ifelse(d$lineage %in% a, d$lineage, 'Other')
d$lineage <- factor(d$lineage, levels = c(a, 'Other'))
colors <- c(grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n), 'gray80')

lineagesPlot <- 
  ggplot(d, aes(dateLabel, fill = lineage)) + 
  theme_bw() +
  geom_bar(stat = 'count') +
  scale_fill_manual(name = 'Lineage', values = colors, drop = FALSE) +
  guides(fill=guide_legend(ncol=2)) +
  labs(x = 'Sample Date', y = 'Genomes') +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16), legend.title=element_text(size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(lineagesPlot, file = file.path(softwareDir, 'summaries/highQualGenomes/lineagesPlot.pdf'), units = 'in', width = 10, height = 7)



# Create variant schematic images and HTML table.
# -------------------------------------------------------------------------------------------------

genomeTable <- bind_rows(lapply(split(representativeSampleSummary_95, 1:nrow(representativeSampleSummary_95)), function(x){
  f <- file.path(softwareDir, 'summaries', 'VSPdata', paste0(x$exp, ifelse(grepl('-', x$exp), '.experiment.RData', '.composite.RData')))
  if(file.exists(f)){
    message('Creating schematic image for: ', f)
    return(variantSchematic(f))
  } else {
    return(tibble()) 
  }
})) %>% dplyr::filter(! grepl('^VSP9', Sample)) %>%  dplyr::arrange(desc(Date))




# Lineage mountain plot.
# Calculate the lineage plot data to determine which lineages should not be considered Other.
# -------------------------------------------------------------------------------------------------

dayBreaks <- '8 days'

d <- genomeMetaData[grepl('random|asymptomatic|hospitalized', 
                          genomeMetaData$rationale, ignore.case = T),] %>%
  dplyr::filter(lineage != 'NA') %>%
  dplyr::filter(ymd(sample_date) >= '2021-02-28') %>%
  dplyr::mutate(date = cut.Date(ymd(sample_date), breaks = dayBreaks)) %>%
  dplyr::group_by(date) %>%
  dplyr::mutate(nGenomes = n_distinct(genome_id)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(lineage = factor(lineage)) %>%
  dplyr::group_by(date, lineage, .drop = FALSE) %>%
  dplyr::summarise(pLineage = n_distinct(genome_id) / nGenomes[1]) %>%
  dplyr::ungroup() %>%
  tidyr::replace_na(list(pLineage = 0)) %>%
  dplyr::arrange(desc(pLineage)) 

# Find the top 8 lineages not in variant of interest list.
o <- unique(as.character(d$lineage))
o <- c(VOIs, o[! o %in% VOIs][1:8])

colors <- c(grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(o)-1), 'gray80')

lineageMountainPlot <- genomeMetaData[grepl('random|asymptomatic|hospitalized', 
                          genomeMetaData$rationale, ignore.case = T),] %>%
     dplyr::filter(lineage != 'NA') %>%
     dplyr::filter(ymd(sample_date) >= '2021-02-28') %>%
     dplyr::arrange(sample_date) %>%
     dplyr::mutate(date = cut.Date(ymd(sample_date), breaks = dayBreaks)) %>%
     dplyr::group_by(date) %>%
     dplyr::mutate(nGenomes = n_distinct(genome_id)) %>%
     dplyr::ungroup() %>%
     dplyr::mutate(lineage = forcats::fct_other(lineage, keep = o)) %>%
     dplyr::group_by(date, lineage, .drop = FALSE) %>%
     dplyr::summarise(pLineage = n_distinct(genome_id) / nGenomes[1], nGenomes = nGenomes[1]) %>%
     dplyr::ungroup() %>%
     tidyr::replace_na(list(pLineage = 0, nGenomes = 0)) %>%
     dplyr::mutate(date2 = ymd(as.character(date))) %>%
     dplyr::group_by(date) %>%
     dplyr::arrange(desc(nGenomes)) %>%
     dplyr::mutate(date2 = paste0(months(date2), ' ', day(date2), '\n', nGenomes[1], ' genomes')) %>%
     dplyr::ungroup() %>%
     dplyr::arrange(date) %>%
     dplyr::mutate(date2 = factor(as.character(date2), levels = unique(date2))) %>%
     ggplot(aes(date2, pLineage, fill = lineage, group = lineage)) +
       scale_fill_manual(values = colors) +
       geom_area(color = "black", size = 0.2) +
       scale_x_discrete(expand = c(0.01,0)) +
       scale_y_continuous(expand = c(0.01,0), labels = scales::percent_format(accuracy = 1)) +
       labs(x = '', y = 'Percent Sampled Lineages') +
       guides(fill=guide_legend(title="Lineages")) +
       theme(axis.text = element_text(size = 12), axis.text.x = element_text(angle = 0, vjust=0),
             panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
             panel.background = element_blank(), axis.line = element_line(colour = "black"))



d <- dplyr::filter(genomeMetaData, lineage != 'NA') %>% dplyr::arrange(sample_date)

# Add a faux row at the start of the data frame to ensure that the date cut factor
# begins on the first of the month.

r <- d[1,]
m <- ymd(r$sample_date)
r$sample_date <- as.character(ymd(paste0(c(year(m), month(m), 1), collapse = '-')))
r$lineage <- 'faux'

d <- bind_rows(r, d) %>%
     dplyr::mutate(dateCut = cut.Date(ymd(sample_date), breaks = dayBreaks)) %>%
     dplyr::filter(lineage != 'faux') %>%
     dplyr::group_by(dateCut) %>%
  dplyr::mutate(nGenomes = n_distinct(genome_id)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(dateCut, lineage, .drop = FALSE) %>%
  dplyr::summarise(pLineage = n_distinct(lab_id) / nGenomes[1], nGenomes = nGenomes[1]) %>%
  dplyr::ungroup() %>%
  tidyr::replace_na(list(pLineage = 0, nGenomes = 0, lineage = 'none')) %>%
  dplyr::mutate(lineage2 = forcats::fct_other(lineage, keep = o)) %>%
  dplyr::mutate(dateCut2 = paste0(month(ymd(as.character(dateCut)), label = TRUE), ' ',
                                  paste0(year(ymd(as.character(dateCut))))))

sampledLineages1 <- ggplot(d, aes(dateCut, pLineage, fill=lineage2)) +
  scale_fill_manual(values = c(grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(o)-1), 'gray80')) +
  geom_col() +
  scale_y_continuous(expand = c(0, 0), labels = scales::percent_format(accuracy = 1)) +
  scale_x_discrete(expand = c(0, 0), 
                   breaks = d$dateCut[! duplicated(d$dateCut2)], 
                   labels = d$dateCut2[! duplicated(d$dateCut2)]) +
  labs(x = '', y = 'Percent Sampled Lineages') +
  guides(fill=guide_legend(title="Lineages")) +
  theme(axis.text = element_text(size = 12), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


sampledLineages2 <- ggplot(dplyr::select(d, dateCut, nGenomes) %>% dplyr::distinct(), aes(dateCut, nGenomes)) +
  geom_col() +
  scale_y_continuous(expand = c(0, 0), trans = 'log2', breaks = c(1,4,16,64,256), labels = c(0,4,16,64,256)) +
  scale_x_discrete(expand = c(0, 0), 
                   breaks = d$dateCut[! duplicated(d$dateCut2)], 
                   labels = d$dateCut2[! duplicated(d$dateCut2)]) +
  labs(x = '', y = expression(Log[2]~Genomes))+
  theme(axis.text = element_text(size = 12), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=14),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

rmarkdown::render('dashboard.Rmd', output_file = 'index.html',  params = list('lastUpdated' = date()))




# Write out high quality genome sequences.
writeXStringSet(concensusSeqs95_5, file = file.path(softwareDir, 'summaries/highQualGenomes/genomes.fasta'))


# Rename the reference genome because by default it is named 'genome' which is needed for assemble.R.
r <- Biostrings::readDNAStringSet(file.path(softwareDir, 'data/references/Wuhan-Hu-1.fasta'))
names(r) <- 'Wuhan-Hu-1'
Biostrings::writeXStringSet(r, file.path(softwareDir, 'summaries/highQualGenomes/referenceGenome.fasta'))


# Shorten sequence ids for trees.
o <- Biostrings::readDNAStringSet(file.path(softwareDir, 'summaries/highQualGenomes/genomes.fasta'))
names(o) <- sapply(names(o), function(x){ paste0(unlist(strsplit(x, '\\|'))[c(1:2,5)], collapse = '|') })
names(o) <- gsub('Saliva_-_Positive_Control', 'Saliva', names(o))
names(o) <- gsub('PennEssentialWorkers_Dec2020', 'PennEssential', names(o))
#names(o) <- gsub('\\|[A-Z\\.\\d]+$', '', names(o), perl = TRUE)
Biostrings::writeXStringSet(o, file.path(softwareDir, 'summaries/highQualGenomes/genomes.fasta2'))


# Hierarchical tree.
system(paste0(mafftPath, ' --namelength ', max(nchar(names(concensusSeqs95_5))), ' --thread 20 --auto --addfragments ',
       file.path(softwareDir, 'summaries/highQualGenomes/genomes.fasta2'), ' ',
       file.path(softwareDir, 'summaries/highQualGenomes/referenceGenome.fasta'), ' > ',
       file.path(softwareDir, 'summaries/highQualGenomes/genomes.mafft')))

v <- ape::read.dna(file.path(softwareDir, "summaries/highQualGenomes/genomes.mafft"), format="fasta")
v_phyDat <- phangorn::phyDat(v, type = "DNA", levels = NULL)
dna_dist <- phangorn::dist.ml(v_phyDat, model="JC69")

dendr    <- ggdendro::dendro_data(hclust(dna_dist, method='average'), type="rectangle") 
segments <- ggdendro::segment(dendr)
labels   <- ggdendro::label(dendr)

concensusSeqPhyloPlot <- 
  ggplot() + 
  geom_segment(data=segments, aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=labels, aes(x=x, y=y, label=label, hjust=0), size=4) +
  coord_flip() + 
  scale_y_reverse(expand=c(0.1, 0)) + 
  labs(x = '', y = 'Distance') +
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

ggsave(concensusSeqPhyloPlot, height = 150, width = 50, units = 'in', limitsize = FALSE, file = file.path(softwareDir, 'summaries/highQualGenomes/hierarchicalPhyloPlot.pdf'))
invisible(file.remove(list.files(file.path(softwareDir, 'summaries/highQualGenomes'), pattern = 'mafft|referenceGenome|fasta2', full.names = TRUE)))


# Write out genome sequences with GISAID style ids.
names(concensusSeqs95_5) <- genomeMetaData$genome_id
writeXStringSet(concensusSeqs95_5, file = file.path(softwareDir, 'summaries/highQualGenomes/genomes.fasta'))


# Genome id check
uniqueVSPids <- unique(sapply(names(concensusSeqs95_5), function(x){ unlist(strsplit(x, '/'))[3] }))
message(length(uniqueVSPids), ' high quality genomes exported.')


# Variant tables.
d <- bind_rows(lapply(list.files(file.path(softwareDir, 'summaries', 'sampleSummaries'), full.names = TRUE, pattern = '.tsv$', recursive = TRUE), function(y){
  t <- read.table(y, sep = '\t', header = TRUE)
  t$Subject <- as.character(t$Subject)
  t$percentRefReadCoverage5 <- as.numeric(sub('%', '', t$percentRefReadCoverage5))
  t$sampleDate2 <- gsub('-', '', t$sampleDate)
  t
})) %>% filter(percentRefReadCoverage5 >= 95 & ! grepl('simulate', sampleType, ignore.case = TRUE))

o <- representativeSampleSummary(d, minPercentRefReadCoverage5 = 95) %>% pull(exp)

r <- bind_rows(lapply(o, function(e){
  f <- file.path(softwareDir, 'summaries', 'VSPdata', paste0(e, ifelse(grepl('-', e), '.experiment.RData', '.composite.RData')))
  if(! file.exists(f)){
    message('Error - could not locate ', f)
    return(tibble())
  }
  load(f)
  if(nrow(opt$variantTableMajor) == 0) return(tibble())
  tibble(experiment = e, 
         position = opt$variantTableMajor$POS,
         variant = paste0(opt$variantTableMajor$REF, opt$variantTableMajor$POS, opt$variantTableMajor$ALT),
         gene = opt$variantTableMajor$genes,
         type = opt$variantTableMajor$type)
}))


r2 <- bind_rows(lapply(split(r, r$position), function(e){
  tibble(position = e$position[1],
         variants = paste0(unique(e$variant), collapse = ', '),
         gene = paste0(unique(e$gene), collapse = ', '),
         types = paste0(unique(e$type), collapse = ', '),
         experiments = n_distinct(e$experiment))
})) %>% arrange(desc(experiments)) 

openxlsx::write.xlsx(r2, file.path(softwareDir, 'summaries/highQualGenomes/SNPsummary.xlsx'))


#--------------------------------------------------------------------------

file.remove(file.path(softwareDir, '/summaries/SARS-CoV-2_reports.zip'))
system(paste0('zip -9 -r ', softwareDir, '/summaries/SARS-CoV-2_reports.zip ', softwareDir, '/summaries/trials/'))


system(paste0('rsync -r --del ', file.path(softwareDir, 'summaries'), 
              ' microb120:/media/lorax/data/SARS-CoV-2/'))

system(paste0('scp ', file.path(softwareDir, 'index.html'), 
              ' microb120:/media/lorax/data/SARS-CoV-2/index.html'))

system('ssh microb120 rm /media/lorax/data/SARS-CoV-2/working')
