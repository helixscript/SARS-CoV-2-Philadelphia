library(parallel)
library(dplyr)
source('lib/lib.R')

CPUs <- 30
overwrite <- FALSE
Rscript <- '/home/opt/R-3.4.0/bin/Rscript'
sequenceDataDir <- '/data/sequencing/Illumina-archive/SARS-CoV-2'
VSPdataDir <- '/home/everett/projects/SARS-CoV-2-Philadelphia/summaries/VSPdata'
sampleData <- '/data/sequencing/Illumina-archive/SARS-CoV-2/samples.tsv'


# Clear previous log files.
invisible(file.remove(list.files('logs', pattern = 'build_VSP_data', full.names = TRUE)))


# Read in sample data table.
samples <- removeProblematicSamples(read.table(sampleData, sep= '\t', header = TRUE, quote = '', stringsAsFactors = FALSE))



cluster <- makeCluster(CPUs)
clusterExport(cluster, c('Rscript', 'overwrite', 'sequenceDataDir'))

# Create a list of all VSP R1 sequencng files 
f1 <- unlist(lapply(strsplit(list.files(sequenceDataDir, recursive = TRUE, pattern = '^VSP', full.names = FALSE), '/'), '[[', 2))
f1 <- unique(unlist(lapply(strsplit(f1, '_'), '[', 1)))

# Create a list of all VSP data objects already created.
f2 <- list.files(VSPdataDir, recursive = TRUE, pattern = '^VSP', full.names = FALSE)
f2 <- unique(unlist(lapply(strsplit(f2, '\\.'), '[', 1)))

# Identify new incoming experiments. 
f3 <- f1[! f1 %in% f2]

VSPs <- unique(sub('\\-\\d+[a-z]?', '', f3))

# Remove patient reports and sample summaries so that they are recreated.
invisible(sapply(VSPs, function(x){
  x <- subset(samples, VSP == x)
  if(nrow(x) == 1){
    if(file.exists(file.path('summaries', 'trials', x$trial_id, paste0(x$patient_id, '.pdf')))){
      file.remove(file.path('summaries', 'trials', x$trial_id, paste0(x$patient_id, '.pdf')))
      
      if(file.exists(file.path('summaries', 'sampleSummaries', x$trial_id, paste0(x$patient_id, '.sampleSummary.tsv')))){
        file.remove(file.path('summaries', 'sampleSummaries', x$trial_id, paste0(x$patient_id, '.sampleSummary.tsv')))
      }
    }
  }
}))


d <- data.frame(vsp = VSPs)
d$s <- ntile(1:nrow(d), CPUs)


invisible(parLapply(cluster, split(d, d$s), function(x){
#invisible(lapply(split(d, d$s), function(x){
  library(dplyr)
  library(stringr)
  
  logFile <- paste0('logs/build_VSP_data.', x$s[1], '.log')
  write(date(), file = logFile, append = FALSE)
  
  invisible(sapply(x$vsp, function(vsp){
    #message('Starting ', vsp, '\n\n')
    
    # All the files to draw R1 and R2 read files from for VSP.
    sampleFiles <- list.files(sequenceDataDir, recursive = TRUE, pattern = paste0('^', vsp, '-'))
  
    # Split experiment files into a list where each key contains the R1 and R2 sequencing files.
    s <- unlist(lapply(strsplit(unlist(lapply(strsplit(sampleFiles, '/'), '[', 2)), '_'), '[', 1))
    
    sampleRunFiles <- split(sampleFiles, s)

    write(paste0('Starting ', vsp, ' with ', length(sampleRunFiles), ' seq experiments.'), file = logFile, append = TRUE)
    
    if(length(sampleRunFiles) > 1){
      outputFile <- paste0('summaries/VSPdata/', vsp, '.composite.RData')
      
      if(overwrite == FALSE & file.exists(outputFile)){
           write(paste0(outputFile , ' exists -- skipping.'), file = logFile, append = TRUE)
      } else {
        if(file.exists(outputFile)) file.remove(outputFile)
        
        comm <- paste0(Rscript, ' assemble.R --outputFile ', outputFile, ' --R1 ', 
                       paste0(sequenceDataDir, '/', sampleFiles[grepl('_R1_', sampleFiles)], collapse = ','), ' --R2 ',
                       paste0(sequenceDataDir, '/', sampleFiles[grepl('_R2_', sampleFiles)], collapse = ','))
                  
        write(paste0(date(), ' (', outputFile, ') ', comm, '\n'), file = logFile, append = TRUE)
        
        tryCatch({
           system(comm)
        }, error = function(e) {
           write(paste0('\n[!] (', outputFile, ') data object creation failed, system error caught: "', e, '".\n'), file = logFile, append = TRUE)
        })
        
        if(! file.exists(outputFile)) write(paste0('\n[!] (', outputFile, ') data object creation failed.'), file = logFile, append = TRUE)
      }
    }
    
    invisible(mapply(function(n, x){
      outputFile <- paste0('summaries/VSPdata/', n, '.experiment.RData')
         
       if(overwrite == FALSE & file.exists(outputFile)){
         write(paste0(outputFile , ' exists -- skipping.'), file = logFile, append = TRUE)
       } else {
         comm <-  paste0(Rscript, ' assemble.R --outputFile ', outputFile, ' --R1 ', 
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
if(file.exists('logs/build_VSP_data.failed')) file.remove('logs/build_VSP_data.failed')
system("find ./logs -name build_VSP_data* | xargs grep failed > logs/build_VSP_data.failed")

