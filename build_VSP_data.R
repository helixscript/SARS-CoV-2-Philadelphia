library(parallel)
library(dplyr)

CPUs <- 30
overwrite <- FALSE
Rscript <- '/home/opt/R-3.4.0/bin/Rscript'

cluster <- makeCluster(CPUs)
clusterExport(cluster, c('Rscript', 'overwrite'))

f <- unlist(lapply(strsplit(list.files('data/sequencing', recursive = TRUE, pattern = 'VSP.+_R1_', full.names = FALSE), '/'), '[[', 2))

d <- data.frame(vsp = unique(unlist(lapply(strsplit(unlist(lapply(strsplit(f, '_'), '[[', 1)), '-'), '[[', 1))))
d$s <- ntile(1:nrow(d), CPUs)


invisible(parLapply(cluster, split(d, d$s), function(x){
#invisible(lapply(split(d, d$s), function(x){
  library(dplyr)
  library(stringr)
  
  logFile <- paste0('logs/build_VSP_data.', x$s[1], '.log')
  write(date(), file = logFile, append = FALSE)
  
  invisible(sapply(x$vsp, function(vsp){

    message(vsp)
    # All the files to draw R1 and R2 read files from for VSP.
    sampleFiles <- list.files('data/sequencing', recursive = TRUE, pattern = paste0('^', vsp, '-'))
  
    # sampleFiles split by sample for the creation of fine grain data objects.
    s <- unlist(lapply(strsplit(unlist(lapply(strsplit(sampleFiles, '/'), '[', 2)), '_'), '[', 1))
    
    sampleRunFiles <- split(sampleFiles, s)

    write(paste0('Starting ', vsp, ' with ', length(sampleRunFiles), ' seq experiments.'), file = logFile, append = TRUE)
    
    if(length(sampleRunFiles) > 1){
      outputFile <- paste0('summaries/VSPdata/', vsp, '.composite.RData')
      
      if(overwrite == FALSE & file.exists(outputFile)){
           write(paste0(outputFile , ' exists -- skipping.'), file = logFile, append = TRUE)
      } else {
        comm <- paste0(Rscript, ' assemble.R --outputFile ', outputFile, ' --R1 ', 
                       paste0('data/sequencing/', sampleFiles[grepl('_R1_', sampleFiles)], collapse = ','), ' --R2 ',
                       paste0('data/sequencing/', sampleFiles[grepl('_R2_', sampleFiles)], collapse = ','))
                  
        write(paste0(date(), ' (', outputFile, ') ', comm, '\n'), file = logFile, append = TRUE)
        
        #browser() 
        
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
                         paste0('data/sequencing/', x[grepl('_R1_', x)], collapse = ','), ' --R2 ',
                         paste0('data/sequencing/', x[grepl('_R2_', x)], collapse = ','))
         
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
system("find ./logs -name 'build_VSP_data*' | xargs grep 'failed' > 'logs/buid_VSP_data.failed'")

