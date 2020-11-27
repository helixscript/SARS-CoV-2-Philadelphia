options(stringsAsFactors = FALSE)
o <- Sys.getenv("PATH")
Sys.setenv(PATH = paste('/home/everett/miniconda3/bin:/home/everett/hisatgenotype:/home/everett/samtools-1.3.1/bin:/home/everett/hisatgenotype/hisat2',o, sep = ":"))
Sys.setenv(PYTHONPATH = '/home/everett/hisatgenotype/hisatgenotype_modules')

f <- list.files('trimmedSeqDataCollations')
d <- data.frame(file = f[!grepl('Control', f)])
d$sample <- unlist(lapply(strsplit(d$file, '\\.'), '[', 1)) 

library(parallel)
cluster <- makeCluster(15)


invisible(parLapply(cluster, split(d, d$sample), function(x){
  comm <- paste0('hisatgenotype --base hla -1 trimmedSeqDataCollations/',
                 x$file[grepl('.R1.', x$file)], ' -2 trimmedSeqDataCollations/', x$file[grepl('.R2.', x$file)],
                 ' -p 3 --out-dir hisatgenotype/', x$sample[1])
  dir.create(paste0('hisatgenotype/', x$sample[1]))
  system(comm)
}))