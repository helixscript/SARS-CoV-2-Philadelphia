library(stringr)
library(dplyr)
f <- unlist(list.files('sequencing', recursive = TRUE, pattern = '*.gz$'))
f <- f[grepl('\\d[ab]_', f)]
d <- tibble(path = f, exp = str_extract(f, 'VSP\\d+\\-\\d'), direction = str_match(f, '_(R\\d)_0')[,2])

r <- bind_rows(lapply(split(d, paste( d$exp, d$direction)), function(x){
       if(as.integer(str_extract(x$exp[1], '\\d+')) >= 9000) return(tibble())
       if(nrow(x) == 1) return(tibble())
  
       d <- file.path('sequencingBackup', str_extract(x$path[1], '[^/]+'))
       if(! dir.exists(d)) dir.create(d)
        
       x$path <- paste0('sequencing/', x$path)
       system(paste('cp', x[1,]$path, x[2,]$path, paste0(d, '/')))
        
       comm <- paste('cat', x[1,]$path, x[2,]$path, '>', sub(paste0(x[1,]$exp, '[ab]'), paste0(x[1,]$exp, 'm'), paste0(x[1,]$path)))
       system(comm)
        
       system(paste('rm', x[1,]$path, x[2,]$path))
        
       tibble(exp = x[1,]$exp, samples = nrow(x), comm = comm)
     }))
