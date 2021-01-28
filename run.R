Rscript <- '/home/opt/R-3.4.0/bin/Rscript'
system(paste0(Rscript, ' build_VSP_data.R'))
system(paste0(Rscript, ' build_reports.R'))
system(paste0(Rscript, ' build_trial_tables.R'))
