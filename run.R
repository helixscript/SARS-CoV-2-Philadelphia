Rscript <- '/home/opt/R-3.4.0/bin/Rscript'
system(paste0(Rscript, ' build_VSP_data.R'))
system(paste0(Rscript, ' build_reports.R'))
system(paste0(Rscript, ' build_consensus_genomes.R'))
system(paste0(Rscript, ' build_variant_tables.R'))
system(paste0(Rscript, ' build_trees.R'))