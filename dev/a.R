library(dplyr)

r <- bind_rows(lapply(readLines('sampleList'), function(x){
  load(paste0('../summaries/VSPdata/', x, '.experiment.RData'))
  if(nrow(opt$variantTableMajor) > 0) opt$variantTableMajor$CHROM <- x
  opt$variantTableMajor
}))

lapply(readLines('sampleList'), function(x){
  file.remove(paste0('../summaries/VSPdata/', x, '.experiment.RData'))
  x2 <- sub('\\-\\d', '', x)
  file.remove(paste0('../summaries/VSPdata/', x2, '.composite.RData'))
})



