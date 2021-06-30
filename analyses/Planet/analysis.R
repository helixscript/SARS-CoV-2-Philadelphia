library(dplyr)

# Read in the Planet VCF files and extract the alternative base frequencies.
d <- bind_rows(lapply(list.files('Planet_vcf', full.names = TRUE), function(x){
  o <- read.table(x, sep = '\t')
  bind_rows(lapply(split(o, 1:nrow(o)), function(x2){
    a <- unlist(strsplit(x2$V9, ':'))
    b <- unlist(strsplit(x2$V10, ':'))
    data.frame(Planet_subject = unlist(strsplit(unlist(strsplit(x, '/'))[2], '_'))[1],
               Planet_pos = x2$V2, Planet_depth = as.integer(stringr::str_extract(x2$V8, '\\d+')),
               Planet_altFreq = as.numeric(b[which(a == "ALT_FREQ")]))
  }))
}))

s <- read.table('samples.tsv', sep = '\t', header = TRUE)

o <- bind_rows(lapply(split(s, 1:nrow(s)), function(x){
  x$Planet_subject <- gsub('_', '-', x$Planet_subject)
  s <- subset(d, Planet_subject == x$Planet_subject)
  load(paste0('../../summaries/VSPdata/', x$VSP, '-1.experiment.RData'))
  left_join(opt$variantTableMajor, subset(d, Planet_subject == x$Planet_subject), by = c('POS' = 'Planet_pos')) %>%
            select(-Planet_subject, -CHROM) %>%
            mutate(altFreqDiff = abs(percentAlt - Planet_altFreq)) %>%
            tibble::add_column(subject = x$Planet_subject[1], .before = 'POS')
}))
  
openxlsx::write.xlsx(o, file = 'Bushman_Planet_variant_table.xlsx')
