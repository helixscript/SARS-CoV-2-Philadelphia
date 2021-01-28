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




# Collate sampe summary table and omit problematic samples and subjects.
summary <- removeProblematicSamples(bind_rows(lapply(list.files('summaries/sampleSummaries', full.names = TRUE), function(x){
  t <- read.table(x, sep = '\t', header = TRUE)
  t$Subject <- as.character(t$Subject)
  t$largestContig <- as.numeric(t$largestContig)
  t$percentRefReadCoverage  <- as.numeric(sub('%', '', t$percentRefReadCoverage))
  t$percentRefReadCoverage5 <- as.numeric(sub('%', '', t$percentRefReadCoverage5))
  t$sampleDate2 <- gsub('-', '', t$sampleDate)
  t
})))


# Identify the best data object for each subject -- typically the composite object.
representativeSampleSummary_90 <- representativeSampleSummary(summary, 90)


# Retrieve VSP objects and create a concatnetated table of variants.
retrieveVariantTables <- function(summary){
  bind_rows(lapply(1:nrow(summary), function(x){
    x <- summary[x,]
    f <- paste0('summaries/VSPdata/', x$exp, '.', ifelse(x$type == 'composite', 'composite', 'experiment'), '.RData')
    if(! file.exists(f)) stop('Can not locate VSP data file -- ', f)
    load(f)
    if(nrow(opt$variantTableMajor) > 0) opt$variantTableMajor$sample <-  paste0(x$Subject, '|', x$sampleType, '|', x$sampleDate2)
    opt$variantTableMajor
  }))
}

variantTable <- retrieveVariantTables(representativeSampleSummary_90)



# Create a table of variants with subject and sample counts.
variantTable$subject <- unlist(lapply(str_extract_all(variantTable$sample, '([^\\|]+)'), '[', 1))
variantSummary <- 
  group_by(variantTable, POS) %>%
  summarise(nSubjects = n_distinct(subject), nSamples = n_distinct(sample), gene = genes[1], mutation = type[1]) %>%
  ungroup() %>% 
  arrange(desc(nSamples))

openxlsx::write.xlsx(variantSummary, file = 'summaries/allGenomes_90_5_variantSummary.xlsx')

v <- group_by(variantTable, subject, sample) %>% mutate(variants = paste0(POS, '|', genes)) %>% ungroup() %>% select(subject, sample, variants)
v2 <- reshape2::dcast(mutate(v, x = 'x'), subject+sample~variants, value.var = 'x')
openxlsx::write.xlsx(v2, file = 'summaries/allGenomes_90_5_variantSummary2.xlsx')



# Order variant table to match dendro plot
variantTable$sample <- factor(as.character(variantTable$sample), levels = unique(dendr$labels$label))
variantTable <- dplyr::arrange(variantTable, POS)
variantTable$genes <- factor(as.character(variantTable$genes), levels = unique(as.character(variantTable$genes)))

variantTable$genes <- as.character(variantTable$genes)
variantTable$genes <- gsub('intergenic', 'Intergenic', variantTable$genes)
variantTable$genes <- factor(as.character(variantTable$genes), levels = unique(as.character(variantTable$genes)))

genomeVariantPlot <- 
  ggplot(variantTable, aes(POS, sample, color = genes)) + 
  theme_bw() +
  labs(x = 'Genome position', y = '') +
  scale_color_manual(name = 'Gene', 
                     values = c('gray80', grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_distinct(variantTable$genes) -1))) +
  geom_point(size = 4) +
  scale_y_discrete(drop=FALSE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.position="bottom") +
  guides(color=guide_legend(nrow=1))

ggsave(genomeVariantPlot, height = 10, width = 12, units = 'in', file = 'summaries/allGenomes_90_5_genomeVariantPlot.pdf', useDingbats = FALSE)


