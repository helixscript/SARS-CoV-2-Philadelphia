library(tidyverse)
source('lib/lib.R')
mafftPath <- '~/ext/mafft/bin/mafft'

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


# Build phylogenetic tree.
system(paste0(mafftPath, '  --thread 10 --auto --addfragments summaries/allGenomes_90_5.fasta data/references/USA-WA1-2020.fasta> summaries/allGenomes_90_5.mafft'))


v <- ape::read.dna("summaries/allGenomes_90_5.mafft", format="fasta")
invisible(file.remove('summaries/allGenomes_90_5.mafft'))
v_phyDat <- phangorn::phyDat(v, type = "DNA", levels = NULL)
dna_dist <- phangorn::dist.ml(v_phyDat, model="JC69")

dendr    <- ggdendro::dendro_data(hclust(dna_dist, method='average'), type="rectangle") 
segments <- ggdendro::segment(dendr)
labels   <- ggdendro::label(dendr)

concensusSeqPhyloPlot <- 
  ggplot() + 
  geom_segment(data=segments, aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=labels, aes(x=x, y=y, label=label, hjust=0), size=3) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

ggsave(concensusSeqPhyloPlot, height = 10, width = 12, units = 'in', file = 'summaries/allGenomes_90_5_hierarchicalPhyloPlot.pdf')


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


