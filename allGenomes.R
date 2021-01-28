library(dplyr)
library(Biostrings)
library(tidyverse)
library(lubridate)
source('lib/lib.R')
mafftPath <- '~/ext/mafft/bin/mafft'


samples <- removeProblematicSamples(read.table('data/samples.tsv', sep= '\t', header = TRUE, quote = '', stringsAsFactors = FALSE))

# Collate sample summary table and omit problematic samples and subjects.
summary <- removeProblematicSamples(bind_rows(lapply(list.files('summaries/sampleSummaries', full.names = TRUE, pattern = '.tsv$', recursive = TRUE), function(x){
  t <- read.table(x, sep = '\t', header = TRUE)
  t$trial_id <- samples[match(sub('\\-\\d+[a-z]?$', '', t$exp, perl = TRUE), samples$VSP), ]$trial_id
  t$Subject <- as.character(t$Subject)
  t$largestContig <- as.numeric(t$largestContig)
  t$percentRefReadCoverage  <- as.numeric(sub('%', '', t$percentRefReadCoverage))
  t$percentRefReadCoverage5 <- as.numeric(sub('%', '', t$percentRefReadCoverage5))
  t$sampleDate2 <- gsub('-', '', t$sampleDate)
  t
}))) %>% filter(! grepl('none|simulated', trial_id, ignore.case = TRUE))


# Identify the best data object for each subject -- typically the composite object.
representativeSampleSummary_90 <- representativeSampleSummary(summary, 90)
concensusSeqs90_5 <- retrieveConcensusSeqs(representativeSampleSummary_90)
writeXStringSet(concensusSeqs90_5, file = 'summaries/highQualGenomes/genomes.fasta')


# Extract lineages from genome ids and write them out to an Excel file.
lineages <- bind_rows(lapply(names(concensusSeqs90_5), function(x){
  x <- unlist(strsplit(x, '\\|'))
  tibble(sample = paste(x[1:4], collapse = '|'), lineage = x[5])
}))

openxlsx::write.xlsx(lineages, file = 'summaries/highQualGenomes/lineages.xlsx')


# Create a longitudinal bar plot of identified lineages.
d <- bind_rows(lapply(1:nrow(lineages), function(x){
       x <- lineages[x,]
       x$date <- unlist(str_split(x$sample, '\\|'))[4]
       o <- ymd(x$date)
       x$days <- as.integer(o)
       x$dateLabel <- paste0(month(o), '/', year(o))
       x
     }))

d <- d[order(d$days),]
d$dateLabel <- factor(d$dateLabel, levels = unique(d$dateLabel))
d$lineage <- factor(d$lineage, levels = unique(d$lineage))

colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_distinct(d$lineage)) 

lineagesPlot <- 
  ggplot(d, aes(dateLabel, fill = lineage)) + 
  theme_bw() +
  geom_bar(stat = 'count') +
  scale_fill_manual(name = 'Lineage', values = colors) +
  guides(fill=guide_legend(ncol=2)) +
  labs(x = 'Sample Date', y = 'Genomes') +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16), legend.title=element_text(size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(lineagesPlot, file = 'summaries/highQualGenomes/lineagesPlot.pdf', units = 'in', width = 10, height = 5)



# Create a raxml tree for all the high quality genomes.

# Rename the refernece genome because by default it is named 'genome' which is needed for assemble.R.
r <- Biostrings::readDNAStringSet('data/references/USA-WA1-2020.fasta')
names(r) <- 'USA-WA1-2020'
Biostrings::writeXStringSet(r, 'summaries/highQualGenomes/referenceGenome.fasta')


# Shorten sequence ids for trees.
o <- Biostrings::readDNAStringSet('summaries/highQualGenomes/genomes.fasta')
names(o) <- gsub('Saliva_-_Positive_Control', 'Saliva', names(o))
names(o) <- gsub('PennEssentialWorkers_Dec2020', 'PennEssential', names(o))
names(o) <- gsub('\\|[A-Z\\.\\d]+$', '', names(o), perl = TRUE)
Biostrings::writeXStringSet(o, 'summaries/highQualGenomes/genomes.fasta2')

system(paste0(mafftPath, ' --phylipout --namelength ', max(nchar(names(concensusSeqs90_5))), ' --thread 10 --auto --addfragments summaries/highQualGenomes/genomes.fasta2 summaries/highQualGenomes/referenceGenome.fasta > summaries/highQualGenomes/genomes.mafft'))
system('raxmlHPC-PTHREADS-SSE3 -s summaries/highQualGenomes/genomes.mafft -m GTRGAMMA -T 25 -n raxmlOut -f a -x 12345 -p 12345 -N autoMRE')

dendr <- ggdendro::dendro_data(phylogram::read.dendrogram('RAxML_bestTree.raxmlOut'), type="rectangle")
segments <- ggdendro::segment(dendr)
labels <- ggdendro::label(dendr)

p <- ggplot() + 
     geom_segment(data=segments, aes(x=x, y=y, xend=xend, yend=yend)) +
     geom_text(data=labels, aes(x=x, y=y, label=label, hjust=0), size=3) +
     coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
     labs(x = '', y = 'Distance') +
     theme(axis.line.y=element_blank(),
           axis.ticks.y=element_blank(),
           axis.text.y=element_blank(),
           axis.title.y=element_blank(),
           panel.background=element_rect(fill="white"),
           panel.grid=element_blank())

invisible(file.remove(list.files(pattern = 'raxmlOut')))

ggsave(p, filename = 'summaries/highQualGenomes/raxmlPhyloPlot.pdf', units = 'in', height = 15, width = 15)


# Hierarchical tree.
system(paste0(mafftPath, ' --namelength ', max(nchar(names(concensusSeqs90_5))), ' --thread 10 --auto --addfragments summaries/highQualGenomes/genomes.fasta2 summaries/highQualGenomes/referenceGenome.fasta > summaries/highQualGenomes/genomes.mafft'))

v <- ape::read.dna("summaries/highQualGenomes/genomes.mafft", format="fasta")
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
  labs(x = '', y = 'Distance') +
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

ggsave(concensusSeqPhyloPlot, height = 12, width = 15, units = 'in', file = 'summaries/highQualGenomes/hierarchicalPhyloPlot.pdf')
invisible(file.remove(list.files('summaries/highQualGenomes', pattern = 'mafft|referenceGenome|fasta2', full.names = TRUE)))


# Variant tables.

d <- bind_rows(lapply(list.files(file.path('summaries', 'sampleSummaries'), full.names = TRUE, pattern = '.tsv$', recursive = TRUE), function(y){
  t <- read.table(y, sep = '\t', header = TRUE)
  t$Subject <- as.character(t$Subject)
  t$percentRefReadCoverage5 <- as.numeric(sub('%', '', t$percentRefReadCoverage5))
  t$sampleDate2 <- gsub('-', '', t$sampleDate)
  t
})) %>% filter(percentRefReadCoverage5 >= 90 & ! grepl('simulate', sampleType, ignore.case = TRUE))

o <- representativeSampleSummary(d, minPercentRefReadCoverage5 = 90) %>% pull(exp)

r <- bind_rows(lapply(o, function(e){
    f <- file.path('summaries', 'VSPdata', paste0(e, ifelse(grepl('-', e), '.experiment.RData', '.composite.RData')))
    if(! file.exists(f)) stop('Error - could not locate ', f)
    load(f)
    if(nrow(opt$variantTableMajor) == 0) return(tibble())
    tibble(experiment = e, 
           position = opt$variantTableMajor$POS,
           variant = paste0(opt$variantTableMajor$REF, opt$variantTableMajor$POS, opt$variantTableMajor$ALT),
           gene = opt$variantTableMajor$genes,
           type = opt$variantTableMajor$type)
  }))
  

r2 <- bind_rows(lapply(split(r, r$position), function(e){
        tibble(position = e$position[1],
               variants = paste0(unique(e$variant), collapse = ', '),
               gene = paste0(unique(e$gene), collapse = ', '),
               types = paste0(unique(e$type), collapse = ', '),
               experiments = n_distinct(e$experiment))
      })) %>% arrange(desc(experiments)) 

openxlsx::write.xlsx(r2, 'summaries/highQualGenomes/variantSummary.xlsx')


