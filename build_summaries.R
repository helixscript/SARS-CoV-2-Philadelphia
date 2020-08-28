library(tidyverse)
library(parallel)
library(Biostrings)
options(stringsAsFactors = FALSE)
overWriteSubjectReports <- FALSE
mafftPath <- '~/ext/mafft/bin/mafft'
CPUs <- 25

# Read in sample data table.
samples <- read.table('data/samples.tsv', sep= '\t', header = TRUE, quote = '', stringsAsFactors = FALSE)
sampleInputs <- read.table('data/sampleInputs.tsv', sep= '\t', header = TRUE, stringsAsFactors = FALSE)


#---

#v <- unique(str_extract(list.files('summaries/VSPdata'), 'VSP\\d+'))
#samples <- filter(samples, VSP %in% v)
#o <- group_by(samples, patient_id) %>% summarise(timePoints = n_distinct(sample_date)) %>% ungroup() %>% filter(timePoints > 1) %>% arrange(desc(timePoints))

samples <- samples[! grepl('MisC', samples$patient_id),]


#---

# Create subject summary reports for all subjects with experimental data (VSP data objects in VSP_data/).
cluster <- makeCluster(CPUs)
clusterExport(cluster, c('samples', 'sampleInputs', 'overWriteSubjectReports'))
if(! dir.exists('summaries/patientReportsData')) dir.create('summaries/patientReportsData')


# Start subject report log.
write(date(), file = 'logs/buid_subject_reports.log', append = FALSE)


# Create a table of patient ids with a chunking vector.
d <- data.frame(patient_id = unique(samples$patient_id))
d$s <- ntile(1:nrow(d), CPUs)


# Build patient reports.
invisible(parLapply(cluster, split(d, d$s), function(p){
#invisible(lapply(split(d, d$s), function(p){  
  library(tidyverse)
  
  invisible(sapply(p$patient_id, function(patientID){
    x <- subset(samples, patient_id == patientID)
  
    if(overWriteSubjectReports == FALSE & file.exists(paste0('summaries/patientReports/', x$patient_id[1], '.pdf'))) return()
    
    message('Starting subject: ', x$patient_id[1], ', samples: ', paste0(x$VSP, collapse = '|'))
    files <- list.files('summaries/VSPdata', pattern = paste0(x$VSP, collapse = '|'), full.names = TRUE)

    if(length(files) == 0){
      write(paste0('[.] No VSP data files for subject ', x$patient_id[1]), file = 'logs/buid_subject_reports.log', append = TRUE)
      return()
    }
  
    dat <- lapply(files, function(f){
           # Here we load VSP data files where there are two types of files, analyses of individual sequencing 
           # experiments and composite analyses where multiple experiments where combined.
             load(f)
             opt$vsp <- str_extract(f, 'VSP\\d+')
             opt$seq_sample <-  str_extract(f, 'VSP\\d+\\-?\\d+?[ab]?')
           
             d <- subset(samples, VSP == opt$vsp)
             opt$sample_id <- d$sample_id[1]
             opt$subject <- d$patient_id[1]
             opt$date <- d$sample_date[1]
             opt$sampleType = d$sample_type[1]
             opt$genomesPerMicroLiter = as.numeric(d$N1_copies_per_ul[1])
             if(! grepl('composite', f)){
               opt$inputGenomes <- opt$genomesPerMicroLiter * subset(sampleInputs, Sample_Name == opt$seq_sample)$uL_Inputs
             } else {
               opt$inputGenomes <- NA
             }
             opt$concensusSeq = as.character(opt$concensusSeq)
             return(opt)
           })
         
    names(dat) <- unlist(lapply(dat, '[[', 'seq_sample'))
         
    save(dat, file = paste0('summaries/patientReportsData/', dat[[1]]$subject, '.RData'))
         
     result = tryCatch({
              rmarkdown::render('report.Rmd',
                                output_file = paste0('summaries/patientReports/', dat[[1]]$subject, '.pdf'),
                                params = list('date'  = format(Sys.time(), "%Y-%m-%d"),
                                              'title' = paste0('COVID-19 subject ', dat[[1]]$subject)))
              }, error = function(e) {
                 write(paste0('[!] Failed to create subject report for ', dat[[1]]$subject), file = 'logs/buid_subject_reports.log', append = TRUE)
              })
    }))
}))


# Collate sampe summary tables.
summary <- bind_rows(lapply(list.files('summaries/sampleSummaries', full.names = TRUE), function(x){
             t <- read.table(x, sep = '\t', header = TRUE)
             t$Subject <- as.character(t$Subject)
             t$largestContig <- as.numeric(t$largestContig)
             t$percentRefReadCoverage <- as.numeric(sub('%', '', t$percentRefReadCoverage))
             t$percentRefReadCoverage5 <- as.numeric(sub('%', '', t$percentRefReadCoverage5))
             t$sampleDate2 <- gsub('-', '', as.character(lubridate::mdy(t$sampleDate)))
             t
}))


summary <- summary[!grepl('MisC', summary$Subject),]

openxlsx::write.xlsx(summary, file = 'summaries/sampleSummary.xlsx')



# Select a sample to represent each subject, sample type, time point combination.
# Select composite samples when available otherwise select the samples with the greated coverage.
# Break ties with read coverage percentages.

representativeSampleSummary <- function(summary, minPercentRefReadCoverage5){
  bind_rows(lapply(split(summary, paste(summary$Subject, summary$sampleType, summary$sampleDate2)), function(x){
    r <- subset(x, type == 'composite')
    if(nrow(r) > 0){
      r <- top_n(r, 1, wt = percentRefReadCoverage5) %>% dplyr::slice(1)
    } else{
      r <- top_n(x, 1, wt = percentRefReadCoverage5) %>% dplyr::slice(1)
    }
    r
  })) %>% dplyr::filter(percentRefReadCoverage5 >= minPercentRefReadCoverage5)
}

representativeSampleSummary_90 <- representativeSampleSummary(summary, 90)
openxlsx::write.xlsx(representativeSampleSummary_90, file = 'summaries/representativeSampleSummary_90.xlsx')


# The reference genome is 29,882 NT long
# 90% of this length is 26,893 NT.
# Select summary records with longer contigs, break ties with read coverage %.

minContigLengthKD <- 26893 / 1000
longContigSampleSummary <- 
  group_by(dplyr::filter(summary, largestContig >= minContigLengthKD), Subject, sampleType, sampleDate2) %>%
  dplyr::top_n(1, wt = largestContig) %>%
  dplyr::arrange(desc(percentRefReadCoverage)) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>% 
  dplyr::arrange(Subject)

openxlsx::write.xlsx(longContigSampleSummary, file = 'summaries/longContigSampleSummary_gte90percentRef.xlsx')



# Retrieve and write out long contings.
contigs <- Reduce('append', lapply(1:nrow(longContigSampleSummary), function(x){
  x <- longContigSampleSummary[x,]
  x$type <- ifelse(x$type == 'composite', x$type, 'experiment')
  load(paste0('summaries/VSPdata/', x$exp, '.', x$type, '.RData'))
  contig <- opt$contigs[order(width(opt$contigs), decreasing = TRUE)][1]
  names(contig) <- paste0(x$Subject, '|', x$sampleType, '|', x$sampleDate2)
  contig
}))

writeXStringSet(contigs, file = 'summaries/contigs90.fasta')



# Retrieve consensus sequences from representative samples.
retrieveConcensusSeqs <- function(summary){
  Reduce('append', lapply(1:nrow(summary), function(x){
    x <- summary[x,]
    f <- paste0('summaries/VSPdata/', x$exp, '.', ifelse(x$type == 'composite', 'composite', 'experiment'), '.RData')
    if(! file.exists(f)) stop('Can not locate VSP data file -- ', f)
    load(f)
    d <- DNAStringSet(opt$concensusSeq)
    names(d) <- paste0(x$Subject, '|', x$sampleType, '|', x$sampleDate2)
    d
  }))
}


concensusSeqs90 <- retrieveConcensusSeqs(representativeSampleSummary_90)
writeXStringSet(concensusSeqs90, file = 'summaries/concensusSeqs90.fasta')




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

# The colated variant table has all variants for each position. 
# Here we reduce it to the dominant variants. 

variantTable <- filter(variantTable,  percentAlt > 0.5 & ! ALT == 'e') %>%
                group_by(sample, POS) %>%
                top_n(1, wt = percentAlt) %>%
                ungroup()


# Create a table of variants with subject and sample counts.
variantTable$subject <- unlist(lapply(str_extract_all(variantTable$sample, '([^\\|]+)'), '[', 1))
variantSummary <- 
  group_by(variantTable, POS) %>%
  summarise(nSubjects = n_distinct(subject), nSamples = n_distinct(sample), gene = genes[1], mutation = type[1]) %>%
  ungroup() %>% 
  arrange(desc(nSamples))

openxlsx::write.xlsx(variantSummary, file = 'summaries/variantSummary.xlsx')

v <- group_by(variantTable, subject, sample) %>% mutate(variants = paste0(POS, '|', genes)) %>% ungroup() %>% select(subject, sample, variants)
v2 <- reshape2::dcast(v, subject+sample~variants, value.var = 'variants')
openxlsx::write.xlsx(v2, file = 'summaries/variantSummary2.xlsx')


# Align the concensus sequences.
if(! file.exists('summaries/concensusSeqs90.mafft')) system(paste0(mafftPath, ' --thread 25 --globalpair --maxiterate 50 summaries/concensusSeqs90.fasta > summaries/concensusSeqs90.mafft'))

# Build phylogenetic tree.
v <- ape::read.dna("summaries/concensusSeqs90.mafft", format="fasta")
v_phyDat <- phangorn::phyDat(v, type = "DNA", levels = NULL)

# Determine best mode and create a distance matrix.
mt <- phangorn::modelTest(v_phyDat)
dna_dist <- phangorn::dist.ml(v_phyDat, model="JC69")

# Create a data frame to plot tree.
dendr <- ggdendro::dendro_data(hclust(dna_dist, method='average'), type="rectangle") 

concensusSeqPhyloPlot <- 
  ggplot() + 
  geom_segment(data=ggdendro::segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  geom_text(data=ggdendro::label(dendr), aes(x=x, y=y, label=label, hjust=0), size=3) +
  coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        panel.background=element_rect(fill="white"),
        panel.grid=element_blank())

ggsave(concensusSeqPhyloPlot, height = 10, width = 12, units = 'in', file = 'summaries/concensusSeqs90PhyloPlot.pdf')


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

ggsave(genomeVariantPlot, height = 10, width = 12, units = 'in', file = 'summaries/genomeVariantPlot.pdf')


# Create a single results download.
file.remove('summaries.zip')
system('zip summaries/summaries.zip summaries/*.xlsx summaries/*.pdf summaries/*.fasta summaries/patientReports/*')



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
      #x$coverage5 <- NA
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
