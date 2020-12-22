library(ShortRead)
library(dplyr)
library(lubridate)
options(stringsAsFactors = FALSE)
states <- readLines('US_state_names.txt')

nextStrainMetaData <- read.table('nextStrainMetaData.tsv', sep = '\t', header = TRUE, quote = '', comment.char = '', fill = TRUE)

gisaid <- readFasta('GISAID/sequences_2020-08-07_17-50.fasta')
gisaid <- gisaid[! duplicated(gisaid@sread)]

g <- readFasta('../../summaries/concensusSeqs95.fasta')

# Remove 228_NP_OP_20200424
g <- g[! grepl('228\\|NP\\-OP\\|20200424', as.character(g@id))]

g <- g[! as.character(g@id) %in% c('CCLB|Vero cells|20200328', 'E6|Vero cells|20200328')]
g <- g[! duplicated(g@sread)]

# Remove overlaps with GISAID and our genomes.
o <- Reduce('append', list(g, gisaid))
gisaid <- gisaid[! as.character(gisaid@id) %in% as.character(o[duplicated(o@sread)]@id)]


metaData <- read.table('GISAID/metadata_2020-08-07_19-18.tsv', header = TRUE, sep = '\t', quote = '', fill = TRUE)

# Remove rows from the meta data for those genomes which were removed due to sequence duplication.
metaData <- subset(metaData, strain %in% as.character(gisaid@id))


metaData$virus <- 'SARS-CoV-2'
metaData <- metaData[grepl('\\d\\d\\d\\d\\-\\d\\d\\-\\d\\d', metaData$date),]
metaData$date2 <- ymd(metaData$date)

metaData$period <- ifelse(metaData$date2 <= ymd('2020-01-31'), 'A', 
                          ifelse(metaData$date2 >= ymd('2020-02-01') & metaData$date2 <= ymd('2020-03-10'), 'B', 'C'))

us <- subset(metaData, metaData$division_exposure %in% states & country_exposure == 'USA')  # "2020-01-25" "2020-01-25" "2020-01-29" "2020-01-29" "2020-02-10"
ny <- subset(us, division_exposure == 'New York')     # "2020-03-02" "2020-03-02" "2020-03-03" "2020-03-03" "2020-03-03"
nj <- subset(us, division_exposure == 'New Jersey')   # "2020-03-03" "2020-03-04" "2020-03-04" "2020-03-10" "2020-03-13"
md <- subset(us, division_exposure == 'Maryland')     #  "2020-03-04" "2020-03-04" "2020-03-04" "2020-03-07" "2020-03-08"
pa <- subset(us, division_exposure == 'Pennsylvania') # "2020-03-05" "2020-03-05" "2020-03-06" "2020-03-07" "2020-03-07"

non_us <- subset(metaData, !strain %in% us$strain)


set.seed(1)
# Limit non_us strains to those in nextStrain since they already filtered their strains to be good rep. strains.
a <- subset(non_us, period == 'A' & strain %in% nextStrainMetaData$Strain)
n <- 34

period_A <- bind_rows(lapply(split(a, a$pangolin_lineage), function(x){
  if(nrow(x) > n){
    return(sample_n(x, n))
  } else {
    return(x)
  }
}))

period_A$source <- 'Non-USA'

#--------------------------------------------------------------------------------------------------

o <- subset(non_us, period == 'B' & strain %in% nextStrainMetaData$Strain)
n <- 6
a <- bind_rows(lapply(split(o, o$pangolin_lineage), function(x){
  if(nrow(x) > n){
    return(sample_n(x, n))
  } else {
    return(x)
  }
}))

o <- subset(us, ! division_exposure %in% c('New York', 'New Jersey', 'Maryland', 'Pennsylvania') & period == 'B')
n <- 13
b <- bind_rows(lapply(split(o, o$pangolin_lineage), function(x){
  if(nrow(x) > n){
    return(sample_n(x, n))
  } else {
    return(x)
  }
}))

c <- subset(us, period == 'B' & division_exposure == 'New York')
d <- subset(us, period == 'B' & division_exposure == 'New Jersey')
e <- subset(us, period == 'B' & division_exposure == 'Maryland')
a$source  <- 'Non-USA'
b$source  <- 'USA-other'
c$source  <- 'NY'
d$source  <- 'NJ'
e$source  <- 'MD'
period_B <- bind_rows(period_A, a, b, c, d, e)  %>% mutate(period = 'B')

#----------------------------------------------------------------------------------


o <- subset(non_us, period == 'C' & strain %in% nextStrainMetaData$Strain)
n <- 2
a <- bind_rows(lapply(split(o, o$pangolin_lineage), function(x){
  if(nrow(x) > n){
    return(sample_n(x, n))
  } else {
    return(x)
  }
}))


o <- subset(us, ! division_exposure %in% c('New York', 'New Jersey', 'Maryland', 'Pennsylvania') & period == 'C')
n <- 2
b <- bind_rows(lapply(split(o, o$pangolin_lineage), function(x){
  if(nrow(x) > n){
    return(sample_n(x, n))
  } else {
    return(x)
  }
}))


o <- subset(us, period == 'C' & division_exposure == 'New York')
n <- 7
c <- bind_rows(lapply(split(o, o$pangolin_lineage), function(x){
  if(nrow(x) > n){
    return(sample_n(x, n))
  } else {
    return(x)
  }
}))

d <- subset(us, period == 'C' & division_exposure == 'New Jersey')

o <- subset(us, period == 'C' & division_exposure == 'Maryland')
n <- 30
e <- bind_rows(lapply(split(o, o$pangolin_lineage), function(x){
  if(nrow(x) > n){
    return(sample_n(x, n))
  } else {
    return(x)
  }
}))

a$source  <- 'Non-USA'
b$source  <- 'USA-other'
c$source  <- 'NY'
d$source  <- 'NJ'
e$source  <- 'MD'
period_C <- bind_rows(period_B, a, b, c, d, e) %>% mutate(period = 'C')

periods <- select(bind_rows(period_A, period_B, period_C), strain, virus, date, pangolin_lineage, GISAID_clade, country_exposure, division_exposure, period, source)
periods <- bind_rows(lapply(split(periods, periods$period), function(x) x[! duplicated(x$strain),]))

system('rm -r nextStrainRun_A/results/* nextStrainRun_A/data/*')
system('rm -r nextStrainRun_B/results/* nextStrainRun_B/data/*')
system('rm -r nextStrainRun_C/results/* nextStrainRun_C/data/*')

writeFasta(gisaid[as.character(gisaid@id) %in% subset(periods, period == 'A')$strain], file = 'nextStrainRun_A/data/sequences.fasta')
write.table(subset(periods, period == 'A'), file = 'nextStrainRun_A/data/metadata.tsv', sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

writeFasta(gisaid[as.character(gisaid@id) %in% subset(periods, period == 'B')$strain], file = 'nextStrainRun_B/data/sequences.fasta')
write.table(subset(periods, period == 'B'), file = 'nextStrainRun_B/data/metadata.tsv', sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)

a <- tibble(strain = as.character(g@id), 
            virus = 'SARS-CoV-2', 
            date = unlist(lapply(strsplit(as.character(g@id), '\\|'), function(x) as.character(ymd(x[3])))),
            pangolin_lineage = 'X',
            GISAID_clade = 'X', 
            country_exposure = 'USA', 
            division_exposure = 'Pennsylvania', 
            period = 'C', 
            source = 'Philadelphia')

writeFasta(Reduce('append', list(g, gisaid[as.character(gisaid@id) %in% subset(periods, period == 'C')$strain])), file = 'nextStrainRun_C/data/sequences.fasta')
write.table(bind_rows(a, subset(periods, period == 'C')), file = 'nextStrainRun_C/data/metadata.tsv', sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)


# Run next Strain sofware...
# conda activate nextstrain
# cd into nextStrainRun dir
# sh ../nextStrainComms.sh


library(ggplot2)

createTree <- function(dir){
  o <- phylogram::read.dendrogram(file.path(dir, 'results/tree.nwk'))
 
  # Remove outliers by looking at the first points in the pointData df below. 
  o <- phylogram::prune(o, 'Malaysia/IMR_WC085/2020')
  o <- phylogram::prune(o, 'USA/WA-UW-4572/2020')
  
  m <- read.table(file.path(dir, 'data/metadata.tsv'), sep = '\t', header = TRUE)
  dendr <- ggdendro::dendro_data(o)
  # browser()
  pointData <- group_by(subset(ggdendro::segment(dendr), x %in% ggdendro::label(dendr)$x), x) %>%
               top_n(-1, wt = yend) %>% 
               slice(1) %>%
               ungroup() %>%
               left_join(select(ggdendro::label(dendr), -y), by = 'x') %>%
               left_join(select(m, strain, source), by = c('label' = 'strain'))

  pointData$source <- factor(as.character(pointData$source), levels = c("Non-USA", "USA-other", "NJ", "NY", "MD", "Philadelphia"))

  a <- subset(pointData, source == 'Philadelphia')
  b <- subset(pointData, source != 'Philadelphia')
  pointData <- bind_rows(b, a)

  ggplot() +
    geom_segment(data=ggdendro::segment(dendr), aes(x=x, y=y, xend=xend, yend=yend), size = 0.2, color = 'gray') + 
    scale_color_manual(values = c('gray50', 'black', 'gold2', 'red2', 'green3', 'dodgerblue'), drop = FALSE) +
    geom_point(data=pointData, aes(x=x, y=yend, color = source), size = 1.5) +
    coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
    theme(axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
}


A <- createTree('nextStrainRun_A')
B <- createTree('nextStrainRun_B')
C <- createTree('nextStrainRun_C')

ggsave(A, filename = 'A.pdf', height = 3,  width = 3, units = 'in', useDingbats = FALSE)
ggsave(B, filename = 'B.pdf', height = 5,  width = 4, units = 'in', useDingbats = FALSE)
ggsave(C, filename = 'C.pdf', height = 10, width = 5, units = 'in', useDingbats = FALSE)
