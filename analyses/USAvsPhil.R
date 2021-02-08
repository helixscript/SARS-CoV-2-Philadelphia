library(tidyverse)
library(lubridate)

options(stringsAsFactors = FALSE)

d <- read.table('../data/nextStrain_metaData_northAmerica_Jan2021.txt', sep = '\t', header = TRUE, fill = TRUE, quote = '')
d <- d[d$Country == 'USA' & d$Country.of.Exposure == 'USA',]
d <- d[,c('Strain', 'Collection.Data', 'PANGO.Lineage')]

d$Collection.Data <- ymd(d$Collection.Data)
d$days <- as.integer(d$Collection.Data)
d$dateLabel <- paste0(month(d$Collection.Data), '/', year(d$Collection.Data))
d <- d[d$dateLabel != 'NA/NA',]
d <- arrange(d, days)  


lineages <- names(sort(table(d[grepl('2021', d$Collection.Data),]$PANGO.Lineage), decreasing = TRUE)[1:13])

# Read in Philadelphia lineages.
p <- openxlsx::read.xlsx('../summaries/highQualGenomes/lineages.xlsx')

p <- bind_rows(lapply(1:nrow(p), function(x){
  x <- p[x,]
  x$date <- unlist(str_split(x$sample, '\\|'))[4]
  o <- ymd(x$date)
  x$days <- as.integer(o)
  x$dateLabel <- paste0(month(o), '/', year(o))
  rename(x, PANGO.Lineage = lineage)
}))

d <- select(d, PANGO.Lineage, dateLabel, days)
p <- select(p, PANGO.Lineage, dateLabel, days)


# Determine which of the more abundant Philadelphia linages are not in the nextStrain 2021 lineages 
# and add 5 of them to the legend.
o <- names(sort(table(p$PANGO.Lineage), decreasing = TRUE))
lineages <- c(lineages, o[! o %in% lineages][1:2])

# Order the lineages to match the most abundant lineages/
o <- names(sort(table(c(d$PANGO.Lineage, p$PANGO.Lineage)), decreasing = TRUE))
lineages <- c( 'Other', o[o %in% lineages])


d$dateLabel <- factor(d$dateLabel, levels = unique(d$dateLabel))
p$dateLabel <- factor(p$dateLabel, levels = unique(d$dateLabel))

d$PANGO.Lineage <- ifelse(d$PANGO.Lineage %in% lineages, d$PANGO.Lineage, 'Other')
d$PANGO.Lineage <- factor(d$PANGO.Lineage, levels = lineages)

p$PANGO.Lineage <- ifelse(p$PANGO.Lineage %in% lineages, p$PANGO.Lineage, 'Other')
p$PANGO.Lineage <- factor(p$PANGO.Lineage, levels = lineages)

colors <- c('gray90', grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(lineages) - 1)) 

colors2 <- c('gray90',  '#3e647d', '#7b92a8', '#82c0e9', '#2d6d66', '#bfa19c', '#008bbc',
            '#97b6b0', '#d7d29e', '#1a476f', '#90353b', '#9c8847', '#938dd2', '#6e8e84',
            '#c10534', '#bab05d')
            
nextStrainPlot <- 
  ggplot(d, aes(dateLabel, fill = PANGO.Lineage)) + 
  theme_bw() +
  geom_bar(stat = 'count') +
  scale_fill_manual(name = 'Lineage', values = colors, drop = FALSE) +
  guides(fill=guide_legend(ncol=2)) +
  labs(x = 'Sample Date', y = 'Genomes') +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16), legend.title=element_text(size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(nextStrainPlot, file = 'nextStrain_USA.pdf', units = 'in', width = 10, height = 10)
nextStrainPlot <- nextStrainPlot + scale_fill_manual(name = 'Lineage', values = colors2, drop = FALSE)
ggsave(nextStrainPlot, file = 'nextStrain_USA_palette2.pdf', units = 'in', width = 10, height = 10)

PhlStrainPlot <- 
  ggplot(p, aes(dateLabel, fill = PANGO.Lineage)) + 
  theme_bw() +
  geom_bar(stat = 'count') +
  scale_fill_manual(name = 'Lineage', values = colors, drop = FALSE) +
  guides(fill=guide_legend(ncol=2)) +
  labs(x = 'Sample Date', y = 'Genomes') +
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16), legend.title=element_text(size=16),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

ggsave(PhlStrainPlot, file = 'Philadelphia_.pdf', units = 'in', width = 10, height = 10)
PhlStrainPlot <- PhlStrainPlot + scale_fill_manual(name = 'Lineage', values = colors2, drop = FALSE)
ggsave(PhlStrainPlot, file = 'Philadelphia_palette2.pdf', units = 'in', width = 10, height = 10)
                      