library(tidyverse)
library(Biostrings)
mafftPath <- '~/ext/mafft/bin/mafft'


samples <- read.table('data/samples.tsv', sep= '\t', header = TRUE, quote = '', stringsAsFactors = FALSE)

VSPs <- c('VSP0526', 'VSP0527', 'VSP0528', 'VSP0529', 'VSP0530', 'VSP0531', 'VSP0532', 'VSP0533',
          'VSP0534', 'VSP0535', 'VSP0536', 'VSP0537', 'VSP0538', 'VSP0539', 'VSP0540', 'VSP0541',
          'VSP0542', 'VSP0543', 'VSP0544', 'VSP0545', 'VSP0546', 'VSP0547', 'VSP9977')

s <- Reduce('append', lapply(VSPs, function(x){
       message(x)
       load(file.path('summaries/VSPdata/', paste0(x, '-1.experiment.RData')))
       if('refGenomePercentCovered_5reads' %in% names(opt)){
         if(opt$refGenomePercentCovered_5reads >= 0.90){
           o <- DNAStringSet(opt$concensusSeq)
           names(o) <- samples[samples$VSP == x,]$patient_id
           return(o)
         } else {
           return(DNAStringSet())           
         }
       } else {
         return(DNAStringSet())  
       }
     }))

writeXStringSet(s, file = 'tmp.fasta')
system(paste0(mafftPath, ' --phylipout --thread 10 --globalpair --maxiterate 5 tmp.fasta > tmp.mafft'))

comm <- 'raxmlHPC-PTHREADS-SSE3 -s tmp.mafft -m GTRGAMMA -T 25 -n raxmlOut -f a -x 12345 -p 12345 -N autoMRE -o "England_EPI_ISL_747519"'

dendr <- ggdendro::dendro_data(phylogram::read.dendrogram('RAxML_bestTree.raxmlOut'), type="rectangle")
segments <- ggdendro::segment(dendr)
labels <- ggdendro::label(dendr)

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


