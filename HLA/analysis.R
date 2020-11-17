library(dplyr)
library(ShortRead)
library(Biostrings)
library(ggplot2)

bwaPath <- '~/ext/bwa'
samtoolsBin <- '~/ext/samtools/bin'
bcftoolsBin <- '~/ext/bcftools/bin'
mafftPath   <- '~/ext/mafft/bin/mafft'

d <- tibble(file = list.files('seqData'),
            sample = unlist(lapply(strsplit(file, '_'), 
                                   function(x) substr(x[1], 1, nchar(x[1])-1))),
            exp = unlist(lapply(strsplit(file, '_'), 
                                function(x) substr(x[1], nchar(x[1]), nchar(x[1])))),
            gene = ifelse(unlist(lapply(strsplit(file, '_'), 
                                       function(x) substr(x[1], nchar(x[1]), nchar(x[1])))) %in% c('A', 'B'), 'HLA-A', 'HLA-B'),
            dir = ifelse(grepl('_R1_', file), 'R1', 'R2')) %>%
     filter(! grepl('Control', sample))

invisible(lapply(split(d, paste0(d$sample, d$exp)), function(x){

  x <- arrange(x, dir)
  
  system(paste0(bwaPath, ' mem -t 30 -M db/', ifelse(x$gene[1] == 'HLA-A', 'HLA-A_cds', 'HLA-B_cds'),
                ' seqData/', x[1,]$file, ' seqData/', x[2,]$file, ' > o.sam'))
  system(paste0(samtoolsBin, '/samtools view -S -b o.sam  > o.bam'))
  invisible(file.remove('o.sam'))
  system(paste0(samtoolsBin, '/samtools view -q 30 -f 0x2 -b o.bam > o.filt.bam'))
  invisible(file.remove('o.bam'))
  system(paste0(samtoolsBin, '/samtools sort -o o.filt.sorted.bam o.filt.bam'))
  invisible(file.remove('o.filt.bam'))
  system(paste0(samtoolsBin, '/samtools index o.filt.sorted.bam'))
  
  write(c('new',
          paste0('genome db/', ifelse(x$gene[1] == 'HLA-A', 'HLA-A_cds', 'HLA-B_cds'), '.genome'),
          'load o.filt.sorted.bam',
          'snapshotDirectory read_pileup_images/',
          paste0('goto ', ifelse(x$gene[1] == 'HLA-A', 'HLA-A_cds:1-1098', 'HLA-B_cds:1-1089')),   
          'sort position',
          'collapse',
          paste0('snapshot ', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.png'),
          'exit'), file = 'batchFile', append = FALSE)
  
  system('xvfb-run --auto-servernum  /home/everett/ext/IGV/igv.sh -b batchFile')
  
  system(paste0(samtoolsBin, '/samtools mpileup --output ', paste0('read_pileup_data/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.pileup'), 
                ' --max-depth 100000 -f db/', ifelse(x$gene[1] == 'HLA-A', 'HLA-A_cds', 'HLA-B_cds'), ' o.filt.sorted.bam'))
  
  system(paste0(bcftoolsBin, '/bcftools mpileup -A --max-depth 100000 -Ou -f db/', ifelse(x$gene[1] == 'HLA-A', 'HLA-A_cds', 'HLA-B_cds'), 
                ' o.filt.sorted.bam |  ', bcftoolsBin,  '/bcftools call -mv -Oz ', 
                '-o o.filt.sorted.vcf.gz'))
  
  system(paste0(bcftoolsBin, "/bcftools filter -i'QUAL>20&DP>50' o.filt.sorted.vcf.gz -O z -o o.filt.sorted.vcf2.gz"))
  
  system(paste0(bcftoolsBin, '/bcftools index o.filt.sorted.vcf2.gz'))
  
  system(paste0('zcat o.filt.sorted.vcf2.gz > variant_calls/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.vcf'))
  
  system(paste0('cat  db/', ifelse(x$gene[1] == 'HLA-A', 'HLA-A_cds', 'HLA-B_cds'), ' | ', bcftoolsBin, '/bcftools consensus ',  
                'o.filt.sorted.vcf2.gz > consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.fasta'))
  
  o <- Biostrings::readDNAStringSet(paste0('consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.fasta'))
  names(o) <- paste0(x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp)
  Biostrings::writeXStringSet(o, file = paste0('consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.fasta'), append = FALSE)
  
  invisible(file.remove(list.files(pattern = '^o\\.|^batch')))
}))


invisible(lapply(split(d, d$sample), function(x){
  a <- readDNAStringSet(paste0('consensus_seqs/', x$sample[1], '_HLA-A_A.fasta'))
  b <- readDNAStringSet(paste0('consensus_seqs/', x$sample[1], '_HLA-B_C.fasta'))
  o <- DNAStringSet(paste0(as.character(a), as.character(b)))
  names(o) <- x$sample[1]
  writeXStringSet(o, file = paste0('tree_sequences/', x$sample[1], '.fasta'), append = FALSE)
}))

system('cat tree_sequences/* > consensusSeqs.fasta')
system(paste0(mafftPath, ' --thread 20 --globalpair --maxiterate 5 consensusSeqs.fasta > consensusSeqs.mafft'))

v <- ape::read.dna("consensusSeqs.mafft", format="fasta")
v_phyDat <- phangorn::phyDat(v, type = "DNA", levels = NULL)

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

ggsave(concensusSeqPhyloPlot, file = 'concensusSeqPhyloPlot.pdf', width=8, units = 'in')



