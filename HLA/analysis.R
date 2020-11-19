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











dbs <- list('A' = 'db/HLA-A_cds_A', 'B' = 'db/HLA-A_cds_B', 'C' = 'db/HLA-B_cds_C', 'D' = 'db/HLA-B_cds_D')

invisible(lapply(split(d, paste0(d$sample, d$exp)), function(x){

  if(file.exists(paste0('touch done/', paste0(x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp)))) return()
  x <- arrange(x, dir)

  system(paste0(bwaPath, ' mem -t 30 -M ', dbs[[x$exp[1]]],
                ' seqData/', x[1,]$file, ' seqData/', x[2,]$file, ' > o.sam'))
  system(paste0(samtoolsBin, '/samtools view -S -b o.sam  > o.bam'))
  invisible(file.remove('o.sam'))
  system(paste0(samtoolsBin, '/samtools view -q 30 -f 0x2 -b o.bam > o.filt.bam'))
  invisible(file.remove('o.bam'))
  system(paste0(samtoolsBin, '/samtools sort -o o.filt.sorted.bam o.filt.bam'))
  invisible(file.remove('o.filt.bam'))
  system(paste0(samtoolsBin, '/samtools index o.filt.sorted.bam'))
  
  
  write(c('new',
          paste0('genome ', dbs[[x$exp[1]]], '.genome'),
          'load o.filt.sorted.bam',
          'snapshotDirectory read_pileup_images/',
          paste0('goto ', sub('db/', '', dbs[[x$exp[1]]]), ':1-', width(readFasta(paste0(dbs[[x$exp[1]]], '.ff')))),
          'sort position',
          'collapse',
          paste0('snapshot ', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.png'),
          'exit'), file = 'batchFile', append = FALSE)
  
  system('xvfb-run --auto-servernum  /home/everett/ext/IGV/igv.sh -b batchFile')
  
  system(paste0(bcftoolsBin, '/bcftools mpileup -A --max-depth 100000 -Ou -f ', paste0(dbs[[x$exp[1]]], '.ff'), 
                ' o.filt.sorted.bam |  ', bcftoolsBin,  '/bcftools call -mv -Oz ', 
                '-o o.filt.sorted.vcf.gz'))
  
  system(paste0(bcftoolsBin, "/bcftools filter -i'QUAL>20&DP>50' o.filt.sorted.vcf.gz -O z -o o.filt.sorted.vcf2.gz"))
  
  z <- system('zcat o.filt.sorted.vcf2.gz', intern = TRUE)
  header <- z[grepl('^#', z)]
  body  <- z[! grepl('^#', z)]
  body1 <- z[! grepl('1/1', z)]
  body2 <- z[grepl('1/1', z)]
  
  write(c(header, body1), file = 'o.filt.sorted.vcf2.single')
  write(c(header, body2), file = 'o.filt.sorted.vcf2.double')
  system('bgzip o.filt.sorted.vcf2.single')
  system('bgzip o.filt.sorted.vcf2.double')
  
  tryCatch({
     system(paste0(bcftoolsBin, '/bcftools index o.filt.sorted.vcf2.gz'))
     system(paste0('zcat o.filt.sorted.vcf2.gz > variant_calls/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.vcf'))
     system(paste0('cat ', dbs[[x$exp[1]]], '.ff | ', bcftoolsBin, '/bcftools consensus ',  
                   'o.filt.sorted.vcf2.gz > consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.fasta'))
     
     o <- Biostrings::readDNAStringSet(paste0('consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.fasta'))
     names(o) <- paste0(x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp)
     Biostrings::writeXStringSet(o, file = paste0('consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.fasta'), append = FALSE)
  })
  
  tryCatch({
    system(paste0(bcftoolsBin, '/bcftools index o.filt.sorted.vcf2.single.gz'))
    system(paste0('zcat o.filt.sorted.vcf2.single.gz > variant_calls/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.vcf.single'))
    system(paste0('cat ', dbs[[x$exp[1]]], '.ff | ', bcftoolsBin, '/bcftools consensus ',  
                  'o.filt.sorted.vcf2.single.gz > consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.single.fasta'))
    
    o <- Biostrings::readDNAStringSet(paste0('consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.single.fasta'))
    names(o) <- paste0(x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp)
    Biostrings::writeXStringSet(o, file = paste0('consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.single.fasta'), append = FALSE)
  })
  
  tryCatch({
    system(paste0(bcftoolsBin, '/bcftools index o.filt.sorted.vcf2.double.gz'))
    system(paste0('zcat o.filt.sorted.vcf2.double.gz > variant_calls/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.vcf.double'))
    system(paste0('cat ', dbs[[x$exp[1]]], '.ff | ', bcftoolsBin, '/bcftools consensus ',  
                  'o.filt.sorted.vcf2.double.gz > consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.double.fasta'))
    
    o <- Biostrings::readDNAStringSet(paste0('consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.double.fasta'))
    names(o) <- paste0(x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp)
    Biostrings::writeXStringSet(o, file = paste0('consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.double.fasta'), append = FALSE)
  })
     
  invisible(file.remove(list.files(pattern = '^o\\.|^batch')))
  system(paste0('touch done/', paste0(x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp)))
}))



f <- tibble(file =list.files('read_pileup_images', pattern = 'HLA-A_A.png|HLA-B_C.png'))
f$sample <- unlist(lapply(strsplit(f$file, '_'), '[', 1))
f$exp <- unlist(lapply(strsplit(f$file, '_'), '[', 2))

invisible(lapply(split(f, f$sample), function(x){
  x <- arrange(x, exp)
  f <- image_read(paste0('read_pileup_images/', x[1,]$file))
  a <- image_crop(f, "4500x45+170+130")
  a <- image_scale(a, '100%x150%')
  a <- image_border(a, "#FFFFFF", "0x15")
  s <- stringr::str_extract(x[1,]$sample, '[A-Z]+')
  x[1,]$sample <- sub(s, paste0(' ', s, ' '), x[1,]$sample)
  a <- image_annotate(a, text = paste0(x[1,]$sample, '   ', x[1,]$exp), location = "+15+85", font = "Arial", size = 12)
  image_write(a, path = "a.png", format = "png")
  
  f <- image_read(paste0('read_pileup_images/', x[2,]$file))
  a <- image_crop(f, "4500x45+170+130")
  a <- image_scale(a, '100%x150%')
  a <- image_border(a, "#FFFFFF", "0x15")
  s <- stringr::str_extract(x[1,]$sample, '[A-Z]+')
  a <- image_annotate(a, text = paste0(x[1,]$sample, '   ', x[2,]$exp), location = "+15+85", font = "Arial", size = 12)
  image_write(a, path = "b.png", format = "png")
  
  system(paste0('convert a.png b.png +append annotated_read_pileup_images/', x[2,]$sample, '.png'))
  invisible(file.remove(c('a.png', 'b.png')))
}))




# 
# invisible(lapply(split(d, d$sample), function(x){
#   a <- readDNAStringSet(paste0('consensus_seqs/', x$sample[1], '_HLA-A_A.fasta'))
#   b <- readDNAStringSet(paste0('consensus_seqs/', x$sample[1], '_HLA-B_C.fasta'))
#   o <- DNAStringSet(paste0(as.character(a), as.character(b)))
#   names(o) <- x$sample[1]
#   writeXStringSet(o, file = paste0('tree_sequences/', x$sample[1], '.fasta'), append = FALSE)
# }))
# 
# system('cat tree_sequences/* > consensusSeqs.fasta')
# system(paste0(mafftPath, ' --thread 20 --globalpair --maxiterate 5 consensusSeqs.fasta > consensusSeqs.mafft'))
# 
# v <- ape::read.dna("consensusSeqs.mafft", format="fasta")
# v_phyDat <- phangorn::phyDat(v, type = "DNA", levels = NULL)
# 
# dna_dist <- phangorn::dist.ml(v_phyDat, model="JC69")
# 
# # Create a data frame to plot tree.
# dendr <- ggdendro::dendro_data(hclust(dna_dist, method='average'), type="rectangle") 
# 
# concensusSeqPhyloPlot <- 
#   ggplot() + 
#   geom_segment(data=ggdendro::segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
#   geom_text(data=ggdendro::label(dendr), aes(x=x, y=y, label=label, hjust=0), size=3) +
#   coord_flip() + scale_y_reverse(expand=c(0.2, 0)) + 
#   theme(axis.line.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.title.y=element_blank(),
#         panel.background=element_rect(fill="white"),
#         panel.grid=element_blank())
# 
# ggsave(concensusSeqPhyloPlot, file = 'concensusSeqPhyloPlot.pdf', width=8, units = 'in')
# 
# 
# 
