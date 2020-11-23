library(dplyr)
library(ShortRead)
library(Biostrings)
library(ggplot2)
library(magick)
tmpFile <- function(){ paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') }

bwaPath <- '~/ext/bwa'
samtoolsBin <- '~/ext/samtools/bin'
bcftoolsBin <- '~/ext/bcftools/bin'
mafftPath   <- '~/ext/mafft/bin/mafft'

d <- tibble(file = list.files('trimmedSeqData'),
            sample = unlist(lapply(strsplit(file, '_'), 
                                   function(x) substr(x[1], 1, nchar(x[1])-1))),
            exp = unlist(lapply(strsplit(file, '_'), 
                                function(x) substr(x[1], nchar(x[1]), nchar(x[1])))),
            gene = ifelse(unlist(lapply(strsplit(file, '_'), 
                                       function(x) substr(x[1], nchar(x[1]), nchar(x[1])))) %in% c('A', 'B'), 'HLA-A', 'HLA-B'),
            dir = ifelse(grepl('_R1_', file), 'R1', 'R2')) %>%
     filter(! grepl('Control', sample))


dbs <- list('A' = 'db/HLA-A.A', 'B' = 'db/HLA-A.B', 'C' = 'db/HLA-B.C', 'D' = 'db/HLA-B.D')

invisible(lapply(split(d, paste0(d$sample, d$exp)), function(x){

  if(file.exists(paste0('touch done/', paste0(x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp)))) return()
  x <- arrange(x, dir)

  system(paste0(bwaPath, ' mem -t 30 -M ', dbs[[x$exp[1]]],
                ' trimmedSeqData/', x[1,]$file, ' trimmedSeqData/', x[2,]$file, ' > o.sam'))
  
  system(paste0('~/ext/samclip --ref ', dbs[[x$exp[1]]], ' o.sam > o.sam2'))
  
  system(paste0(samtoolsBin, '/samtools view -S -b o.sam2  > o.bam'))
  invisible(file.remove(c('o.sam', 'o.sam2')))
  
  system(paste0(samtoolsBin, '/samtools view -q 30 -f 0x2 -b o.bam > o.filt.bam'))
  invisible(file.remove('o.bam'))
  system(paste0(samtoolsBin, '/samtools sort -o o.filt.sorted.bam o.filt.bam'))
  invisible(file.remove('o.filt.bam'))
  system(paste0(samtoolsBin, '/samtools index o.filt.sorted.bam'))
  
  write(c('new',
          paste0('genome ', dbs[[x$exp[1]]], '.genome'),
          'load o.filt.sorted.bam',
          'snapshotDirectory read_pileup_images/',
          paste0('goto ', sub('db/', '', dbs[[x$exp[1]]]), ':1-', width(readFasta(dbs[[x$exp[1]]]))+10),
          'sort position',
          'collapse',
          paste0('snapshot ', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.png'),
          'exit'), file = 'batchFile', append = FALSE)
  
  system('xvfb-run --auto-servernum  /home/everett/ext/IGV/igv.sh -b batchFile')
  
  system(paste0(bcftoolsBin, '/bcftools mpileup -A --max-depth 100000 -Ou -f ', dbs[[x$exp[1]]], 
                ' o.filt.sorted.bam |  ', bcftoolsBin,  '/bcftools call -mv -Oz ', 
                '-o o.filt.sorted.vcf.gz'))
  
  system(paste0(bcftoolsBin, "/bcftools filter -i'QUAL>50&DP>100' o.filt.sorted.vcf.gz -O z -o o.filt.sorted.vcf2.gz"))
  
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
     system(paste0('cat ', dbs[[x$exp[1]]], ' | ', bcftoolsBin, '/bcftools consensus ',  
                   'o.filt.sorted.vcf2.gz > consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.fasta'))
     
     o <- Biostrings::readDNAStringSet(paste0('consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.fasta'))
     names(o) <- paste0(x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp)
     Biostrings::writeXStringSet(o, file = paste0('consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.fasta'), append = FALSE)
  })
  
  tryCatch({
    system(paste0(bcftoolsBin, '/bcftools index o.filt.sorted.vcf2.single.gz'))
    system(paste0('zcat o.filt.sorted.vcf2.single.gz > variant_calls/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.vcf.single'))
    system(paste0('cat ', dbs[[x$exp[1]]], ' | ', bcftoolsBin, '/bcftools consensus ',  
                  'o.filt.sorted.vcf2.single.gz > consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.single.fasta'))
    
    o <- Biostrings::readDNAStringSet(paste0('consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.single.fasta'))
    names(o) <- paste0(x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp)
    Biostrings::writeXStringSet(o, file = paste0('consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.single.fasta'), append = FALSE)
  })
  
  tryCatch({
    system(paste0(bcftoolsBin, '/bcftools index o.filt.sorted.vcf2.double.gz'))
    system(paste0('zcat o.filt.sorted.vcf2.double.gz > variant_calls/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.vcf.double'))
    system(paste0('cat ', dbs[[x$exp[1]]], ' | ', bcftoolsBin, '/bcftools consensus ',  
                  'o.filt.sorted.vcf2.double.gz > consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.double.fasta'))
    
    o <- Biostrings::readDNAStringSet(paste0('consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.double.fasta'))
    names(o) <- paste0(x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp)
    Biostrings::writeXStringSet(o, file = paste0('consensus_seqs/', x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp, '.double.fasta'), append = FALSE)
  })
     
  invisible(file.remove(list.files(pattern = '^o\\.|^batch')))
  system(paste0('touch done/', paste0(x[1,]$sample, '_', x[1,]$gene, '_', x[1,]$exp)))
}))


# Create HLA-A and HLA-B pileup images by extracting regions from the IGV generated snapshots.
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



buildTree <- function(x){
  f <- tmpFile()
  writeXStringSet(x, file = f)
  system(paste0(mafftPath, ' --thread 20 --globalpair --maxiterate 5 ', f, ' > ', f, '.mafft'))

  v <- ape::read.dna(paste0(f, '.mafft'), format="fasta")
  invisible(file.remove(c(f, paste0(f, '.mafft'))))
  v_phyDat <- phangorn::phyDat(v, type = "DNA", levels = NULL)
  dna_dist <- phangorn::dist.ml(v_phyDat, model="JC69")

  dendr <- ggdendro::dendro_data(hclust(dna_dist, method='average'), type="rectangle")

  r <- list()
  r$seqs <- x
  r$segments<- ggdendro::segment(dendr)
  r$labels <- ggdendro::label(dendr)
  r$labels$label2 <- stringr::str_extract(r$labels$label, '\\d+')
  r$labels$label2 <- factor(r$labels$label2)

  r$plot <-
    ggplot() +
    geom_segment(data=r$segments, aes(x=x, y=y, xend=xend, yend=yend)) +
    scale_color_manual(values=c(grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(levels(r$labels$label2))))) +
    geom_text(data=r$labels, aes(x=x, y=y, label=label, color=label2, hjust=0), size=3) +
    labs(x = '', y = '')+
    coord_flip() + scale_y_reverse(expand=c(0.2, 0)) +
    theme(legend.position = "none",
          axis.line.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          panel.background=element_rect(fill="white"),
          panel.grid=element_blank())

  r
}

variantGrid <- function(files, labels){
  
  a <- tibble(file = files,
              sample = unlist(lapply(strsplit(file, '_'), '[', 1)))

  b <- bind_rows(lapply(files, function(x){
    o <- readLines(paste0('variant_calls/', x))
    o <- o[! grepl('^#', o)]
    d <-bind_rows(lapply(o, function(a){
      a <- unlist(strsplit(a, '\t'))
      tibble(file = x, exp = a[1], pos = a[2], v = paste0(a[4], ' > ', paste0(sort(unlist(strsplit(a[5], ','))), collapse = ',')))
    }))
  }))
  
  o <- left_join(a, b, by = 'file')
  o <- o[! is.na(o$exp),]
  o$pos <- as.integer(o$pos)
  o$pos <- factor(o$pos, levels = unique(sort(o$pos)))
  o$v <- factor(o$v)
  o$sample <- factor(o$sample, levels = labels)
  
  ggplot(o, aes(sample, pos, fill = v)) + 
    coord_flip() +
    geom_tile(size=0.25, color = 'black') +
    scale_fill_manual(name = 'Base change', values=c(grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(length(levels(o$v))))) +
    facet_grid(exp~.) +
    scale_y_discrete(expand=c(0,0)) +
    scale_x_discrete(expand=c(0,0)) +
    theme(legend.position = 'none', panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          text = element_text(size=14))
}


# Create collections of concensus sequences (references with variants applied) for tree building.
HLA_AB.double <- Reduce('append', lapply(split(d, d$sample), function(x){
  a <- readDNAStringSet(paste0('consensus_seqs/', x$sample[1], '_HLA-A_A.double.fasta'))
  b <- readDNAStringSet(paste0('consensus_seqs/', x$sample[1], '_HLA-B_C.double.fasta'))
  o <- DNAStringSet(paste0(as.character(a), as.character(b)))
  names(o) <- x$sample[1]
  o
}))

HLA_AB.double.tree <- buildTree(HLA_AB.double)
ggsave(HLA_AB.all.tree$plot, width = 8, units = 'in', file = 'trees_and_grids/allSubjects_HLA-A_HLA-B_bothAlleles_tree.pdf')
labels <- as.character(HLA_AB.double.tree$labels$label)
files <- list.files('variant_calls', pattern = 'HLA-A_A.vcf.double$|HLA-B_C.vcf.double$', full.names = FALSE)
v <- variantGrid(files, labels)
ggsave(v, file = 'trees_and_grids/allSubjects_HLA-A_HLA-B_bothAlleles_grid.pdf', width = 12, height = 6, units = 'in')




HLA_B.double <- Reduce('append', lapply(split(d, d$sample), function(x){
  o <- DNAStringSet(readDNAStringSet(paste0('consensus_seqs/', x$sample[1], '_HLA-B_C.double.fasta')))
  names(o) <- x$sample[1]
  o
}))

HLA_B.double.tree <- buildTree(HLA_B.double)
ggsave(HLA_B.double.tree$plot, width = 8, units = 'in', file = 'trees_and_grids/allSubjects_HLA-B_bothAlleles_tree.pdf')
labels <- as.character(HLA_B.double.tree$labels$label)
files <- list.files('variant_calls', pattern = 'HLA-B_C.vcf.double$', full.names = FALSE)
v <- variantGrid(files, labels)
ggsave(v, file = 'trees_and_grids/allSubjects_HLA-B_bothAlleles_grid.pdf', width = 12, height = 6, units = 'in')


HLA_AB.all <- Reduce('append', lapply(split(d, d$sample), function(x){
  a <- readDNAStringSet(paste0('consensus_seqs/', x$sample[1], '_HLA-A_A.fasta'))
  b <- readDNAStringSet(paste0('consensus_seqs/', x$sample[1], '_HLA-B_C.fasta'))
  o <- DNAStringSet(paste0(as.character(a), as.character(b)))
  names(o) <- x$sample[1]
  o
}))

HLA_AB.all.tree <- buildTree(HLA_AB.all)
ggsave(HLA_AB.all.tree$plot, width = 8, units = 'in', file = 'trees_and_grids/allSubjects_HLA-A_HLA-B_tree.pdf')
labels <- as.character(HLA_AB.all.tree$labels$label)
files <- list.files('variant_calls', pattern = 'HLA-A_A.vcf$|HLA-B_C.vcf$', full.names = FALSE)
v <- variantGrid(files, labels)
ggsave(v, file = 'trees_and_grids/allSubjects_HLA-A_HLA-B_grid.pdf', width = 30, height = 5, units = 'in')


# Subject specific grids.
d$subject <- stringr::str_extract(d$file, '^\\d+')
invisible(lapply(split(d, d$subject), function(x){
  files <- list.files('variant_calls', pattern = 'HLA-A_A.vcf$|HLA-B_C.vcf$', full.names = FALSE)
  files <- files[grepl(paste0(x$sample, collapse = '|'), files)]
  v <- variantGrid(files, unique(x$sample))
  h <- c(2.5, 3.0, 3.5, 4, 4)
  ggsave(v, file = paste('trees_and_grids/subject_', x$subject[1], '_HLA-A_HLA-B_grid.pdf'), height = h[[n_distinct(x$sample)]], width = 20, units = 'in')
}))
