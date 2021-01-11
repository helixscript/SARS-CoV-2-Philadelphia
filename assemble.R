library(ShortRead)
library(tidyverse)
library(optparse)
library(genbankr)
library(GenomicRanges)
source('lib/assemble.lib.R')

option_list = list(
  make_option(c("--outputFile"), type="character", default='save.RData', help="path to output RData file", metavar="character"),
  make_option(c("--R1"), type="character", default=NULL, help="comma delimited list of R1 fastq files", metavar="character"),
  make_option(c("--R2"), type="character", default=NULL, help="comma delimited list of R2 fastq files", metavar="character"),
  make_option(c("--refGenomeFasta"), type="character", default='data/references/USA-WA1-2020.fasta', help="reference genome FASTA path", metavar="character"),   
  make_option(c("--refGenomeBWA"), type="character", default='data/references/USA-WA1-2020.fasta', help="reference genome BWA db path", metavar="character"),   
  make_option(c("--refGenomeGenBank"), type="character", default='data/references/USA-WA1-2020.gb', help="reference genome BWA db path", metavar="character"), 
  make_option(c("--minVariantPhredScore"), type="integer", default=20, help="minimum PHRED score allowed for called varinats", metavar="character"),
  make_option(c("--bwaPath"), type="character", default='~/ext/bwa', help="path to bwa binary", metavar="character"), 
  make_option(c("--megahitPath"), type="character", default='~/ext/megahit/bin/megahit', help="path to megahit binary", metavar="character"), 
  make_option(c("--minBWAmappingScore"), type="integer", default=30, help="minimum BWA mapping score", metavar="character"), 
  make_option(c("--samtoolsBin"), type="character", default='~/ext/samtools/bin', help="path to samtools bin", metavar="character"), 
  make_option(c("--bcftoolsBin"), type="character", default='~/ext/bcftools/bin',  help="path to bcftools bin", metavar="character"))

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

opt$errorCode <- 0
opt$errorMessage <- NA
opt$workDir <- tmpFile()


# Define data structures which must be present even if empty for downstream analyses.
opt$variantTable      <- data.frame()
opt$variantTableMajor <- data.frame()
opt$contigs           <- Biostrings::DNAStringSet()

# Testing data
# opt$R1 <- 'data/sequencing/simulated/VSP9977-1_S1_L001_R1_001.fastq.gz'
# opt$R2 <- 'data/sequencing/simulated/VSP9977-1_S1_L001_R2_001.fastq.gz'

if(! 'outputFile' %in% names(opt)) stop('--workDir must be defined.')
if(! 'workDir' %in% names(opt)) stop('--workDir must be defined.')
if(dir.exists(opt$workDir)) stop('Error -- output directory already exists')
dir.create(opt$workDir)
if(! dir.exists(opt$workDir)) stop('Error -- could not create the output directory.')

if(! 'R1' %in% names(opt)) stop('--R1 must be defined.')
if(! 'R2' %in% names(opt)) stop('--R2 must be defined.')

R1s <- unlist(strsplit(opt$R1, ','));  if(! all(file.exists(R1s))) stop('All the R1 files could not be found.')
R2s <- unlist(strsplit(opt$R2, ','));  if(! all(file.exists(R2s))) stop('All the R1 files could not be found.')

t1 <- paste0(opt$workDir, '/x')

# Combine R1 and R2 data files in composite R1 and R2 files.
system(paste0('cat ', paste0(R1s, collapse = ' '), ' > ', t1, '_R1.fastq'))
system(paste0('cat ', paste0(R2s, collapse = ' '), ' > ', t1, '_R2.fastq'))
  

# Quality trim reads and create trimmed FASTA files.
r <- prepareTrimmedReads(readFastq(paste0(t1, '_R1.fastq')), readFastq(paste0(t1, '_R2.fastq')))
writeFasta(r[[1]], file = paste0(t1, '_R1.trimmed.fasta'))
writeFasta(r[[2]], file = paste0(t1, '_R2.trimmed.fasta'))


# Align trimmed reads to the reference genome.
system(paste0(opt$bwaPath, ' mem -M ', opt$refGenomeBWA, ' ',  paste0(t1, '_R1.trimmed.fasta'), ' ', 
              paste0(t1, '_R2.trimmed.fasta'),  ' > ', paste0(t1, '_genome.sam')))
system(paste0(opt$samtoolsBin, '/samtools view -S -b ', paste0(t1, '_genome.sam'), ' > ', paste0(t1, '_genome.bam')))
invisible(file.remove(paste0(t1, '_genome.sam')))


# Remove read pairs with mapping qualities below the provided min score and read pairs not properly mapped.
# system(paste0(opt$samtoolsBin, '/samtools view -q ', opt$minBWAmappingScore, ' -f 0x2 -b ', 
#               t1, '_genome.bam > ', t1, '_genome.filt.bam'))

system(paste0(opt$samtoolsBin, '/samtools view -q ', opt$minBWAmappingScore, ' -b ', t1, '_genome.bam > ', t1, '_genome.filt.bam'))


# Retrieve a list of aligned reads.
alignedReadsIDs <- system(paste0(opt$samtoolsBin, '/samtools view ', t1, '_genome.filt.bam | cut  -f 1 | uniq'), intern = TRUE)

if(length(alignedReadsIDs) == 0){
  opt$errorCode <- 1
  opt$errorMessage <- 'No quality trimmed readsa aligned to the reference genome.'
  save(opt, file = opt$outputFile)
  unlink(opt$workDir, recursive = TRUE)
  stop()
}

# Build contigs with reads that aligned to the reference genome in the proper orientation and sufficient mapping scores.
writeFasta(r[[1]][names(r[[1]]) %in% alignedReadsIDs], file = paste0(t1, '_R1.trimmed.genomeAligned.fasta'))
writeFasta(r[[2]][names(r[[2]]) %in% alignedReadsIDs], file = paste0(t1, '_R2.trimmed.genomeAligned.fasta'))

opt$contigs <- megaHitContigs(paste0(t1, '_R1.trimmed.genomeAligned.fasta'), paste0(t1, '_R2.trimmed.genomeAligned.fasta'), 
                              workDir = paste0(t1, '_megahit'), megahit.path = opt$megahitPath)


refGenomeLength <- width(readFasta(opt$refGenomeFasta))

# Align contigs to reference genome and only retain those that map well (alignment flag == 0 or 16).
writeFasta(opt$contigs, file = paste0(t1, '_contigs.fasta'))
system(paste0(opt$bwaPath, ' mem -M ', opt$refGenomeBWA, ' ',  paste0(t1, '_contigs.fasta'), ' > ', paste0(t1, '_contigs.sam')))
sam <- readLines(paste0(t1, '_contigs.sam'))
samLines <- sam[! grepl('^@', sam)]

if(length(samLines) != 0){
  sam <- subset(read.table(textConnection(samLines), sep = '\t', header = FALSE, fill = TRUE), V2 %in% c(0, 16))
  sam$length <- sapply(as.character(sam$V13), samMD2length)
  opt$contigs <- opt$contigs[names(opt$contigs) %in% sam$V1]

  if(length(opt$contigs) > 0){
    d <- group_by(sam, V1) %>% top_n(1, length) %>% dplyr::slice(1) %>% 
         summarise(start = V4, length = length, editDist = as.integer(str_extract(V12, '\\d+'))) %>% ungroup()
  
    opt$contigStartPos   <- d[match(names(opt$contigs), d$V1),]$start
    opt$contigEndPos     <- opt$contigStartPos + d[match(names(opt$contigs), d$V1),]$length
    opt$contigsEditDists <- d[match(names(opt$contigs), d$V1),]$editDist
    names(opt$contigs)   <- paste0(names(opt$contigs), ' [', opt$contigsEditDists, ']')
    
    # Check to make sure that the assembler did not create something much longer than the reference.
    opt$contigs <- opt$contigs[! opt$contigEndPos > refGenomeLength + 50]
    
  } else {
    opt$contigStartPos <- NA
    opt$contigEndPos <- NA
    opt$contigsEditDists <- NA
  }
} else 
{
  opt$contigs <- Biostrings::DNAStringSet()
  opt$contigStartPos <- NA
  opt$contigEndPos <- NA
  opt$contigsEditDists <- NA
}


system(paste0(opt$samtoolsBin, '/samtools sort -o ', paste0(t1, '_genome.filt.sorted.bam'), ' ', paste0(t1, '_genome.filt.bam')))
system(paste0(opt$samtoolsBin, '/samtools index ', paste0(t1, '_genome.filt.sorted.bam')))


# Save bam and bam index files for downstream analyses.
system(paste0('cp ', t1, '_genome.filt.sorted.bam ', sub('.RData', '.bam', sub('VSPdata', 'VSPalignments', opt$outputFile))))
system(paste0('cp ', t1, '_genome.filt.sorted.bam.bai ', sub('.RData', '.bam.bai', sub('VSPdata', 'VSPalignments', opt$outputFile))))


# Create pileup data file for determining depth at specific locations.
system(paste0(opt$samtoolsBin, '/samtools mpileup -A -a --output ', paste0(t1, '.pileup'), ' --max-depth 100000 -f ', 
              opt$refGenomeFasta, ' ', paste0(t1, '_genome.filt.sorted.bam')))


# Determine the percentage of the reference genome covered in the pileup data.
opt$pileupData <- tryCatch({
                              read.table(paste0(t1, '.pileup'), sep = '\t', header = FALSE, quote = '')[,1:5]
                           }, error = function(e) {
                              return(data.frame())
                           })



# Pileup format reports the number of read pairs (column 4) while VCF format (DP) 
# reports the number of reads which appears to report 2x the pileup format value. 
# Confirmed by looking at pileup in IGV.


if(nrow(opt$pileupData) > 0){
  refGenomeLength <- nchar(as.character(readFasta(opt$refGenomeFasta)@sread))
  opt$refGenomePercentCovered <- nrow(opt$pileupData) / refGenomeLength
  opt$refGenomePercentCovered_5reads <- nrow(subset(opt$pileupData,  V4 >= 5))  / refGenomeLength
  opt$refGenomePercentCovered_10reads <- nrow(subset(opt$pileupData, V4 >= 10)) / refGenomeLength
  opt$refGenomePercentCovered_25reads <- nrow(subset(opt$pileupData, V4 >= 25)) / refGenomeLength
  opt$refGenomePercentCovered_50reads <- nrow(subset(opt$pileupData, V4 >= 50)) / refGenomeLength
  
  # If pileup data could be created then we can try to call variants.
  system(paste0(opt$bcftoolsBin, '/bcftools mpileup -A --max-depth 100000 -Ou -f ', opt$refGenomeFasta, ' ',
                paste0(t1, '_genome.filt.sorted.bam'), ' |  ', opt$bcftoolsBin,  '/bcftools call -mv -Oz ', 
                '-o ', paste0(t1, '.vcf.gz')))
  
  
  # Read in the variant table created by bcf tools. 
  # We use tryCatch() here because the table may be empty only containing the header information.
  opt$variantTable <- tryCatch({
      system(paste0(opt$bcftoolsBin, "/bcftools filter -i'QUAL>", opt$minVariantPhredScore, "' ", 
                    paste0(t1, '.vcf.gz'), " -O z -o ", paste0(t1, '.filt.vcf.gz')))
      
      system(paste0(opt$bcftoolsBin, '/bcftools index ', t1, '.filt.vcf.gz'))
      
      x <- read.table(paste0(t1, '.filt.vcf.gz'), sep = '\t', header = FALSE, comment.char = '#')
      names(x) <- c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'OTHER')
      x <- x[, !is.na(names(x))]
      x
    },  error=function(cond) {
      return(data.frame()) 
    })
  
  if(nrow(opt$variantTable) > 0){
    opt$variantTable <- tryCatch({
      # Here we parse the pileup data to create a more informative alt call for variants.
      x <- bind_rows(lapply(1:nrow(opt$variantTable), function(i){
                  x <- opt$variantTable[i,]
                  p <- parsePileUpString(subset(opt$pileupData, V2 == x$POS)$V5)
                  
                  # Expand the variant call to include the different possibilities.
                  x <- x[rep(1, length(p)),]
                  x$ALT <- names(p)
                  x$percentAlt <- p
                  x <- subset(x, percentAlt > 0)
                  x$reads <- subset(opt$pileupData, V2 == x$POS[1])$V4
                  
                  dplyr::select(x, -ID, -FILTER, -INFO, -FORMAT, -OTHER)
                }))
      x
    }, error=function(cond) {
      stop('Error parsing variant occurrences.')
    })
    
    # Select the major variant for each position.
    opt$variantTableMajor <- dplyr::filter(opt$variantTable,  percentAlt >= 0.5 & ! ALT == 'e') %>%
                             dplyr::group_by(POS) %>%
                             dplyr::top_n(1, wt = percentAlt) %>%
                             dplyr::slice(1) %>%
                             dplyr::ungroup() %>%
                             dplyr::filter(reads >= 5)
    
  }
} else {
  opt$refGenomePercentCovered <- 0
  opt$refGenomePercentCovered_5reads  <- 0
  opt$refGenomePercentCovered_10reads <- 0
  opt$refGenomePercentCovered_25reads <- 0
  opt$refGenomePercentCovered_50reads <- 0
  
  opt$errorCode <- 5
  opt$errorMessage <- 'No pileup or variant data available.'
}


# Determine the result of the variants on the AA sequence of viral proteins.
# This approach assumes all orfs are non-overlapping pos. strand sequenes.
if(nrow(opt$variantTableMajor) > 0){

  # Here we copy the variant vcf file and selectively remove calls 
  # which are not found in our filtered variantTable. This is done 
  # to preserve the additional comment lines which appear to be necessary. 

  tryCatch({
    system(paste0('cp ', t1, '.filt.vcf.gz ', t1, '.filt.vcf.copy.gz'))
    system(paste0('gunzip ', t1, '.filt.vcf.copy.gz'))
    
    o <- lapply(readLines(paste0(t1, '.filt.vcf.copy')), function(x){
      if(grepl('^#', x)){
        return(x)
      } else {
        if(as.integer(unlist(strsplit(x, '\t'))[2]) %in% opt$variantTableMajor$POS){
          return(x)
        } else {
          return(NULL)
        }
      }
    })
    
    write(unlist(o[! sapply(o, is.null)]), file = paste0(t1, '.filt.vcf.copy2'))
    system(paste0('bgzip ', t1, '.filt.vcf.copy2'))
    system(paste0(opt$bcftoolsBin, '/bcftools index ', t1, '.filt.vcf.copy2.gz'))
    system(paste0('cat  ', opt$refGenomeFasta, ' | ', opt$bcftoolsBin, '/bcftools consensus ',  
                  t1, '.filt.vcf.copy2.gz > ', t1, '.consensus.fasta'))
  
    opt$concensusSeq <- as.character(readFasta(paste0(t1, '.consensus.fasta'))@sread)
  }, error=function(cond) {
    stop('Error creating concensus sequence.')
  })
  
  gb <- readGenBank(opt$refGenomeGenBank)
  
  cds <- gb@cds
  seqlevels(cds) <- 'genome'
  seqnames(cds)  <- 'genome'
  
  # Calculate how the shift left or right caused by deletions and insertions.
  opt$variantTableMajor$shift <- ifelse(grepl('del', opt$variantTableMajor$ALT), (nchar(opt$variantTableMajor$ALT)-3)*-1, 0)
  opt$variantTableMajor$shift <- ifelse(grepl('ins', opt$variantTableMajor$ALT), (nchar(opt$variantTableMajor$ALT)-3), opt$variantTableMajor$shift)
  
  # Remove variant positions flanking indels since they appear to be artifacts. 
  artifacts <- c(opt$variantTableMajor[grep('ins|del', opt$variantTableMajor$ALT),]$POS + abs(opt$variantTableMajor[grep('ins|del', opt$variantTableMajor$ALT),]$shift),
                 opt$variantTableMajor[grep('ins|del', opt$variantTableMajor$ALT),]$POS -1)
  
  opt$variantTableMajor <- opt$variantTableMajor[! opt$variantTableMajor$POS %in% artifacts,]
  
  opt$variantTableMajor <- bind_rows(lapply(split(opt$variantTableMajor, 1:nrow(opt$variantTableMajor)), function(x){
    
    
    
    # Determine the offset of this position in the concensus sequence because it may not be the same length
    # if indels have been applied. Here we sum the indel shifts before this variant call.
    offset <- sum(opt$variantTableMajor[1:grep(x$POS, opt$variantTableMajor$POS),]$shift)
    
    cds2 <- cds
    start(cds2) <- start(cds2) + offset
    end(cds2) <- end(cds2) + offset
    
    v <- GRanges(seqnames = 'genome', ranges = IRanges(x$POS, end = x$POS), strand = '+')
    o <- GenomicRanges::findOverlaps(v, cds2)
    
    if(length(o) == 0){
      x$genes <- 'intergenic'
      x$type <- ' '
    } else {
      
      # Define the gene the variant is within.
      hit <- cds2[subjectHits(o)]
      x$genes <- paste0(hit$gene, collapse = ', ')
      
      # Retrieve the amino acid sequence of the gene the variant is within.
      orf  <- as.character(translate(DNAString(substr(as.character(readFasta(opt$refGenomeFasta)@sread), start(hit), end(hit)))))
      
      orf2 <- as.character(translate(DNAString(substr(opt$concensusSeq, start(hit), end(hit)))))
      
    
      # Determine the offset of this position in the concensus sequence because it may not be the same length
      # if indels have been applied. Here we sum the indel shifts before this variant call.
      # offset <- sum(opt$variantTableMajor[1:grep(x$POS, opt$variantTableMajor$POS),]$shift)
      
      #              1   2   3   4   5   6   7   8
      # 123 456 789 012 345 678 901 234 567 890 123
      # ATG CAT TGA ATG GGC TTA CGA GCT TAA GTA TAG
      #             ^             x  21-10 + 2 = 13/3 = 4.3 ~ 4
      #                          x   20-10 + 2 = 12/3 = 4.0 = 4
      #                         x    19-10 + 2 = 11/3 = 3.6 ~ 4
      #                                 x   25-10 + 2 = 17/3 = 5.6 ~ 6
      #                                  x  26-10 + 2 = 18/3 = 6.0 = 6
      #                                   x 27-10 + 2 = 19/3 = 6.3 ~ 6
      
      aa <- round(((x$POS - start(hit)) + 2)/3)
      orf_aa <- substr(orf, aa, aa)
      
      aa2 <- round((((x$POS + offset) - start(hit)) + 2)/3)
      orf2_aa <- substr(orf2, aa2, aa2)
      
      maxALTchars <- max(nchar(unlist(strsplit(as.character(x$ALT), ','))))
      
      if(nchar(as.character(x$REF)) == 1 & nchar(as.character(x$ALT)) > 1 & maxALTchars == 1){
        x$type <- paste0(x$POS, '_mixedPop')
      } else if (grepl('ins', as.character(x$ALT))){  
        x$type <- paste0('ins ', nchar(x$ALT)-3)
      } else if (grepl('del', as.character(x$ALT))){
        x$type <- paste0('del ', nchar(x$ALT)-3)
      } else if (orf_aa != orf2_aa){
        x$type <- paste0(orf_aa, aa2, orf2_aa)
      } else {
        x$type <- 'silent'
      }
    }
    x
  }))
  
} else {
  # There were no variants called so we report the reference as the concensus sequence.
  opt$concensusSeq <- as.character(readFasta(opt$refGenomeFasta)@sread)
}


# BCFtools calls indels bythe base preceding the modification.
# Here we find deletion calls and increment the position by one to mark the first deleted base.
i <- opt$variantTable$POS %in% opt$variantTable[grep('del', opt$variantTable$ALT),]$POS
opt$variantTable[i,]$POS <- opt$variantTable[i,]$POS + 1

i <- opt$variantTableMajor$POS %in% opt$variantTableMajor[grep('del', opt$variantTableMajor$ALT),]$POS
opt$variantTableMajor[i,]$POS <- opt$variantTableMajor[i,]$POS + 1



save(opt, file = opt$outputFile)
unlink(opt$workDir, recursive = TRUE)
