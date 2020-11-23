R1 <- 'trimmedSeqData/228NPOP20200529B_S24_L001_R1_001.fastq.gz'
R2 <- 'trimmedSeqData/228NPOP20200529B_S24_L001_R2_001.fastq.gz'
reference <- 'db/HLA-A.B'
#reference <- 'db/HLA-A.ccds'

bwaPath <- '~/ext/bwa'
samtoolsBin <- '~/ext/samtools/bin'
bcftoolsBin <- '~/ext/bcftools/bin'

system(paste0(bwaPath, ' mem -t 30 -M ', reference, ' ', R1, ' ', R2, ' > o.sam'))
system(paste0('~/ext/samclip --ref ', reference, ' o.sam > o.sam2'))
system(paste0(samtoolsBin, '/samtools view -S -b o.sam2  > o.bam'))
invisible(file.remove(c('o.sam', 'o.sam2')))
system(paste0(samtoolsBin, '/samtools view -q 30 -f 0x2 -b o.bam > o.filt.bam'))
invisible(file.remove('o.bam'))
system(paste0(samtoolsBin, '/samtools sort -o o.filt.sorted.bam o.filt.bam'))
invisible(file.remove('o.filt.bam'))
system(paste0(samtoolsBin, '/samtools index o.filt.sorted.bam'))

