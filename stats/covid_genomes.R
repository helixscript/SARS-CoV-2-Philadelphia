library(tidyverse)
library(vegan)
library(reshape2)

load(file = "/home/everett/projects/SARS-CoV-2-Philadelphia/stats/treeData.RData")
dna_dist <- dist_setNames(dna_dist, gsub(" ","_",gsub("-","_",gsub("\\|","_",labels(dna_dist)))))
dendr$labels$label <- gsub(" ","_",gsub("-","_",gsub("\\|","_",dendr$labels$label)))

which(unique(nextClade$seqName) %in% labels(dna_dist))

metadata <- read.table("/home/everett/projects/SARS-CoV-2-Philadelphia/stats/metadata.csv",
                       sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(zip=if_else(Zip.code<10000,paste0("0",as.character(Zip.code)),as.character(Zip.code)),
         Outcome=toupper(Outcome))

dnaDistDf <- melt(as.matrix(dna_dist), varnames =c("row","col")) %>%
  mutate(rowP=substr(row,1,3),
         colP=substr(col,1,3)) %>%
  filter(rowP != colP) %>%
  group_by(rowP, colP) %>%
  summarize(value=min(value)) %>%
  spread(colP, value) %>%
  filter(rowP != "CCL" & rowP != "E6_") %>%
  select(-CCL, -E6_) %>%
  as.data.frame()
rownames(dnaDistDf) <- dnaDistDf$rowP
patient_dna_dist <- as.dist(as.matrix(dnaDistDf %>% select(-rowP)))

treeZips <- metadata$zip[match(as.integer(labels(patient_dna_dist)), metadata$Subject.ID)]
adonis(patient_dna_dist ~ treeZips, method = "bray", perm = 99999)

treeOutcomes <- metadata$Outcome[match(as.integer(labels(patient_dna_dist)), metadata$Subject.ID)]
adonis(patient_dna_dist ~ treeOutcomes, method = "bray", perm = 99999)

treeWHO <- metadata$Max.WHO.score[match(as.integer(labels(patient_dna_dist)), metadata$Subject.ID)]
adonis(patient_dna_dist ~ treeWHO, method = "bray", perm = 99999)

treeAge <- metadata$Age[match(as.integer(labels(patient_dna_dist)), metadata$Subject.ID)]
adonis(patient_dna_dist ~ treeAge, method = "bray", perm = 99999)

nextClade <- read.table("/home/everett/projects/SARS-CoV-2-Philadelphia/stats/nextclade.csv",
                        sep = ";", header = TRUE)

nextCladeSnps <- nextClade %>%
  select(seqName, substitutions) %>%
  separate(substitutions, paste0("mut",seq(101,117)), sep = ",") %>%
  gather("sname","substitution",2:18) %>%
  select(seqName, substitution) %>%
  filter(!is.na(substitution)) %>%
  mutate(t=TRUE) %>%
  spread(substitution, t, fill = FALSE) %>%
  mutate(patient=substr(seqName,1,3)) %>%
  filter(!patient %in% c("CCL","E6_")) %>%
  select(-seqName) %>%
  group_by(patient) %>%
  summarize_all(any) %>%
  left_join(metadata %>% mutate(patient=as.character(Subject.ID), live=Outcome=="LIVE ") %>%
                                   select(patient, live), by="patient") %>%
  as.data.frame()

snps <-
  sapply(colnames(nextCladeSnps)[2:(ncol(nextCladeSnps)-1)], function(col) {
    fishtest <-
      fisher.test(matrix(c(sum(nextCladeSnps[,col] & nextCladeSnps$live),
                           sum(nextCladeSnps[,col] & !nextCladeSnps$live),
                           sum(!nextCladeSnps[,col] & nextCladeSnps$live),
                           sum(!nextCladeSnps[,col] & !nextCladeSnps$live)),
                         ncol=2))
    return(fishtest$p.value)
  })

sig_snps <- snps[snps==min(snps)]

col <- names(sig_snps)[1]
matrix(c(sum(nextCladeSnps[,col] & nextCladeSnps$live),
         sum(nextCladeSnps[,col] & !nextCladeSnps$live),
         sum(!nextCladeSnps[,col] & nextCladeSnps$live),
         sum(!nextCladeSnps[,col] & !nextCladeSnps$live)),
       ncol=2)

genes <- read.table("/home/everett/projects/SARS-CoV-2-Philadelphia/stats/genes.csv",
                    sep = ",", header = TRUE) %>%
  select(gene, start) %>%
  rename(geneStart=start)

nextCladeAA <- nextClade %>%
  select(seqName, aminoacidChanges) %>%
  separate(aminoacidChanges, paste0("mut",seq(101,117)), sep = ",") %>%
  gather("aname","AAchange",2:18) %>%
  select(seqName, AAchange) %>%
  filter(!is.na(AAchange)) %>%
  separate(AAchange, c("gene", "AA"), sep = ":") %>%
  mutate(AAposition=as.integer(substr(AA,2,nchar(AA)-1))) %>%
  left_join(genes, by="gene") %>%
  mutate(NTposition=geneStart + 3 * AAposition)

nextCladeSubs <- nextClade %>%
  select(seqName, substitutions) %>%
  separate(substitutions, paste0("mut",seq(101,117)), sep = ",") %>%
  gather("sname","substitution",2:18) %>%
  select(seqName, substitution) %>%
  filter(!is.na(substitution)) %>%
  mutate(position=as.integer(substr(substitution,2,nchar(substitution)-1)))

subsList <- nextCladeSubs %>%
  arrange(position) %>%
  pull(position)
AAList <- c(nextCladeAA$NTposition,nextCladeAA$NTposition-1,nextCladeAA$NTposition-2)

synonymousSubs <- nextCladeSubs %>%
  filter(!position %in% AAList) %>%
  pull(substitution) %>%
  unique()

nextCladeSubs$synonymous <- nextCladeSubs$substitution %in% synonymousSubs

nextCladeSubs %>% filter(substitution %in% names(sig_snps)) %>%
  arrange(position, seqName)

nextCladeAA %>% filter(NTposition==18999)
