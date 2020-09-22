library(tidyverse)
library(ggplot2)
library(vegan)
library(reshape2)

nextClade <- read.table("/home/kevin/projects/covid/data/nextclade.csv",
                        sep = ";", header = TRUE)

genes <- read.table("/home/kevin/projects/covid/data/genes.csv",
                    sep = ",", header = TRUE) %>%
  select(gene, start) %>%
  rename(geneStart=start)

nextCladeSubs <- nextClade %>%
  select(seqName, substitutions) %>%
  separate(substitutions, paste0("mut",seq(101,117)), sep = ",") %>%
  gather("sname","substitution",2:18) %>%
  select(seqName, substitution) %>%
  filter(!is.na(substitution)) %>%
  mutate(position=as.integer(substr(substitution,2,nchar(substitution)-1)))


nextCladeSubs %>%
  select(seqName,substitution) %>%
  mutate(t=TRUE) %>%
  filter(substitution %in% c("C241T","C3037T","A23403G","G25563T")) %>%
  spread(substitution,t) %>%
  arrange(G25563T, seqName)

nextCladeSubs %>%
  filter(!seqName %in% (nextCladeSubs %>%
                          select(seqName,substitution) %>%
                          mutate(t=TRUE) %>%
                          filter(substitution %in% c("C241T","C3037T","A23403G","G25563T")) %>%
                          pull(seqName)))

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

subsList <- nextCladeSubs %>%
  arrange(position) %>%
  pull(position)
AAList <- c(nextCladeAA$NTposition,nextCladeAA$NTposition-1,nextCladeAA$NTposition-2)

subsList[!subsList %in% AAList]
length(subsList[!subsList %in% AAList])

synonymousSubs <- nextCladeSubs %>%
  filter(!position %in% AAList) %>%
  pull(substitution) %>%
  unique()

nextCladeSubs$synonymous <- nextCladeSubs$substitution %in% synonymousSubs

write.table(nextCladeSubs, file = "/home/kevin/projects/covid/data/scott_subs.csv",
            sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
write.table(nextCladeAA, file = "/home/kevin/projects/covid/data/scott_AA.csv",
            sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

load(file = "/home/kevin/projects/covid/data/treeData.RData")
dna_dist <- dist_setNames(dna_dist, gsub(" ","_",gsub("-","_",gsub("\\|","_",labels(dna_dist)))))
dendr$labels$label <- gsub(" ","_",gsub("-","_",gsub("\\|","_",dendr$labels$label)))

which(unique(nextClade$seqName) %in% labels(dna_dist))

metadata <- read.table("/home/kevin/projects/covid/data/metadata.csv",
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
adonis(patient_dna_dist ~ treeZips, method = "bray", perm = 999)

treeOutcomes <- metadata$Outcome[match(as.integer(labels(patient_dna_dist)), metadata$Subject.ID)]
adonis(patient_dna_dist ~ treeOutcomes, method = "bray", perm = 999)



dnaDistDf <- melt(as.matrix(dna_dist), varnames =c("row","col")) %>%
  mutate(rowP=substr(row,1,3),
         colP=substr(col,1,3)) %>%
  filter(rowP != colP) %>%
  filter(rowP != "CCL" & rowP != "E6_") %>%
  select(-rowP, -colP) %>%
  spread(col, value) %>%
  select(-CCLB_Vero_cells_20200328, -E6_Vero_cells_20200328) %>%
  as.data.frame()
rownames(dnaDistDf) <- dnaDistDf$row
patient_dna_dist <- as.dist(as.matrix(dnaDistDf %>% select(-row)))

treeZips <- metadata$zip[match(as.integer(substr(labels(patient_dna_dist),1,3)), metadata$Subject.ID)]
adonis(patient_dna_dist ~ treeZips, method = "bray", perm = 999)

treeOutcomes <- metadata$Outcome[match(as.integer(labels(patient_dna_dist)), metadata$Subject.ID)]
adonis(patient_dna_dist ~ treeOutcomes, method = "bray", perm = 999)




