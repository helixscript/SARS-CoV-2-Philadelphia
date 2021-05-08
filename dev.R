
dayBreaks <- '8 days'
minSlopeTimePoints <- 3

d <- genomeMetaData[grepl('random|asymptomatic|hospitalized', 
                               genomeMetaData$rationale, ignore.case = T),] %>%
  dplyr::filter(lineage != 'NA') %>%
  dplyr::filter(ymd(sample_date) >= '2021-02-28') %>%
  dplyr::arrange(sample_date) %>%
  dplyr::mutate(date = cut.Date(ymd(sample_date), breaks = dayBreaks)) 

cluster <- makeCluster(10)
clusterExport(cluster, 'representativeSampleSummary_95')
r <- bind_rows(parLapply(cluster, split(d, paste0(d$date, d$lab_id)), function(x){
        load(subset(representativeSampleSummary_95, VSP == x$lab_id)$dataFile)
        data.frame(date = x$date, VSP = x$lab_id[1], mutation = paste(opt$variantTableMajor$POS, opt$variantTableMajor$genes, opt$variantTableMajor$type))
     }))

r2 <- bind_rows(lapply(split(r, r$date), function(x){
        data.frame(date = x$date[1], f = table(x$mutation) / n_distinct(x$VSP)) 
      }))


r3 <- bind_rows(lapply(split(r2, r2$f.Var1), function(x){
        x <- dplyr::arrange(x, date)
        if(nrow(x) >= minSlopeTimePoints){
          m <- (x[nrow(x),]$f.Freq - x[1,]$f.Freq) / as.integer(ymd(x[nrow(x),]$date) - ymd(x[1,]$date))
          return(tibble(mutation = x$f.Var1[1], m = m))
        } else {
         return(tibble())
        }
       })) %>% arrange(desc(m))

colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(21)

ggplot(subset(r2, f.Var1 %in% r3[1:21,]$mutation), aes(date, f.Freq, group = f.Var1, color = f.Var1)) +
 theme_bw()+
   scale_color_manual(name = '', values = colors) +
   geom_point(size = 2) +
  geom_line() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = 'Date', y = 'Sample frequency') +
  guides(color = guide_legend(ncol=3)) +
  theme(legend.position="bottom",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))


o <- subset(d, date == levels(d$date)[length(levels(d$date))])
o <- left_join(o, dplyr::select(representativeSampleSummary_95, VSP, dataFile), by = c('lab_id' = 'VSP'))

r <- data.frame(table(unlist(lapply(split(o, 1:nrow(o)), function(x){
  load(x$dataFile)
  paste(opt$variantTableMajor$POS, opt$variantTableMajor$genes, opt$variantTableMajor$type) 
})))) %>%
  mutate(p = Freq / nrow(o)) %>%
  arrange(desc(p))

r
r$p <= r$Freq / nro