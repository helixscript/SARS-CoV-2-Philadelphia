dashBoardPlots <- function(file){
  library(readxl)
  library(plyr)
  library(dplyr)
  library(lubridate)
  library(scales)
  library(ggplot2)
  library(dplyr)
  library(viridis)
  library(hrbrthemes)
  
  start_date <- "2020-01-01"
  lineage_important <- c("B.1.526", "B.1.526.1", "B.1.525", "P.2", "B.1.1.7", "P.1", "B.1.351", "B.1.427", "B.1.429", "B.1.311",  "B.1", "B.1.2", "B.1.243")
  
  # R Script _____________________________________________________
  
  # Read the excel file
  df_original <- read_excel(file)
  
  # Make a new dataframe with only the necessary columns
  df1 <- df_original[,c("sample_date", "lineage", "rationale")]
  
  # Remove missing dates from df1
  df2 <- na.omit(df1)
  
  # Make the date column a date object
  df2$sample_date <- as.Date(df2$sample_date, "%Y-%m-%d")
  
  # Make a dataframe of only the random samples
  df3a <- df2[df2$rationale == 'Vaccine breakthrough',]
  df3b <- df2[df2$rationale == 'S drop out; Vaccine Breakthrough',]
  df3c <- df2[df2$rationale == 'Random; Vaccine Breakthrough',]
  df3d <- df2[df2$rationale == 'S drop out; Vaccine Breakthrough',]
  df3e <- df2[df2$rationale == 'Random; Vaccine Breakthrough, S drop',]
  df4 <- rbind(df3a, df3b, df3c, df3d, df3e)
  
  # Summarize the data
  df4$sample_date_week <- floor_date(df4$sample_date, "month") # Adds a column with date for the week only
  df5 <- ddply(df4, c("sample_date_week", "lineage"), summarise, n = length(lineage)) # Makes a new dataframe with the counts of each lineage for each of the week categories
  df6 <- na.omit(df5) # Remove any data that does not have an associated date
  
  # Make a list of the unique dates and lineages
  lineage_unique <- unique(df6[c("lineage")])
  date_week_unique <- unique(df6[c('sample_date_week')])
  
  # Add important lineages to unique lineages. This makes sure the color assigned to lineages is consistent.
  for(li in 1:length(lineage_important)){
    lineage_unique <- lineage_unique %>% add_row(lineage = lineage_important[li])
  }
  lineage_unique <- unique(lineage_unique[c("lineage")])
  
  # Make the columns of the graph with the repeating lineages for each date
  lineage_repeat <- lineage_unique
  for(d in 2:lengths(date_week_unique)){
    lineage_repeat <- Map(c, lineage_repeat, lineage_unique)
  }
  
  # Make the column of the graph with the repeating date, where there are enough dates for each lineage.
  date_repeat <- list(sample_date = date_week_unique[1,1])
  for(d in 1:lengths(date_week_unique)){
    for(l in 1:lengths(lineage_unique)){
      date_repeat <-Map(c, date_repeat, date_week_unique[d,1])
    }
  } 
  lengths(date_repeat)
  
  # Remove the initial repeated date from the creation of this new list
  date_repeat <- date_repeat$sample_date[-2]
  date_repeat <- list(sample_date = date_repeat)
  lengths(date_repeat)
  
  # Combine the repeated dates and lineages into a single data frame
  df7 <- c(date_repeat, lineage_repeat)
  df7 <- data.frame(df7)
  
  # Make a list of zero counts to be temporary placeholders.  
  count_zero <- rep(c(0),each=length(df7$sample_date))
  count_zero <- list(counts = count_zero)
  df8 <- c(date_repeat, lineage_repeat, count_zero)
  df8 <- data.frame(df8)
  
  # Loop through the data to determine which dates have counts for each lineage and replace zero placeholders with appropriate counts
  short_length <- length(df6$sample_date_week)
  long_length <- length(df8$sample_date)
  for(s in 1:short_length){
    temp_short_date <- df6[s,1]
    temp_short_lineage <- df6[s,2]
    for(l in 1:long_length){
      temp_long_date <- df8[l,1]
      temp_long_lineage <- df8[l,2]
      if(temp_short_date == temp_long_date && temp_short_lineage == temp_long_lineage){
        df8[l,3] <- df6[s,3] 
      }
    }
  }
  
  
  # Apply the date cut off to include reads only as old as the start_date specified
  df9 <- df8 %>%
    select(sample_date, lineage, counts) %>%
    filter(sample_date >= as.Date(start_date))
  
  # Assign data frame
  data <- df9
  
  # Convert variant counts to percentages
  data2 <- data  %>%
    group_by(sample_date, lineage) %>%
    summarise(n = sum(counts)) %>%
    mutate(percentage =100 * n / sum(n),
           date=as.Date(sample_date, "%m/%d/%Y")) %>%
    ungroup() %>% as.data.frame()
  
  # Determine which lineages to drop
  lineage_unique2 <- c("")
  for(lu in 1:lengths(lineage_unique)){
    lineage_unique2 <- c(lineage_unique2, lineage_unique[lu,1])
  }
  lineage_unique2 <- c(lineage_unique2[2:length(lineage_unique2)])
  lineage_drop <- setdiff(lineage_unique2, lineage_important)
  #lineage_drop <- c("NA")
  
  # Lump variants into 'Other' category
  df9 %>%
    as_tibble() %>%
    mutate(date = as.Date(sample_date, "%Y-%m-%d")) %>%
    # Lumped "Other" Variants chosen at random here
    mutate(Lineage = factor(lineage),
           Lineage = forcats::fct_other(f = Lineage, drop = lineage_drop, other_level = "Other")) %>%
    group_by(date, Lineage) %>%
    summarise(Lineage_count = sum(counts, na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(date) %>%
    mutate(Lineage_perc = 100 * Lineage_count / sum(Lineage_count)) %>%
    ungroup() %>%
    identity() -> data3
  
  # Add specimen totals
  data3 %>%
    group_by(date) %>%
    summarise(total_seqs = sum(Lineage_count, na.rm = TRUE)) %>%
    ggplot(data = ., aes(x = date, y = total_seqs)) +
    geom_point() +
    geom_line() +
    theme_ipsum() +
    labs(x = "Date", y = "Genomes Sequenced")
  
  # Add colloquial names
  data4 <- data3
  data4 <- data.frame(lapply(data4, function(x) {gsub("B.1.526.1", "temp_lineage_01", x)}))
  data4 <- data.frame(lapply(data4, function(x) {gsub("B.1.526", "B.1.526 (NY)", x)}))
  data4 <- data.frame(lapply(data4, function(x) {gsub("B.1.1.7", "B.1.1.7 (UK)", x)}))
  data4 <- data.frame(lapply(data4, function(x) {gsub("P.1", "P.1 (Brazil)", x)}))
  data4 <- data.frame(lapply(data4, function(x) {gsub("B.1.351", "B.1.351 (South Africa)", x)}))
  data4 <- data.frame(lapply(data4, function(x) {gsub("B.1.427", "B.1.427 (California)", x)}))
  data4 <- data.frame(lapply(data4, function(x) {gsub("B.1.429", "B.1.429 (California)", x)}))
  data4 <- data.frame(lapply(data4, function(x) {gsub("temp_lineage_01", "B.1.526.1 (NY)", x)}))
  data4$date <- as.Date(data4$date)
  data4$Lineage <- as.factor(data4$Lineage)
  data4$Lineage_count <- as.numeric(data4$Lineage_count)
  data4$Lineage_perc <- as.numeric(data4$Lineage_perc)
  data4$Lineage <- factor(data4$Lineage, levels = c("B.1", "B.1.1.7 (UK)", "B.1.2", "B.1.243", "B.1.311", "B.1.351 (South Africa)", "B.1.427 (California)", "B.1.429 (California)", "B.1.525", "B.1.526 (NY)", "B.1.526.1 (NY)", "P.1 (Brazil)", "P.2", "Other"))
  
  # Plot lineages using bar graph
  
  p1 <- ggplot(data4, aes(fill=Lineage, y=Lineage_perc, x=date)) +
    geom_bar(position="stack", stat="identity") +
    theme_ipsum() +
    ggtitle("Philadelphia Vaccine Breakthrough SARS-CoV-2 Variants") +
    labs(x = "Date", y = "Percentage") +
    scale_x_date(date_labels = "%b %Y", breaks = scales::pretty_breaks(n = 5)) +
    scale_fill_manual(values = c("yellow2", "palegreen", "lightseagreen", "skyblue", "steelblue1", "royalblue", "navy", "purple", "hotpink", "indianred", "red", "orange","lightgoldenrod1", "honeydew4"))
  
  # Plot sequencing counts using bar graph
  p2 <- data4 %>%
    group_by(date) %>%
    summarise(total_seqs = sum(Lineage_count, na.rm = TRUE)) %>%
    ggplot(data = ., aes(x = date, y = total_seqs)) +
    geom_bar(position="stack", stat="identity", width = 5) +
    theme_ipsum() +
    labs(x = "Date", y = "Genomes Sequenced") +
    scale_x_date(date_labels = "%b %Y", breaks = scales::pretty_breaks(n = 3))
  
  # Combine lineages and sequence counts graphs
  library(patchwork)
  p1 + p2 + plot_layout(ncol = 1, heights = c(1,0.4), guides = 'collect') 
  
  
  
  list('mountainPlot' = p1 + p2 + plot_layout(ncol = 1, heights = c(1,0.4)))
}



tmpFile <- function(){ paste0('tmp.', paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = '')) }

removeProblematicSamples <- function(x){
  if('exp' %in% names(x))        x <- dplyr::filter(x, ! grepl('VSP0069', exp))
  if('Subject' %in% names(x))    x <- dplyr::filter(x, Subject != '237')
  if('patient_id' %in% names(x)) x <- dplyr::filter(x, patient_id != '237')
  if('patient_id' %in% names(x)) x <- dplyr::filter(x, ! grepl('MisC', patient_id, ignore.case = TRUE))
  if('trial_id' %in% names(x))   x <- dplyr::filter(x, ! grepl('ODoherty', trial_id, ignore.case = TRUE))
  
  # Use Testing date rather than collection date where available.
  if('sampleCollection_date' %in% names(x) & 'Testing_Date' %in% names(x)){
    x$sampleCollection_date <- mdy(x$sampleCollection_date)
    x$Testing_Date <- mdy(x$Testing_Date)
    x[! is.na(x$Testing_Date),]$sampleCollection_date <- x[! is.na(x$Testing_Date),]$Testing_Date
  }
  
  x
}


variantSchematic <- function(VSPobjPath){
  load(VSPobjPath)
  d <- data.frame(x1 = 0, x2 = nchar(opt$concensusSeq), y1 = 0, y2 = 1)
  m <- data.frame(x1 = opt$variantTableMajor$POS, x2 = opt$variantTableMajor$POS, 
                  y1 = 0, y2 = 1, ALT = as.character(opt$variantTableMajor$ALT))

  m$ALT2 <- ifelse(m$ALT == 'A', 'A', 
                   ifelse(m$ALT == 'T', 'T',
                          ifelse(m$ALT == 'C', 'C',
                                 ifelse(m$ALT == 'G', 'G',
                                        ifelse(grepl('del', m$ALT, ignore.case = T), 'Del',
                                               ifelse(grepl('ins', m$ALT, ignore.case = T), 'Ins', 'Other'))))))

  m$ALT2 <- factor(as.character(m$ALT2), levels = c('A', 'T', 'C', 'G', 'Del', 'Ins', 'Other'))

  p <- ggplot() +
    scale_color_manual(values = c('red', 'blue', 'green2', 'gold2', 'black', 'gray50')) +
    geom_rect(data = d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2), color="black", fill = 'white') +
    geom_segment(data = m, aes(x = x1, y = y1, xend = x2, yend = y2, color = ALT2), size = 1) +
    theme(axis.line=element_blank(),axis.text.x=element_blank(),
          axis.text.y=element_blank(),axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),legend.position="none",
          panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),plot.background=element_blank())

  VSPid <- stringr::str_extract(VSPobjPath, 'VSP\\d+')
  
  ggsave(p, file = file.path('summaries', 'images', paste0(VSPid, '.png')), units = 'in', height = 0.5, width = 3)
  
  o <- subset(samples, VSP == VSPid)
  if(nrow(o) == 1){
    return(tibble('Sample' = VSPid, Date = o$sampleCollection_date, 
                  'Num. mutations' = nrow(opt$variantTableMajor), lineage = as.character(opt$pangolinAssignment), class = o$Rationale,
                  'Analysis report' = paste0('<a href="', file.path('summaries', 'trials', o$trial_id, paste0(o$patient_id, '.pdf')), '">link</a>'),
                  'Mutation positions' = paste0("<img src='",  paste0('summaries/images/', VSPid, '.png'), "'>")))
  } else {
    return(tibble())
  }
}



# Select a sample to represent each subject, sample type, time point combination.
# Select composite samples when available otherwise select the samples with the greated coverage.
# Break ties with read coverage percentages.

representativeSampleSummary <- function(summary, minPercentRefReadCoverage5){
  bind_rows(lapply(split(summary, paste(summary$Subject, summary$sampleType, summary$sampleDate2)), function(x){
    r <- subset(x, type == 'composite')
    if(nrow(r) > 0){
      r <- top_n(r, 1, wt = percentRefReadCoverage5) %>% dplyr::slice(1)
    } else{
      r <- top_n(x, 1, wt = percentRefReadCoverage5) %>% dplyr::slice(1)
    }
    r
  })) %>% dplyr::filter(percentRefReadCoverage5 >= minPercentRefReadCoverage5)
}


# Retrieve consensus sequences from representative samples.
retrieveConcensusSeqs <- function(summary){
  Reduce('append', lapply(1:nrow(summary), function(x){
    x <- summary[x,]
    f <- paste0(softwareDir, '/summaries/VSPdata/', x$exp, '.', ifelse(x$type == 'composite', 'composite', 'experiment'), '.RData')
    if(! file.exists(f)){
      #stop('Can not locate VSP data file -- ', f)
      message('Can not locate VSP data file -- ', f)
      return(DNAStringSet())
    }
    load(f)
    d <- DNAStringSet(opt$concensusSeq)
    
    names(d) <- gsub('\\s+', '_', paste0(x$trial_id, '|', x$Subject, '|', x$sampleType, '|', 
                                         x$sampleDate2, '|', stringr::str_extract(x$exp, '^VSP\\d+'), '|', 
                                         ceiling(mean(opt$pileupData$V4)), 'x mean coverage|', as.character(opt$pangolinAssignment)))
    d
  }))
}