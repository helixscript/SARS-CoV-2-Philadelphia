library(lubridate)

samples <- readr::read_delim('/data/SARS-CoV-2/samples.tsv', '\t')

d1 <- mdy(samples$sampleCollection_date)
d2 <- mdy(samples$Testing_Date)


i <- unique(c(which(abs(d1 - d2) > 30, d1 - d2 < 0)))
           
s <- samples[i, c("sample_id", "patient_id", "trial_id", "sampleCollection_date", "sample_type","VSP", "Testing_Date")] 
s$delta_date <- mdy(s$sampleCollection_date) - mdy(s$Testing_Date)

s$sequenced <- s$VSP %in% unique(stringr::str_extract(list.files('summaries/VSPdata'), 'VSP\\d+'))
openxlsx::write.xlsx(s, file = 'dateTrouble.xlsx')