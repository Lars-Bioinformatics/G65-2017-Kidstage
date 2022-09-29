library(tidyverse)

setwd("/work/G65-2017-Kidstage/Connor/MRDetect_result/variants_final/mrdetect/")

# patients = paste0("G65-P",1:2,"-")
patients = c("G65-P1-6A6A","G65-P10-BA34","G65-P11-DCC0","G65-P12-BDCB","G65-P13-F083","G65-P14-9201",
             "G65-P18-8328","G65-P19-B0D7-1","G65-P19-B0D7-2","G65-P19-B0D7-3","G65-P19-B0D7-4",
             "G65-P2-C99C","G65-P20-8E72","G65-P21-3822","G65-P22-3C23","G65-P23-EBBE","G65-P25-F720",
             "G65-P26-6930","G65-P28-42A6","G65-P28-42A6","G65-P29-6AE9","G65-P30-D09E","G65-P31-45BB",
             "G65-P32-9291-1","G65-P32-9291-2","G65-P32-9291-3","G65-P32-9291-4","G65-P34-7603",
             "G65-P35-7F49","G65-P36-3A8B","G65-P36-3A8B","G65-P37-9CAE","G65-P37-9CAE","G65-P39-7D1F-1",
             "G65-P39-7D1F-2","G65-P39-7D1F-3","G65-P4-413E","G65-P5-F139","G65-P7-C59D","G65-P8-E0B8")
patient = patients[1]
for (patient in patients){
  files = list.files(pattern="_RESULT", recursive=T)
  data = map_dfr(.x = files, .f = function(file){read.csv(file, header = T, nrows = 1)})
  data = setNames(data,c(names(data)[2:6]))
  data = data[,1:5]
  case = data %>% filter(str_detect(row.names(data),patient))
  controls = data %>% filter(!str_detect(row.names(data),patient))
  z_score = case$sites.detected - mean(controls$sites.detected) / sd(controls$sites.detected)
  print(paste0(patient,": ",z_score, " - sites detected: ", case$sites.detected))
}
