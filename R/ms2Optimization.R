# setwd("E:/project/MSMS optimization/QC10_HILIC")
# # dataProcessing(path = ".", polarity = "positive", 
# #                ppm = 25, peakwidth = c(5, 30), noise = 1000)
# 
# load("Result/raw_data")
# peak_table <- readr::read_csv("Result/Peak.table.csv")
# 
# library(tidyverse)
# 
# peak_table <- peak_table %>% 
# filter(., !is.na(QC10_pHILIC.1)) %>% 
#   filter(., rank(desc(QC10_pHILIC.1)) <= 50)
# write.csv(peak_table, "Peak_table.csv")
# colnames(peak_table)
# 
# mz <- peak_table$mzmed
# rt <- peak_table$rtmed
# 
# 
# # temp.data <- data.frame(mz, rt, intensity = log(QC10_pHILIC.1, 10),
# #                         stringsAsFactors = FALSE)
# # 
# # plot <- ggplot(temp.data, aes(x = rt, y = mz, colour = intensity)) +
# #   geom_point()
# # 
# # library(rayshader)
# # plot_gg(plot,multicore=TRUE,width=5,height=5,scale=250)
# 
# 
# mz_range <- lapply(mz, function(x){
#   c(x - 12.5 * ifelse(x < 400, 400, x) / 10^6,
#   x + 12.5 * ifelse(x < 400, 400, x) / 10^6)
# })
# 
# rt_range <- lapply(rt, function(x){
#   c(x - 100, x + 100)
# })
# 
# for(i in 1:50){
#   cat(i, " ")
#   temp <- chromatogram(raw_data, mz = mz_range[[i]], 
#                        rt = rt_range[[i]])
#   png(filename = paste(i,"png", sep = "."), width = 8, height = 6, 
#       units = "in", res = 300)
#   plot(temp)
#   dev.off()
# }
# 
# 
# 
# setwd("E:/project/MSMS optimization/QC10_RPLC")
# dataProcessing(path = ".", polarity = "positive", 
#                ppm = 25, peakwidth = c(5, 30),
#                noise = 1000)
# 
# load("Result/raw_data")
# peak_table <- readr::read_csv("Result/Peak.table.csv")
# 
# library(tidyverse)
# 
# peak_table <- peak_table %>% 
#   filter(., !is.na(QC10_pRPLC.1)) %>% 
#   filter(., rank(desc(QC10_pRPLC.1)) <= 50)
# write.csv(peak_table, "Peak_table.csv")
# colnames(peak_table)
# 
# mz <- peak_table$mzmed
# rt <- peak_table$rtmed
# 
# mz_range <- lapply(mz, function(x){
#   c(x - 12.5 * ifelse(x < 400, 400, x) / 10^6,
#     x + 12.5 * ifelse(x < 400, 400, x) / 10^6)
# })
# 
# rt_range <- lapply(rt, function(x){
#   c(x - 100, x + 100)
# })
# 
# for(i in 1:50){
#   cat(i, " ")
#   temp <- chromatogram(raw_data, mz = mz_range[[i]], 
#                        rt = rt_range[[i]])
#   png(filename = paste(i,"png", sep = "."), width = 8, height = 6, 
#       units = "in", res = 300)
#   plot(temp)
#   dev.off()
# }
# 
# 
# load("E:/project/MSMS optimization/QC10_HILIC/Result/raw_data")
# 
# 
# 
# 
# 
# match.result <- as.data.frame(match.result)
# 
# library(tidyverse)
# 
# library(plyr)
# temp.data <- plyr::dlply(.data = match.result, .variables = .(Index1))
# 
# count <- unlist(lapply(temp.data, nrow))
# count <- unname(count)
# count <- data.frame(count)
# library(ggplot2)
# ggplot(data = count) +
#   geom_bar(mapping = aes(x = count, fill = factor(count))) +
#   annotate(geom = "text", x = c(1,2,3), y = c(table(count$count))/2, label = c(table(count$count)),
#            color = "white", size = 10) + theme_bw()+
#   ggthemes::theme_economist()
# plot(unlist(lapply(temp.data, nrow)))
# 
# 
# 
# 
# getwd()
# 
# setwd("E:/project/MSMS optimization")
# peak.table.rplc <- readr::read_csv("QC10_RPLC/Peak_table.csv")
# rt1 <- peak.table.rplc$`Peak width`
# 
# peak.table.hilic <- readr::read_csv("QC10_HILIC/Peak_table.csv")
# rt2 <- peak.table.hilic$`Peak width`
# 
# 
# rt <- rbind(data.frame(RT = rt1, LC = "RPLC", stringsAsFactors = FALSE),
#             data.frame(RT = rt2, LC = "HILIC", stringsAsFactors = FALSE))
# 
# 
# library(ggplot2)
# 
# rt <- filter(.data = rt, !is.na(RT))
# plot <- ggplot(data = rt, mapping = aes(x = LC, y = RT)) + 
# geom_boxplot(fill = mypal[7]) +
#   theme_bw() +
#   labs(x = "LC", y = "Peak width (second)")+
#   annotate(geom = "text", x = c(1,2), y = c(43,22), 
#            label = c("34.6 second", "18.4 second"), col = "white", size = 5)
#   
# library(ggsci)
# mypal = pal_npg("nrc", alpha = 0.7)(9)
# mypal
# library("scales")
# show_col(mypal)
# 
# plot <- plot + 
#   scale_fill_manual(values = mypal[3:4])
# 
# mean(rt1, na.rm = TRUE)
# mean(rt2, na.rm = TRUE)
# export::graph2ppt(x = plot, file = "plot")
