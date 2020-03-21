# setwd("E:/project/hPOP/hPOP-HILIC-POS-MS2/MS1")
# processData(
#   path = ".",
#   polarity = "negative",
#   ppm = 25,
#   peakwidth = c(5, 30),
#   noise = 5000,
#   threads = 6
# )
# 
# load("Result/raw_data")
# peak_table <- readr::read_csv("Result/Peak.table.csv")
# 
# #####mz distribution
# library(tidyverse)
# 
# mz <- peak_table$mzmed
# mz <- sort(mz)
# mz <- as.data.frame(mz)
# 
# cut.point <- round(seq(1, nrow(mz), length.out = 5))
# cut.mz <- mz[cut.point, ]
# 
# 
# round(nrow(mz) / 4)
# 
# library(ggplot2)
# 
# plot <-
#   ggplot(mz, aes(x = mz)) +
#   geom_histogram(binwidth = 1) +
#   theme_bw() +
#   geom_vline(xintercept = cut.mz, colour = "red") +
#   annotate(
#     geom = "text",
#     x = cut.mz,
#     y = 20,
#     label = round(cut.mz, 4),
#     colour = "red"
#   )
# 
# export::graph2ppt(plot, "segment", width = 8, height = 6)
# 
# 
# 
# 
# peak_table <- peak_table %>%
#   filter(., !is.na(QC_N.1)) %>%
#   filter(., rank(desc(QC_N.1)) <= 50)
# # write.csv(peak_table, "Peak_table.csv")
# # colnames(peak_table)
# 
# mz <- peak_table$mzmed
# rt <- peak_table$rtmed
# 
# mz_range <- lapply(mz, function(x) {
#   c(x - 12.5 * ifelse(x < 400, 400, x) / 10 ^ 6,
#     x + 12.5 * ifelse(x < 400, 400, x) / 10 ^ 6)
# })
# 
# rt_range <- lapply(rt, function(x) {
#   c(x - 100, x + 100)
# })
# 
# for (i in 1:50) {
#   cat(i, " ")
#   temp <- chromatogram(raw_data, mz = mz_range[[i]],
#                        rt = rt_range[[i]])
#   png(
#     filename = paste(i, "png", sep = "."),
#     width = 8,
#     height = 6,
#     units = "in",
#     res = 300
#   )
#   plot(temp)
#   dev.off()
# }
# 
# 
# 
# 
# 
# 
# setwd("E:/project/hPOP/hPOP-HILIC-POS-MS2/NCE25/")
# ms2.file <- grep("mzXML", dir(), value = TRUE)
# ms2.data <- readMZXML(file = ms2.file, threads = 4)
# 
# 
# ms1.info <- lapply(ms2.data, function(x) {
#   x[[1]]
# })
# 
# ms1.info <- do.call(what = rbind, args = ms1.info)
# ms1.info <- as.data.frame(ms1.info)
# rownames(ms1.info) <- NULL
# 
# duplicated.name <- unique(ms1.info$name[duplicated(ms1.info$name)])
# if (length(duplicated.name) > 0) {
#   lapply(duplicated.name, function(x) {
#     ms1.info$name[which(ms1.info$name == x)] <-
#       paste(x, c(1:sum(ms1.info$name == x)), sep = "_")
#   })
# }
# 
# 
# ms1.data.hilic.pos <- readr::read_csv(file = "Peak.table.csv",
#                                       col_types = readr::cols())
# colnames(ms1.data)[1:3] <- c("name", "mz", "rt")
# match.result.hilic.pos <-
#   SXTMTmatch(
#     data1 = ms1.data.hilic.pos[, c(2, 3)],
#     data2 = ms1.info[, c(2, 3)],
#     mz.tol = 25,
#     rt.tol = 10,
#     rt.error.type = "abs"
#   )
# 
# 
# match.result.hilic.pos <- as.data.frame(match.result.hilic.pos)
# 
# 
# 
# 
# 
# library(tidyverse)
# 
# library(plyr)
# temp.data1 <-
#   plyr::dlply(.data = match.result.hilic.neg, .variables = .(Index1))
# temp.data2 <-
#   plyr::dlply(.data = match.result.hilic.pos, .variables = .(Index1))
# temp.data3 <-
#   plyr::dlply(.data = match.result.rplc.pos, .variables = .(Index1))
# temp.data4 <-
#   plyr::dlply(.data = match.result.rplc.neg, .variables = .(Index1))
# 
# count1 <- unlist(lapply(temp.data1, nrow))
# count1 <- unname(count1)
# count1 <- data.frame(count1)
# 
# count2 <- unlist(lapply(temp.data2, nrow))
# count2 <- unname(count2)
# count2 <- data.frame(count2)
# 
# count3 <- unlist(lapply(temp.data3, nrow))
# count3 <- unname(count3)
# count3 <- data.frame(count3)
# 
# count4 <- unlist(lapply(temp.data4, nrow))
# count4 <- unname(count4)
# count4 <- data.frame(count4)
# 
# colnames(count1) <-
#   colnames(count2) <- colnames(count3) <- colnames(count4) <-
#   "count"
# 
# count <- rbind(
#   data.frame(count1, "class" = "HILIC_NEG"),
#   data.frame(count2, "class" = "HILIC.POS"),
#   data.frame(count3, "class" = "RPLC_NEG"),
#   data.frame(count4, "class" = "RPLC_POS")
# )
# 
# library(ggplot2)
# plot1 <- ggplot(data = count) +
#   geom_bar(mapping = aes(x = class, fill = factor(count)),
#            position = "fill") +
#   scale_x_discrete(name = "") +
#   scale_y_continuous(name = "Percentage") +
#   # annotate(geom = "text", x = c(1,2,3), y = c(table(count$count))/2, label = c(table(count$count)),
#   #          color = "white", size = 10) + theme_bw()+
#   theme_bw() +
#   guides(fill = guide_legend(title = "MS2 spectra number/peak"))
# 
# 
# temp.data <- data.frame(
#   "HILIC_NEG" = length(unique(match.result.hilic.neg$Index1)) / nrow(ms1.data.hilic.neg),
#   "HILIC_POS" = length(unique(match.result.hilic.pos$Index1)) / nrow(ms1.data.hilic.pos),
#   "RPLC_NEG" = length(unique(match.result.rplc.neg$Index1)) / nrow(ms1.data.rplc.neg),
#   "RPLC_POS" = length(unique(match.result.rplc.pos$Index1)) / nrow(ms1.data.rplc.pos) ,
#   stringsAsFactors = FALSE
# )
# 
# 
# temp.data <- t(temp.data)
# 
# colnames(temp.data) <- "coverage"
# temp.data <- as.data.frame(temp.data)
# 
# temp.data <- data.frame("Class" = rownames(temp.data),
#                         temp.data,
#                         stringsAsFactors = FALSE)
# 
# temp.data[, 2] <- temp.data[, 2] * 100
# 
# plot2 <-
#   ggplot(temp.data) +
#   geom_bar(aes(y = coverage, x = Class, fill = Class), stat = "identity") +
#   theme_bw()
# 
# 
# export::graph2ppt(plot1, "plot1", width = 8, height = 6)
# export::graph2ppt(plot2, "plot2", width = 8, height = 6)
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
# rt <-
#   rbind(
#     data.frame(
#       RT = rt1,
#       LC = "RPLC",
#       stringsAsFactors = FALSE
#     ),
#     data.frame(
#       RT = rt2,
#       LC = "HILIC",
#       stringsAsFactors = FALSE
#     )
#   )
# 
# 
# library(ggplot2)
# 
# rt <- filter(.data = rt, !is.na(RT))
# plot <- ggplot(data = rt, mapping = aes(x = LC, y = RT)) +
#   geom_boxplot(fill = mypal[7]) +
#   theme_bw() +
#   labs(x = "LC", y = "Peak width (second)") +
#   annotate(
#     geom = "text",
#     x = c(1, 2),
#     y = c(43, 22),
#     label = c("34.6 second", "18.4 second"),
#     col = "white",
#     size = 5
#   )
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
# 
# 
# 
# 
# 
# 
# 
# library(tidyverse)
# peak.table <- readr::read_csv("Peak.table.csv")
# 
# mz <- peak.table$mzmed
# mz <- sort(mz)
# mz <- as.data.frame(mz)
# 
# cut.point <- round(seq(1, nrow(mz), length.out = 5))
# cut.mz <- mz[cut.point, ]
# 
# 
# round(nrow(mz) / 4)
# 
# library(ggplot2)
# 
# ggplot(mz, aes(x = mz)) +
#   geom_histogram(binwidth = 1) +
#   theme_bw() +
#   geom_vline(xintercept = cut.mz, colour = "red") +
#   annotate(
#     geom = "text",
#     x = cut.mz,
#     y = 20,
#     label = round(cut.mz, 4),
#     colour = "red"
#   )
