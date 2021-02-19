#' @title processData
#' @description Raw MS data processing using xcms.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param path Work directory.
#' @param polarity The polarity of data, "positive"or "negative".
#' @param ppm see xcms.
#' @param peakwidth See xcms.
#' @param snthresh See xcms.
#' @param prefilter See xcms.
#' @param fitgauss see xcms.
#' @param integrate see xcms.
#' @param mzdiff see xcms.
#' @param noise See xcms.
#' @param threads Number of threads.
#' @param binSize see xcms.
#' @param bw see xcms.
#' @param output.tic Output TIC plot or not.
#' @param output.bpc Output BPC plot or not.
#' @param output.rt.correction.plot Output rt correction plot or not.
#' @param min.fraction See xcms.
#' @param fill.peaks Fill peaks NA or not.
#' @param output.peak.eic Output some peaks EIC or not. 
#' If you want to output this, please place the 'is.xlsx' in the folder. And three columns, name, mz and rt.
#' @param is.table International standard table name.
#' @param is.mz.tolerance mz tolerance for internal standards. Default is 25 ppm.
#' @param is.rt.tolerance RT tolerance for internal standards. Default is 50 s..
#' @param group.for.figure Which group you want to use to output TIC and BPC and EIC. Default is QC.
#' @return Peak table.
#' @export
#' @import xcms 
#' @importFrom xcms CentWaveParam findChromPeaks adjustRtime ObiwarpParam chromatogram PeakDensityParam groupChromPeaks featureChromatograms groupnames featureDefinitions featureValues fillChromPeaks
#' @importFrom Biobase featureData
#' @import MSnbase
#' @import mzR
#' @import ggsci
#' @import grDevices


# #debug
# sxtTools::setwd_project()
# setwd("test_data/mzxml/")
# # rm(list = ls())
# library(xcms)
# library(MSnbase)
# library(mzR)
# library(tidyverse)
# 
# 
# 
# 
# processData(path = ".",
#             polarity = "positive",
#             ppm = 15,
#             peakwidth = c(5, 30),
#             snthresh = 10,
#             prefilter = c(3, 500),
#             fitgauss = FALSE,
#             integrate = 2,
#             mzdiff = 0.01,
#             noise = 500,
#             threads = 20,
#             binSize = 0.025,
#             bw = 5,
#             output.tic = FALSE,
#             output.bpc = FALSE,
#             output.rt.correction.plot = FALSE,
#             min.fraction = 0.5,
#             fill.peaks = FALSE,
#             output.peak.eic = TRUE,
#             is.table = "is.xlsx",
#             group.for.figure = "QC"
# )

# peak_table1 <- readr::read_csv("Result/Peak_table.csv")
# peak_table2 <- readr::read_tsv("results/result.tsv")
# 
# load("Result/intermediate_data/xdata3")
# 
# plot(peak_table1$mzmax - peak_table1$mzmin)
# 
# plot(peak_table2$mzmax - peak_table2$mzmin)
# 
# 
# plot(peak_table1$rtmax - peak_table1$rtmin)
# 
# plot(c(peak_table2$rtmax - peak_table2$rtmin)*60)
# 
# 
# is_table <- readxl::read_xlsx("is.xlsx")
# is_table$mz
# 
# temp.idx <- which(abs(peak_table1$mzmed - is_table$mz[1])*10^6/is_table$mz[1] < 25 &
#                     abs(peak_table1$rtmed - 123) < 50)
# 
# 
# temp.idx
# peak_table1[temp.idx,]$mzmed
# peak_table1[temp.idx,]$rtmed
# peak_table1[temp.idx,c("QC1.10", "QC2.3", "QCU1", "QCU17")]
# 
# temp.idx
# outputFeatureEIC(object = xdata3,
#                  feature.index = 11970,
#                  interactive.plot = TRUE)
# 
# test <- extractPeaks(path = "extractPeak",
#                                ppm = 15, 
#                                threads = 10,
#                                is.table = "is.xlsx")
# metflow2::showPeak(object = test, peak.index = 6, alpha = 0)

setGeneric(name = "processData", 
           def = function(
             path = ".",
             polarity = c("positive", "negative"),
             ppm = 15,
             peakwidth = c(5, 30),
             snthresh = 10,
             prefilter = c(3, 500),
             fitgauss = FALSE,
             integrate = 2,
             mzdiff = 0.01,
             noise = 500,
             threads = 6,
             binSize = 0.025,
             bw = 5,
             output.tic = TRUE,
             output.bpc = TRUE,
             output.rt.correction.plot = TRUE,
             min.fraction = 0.5,
             fill.peaks = FALSE,
             output.peak.eic = TRUE,
             is.table = "is.table.xlsx",
             is.mz.tolerance = 25,
             is.rt.tolerance = 50,
             group.for.figure = "QC"
           ){
             polarity <- match.arg(polarity)
             output_path <- file.path(path, "Result")
             dir.create(output_path, showWarnings = FALSE)
             intermediate_data_path <- file.path(output_path, "intermediate_data")
             dir.create(intermediate_data_path, showWarnings = FALSE)
            
             ##paramters
             parameters <- list(
               path = path,
               polarity = polarity,
               ppm = ppm,
               peakwidth = peakwidth,
               snthresh = snthresh,
               prefilter = prefilter,
               fitgauss = fitgauss,
               integrate = integrate,
               mzdiff = mzdiff,
               noise = noise,
               threads = threads,
               binSize = binSize,
               bw = bw,
               output.tic = output.tic,
               output.bpc = output.bpc,
               output.rt.correction.plot = output.rt.correction.plot,
               min.fraction = min.fraction,
               fill.peaks = fill.peaks
             )
             
             save(parameters, file = file.path(intermediate_data_path, "parameters"))
             
             ##------------------------------------------------------------------------------------
             #peak detection
             
             f.in <- list.files(path = path,
                                pattern = '\\.(mz[X]{0,1}ML|cdf)',
                                recursive = TRUE)
             
             sample_group <-
               unlist(lapply(stringr::str_split(string = f.in, pattern = "/"), function(x) {
                 x[1]
               }))
             
             sample_group[grep("\\.(mz[X]{0,1}ML|cdf)", sample_group)] <- "Group0"
             
             if(!group.for.figure %in% sample_group){
               group.for.figure2 <- 
                 plyr::count(sample_group) %>% 
                 dplyr::filter(freq == min(freq) & stringr::str_to_lower(x) != "blank") %>% 
                 pull(x) %>% 
                 `[`(1) %>% 
                 as.character()
               cat(crayon::yellow(group.for.figure, "is not in you directory, so set is as ", 
                                  group.for.figure2), "\n")
               group.for.figure <- group.for.figure2 
             }
             
             pd <-
               data.frame(
                 sample_name = sub(
                   basename(f.in),
                   pattern = ".mzXML",
                   replacement = "",
                   fixed = TRUE
                 ),
                 sample_group = sample_group,
                 stringsAsFactors = FALSE
               )
             
             ## Define colors for different groups
             group_colors <-
               paste0(RColorBrewer::brewer.pal(9, "Set1")[1:length(unique(sample_group))], "60")
             # group_colors <-
             #   grDevices::colorRampPalette(ggsci::pal_npg()(10))(100)[1:length(unique(sample_group))]
             names(group_colors) <- unique(sample_group)
             
             # requireNamespace("xcms")
             cat(crayon::green("Reading raw data, it will take a while...\n"))
             
             if (any(dir(intermediate_data_path) == "raw_data")) {
               cat(crayon::yellow("Use old saved data in Result.\n"))
               load(file.path(intermediate_data_path, "raw_data"))
             } else{
               raw_data <- MSnbase::readMSData(
                 files = f.in,
                 pdata = new("NAnnotatedDataFrame", pd),
                 mode = "onDisk",
                 verbose = TRUE
               )
               save(raw_data,
                    file = file.path(intermediate_data_path, "raw_data")
                    # compress = "xz"
               )
             }
             
             cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
             
             #----------------------------------------------------------------------------
             cat(crayon::green("Peak detecting...\n"))
             ###peak detection
             cwp <- suppressMessages(
               xcms::CentWaveParam(
               ppm = ppm,
               prefilter = prefilter,
               integrate = integrate,
               peakwidth = peakwidth,
               snthresh = snthresh,
               mzdiff = mzdiff,
               noise = noise, 
               fitgauss = fitgauss,
             )
             )
             
             if (any(dir(intermediate_data_path) == "xdata")) {
               cat(crayon::yellow("Use old saved data in Result.\n"))
               load(file.path(intermediate_data_path, "xdata"))
             } else{
               xdata <- 
                 suppressMessages(
                   try(xcms::findChromPeaks(
                     raw_data,
                     param = cwp,
                     BPPARAM = BiocParallel::SnowParam(workers = threads,
                                                       progressbar = TRUE)
                   ), silent = FALSE
                   )
                 )

               
               if(class(xdata) == "try-error"){
                 stop("Error in xcms::findChromPeaks.\n")
               }
               
               save(xdata,
                    file = file.path(intermediate_data_path, "xdata")
                    # compress = "xz"
               )
             }
             
             rm(list = "raw_data")
             cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
             
             #--------------------------------------------------------------------------------
             #retention time correction
             #Alignment
             cat(crayon::green("Correcting rentention time...\n "))
             
             if (any(dir(intermediate_data_path) == "xdata2")) {
               cat(crayon::yellow("Use old saved data in Result.\n"))
               load(file.path(intermediate_data_path, "xdata2"))
             } else{
               xdata2 <- try(xcms::adjustRtime(xdata,
                                               param = xcms::ObiwarpParam(binSize = 0.5)),
                             silent = FALSE)
               save(xdata2,
                    file = file.path(intermediate_data_path, "xdata2")
                    # compress = "xz"
               )
             }
             
             cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
             
             if (class(xdata2) == "try-error") {
               xdata2 <- xdata
               save(xdata2,
                    file = file.path(intermediate_data_path, "xdata2")
                    # compress = "xz"
               )
             } else{
               ## Plot also the difference of adjusted to raw retention time.
               if (output.rt.correction.plot) {
                 cat(crayon::green("Drawing RT correction plot..."))
                 rt.correction.plot <- plotAdjustedRT(object = xdata2)
                 # save(
                 #   rt.correction.plot,
                 #   file = file.path(intermediate_data_path, "rt.correction.plot"),
                 #   compress = "xz"
                 # )
                 ggplot2::ggsave(
                   filename = file.path(output_path, "RT correction plot.png"),
                   plot = rt.correction.plot,
                   width = 20,
                   height = 7
                 )
                 rm(list = c("rt.correction.plot"))
                 cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
               }
             }
             
             rm(list = "xdata")
             
             ###TIC
             if (output.tic) {
               cat(crayon::green("Drawing TIC plot..."))
               tic.plot <- xcms::chromatogram(object = xdata2,
                                              aggregationFun = "sum",
                                              BPPARAM =
                                                BiocParallel::SnowParam(workers = threads,
                                                                        progressbar = TRUE)
               )
               
               ## Plot all chromatograms.
               # save(tic.plot,
               #      file = file.path(intermediate_data_path, "tic.plot"),
               #      compress = "xz")
               plot <- chromatogramPlot(object = tic.plot,
                                        title = "TIC", 
                                        group.for.figure = group.for.figure)
               
               if(!is.null(plot)){
                 ggplot2::ggsave(
                   filename = file.path(output_path, "TIC.png"),
                   plot = plot,
                   width = 20,
                   height = 7
                 )  
               }
               
               rm(list = c("plot", "tic.plot"))
               cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
             }
             
             ###BPC
             if (output.bpc) {
               cat(crayon::green("Drawing BPC plot..."))
               bpc.plot <- xcms::chromatogram(object = xdata2,
                                              aggregationFun = "max",
                                              BPPARAM =
                                                BiocParallel::SnowParam(workers = threads,
                                                                        progressbar = TRUE))
               
               ## Plot all chromatograms.
               # save(bpc.plot,
               #      file = file.path(intermediate_data_path, "bpc.plot"),
               #      compress = "xz")
               
               plot <- chromatogramPlot(object = bpc.plot, title = "BPC", 
                                        group.for.figure = group.for.figure)
               if(!is.null(plot)){
                 ggplot2::ggsave(
                   filename = file.path(output_path, "BPC.png"),
                   plot = plot,
                   width = 20,
                   height = 7
                 )      
               }
               
               rm(list = c("plot", "bpc.plot"))
               cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
             }
             
             
             #-------------------------------------------------------------------------------------
             ## Perform the correspondence
             cat(crayon::green("Grouping peaks across samples...\n"))
             
             if (any(dir(intermediate_data_path) == "xdata3")) {
               cat(crayon::yellow("Use old saved data in Result.\n"))
               load(file.path(intermediate_data_path, "xdata3"))
             } else{
               pdp <- xcms::PeakDensityParam(
                 sampleGroups = xdata2$sample_group,
                 minFraction = min.fraction,
                 bw = bw, 
                 binSize = binSize, 
                 minSamples = 1, 
                 maxFeatures = 100
               )
               
               xdata3 <- xcms::groupChromPeaks(xdata2, param = pdp)
               save(xdata3,
                    file = file.path(intermediate_data_path, "xdata3")
                    # compress = "xz"
               )
             }
             
             cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
             rm(list = "xdata2")
             
             if (fill.peaks) {
               ## Filling missing peaks using default settings. Alternatively we could
               ## pass a FillChromPeaksParam object to the method.
               xdata3 <- xcms::fillChromPeaks(xdata3)
               save(xdata3,
                    file = file.path(intermediate_data_path, "xdata3")
                    # compress = "xz"
               )
             }
             
             cat(crayon::green("Outputting peak table...\n"))
             ##output peak table
             values <- xcms::featureValues(xdata3, value = "into")
             definition <- xcms::featureDefinitions(object = xdata3)
             definition <- definition[, -ncol(definition)]
             definition <- 
               definition[names(definition) != "peakidx"]
             definition <- 
               definition@listData %>% 
               do.call(cbind, .) %>% 
               as.data.frame()
             
             peak_name <- xcms::groupnames(xdata3)
             peak_name <- paste(peak_name, ifelse(polarity == "positive", "POS", "NEG"), sep = "_")
             
             colnames(values) <-
               stringr::str_replace(
                 string = colnames(values),
                 pattern = "\\.mz[X]{0,1}ML",
                 replacement = ""
               )
             
             peak_table <- data.frame(peak.name = peak_name,
                                      definition,
                                      values,
                                      stringsAsFactors = FALSE)
             rownames(peak_table) <- NULL
             
             colnames(peak_table) <-
               stringr::str_replace(
                 string = colnames(peak_table),
                 pattern = "\\.mz[X]{0,1}ML",
                 replacement = ""
               )
            
             readr::write_csv(peak_table, path = file.path(output_path, "Peak_table.csv"))
              
             peak_table_for_cleaning <-
               definition %>%
               dplyr::select(-c("mzmin", 'mzmax', 'rtmin', 'rtmax', 'npeaks', unique(sample_group))) %>% 
               dplyr::rename(mz = mzmed, rt = rtmed) %>% 
               data.frame(name = peak_name, ., values, stringsAsFactors = FALSE)
             
            
             readr::write_csv(peak_table_for_cleaning, path = file.path(output_path, "Peak_table_for_cleaning.csv"))
             cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
             
             rm(list = c("peak_table", "peak_table_for_cleaning"))
             cat(crayon::red("OK\n"))
             ##-------------------------------------------------------------------------------------
             ##output EIC of all peaks
             
             is_table <- 
               try(readxl::read_xlsx(file.path(path, is.table)), silent = FALSE)
             
             if(class(is_table) != "try-error")  {
               data1 <- as.matrix(is_table[,c(2,3)])
               data2 <- as.matrix(definition[,c("mzmed", "rtmed")])
               
               match_result <- 
                 SXTMTmatch(data1 = data1, 
                            data2 = data2,
                            mz.tolerance = is.mz.tolerance,
                            rt.tolerance = is.rt.tolerance)  
               if(is.null(match_result) | nrow(match_result) == 0){
                 output.peak.eic <- FALSE
               }
             }else{
               output.peak.eic <- FALSE
             }
             
             
             if(output.peak.eic){
               cat(crayon::green("Outputting peak EICs...\n"))
               feature_EIC_path <- file.path(output_path, "feature_EIC")
               dir.create(feature_EIC_path, showWarnings = FALSE)
               
               
               temp_fun <- function(idx = 100,
                                    feature_eic_data,
                                    path = ".",
                                    peak.name,
                                    metabolite.name) {
                 
                 # suppressMessages(require(magrittr))
                 peak.name <- peak.name[idx]
                 plot.name <- paste("",peak.name, sep = "")
                 metabolite.name <- metabolite.name[idx]
                 feature_eic_data <- 
                   feature_eic_data[[idx]]  
                 
                 rt_range <- c(min(feature_eic_data$rt, na.rm = TRUE), 
                               max(feature_eic_data$rt, na.rm = TRUE))
                 
                 if(nrow(feature_eic_data) != 0){
                   plot <- 
                     ggplot2::ggplot(feature_eic_data, ggplot2::aes(rt, intensity, group = sample_name)) +
                     ggplot2::geom_line(ggplot2::aes(color = sample_name)) +
                     ggsci::scale_color_lancet() +
                     ggplot2::labs(x = "Retention time", 
                                   title = paste("RT range:", rt_range[1], rt_range[2], metabolite.name, sep = "_")) +
                     ggplot2::theme_bw()
                   
                   ggplot2::ggsave(plot, 
                                   file = file.path(path, paste(plot.name, "png", sep = ".")),
                                   width = 10, 
                                   height = 6) 
                 }
               }
               
               
               index2 <- sort(unique(match_result[, 2]))
               metabolite_name <- is_table$name[match_result[match(index2, match_result[,2]), 1]]
               feature_eic <-
                 xcms::featureChromatograms(
                   x = xdata3,
                   features = index2,
                   expandRt = 0,
                   BPPARAM =
                     BiocParallel::SnowParam(workers = threads,
                                             progressbar = TRUE)
                 )
               
               
               
               feature_eic_data <- feature_eic@.Data %>% 
                 pbapply::pbapply(1, function(y){
                   y <- lapply(y, function(x){
                     if(class(x) == "XChromatogram"){
                       if(nrow(x@chromPeaks) == 0){
                         data.frame(rt.med = NA,
                                    rt.min = NA,
                                    rt.max = NA,
                                    rt = NA, 
                                    min.intensity = 0,
                                    max.intensity = NA,
                                    intensity = NA,
                                    stringsAsFactors = FALSE) 
                       }else{
                         if(nrow(x@chromPeaks) > 1){
                           x@chromPeaks <- 
                             tibble::as_tibble(x@chromPeaks) %>%
                             dplyr::filter(maxo == max(maxo)) %>% 
                             as.matrix()
                         }
                         data.frame(rt.med = x@chromPeaks[,4],
                                    rt.min = x@chromPeaks[,5],
                                    rt.max = x@chromPeaks[,6],
                                    rt = x@rtime, 
                                    min.intensity = 0,
                                    max.intensity = x@chromPeaks[,"maxo"],
                                    intensity = x@intensity,
                                    stringsAsFactors = FALSE)  
                       } 
                     }else{
                       
                       
                     }
                   }
                   )
                   y <- 
                     mapply(function(y, sample.group, sample.name){
                       data.frame(y, 
                                  sample_group = sample.group,
                                  sample_name = sample.name,
                                  stringsAsFactors = FALSE) %>% 
                         list()
                     },
                     y = y,
                     sample.group = feature_eic@phenoData@data$sample_group,
                     sample.name = feature_eic@phenoData@data$sample_name
                     )
                   
                   y <- do.call(rbind, y)
                   y
                   
                 })
               
               
               feature_eic_data <- 
                 feature_eic_data %>% 
                 lapply(function(x){
                   x <- 
                     x %>% 
                     dplyr::filter(sample_group %in% group.for.figure)
                   
                   if(length(unique(x$sample_name)) > 8){
                     idx <- which(x$sample_name %in% sort(sample(unique(x$sample_name), 18))) %>% 
                       sort()
                     x <- x[idx, , drop = FALSE]
                   }
                   x
                 })
               
               
               BiocParallel::bplapply(1:length(index2), 
                                      FUN = temp_fun,
                                      BPPARAM = BiocParallel::SnowParam(workers = threads,
                                                                        progressbar = TRUE),
                                      feature_eic_data = feature_eic_data, 
                                      path = feature_EIC_path, 
                                      peak.name = peak_name[index2],
                                      metabolite.name = metabolite_name
               )
               
             }
             
             rm(list = "xdata3")
             
             cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
             cat(crayon::bgRed(clisymbols::symbol$tick ,"All done!\n"))
  
})




