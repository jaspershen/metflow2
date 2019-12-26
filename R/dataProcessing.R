#' @title processData
#' @description Raw MS data processing using xcms.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param path Work directory.
#' @param polarity The polarity of data, "positive"or "negative".
#' @param ppm see xcms.
#' @param peakwidth See xcms.
#' @param snthresh See xcms.
#' @param mzdiff See xcms.
#' @param noise See xcms.
#' @param threads Number of threads.
#' @param output.tic Output TIC plot or not.
#' @param output.bpc Output BPC plot or not.
#' @param output.rt.correction.plot Output rt correction plot or not.
#' @param min.fraction See xcms.
#' @param fill.peaks Fill peaks NA or not.
#' @return Peak table.
#' @export
#' @import xcms 
#' @import MSnbase
#' @import mzR


##debug
# sxtTools::setwd_project()
# setwd("demo_data/POS/")
# path <- "."
# polarity <- "positive"
# ppm <- 25
# peakwidth = c(5, 30)
# snthresh = 10
# mzdiff = -0.001
# noise = 500
# threads = 4
# output.tic = FALSE
# output.bpc = FALSE
# output.rt.correction.plot = FALSE
# min.fraction = 0
# fill.peaks = FALSE



# processData()


processData <- function(path = ".",
                        polarity = c("positive", "negative"),
                        ppm = 25,
                        peakwidth = c(5, 30),
                        snthresh = 10,
                        mzdiff = -0.001,
                        noise = 500,
                        threads = 4,
                        output.tic = TRUE,
                        output.bpc = TRUE,
                        output.rt.correction.plot = TRUE,
                        min.fraction = 0,
                        fill.peaks = FALSE) {
  
  output.path <- file.path(path, "Result")
  dir.create(output.path)
  ##paramters
  parameters <- list(
    path = path,
    polarity = polarity,
    ppm = ppm,
    peakwidth = peakwidth,
    snthresh = snthresh,
    mzdiff = mzdiff,
    noise = noise,
    threads = threads,
    output.tic = output.tic,
    output.bpc = output.bpc,
    output.rt.correction.plot = output.rt.correction.plot,
    min.fraction = min.fraction,
    fill.peaks = fill.peaks
  )
  
  save(parameters, file = file.path(output.path, "parameters"))
  ##peak detection
  
  f.in <- list.files(path = path,
                     pattern = '\\.(mz[X]{0,1}ML|cdf)',
                     recursive = TRUE)
  
  sample_group <-
    unlist(lapply(stringr::str_split(string = f.in, pattern = "/"), function(x) {
      x[1]
    }))
  
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
  
  
  # requireNamespace("xcms")
  cat(crayon::green("Reading raw data, it will take a while...\n"))
  if (any(dir(file.path(path, "Result")) == "raw_data")) {
    load(file.path(path, "Result/raw_data"))
  } else{
    raw_data <- MSnbase::readMSData(
      files = f.in,
      pdata = new("NAnnotatedDataFrame", pd),
      mode = "onDisk",
      verbose = TRUE
    )
    save(raw_data,
         file = file.path(output.path, "raw_data"),
         compress = "xz")
  }
  
  cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
  
  
  #----------------------------------------------------------------------------
  #fix bug
  ## Define the rt and m/z range of the peak area
  # rtr <- c(100, 140)
  # mzr <- c(207.1415 - 25 *  207.1415 / 10 ^ 6, 207.1415 + 25 *  207.1415 / 10 ^ 6)
  # ## extract the chromatogram
  # chr_raw <- xcms::chromatogram(object = raw_data, mz = mzr, rt = rtr)
  # plot(chr_raw, col = group_colors[chr_raw$sample_group])
  # 
  # raw_data %>%
  #   filterRt(rt = rtr) %>%
  #   filterMz(mz = mzr) %>%
  #   plot(type = "XIC")
  # 
  # xchr <- xcms::findChromPeaks(chr_raw, param = CentWaveParam(snthresh = 2))
  # 
  # head(chromPeaks(xchr))
  # chromPeakData(xchr)
  # 
  # sample_colors <- group_colors[xchr$sample_group]
  # layout(1)
  # plot(xchr, 
  #      col = sample_colors,
  #      # peakBg = sample_colors[chromPeaks(xchr)[, "column"]],
  #      peakBg = sample_colors
  #      )
  
  
  
  #----------------------------------------------------------------------------
  cat(crayon::green("Peak detecting...\n"))
  ###peak detection
  cwp <- xcms::CentWaveParam(
    ppm = ppm,
    peakwidth = peakwidth,
    snthresh = snthresh,
    mzdiff = mzdiff,
    noise = noise
  )
  
  if (any(dir(file.path(path, "Result")) == "xdata")) {
    load(file.path(path, "Result/xdata"))
  } else{
    xdata <- xcms::findChromPeaks(
      raw_data,
      param = cwp,
      BPPARAM = BiocParallel::SnowParam(workers = threads,
                                        progressbar = TRUE)
    )
    
    save(xdata,
         file = file.path(output.path, "xdata"),
         compress = "xz")
  }
  
  cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
  
  cat(crayon::green("Correcting rentention time...\n "))
  
  if (any(dir(file.path(path, "Result")) == "xdata2")) {
    load(file.path(path, "Result/xdata2"))
  } else{
    xdata2 <- try(xcms::adjustRtime(xdata,
                                    param = xcms::ObiwarpParam(binSize = 0.6)),
                  silent = TRUE)
  }
  
  cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
  
  if (class(xdata2) == "try-error") {
    xdata2 <- xdata
  } else{
    ## Plot also the difference of adjusted to raw retention time.
    if (output.rt.correction.plot) {
      cat(crayon::green("Drawing RT correction plot..."))
      rt.correction.plot <- plotAdjustedRT(object = xdata2)
      save(
        rt.correction.plot,
        file = file.path(output.path, "rt.correction.plot"),
        compress = "xz"
      )
      ggplot2::ggsave(
        filename = file.path(output.path, "RT correction plot.png"),
        plot = rt.correction.plot,
        width = 20,
        height = 7
      )
      rm(list = c("rt.correction.plot"))
      cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
    }
  }
  
  save(xdata2,
       file = file.path(output.path, "xdata2"),
       compress = "xz")
  
  ###TIC
  if (output.tic) {
    cat(crayon::green("Drawing TIC plot..."))
    tic.plot <- xcms::chromatogram(object = xdata2,
                                   aggregationFun = "sum")
    ## Define colors for different groups
    group_colors <-
      paste0(RColorBrewer::brewer.pal(9, "Set1")[1:length(unique(sample_group))], "60")
    names(group_colors) <- unique(sample_group)
    ## Plot all chromatograms.
    save(tic.plot,
         file = file.path(output.path, "tic.plot"),
         compress = "xz")
    plot <- chromatogramPlot(object = tic.plot, title = "TIC")
    ggplot2::ggsave(
      filename = file.path(output.path, "TIC.png"),
      plot = plot,
      width = 20,
      height = 7
    )
    rm(list = c("plot", "tic.plot"))
    cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
  }
  
  ###BPC
  if (output.bpc) {
    cat(crayon::green("Drawing BPC plot..."))
    bpc.plot <- xcms::chromatogram(object = xdata2,
                                   aggregationFun = "max")
    ## Define colors for different groups
    group_colors <-
      paste0(RColorBrewer::brewer.pal(9, "Set1")[1:length(unique(sample_group))], "60")
    names(group_colors) <- unique(sample_group)
    ## Plot all chromatograms.
    save(bpc.plot,
         file = file.path(output.path, "bpc.plot"),
         compress = "xz")
    plot <- chromatogramPlot(object = bpc.plot, title = "BPC")
    ggplot2::ggsave(
      filename = file.path(output.path, "BPC.png"),
      plot = plot,
      width = 20,
      height = 7
    )
    rm(list = c("plot", "bpc.plot"))
    cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
  }
  
  ## Perform the correspondence
  cat(crayon::green("Grouping peaks across samples...\n"))
  pdp <- xcms::PeakDensityParam(
    sampleGroups = xdata2$sample_group,
    minFraction = min.fraction,
    bw = 20
  )
  xdata2 <- xcms::groupChromPeaks(xdata2, param = pdp)
  cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
  
  if (fill.peaks) {
    ## Filling missing peaks using default settings. Alternatively we could
    ## pass a FillChromPeaksParam object to the method.
    xdata2 <- xcms::fillChromPeaks(xdata2)
  }
  
  cat(crayon::green("Outputting peak table..."))
  ##output peak table
  values <- xcms::featureValues(xdata2, value = "into")
  definition <- xcms::featureDefinitions(object = xdata2)
  definition <- definition[, -ncol(definition)]
  peak.name <- xcms::groupnames(xdata2)
  
  peak.table <- data.frame(peak.name = peak.name,
                           definition,
                           values,
                           stringsAsFactors = FALSE)
  rownames(peak.table) <- NULL
  colnames(peak.table) <-
    stringr::str_replace(
      string = colnames(peak.table),
      pattern = "\\.mz[X]{0,1}ML",
      replacement = ""
    )
  readr::write_csv(peak.table, path = file.path(output.path, "Peak.table.csv"))
  cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
  cat(crayon::bgGreen(clisymbols::symbol$tick ,"All is done!\n"))
}


#' @title chromatogramPlot
#' @description Draw TIC or BPC.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param object Object for tic.plot or bpc.plot.
#' @param title.size Font size of title..
#' @param lab.size Font size of lab title.
#' @param axis.text.size Font size of axis text.
#' @param alpha alpha.
#' @param title Title of the plot..
#' @param interactive interactive.
#' @return A ggplot2 object.
#' @export

setGeneric(
  name = "chromatogramPlot",
  def = function(object,
                 title.size = 15,
                 lab.size = 15,
                 axis.text.size = 15,
                 alpha = 0.5,
                 title = "",
                 interactive = FALSE) {
    options(warn = -1)
    info <- object@phenoData@data
    data <- object@.Data
    rm(list = c("object"))
    data <- apply(data, 2, function(x) {
      x <- x[[1]]
      x <-
        data.frame(
          "mz" = x@rtime,
          "intensity" = x@intensity,
          stringsAsFactors = FALSE
        )
      list(x)
    })
    
    data <- lapply(data, function(x) {
      x[[1]]
    })
    
    data <- mapply(
      FUN = function(x, y, z) {
        x <- data.frame(
          x,
          "group" = y,
          "sample" = z,
          stringsAsFactors = FALSE
        )
        list(x)
      },
      x = data,
      y = info[, 2],
      z = info[, 1]
    )
    
    # data <- lapply(data, function(x){
    #   x <- plyr::dlply(.data = x, .variables = plyr::.(sample))
    # })
    
    data <- do.call(rbind, args = data)
    
    # data <- plyr::dlply(.data = data, .variables = plyr::.(sample))
    
    plot <-
      ggplot2::ggplot(data = data,
                      ggplot2::aes(x = mz, y = intensity)) +
      ggplot2::geom_line(
        data = data,
        mapping = ggplot2::aes(colour = group, group = sample),
        alpha = alpha
      ) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "Retention time (RT, second)", y = "Intensity", title = title) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          color = "black",
          size = title.size,
          face = "plain",
          hjust = 0.5
        ),
        axis.title = ggplot2::element_text(
          color = "black",
          size = lab.size,
          face = "plain"
        ),
        axis.text = ggplot2::element_text(
          color = "black",
          size = axis.text.size,
          face = "plain"
        )
      )
    
    if (interactive) {
      plot <- plotly::ggplotly(plot)
    }
    
    return(plot)
    
  }
)


setGeneric(
  name = "plotAdjustedRT",
  def = function(object,
                 title.size = 15,
                 lab.size = 15,
                 axis.text.size = 15) {
    diffRt <- xcms::rtime(object, adjusted = TRUE) - xcms::rtime(object,
                                                                 adjusted = FALSE)
    diffRt <- split(diffRt, MSnbase::fromFile(object))
    xRt <- xcms::rtime(object,
                       adjusted = TRUE,
                       bySample = TRUE)
    
    sample_name <- object@phenoData@data$sample_name
    sample_group <- object@phenoData@data$sample_group
    
    diffRt <- mapply(
      FUN = function(x, y) {
        list(data.frame(x, y, stringsAsFactors = FALSE))
      },
      x = diffRt,
      y = sample_name
    )
    
    xRt <- mapply(
      FUN = function(x, y) {
        list(data.frame(x, y, stringsAsFactors = FALSE))
      },
      x = xRt,
      y = sample_name
    )
    
    diffRt <- do.call(what = rbind, args = diffRt)
    xRt <- do.call(rbind, xRt)
    
    temp.data <-
      data.frame(xRt, diffRt, stringsAsFactors = FALSE)
    
    colnames(temp.data) <-
      c("rt", "sample.name", "diffRT", "sample.name2")
    rm(list = c("object", "xRt", "diffRt"))
    
    plot <-
      ggplot2::ggplot(data = temp.data, ggplot2::aes(x = rt, y = diffRT)) +
      ggplot2::geom_line(data = temp.data, ggplot2::aes(color = sample.name)) +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "Retention time (second)", y = "RT deviation (second)") +
      ggplot2::theme(
        legend.position = "none",
        axis.title = ggplot2::element_text(
          color = "black",
          size = lab.size,
          face = "plain"
        ),
        axis.text = ggplot2::element_text(
          color = "black",
          size = axis.text.size,
          face = "plain"
        )
      )
    plot
  }
)



setGeneric(
  name = "extractEIC",
  def = function(object,
                 mz.range,
                 rt.range,
                 title.size = 15,
                 lab.size = 15,
                 axis.text.size = 15,
                 alpha = 0.5,
                 title = "") {
    data <- data.frame(
      rt = object@.Data[1, 1][[1]]@rtime,
      intensity = object@.Data[1, 1][[1]]@intensity,
      stringsAsFactors = FALSE
    )
    
    # fit <- MASS::fitdistr(data$intensity, densfun = "normal")
    # temp.data <- rnorm(n = nrow(data), mean = fit$sd[1], sd = fit$sd[2])
    
    plot <-
      ggplot2::ggplot(data = data, ggplot2::aes(x = rt, y = intensity)) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::theme_bw() +
      ggplot2::labs(x = "Retention time (second)", y = "Intensity") +
      ggplot2::theme(
        legend.position = "none",
        axis.title = ggplot2::element_text(
          color = "black",
          size = lab.size,
          face = "plain"
        ),
        axis.text = ggplot2::element_text(
          color = "black",
          size = axis.text.size,
          face = "plain"
        )
      )
    plot
  }
)


#' @title extractPeaks
#' @description From mzXML data extract peaks according to IS table.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param path Work directory.
#' @param ppm see xcms.
#' @param threads Number of threads.
#' @param is.table Peak table. Two columns, column 1 is name of peak, column 2 is m/z of peaks.
#' @return Result contains EIC of peaks.
#' @export
#' @import xcms 
#' @import MSnbase
#' @importFrom MSnbase selectFeatureData
#' @import mzR
#' @import stringr
#' @import tidyverse



extractPeaks <- function(path = ".",
                         ppm = 5,
                         threads = 4,
                         is.table = "is.xlsx"
                         # rt.expand = 1
                         ) {
  output.path <- path
  # dir.create(output.path)
  ##peak detection
  
  f.in <- list.files(path = path,
                     pattern = '\\.(mz[X]{0,1}ML|cdf)',
                     recursive = TRUE)
  
  sample_group <-
    unlist(lapply(stringr::str_split(string = f.in, pattern = "/"), function(x) {
      x[1]
    }))
  
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
  
  # requireNamespace("xcms")
  cat(crayon::green("Reading raw data, it will take a while...\n"))
  
  if (any(dir(path) == "raw_data")) {
    load(file.path(path, "raw_data"))
  } else{
    raw_data <- MSnbase::readMSData(
      files = f.in,
      pdata = new("NAnnotatedDataFrame", pd),
      mode = "onDisk",
      verbose = TRUE
    )
    
    save(raw_data,
         file = file.path(output.path, "raw_data"),
         compress = "xz")
  }
  
  cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
  
  is.table <- readxl::read_xlsx(file.path(path, is.table))
  
  mz <-
    is.table %>%
    dplyr::pull(2)
  
  mz <- as.numeric(mz)
  
  mz_range <-
    lapply(mz, function(x) {
      c(x - ppm * x / 10 ^ 6, ppm * x / 10 ^ 6 + x)
    })
  
  mz_range <- do.call(rbind, mz_range)
  
  cat(crayon::green("Extracting peaks, it will take a while..."))
  peak_data <- xcms::chromatogram(object = raw_data,
                                  mz = mz_range)
  cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
  ##geth the extract peaks table of IS
  # cat(crayon::green("Generating peak table, it will take a while...\n"))
  # peak_table <- xcms::findChromPeaks(object = peak_data, 
  #                      param = xcms::CentWaveParam(peakwidth = c(5, 60), 
  #                                                  snthresh = 2,
  #                                                  mzdiff = -0.001,
  #                                                  noise = 0),
  #                      BPPARAM = BiocParallel::SnowParam(workers = threads,
  #                                                        progressbar = TRUE))
  # 
  # peak_table2 <- try(xcms::adjustRtime(peak_table,
  #                                 param = xcms::ObiwarpParam(binSize = 0.6)),
  #               silent = TRUE)
  # 
  # if(class(peak_table2) == 'try-error'){
  #   peak_table2 <- peak_table 
  # }
  # 
  # pdp <- xcms::PeakDensityParam(
  #   sampleGroups = peak_table2$sample_group,
  #   minFraction = 0,
  #   bw = 20
  # )
  # 
  # peak_table3 <- xcms::groupChromPeaks(peak_table2, param = pdp)
  # cat(crayon::green("Outputting peak table...\n"))
  # ##output peak table
  # values <- xcms::featureValues(peak_table3, value = "into")
  # definition <- xcms::featureDefinitions(object = peak_table3)
  # definition <- definition[, -ncol(definition)]
  # peak.name <- xcms::groupnames(peak_table3)
  # 
  # peak.table <- data.frame(peak.name = peak.name,
  #                          definition,
  #                          values,
  #                          stringsAsFactors = FALSE)
  # rownames(peak.table) <- NULL
  # colnames(peak.table) <-
  #   stringr::str_replace(
  #     string = colnames(peak.table),
  #     pattern = "\\.mz[X]{0,1}ML",
  #     replacement = ""
  #   )
  # readr::write_csv(peak.table, path = file.path(output.path, "Peak.table.csv"))
  
  save(peak_data, file = file.path(output.path, "peak_data"))
  return(peak_data)
}


#' @title showPeak
#' @description Show the peaks from result from extractPeaks function.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param object Object from extractPeaks.
#' @param peak.index Which peak to show. Index.
#' @param title.size Title size.
#' @param lab.size Lab titile size.
#' @param axis.text.size Text size of axis.
#' @param alpha alpha.
#' @param title Title of the plot.
#' @param interactive Interactive or not.

#' @return Result contains EIC of peaks.
#' @export
#' @import xcms 
#' @import MSnbase
#' @import stringr
#' @import tidyverse

setGeneric(
  name = "showPeak",
  def = function(object,
                 peak.index = 1,
                 title.size = 15,
                 lab.size = 15,
                 axis.text.size = 15,
                 alpha = 0.5,
                 title = "",
                 interactive = TRUE) {
    options(warn = -1)
    info <- object@phenoData@data
    data <- object@.Data
    rm(list = c("object"))
    if (peak.index > nrow(data)) {
      peak.index <- nrow(data)
      cat("peak.index is ", nrow(data), '\n')
    }
    data <- apply(data, 2, function(x) {
      x <- x[[peak.index]]
      x <-
        data.frame(
          "rt" = x@rtime,
          "intensity" = x@intensity,
          stringsAsFactors = FALSE
        )
      list(x)
    })
    
    data <- lapply(data, function(x) {
      x[[1]]
    })
    
    data <- mapply(
      FUN = function(x, y, z) {
        x <- data.frame(
          x,
          "group" = y,
          "sample" = z,
          stringsAsFactors = FALSE
        )
        list(x)
      },
      x = data,
      y = info[, 2],
      z = info[, 1]
    )
    
    data <- do.call(rbind, args = data)
    data$intensity[is.na(data$intensity)] <- 0
    
    plot <-
      ggplot2::ggplot(data = data,
                      ggplot2::aes(x = rt, y = intensity)) +
      ggplot2::geom_line(
        data = data,
        mapping = ggplot2::aes(colour = group)
      ) +
      ggplot2::geom_area( mapping = ggplot2::aes(fill = group),
                          alpha = alpha) +
      # ggsci::scale_color_tron(alpha = alpha) +
      # ggsci::scale_fill_tron(alpha = alpha) +
      # ggplot2::scale_y_continuous(breaks = intensity,
      #   labels = ecoflux::scientific_10x(values = intensity, digits = 2)) + 
      # ggplot2::scale_y_continuous(breaks = intensity,
      #                             labels = scales::math_format(10^.intensity)) + 
      
      ggplot2::theme_bw() +
      ggplot2::labs(x = "Retention time (RT, second)", 
                    y = "Intensity", title = title) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          color = "black",
          size = title.size,
          face = "plain",
          hjust = 0.5
        ),
        axis.title = ggplot2::element_text(
          color = "black",
          size = lab.size,
          face = "plain"
        ),
        axis.text = ggplot2::element_text(
          color = "black",
          size = axis.text.size,
          face = "plain"
        )
      )
    
    if (interactive) {
      plot <- plotly::ggplotly(plot)
    }
    return(plot)
  }
)
