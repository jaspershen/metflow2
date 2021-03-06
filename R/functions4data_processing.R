#' @title plot_chromatogram
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
#' @param group.for.figure What groups to show EIC.
#' @return A ggplot2 object.
#' @export

plot_chromatogram = function(
  object,
  title.size = 13,
  lab.size = 13,
  axis.text.size = 12,
  alpha = 0.5,
  title = "",
  interactive = FALSE,
  group.for.figure = "QC"
){
  options(warn = -1)
  info <- object@phenoData@data
  data <- object@.Data
  data <-
    data[1, which(info$sample_group %in% group.for.figure), drop = FALSE]
  info <-
    info[info$sample_group %in% group.for.figure, , drop = FALSE]
  
  if (nrow(info) == 0) {
    return(NULL)
  }
  
  if (nrow(info) > 10) {
    idx <- sort(sample(1:nrow(info), 10))
    info <- info[idx, , drop = FALSE]
    data <- data[, idx, drop = FALSE]
  }
  
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
    ggplot2::geom_line(data = data,
                       mapping = ggplot2::aes(color = sample, group = sample),
                       alpha = alpha) +
    ggsci::scale_color_aaas() +
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











#' @title plotAdjustedRT
#' @description plotAdjustedRT
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param object object
#' @param title.size title.size
#' @param lab.size lab.size
#' @param axis.text.size axis.text.size
#' @return A ggplot2 object.

plotAdjustedRT = function(
  object,
  title.size = 15,
  lab.size = 15,
  axis.text.size = 15
){
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





#' @title extractEIC
#' @description extractEIC
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param object object
#' @param mz.range mz.range
#' @param rt.range rt.range
#' @param title.size title.size
#' @param lab.size lab.size
#' @param axis.text.size axis.text.size
#' @param alpha alpha
#' @param title title
#' @return result

extractEIC = function(object,
                      mz.range,
                      rt.range,
                      title.size = 15,
                      lab.size = 15,
                      axis.text.size = 15,
                      alpha = 0.5,
                      title = ""){
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



#' @title extractPeaks
#' @description From mzXML data extract peaks according to IS table.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param path Work directory.
#' @param ppm see xcms.
#' @param threads Number of threads.
#' @param is.table Peak table. Two columns, column 1 is name of peak, column 2 is m/z of peaks.
#' @param mz mz
#' @param rt rt
#' @param rt.tolerance Rt tolerance.
#' @return Result contains EIC of peaks.
#' @export

extractPeaks = function(path = ".",
                        ppm = 15,
                        threads = 4,
                        is.table = "is.xlsx",
                        mz = NULL,
                        rt = NULL,
                        rt.tolerance = 40){
  options(warn = -1)
  output_path <- path
  dir.create(output_path, showWarnings = FALSE)
  ##peak detection
  
  f.in <- list.files(
    path = path,
    pattern = '\\.(mz[X]{0,1}ML|cdf)',
    recursive = TRUE,
    full.names = TRUE
  )
  
  sample_group <-
    BiocGenerics::basename(f.in) %>%
    stringr::str_replace("\\.(mz[X]{0,1}ML|cdf)", "")
  
  pd <-
    data.frame(
      basename(f.in),
      sample_group = sample_group,
      stringsAsFactors = FALSE)
  
  cat(crayon::green("Reading raw data, it will take a while...\n"))
  
  if (any(dir(path) == "raw_data")) {
    cat(crayon::yellow("Use old data.\n"))
    load(file.path(path, "raw_data"))
  } else{
    raw_data <- MSnbase::readMSData(
      files = f.in,
      pdata = new("NAnnotatedDataFrame", pd),
      mode = "onDisk",
      verbose = TRUE
    )
    
    save(raw_data,
         file = file.path(output_path, "raw_data"),
         compress = "xz")
  }
  
  cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
  
  is.table <-
    try(readxl::read_xlsx(file.path(path, is.table)), silent = TRUE)
  
  if (!is.null(mz) & !is.null(rt)) {
    if (length(mz) != length(rt)) {
      cat(crayon::yellow("Lenght of mz and rt you provied are different.\n"))
    }
    is.table <- data.frame(mz = as.numeric(mz),
                           rt = as.numeric(rt),
                           stringsAsFactors = FALSE)
    is.table$name <- paste("feature", 1:nrow(is.table), sep = "_")
    
    is.table <-
      is.table %>%
      dplyr::select(name, mz, rt)
  }
  
  if (class(is.table)[1] == "try-error") {
    stop(crayon::red('Please provide right is table or mz and rt.\n'))
  }
  
  mz <-
    is.table %>%
    dplyr::pull(2) %>% 
    as.numeric()
  
  mz_range <-
    lapply(mz, function(x) {
      c(x - ppm * x / 10 ^ 6, ppm * x / 10 ^ 6 + x)
    }) %>% 
    do.call(rbind, .) %>% 
    as.data.frame()
  
  if (any(colnames(is.table) == "rt")) {
    rt <-
      is.table %>%
      dplyr::pull(3) %>%
      as.numeric()
    
    rt_range <-
      lapply(rt, function(x) {
        c(x - rt.tolerance, x + rt.tolerance)
      }) %>%
      do.call(rbind, .)
  } else{
    rt_range <- NA
  }
  
  cat(crayon::green("Extracting peaks, it will take a while..."))
  if (!is.na(rt_range)) {
    peak_data <- xcms::chromatogram(object = raw_data,
                                    mz = mz_range,
                                    rt = rt_range)
  } else{
    peak_data <- xcms::chromatogram(object = raw_data,
                                    mz = mz_range)
  }
  cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
  
  save(peak_data, file = file.path(output_path, "peak_data"))
  return(peak_data)
}


#' @title show_peak
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

show_peak = function(object,
                    peak.index = 1,
                    title.size = 15,
                    lab.size = 15,
                    axis.text.size = 15,
                    alpha = 0.5,
                    title = "",
                    interactive = TRUE){
  
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
    ggplot2::geom_line(data = data,
                       mapping = ggplot2::aes(colour = group)) +
    ggplot2::geom_area(mapping = ggplot2::aes(fill = group),
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














#' @title plot_tic
#' @description Extract TIC or BPC from in a RT range from mzXML format data.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param file name of mzXML data.
#' @param path work directory.
#' @param type tic or bpc.
#' @param threads thread number
#' @param legend legend or not.
#' @return A object for plot_chromatogram() function.
#' @export

plot_tic = function(file,
                    path = ".",
                    type = c("tic", "bpc"),
                    threads = 4,
                    legend = FALSE){
  type <- match.arg(type)
  sample_group <-
    rep("QC", length(file))
  
  pd <-
    data.frame(
      sample_name = sub(
        basename(file),
        pattern = ".mzXML",
        replacement = "",
        fixed = TRUE
      ),
      sample_group = sample_group,
      stringsAsFactors = FALSE
    )
  
  
  cat(crayon::green("Reading raw data, it will take a while...\n"))
  
  raw_data <- MSnbase::readMSData(
    files = file.path(path, file),
    pdata = new("NAnnotatedDataFrame", pd),
    mode = "onDisk",
    verbose = TRUE
  )
  
  
  cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
  
  #retention time correction
  #Alignment
  cat(crayon::green("Correcting rentention time...\n "))
  
  xdata <- try(xcms::adjustRtime(raw_data,
                                 param = xcms::ObiwarpParam(binSize = 0.5)),
               silent = TRUE)
  
  cat(crayon::red(clisymbols::symbol$tick, "OK\n"))
  
  if (class(xdata) == "try-error") {
    xdata <- raw_data
  }
  
  rm(list = "raw_data")
  
  if(tinyTools::get_os() == "windows"){
    bpparam =
      BiocParallel::SnowParam(workers = threads,
                              progressbar = TRUE)
  }else{
    bpparam = BiocParallel::MulticoreParam(workers = threads,
                                           progressbar = TRUE)
  }
  
  tic.plot <- xcms::chromatogram(
    object = xdata,
    aggregationFun = ifelse(type == "tic", "sum", "max"),
    BPPARAM = bpparam
  )
  
  plot <- plot_chromatogram(object = tic.plot,
                            title = "",
                            alpha = 1,
                            group.for.figure = "QC")
  plot +
    ggplot2::theme(legend.position = "none")
  
}







