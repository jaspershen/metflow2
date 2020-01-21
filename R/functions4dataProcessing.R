
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
#' @param group.for.figure What groups to show EIC.
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
                 interactive = FALSE,
                 group.for.figure = "QC") {
    options(warn = -1)
    info <- object@phenoData@data
    data <- object@.Data
    data <- data[1,which(info$sample_group %in% group.for.figure), drop = FALSE]
    info <- info[info$sample_group %in% group.for.figure,,drop = FALSE]
    
    if(nrow(info) == 0){
      return(NULL)
    }
    
    if(nrow(info) > 10){
      idx <- sort(sample(1:nrow(info), 10))
      info <- info[idx, , drop = FALSE]
      data <- data[,idx,drop = FALSE]
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
      ggplot2::geom_line(
        data = data,
        mapping = ggplot2::aes(colour = sample, group = sample)
      ) +
      ggsci::scale_fill_lancet() +
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
                         ppm = 15,
                         threads = 4,
                         is.table = "is.xlsx"
                         # rt.expand = 1
) {
  output_path <- path
  # dir.create(output_path)
  ##peak detection
  
  f.in <- list.files(path = path,
                     pattern = '\\.(mz[X]{0,1}ML|cdf)',
                     recursive = TRUE, 
                     full.names = TRUE)
  
  sample_group <-
    BiocGenerics::basename(f.in) %>% 
    stringr::str_replace("\\.(mz[X]{0,1}ML|cdf)", "")
  
  pd <-
    data.frame(
      # sample_name = sub(
      basename(f.in),
      #   pattern = ".mzXML",
      #   replacement = "",
      #   fixed = TRUE
      # ),
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
         file = file.path(output_path, "raw_data"),
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
  
  
  save(peak_data, file = file.path(output_path, "peak_data"))
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





##-------------------------------------------------------------------------
#' @title outputFeatureEIC
#' @description Output EICs of all peaks in peak table.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param object xdata class object from dataProcessing.
#' @param feature.index Which feature to output. Index.
#' @param polygon Add polygon or not.
#' @param save.plot Save plot or not.
#' @param path Work directory.
#' @param file.type pdf or png.
#' @param plot.name The name of plot.
#' @param width Width.
#' @param height Height.
#' @param interactive.plot Interactive or not.
#' @param return.plot Return plot or not.
#' @param message Show message or not.
#' @return ggplot object.
#' @export
#' @import xcms 
#' @import tidyverse
#' @import plotly


setGeneric(
  name = "outputFeatureEIC",
  def = function(object,
                 feature.index = 1,
                 polygon = TRUE,
                 save.plot = TRUE,
                 path = ".",
                 file.type = c("pdf", "png"),
                 plot.name = "peak_plot",
                 width = 7,
                 height = 7,
                 interactive.plot = FALSE,
                 return.plot = TRUE,
                 message = TRUE) {
    file.type <- match.arg(file.type)
    if (message) {
      cat(crayon::green("Extracting EICs..."))
    }
    
    feature_eic <- xcms::featureChromatograms(x = object,
                                              features = feature.index)
    if (message) {
      cat(crayon::red("OK\n"))
    }
    
    feature_eic_data <-
      object@.Data[1, ] %>%
      lapply(function(x) {
        if (nrow(x@chromPeaks) == 0) {
          data.frame(
            rt.med = NA,
            rt.min = NA,
            rt.max = NA,
            rt = NA,
            min.intensity = 0,
            max.intensity = NA,
            intensity = NA,
            stringsAsFactors = FALSE
          )
        } else{
          if (nrow(x@chromPeaks) > 1) {
            x@chromPeaks <-
              as_tibble(x@chromPeaks) %>%
              filter(maxo == max(maxo)) %>%
              as.matrix()
          }
          data.frame(
            rt.med = x@chromPeaks[, 4],
            rt.min = x@chromPeaks[, 5],
            rt.max = x@chromPeaks[, 6],
            rt = x@rtime,
            min.intensity = 0,
            max.intensity = x@chromPeaks[, "maxo"],
            intensity = x@intensity,
            stringsAsFactors = FALSE
          )
        }
      })
    
    feature_eic_data <-
      mapply(
        function(x, y, z) {
          data.frame(
            x,
            sample_group = y,
            sample_name = z,
            stringsAsFactors = FALSE
          ) %>%
            list()
        },
        x = feature_eic_data,
        y = object@phenoData@data$sample_group,
        z = object@phenoData@data$sample_name
      )
    
    
    feature_eic_data <-
      do.call(rbind, feature_eic_data)
    
    plot <-
      feature_eic_data %>%
      ggplot(aes(rt, intensity, group = sample_name)) +
      geom_line(aes(color = sample_group)) +
      ggsci::scale_color_lancet() +
      labs(x = "Retention time") +
      theme_bw()
    
    if (polygon) {
      plot2 <-
        plot +
        new_scale_color() +
        geom_rect(
          aes(
            xmin = rt.min,
            xmax = rt.max,
            ymin = min.intensity,
            ymax = max.intensity,
            color = sample_group
          ),
          fill = NA,
          linetype = 2
        ) +
        ggsci::scale_color_lancet(alpha = 0.3)
    }
    
    if (save.plot) {
      if (polygon) {
        if (file.type == "pdf") {
          ggplot2::ggsave(
            plot2,
            file = file.path(path, paste(plot.name, "pdf", sep = ".")),
            width = width,
            height = height
          )
        } else{
          ggplot2::ggsave(
            plot2,
            file = file.path(path, paste(plot.name, "png", sep = ".")),
            width = width,
            height = height
          )
        }
        
      } else{
        if (file.type == "pdf") {
          ggplot2::ggsave(
            plot,
            file = file.path(path, paste(plot.name, "pdf", sep = ".")),
            width = width,
            height = height
          )
        } else{
          ggplot2::ggsave(
            plot,
            file = file.path(path, paste(plot.name, "png", sep = ".")),
            width = width,
            height = height
          )
        }
      }
    }
    
    if (interactive.plot) {
      plot <-
        plotly::ggplotly(plot)
    }
    
    if (return.plot) {
      return(plot)
    }
    
  }
)



new_scale <- function(new_aes) {
  structure(ggplot2::standardise_aes_names(new_aes), class = "new_aes")
}

new_scale_color <- function() {
  new_scale("colour")
}


