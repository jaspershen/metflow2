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

chromatogramPlot = function(
  object,
  title.size = 15,
  lab.size = 15,
  axis.text.size = 15,
  alpha = 0.5,
  title = "",
  interactive = FALSE,
  group.for.figure = "QC"
){
  
  cat(crayon::yellow("chromatogramPlot() is deprecated, please use the plot_chromatogram().\n"))
  
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
                       mapping = ggplot2::aes(colour = sample, 
                                              group = sample)) +
    # ggsci::scale_fill_lancet() +
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

showPeak = function(object,
                    peak.index = 1,
                    title.size = 15,
                    lab.size = 15,
                    axis.text.size = 15,
                    alpha = 0.5,
                    title = "",
                    interactive = TRUE){
  
  cat(crayon::yellow("showPeak() is deprecated, please use the show_peak().\n"))
  
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

outputFeatureEIC = function(object,
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
                            message = TRUE){
  
  cat(crayon::yellow("outputFeatureEIC() is deprecated.\n"))
  
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
    object@.Data[1,] %>%
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
    # ggsci::scale_color_lancet() +
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
      ) 
    # ggsci::scale_color_lancet(alpha = 0.3)
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



#' @title new_scale
#' @description new_scale
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param new_aes new_aes
#' @return result

new_scale <- function(new_aes) {
  structure(ggplot2::standardise_aes_names(new_aes), class = "new_aes")
}

#' @title new_scale_color
#' @description new_scale_color
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @return result
new_scale_color <- function() {
  new_scale("colour")
}











#' @title getWorklist
#' @description Generate sample and worklist..
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param table.name Sample name.
#' @param instrument instrument.
#' @param each.mode.number each.mode.number.
#' @param randommethod See xcms.
#' @param samplenumber See xcms.
#' @param QCstep See xcms.
#' @param conditionQCnumber See xcms.
#' @param qc.index.from Number of threads.
#' @param dir dir.
#' @param method.path method.path.
#' @param ms1.method.pos ms1.method.pos.
#' @param ms1.method.neg ms1.method.neg.
#' @param ms2.method.pos ms2.method.pos.
#' @param ms2.method.neg ms2.method.neg.
#' @param path The working directory.
#' @return Peak table.
#' @export

getWorklist = function(table.name = "batch.xlsx",
                       instrument = c("Thermo", "Agilent", "AB"),
                       each.mode.number = 32,
                       randommethod = c("no", "position", "injection"),
                       samplenumber = NULL,
                       QCstep = 8,
                       conditionQCnumber = 8,
                       qc.index.from = 1,
                       dir = "D:\\Liang\\data\\PS4U\\HILIC\\batch3\\",
                       method.path = "D:\\Liang\\Method\\urine\\HILIC\\",
                       ms1.method.pos = "ZIC-HILIC_MS_pos",
                       ms1.method.neg = "ZIC-HILIC_MS_neg",
                       ms2.method.pos = c(
                         "ZIC-HILIC_MSMS_pos_NCE25",
                         "ZIC-HILIC_MSMS_pos_NCE25",
                         "ZIC-HILIC_MSMS_pos_NCE25",
                         "ZIC-HILIC_MSMS_pos_NCE25"
                       ),
                       ms2.method.neg = c(
                         "ZIC-HILIC_MSMS_neg_NCE25",
                         "ZIC-HILIC_MSMS_neg_NCE25",
                         "ZIC-HILIC_MSMS_neg_NCE25",
                         "ZIC-HILIC_MSMS_neg_NCE25"
                       ),
                       path = ".") {
  
  cat(crayon::yellow("getWorklist() is deprecated, please use the get_worklist().\n"))
  
  dir.create(path = path, showWarnings = FALSE)
  instrument <- match.arg(instrument)
  randommethod <- match.arg(randommethod)
  options(warn = -1)
  file <- dir(path)
  if(all(file != table.name)){
    stop("No ", table.name, " in ", path)
  }
  if (instrument == "Thermo") {
    batch <- readxl::read_excel(file.path(path, table.name))
    batch <- batch[, 1]
    if (randommethod == "no") {
      ###add position
      position <-
        unlist(lapply(c("B", "G", "R"), function(y) {
          paste(y,
                unlist(lapply(LETTERS[1:5], function(x) {
                  paste(x, 1:8, sep = "")
                })),
                sep = "")
        }))
      
      position <- position[-c(1:8)]
      
      condition.qc <-
        matrix(rep(
          c("Condition_QC", 'Condition_QC', "BA1"),
          conditionQCnumber
        ),
        ncol = 3,
        byrow = TRUE)
      
      condition.qc <-
        as.data.frame(condition.qc, stringsAsFactors = FALSE)
      
      colnames(condition.qc) <-
        c("File.Name", "Sample.ID", 'Position')
      
      blank.qc <-
        data.frame(
          "File.Name" = c("Blank", "QC"),
          "Sample.ID" = c("Blank", "QC"),
          "Position" = c("BA2", "BA1"),
          stringsAsFactors = FALSE
        )
      
      blank <-
        data.frame(
          "File.Name" = c("Blank", "Blank", "Blank"),
          "Sample.ID" = c("Blank", "Blank", "Blank"),
          "Position" = c("BA2", "BA2", "BA2"),
          stringsAsFactors = FALSE
        )
      
      dilution.qc <-
        data.frame(
          "File.Name" = c(
            "DL_QC1_1",
            "DL_QC2_1",
            "DL_QC4_1",
            "DL_QC8_1",
            "DL_QC1_2",
            "DL_QC2_2",
            "DL_QC4_2",
            "DL_QC8_2"
          ),
          "Sample.ID" = c(
            "DL_QC1_1",
            "DL_QC2_1",
            "DL_QC4_1",
            "DL_QC8_1",
            "DL_QC1_2",
            "DL_QC2_2",
            "DL_QC4_2",
            "DL_QC8_2"
          ),
          "Position" = c("BA3", "BA4", "BA5", "BA6", "BA3", "BA4", "BA5", "BA6"),
          stringsAsFactors = FALSE
        )
      
      ms2.qc <-
        data.frame(
          "File.Name" = c(
            "QC_MS2_NCE25_1",
            "QC_MS2_NCE25_2",
            "QC_MS2_NCE25_3",
            "QC_MS2_NCE25_4",
            "QC_MS2_NCE50_1",
            "QC_MS2_NCE50_2",
            "QC_MS2_NCE50_3",
            "QC_MS2_NCE50_4"
          ),
          "Sample.ID" = c(
            "QC_MS2_NCE25_1",
            "QC_MS2_NCE25_2",
            "QC_MS2_NCE25_3",
            "QC_MS2_NCE25_4",
            "QC_MS2_NCE50_1",
            "QC_MS2_NCE50_2",
            "QC_MS2_NCE50_3",
            "QC_MS2_NCE50_4"
          ),
          "Position" = c("BA1", "BA1", "BA1", "BA1", "BA1", "BA1", "BA1", "BA1"),
          stringsAsFactors = FALSE
        )
      
      batch <- cbind(batch, batch)
      batch <- data.frame(batch,
                          "Position" = rep(position, nrow(batch))[1:nrow(batch)],
                          stringsAsFactors = FALSE)
      colnames(batch) <-
        c("File.Name", "Sample.ID", "Position")
      temp.class <- sort(rep(1:nrow(batch), QCstep))
      temp.class <- temp.class[1:nrow(batch)]
      batch <-
        data.frame(temp.class, batch, stringsAsFactors = FALSE)
      batch <-
        plyr::dlply(.data = batch,
                    .variables = plyr::.(temp.class))
      
      batch <- lapply(batch, function(x) {
        rbind(blank.qc, x[, -1, drop = FALSE])
      })
      
      batch <- do.call(rbind, batch)
      batch <- rbind(batch, blank.qc)
      rownames(batch) <- NULL
      
      batch.pos <- data.frame(
        batch,
        Path = paste(dir, "POS", sep = ""),
        Instrument.Method = paste(method.path, ms1.method.pos, sep = "\\")
      )
      
      batch.neg <- data.frame(
        batch,
        Path = paste(dir, "NEG", sep = ""),
        Instrument.Method = paste(method.path, ms1.method.neg, sep = "\\")
      )
      
      batch.pos <-
        batch.pos[, c("File.Name",
                      "Sample.ID",
                      "Path",
                      "Instrument.Method",
                      "Position")]
      batch.neg <-
        batch.neg[, c("File.Name",
                      "Sample.ID",
                      "Path",
                      "Instrument.Method",
                      "Position")]
      
      condition.qc.pos <- data.frame(
        condition.qc,
        "Path" = paste(dir, 'POS', sep = ""),
        Instrument.Method = paste(method.path, ms1.method.pos, sep = "\\"),
        stringsAsFactors = FALSE
      )
      
      condition.qc.neg <- data.frame(
        condition.qc,
        "Path" = paste(dir, 'NEG', sep = ""),
        Instrument.Method = paste(method.path, ms1.method.neg, sep = "\\"),
        stringsAsFactors = FALSE
      )
      
      blank.pos <- data.frame(
        blank,
        "Path" = paste(dir, 'POS', sep = ""),
        Instrument.Method = paste(method.path, ms1.method.pos, sep = "\\"),
        stringsAsFactors = FALSE
      )
      
      blank.neg <- data.frame(
        blank,
        "Path" = paste(dir, 'NEG', sep = ""),
        Instrument.Method = paste(method.path, ms1.method.neg, sep = "\\"),
        stringsAsFactors = FALSE
      )
      
      
      dilution.qc.pos <- data.frame(
        dilution.qc,
        "Path" = paste(dir, 'POS', sep = ""),
        Instrument.Method = paste(method.path, ms1.method.pos, sep = "\\"),
        stringsAsFactors = FALSE
      )
      
      dilution.qc.neg <- data.frame(
        dilution.qc,
        "Path" = paste(dir, 'NEG', sep = ""),
        Instrument.Method = paste(method.path, ms1.method.neg, sep = "\\"),
        stringsAsFactors = FALSE
      )
      
      
      ms2.qc.pos <- data.frame(
        ms2.qc,
        "Path" = paste(dir, 'POS', sep = ""),
        Instrument.Method = paste(method.path, ms2.method.pos, sep = "\\"),
        stringsAsFactors = FALSE
      )
      
      ms2.qc.neg <- data.frame(
        ms2.qc,
        "Path" = paste(dir, 'NEG', sep = ""),
        Instrument.Method = paste(method.path, ms2.method.neg, sep = "\\"),
        stringsAsFactors = FALSE
      )
      
      condition.qc.pos <- condition.qc.pos[, colnames(batch.pos)]
      condition.qc.neg <- condition.qc.neg[, colnames(batch.neg)]
      
      blank.pos <- blank.pos[, colnames(batch.pos)]
      blank.neg <- blank.neg[, colnames(batch.neg)]
      
      dilution.qc.pos <- dilution.qc.pos[, colnames(batch.pos)]
      dilution.qc.neg <- dilution.qc.neg[, colnames(batch.neg)]
      
      ms2.qc.pos <- ms2.qc.pos[, colnames(batch.pos)]
      ms2.qc.neg <- ms2.qc.neg[, colnames(batch.neg)]
      
      batch.pos <-
        rbind(condition.qc.pos,
              dilution.qc.pos,
              ms2.qc.pos,
              batch.pos,
              blank.pos)
      batch.neg <-
        rbind(condition.qc.neg,
              dilution.qc.neg,
              ms2.qc.neg,
              batch.neg,
              blank.neg)
      
      ###rename
      batch.pos[which(batch.pos[, 1] == "Condition_QC"), 1] <-
        paste(batch.pos[which(batch.pos[, 1] == "Condition_QC"), 1],
              1:length(which(batch.pos[, 1] == "Condition_QC")), sep = "_")
      
      batch.pos[which(batch.pos[, 1] == "QC"), 1] <-
        paste("QC", qc.index.from:(sum(batch.pos[, 1] == "QC") + qc.index.from -
                                     1),
              sep = "_")
      
      batch.pos[which(batch.pos[, 1] == "Blank"), 1] <-
        paste(batch.pos[which(batch.pos[, 1] == "Blank"), 1],
              1:length(which(batch.pos[, 1] == "Blank")), sep = "_")
      
      batch.pos[, 2] <- batch.pos[, 1]
      
      
      batch.neg[which(batch.neg[, 1] == "Condition_QC"), 1] <-
        paste(batch.neg[which(batch.neg[, 1] == "Condition_QC"), 1],
              1:length(which(batch.neg[, 1] == "Condition_QC")), sep = "_")
      
      batch.neg[which(batch.neg[, 1] == "QC"), 1] <-
        paste("QC", qc.index.from:(sum(batch.neg[, 1] == "QC") + qc.index.from -
                                     1),
              sep = "_")
      
      batch.neg[which(batch.neg[, 1] == "Blank"), 1] <-
        paste(batch.neg[which(batch.neg[, 1] == "Blank"), 1],
              1:length(which(batch.neg[, 1] == "Blank")), sep = "_")
      
      batch.neg[, 2] <- batch.neg[, 1]
      
      write.csv(batch.pos, file.path(path, "worklist.pos.csv"), row.names = FALSE)
      write.csv(batch.neg, file.path(path, "worklist.neg.csv"), row.names = FALSE)
      
      
      ###pos and neg
      idx_blank <-
        which(stringr::str_detect(batch.pos$File.Name, "^Blank_[0-9]{1,3}"))
      idx_qc <-
        which(stringr::str_detect(batch.pos$File.Name, "^QC_[0-9]{1,3}"))
      
      # batch.pos$File.Name[idx_blank]
      # batch.pos$File.Name[idx_qc]
      
      ##for positive mode
      batch.pos.head <-
        batch.pos[1:(idx_blank[1] - 1),]
      
      batch.pos.tail <-
        batch.pos[(tail(idx_qc, 1) + 1):nrow(batch.pos),]
      
      batch.pos.middle <-
        batch.pos[idx_blank[1]:tail(idx_qc, 1),]
      
      temp_class <- rep(1:(QCstep + 2))
      
      idx1 <- grep("^Blank", batch.pos.middle$File.Name)
      idx1 <- idx1[-length(idx1)]
      
      
      idx1 <-
        c(idx1[seq(1, length(idx1), by = round(each.mode.number / QCstep))], tail(idx1, 1)) %>%
        unique
      
      idx2 <- c(idx1[-1] - 1, nrow(batch.pos.middle))
      
      
      batch.pos.middle <-
        mapply(function(x, y) {
          list(batch.pos.middle[x:y,])
        },
        x = idx1,
        y = idx2)
      
      
      
      
      batch.pos.middle <-
        lapply(batch.pos.middle, function(x) {
          temp.idx <- grep("Blank", x$File.Name)
          temp.idx <- temp.idx[-1]
          x <- x[-temp.idx, , drop = FALSE]
          x
        })
      
      batch.pos.middle <-
        lapply(batch.pos.middle, function(x) {
          if (length(grep("QC", tail(x$File.Name))) == 0) {
            x <- rbind(x, x[grep("QC", x$File.Name)[1], ])
            x
          } else{
            x
          }
        })
      
      
      #rename
      
      
      
      batch.neg.head <-
        batch.neg[1:(idx_blank[1] - 1),]
      
      batch.neg.tail <-
        batch.neg[(tail(idx_qc, 1) + 1):nrow(batch.neg),]
      
      batch.neg.middle <-
        batch.neg[idx_blank[1]:tail(idx_qc, 1),]
      
      temp_class <- rep(1:(QCstep + 2))
      
      idx1 <- grep("^Blank", batch.neg.middle$File.Name)
      idx1 <- idx1[-length(idx1)]
      
      idx1 <-
        c(idx1[seq(1, length(idx1), by = round(each.mode.number / QCstep))], tail(idx1, 1)) %>%
        unique
      
      idx2 <- c(idx1[-1] - 1, nrow(batch.neg.middle))
      
      
      batch.neg.middle <-
        mapply(function(x, y) {
          list(batch.neg.middle[x:y,])
        },
        x = idx1,
        y = idx2)
      
      batch.neg.middle <-
        lapply(batch.neg.middle, function(x) {
          temp.idx <- grep("Blank", x$File.Name)
          temp.idx <- temp.idx[-1]
          x <- x[-temp.idx, , drop = FALSE]
          x
        })
      
      batch.neg.middle <-
        lapply(batch.neg.middle, function(x) {
          if (length(grep("QC", tail(x$File.Name))) == 0) {
            x <- rbind(x, x[grep("QC", x$File.Name)[1], ])
            x
          } else{
            x
          }
        })
      
      batch.neg.head <-
        batch.neg.head[-grep("Condition_QC", batch.neg.head$File.Name),]
      
      
      batch.pos.middle[[1]] <-
        rbind(batch.pos.head, batch.pos.middle[[1]])
      
      batch.neg.middle[[1]] <-
        rbind(batch.neg.head, batch.neg.middle[[1]])
      
      batch.neg.middle[[length(batch.neg.middle)]] <-
        rbind(batch.neg.middle[[length(batch.neg.middle)]],
              batch.neg.tail)
      
      batch <- mapply(function(x, y) {
        list(rbind(x, y))
      },
      x = batch.pos.middle,
      y = batch.neg.middle)
      
      batch <- do.call(rbind, batch)
      
      
      # batch$File.Name
      
      mode <- rep(NA, nrow(batch))
      mode[grep("POS", batch$Path)] <- "POS"
      mode[grep("NEG", batch$Path)] <- "NEG"
      
      blank_pos <-
        which(
          mode == "POS" &
            stringr::str_detect(string = batch$File.Name, pattern = "^Blank_")
        )
      
      blank_neg <-
        which(
          mode == "NEG" &
            stringr::str_detect(string = batch$File.Name, pattern = "^Blank_")
        )
      
      
      QC_pos <-
        which(
          mode == "POS" &
            stringr::str_detect(string = batch$File.Name, pattern = "^QC_[0-9]{1,3}")
        )
      
      QC_neg <-
        which(
          mode == "NEG" &
            stringr::str_detect(string = batch$File.Name, pattern = "^QC_[0-9]{1,3}")
        )
      
      batch$File.Name[blank_pos] <-
        batch$Sample.ID[blank_pos] <-
        paste("Blank", 1:length(blank_pos), sep = "_")
      
      batch$File.Name[blank_neg] <-
        batch$Sample.ID[blank_neg] <-
        paste("Blank", 1:length(blank_neg), sep = "_")
      
      batch$File.Name[QC_pos] <-
        batch$Sample.ID[QC_pos] <-
        paste("QC", qc.index.from:(length(QC_pos) + qc.index.from - 1), sep = "_")
      
      batch$File.Name[QC_neg] <-
        batch$Sample.ID[QC_neg] <-
        paste("QC", qc.index.from:(length(QC_pos) + qc.index.from - 1), sep = "_")
      
      
      
      write.csv(batch, file.path(path, "worklist.csv"), row.names = FALSE)
      
    }
  }
}

















#' @title Impute MV in data.
#' @description Impute MV in data.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object A metflowClass object.
#' @param method Imputation method. It
#' contains "knn", "rf" (missForest), "mean", "median", "zero", "minium",
#' "bpca" (BPCA), "svd" (SVD) and "ppca" (PPCA). Default is "knn".
#' The detial of
#' this method can be find in detail and reference paperes.
#' @param k See ?impute.knn
#' @param rowmax See ?impute.knn
#' @param colmax See ?impute.knn
#' @param maxp See ?impute.knn
#' @param rng.seed See ?impute.knn
#' @param maxiter See ?missForest
#' @param ntree See ?missForest
#' @param decreasing See ?missForest
#' @param nPcs See ?bpca
#' @param maxSteps See ?bpca
#' @param threshold See ?bpca
#' @param ... Other arguments.
#' @return A new metflowClass object.
#' @export

imputeMV = function(object,
                    method = c("knn",
                               "rf",
                               "mean",
                               "median",
                               "zero",
                               "minimum",
                               "bpca",
                               "svdImpute",
                               "ppca"),
                    k = 10,
                    rowmax = 0.5,
                    colmax = 0.8,
                    maxp = 1500,
                    rng.seed = 362436069,
                    # missForest parameters
                    maxiter = 10,
                    ntree = 100,
                    decreasing = FALSE,
                    #BPCA PPCA, and SVD parameters
                    nPcs = 2,
                    maxSteps = 100,
                    threshold = 1e-04,
                    ...){
  cat(crayon::yellow("`imputeMV()` is deprecated, please use `impute_mv()`"))
  
  options(warn = -1)
  if (class(object) != "metflowClass") {
    stop("Only the metflowClass is supported!\n")
  }
  
  #### MV imputation
  ms1_data <- object@ms1.data
  if (length(ms1_data) > 1) {
    stop("Please align batch first.\n")
  }
  ms1_data <- ms1_data[[1]]
  qc_data <- getData(object = object, slot = "QC")
  subject_data <-
    getData(object = object, slot = "Subject")
  subject_qc_data <- cbind(qc_data, subject_data)
  subject_qc_data <- tibble::as_tibble(subject_qc_data)
  
  if (sum(is.na(subject_qc_data)) == 0) {
    cat("No missing values.\n")
    invisible(object)
  }
  
  subject_qc_data <- sxtMVimputation(
    data = subject_qc_data,
    method = method,
    # knn parameters
    k = k,
    rowmax = rowmax,
    colmax = colmax,
    maxp = maxp,
    rng.seed = rng.seed,
    # missForest parameters
    maxiter = maxiter,
    ntree = ntree,
    decreasing = decreasing,
    #BPCA PPCA, and SVD parameters
    nPcs = nPcs,
    maxSteps = maxSteps,
    threshold = threshold,
    ...
  )
  
  object@process.info$imputeMV <- list()
  object@process.info$imputeMV$method = method
  object@process.info$imputeMV$k = k
  object@process.info$imputeMV$rowmax = rowmax
  object@process.info$imputeMV$colmax = colmax
  object@process.info$imputeMV$maxp = maxp
  object@process.info$imputeMV$rng.seed = rng.seed
  object@process.info$imputeMV$maxiter = maxiter
  object@process.info$imputeMV$ntree = ntree
  object@process.info$imputeMV$decreasing = decreasing
  object@process.info$imputeMV$nPcs = nPcs
  object@process.info$imputeMV$maxSteps = maxSteps
  object@process.info$imputeMV$threshold = threshold
  object@process.info$imputeMV$maxSteps = maxSteps
  
  sample_info <- object@sample.info
  subject_qc_name <-
    dplyr::filter(.data = sample_info, class %in% c("Subject", "QC")) %>%
    dplyr::pull(., sample.name)
  
  subject_qc_data <-
    subject_qc_data[, match(subject_qc_name, colnames(subject_qc_data))]
  
  ms1_data[, match(subject_qc_name, colnames(ms1_data))] <-
    subject_qc_data
  ms1_data <- list(ms1_data)
  object@ms1.data <- ms1_data
  invisible(object)
}











#' @title creatMetflowObject
#' @description Creat metflowClass object.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param ms1.data MS1 peak table name.
#' @param sample.information Sample information name.
#' @param path Work directory.
#' @return A metflowClass object.
#' @export

creatMetflowObject = function(ms1.data,
                              sample.information,
                              path = ".") {
  
  cat(crayon::yellow("`creatMetflowObject()` is deprecated, please use `create_metflow_object()`"))
  
  check.result <- check_data(data = ms1.data,
                             sample.info = sample.information,
                             path = path)
  
  if (any(check.result$`Check result` != "Valid")) {
    stop("Error is in you data. Please check it.\n")
  }
  
  ms1.data.name <- ms1.data
  sample.information.name <- sample.information
  
  ms1.data <- pbapply::pblapply(ms1.data, function(x) {
    readr::read_csv(file = file.path(path, x),
                    col_types = readr::cols())
  })
  
  sample.info <-
    readr::read_csv(file.path(path, sample.information),
                    col_types = readr::cols())
  
  object <- new(
    Class = "metflowClass",
    ms1.data = ms1.data,
    sample.info = sample.info,
    version = "0.1.0"
  )
  invisible(object)
}








#' @title getData
#' @description Get data from metflowClass object.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object A metflowClass object.
#' @param slot Class of data.
#' @param silence.deprecated Silence deprecated information or not.
#' @return A data frame.
#' @export

getData = function(
  object,
  slot = c("Subject", "QC", "QC.DL", "Blank", "Tags"),
  silence.deprecated = TRUE
){
  if (!silence.deprecated) {
    cat(crayon::yellow("`getData()` is deprecated, please use `get_data()`"))
  }
  
  # requireNamespace("magrittr")
  if (class(object) != "metflowClass") {
    stop("Only the metflowClass is supported!\n")
  }
  
  if (length(object@ms1.data) > 1) {
    stop("Plase align batch first.\n")
  }
  
  # slot <- stringr::str_to_title(slot)
  slot <- match.arg(slot)
  
  if (slot == "Tags") {
    result <- object@ms1.data[[1]] %>%
      dplyr::select(., -one_of(object@sample.info$sample.name))
    return(result)
  }
  
  result <-
    try(dplyr::filter(.data = object@sample.info, class == slot)$sample.name %>%
          dplyr::select(.data = object@ms1.data[[1]], .))
  
  if (ncol(result) == 0) {
    return(NULL)
  }
  return(result)
}




#' @title getMVplot4sample
#' @description get MV plot of subject samples.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object A metflowClass object.
#' @return A ggplot2 object.
#' @export

getMVplot4sample = function(object){
  if (!silence.deprecated) {
    cat(
      crayon::yellow(
        "`getMVplot4sample()` is deprecated, please use `get_mv_plot_samples()`"
      )
    )
  }
  
  
  if (class(object) != "metflowClass") {
    stop("Only the metflowClass is supported!\n")
  }
  plot <- try(object@process.info$filterSample$plot)
  if (class(plot)[1] == "try-error") {
    return(NULL)
  } else{
    plot
  }
}
