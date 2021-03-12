#' @title output_eic
#' @description Output EIC for some peaks in some samples.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param path Work directory.
#' @param query_sample_name sample names
#' @param query_peak_name Feature names
#' @param polarity The polarity of data, "positive"or "negative".
#' @param threads Number of threads.
#' @return EICs.
#' @export


# tinyTools::setwd_project()
# setwd("example/POS/")
# 
# peak_table = readr::read_csv("Result/Peak_table_for_cleaning.csv")
# 
# query_sample_name = colnames(peak_table)[c(4,5)]
# query_peak_name = peak_table$name[sample(1:10000, 20)]
# 
# output_eic(
#   path = ".",
#   query_sample_name = query_sample_name,
#   query_peak_name = query_peak_name,
#   polarity = "positive",
#   threads = 5
# )

output_eic = function(path = ".",
                      query_sample_name,
                      query_peak_name,
                      polarity = c("positive", "negative"),
                      threads = 6) {
  # browser()
  if(missing(query_sample_name)){
    stop("Please provide sample names.\n")
  }
  
  if(missing(query_peak_name)){
    stop("Please provide feature names.\n")
  }
  
  polarity <- match.arg(polarity)
  output_path <- file.path(path, "Result")
  dir.create(output_path, showWarnings = FALSE)
  intermediate_data_path <-
    file.path(output_path, "intermediate_data")
  dir.create(intermediate_data_path, showWarnings = FALSE)
  
  f.in <- list.files(path = path,
                     pattern = '\\.(mz[X]{0,1}ML|cdf)',
                     recursive = TRUE)
  
  sample_group <-
    unlist(lapply(stringr::str_split(string = f.in, pattern = "/"), function(x) {
      x[1]
    }))
  
  sample_group[grep("\\.(mz[X]{0,1}ML|cdf)", sample_group)] <-
    "Group0"

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
  
  names(group_colors) <- unique(sample_group)
  
  if (all(dir(intermediate_data_path) != "xdata3")) {
    stop("Please run processData() or process_data() first.\n")
  } else{
    load(file.path(intermediate_data_path, "xdata3"))
  }

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
  peak_name <-
    paste(peak_name, ifelse(polarity == "positive", "POS", "NEG"), sep = "_")
  
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
  
  ##-------------------------------------------------------------------------------------
  ##output EIC of all peaks
  query_sample_name2 = 
    intersect(query_sample_name, colnames(peak_table))
  
  if(length(query_sample_name2) == 0){
    stop("All the samples you provided are not in this project.\n")
  }
  
  query_peak_name2 = 
    intersect(query_peak_name, peak_table$peak.name)
  
  if(length(query_peak_name2) == 0){
    stop("All the features you provided are not in this project.\n")
  }
  
  is_table <-
    peak_table[,c("peak.name", "mzmed", "rtmed")] %>% 
    dplyr::filter(peak.name %in% query_peak_name2) %>% 
    dplyr::rename(name = peak.name, mz = mzmed, rt = rtmed)
  
  # if (output.peak.eic) {
  options(warn = -1)  
  cat(crayon::green("Outputting peak EICs...\n"))
    name = paste("feature_EIC", as.character(Sys.time()), sep = "_") %>% 
      stringr::str_replace_all("-", "_") %>% 
      stringr::str_replace_all(":", "_") %>% 
      stringr::str_replace_all(" ", "_")
    feature_EIC_path <- file.path(output_path, name)
    dir.create(feature_EIC_path, showWarnings = FALSE)
    
    temp_fun <- function(idx = 1,
                         feature_eic_data,
                         path = ".",
                         peak.name,
                         metabolite.name) {
      
      peak.name <- peak.name[idx]
      plot.name <- paste("", peak.name, sep = "")
      metabolite.name <- metabolite.name[idx]
      feature_eic_data <-
        feature_eic_data[[idx]]
      
      rt_range <-
        c(min(feature_eic_data$rt, na.rm = TRUE),
          max(feature_eic_data$rt, na.rm = TRUE))
      
      if (nrow(feature_eic_data) != 0) {
        plot <-
          ggplot2::ggplot(feature_eic_data,
                          ggplot2::aes(rt, intensity, group = sample_name)) +
          ggplot2::geom_line(ggplot2::aes(color = sample_name)) +
          ggplot2::geom_point(ggplot2::aes(color = sample_name)) +
          # ggsci::scale_color_lancet() +
          ggplot2::labs(
            x = "Retention time",
            title = paste(
              "RT range:",
              rt_range[1],
              rt_range[2],
              metabolite.name,
              sep = "_"
            )
          ) +
          ggplot2::theme_bw()
        
        ggplot2::ggsave(
          plot,
          file = file.path(path, paste(plot.name, "pdf", sep = ".")),
          width = 10,
          height = 6
        )
      }
    }
    
    metabolite_name <-
      is_table$name
    
    if (tinyTools::get_os() == "windows") {
      bpparam =
        BiocParallel::SnowParam(workers = threads,
                                progressbar = TRUE)
    } else{
      bpparam = BiocParallel::MulticoreParam(workers = threads,
                                             progressbar = TRUE)
    }
    
    index2 = match(is_table$name, peak_table$peak.name)
    
    feature_eic <-
      xcms::featureChromatograms(
        x = xdata3,
        features = index2,
        expandRt = 0,
        BPPARAM = bpparam
      )
    
    feature_eic_data <- feature_eic@.Data %>%
      pbapply::pbapply(1, function(y) {
        y <- lapply(y, function(x) {
          if (class(x) == "XChromatogram") {
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
                  tibble::as_tibble(x@chromPeaks) %>%
                  dplyr::filter(maxo == max(maxo)) %>%
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
          } else{
            
          }
        })
        y <-
          mapply(
            function(y, sample.group, sample.name) {
              data.frame(
                y,
                sample_group = sample.group,
                sample_name = sample.name,
                stringsAsFactors = FALSE
              ) %>%
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
      lapply(function(x) {
        x <-
          x %>%
          dplyr::filter(sample_name %in% query_sample_name2)
        
        if (length(unique(x$sample_name)) > 8) {
          idx <-
            which(x$sample_name %in% sort(sample(unique(
              x$sample_name
            ), 18))) %>%
            sort()
          x <- x[idx, , drop = FALSE]
        }
        x
      })
    
    if (tinyTools::get_os() == "windows") {
      bpparam =
        BiocParallel::SnowParam(workers = threads,
                                progressbar = TRUE)
    } else{
      bpparam = BiocParallel::MulticoreParam(workers = threads,
                                             progressbar = TRUE)
    }
    
    BiocParallel::bplapply(
      1:length(index2),
      FUN = temp_fun,
      BPPARAM = bpparam,
      feature_eic_data = feature_eic_data,
      path = feature_EIC_path,
      peak.name = peak_name[index2],
      metabolite.name = metabolite_name
    )
  # }
  rm(list = "xdata3")
  cat(crayon::red(clisymbols::symbol$tick, "Done\n"))
  }