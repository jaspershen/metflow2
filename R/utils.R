#' @title calFC
#' @description Calculate fold change for metflowClass.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object metflowClass object
#' @param control.group Control group name.
#' @param case.group Case group name.
#' @param type type
#' @return Fold change.
#' @export

calFC = function(object,
                 control.group,
                 case.group,
                 type = c("median", "mean")){
  type <- match.arg(type)
  if (missing(control.group) | missing(case.group)) {
    stop("Please set control.group or case.group.\n")
  }
  
  if (!all(c(control.group, case.group) %in% object@sample.info$group)) {
    stop("Please make sure that the control and case group are in your sample information.\n")
  }
  
  ms1_data <- object@ms1.data
  if (length(ms1_data) > 1) {
    stop("Please align batches first.\n")
  }
  
  ms1_data <- ms1_data[[1]]
  sample_info <- object@sample.info
  
  control_data <-
    sample_info %>%
    filter(., group == control.group) %>%
    dplyr::pull(., sample.name) %>%
    match(., colnames(ms1_data)) %>%
    ms1_data[, .]
  
  case_data <-
    sample_info %>%
    filter(., group == case.group) %>%
    dplyr::pull(., sample.name) %>%
    match(., colnames(ms1_data)) %>%
    ms1_data[, .]
  
  if (sum(is.na(control_data)) != 0 |
      sum(is.na(case_data)) != 0) {
    stop("Please impute MV first.\n")
  }
  
  fc <- apply(case_data, 1, function(x) {
    switch(type,
           mean = mean(x),
           median = median(x))
  }) / apply(control_data, 1, function(x) {
    switch(type,
           mean = mean(x),
           median = median(x))
  })
  
  name <- object %>%
    getData(., slot = "Tags") %>%
    dplyr::pull(., name)
  getData(object = object, slot = "Tags")
  fc <- data.frame(index = 1:length(fc),
                   name,
                   fc,
                   stringsAsFactors = FALSE)
  rownames(fc) <- NULL
  invisible(fc)
}



#' @title SXTMTmatch
#' @description Match two data according to mz and RT.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param data1 First data for matching, first column must be mz
#' and seconod column must be rt.
#' @param data2 Second data for matching, first column must be mz
#' and seconod column must be rt.
#' @param mz.tolerance mz tolerance for ms1 and ms2 data matching.
#' @param rt.tolerance RT tolerance for ms1 and ms2 data matching.
#' @return Return a result which give the matching result of data1 and database.


SXTMTmatch = function(data1,
                      data2,
                      mz.tolerance = 25,
                      rt.tolerance = 180) {
  if (nrow(data1) == 0 | nrow(data2) == 0) {
    result <- NULL
    return(result)
  }
  mz1 <- as.numeric(data1[, 1])
  rt1 <- as.numeric(data1[, 2])
  
  mz2 <- as.numeric(data2[, 1])
  rt2 <- as.numeric(data2[, 2])
  
  result <- NULL
  cat("finished: %")
  cat("\n")
  for (i in seq_along(mz1)) {
    mz.error <- abs(mz1[i] - mz2) * 10 ^ 6 / mz1[i]
    rt.error <- abs(rt1[i] - rt2)
    j <- which(mz.error <= mz.tolerance & rt.error <= rt.tolerance)
    if (length(j) != 0) {
      result1 <-
        cbind(i, j, mz1[i], mz2[j], mz.error[j], rt1[i], rt2[j], rt.error[j])
      result <- rbind(result, result1)
    }
    
    count <- floor((length(mz1)) * c(seq(0, 1, 0.01)))
    if (any(i == count)) {
      cat(ceiling (i * 100 / length(mz1)))
      cat(" ")
    }
    
  }
  cat("\n")
  if (is.null(result)) {
    cat("There are not any peak be matched\n,
        please change the mz or rt tolerance and try again")
    cat("\n")
  }
  else {
    number1 <- length(unique(result[, 1]))
    number2 <- length(unique(result[, 2]))
    cat(
      paste(
        "There are",
        number1,
        "peaks in data1, and",
        number2,
        "peaks in data2 are matched"
      )
    )
    cat("\n")
    colnames(result) <-
      c("Index1",
        "Index2",
        "mz1",
        "mz2",
        "mz error",
        "rt1",
        "rt2",
        "rt error")
    return(result)
  }
}





#' @title calP
#' @description Calculate p value and AUC for each feature.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object object
#' @param control.group control.group.
#' @param case.group case.group.
#' @param test.method test.method.
#' @param adjust.method adjust.method.
#' @param log.scale log.scale.
#' @param ... other arguments.
#' @return p values.
#' @export

calP = function(object,
                control.group,
                case.group,
                test.method = c("t", "wilcox"),
                adjust.method = c("holm",
                                  "hochberg",
                                  "hommel",
                                  "bonferroni",
                                  "BH",
                                  "BY",
                                  "fdr",
                                  "none"),
                log.scale = FALSE,
                ...) {
  test.method <- match.arg(test.method)
  adjust.method <- match.arg(adjust.method)
  if(missing(control.group) | missing(case.group)){
    stop("Please set control.group or case.group.\n")
  }
  
  if(!all(c(control.group, case.group) %in% object@sample.info$group)){
    stop("Please make sure that the control and case group are in your sample information.\n")
  }
  
  ms1_data <- object@ms1.data
  if(length(ms1_data) > 1){
    stop("Please align batches first.\n")
  }
  
  ms1_data <- ms1_data[[1]]
  sample_info <- object@sample.info
  
  control_data <-
    sample_info %>%
    filter(., group == control.group) %>%
    dplyr::pull(., sample.name) %>%
    match(., colnames(ms1_data)) %>%
    ms1_data[,.]
  
  case_data <-
    sample_info %>%
    filter(., group == case.group) %>%
    dplyr::pull(., sample.name) %>%
    match(., colnames(ms1_data)) %>%
    ms1_data[,.]
  
  if(sum(is.na(control_data)) != 0 | sum(is.na(case_data)) != 0){
    stop("Please impute MV first.\n")
  }
  
  if(log.scale == TRUE){
    control_data <- log(control_data + 1, 10)
    case_data <- log(case_data + 1, 10)
  }
  
  control_data <- split(control_data, seq(nrow(control_data)))
  case_data <- split(case_data, seq(nrow(case_data)))
  pbapply::pboptions(type = "timer")
  
  p_value <- pbapply::pbmapply(function(x, y){
    temp_p <- try(
      switch(test.method,
             t = t.test(as.numeric(x), as.numeric(y), ...),
             wilcox = wilcox.test(as.numeric(x), as.numeric(y), ...))
    )
    temp_p$p.value
  },
  x = control_data,
  y = case_data)
  
  p_value <-
    object %>%
    getData(., "Tags") %>%
    dplyr::pull(., name) %>%
    data.frame(index = 1:length(p_value), ., p_value)
  colnames(p_value)[2] <- "name"
  
  if (adjust.method != "none") {
    p_value <- p_value %>%
      dplyr::mutate(., p_value_adjust = p.adjust(p_value,
                                                 method = adjust.method))
  }
  return(p_value)
}



#' @title VolcanoPlot
#' @description Draw volcano plot.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param p.value p.value.
#' @param fc fc.
#' @param p.cutoff p.cutoff.
#' @param fc.cutoff fc.cutoff.
#' @param xlab xlab.
#' @param ylab ylab.
#' @return volcano plot.
#' @export

VolcanoPlot = function(p.value,
                       fc,
                       p.cutoff = 0.05,
                       fc.cutoff = 1.3,
                       xlab = "Log2 (Fold change)",
                       ylab = "-Log10 (P-value)"){
  p.value <- as.numeric(p.value)
  fc <- as.numeric(fc)
  
  marker <- rep("No", length(p.value))
  marker[p.value < p.cutoff &
           fc > fc.cutoff] <- "Increase"
  marker[p.value < p.cutoff &
           fc < 1 / fc.cutoff] <- "Decrease"
  
  data <-
    data.frame(p.value, fc, marker, stringsAsFactors = FALSE)
  
  # requireNamespace("ggplot2")
  plot <- ggplot(data, aes(log(fc, 2), -log(p.value, 10),
                           colour = marker)) +
    geom_point() +
    scale_colour_manual(
      guide = guide_legend(title = NULL),
      values = c(
        'No' = "grey",
        "Increase" = "#ED00007F",
        'Decrease' = "#00468B7F"
      )
    ) +
    theme_bw() +
    geom_hline(yintercept = -log(p.cutoff, 10),
               linetype = 2) +
    geom_vline(xintercept = c(log(fc.cutoff, 2), -log(fc.cutoff, 2)), linetype = 2) +
    labs(x = xlab, y = ylab) +
    theme(
      axis.title = element_text(size = 15),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 15),
      legend.text = element_text(size = 12),
      legend.position = c(1, 1),
      legend.justification = c(1, 1),
      legend.background = element_rect(fill = NULL, colour = "black")
    )
  plot
  
}




#' @title sxtScale
#' @description sxtScale
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param df df
#' @param method method
#' @return result

sxtScale = function(df,
                    method = c("no", "auto", "pareto", "center")) {
  
  method <- match.arg(method)
  if(method == "no"){
    return(as.data.frame(df))
  }
  
  if (method == "auto") {
    df <- apply(df, 2, function(x) {
      (x - mean(x)) / sd(x)
    })
    return(as.data.frame(df))
  }
  
  if (method == "pareto") {
    df <- apply(df, 2, function(x) {
      (x - mean(x, na.rm = TRUE)) / sqrt(sd(x, na.rm = TRUE))
    })
    return(as.data.frame(df))
  }
  
  if (method == "center") {
    df <- apply(df, 2, function(x) {
      x - mean(x, na.rm = TRUE)
    })
    return(as.data.frame(df))
  }
}



#' @title showPeakPlot
#' @description Show peak plot.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object1 metflowClass object 1.
#' @param object2 metflowClass object 2.
#' @param peak.index Peak index
#' @return Peak plot (ggplot2 object).
#' @export

showPeakPlot = function(object1, object2, peak.index = 1){
  object <- list(object1, object2)
  if (lapply(object, function(x) {
    class(x) != "metflowClass"
  }) %>%
  unlist() %>%
  any()) {
    stop("Object must be metflowClass object.\n")
  }
  
  if(
    lapply(object, function(x){
      nrow(x@ms1.data[[1]])
    }) %>% 
    unlist() %>% 
    unique() %>% 
    length() > 1  
  ){
    stop("Objects must have same peak number.\n")
  } 
  
  peak_plots <- 
    lapply(object, function(x){
      ms1_data <- x@ms1.data[[1]]
      sample_info <- x@sample.info
      sample_info <- 
        sample_info %>% 
        filter(class == "Subject" | class == "QC")
      ms1_data <- 
        ms1_data %>% 
        select(one_of(sample_info$sample.name))
      if(peak.index > nrow(ms1_data)){
        peak.index <- nrow(ms1_data)
        cat("peak.index is bigger than the nrow of data.\n")
      }
      
      peak_data <- ms1_data[peak.index,]
      peak_data <- 
        peak_data[,match(sample_info$sample.name, colnames(peak_data))]  
      
      peak_data <-
        data.frame(
          'injection.order' = sample_info$injection.order,
          'class' = sample_info$class,
          'batch' = sample_info$batch,
          'int' = unlist(peak_data),
          stringsAsFactors = FALSE
        ) %>% 
        arrange(., injection.order)
      
      v_lines <- 
        which(!duplicated(sample_info$batch))[-1]
      peak_plot <-
        ggplot(peak_data, aes(injection.order, int)) +
        geom_vline(xintercept = v_lines, linetype = 2) +
        geom_point(aes(colour = class, shape = factor(batch))) +
        theme_bw() +
        labs(x = "Injection order", 
             y = "Intensity", 
             title = x@ms1.data[[1]]$name[peak.index]) +
        scale_colour_manual(values = c("#E18727FF", "#0072B5FF")) +
        guides(colour = guide_legend(title = "Class",
                                     override.aes = list(size = 3)),
               shape = guide_legend(title = "Batch", 
                                    override.aes = list(size = 3))) +
        theme(axis.title = element_text(size = 15),
              axis.text = element_text(size = 12),
              legend.title = element_text(size = 15),
              legend.text = element_text(size = 12))
      
      peak_plot
    })
  
  plot <- 
    peak_plots[[1]] + peak_plots[[2]] + patchwork::plot_layout(ncol = 1)
  
  plot
}

