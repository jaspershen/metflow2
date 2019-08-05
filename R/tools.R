#' @title calFC
#' @description Calculate fold change for metflowClass.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object metflowClass object
#' @param control.group Control group name.
#' @param case.group Case group name.
#' @return Fold change.
#' @export

calFC <- function(object,
                  control.group,
                  case.group,
                  type = c("median", "mean")) {
  # requireNamespace("tidyverse")
  type <- match.arg(type)
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
  
  fc <- apply(case_data, 1, function(x){
    switch(type, 
           mean = mean(x),
           median = median(x))
  })/apply(control_data, 1, function(x){
    switch(type, 
           mean = mean(x),
           median = median(x))
  })
  
  name <- object %>% 
    getData(., slot = "Tags") %>% 
    dplyr::pull(., name)
    getData(object = object, slot = "Tags")
  fc <- data.frame(index = 1:length(fc), name, fc,
                   stringsAsFactors = FALSE)
  rownames(fc) <- NULL
  invisible(fc)
}


#' @title HeatMap
#' @description Heat map
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param log.scale log transformation or not.
#' @param color Color list for sample group.
#' @param variable "all" or "marker" for heatmap.
#' @param Group group for heatmap.
#' @param scale.method scale method.
#' @param show_rownames Default is FALSE.
#' @param show_colnames Default is FALSE.
#' @param path Work directory.
#' @param width Plot width
#' @param height Plot height.
#' @param border_color Default is NA.
#' @param fontsize_row Default is 10.
#' @param cluster_rows Default is TRUE.
#' @param cluster_cols Default is TURE.
#' @param clustering_method Default is"ward.D",
#' @param ... other parameters for pheatmap.
#' @return A heatmap plot.
#' @seealso \code{\link[pheatmap]{pheatmap}}

HeatMap <- function(MetFlowData,
                    log.scale = FALSE,
                    color = c("palegreen",
                              "firebrick1",
                              "royalblue",
                              "yellow",
                              "black",
                              "cyan",
                              "gray48"),
                    variable = "all",
                    Group = c("control", "case"),
                    scale.method = "auto",
                    show_rownames = FALSE,
                    show_colnames = FALSE,
                    path = ".",
                    width = 7,
                    height = 7,
                    border_color = NA,
                    fontsize_row = 10,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE,
                    clustering_method = "ward.D",
                    ...) {
  if (path != ".") {
    dir.create(path)
  }

  subject <- MetFlowData@subject
  tags <- MetFlowData@tags
  subject.info <- MetFlowData@subject.info
  group <- subject.info[, "group"]

  idx <- which(group %in% Group)
  subject.info <- subject.info[idx, ]
  subject <- subject[, idx]
  group <- subject.info[, "group"]
  ## data organization
  if (variable == "all") {
    data <- t(subject)
  } else{
    if (all(colnames(tags) != "is.marker")) {
      stop("Please select marker first.")
    }
    is.marker <- tags[, "is.marker"]
    var.index <- which(is.marker == "yes")
    data <- t(subject[var.index, ])
  }

  ##log transformation
  if (log.scale == FALSE) {
    data <- data
  }

  if (log.scale == "e") {
    data <- log(data + 1)
  }

  if (log.scale != FALSE & log.scale != "e") {
    data <- log(data + 1, as.numeric(log.scale))
  }

  data1 <- SXTscale(data, method = scale.method)
  data1.range <- abs(range(data1))
  dif <- data1.range[1] - data1.range[2]
  if (dif < 0) {
    data1[data1 > data1.range[1]] <- data1.range[1]
  }
  if (dif > 0) {
    data1[data1 < -1 * data1.range[2]] <- -1 * data1.range[2]
  }

  annotation_col <- data.frame(Group = factor(c(group)))

  rownames(annotation_col) <- rownames(data)

  # Specify colors
  ann_col <- NULL
  for (i in seq_along(Group)) {
    ann_col[i] <- color[i]
  }

  ann_colors = list(Group = ann_col)
  names(ann_colors[[1]]) <- Group

  pdf(file.path(path, "heatmap.pdf"),
      width = width,
      height = height)
  par(mar = c(5,5,4,2))
  pheatmap::pheatmap(
    t(data1),
    color = colorRampPalette(c("green", "black", "red"))(1000),
    scale = "none",
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    border_color = border_color,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    fontsize_row = fontsize_row,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    clustering_method = clustering_method,
    ...
  )
  dev.off()
}


SXTdummy <- function (Y) {
  dummy <- matrix(0, nrow = length(Y), ncol = length(table(Y)))
  for (i in seq_along(Y)) {
    for (j in 1:ncol(dummy)) {
      if (Y[i] == names(table(Y))[j])
        dummy[i, j] = 1
    }
  }
  return(dummy)
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

SXTMTmatch <- function(data1,
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





#' @title SXTvip
#' @description Get VIP from pls object.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object PLS object.
#' @return Return VIP.
#' @export
#' @examples
#' library(pls)
#' x <- matrix(rnorm(1000),nrow = 10,ncol = 100)
#' y <- rep(0:1,5)
#' res <- plsr(y~x, method = "oscorespls")
#' SXTvip(res)


SXTvip <- function(object) {
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.
         Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}




#' @title UnivariateTest
#' @description Calculate p value and AUC for each feature.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param test.method Which test you want to use? "t" means stutent t test
#' and "wilcox" mean wilcoxon test.
#' @param adjust.method p value correction method. See p.adjust function.
#' @param log.scale Data transformation method, defaulst is FALSE.
#' @param class Class used to do test.
#' @return MetFlowData which has been added p and AUC information in tags.
#' @seealso The details of univariate test can
#' be found in \code{\link[stats]{t.test}},
#' \code{\link[stats]{p.adjust}} and \code{\link[stats]{wilcox.test}}.
#' @examples
#' @export

setGeneric(
  name = "calP",
  def = function(object,
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
)



VIP <- function(MetFlowData,
                #used data
                log.scale = 10,
                scalemethod="auto",
                plsmethod = "plsreg",
                path = NULL)
  #parameter setting
{
  #
  options(warn = -1)

  if(is.null(path)) {path <- getwd()
  }else{
    dir.create(path)
  }

  subject <- MetFlowData@subject
  qc <- MetFlowData@qc
  subject.info <- MetFlowData@subject.info
  group <- subject.info[,"group"]
  group.unique <- sort(unique(group))
  subject.name <- subject.info[,1]

  if (is.null(qc)) {QC <- FALSE}

  info <- list()
  for (i in seq_along(group.unique)) {
    info[[i]] <- subject.name[which(group == group.unique[i])]
  }

  names(info) <- group.unique

  #load needed packages
  need.packages1 <- c("pls","plsdepot")

  packages <- library()[[2]][,1]
  for (i in need.packages1) {
    if (!any(packages == i)) {install.packages(i)}
  }

  int <- t(subject)
  index <- NULL
  for (i in seq_along(info)) {
    index1 <- as.character(info[[i]])
    index <- c(index,index1)
  }
  if (length(which(index==""))!=0)  {index<-index[-which(index=="")]}

  index <- index[!is.na(index)]
  index <- match(index, rownames(int))
  index <- index[!is.na(index)]
  int <- int[index, ]

  #######
  name <- rownames(int)
  #
  Y <- NULL
  label <- list()
  for (i in seq_along( info )) {
    label[[i]] <- match(as.character(info[[i]]),name)
    label[[i]] <- label[[i]][!is.na(label[[i]])]
    Y[label[[i]]] <- i-1
  }
  #
  int <- log(int+1, log.scale)
  int.scale <- SXTscale(int,method = scalemethod)
  # int.Y<-SXTscale(Y,method=scalemethod)
  int.Y <- Y
  ncompa <- nrow(int) - 1

  if (plsmethod == "plsr") {
    pls1 <- plsr(int.Y~int.scale,scale = FALSE,
                 validation = "CV",ncomp = ncompa,method = "oscorespls")
    save(pls1,file = "pls1")

    #########select the number of compents#################
    msep <- MSEP(pls1)
    save(msep,file = file.path(path, "msep"))
    msep <- msep$val[,,]

    yesnot <- "y"
    while (yesnot == "y"|yesnot == "") {
      comps.number <-readline("How many comps do you want to see? ")
      while (!exists("comps.number")|comps.number=="")
      {cat("You must give a comps number to continute.\n")
        comps.number <- readline("How many comps do you want to see? ")}
      comps.number <- as.numeric(comps.number)
      plot(x = c(1:comps.number),y=msep[1,2:(comps.number+1)],
           type="b",col="firebrick1",pch=20,
           xlab = "ncomp",ylab = "MSEP",cex.lab=1.3,cex.axis=1.3)
      points(x=c(1:comps.number),y=msep[2,2:(comps.number+1)],type="b",pch=2)
      legend("top",legend = c("CV","adjCV"),
             col = c("firebrick1","black"),pch=c(20,2),
             bty = "n",cex = 1.3,pt.cex = 1.3)
      yesnot <- readline("Do you want to see the next plot? (y/n)")
    }

    pdf(file.path(path, "MSEP plot.pdf"))
    plot(x=c(1:comps.number),y=msep[1,2:(comps.number+1)],
         type="b",col="firebrick1",pch=20,
         xlab="ncomp",ylab="MSEP",cex.lab=1.3,cex.axis=1.3)
    points(x=c(1:comps.number),y=msep[2,2:(comps.number+1)],type="b",pch=2)
    legend("top",legend = c("CV","adjCV"),
           col = c("firebrick1","black"),pch=c(20,2),
           bty = "n",cex = 1.3,pt.cex = 1.3)
    dev.off()

    number<-readline("Please type number and
                     press Enter  to continute:  ")
    while (!exists("number")|number=="") {cat("You must give a number to continute.\n")
      number<-
        readline("Please type comps number value and press Enter  to continute: ")}
    number<-as.numeric(number)

    ##################construct final pls model###################
    pls2 <- plsr(int.Y~int.scale, scale = FALSE,validation = "CV",
                 ncomp = number,method = "oscorespls")
    save(pls2, file = file.path(path, "pls2"))
    vip <- SXTvip(pls2)
    vip <- apply(vip, 1, mean)
  }

  else {
    #
    # require(SXTdummy)
    dummy <- SXTdummy(Y)
    # int.dummy<-SXTscale(dummy,method=scalemethod)
    int.dummy <- dummy
    # ncompa = nrow(int.scale) - 1
    ncompa <- min(nrow(int), ncol(int))
    pls1 <- plsdepot::plsreg1(int.scale,Y, comps = ncompa)
    save(pls1,file = file.path(path, "pls1"))
    #########select the number of compents#################
    Q2cum <- pls1$Q2[,5]
    Q2cum[is.nan(Q2cum)] <- 1
    yesnot <- "y"
    while (yesnot=="y"|yesnot=="") {
      comps.number<-readline("How many comps do you want to see? ")
      while (!exists("comps.number")|comps.number=="") {cat("You must give a comps number to continute.\n")
        comps.number<-readline("How many comps do you want to see? ")}
      comps.number<-as.numeric(comps.number)
      barplot(Q2cum[1:comps.number],xlab="ncomp",
              ylab="Q2cum",cex.lab=1.3,cex.axis=1.3)
      a <- barplot(Q2cum[1:comps.number],xlab="ncomp",
                   ylab="Q2cum",cex.lab=1.3,cex.axis=1.3)
      abline(h=0)
      points(a,Q2cum[1:comps.number],type="b",col="red",pch=20,cex=2)
      yesnot <- readline("Do you want to see the next plot? (y/n)")
    }

    number <- readline("Please type number and press Enter  to continute:  ")
    while (!exists("number")|number=="") {cat("You must give a number to continute.\n")
      number <-
        readline("Please type comps number value and press Enter  to continute: ")}
    number <- as.numeric(number)

    ##################construct final pls model###################
    cat(paste("Construct PLS model with all peaks using",number,"comps ...","\n"))
    pls2 <- plsdepot::plsreg1(int.scale,Y,comps = number)
    pls.temp <- plsdepot::plsreg2(int.scale,int.dummy, comps = number)
    vip <- pls.temp$VIP
    vip <- apply(vip, 1, mean)
  }

  tags <- MetFlowData@tags
  tags <- data.frame(tags, vip)
  MetFlowData@tags <- as.matrix(tags)
  return(MetFlowData)
  }

#' @title VolcanoPlot
#' @description Draw volcano plot.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param MetFlowData MetFlowData.
#' @param x x axis factor.
#' @param y y axis factor.
#' @param z z axis factor.
#' @param col Colour for markers and non-markers.
#' @param foldchange.cutoff Fold chagne cutoff.
#' @param p.cutoff p value cutoff.
#' @param vip.cutoff VIP value cutoff.
#' @param path The directory you want to write results.
#' @return volcano plot.
#' @examples
#' \donttest{
#' }

###volcanoplot
setGeneric(
  name = "volcanoPlot",
  def = function(p.value,
                 fc,
                 p.cutoff = 0.05,
                 fc.cutoff = 1.3,
                 xlab = "Log2 (Fold change)",
                 ylab = "-Log10 (P-value)") {
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
)




ChangeSampleName <- function(data = "data.csv",
                             sample.information = "sample.information.csv",
                             polarity = "positive",
                             posfix = NULL,
                             ordered.qc = FALSE,
                             output = TRUE,
                             path = ".") {
  #
  if (path != ".") {
    dir.create(path)
  }
  # data <-
  #   read.csv(file.path(path, data),
  #            stringsAsFactors = FALSE,
  #            check.names = FALSE)

  data <-
    readr::read_csv(file.path(path, data),
                    col_types = readr::cols(),
                    progress = FALSE)

  data <- as.data.frame(data)

  if (sum(duplicated(colnames(data))) > 0) {
    stop("There are duplicated samples (names) in you data.")
  }

  sample.information <-
    read.csv(
      file.path(path, sample.information),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  sample.information <-
    sample.information[!is.na(sample.information[, 1]), ]

  ## sort sample information according to sample order
  sample.information <-
    sample.information[order(as.numeric(sample.information[, 2])), ]
  write.csv(sample.information,
            file.path(path, "sample.information.csv"),
            row.names = FALSE)

  sample.name <- as.character(sample.information[, 1])
  injection.order <- as.numeric(sample.information[, 2])
  class <- sample.information[, 3]
  batch <- sample.information[, 4]

  ## any sample in sample information are not in data?
  index <- match(sample.name, colnames(data))
  sample.not.in.data.index <- which(is.na(index))
  if (length(sample.not.in.data.index) != 0)
  {
    stop(paste(
      paste(sample.name[sample.not.in.data.index], collapse = " "),
      "in sample information are not found in data."
    ))
  }

  ## sort sample in data according to sample order
  sample <- data[, index]
  tags <- data[, -index]
  data <- cbind(tags, sample)
  index <- match(sample.name, colnames(data))

  if (!is.null(posfix)) {
    if (polarity == "positive") {
      sample.name <-
        substr(sample.name,
               start = 1,
               stop = unlist(gregexpr(".POS", sample.name)) - 1)
    }
    if (polarity == "negative") {
      sample.name <-
        substr(sample.name,
               start = 1,
               stop = unlist(gregexpr(".NEG", sample.name)) - 1)
    }
  }

  subject.index <- grep("Subject", class)
  qc.index <- grep("QC", class)

  subject.name <- sample.name[subject.index]
  qc.name <- sample.name[qc.index]

  subject.order <- injection.order[subject.index]
  qc.order <- injection.order[qc.index]

  subject.name <-
    paste(paste("Sample", subject.order, sep = ""), subject.name, sep = "_")
  if (ordered.qc) {
    qc.name <- paste(paste("Sample", qc.order, sep = ""), qc.name, sep = "_")
  }
  else {
    qc.name <- paste("QC", rank(qc.order), sep = "")
    qc.name <-
      paste(paste("Sample", qc.order, sep = ""), qc.name, sep = "_")
  }

  sample.name[subject.index] <- subject.name
  sample.name[qc.index] <- qc.name
  if (polarity == "positive") {
    sample.name <- paste(sample.name, "POS", sep = "_")
  }
  if (polarity == "negative") {
    sample.name <- paste(sample.name, "NEG", sep = "_")
  }
  if (polarity == "none") {
    sample.name <- paste(sample.name, "NONE", sep = "_")
  }

  colnames(data)[index] <- sample.name
  sample.information[, 1] <- sample.name
  write.csv(sample.information,
            file.path(path, "sample.information1.csv"),
            row.names = FALSE)
  if (output) {
    write.csv(data, file.path(path, "data1.csv"), row.names = FALSE)
  }
  return(data)
}










#' @title ChangeWorklist
#' @description Change date in worklist from GetWorklist.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param date Date you want to change in you worklist. Default is today.
#' @return Return a new worklist.
#' @examples
#' \donttest{
#' #demo data
#' batch.design <- paste("A",c(1:100), sep = "")
#' ##create a folder for demo
#' dir.create("Demo")
#' setwd("Demo")
#' write.csv(batch.design, "batch.design.csv", row.names = FALSE)
#' #get worklist
#' GetWorklist()
#' ##remove other information
#' file.remove(dir()[dir()!="worklist POS.csv"])
#' #run ChangeWorklist
#' ChangeWorklist()
#' }


ChangeWorklist <-
  function(date = gsub("-", "", as.character(Sys.Date()))) {
    #used to change the date in worklist
    options(warn = -1)
    #change the date
    file <- dir()
    file <- file[!file.info(file)$isdir]

    worklist <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)



    Data.File <- worklist[, grep("Data.File", colnames(worklist))]
    date.old <- substr(file, 1, 8)
    Data.File <- gsub(date.old, date, Data.File)
    worklist[, grep("Data.File", colnames(worklist))] <- Data.File
    worklistname <- paste(date, substr(file, 9, (nchar(file) - 5)), sep =
                            "")

    write.csv(worklist, sprintf('%s.csv', worklistname), row.names = FALSE)
    cat("The date has been changed.\n")
  }





#' Generate Worklist for data acquisition.
#'
#' @title GetWorklist
#' @description Generate Worklist for data acquisition.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param x batch.design file.
#' @param instrument Which instrument you use?
#' "Agilent" or "AB", default is "Agilent".
#' @param name The name of worklist.
#' @param randommethod Which random method you want to use?
#' "no", "position" or "injection". Default is "no".
#' @param samplenumber Sample number.
#' @param replication Replication times.
#' @param QCstep QC step.
#' @param conditionQCnumber Condition QC number.
#' @param testmixstep Test mixture step.
#' @param injectionfrom Injection order from which? Default is 1.
#' @param user Default is "other".
#' @param dir Directory.
#' @return New worklist.
#' @examples
#' \donttest{
#' #demo data
#' batch.design <- paste("A",c(1:100), sep = "")
#' ##create a folder for demo
#' dir.create("Demo")
#' setwd("Demo")
#' write.csv(batch.design, "batch.design.csv", row.names = FALSE)
#'
#' #run ChangeWorklist
#' GetWorklist()
#' }

GetWorklist <- function(x = NULL,
                        instrument = "Agilent",
                        name = "worklist",
                        randommethod = "no",
                        samplenumber = NULL,
                        replication = 1,
                        QCstep = 8,
                        conditionQCnumber = 10,
                        testmixstep = 0,
                        injectionfrom = 1,
                        user = "other",
                        dir = "D:\\MassHunter\\Data\\SXT\\") {
  #names is the name of the folder,plates is the used plates,
  #if AB,dir is ""

  options(warn = -1)
  file <- dir()

  if (instrument == "AB") {dir = ""}
  if (is.null(x)) {
    x <- read.csv(file[file == "batch.design.csv"],
                  check.names = FALSE,
                  stringsAsFactors = FALSE)
  }
  # --------------------------------------------------------------------------
  options(warn = 0)
  x <- as.character(x[, 1])
  na.number <- sum(is.na(x))
  x <- x[!is.na(x)]
  space.number <- sum(x == "")
  x <- x[x != ""]

  cat(
    paste(
      "\nThere are",
      na.number,
      "NAs and",
      space.number,
      "spaces in your batch design, please confirm there are no error.\n"
    )
  )

  # ---------------------------------------------------------------------------
  if (is.null(samplenumber)) {samplenumber <- length(x)}
  else {
    if (samplenumber > length(x)) {
      samplenumber <- samplenumber
      warning("The sample number you set is larger than
              the sample in your batch design.\n")
    }
    else {
      samplenumber <- samplenumber
    }

    }
  x <- x[1:samplenumber]
  # ------------------------------------------------------------------------
  options(warn = -1)
  plate1 <- rep(1, 51)
  plate2 <- rep(2, 51)
  plate3 <- rep(3, 51)
  plate4 <- rep(4, 51)
  plate5 <- rep(5, 51)
  plate6 <- rep(6, 51)

  plate <- rep(c(plate1, plate2), 6)
  plate <- plate[1:(samplenumber * replication)]
  real.plate <-
    rep(c(plate1, plate2, plate3, plate4, plate5, plate6), 6)
  real.plate <- real.plate[1:(samplenumber * replication)]
  vial.position <- rep(c(1:51), 6)
  vial.position <- vial.position[1:(samplenumber * replication)]

  sub.position <- list()
  for (i in 1:6) {
    sub.position[[i]] <- paste(LETTERS[i], c(1:9), sep = "")
  }
  sub.position <- unlist(sub.position)

  real.position <- list()
  for (i in 1:12) {
    real.position[[i]] <-
      paste(paste("P", i, sep = ""), sub.position, sep = "-")
  }

  real.position <- unlist(real.position)

  position <- rep(real.position[1:108], 6)

  position <- position[1:(samplenumber * replication)]
  real.position <- real.position[1:(samplenumber * replication)]
  sub.position96 <-
    c(
      paste("A", c(1:12), sep = ""),
      paste("B", c(1:12), sep = ""),
      paste("C", c(1:12), sep = ""),
      paste("D", c(1:12), sep = ""),
      paste("E", c(1:12), sep = ""),
      paste("F", c(1:12), sep = ""),
      paste("G", c(1:12), sep = ""),
      paste("H", c(1:12), sep = "")
    )
  real.position96 <-
    rep(c(
      paste("P1", sub.position96, sep = "-"),
      paste("P2", sub.position96, sep = "-"),
      paste("P3", sub.position96, sep = "-"),
      paste("P4", sub.position96, sep = "-"),
      paste("P5", sub.position96, sep = "-"),
      paste("P6", sub.position96, sep = "-")
    ), 4)
  real.position96 <- real.position96[1:(samplenumber * replication)]

  position96 <-
    rep(c(
      paste("P1", sub.position96, sep = "-"),
      paste("P2", sub.position96, sep = "-")
    ), 8)
  position96 <- position96[1:(samplenumber * replication)]

  # -------------------------------------------------------------------------

  ###repeat samples?
  if (replication == 1) {
    x <- x
  } else{
    x2 <- NULL
    for (i in seq_len(replication)) {
      x1 <- paste(x, i, sep = "_")
      x2 <- cbind(x2, x1)
    }
    x <- x2
  }

  x <- as.character(x)
  #random position or random injection order or no
  if (instrument == "Agilent") {
    if (randommethod == "position") {
      random.order <- sample(1:(samplenumber * replication))
      x <- data.frame(random.order, x)
      x <- x[order(as.numeric(x[, 1])),]
      x <- as.character(x[, -1])
      x <-
        cbind(x, position, real.position, position96, real.position96)
      colnames(x) <-
        c(
          "Sample Name",
          "Position in 54 plate",
          "Real position in 54 plate",
          "Position in 96 plate",
          "Real position in 96 plate"
        )
      write.csv(x, sprintf("%s sample info.csv", name))
      x <- x[, -c(3, 4, 5)]
    }
    if (randommethod == "injection") {
      if (length(x) > 108)
      {
        warning("The sample number is larger than 108,
                injection order random is not commended.")
      }
      x <-
        cbind(x, position, real.position, position96, real.position96)
      colnames(x) <-
        c(
          "Sample Name",
          "Position in 54 plate",
          "Real position in 54 plate",
          "Position in 96 plate",
          "Real position in 96 plate"
        )
      write.csv(x, sprintf("%s sample info.csv", name))
      random.order <- sample(1:(samplenumber * replication))
      x <- data.frame(random.order, x)
      x <- x[order(x[, 1]),]
      x <- x[, -c(1, 4, 5, 6)]
      }
    if (randommethod == "no") {
      x <- cbind(x, position, real.position, position96, real.position96)
      colnames(x) <-
        c(
          "Sample Name",
          "Position in 54 plate",
          "Real position in 54 plate",
          "Position in 96 plate",
          "Real position in 96 plate"
        )
      write.csv(x, sprintf("%s sample info.csv", name))
      x <- x[, -c(3, 4, 5)]
    }
  }
  ##AB instrument
  if (instrument == "AB") {
    if (randommethod == "position") {
      random.order <- sample(1:(samplenumber * replication))
      x <- data.frame(random.order, x)
      x <- x[order(as.numeric(x[, 1])),]
      x <- x[, -1]
      x <-
        cbind(x,
              plate,
              vial.position,
              real.plate,
              position96,
              real.position96)
      colnames(x) <-
        c(
          "Sample Name",
          "Plate",
          "Position in 54 plate",
          "Real plate of 54 plate",
          "Position in 96 plate",
          "Real position 96 plate"
        )
      write.csv(x, "sample info.csv")
      x <- x[, -c(4, 5, 6)]
    }

    if (randommethod == "injection") {
      x <-
        cbind(x,
              plate,
              vial.position,
              real.plate,
              position96,
              real.position96)
      colnames(x) <-
        c(
          "Sample Name",
          "Plate",
          "Position in 54 plate",
          "Real Plate of 54 plate",
          "Position in 96 plate",
          "Real position in 96 plate"
        )
      write.csv(x, "sample info.csv")
      random.order <- sample(1:(samplenumber * replication))
      x <- cbind(random.order, x)
      x <- x[order(as.numeric(x[, 1])), ]
      x <- x[, -c(1, 5, 6, 7)]
    }

    if (randommethod == "no") {
      x <-
        data.frame(x,
                   plate,
                   vial.position,
                   real.plate,
                   position96,
                   real.position96)
      colnames(x) <-
        c(
          "Sample Name",
          "Plate",
          "Position in 54 plate",
          "Real Plate of 54 plate",
          "Position in 96 plate",
          "Real position in 96 plate"
        )
      write.csv(x, "sample info.csv")
      x <- x[, -c(4, 5, 6)]
    }
  }
  #now x column 1 is Sample Name, column 2 is Sample Position
  if (instrument == "Agilent") {
    Blank <- c("Blank", "Vial1")
    Test.mix <- matrix(c("Test_mix", "Vial2"), ncol = 2)
    QC <- c("QC", "Vial3")
    Blank.QC <- rbind(Blank, QC)
  }
  else {
    Blank <- c("Blank", "1", "52")
    Test.mix <- c("Test_mix", "1", "53")
    QC <- c("QC", "1", "54")
    Blank.QC <- rbind(Blank, QC)
  }

  #insert Blank and QC
  x <-
    lapply(seq(1, nrow(x), by = QCstep), function(y)
      if (y + QCstep - 1 <= (samplenumber * replication)) {
        x[y:(y + QCstep - 1),]
      }
      else {
        x[y:nrow(x),]
      })
  colnames(Blank.QC) <- colnames(x[[1]])
  x <- lapply(x, function(y)
    rbind(Blank.QC, y))
  ###
  x2 <- NULL
  for (i in seq_along(x)) {
    x1 <- x[[i]]
    x2 <- rbind(x2, x1)
  }
  x <- x2
  x <- x[-1, ]

  x <- rbind(x, Blank.QC)
  x <- x[-(nrow(x) - 1), ]
  #insert Test.mix
  if (testmixstep == 0) {
    x = x
  }
  else {
    x <-
      lapply(seq(1, nrow(x), by = testmixstep), function(y)
        if (y + testmixstep - 1 <= nrow(x)) {
          x[y:(y + testmixstep - 1),]
        } else {
          x[y:nrow(x),]
        })

    colnames(Test.mix) <- colnames(x[[1]])
    x <- lapply(x, function(y)
      rbind(Test.mix, y))

    x3 <- NULL
    for (i in seq_along(x)) {
      x1 <- x[[i]]
      x3 <- rbind(x3, x1)
    }
    x <- x3
    x <- rbind(x, Test.mix)
  }

  if (instrument == "Agilent") {
    colnames(x) <- c('Sample.Name', "Sample.Position")
    x[, 1] <- as.character(x[, 1])
    x[, 2] <- as.character(x[, 2])
  }
  if (instrument == "AB") {
    colnames(x) <- c('Sample.Name', "Plate.Position", "Vial.Postion")
    x[, 1] <- as.character(x[, 1])
    x[, 2] <- as.character(x[, 2])
    x[, 3] <- as.character(x[, 3])
  }

  if (instrument == "Agilent")
  {
    temp1 <- matrix(rep(Blank, 3), ncol = 2, byrow = TRUE)
    temp2 <- matrix(rep(QC, conditionQCnumber),
                    ncol = 2,
                    byrow = TRUE)
    temp3 <- matrix(rep(Blank, 3), ncol = 2, byrow = TRUE)
    colnames(temp1) <-
      colnames(temp2) <- colnames(temp3) <- colnames(x)
    x <- rbind(temp1, temp2, x, temp3)
    x[, 1] <- as.character(x[, 1])
    x[, 2] <- as.character(x[, 2])
  }
  if (instrument == "AB") {
    temp1 <- matrix(rep(Blank, 3), ncol = 3, byrow = TRUE)
    temp2 <- matrix(rep(QC, conditionQCnumber),
                    ncol = 3,
                    byrow = TRUE)
    temp3 <- matrix(rep(Blank, 3), ncol = 3, byrow = TRUE)
    colnames(temp1) <-
      colnames(temp2) <- colnames(temp3) <- colnames(x)
    x <- rbind(temp1, temp2, x, temp3)
    x[, 1] <- as.character(x[, 1])
    x[, 2] <- as.character(x[, 2])
    x[, 3] <- as.character(x[, 3])
  }

  Blank.number <- length(grep("Blank", x[, 1]))
  x[, 1][grep("Blank", x[, 1])] <-
    paste("Blank", c(1:Blank.number), sep = "")

  Test.mix.number <- length(grep("Test_mix", x[, 1]))
  x[, 1][grep("Test_mix", x[, 1])] <-
    paste("Test_mix", c(1:Test.mix.number), sep = "")

  QC.number <- length(grep("QC", x[, 1]))
  x[, 1][grep("QC", x[, 1])][1:conditionQCnumber] <-
    paste("Condition_QC", c(1:conditionQCnumber), sep = "")
  x[, 1][grep("QC", x[, 1])][conditionQCnumber + 1:QC.number] <-
    paste("QC", c(1:(QC.number - conditionQCnumber)), sep = "")
  first <- which(x[, 1] == "QC1")
  last <-
    which(x[, 1] == sprintf("QC%s", QC.number - conditionQCnumber))

  before.info <- x[1:(first - 1),]
  Data.File1 <- before.info[, 1]

  after.info <- x[(last + 1):nrow(x),]
  Data.File5 <- after.info[, 1]

  middle.info <- x[first:last, ]
  middle.info <- cbind(middle.info, c(1:nrow(middle.info)))
  Sample.QC <-
    middle.info[setdiff(1:nrow(middle.info),grep("Blank", middle.info[, 1])), ]
  #remove Blank in middle.info
  Sample.QC <-
    Sample.QC[setdiff(1:nrow(Sample.QC), grep("Test_mix", Sample.QC[, 1])), ]
  #remove Test.mix in middle.info

  middle.blank <-
    middle.info[grep("Blank", middle.info[, 1]), , drop = FALSE]
  #column 3 is number
  middle.testmix <-
    middle.info[grep("Test_mix", middle.info[, 1]), , drop = FALSE]

  Data.File2 <-
    paste("Sample", c(injectionfrom:(nrow(Sample.QC) + injectionfrom - 1)),
          sep = "")

  Data.File2 <- paste(Data.File2, Sample.QC[, 1], sep = "_")

  if (user == "other") {
    Data.File2 <- Sample.QC[, 1]
  }

  if (instrument == "Agilent") {
    Data.File2 <- cbind(Data.File2, Sample.QC[, 3])
    Data.File3 <- middle.blank[, c(1, 3)]
    Data.File4 <- middle.testmix[, c(1, 3)]
  }
  if (instrument == "AB") {
    Data.File2 <- cbind(Data.File2, Sample.QC[, 4])
    Data.File3 <- middle.blank[, c(1, 4)]
    Data.File4 <- middle.testmix[, c(1, 4)]
  }

  colnames(Data.File4) <- colnames(Data.File3) <-
    colnames(Data.File2) <- paste("test", c(1:ncol(Data.File2)))
  Data.File2 <- rbind(Data.File2, Data.File3, Data.File4)
  Data.File2 <- Data.File2[order(as.numeric(Data.File2[, 2])),]
  Data.File2 <- Data.File2[, -2]

  Data.File <- c(Data.File1, as.character(Data.File2), Data.File5)
  name.POS <- paste(name, "POS")
  name.NEG <- paste(name, "NEG")
  Data.File.POS <- paste(dir, name.POS, "\\", Data.File, sep = "")
  Data.File.NEG <- paste(dir, name.NEG, "\\", Data.File, sep = "")
  Data.File.POS <- paste(Data.File.POS, "POS", sep = "_")
  Data.File.NEG <- paste(Data.File.NEG, "NEG", sep = "_")

  if (instrument == "Agilent")
  {
    Data.File.POS <- paste(Data.File.POS, "d", sep = ".")
    Data.File.NEG <- paste(Data.File.NEG, "d", sep = ".")
  }

  x.POS <- cbind(x, Data.File.POS)
  x.NEG <- cbind(x, Data.File.NEG)

  write.csv(x.POS, sprintf("%s POS.csv", name), row.names = FALSE)
  write.csv(x.NEG, sprintf("%s NEG.csv", name), row.names = FALSE)
  cat("Worklist is generated.\n")
  return(TRUE)
  }


setGeneric(name = "SXTMTmatch2",
           def = function(data1,
                          data2,
                          mz.tol,
                          #rt.tol is relative
                          rt.tol = 30,
                          rt.error.type = c("relative", "abs")){
             rt.error.type <- match.arg(rt.error.type)
             #
             if (nrow(data1) == 0 | nrow(data2) == 0) {
               result <- NULL
               return(result)
             }
             # mz1 <- as.numeric(data1[, 1])
             # rt1 <- as.numeric(data1[, 2])
             info1 <- data1[,c(1,2)]
             info1 <- apply(info1, 1, list)

             mz2 <- as.numeric(data2[, 1])
             rt2 <- as.numeric(data2[, 2])

             result <- pbapply::pblapply(info1, function(x) {
               temp.mz1 <- x[[1]][[1]]
               temp.rt1 <- x[[1]][[2]]
               mz.error <- abs(temp.mz1 - mz2) * 10 ^ 6 / temp.mz1
               if(rt.error.type == "relative"){
                 rt.error <- abs(temp.rt1 - rt2) * 100 / temp.rt1
               }else{
                 rt.error <- abs(temp.rt1 - rt2)
               }

               j <- which(mz.error <= mz.tol & rt.error <= rt.tol)
               if(length(j) == 0){
                 matrix(NA, ncol = 7)
               }else{
                 cbind(j, temp.mz1, mz2[j], mz.error[j], temp.rt1, rt2[j], rt.error[j])
               }
             })

             if(length(result) == 1){
               result <- cbind(1,result[[1]])
             }else{
               result <- mapply(function(x,y){list(cbind(x,y))},
                                x <- 1:length(info1),
                                y = result)
               result <- do.call(rbind, result)
             }

             result <- matrix(result[which(!apply(result,1,function(x) any(is.na(x)))),], ncol = 8)
             if(nrow(result) == 0) return(NULL)
             colnames(result) <-
               c("Index1",
                 "Index2",
                 "mz1",
                 "mz2",
                 "mz error",
                 "rt1",
                 "rt2",
                 "rt error")
             result <- result
           })




setGeneric(name = "sxtScale", 
           def = function(df,
                         method = c("no", "auto", "pareto", "center")){
             
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
           })


