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


# #' @title HeatMap
# #' @description Heat map
# #' @author Xiaotao Shen
# #' \email{shenxt@@sioc.ac.cn}
# #' @param MetFlowData MetFlowData.
# #' @param log.scale log transformation or not.
# #' @param color Color list for sample group.
# #' @param variable "all" or "marker" for heatmap.
# #' @param Group group for heatmap.
# #' @param scale.method scale method.
# #' @param show_rownames Default is FALSE.
# #' @param show_colnames Default is FALSE.
# #' @param path Work directory.
# #' @param width Plot width
# #' @param height Plot height.
# #' @param border_color Default is NA.
# #' @param fontsize_row Default is 10.
# #' @param cluster_rows Default is TRUE.
# #' @param cluster_cols Default is TURE.
# #' @param clustering_method Default is"ward.D",
# #' @param ... other parameters for pheatmap.
# #' @return A heatmap plot.
# #' @seealso \code{\link[pheatmap]{pheatmap}}
#' 
#' HeatMap <- function(MetFlowData,
#'                     log.scale = FALSE,
#'                     color = c("palegreen",
#'                               "firebrick1",
#'                               "royalblue",
#'                               "yellow",
#'                               "black",
#'                               "cyan",
#'                               "gray48"),
#'                     variable = "all",
#'                     Group = c("control", "case"),
#'                     scale.method = "auto",
#'                     show_rownames = FALSE,
#'                     show_colnames = FALSE,
#'                     path = ".",
#'                     width = 7,
#'                     height = 7,
#'                     border_color = NA,
#'                     fontsize_row = 10,
#'                     cluster_rows = TRUE,
#'                     cluster_cols = TRUE,
#'                     clustering_method = "ward.D",
#'                     ...) {
#'   if (path != ".") {
#'     dir.create(path)
#'   }
#' 
#'   subject <- MetFlowData@subject
#'   tags <- MetFlowData@tags
#'   subject.info <- MetFlowData@subject.info
#'   group <- subject.info[, "group"]
#' 
#'   idx <- which(group %in% Group)
#'   subject.info <- subject.info[idx, ]
#'   subject <- subject[, idx]
#'   group <- subject.info[, "group"]
#'   ## data organization
#'   if (variable == "all") {
#'     data <- t(subject)
#'   } else{
#'     if (all(colnames(tags) != "is.marker")) {
#'       stop("Please select marker first.")
#'     }
#'     is.marker <- tags[, "is.marker"]
#'     var.index <- which(is.marker == "yes")
#'     data <- t(subject[var.index, ])
#'   }
#' 
#'   ##log transformation
#'   if (log.scale == FALSE) {
#'     data <- data
#'   }
#' 
#'   if (log.scale == "e") {
#'     data <- log(data + 1)
#'   }
#' 
#'   if (log.scale != FALSE & log.scale != "e") {
#'     data <- log(data + 1, as.numeric(log.scale))
#'   }
#' 
#'   data1 <- SXTscale(data, method = scale.method)
#'   data1.range <- abs(range(data1))
#'   dif <- data1.range[1] - data1.range[2]
#'   if (dif < 0) {
#'     data1[data1 > data1.range[1]] <- data1.range[1]
#'   }
#'   if (dif > 0) {
#'     data1[data1 < -1 * data1.range[2]] <- -1 * data1.range[2]
#'   }
#' 
#'   annotation_col <- data.frame(Group = factor(c(group)))
#' 
#'   rownames(annotation_col) <- rownames(data)
#' 
#'   # Specify colors
#'   ann_col <- NULL
#'   for (i in seq_along(Group)) {
#'     ann_col[i] <- color[i]
#'   }
#' 
#'   ann_colors = list(Group = ann_col)
#'   names(ann_colors[[1]]) <- Group
#' 
#'   pdf(file.path(path, "heatmap.pdf"),
#'       width = width,
#'       height = height)
#'   par(mar = c(5,5,4,2))
#'   pheatmap::pheatmap(
#'     t(data1),
#'     color = colorRampPalette(c("green", "black", "red"))(1000),
#'     scale = "none",
#'     show_rownames = show_rownames,
#'     show_colnames = show_colnames,
#'     border_color = border_color,
#'     annotation_col = annotation_col,
#'     annotation_colors = ann_colors,
#'     fontsize_row = fontsize_row,
#'     cluster_rows = cluster_rows,
#'     cluster_cols = cluster_cols,
#'     clustering_method = clustering_method,
#'     ...
#'   )
#'   dev.off()
#' }


# SXTdummy <- function(Y) {
#   dummy <- matrix(0, nrow = length(Y), ncol = length(table(Y)))
#   for (i in seq_along(Y)) {
#     for (j in 1:ncol(dummy)) {
#       if (Y[i] == names(table(Y))[j])
#         dummy[i, j] = 1
#     }
#   }
#   return(dummy)
# }


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





#' #' @title SXTvip
#' #' @description Get VIP from pls object.
#' #' @author Xiaotao Shen
#' #' \email{shenxt@@sioc.ac.cn}
#' #' @param object PLS object.
#' #' @return Return VIP.
#' #' @export
#' #' @examples
#' #' library(pls)
#' #' x <- matrix(rnorm(1000),nrow = 10,ncol = 100)
#' #' y <- rep(0:1,5)
#' #' res <- plsr(y~x, method = "oscorespls")
#' #' SXTvip(res)
#' 
#' 
#' SXTvip <- function(object) {
#'   if (object$method != "oscorespls")
#'     stop("Only implemented for orthogonal scores algorithm.
#'          Refit with 'method = \"oscorespls\"'")
#'   if (nrow(object$Yloadings) > 1)
#'     stop("Only implemented for single-response models")
#'   SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
#'   Wnorm2 <- colSums(object$loading.weights^2)
#'   SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
#'   sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
#' }




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


# setGeneric(name = "calP",
#            def = function(object,
#                           control.group,
#                           case.group,
#                           test.method = c("t", "wilcox"),
#                           adjust.method = c("holm",
#                                             "hochberg",
#                                             "hommel",
#                                             "bonferroni",
#                                             "BH",
#                                             "BY",
#                                             "fdr",
#                                             "none"),
#                           log.scale = FALSE,
#                           ...)
#            {
#              standardGeneric("calP")
#            }
# )
# 
# setMethod(f = "calP", 
#           signature(object = "metflowClass"), 
#           function(object,
#                    control.group,
#                    case.group,
#                    test.method = c("t", "wilcox"),
#                    adjust.method = c("holm",
#                                      "hochberg",
#                                      "hommel",
#                                      "bonferroni",
#                                      "BH",
#                                      "BY",
#                                      "fdr",
#                                      "none"),
#                    log.scale = FALSE,
#                    ...) {
#             test.method <- match.arg(test.method)
#             adjust.method <- match.arg(adjust.method)
#             if(missing(control.group) | missing(case.group)){
#               stop("Please set control.group or case.group.\n")
#             }
#             
#             if(!all(c(control.group, case.group) %in% object@sample.info$group)){
#               stop("Please make sure that the control and case group are in your sample information.\n")
#             }
#             
#             ms1_data <- object@ms1.data
#             if(length(ms1_data) > 1){
#               stop("Please align batches first.\n")  
#             }
#             
#             ms1_data <- ms1_data[[1]]
#             sample_info <- object@sample.info
#             
#             control_data <- 
#               sample_info %>% 
#               filter(., group == control.group) %>% 
#               dplyr::pull(., sample.name) %>% 
#               match(., colnames(ms1_data)) %>% 
#               ms1_data[,.]
#             
#             case_data <- 
#               sample_info %>% 
#               filter(., group == case.group) %>% 
#               dplyr::pull(., sample.name) %>% 
#               match(., colnames(ms1_data)) %>% 
#               ms1_data[,.]
#             
#             if(sum(is.na(control_data)) != 0 | sum(is.na(case_data)) != 0){
#               stop("Please impute MV first.\n")
#             }
#             
#             if(log.scale == TRUE){
#               control_data <- log(control_data + 1, 10)
#               case_data <- log(case_data + 1, 10)
#             }
#             
#             control_data <- split(control_data, seq(nrow(control_data)))
#             case_data <- split(case_data, seq(nrow(case_data)))
#             pbapply::pboptions(type = "timer")
#             
#             p_value <- pbapply::pbmapply(function(x, y){
#               temp_p <- try(
#                 switch(test.method,
#                        t = t.test(as.numeric(x), as.numeric(y), ...),
#                        wilcox = wilcox.test(as.numeric(x), as.numeric(y), ...))
#               )
#               temp_p$p.value
#             },
#             x = control_data,
#             y = case_data)
#             
#             p_value <- 
#               object %>% 
#               getData(., "Tags") %>% 
#               dplyr::pull(., name) %>% 
#               data.frame(index = 1:length(p_value), ., p_value)
#             colnames(p_value)[2] <- "name"
#             
#             if (adjust.method != "none") {
#               p_value <- p_value %>%
#                 dplyr::mutate(., p_value_adjust = p.adjust(p_value,
#                                                            method = adjust.method))
#             }
#             return(p_value)
# })


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



# VIP <- function(MetFlowData,
#                 #used data
#                 log.scale = 10,
#                 scalemethod="auto",
#                 plsmethod = "plsreg",
#                 path = NULL)
#   #parameter setting
# {
#   #
#   options(warn = -1)
# 
#   if(is.null(path)) {path <- getwd()
#   }else{
#     dir.create(path)
#   }
# 
#   subject <- MetFlowData@subject
#   qc <- MetFlowData@qc
#   subject.info <- MetFlowData@subject.info
#   group <- subject.info[,"group"]
#   group.unique <- sort(unique(group))
#   subject.name <- subject.info[,1]
# 
#   if (is.null(qc)) {QC <- FALSE}
# 
#   info <- list()
#   for (i in seq_along(group.unique)) {
#     info[[i]] <- subject.name[which(group == group.unique[i])]
#   }
# 
#   names(info) <- group.unique
# 
#   #load needed packages
#   need.packages1 <- c("pls","plsdepot")
# 
#   packages <- library()[[2]][,1]
#   for (i in need.packages1) {
#     if (!any(packages == i)) {install.packages(i)}
#   }
# 
#   int <- t(subject)
#   index <- NULL
#   for (i in seq_along(info)) {
#     index1 <- as.character(info[[i]])
#     index <- c(index,index1)
#   }
#   if (length(which(index==""))!=0)  {index<-index[-which(index=="")]}
# 
#   index <- index[!is.na(index)]
#   index <- match(index, rownames(int))
#   index <- index[!is.na(index)]
#   int <- int[index, ]
# 
#   #######
#   name <- rownames(int)
#   #
#   Y <- NULL
#   label <- list()
#   for (i in seq_along( info )) {
#     label[[i]] <- match(as.character(info[[i]]),name)
#     label[[i]] <- label[[i]][!is.na(label[[i]])]
#     Y[label[[i]]] <- i-1
#   }
#   #
#   int <- log(int+1, log.scale)
#   int.scale <- SXTscale(int,method = scalemethod)
#   # int.Y<-SXTscale(Y,method=scalemethod)
#   int.Y <- Y
#   ncompa <- nrow(int) - 1
# 
#   if (plsmethod == "plsr") {
#     pls1 <- plsr(int.Y~int.scale,scale = FALSE,
#                  validation = "CV",ncomp = ncompa,method = "oscorespls")
#     save(pls1,file = "pls1")
# 
#     #########select the number of compents#################
#     msep <- MSEP(pls1)
#     save(msep,file = file.path(path, "msep"))
#     msep <- msep$val[,,]
# 
#     yesnot <- "y"
#     while (yesnot == "y"|yesnot == "") {
#       comps.number <-readline("How many comps do you want to see? ")
#       while (!exists("comps.number")|comps.number=="")
#       {cat("You must give a comps number to continute.\n")
#         comps.number <- readline("How many comps do you want to see? ")}
#       comps.number <- as.numeric(comps.number)
#       plot(x = c(1:comps.number),y=msep[1,2:(comps.number+1)],
#            type="b",col="firebrick1",pch=20,
#            xlab = "ncomp",ylab = "MSEP",cex.lab=1.3,cex.axis=1.3)
#       points(x=c(1:comps.number),y=msep[2,2:(comps.number+1)],type="b",pch=2)
#       legend("top",legend = c("CV","adjCV"),
#              col = c("firebrick1","black"),pch=c(20,2),
#              bty = "n",cex = 1.3,pt.cex = 1.3)
#       yesnot <- readline("Do you want to see the next plot? (y/n)")
#     }
# 
#     pdf(file.path(path, "MSEP plot.pdf"))
#     plot(x=c(1:comps.number),y=msep[1,2:(comps.number+1)],
#          type="b",col="firebrick1",pch=20,
#          xlab="ncomp",ylab="MSEP",cex.lab=1.3,cex.axis=1.3)
#     points(x=c(1:comps.number),y=msep[2,2:(comps.number+1)],type="b",pch=2)
#     legend("top",legend = c("CV","adjCV"),
#            col = c("firebrick1","black"),pch=c(20,2),
#            bty = "n",cex = 1.3,pt.cex = 1.3)
#     dev.off()
# 
#     number<-readline("Please type number and
#                      press Enter  to continute:  ")
#     while (!exists("number")|number=="") {cat("You must give a number to continute.\n")
#       number<-
#         readline("Please type comps number value and press Enter  to continute: ")}
#     number<-as.numeric(number)
# 
#     ##################construct final pls model###################
#     pls2 <- plsr(int.Y~int.scale, scale = FALSE,validation = "CV",
#                  ncomp = number,method = "oscorespls")
#     save(pls2, file = file.path(path, "pls2"))
#     vip <- SXTvip(pls2)
#     vip <- apply(vip, 1, mean)
#   }
# 
#   else {
#     #
#     # require(SXTdummy)
#     dummy <- SXTdummy(Y)
#     # int.dummy<-SXTscale(dummy,method=scalemethod)
#     int.dummy <- dummy
#     # ncompa = nrow(int.scale) - 1
#     ncompa <- min(nrow(int), ncol(int))
#     pls1 <- plsdepot::plsreg1(int.scale,Y, comps = ncompa)
#     save(pls1,file = file.path(path, "pls1"))
#     #########select the number of compents#################
#     Q2cum <- pls1$Q2[,5]
#     Q2cum[is.nan(Q2cum)] <- 1
#     yesnot <- "y"
#     while (yesnot=="y"|yesnot=="") {
#       comps.number<-readline("How many comps do you want to see? ")
#       while (!exists("comps.number")|comps.number=="") {cat("You must give a comps number to continute.\n")
#         comps.number<-readline("How many comps do you want to see? ")}
#       comps.number<-as.numeric(comps.number)
#       barplot(Q2cum[1:comps.number],xlab="ncomp",
#               ylab="Q2cum",cex.lab=1.3,cex.axis=1.3)
#       a <- barplot(Q2cum[1:comps.number],xlab="ncomp",
#                    ylab="Q2cum",cex.lab=1.3,cex.axis=1.3)
#       abline(h=0)
#       points(a,Q2cum[1:comps.number],type="b",col="red",pch=20,cex=2)
#       yesnot <- readline("Do you want to see the next plot? (y/n)")
#     }
# 
#     number <- readline("Please type number and press Enter  to continute:  ")
#     while (!exists("number")|number=="") {cat("You must give a number to continute.\n")
#       number <-
#         readline("Please type comps number value and press Enter  to continute: ")}
#     number <- as.numeric(number)
# 
#     ##################construct final pls model###################
#     cat(paste("Construct PLS model with all peaks using",number,"comps ...","\n"))
#     pls2 <- plsdepot::plsreg1(int.scale,Y,comps = number)
#     pls.temp <- plsdepot::plsreg2(int.scale,int.dummy, comps = number)
#     vip <- pls.temp$VIP
#     vip <- apply(vip, 1, mean)
#   }
# 
#   tags <- MetFlowData@tags
#   tags <- data.frame(tags, vip)
#   MetFlowData@tags <- as.matrix(tags)
#   return(MetFlowData)
#   }

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




# ChangeSampleName <- function(data = "data.csv",
#                              sample.information = "sample.information.csv",
#                              polarity = "positive",
#                              posfix = NULL,
#                              ordered.qc = FALSE,
#                              output = TRUE,
#                              path = ".") {
#   #
#   if (path != ".") {
#     dir.create(path)
#   }
#   # data <-
#   #   read.csv(file.path(path, data),
#   #            stringsAsFactors = FALSE,
#   #            check.names = FALSE)
# 
#   data <-
#     readr::read_csv(file.path(path, data),
#                     col_types = readr::cols(),
#                     progress = FALSE)
# 
#   data <- as.data.frame(data)
# 
#   if (sum(duplicated(colnames(data))) > 0) {
#     stop("There are duplicated samples (names) in you data.")
#   }
# 
#   sample.information <-
#     read.csv(
#       file.path(path, sample.information),
#       stringsAsFactors = FALSE,
#       check.names = FALSE
#     )
#   sample.information <-
#     sample.information[!is.na(sample.information[, 1]), ]
# 
#   ## sort sample information according to sample order
#   sample.information <-
#     sample.information[order(as.numeric(sample.information[, 2])), ]
#   write.csv(sample.information,
#             file.path(path, "sample.information.csv"),
#             row.names = FALSE)
# 
#   sample.name <- as.character(sample.information[, 1])
#   injection.order <- as.numeric(sample.information[, 2])
#   class <- sample.information[, 3]
#   batch <- sample.information[, 4]
# 
#   ## any sample in sample information are not in data?
#   index <- match(sample.name, colnames(data))
#   sample.not.in.data.index <- which(is.na(index))
#   if (length(sample.not.in.data.index) != 0)
#   {
#     stop(paste(
#       paste(sample.name[sample.not.in.data.index], collapse = " "),
#       "in sample information are not found in data."
#     ))
#   }
# 
#   ## sort sample in data according to sample order
#   sample <- data[, index]
#   tags <- data[, -index]
#   data <- cbind(tags, sample)
#   index <- match(sample.name, colnames(data))
# 
#   if (!is.null(posfix)) {
#     if (polarity == "positive") {
#       sample.name <-
#         substr(sample.name,
#                start = 1,
#                stop = unlist(gregexpr(".POS", sample.name)) - 1)
#     }
#     if (polarity == "negative") {
#       sample.name <-
#         substr(sample.name,
#                start = 1,
#                stop = unlist(gregexpr(".NEG", sample.name)) - 1)
#     }
#   }
# 
#   subject.index <- grep("Subject", class)
#   qc.index <- grep("QC", class)
# 
#   subject.name <- sample.name[subject.index]
#   qc.name <- sample.name[qc.index]
# 
#   subject.order <- injection.order[subject.index]
#   qc.order <- injection.order[qc.index]
# 
#   subject.name <-
#     paste(paste("Sample", subject.order, sep = ""), subject.name, sep = "_")
#   if (ordered.qc) {
#     qc.name <- paste(paste("Sample", qc.order, sep = ""), qc.name, sep = "_")
#   }
#   else {
#     qc.name <- paste("QC", rank(qc.order), sep = "")
#     qc.name <-
#       paste(paste("Sample", qc.order, sep = ""), qc.name, sep = "_")
#   }
# 
#   sample.name[subject.index] <- subject.name
#   sample.name[qc.index] <- qc.name
#   if (polarity == "positive") {
#     sample.name <- paste(sample.name, "POS", sep = "_")
#   }
#   if (polarity == "negative") {
#     sample.name <- paste(sample.name, "NEG", sep = "_")
#   }
#   if (polarity == "none") {
#     sample.name <- paste(sample.name, "NONE", sep = "_")
#   }
# 
#   colnames(data)[index] <- sample.name
#   sample.information[, 1] <- sample.name
#   write.csv(sample.information,
#             file.path(path, "sample.information1.csv"),
#             row.names = FALSE)
#   if (output) {
#     write.csv(data, file.path(path, "data1.csv"), row.names = FALSE)
#   }
#   return(data)
# }





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



#' @title showPeakPlot
#' @description Show peak plot.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object1 metflowClass object 1.
#' @param object2 metflowClass object 2.
#' @param peak.index Peak index
#' @return Peak plot (ggplot2 object).
#' @export
#' @import tidyverse
#' @import patchwork

setGeneric(
  name = "showPeakPlot",
  def = function(object1, object2, peak.index = 1) {
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
)








