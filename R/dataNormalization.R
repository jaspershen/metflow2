#' @title normalizeData
#' @description Data normalization.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object A metflowClass object.
#' @param method Normalization method, mean, median, total svr or loess,
#' default is svr. Please see the details.
#' @param keep.scale Remain scale or not. Default is TRUE.
#' @param begin 0.5
#' @param end 1
#' @param step 0.2
#' @param threads 4
#' @export
#' @import tidyverse
#' @return A new metflowClass object.
#' @examples

setGeneric(
  name = "normalizeData",
  def = function(object,
                 method = c("svr", "total", "median", "mean", "pqn", "loess"),
                 keep.scale = TRUE,
                 begin = 0.5,
                 end = 1, 
                 step = 0.2, 
                 threads = 4) {
    method <- match.arg(method)
    if (class(object) != "metflowClass") {
      stop("Only the metflowClass is supported!\n")
    }
    
    ms1_data <- object@ms1.data
    if(length(ms1_data) > 1){
      stop("Please align batch first.\n")
    }
    
    ms1_data <- ms1_data[[1]]
    
    qc_data <- getData(object = object, slot = "QC")
    subject_data <- getData(object = object, slot = "Subject")
    
    if(sum(is.na(qc_data)) +  sum(is.na(subject_data)) > 0){
      stop("Please impute MV first.\n")
    }
    
    if(is.null(qc_data)){
      if(method %in% c("svr", "loess")){
        stop("No QC samples in your data, please change other method.\n")
      }
    }
    
    ##sample-wise methods
    if(method %in% c("total", "median", "mean")){
      subject_data <- apply(subject_data, 2, function(x){
        x <- as.numeric(x)
        switch(method,
               total = {x/sum(x)},
               median = {x/median(x)},
               mean = {x/mean(x)}
        )
      })
      subject_data <- as.data.frame(subject_data)
      
      if(!is.null(qc_data)){
        qc_data <- apply(qc_data, 2, function(x){
          x <- as.numeric(x)
          switch(method,
                 total = {x/sum(x)},
                 median = {x/median(x)},
                 mean = {x/mean(x)}
          )
        })
        qc_data <- as.data.frame(qc_data)
      }
    }
    
    ##pqn (Probabilistic Quotient Normalization) method
    if(method == "pqn"){
      subject_data <- KODAMA::normalization(Xtrain = subject_data, 
                                            method = "pqn")$newXtrain
      if(!is.null(qc_data)){
        qc_data <- KODAMA::normalization(Xtrain = qc_data, 
                                         method = "pqn")$newXtrain
      }
    }
    
    
    if(method == "loess"){
      sample_info <- object@sample.info
      sample_info <- 
        sample_info %>% 
        filter(class %in% c('QC', 'Subject'))
      
      ms1_data <- object@ms1.data[[1]]
      ms1_data <- 
        ms1_data %>% 
        select(one_of(sample_info$sample.name))
      
      ###split data according to batch
      ##sample_info is a list
      sample_info <- 
        plyr::dlply(sample_info, .variables = .(batch))
      
      subject_data <- 
        lapply(sample_info, function(x){
          temp_subject_data <- 
            ms1_data %>% 
            select(one_of(x$sample.name[x$class == "Subject"]))
        })
      
      qc_data <- 
        lapply(sample_info, function(x){
          temp_subject_data <- 
            ms1_data %>% 
            select(one_of(x$sample.name[x$class == "QC"]))
        })
      
      subject_order <- 
        lapply(subject_data, function(x){
          object@sample.info$injection.order[match(colnames(x), object@sample.info$sample.name)]
        })
      
      qc_order <- 
        lapply(qc_data, function(x){
          object@sample.info$injection.order[match(colnames(x), 
                                                   object@sample.info$sample.name)]
        })
      
      ####begin data normalization
      qc_subject_data <- 
        vector(mode = "list", length = length(subject_data))
      
      for(batch_idx in 1:length(subject_data)){
        cat("Batch", batch_idx, "...", "\n")
        qc_subject_data[[batch_idx]] <-
          loessNor(
            subject_data = subject_data[[batch_idx]],
            qc_data = qc_data[[batch_idx]],
            subject_order = subject_order[[batch_idx]],
            qc_order = qc_order[[batch_idx]],
            optimization = TRUE,
            path = ".", 
            begin = begin, 
            end = end,
            step = step, 
            threads = threads
          )
        cat("\n")
      }
      
      qc_data <- 
        lapply(qc_subject_data, function(x){
          x[[1]]
        }) %>% 
        do.call(cbind, .)
      
      subject_data <- 
        lapply(qc_subject_data, function(x){
          x[[2]]
        }) %>% 
        do.call(cbind, .)
      
    }
    
    object@process.info$normalizeData <- list()
    object@process.info$normalizeData$method <- method
    object@process.info$normalizeData$keep.scale <- keep.scale
    object@process.info$normalizeData$begin <- begin
    object@process.info$normalizeData$end <- end
    object@process.info$normalizeData$step <- step
    
    sample_info <- object@sample.info
    subject_qc_data <- cbind(qc_data, subject_data)
    
    subject_qc_name <- dplyr::filter(.data = sample_info, class %in% c("Subject", "QC")) %>% 
      dplyr::pull(., sample.name)
    
    subject_qc_data <- subject_qc_data[, match(subject_qc_name, colnames(subject_qc_data))]
    ms1_data <- 
      object@ms1.data[[1]]
    ms1_data[,match(subject_qc_name, colnames(ms1_data))] <- subject_qc_data
    ms1_data <- list(ms1_data)
    object@ms1.data <- ms1_data
    invisible(object)
  }
)


# ##############svr normalization function
# SXTsvrNor <- function(sample,
#                       QC,
#                       tags,
#                       sample.order,
#                       QC.order,
#                       #used data
#                       multiple = 5,
#                       rerun = TRUE,
#                       peakplot = TRUE,
#                       path = NULL,
#                       datastyle = "tof",
#                       dimension1 = TRUE,
#                       threads = 1
#                       #parameters setting
# ) {
# 
#   options(warn = -1)
#   ######is there the e1071?
#   if (is.null(path)) {
#     path <- getwd()
#   } else{
#     dir.create(path)
#   }
# 
#   path1 <- file.path(path, "svr normalization result")
# 
#   dir.create(path1)
# 
#   if (!rerun) {
#     cat("Use previous normalization data\n")
#     # Sys.sleep(1)
#     load(file.path(path1, "normalization file"))
#   } else {
#     # library(snow)
#     # library(wordcloud)
# 
#     ichunks <- split((1:ncol(sample)), 1:threads)
#     svr.data <- BiocParallel::bplapply(ichunks,
#                                        FUN = svr.function,
#                                        BPPARAM = BiocParallel::SnowParam(workers = threads,
#                                                                          progressbar = TRUE),
#                                        sample = sample,
#                                        QC = QC,
#                                        sample.order = sample.order,
#                                        QC.order = QC.order,
#                                        multiple = multiple)
# 
#     sample.nor <- lapply(svr.data, function(x) {
#       x[[1]]
#     })
# 
#     QC.nor <- lapply(svr.data, function(x) {
#       x[[2]]
#     })
# 
#     index <- lapply(svr.data, function(x) {
#       x[[3]]
#     })
# 
# 
#     sample.nor <- do.call(cbind, sample.nor)
#     QC.nor <- do.call(cbind, QC.nor)
# 
#     index <- unlist(index)
# 
#     sample.nor <- sample.nor[,order(index)]
#     QC.nor <- QC.nor[,order(index)]
# 
#     QC.median <- apply(QC, 2, median)
#     if (dimension1) {
#       QC.nor <- t(t(QC.nor) * QC.median)
#       sample.nor <- t(t(sample.nor) * QC.median)
#     }
# 
#     # if (datastyle == "tof") {
#     #   colnames(QC.nor) <- colnames(sample.nor) <- tags["name", ]
#     # }
#     # if (datastyle == "mrm") {
#     #   colnames(QC.nor) <- colnames(sample.nor) <- tags["name", ]
#     # }
# 
#     save(QC.nor, sample.nor, file = file.path(path1, "normalization file"))
#   }
# 
#   rsd <- function(x) {
#     x <- sd(x) * 100 / mean(x)
#   }
# 
#   #following objects are the rsd of sample
#   #and QC before and after normalization
#   sample.rsd <- apply(sample, 2, rsd)
#   sample.nor.rsd <- apply(sample.nor, 2, rsd)
#   QC.rsd <- apply(QC, 2, rsd)
#   QC.nor.rsd <- apply(QC.nor, 2, rsd)
# 
# 
#   #sample.no.nor is the no normalization data added rsd information
#   #sample.svr is the normalization data added rsd information
# 
# 
#   sample.no.nor <- rbind(tags, sample.rsd, QC.rsd, sample, QC)
#   sample.svr <-
#     rbind(tags, sample.nor.rsd, QC.nor.rsd, sample.nor, QC.nor)
# 
#   save(sample.nor,
#        QC.nor,
#        tags,
#        sample.order,
#        QC.order,
#        file = file.path(path1, "data svr nor"))
#   write.csv(t(sample.svr), file.path(path1, "data svr nor.csv"))
# 
#   #generate all peaks plot
# 
#   if (peakplot) {
#     path2 <- file.path(path1, "peak plot")
#     dir.create(path2)
# 
#     cl <- snow::makeCluster(threads, type = "SOCK")
#     nc <- length(cl)
#     options(warn = -1)
#     ichunks <- split((1:ncol(sample)), 1:threads)
# 
#     if (datastyle == "tof")
#     {
#       snow::clusterApply(
#         cl,
#         x = ichunks,
#         fun = peakplot5,
#         sample = sample,
#         sample.nor = sample.nor,
#         QC = QC,
#         QC.nor = QC.nor,
#         sample.order = sample.order,
#         QC.order = QC.order,
#         tags = tags,
#         path = path2,
#         sample.rsd = sample.rsd,
#         QC.rsd = QC.rsd,
#         sample.nor.rsd = sample.nor.rsd,
#         QC.nor.rsd = QC.nor.rsd
#       )
#     }
#     else {
#       snow::clusterApply(
#         cl,
#         x = ichunks,
#         fun = peakplot6,
#         sample = sample,
#         sample.nor = sample.nor,
#         QC = QC,
#         QC.nor = QC.nor,
#         sample.order = sample.order,
#         QC.order = QC.order,
#         tags = tags,
#         path = path2,
#         sample.rsd = sample.rsd,
#         QC.rsd = QC.rsd,
#         sample.nor.rsd = sample.nor.rsd,
#         QC.nor.rsd = QC.nor.rsd
#       )
#     }
#   }
# 
# 
#   ##generate some statistics information
# 
#   compare.rsd(
#     sample.rsd = sample.rsd,
#     sample.nor.rsd = sample.nor.rsd,
#     QC.rsd = QC.rsd,
#     QC.nor.rsd =
#       QC.nor.rsd,
#     path = path1
#   )
#   options(warn = 0)
#   cat("SVR normalization is done\n")
# }



####LOESS normalization function
loessNor <- function(subject_data,
                     qc_data,
                     subject_order,
                     qc_order,
                     optimization = TRUE,
                     begin = 0.5,
                     end = 1,
                     step = 0.2,
                     path = ".",
                     threads = 4
) {
  cat("LOESS normalization...\n")
  
  temp.fun <- 
    function(idx,
             qc_data, 
             qc_order,
             subject_data,
             subject_order,
             optimization = TRUE,
             begin,
             end,
             step,
             cvMSE){
      if (optimization) {
        para <- cvMSE(
          unlist(qc_data[idx, ]),
          qc_order,
          begin1 = begin,
          end1 = end,
          step1 = step
        )
        
        loess.reg <-
          loess(unlist(qc_data[idx, ]) ~ qc_order,
                span = para[2],
                degree = para[1])
      }
      else {
        loess.reg <- loess(unlist(qc_data[idx, ]) ~ qc_order)
      }
      
      qc_data_pred <-
        summary(loess.reg)$fitted
      qc_nor1 <-
        unlist(qc_data[idx, ]) / qc_data_pred
      #if the predict value is 0, then set the ratio to 0
      qc_nor1[is.nan(unlist(qc_nor1))] <- 0
      qc_nor1[is.infinite(unlist(qc_nor1))] <- 0
      
      subject_data_pred <-
        predict(loess.reg, data.frame(qc_order = c(subject_order)))
      
      subject_nor1 <- unlist(subject_data[idx, ]) / subject_data_pred
      
      subject_nor1[is.nan(unlist(subject_nor1))] <-
        0
      subject_nor1[is.infinite(unlist(subject_nor1))] <-
        0
      subject_nor1[is.na(unlist(subject_nor1))] <-
        0
      subject_nor1[which(unlist(subject_nor1) < 0)] <-
        0
      
      return_result <- list(qc_nor1, subject_nor1)
      return(return_result)
    }
  
  peak_index <- 1:nrow(qc_data)
  
  data_nor <- 
    BiocParallel::bplapply(peak_index, 
                           FUN = temp.fun,
                           BPPARAM = BiocParallel::SnowParam(workers = threads,
                                                             progressbar = TRUE), 
                           qc_data = qc_data,
                           qc_order = qc_order,
                           subject_data = subject_data,
                           subject_order = subject_order,
                           optimization = optimization,
                           begin = begin,
                           end = end,
                           step = step,
                           cvMSE
    )
  
  qc_data_nor <- 
    lapply(data_nor, function(x){
      x[[1]]
    })
  
  qc_data_nor <- 
    do.call(rbind, qc_data_nor)
  
  subject_data_nor <- 
    lapply(data_nor, function(x){
      x[[2]]
    })
  subject_data_nor <- 
    do.call(rbind, subject_data_nor)
  
  qc_median <- apply(qc_data, 1, median)
  
  qc_data_nor <- qc_median * qc_data_nor
  subject_data_nor <- qc_median * subject_data_nor
  
  return_result <- list(qc_data_nor, subject_data_nor)
  return(return_result)
  cat("\n")
  cat("LOESS normalization is done\n")
}


#cvMSE is loess parameter optimization function
cvMSE <- function(qc, QC.order, begin1, end1, step1) {
  mse <- NULL
  nmse <- NULL
  cvmse <- NULL
  cvmse2 <- NULL
  
  para <- seq(begin1, end1, by = step1)
  for (i in 1:2) {
    for (j in para) {
      for (k in 2:(length(qc) - 1)) {
        loess.reg <- loess(qc[-k] ~ QC.order[-k], span = j, degree = i)
        predict.qc <- predict(loess.reg, QC.order[k])
        mse[k] <- (qc[k] - predict.qc) ^ 2
        nmse[k] <- (qc[k] - mean(qc)) ^ 2
      }
      cvmse1 <-
        rbind(j, mean(mse, na.rm = TRUE) / mean(nmse, na.rm = TRUE))
      cvmse2 <- cbind(cvmse2, cvmse1)
      mse <- NULL
      nmse <- NULL
    }
    
    cvmse3 <- rbind(i, cvmse2)
    cvmse <- cbind(cvmse, cvmse3)
    cvmse3 <- NULL
    cvmse2 <- NULL
  }
  return(cvmse[, which.min(cvmse[3,])])
}





