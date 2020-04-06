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
    # browser()
    method <- match.arg(method)
    
    if (class(object) != "metflowClass") {
      stop("Only the metflowClass is supported!\n")
    }
    
    ms1_data <- object@ms1.data
    
    if(length(ms1_data) > 1){
      stop("Please align batch first.\n")
    }
    
    ms1_data <- ms1_data[[1]]
    
    if(method == "svr" | method == "loess"){
      if(all(unique(object@sample.info$class) != "QC")){
        stop("No qc samples in your data, svr and loess is not available.\n")
      }
    }
    
    qc_data <- get_data(object = object, slot = "QC")
    subject_data <- get_data(object = object, slot = "Subject")
    
    if(sum(is.na(qc_data)) + sum(is.na(subject_data)) > 0){
      stop("Please impute MV first.\n")
    }
    
    if(is.null(qc_data)){
      if(method %in% c("svr", "loess")){
        stop("No qc samples in your data, please change other method.\n")
      }
    }
    
    ##sample-wise methods
    if(method %in% c("total", "median", "mean")){
      
      object@process.info$normalizeData <- list()
      object@process.info$normalizeData$method <- method
      object@process.info$normalizeData$keep.scale <- keep.scale
      
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
      
      object@process.info$normalizeData <- list()
      object@process.info$normalizeData$method <- method
      object@process.info$normalizeData$keep.scale <- keep.scale
      
      subject_data <- KODAMA::normalization(Xtrain = subject_data, 
                                            method = "pqn")$newXtrain
      if(!is.null(qc_data)){
        qc_data <- KODAMA::normalization(Xtrain = qc_data, 
                                         method = "pqn")$newXtrain
      }
    }
    
    
    if(method == "loess") {
      
      object@process.info$normalizeData <- list()
      object@process.info$normalizeData$method <- method
      object@process.info$normalizeData$keep.scale <- keep.scale
      object@process.info$normalizeData$begin <- begin
      object@process.info$normalizeData$end <- end
      object@process.info$normalizeData$step <- step
      
      sample_info <- object@sample.info
      sample_info <- 
        sample_info %>% 
        dplyr::filter(class %in% c('QC', 'Subject'))
      
      ms1_data <- object@ms1.data[[1]]
      
      ms1_data <- 
        ms1_data %>% 
        select(one_of(sample_info$sample.name))
      
      ###split data according to batch
      ##sample_info is a list
      sample_info <- 
        plyr::dlply(sample_info, .variables = plyr::.(batch))
      
      subject_data <- 
        lapply(sample_info, function(x){
          temp_subject_data <- 
            ms1_data %>% 
            dplyr::select(dplyr::one_of(x$sample.name[x$class == "Subject"]))
        })
      
      qc_data <- 
        lapply(sample_info, function(x){
          temp_subject_data <- 
            ms1_data %>% 
            dplyr::select(dplyr::one_of(x$sample.name[x$class == "QC"]))
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
        cat(crayon::yellow("Batch", batch_idx, "...", "\n"))
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
    
    
    if(method == "svr"){
      
      object@process.info$normalizeData <- list()
      object@process.info$normalizeData$method <- method
      
      sample_info <- object@sample.info
      sample_info <- 
        sample_info %>% 
        dplyr::filter(class %in% c('QC', 'Subject'))
      
      ms1_data <- object@ms1.data[[1]]
      
      ms1_data <- 
        ms1_data %>% 
        select(one_of(sample_info$sample.name))
      
      ###split data according to batch
      ##sample_info is a list
      sample_info <- 
        plyr::dlply(sample_info, .variables = plyr::.(batch))
      
      subject_data <- 
        lapply(sample_info, function(x){
          temp_subject_data <- 
            ms1_data %>% 
            dplyr::select(dplyr::one_of(x$sample.name[x$class == "Subject"]))
        })
      
      qc_data <- 
        lapply(sample_info, function(x){
          temp_subject_data <- 
            ms1_data %>% 
            dplyr::select(dplyr::one_of(x$sample.name[x$class == "QC"]))
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
        cat(crayon::yellow("Batch", batch_idx, "...", "\n"))
        
        qc_subject_data[[batch_idx]] <-
          svrNor(
            sample = subject_data[[batch_idx]],
            qc = qc_data[[batch_idx]],
            sample.order = subject_order[[batch_idx]],
            qc.order = qc_order[[batch_idx]],
            path = ".", 
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
    
    rownames(subject_data) <-
      rownames(qc_data) <-
      object@ms1.data[[1]]$name
  
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
svrNor <- function(sample,
                   qc,
                   sample.order,
                   qc.order,
                   #used data
                   multiple = 5,
                   path = ".",
                   dimension1 = TRUE,
                   threads = 3){
  options(warn = -1)
    sample <- 
      sample %>% 
      as.data.frame() %>% 
      t() %>% 
      as.data.frame()
    
    qc <-
      qc %>%
      as.data.frame() %>%
      t() %>%
      as.data.frame()
    
    ichunks <- split((1:ncol(sample)), 1:threads)
    
    svr.data <- BiocParallel::bplapply(
      ichunks,
      FUN = svr_function,
      BPPARAM = BiocParallel::SnowParam(workers = threads,
                                        progressbar = TRUE),
      sample = sample,
      qc = qc,
      sample.order = sample.order,
      qc.order = qc.order,
      multiple = multiple
    )
    
    sample.nor <- lapply(svr.data, function(x) {
      x[[1]]
    })

    qc.nor <- lapply(svr.data, function(x) {
      x[[2]]
    })

    index <- lapply(svr.data, function(x) {
      x[[3]]
    })


    sample.nor <- do.call(cbind, sample.nor)
    qc.nor <- do.call(cbind, qc.nor)

    index <- unlist(index)

    sample.nor <- sample.nor[,order(index)]
    qc.nor <- qc.nor[,order(index)]

    qc.median <- apply(qc, 2, median)
    
    if (dimension1) {
      qc.nor <- t(t(qc.nor) * qc.median)
      sample.nor <- t(t(sample.nor) * qc.median)
    }

  ##generate some statistics information
  cat(crayon::yellow("SVR normalization is done\n"))
  qc.nor <- 
    qc.nor %>% 
    t() %>% 
    as.data.frame()
  
  sample.nor <- 
    sample.nor %>% 
    t() %>% 
    as.data.frame()
  
  return_result <- list(qc.nor, sample.nor)
  return(return_result)
}

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
  cat(crayon::green("LOESS normalization...\n"))
  
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
  cat("\n")
  cat(crayon::green("LOESS normalization is done\n"))
  return(return_result)
}


#cvMSE is loess parameter optimization function
cvMSE <- function(qc, qc.order, begin1, end1, step1) {
  mse <- NULL
  nmse <- NULL
  cvmse <- NULL
  cvmse2 <- NULL
  
  para <- seq(begin1, end1, by = step1)
  for (i in 1:2) {
    for (j in para) {
      for (k in 2:(length(qc) - 1)) {
        loess.reg <- loess(qc[-k] ~ qc.order[-k], span = j, degree = i)
        predict.qc <- predict(loess.reg, qc.order[k])
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





#' @title svr_function
#' @description svr_function
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param index index
#' @param sample sample
#' @param qc qc
#' @param sample.order sample.order
#' @param qc.order qc.order
#' @param multiple multiple
#' @import crayon
#' @import tidyverse
#' @import BiocParallel
#' @import e1071
#' @importFrom magrittr %>%

setGeneric(
  name = "svr_function",
  def = function(index,
                 sample,
                 qc,
                 sample.order,
                 qc.order,
                 multiple) {
    # library(e1071)
    colnames(sample) <- colnames(qc)
    sample <- sample[, index, drop = FALSE]
    qc <- qc[, index, drop = FALSE]
    # cat("SVR normalization is finished: %\n")
    data.order <- c(sample.order, qc.order)
    
    data.nor <- lapply(c(1:ncol(sample)), function(i) {
      if (multiple != 1) {
        correlation <-
          abs(cor(x = rbind(sample, qc)[, i], y = rbind(sample, qc))[1, ])
        cor.peak <-
          match(names(sort(correlation, decreasing = TRUE)[1:6][-1]),
                names(correlation))
        rm(list = "correlation")
        svr.reg <- e1071::svm(qc[, cor.peak], qc[, i])
      } else{
        svr.reg <- e1071::svm(unlist(qc[, i]) ~ qc.order)
      }
      
      predict.qc <- summary(svr.reg)$fitted
      qc.nor1 <- qc[, i] / predict.qc
      
      #if the predict value is 0, then set the ratio to 0
      qc.nor1[is.nan(unlist(qc.nor1))] <- 0
      qc.nor1[is.infinite(unlist(qc.nor1))] <- 0
      qc.nor1[is.na(unlist(qc.nor1))] <- 0
      qc.nor1[which(unlist(qc.nor1) < 0)] <- 0
      
      if (multiple != 1) {
        predict.sample <- predict(svr.reg, sample[, cor.peak])
      } else{
        predict.sample <-
          predict(svr.reg, data.frame(qc.order = c(sample.order)))
      }
      
      sample.nor1 <- sample[, i] / predict.sample
      sample.nor1[is.nan(unlist(sample.nor1))] <- 0
      sample.nor1[is.infinite(unlist(sample.nor1))] <- 0
      sample.nor1[is.na(unlist(sample.nor1))] <- 0
      sample.nor1[which(unlist(sample.nor1) < 0)] <- 0
      
      return(list(sample.nor1, qc.nor1))
      
    })
    
    sample.nor <- lapply(data.nor, function(x)
      x[[1]])
    qc.nor <- lapply(data.nor, function(x)
      x[[2]])
    rm(list = "data.nor")
    sample.nor <- t(do.call(rbind, sample.nor))
    qc.nor <- t(do.call(rbind, qc.nor))
    
    colnames(sample.nor) <-
      colnames(qc.nor) <- colnames(sample)
    rm(list = c("sample", "qc"))
    
    svr.data <-
      list(sample.nor = sample.nor,
           qc.nor = qc.nor,
           index = index)
    rm(list = c("sample.nor", "qc.nor"))
    return(svr.data)
    
  }
)
