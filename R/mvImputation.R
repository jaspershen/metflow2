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

impute_mv = function(object,
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
                     ...) {
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
  
  qc_data <- get_data(object = object, slot = "QC")
  
  subject_data <-
    get_data(object = object, slot = "Subject")
  
  ##remove the peaks are NAs in all qc or objects
  if (!is.null(qc_data)) {
    idx1 <- apply(qc_data, 1, function(x) {
      sum(is.na(x)) / ncol(qc_data)
    }) %>%
      `==`(1) %>%
      which()
  } else{
    idx1 <- NULL
  }
  
  if (!is.null(subject_data)) {
    idx2 <- apply(subject_data, 1, function(x) {
      sum(is.na(x)) / ncol(subject_data)
    }) %>%
      `==`(1) %>%
      which()
  } else{
    idx2 <- NULL
  }
  
  remove_idx <-
    intersect(idx1, idx2)
  
  if (length(remove_idx) > 0) {
    qc_data <- qc_data[-remove_idx, , drop = FALSE]
    subject_data <- subject_data[-remove_idx, , drop = FALSE]
    object@ms1.data[[1]] <-
      object@ms1.data[[1]][-remove_idx, , drop = FALSE]
    ms1_data <- object@ms1.data[[1]]
  }
  
  if (is.null(qc_data)) {
    subject_qc_data <- subject_data
  }
  
  if (is.null(subject_data)) {
    subject_qc_data <- qc_data
  }
  
  if (!is.null(qc_data) & !is.null(subject_data)) {
    subject_qc_data <- cbind(qc_data, subject_data)
  }
  
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




#' @title sxtMVimputation
#' @description sxtMVimputation
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param data data
#' @param method method
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
#' @return result

sxtMVimputation = function(data,
                           method = c("knn",
                                      "rf",
                                      "mean",
                                      "median",
                                      "zero",
                                      "minimum",
                                      "bpca",
                                      "svdImpute",
                                      "ppca"),
                           ## knn parameters
                           k = 10,
                           rowmax = 0.5,
                           colmax = 0.8,
                           maxp = 1500,
                           rng.seed = 362436069,
                           ## missForest parameters
                           maxiter = 10,
                           ntree = 100,
                           decreasing = FALSE,
                           nPcs = 2,
                           maxSteps = 100,
                           threshold = 1e-04,
                           ...) {
  method = match.arg(method)
  ## KNN method
  if (method == "knn") {
    # library(impute)
    if (exists(".Random.seed"))
      rm(.Random.seed)
    data.knn <- impute::impute.knn(
      as.matrix(t(data)),
      k = k,
      rowmax = rowmax,
      colmax = colmax,
      maxp = maxp,
      rng.seed = rng.seed
    )
    data.knn <- tibble::as_tibble(t(data.knn[["data"]]))
    return(data.knn)
  }
  
  if (method == "rf") {
    data.rf <- missForest::missForest(
      t(data),
      maxiter = maxiter,
      ntree = ntree,
      decreasing = decreasing,
      ...
    )
    data.rf <- tibble::as_tibble(t(data.rf$ximp))
    return(data.rf)
  }
  
  ## mean imputation
  if (method %in% c("mean", "median", "zero", "minimum")) {
    data.result <- apply(data, 1, function(x) {
      x <- as.numeric(x)
      if (all(!is.na(x))) {
        return(x)
      } else{
        x[is.na(x)] <-
          switch(
            EXPR = method,
            mean = mean(x, na.rm = TRUE),
            median = median(x, na.rm = TRUE),
            zero = 0,
            minimum = min(x, na.rm = TRUE)
          )
        x
      }
    })
    rownames(data.result) <- colnames(data)
    return(tibble::as_tibble(t(data.result)))
  }
  
  ##BPCA, SVD and ppca
  if (method %in% c("bpca", "svdImpute", "ppca")) {
    data.result <- pcaMethods::pca(
      t(data),
      method = method,
      nPcs = nPcs,
      maxSteps = maxSteps,
      threshold = threshold
    )
    data.result <- t(pcaMethods::completeObs(data.bpca))
    return(data.result)
  }
}
