#' @title align_batch
#' @description Align different batch peaks tables.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param object A metflowClass object.
#' @param combine.mz.tol m/z tolerance for batch alignment, default is 25 ppm.
#' @param combine.rt.tol RT tolerance for batch alignment, default is 30 seconds.
#' @param use.int.tol Whether use intensity match for batch aglignment.
#' @return A new metflowClass object.
#' @export

align_batch = function(
  object,
  combine.mz.tol = 25,
  combine.rt.tol = 30,
  use.int.tol = FALSE
){
  if (class(object) != "metflowClass") {
    stop("Only for metflowClass object\n")
  }
  
  ms1_data <- object@ms1.data
  if (length(ms1_data) == 1) {
    return(object)
  }
  
  cat("Rough aligning...\n")
  roughMatchResult <- roughAlign(
    peak.table = ms1_data,
    combine.mz.tol = combine.mz.tol,
    combine.rt.tol = combine.rt.tol
  )
  
  cat("Accurate aligning...\n")
  accurateMatchResult <-
    accurateAlign(
      peak.table = ms1_data,
      simple.data = roughMatchResult,
      use.int.tol = use.int.tol
    )
  object@ms1.data <- list(accurateMatchResult)
  
  object@process.info$alignBatch <- list()
  object@process.info$alignBatch$combine.mz.tol <-
    combine.mz.tol
  object@process.info$alignBatch$combine.rt.tol <-
    combine.rt.tol
  
  invisible(object)
}


#' @title alignBatch
#' @description Align different batch peaks tables.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param object A metflowClass object.
#' @param combine.mz.tol m/z tolerance for batch alignment, default is 25 ppm.
#' @param combine.rt.tol RT tolerance for batch alignment, default is 30 seconds.
#' @param use.int.tol Whether use intensity match for batch aglignment.
#' @return A new metflowClass object.
#' @export

alignBatch = function(
  object,
  combine.mz.tol = 25,
  combine.rt.tol = 30,
  use.int.tol = FALSE
){
  
  if(!silence.deprecated){
    cat(crayon::yellow("`alignBatch()` is deprecated, please use `align_batch()`"))
  }
  
  if (class(object) != "metflowClass") {
    stop("Only for metflowClass object\n")
  }
  
  ms1_data <- object@ms1.data
  
  if (length(ms1_data) == 1) {
    return(object)
  }
  
  cat("Rough aligning...\n")
  roughMatchResult <- roughAlign(
    peak.table = ms1_data,
    combine.mz.tol = combine.mz.tol,
    combine.rt.tol = combine.rt.tol
  )
  
  cat("Accurate aligning...\n")
  accurateMatchResult <-
    accurateAlign(
      peak.table = ms1_data,
      simple.data = roughMatchResult,
      use.int.tol = use.int.tol
    )
  object@ms1.data <- list(accurateMatchResult)
  
  object@process.info$alignBatch <- list()
  object@process.info$alignBatch$combine.mz.tol <-
    combine.mz.tol
  object@process.info$alignBatch$combine.rt.tol <-
    combine.rt.tol
  
  invisible(object)
}




#' @title roughAlign
#' @description roughAlign
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param peak.table peak.table
#' @param combine.mz.tol combine.mz.tol
#' @param combine.rt.tol combine.rt.tol
#' @return result

roughAlign <- function(peak.table,
                       combine.mz.tol = 25,
                       combine.rt.tol = 30) {
  if (length(peak.table) == 1)
    return(peak.table[[1]])
  
  batch1 <- peak.table[[1]]
  batch2 <- peak.table[[2]]
  
  ###generate sinplifued datasets
  simple.batch1 <- simplyData(
    data = batch1,
    combine.mz.tol = combine.mz.tol,
    combine.rt.tol = combine.rt.tol
  )
  
  simple.batch2 <- simplyData(
    data = batch2,
    combine.mz.tol = combine.mz.tol,
    combine.rt.tol = combine.rt.tol
  )
  
  rm(list = c("batch1", "batch2"))
  
  ### rough matching
  data1 <- simple.batch1[, c(2:3)]
  data2 <- simple.batch2[, c(2:3)]
  
  match.result <- SXTMTmatch2(
    data1 = as.matrix(data1),
    data2 = as.matrix(data2),
    mz.tol = combine.mz.tol,
    rt.tol = combine.rt.tol,
    rt.error.type = "abs"
  )
  
  rm(list = c("data1", "data2"))
  
  simple.batch1 <- simple.batch1[match.result[, 1],]
  simple.batch2 <- simple.batch2[match.result[, 2],]
  
  rm("match.result")
  
  
  simple.data <- list(simple.batch1, simple.batch2)
  return(simple.data)
}


#' @title simplyData
#' @description simplyData
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param data data
#' @param combine.mz.tol combine.mz.tol
#' @param combine.rt.tol combine.rt.tol
#' @return result

simplyData <- function(data,
                       combine.mz.tol = 5,
                       combine.rt.tol = 30) {
  data <- data[order(data$mz),]
  name <- data$name
  mz <- data$mz
  rt <- data$rt
  int <- apply(data[, -c(1:3)], 1, function(x)
    mean(x, na.rm = TRUE))
  names(mz) <- names(rt) <- names(int) <- name
  
  #group peaks according to mz and RT
  
  all.index <- 1:length(name)
  remain.idx <- vector(mode = "list", length = length(name))
  for (idx in 1:length(name)) {
    if (all(idx != all.index))
      next()
    temp.mz <- mz[idx]
    temp.rt <- rt[idx]
    mz.left <- temp.mz - combine.mz.tol * temp.mz / 10 ^ 6
    mz.right <- temp.mz + combine.mz.tol * temp.mz / 10 ^ 6
    rt.left <- temp.rt - combine.rt.tol
    rt.right <- temp.rt + combine.rt.tol
    
    temp.idx <- which(mz > mz.left & mz <=  mz.right &
                        rt > rt.left & rt <= rt.right)
    if (length(temp.idx) == 1) {
      remain.idx[[idx]] <- temp.idx
    } else{
      max.idx <- temp.idx[which.max(int[temp.idx])]
      remain.idx[[idx]] <- max.idx
      all.index <- setdiff(all.index, setdiff(temp.idx, max.idx))
    }
    
  }
  
  
  remain.idx <- unique(unlist(remain.idx))
  
  simple.data <- data[remain.idx,]
  rm(list = c("data", "mz", "rt", "name", "int"))
  return(simple.data)
}


#' @title baMZplot
#' @description baMZplot
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param simple.data data
#' @return result

baMZplot <- function(simple.data) {
  my.theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18)
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 15),
      axis.text.y = ggplot2::element_text(size = 15)
    ) +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 12)) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 10))
  
  mz.error <-
    (simple.data[[2]]$mz - simple.data[[1]]$mz) * 10 ^ 6 / simple.data[[2]]$mz
  mz1 <- simple.data[[1]]$mz
  mz.error.sd <- sd(abs(mz.error))
  mz.error.sd <-
    paste("m/z error standard:", round(mz.error.sd, 2), "ppm")
  temp.data <- data.frame(mz1, mz.error, stringsAsFactors = FALSE)
  mz.plot <- ggplot2::ggplot(data = temp.data,
                             ggplot2::aes(x = mz1, y = mz.error)) +
    ggplot2::geom_point() +
    my.theme +
    ggplot2::labs(x = "m/z (Batch 1)",
                  y = "m/z deviation (ppm)") +
    ggplot2::ggtitle(paste("m/z vs. m/z deviation; ", mz.error.sd)) +
    ggplot2::geom_hline(aes(yintercept = 0)) +
    ggplot2::annotate(
      geom = "text",
      x = -Inf,
      y = Inf,
      hjust = -0.2,
      vjust = 2,
      label = mz.error.sd
    )
  return(mz.plot)
}


#' @title baRTplot
#' @description baRTplot
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param simple.data simple.data
#' @return result
baRTplot <- function(simple.data) {
  my.theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18)
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 15),
      axis.text.y = ggplot2::element_text(size = 15)
    ) +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 12)) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 10))
  
  rt.error <- simple.data[[2]]$rt - simple.data[[1]]$rt
  rt1 <- simple.data[[1]]$rt
  temp.data <- data.frame(rt1, rt.error, stringsAsFactors = FALSE)
  rt.error.sd <- sd(abs(rt.error))
  rt.error.sd <-
    paste("RT error standard:", round(rt.error.sd, 2), "second")
  rt.plot <- ggplot2::ggplot(data = temp.data,
                             ggplot2::aes(x = rt1, y = rt.error)) +
    ggplot2::geom_point() +
    my.theme +
    ggplot2::labs(x = "Retention time (s, Batch 1)",
                  y = "RT deviation (second)") +
    ggplot2::ggtitle(paste("RT vs. RT deviation;", rt.error.sd)) +
    ggplot2::geom_hline(aes(yintercept = 0)) +
    ggplot2::annotate(
      geom = "text",
      x = -Inf,
      y = Inf,
      hjust = -0.2,
      vjust = 2,
      label = rt.error.sd
    )
  
  return(rt.plot)
}


#' @title baINTplot
#' @description baINTplot
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param simple.data simple.data
#' @return result

baINTplot <- function(simple.data) {
  my.theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18)
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 15),
      axis.text.y = ggplot2::element_text(size = 15)
    ) +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 12)) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 10))
  
  int1 <-
    log(apply(simple.data[[1]][, -c(1:3)], 1, function(x)
      mean(x, na.rm = TRUE)) + 1, 10)
  int2 <-
    log(apply(simple.data[[2]][, -c(1:3)], 1, function(x)
      mean(x, na.rm = TRUE)) + 1, 10)
  int.error <- int2 - int1
  
  int.error.sd <- sd(abs(int.error))
  int.error.sd <-
    paste("Log10int error standard:", round(int.error.sd, 2))
  temp.data <- data.frame(int1, int.error, stringsAsFactors = FALSE)
  int.plot <- ggplot2::ggplot(data = temp.data,
                              ggplot2::aes(x = int1, y = int.error)) +
    ggplot2::geom_point() +
    my.theme +
    ggplot2::labs(x = "Log10intensity (mean, Batch 1)",
                  y = "Log10intensity deviation (second)") +
    ggplot2::ggtitle(paste("Log10int vs. Log10int deviation;", int.error.sd)) +
    ggplot2::geom_hline(aes(yintercept = 0)) +
    ggplot2::annotate(
      geom = "text",
      x = -Inf,
      y = Inf,
      hjust = -0.2,
      vjust = 2,
      label = int.error.sd
    )
  return(int.plot)
}


#' @title peak.table
#' @description peak.table
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param peak.table peak.table
#' @param simple.data simple.data
#' @param use.int.tol use.int.tol
#' @return result

accurateAlign <- function(peak.table,
                          simple.data,
                          use.int.tol) {
  ##retrieve mz ,RT and int sd
  if (length(peak.table) == 1)
    return(peak.table[[1]])
  mz.error <-
    (simple.data[[2]]$mz - simple.data[[1]]$mz) * 10 ^ 6 / simple.data[[2]]$mz
  mz.error.sd <- sd(abs(mz.error))
  
  rt.error <- simple.data[[2]]$rt - simple.data[[1]]$rt
  rt.error.sd <- sd(abs(rt.error))
  
  int1 <-
    log(apply(simple.data[[1]][, -c(1:3)], 1, function(x)
      mean(x, na.rm = TRUE)) + 1, 10)
  int2 <-
    log(apply(simple.data[[2]][, -c(1:3)], 1, function(x)
      mean(x, na.rm = TRUE)) + 1, 10)
  int.error <- int2 - int1
  
  int.error.sd <- sd(abs(int.error))
  
  if (!use.int.tol) {
    int.error.sd <- int.error.sd * 1000000
  }
  
  
  ###begin alignment
  ref.batch <- peak.table[[1]]
  for (i in 2:length(peak.table)) {
    cor.batch <- peak.table[[i]]
    new.batch <- align2Batch(
      batch1 = ref.batch,
      batch2 = cor.batch,
      mz.error.sd = mz.error.sd,
      rt.error.sd = rt.error.sd,
      int.error.sd = int.error.sd,
      fold = 4,
      mz.weight = 0.4,
      rt.weight = 0.4,
      int.weight = 0.2
    )
    ref.batch <- new.batch
    rm("new.batch")
  }
  
  rm(list = c("peak.table"))
  colnames(ref.batch)[1:3] <- c("name", "mz", "rt")
  return(ref.batch)
}


#' @title align2Batch
#' @description peak.table
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param batch1 batch1
#' @param batch2 batch2
#' @param mz.error.sd mz.error.sd
#' @param rt.error.sd rt.error.sd
#' @param int.error.sd int.error.sd
#' @param fold fold
#' @param mz.weight mz.weight
#' @param rt.weight rt.weight
#' @param int.weight int.weight
#' @return result

align2Batch <- function(batch1,
                        batch2,
                        mz.error.sd,
                        rt.error.sd,
                        int.error.sd,
                        fold = 4,
                        mz.weight = 0.4,
                        rt.weight = 0.4,
                        int.weight = 0.2) {
  data1 <- batch1[, c("mz", "rt")]
  int1 <-
    log(apply(batch1[, -c(1:3)], 1, function(x)
      mean(x, na.rm = TRUE)) + 1, 10)
  data1 <- data.frame(data1, int1, stringsAsFactors = FALSE)
  colnames(data1)[3] <- "int"
  rownames(data1) <- batch1$name
  
  data2 <- batch2[, c("mz", "rt")]
  int2 <-
    log(apply(batch2[, -c(1:3)], 1, function(x)
      mean(x, na.rm = TRUE)) + 1, 10)
  data2 <- data.frame(data2, int2, stringsAsFactors = FALSE)
  colnames(data2)[3] <- "int"
  rownames(data2) <- batch2$name
  
  ##batch1 match batch2
  match.result1 <- MRImatch(
    data1 = data1,
    data2 = data2,
    mz.tol = mz.error.sd * fold,
    rt.tol = rt.error.sd * fold,
    rt.error.type = "abs",
    int.tol = int.error.sd * fold
  )
  unique.index1 <- unique(match.result1[, 1])
  
  remain.idx <- unlist(lapply(unique.index1, function(x) {
    temp.idx <- which(x == match.result1[, 1])
    if (length(temp.idx) == 1)
      return(temp.idx)
    temp.result <- match.result1[temp.idx,]
    ##score
    temp.mz.error <- temp.result[, "mz.error"]
    temp.rt.error <- temp.result[, "rt.error"]
    temp.int.error <- temp.result[, "int.error"]
    temp.mz.score <- sapply(temp.mz.error, function(x) {
      matchScore(error = x, sd = mz.error.sd)
    })
    
    temp.rt.score <- sapply(temp.rt.error, function(x) {
      matchScore(error = x, sd = rt.error.sd)
    })
    
    temp.int.score <- sapply(temp.int.error, function(x) {
      matchScore(error = x, sd = int.error.sd)
    })
    
    temp.score <-
      mz.weight * temp.mz.score + rt.weight * temp.rt.score + int.weight * temp.int.score
    return(temp.idx[which.max(temp.score)])
  }))
  
  
  match.result1 <- match.result1[remain.idx,]
  
  
  ##batch2 match batch1
  match.result2 <- MRImatch(
    data1 = data2,
    data2 = data1,
    mz.tol = mz.error.sd * fold,
    rt.tol = rt.error.sd * fold,
    rt.error.type = "abs",
    int.tol = int.error.sd * fold
  )
  unique.index1 <- unique(match.result2[, 1])
  
  remain.idx <- unlist(lapply(unique.index1, function(x) {
    temp.idx <- which(x == match.result2[, 1])
    if (length(temp.idx) == 1)
      return(temp.idx)
    temp.result <- match.result2[temp.idx,]
    ##score
    temp.mz.error <- temp.result[, "mz.error"]
    temp.rt.error <- temp.result[, "rt.error"]
    temp.int.error <- temp.result[, "int.error"]
    temp.mz.score <- sapply(temp.mz.error, function(x) {
      matchScore(error = x, sd = mz.error.sd)
    })
    
    temp.rt.score <- sapply(temp.rt.error, function(x) {
      matchScore(error = x, sd = rt.error.sd)
    })
    
    temp.int.score <- sapply(temp.int.error, function(x) {
      matchScore(error = x, sd = int.error.sd)
    })
    
    temp.score <-
      mz.weight * temp.mz.score + rt.weight * temp.rt.score + int.weight * temp.int.score
    return(temp.idx[which.max(temp.score)])
  }))
  
  
  match.result2 <- match.result2[remain.idx,]
  
  name1 <- paste(match.result1[, 1], match.result1[, 2], sep = "_")
  name2 <- paste(match.result2[, 2], match.result2[, 1], sep = "_")
  
  intersect.name <- intersect(name1, name2)
  
  temp.index <- which(name1 %in% intersect.name)
  
  match.result <- match.result1[temp.index,]
  
  ##combine two batch data
  batch1 <- batch1[match.result[, 1],]
  batch2 <- batch2[match.result[, 2],]
  
  ##new name, mz and RT
  new.mz <-
    data.frame(batch1[, 2], batch2[, 2], stringsAsFactors = FALSE)
  new.mz <- apply(new.mz, 1, mean)
  
  new.rt <-
    data.frame(batch1[, 3], batch2[, 3], stringsAsFactors = FALSE)
  new.rt <- apply(new.rt, 1, mean)
  
  new.name <-
    paste("M", round(new.mz), "T", round(new.rt), sep = "")
  new.name <- reName(new.name)
  
  
  return.data <- data.frame(new.name, new.mz, new.rt,
                            batch1[, -c(1:3)], batch2[, -c(1:3)], stringsAsFactors = FALSE)
  rownames(return.data) <- new.name
  rm(
    list = c(
      "batch1",
      "batch2",
      "new.mz",
      "new.rt",
      "new.name",
      "match.result",
      "match.result1",
      "match.result2"
    )
  )
  return(return.data)
}



#' @title matchScore
#' @description matchScore
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param error error
#' @param sd sd
#' @return result
matchScore <- function(error, sd) {
  (sd / error) ^ 2
}


#' @title reName
#' @description reName
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param name name
#' @return result

reName <- function(name) {
  temp.name <- unique(name)
  lapply(temp.name, function(x) {
    temp.idx <- which(x == name)
    if (length(temp.idx) > 1) {
      paste(name[temp.idx], 1:length(temp.idx), sep = "_")
    }
  })
}


# reName <- function(name) {
#   temp.name <- unique(name)
#   for (i in 1:length(temp.name)) {
#     temp.idx <- which(temp.name[i] == name)
#     if (length(temp.idx) > 1) {
#       name[temp.idx] <-
#         paste(name[temp.idx], 1:length(temp.idx), sep = "_")
#     }
#   }
#   name
# }


#' @title getBatchAlignmentInfo
#' @description getBatchAlignmentInfo
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param raw.data raw.data
#' @param rough.align.data rough.align.data
#' @param accurate.align.data accurate.align.data
#' @return result


getBatchAlignmentInfo <- function(raw.data,
                                  rough.align.data,
                                  accurate.align.data) {
  if (length(raw.data) == 1)
    return("Only one batch data.")
  mz.error <-
    (rough.align.data[[2]]$mz - rough.align.data[[1]]$mz) * 10 ^ 6 / rough.align.data[[2]]$mz
  mz.error.sd <- sd(abs(mz.error))
  
  rt.error <- rough.align.data[[2]]$rt - rough.align.data[[1]]$rt
  rt.error.sd <- sd(abs(rt.error))
  
  int1 <-
    log(apply(rough.align.data[[1]][, -c(1:3)], 1, function(x)
      mean(x, na.rm = TRUE)) + 1, 10)
  int2 <-
    log(apply(rough.align.data[[2]][, -c(1:3)], 1, function(x)
      mean(x, na.rm = TRUE)) + 1, 10)
  int.error <- int2 - int1
  int.error.sd <- sd(abs(int.error))
  
  parameter.info <-
    paste(
      "The standard for m/z, RT and log10intensity errors are",
      round(mz.error.sd, 2),
      ",",
      round(rt.error.sd, 2),
      ",",
      'and',
      round(int.error.sd, 2),
      ",respectively. So the tolerances for m/z, RT and log10intensity are",
      4 * round(mz.error.sd, 2),
      ",",
      4 * round(rt.error.sd, 2),
      ",",
      'and',
      4 * round(int.error.sd, 2),
      ",respectively."
    )
  
  peak.info <-
    paste(
      "There are ",
      length(raw.data),
      " batches. And the peak number are ",
      paste(unlist(lapply(raw.data, nrow)), collapse = ","),
      " respectively. ",
      "After batch alignment, the peak number is ",
      nrow(accurate.align.data),
      ".",
      sep = ""
    )
  
  return(paste(parameter.info, peak.info))
}



#' @title MRImatch
#' @description MRImatch
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param data1 data1
#' @param data2 data2
#' @param mz.tol mz.tol
#' @param rt.tol rt.tol
#' @param rt.error.type rt.error.type
#' @param int.tol int.tol
#' @return result

MRImatch = function(data1,
                    data2,
                    mz.tol,
                    #rt.tol is relative
                    rt.tol = 30,
                    rt.error.type = c("relative", "abs"),
                    int.tol = 1){
  rt.error.type <- match.arg(rt.error.type)
  #
  if (nrow(data1) == 0 | nrow(data2) == 0) {
    result <- NULL
    return(result)
  }
  # mz1 <- as.numeric(data1[, 1])
  # rt1 <- as.numeric(data1[, 2])
  info1 <- data1[, c(1, 2, 3)]
  info1 <- apply(info1, 1, list)
  
  mz2 <- as.numeric(data2[, 1])
  rt2 <- as.numeric(data2[, 2])
  int2 <- as.numeric(data2[, 3])
  
  result <- pbapply::pblapply(info1, function(x) {
    temp.mz1 <- x[[1]][[1]]
    temp.rt1 <- x[[1]][[2]]
    temp.int1 <- x[[1]][[3]]
    mz.error <- abs(temp.mz1 - mz2) * 10 ^ 6 / temp.mz1
    if (rt.error.type == "relative") {
      rt.error <- abs(temp.rt1 - rt2) * 100 / temp.rt1
    } else{
      rt.error <- abs(temp.rt1 - rt2)
    }
    
    int.error <- abs(temp.int1 - int2)
    
    j <-
      which(mz.error <= mz.tol &
              rt.error <= rt.tol & int.error <= int.tol)
    if (length(j) == 0) {
      matrix(NA, ncol = 10)
    } else{
      cbind(
        j,
        temp.mz1,
        mz2[j],
        mz.error[j],
        temp.rt1,
        rt2[j],
        rt.error[j],
        temp.int1,
        int2[j],
        int.error[j]
      )
    }
  })
  
  if (length(result) == 1) {
    result <- cbind(1, result[[1]])
  } else{
    result <- mapply(function(x, y) {
      list(cbind(x, y))
    },
    x <- 1:length(info1),
    y = result)
    result <- do.call(rbind, result)
  }
  
  result <-
    matrix(result[which(!apply(result, 1, function(x)
      any(is.na(x)))),], ncol = 11)
  if (nrow(result) == 0)
    return(NULL)
  colnames(result) <-
    c(
      "Index1",
      "Index2",
      "mz1",
      "mz2",
      "mz.error",
      "rt1",
      "rt2",
      "rt.error",
      "int1",
      "int2",
      "int.error"
    )
  result <- result 
}



#' @title SXTMTmatch2
#' @description SXTMTmatch2
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param data1 data1
#' @param data2 data2
#' @param mz.tol mz.tol
#' @param rt.tol rt.tol
#' @param rt.error.type rt.error.type
#' @return result

SXTMTmatch2 = function(data1,
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
  info1 <- data1[, c(1, 2)]
  info1 <- apply(info1, 1, list)
  
  mz2 <- as.numeric(data2[, 1])
  rt2 <- as.numeric(data2[, 2])
  
  result <- pbapply::pblapply(info1, function(x) {
    temp.mz1 <- x[[1]][[1]]
    temp.rt1 <- x[[1]][[2]]
    mz.error <- abs(temp.mz1 - mz2) * 10 ^ 6 / temp.mz1
    if (rt.error.type == "relative") {
      rt.error <- abs(temp.rt1 - rt2) * 100 / temp.rt1
    } else{
      rt.error <- abs(temp.rt1 - rt2)
    }
    
    j <- which(mz.error <= mz.tol & rt.error <= rt.tol)
    if (length(j) == 0) {
      matrix(NA, ncol = 7)
    } else{
      cbind(j, temp.mz1, mz2[j], mz.error[j], temp.rt1, rt2[j], rt.error[j])
    }
  })
  
  if (length(result) == 1) {
    result <- cbind(1, result[[1]])
  } else{
    result <- mapply(function(x, y) {
      list(cbind(x, y))
    },
    x <- 1:length(info1),
    y = result)
    result <- do.call(rbind, result)
  }
  
  result <-
    matrix(result[which(!apply(result, 1, function(x)
      any(is.na(x)))), ], ncol = 8)
  if (nrow(result) == 0)
    return(NULL)
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
}

