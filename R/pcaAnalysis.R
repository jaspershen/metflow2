#' @title PCAanalysis
#' @description PCA analysis for MetFlowData.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object A metflowClass object.
#' @param scale.method Scale method.
#' @param slot Class of data.
#' @return A ggplot object.

doPCA <- function(object,
                  scale.method = c("no", "auto", "pareto", "center"),
                  slot = c("QC", "Subject")) {
  
  # requireNamespace("tidyverse")
  # requireNamespace("dplyr")
  if(class(object) != "metflowClass"){
    stop("Only the metflowClass is supported!\n")
  }
  
  if(length(object@ms1.data) != 1){
    stop("Please align batches first.\n")
  }
  
  if(any(!slot %in% c("QC", "Subject"))){
    stop("Slot can only be QC and Subject.\n")
  }
  
  tag <- getData(object = object, slot = "Tags")
  name <- dplyr::pull(.data = tag, name)
  
  data <- lapply(slot, function(x){
    x <- getData(object = object, slot = x)
    if(is.null(x)){
      return(x)
    }
    rownames(x) <- name
    as.data.frame(t(x))
  })
  
  remain_idx <- which(!unlist(lapply(data, is.null)))
  slot <- slot[remain_idx]
  data <- data[remain_idx]
  
  class <- mapply(function(x, y){
  rep(y, nrow(x))
  },
  x = data,
  y = slot)
  
  class <- unlist(class)
  
  data <- do.call(rbind, data)
  if(sum(is.na(data)) != 0){
    stop("Please impute MV first.\n")
  }
  
  data <- sxtScale(df = data, method = scale.method)
  
  data <- data.frame(data, class, stringsAsFactors = FALSE)
  data <- data.frame(data, name = rownames(data), 
                     stringsAsFactors = FALSE)
  
  pca_object <- prcomp(select(data, -one_of(c('name', 'class'))))
  
  plot <- ggplot2::autoplot(pca_object, data = data, colour = 'class',
                   frame = TRUE, frame.type = "norm") +
    # scale_y_continuous(limits = c(0, 0.015)) +
    # scale_x_continuous(limits = c(-0.1, 0.1)) +
    theme_bw() +
    scale_colour_manual(values = c("#E64B35FF", "#4DBBD5FF")) +
    theme(axis.title = element_text(size = 15),
          axis.text = element_text(size = 12), 
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          strip.background = element_rect(fill = "#0099B47F"),
          strip.text = element_text(color = "white", size = 15))
  plot
}


SXTpca <- function(subject = NULL,
                   qc = NULL,
                   info = NULL,
                   #used data
                   QC = FALSE,
                   scale.method = "auto",
                   path = ".") {
  ##remove the peaks with zero ratio > 50
  temp.idx1 <- apply(subject, 1, function(x) sum(x == 0)/ncol(subject)) > 0.5
  temp.idx1 <- which(temp.idx1)
  if(length(temp.idx1) > 0){
    subject <- subject[-temp.idx1, ]
    if(!is.null(qc)) {
      qc <- qc[-temp.idx1,]
    }
  }

  if (path != ".") {
    dir.create(path)
  }
  if (any(is.na(subject)) |
      any(is.na(qc)))
    stop("Please impute MV in subject or QC samples.")
  if (is.null(subject))
    stop("Subject sample is NULL")
  if (!is.null(qc)) {
    if (nrow(subject) != nrow(qc))
      stop("ThSe row number of Subject and QC must same")
  }
  if (is.null(qc) & QC)
    stop("QC shoud be FALSE because qc is NULL")
  if (is.null(info))
    stop("Info must not be NULL")

  #select the subject in info and need QC or not
  index <- NULL
  for (i in seq_along(info)) {
    index1 <- as.character(info[[i]])
    index <- c(index, index1)
  }

  if (length(which(index == "")) != 0)  {
    index <- index[-which(index == "")]
  }

  index <- index[!is.na(index)]
  index <- match(index, colnames(subject))
  index <- index[!is.na(index)]
  subject <- subject[, index]


  ##discard the subject's name who is not in the subject data
  for (i in seq_along(info)) {
    idx <- as.character(info[[i]])
    idx <- match(idx, colnames(subject))
    idx <- idx[!is.na(idx)]
    info[[i]] <- colnames(subject)[idx]
  }

  if (QC) {
    int <- cbind(subject, qc)
  } else {
    int <- subject
  }

  ifelse(QC, int <- cbind(subject, qc) , int <- subject)
  name <- colnames(int)
  #
  q <- grep("QC", name)

  if (scale.method == "auto") {
    int <- apply(int, 1, function(x) {
      (x - mean(x)) / sd(x)
    })
    int <- t(int)
  }
  if (scale.method == "pareto") {
    int <- apply(int, 1, function(x) {
      (x - mean(x, na.rm = TRUE)) / sqrt(sd(x, na.rm = TRUE))
    })
    int <- t(int)
  }

  if (scale.method == "no") {
    int <- int
  }
  if (scale.method == "center") {
    int <- apply(int, 1, function(x) {
      (x - mean(x, na.rm = TRUE))
      int <- t(int)
    })
  }



  ## PCA analysis
  int.pca <-
    prcomp(t(data.frame(int)),
           retx = TRUE,
           center = FALSE,
           scale = FALSE)
  sample.pca <- int.pca

  SXTpcaData <- list(
    sample.pca = sample.pca,
    subject = subject,
    qc = qc,
    info = info,
    QC = QC,
    scale.method = scale.method
  )
  class(SXTpcaData) <- "SXTpcaData"
  return(SXTpcaData)
}

