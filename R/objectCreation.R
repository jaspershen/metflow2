#' @title creatMetflowObject
#' @description Creat metflowClass object.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param ms1.data MS1 peak table name.
#' @param sample.information Sample information name.
#' @param path Work directory.
#' @return A metflowClass object.
#' @export

setGeneric(
  name = "creatMetflowObject",
  def = function(ms1.data,
                 sample.information,
                 path = ".") {
    check.result <- checkData(data = ms1.data,
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
      version = "0.0.5"
    )
    invisible(object)
  }
)

##S4 class for function metIdentification
setClass(
  Class = "metflowClass",
  representation(
    ms1.data = "list",
    sample.info = "data.frame",
    process.info = "list",
    version = "character"
  )
)


setMethod(
  f = "show",
  signature = "metflowClass",
  definition = function(object) {
    # requireNamespace("magrittr")
    cat(paste(rep("-", 20), collapse = ""), "\n")
    cat("metflow2 version:", object@version, "\n")
    cat(paste(rep("-", 20), collapse = ""), "\n")
    cat("MS1 data\n")
    cat("There are",
        length(object@ms1.data),
        "peak tables in your MS1 data.\n")
    info <- lapply(object@ms1.data, dim) %>%
      do.call(rbind, .)
    
    colnames(info) <- c("Peak.number", "Column.number")
    rownames(info) <- paste("Batch", 1:nrow(info), sep = "")
    print(info)
    cat(paste(rep("-", 20), collapse = ""), "\n")
    cat("There are",
        nrow(object@sample.info),
        "samples in your MS1 data.\n")
    class_info <-
      as.data.frame(table(object@sample.info$class), stringsAsFactors = FALSE)
    colnames(class_info) <- c("Class", "Number")
    print(class_info)
    cat(paste(rep("-", 20), collapse = ""), "\n")
    group_info <-
      as.data.frame(table(object@sample.info$group), stringsAsFactors = FALSE)
    colnames(group_info) <- c("Group", "Number")
    print(group_info)
    cat(paste(rep("-", 20), collapse = ""), "\n")
    cat("Processing\n")
    if (.hasSlot(object = object, name = "process.info")) {
      process.info <- object@process.info
      mapply(function(x, y) {
        cat(x, paste(rep("-", 10), collapse = ""), "\n")
        y <- y[which(names(y) != "plot")]
        y <- data.frame(names(unlist(y)), unlist(y))
        colnames(y) <- c("Parameter", "Value")
        rownames(y) <- NULL
        print(y)
      },
      x = names(process.info),
      y = process.info)
    } else{
      cat("There are no processing for your data.\n")
    }
  }
)


#' @title getData
#' @description Get data from metflowClass object.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object A metflowClass object.
#' @param slot Class of data.
#' @return A data frame.
#' @export
setGeneric(
  name = "getData",
  def = function(object,
                 slot = c("Subject", "QC", "QC.DL", "Blank", "Tags")) {
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
)


#' @title getMVplot4sample
#' @description get MV plot of subject samples.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object A metflowClass object.
#' @return A ggplot2 object.
#' @export
setGeneric(
  name = "getMVplot4sample",
  def = function(object) {
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
)

#' @title calRSD
#' @description Calculate RSD of peaks.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object A metflowClass object.
#' @param slot Class of data.
#' @return A data frame with RSD.
#' @export
setGeneric(
  name = "calRSD",
  def = function(object, slot = c("Subject", "QC")) {
    slot <- match.arg(slot)
    if (class(object) != "metflowClass") {
      stop("Only the metflowClass is supported!\n")
    }
    data <- getData(object = object, slot = slot)
    if (sum(is.na(data)) != 0) {
      stop("Please impute MV first!\n")
    }
    
    rsd <- apply(data, 1, function(x) {
      x <- as.numeric(x)
      sd(x) * 100 / mean(x)
    })
    
    rsd <- data.frame(
      index = 1:length(rsd),
      name = object@ms1.data[[1]]$name,
      rsd,
      stringsAsFactors = FALSE
    )
    invisible(rsd)
  }
)


#' @title getParams
#' @description Get parameters from a metflowClass object.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object A metflowClass object.
#' @return A data frame of parameters.
#' @export
setGeneric(
  name = "getParams",
  def = function(object) {
    # requireNamespace("tidyverse")
    if (class(object) != "metflowClass") {
      stop("Only the metflowClass is supported!\n")
    }
    process_info <- object@process.info
    if (length(process_info) == 0) {
      cat("No process for this dataset.\n")
      return(NULL)
    }
    
    process_info <-
      lapply(process_info, function(x) {
        x <- x[which(names(x) != "plot")]
        x <-
          x %>%
          unlist(.) %>%
          data.frame(., stringsAsFactors = FALSE) %>%
          data.frame(rownames(.), ., stringsAsFactors = FALSE)
        
        rownames(x) <- NULL
        colnames(x) <- c("Parameter", "Value")
        x
      })
    process_info
  }
)

