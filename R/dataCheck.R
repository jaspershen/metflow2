#' @title checkData
#' @description Check data format.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param data MS1 data tables.
#' @param sample.info Name of sample information.
#' @param path Work directory.
#' @return Notice of data checking.
#' @export

setGeneric(
  name = "checkData",
  def = function(data = c("batch1.csv", "batch2.csv"),
                 sample.info = "sample.information.csv",
                 path = ".") {
    
    sample.info.record <- NULL
    cat("Reading data...\n")
    cat("--------------------------------------------------------------\n")
    data <- lapply(data, function(x) {
      as.data.frame(readr::read_csv(
        file.path(path, x),
        col_types = readr::cols(),
        progress = TRUE
      ))
    })

    sample.info <-
      readr::read_csv(file.path(path, sample.info),
                      col_types = readr::cols(),
                      progress = TRUE)
    
    sample.info <- as.data.frame(sample.info)
             #####

    ##check data, component
    data.record <- lapply(seq_along(data), function(x) {
      temp.data.record <- NULL
      data.col.name <- colnames(data[[x]])
      if (all(data.col.name != "name")) {
        cat("Error: No name in data.\n")
        temp.data.record <- c(temp.data.record, "Error")
      } else{
        # cat("OK: The first column of data is name.\n")
        temp.data.record <- c(temp.data.record, "OK")
      }

      if (all(data.col.name != "mz")) {
        cat("Error: No mz in data.\n")
        temp.data.record <- c(temp.data.record, "Error")
      } else{
        # cat("OK: The second column of data is mz.\n")
        temp.data.record <- c(temp.data.record, "OK")
      }

      if (all(data.col.name != "rt")) {
        cat("Error: No rt in data.\n")
        temp.data.record <- c(temp.data.record, "Error")
      } else{
        # cat("OK: The third column of data is not rt.\n")
        temp.data.record <- c(temp.data.record, "OK")
      }
      temp.data.record
    })

    cat("--------------------------------------------------------------\n")
    ##check sample.info
    if (sum(is.na(sample.info)) > 0) {
      cat("Error: There are", sum(is.na(sample.info)), "NAs in you sample.info.\n")
      sample.info.record <- c(sample.info.record, "Error")
    } else{
      # cat("OK: There are no NAs in you sample.info.\n")
      sample.info.record <- c(sample.info.record, "OK")
    }

    if (ifelse(is.na(sum(sample.info == "") > 0), FALSE, sum(sample.info == "") > 0)) {
      cat("Error: There are",
          sum(sample.info == ""),
          "spaces in you sample.info.\n")
      sample.info.record <- c(sample.info.record, "Error")
    } else{
      sample.info.record <- c(sample.info.record, "OK")
    }

    if (colnames(sample.info)[1] != "sample.name") {
      cat("Error: The column 1 of sample information must be sample.name")
      sample.info.record <- c(sample.info.record, "Error")
    } else{
      # cat("OK: There are no spaces in you sample.info.\n")
      sample.info.record <- c(sample.info.record, "OK")
    }

    if (colnames(sample.info)[2] != "injection.order") {
      cat("Error: The column 2 of sample information must be injection.order")
      sample.info.record <- c(sample.info.record, "Error")
    } else{
      # cat("OK: There are no spaces in you sample.info.\n")
      sample.info.record <- c(sample.info.record, "OK")
    }

    if (colnames(sample.info)[3] != "class") {
      cat("Error: The column 3 of sample information must be class")
      sample.info.record <- c(sample.info.record, "Error")
    } else{
      # cat("OK: There are no spaces in you sample.info.\n")
      sample.info.record <- c(sample.info.record, "OK")
    }

    if (colnames(sample.info)[4] != "batch") {
      cat("Error: The column 4 of sample information must be batch")
      sample.info.record <- c(sample.info.record, "Error")
    } else{
      # cat("OK: There are no spaces in you sample.info.\n")
      sample.info.record <- c(sample.info.record, "OK")
      if (any(table(sample.info[, "batch"]) <= 5)) {
        cat("There are some batch only have less than 5 samples. Please check it.")
        sample.info.record <-
          c(sample.info.record, "Error")
      }
    }

    if (colnames(sample.info)[5] != "group") {
      cat("Error: The column 5 of sample information must be group")
      sample.info.record <- c(sample.info.record, "Error")
    } else{
      # cat("OK: There are no spaces in you sample.info.\n")
      sample.info.record <- c(sample.info.record, "OK")
    }

    class <-  unique(as.character(sample.info[, 3]))
    if (all(c("QC", "Subject") %in% class)) {
      sample.info.record <- c(sample.info.record, "OK")
    } else{
      cat("Error: The class of sample.information must be QC and Subject.\n")
      sample.info.record <- c(sample.info.record, "Error")
    }

    sample.idx <- match(sample.info$sample.name,
                        unique(unlist(lapply(data, function(x)
                          colnames(x)))))

    if (any(is.na(sample.idx))) {
      cat("Error: The sample names in sample.inforamtion and data are not same.\n")
      sample.info.record <- c(sample.info.record, "Error")
    } else{
      sample.info.record <- c(sample.info.record, "OK")
    }

    cat("--------------------------------------------------------------\n")
    cat("Summary:\n")

    stat <- lapply(c(data.record, list(sample.info.record)),
                   function(x) {
                     return(c(
                       ifelse(all(x != "Error"), "Valid", "Invalid"),
                       sum(x == "OK"),
                       sum(x == "Warning"),
                       sum(x == "Error")
                     ))
                   })
    stat <- do.call(rbind, stat)
    stat <- as.data.frame(stat)
    colnames(stat) <-
      c("Check result", "OK", "Warning", "Error")
    rownames(stat)[nrow(stat)] <- "sample.info"
    rownames(stat)[1:(nrow(stat) - 1)] <-
      paste("batch", 1:(nrow(stat) - 1), sep = "")

    print(stat, quote = FALSE)
    cat("\n")
    cat("\n")
    cat("data:\n")
    for (idx in 1:length(data.record)) {
      if (all(data.record[idx] != "Error")) {
        cat(paste("Batch", idx), "is valid.\n")
      } else{
        if (sum(data.record[idx] == "Warning") > 0) {
          cat(
            "There",
            ifelse(sum(data.record[idx] == "Warning") > 1, "are", "is"),
            sum(data.record[idx] == "Warning"),
            ifelse(sum(data.record[idx] == "Warning") > 1, "Warnings", "Warning"),
            "in your",
            paste("Batch", idx),
            ".Please check it according to the information.\n"
          )
        }
        
        if (sum(data.record[idx] == "Error") > 0) {
          cat(
            "There",
            ifelse(sum(data.record[idx] == "Error") > 1, "are", "is"),
            sum(data.record[idx] == "Error"),
            ifelse(sum(data.record[idx] == "Error") > 1, "Errors", "Error"),
            "in your",
            paste("Batch", idx),
            ".Please check it according to the information.\n"
          )
        }
      }
    }

    cat("\n")
    cat("sample.info:\n")
    if (all(sample.info.record != "Error")) {
      cat("sample.info is valid.\n")
    } else{
      if (sum(sample.info.record == "Warning") > 0) {
        cat(
          "There",
          ifelse(sum(sample.info.record == "Warning") > 1, "are", "is"),
          sum(sample.info.record == "Warning"),
          ifelse(
            sum(sample.info.record == "Warning") > 1,
            "Warnings",
            "Warning"
          ),
          "in your sample.info. Please check it according to the information.\n"
        )
      }

      if (sum(sample.info.record == "Error") > 0) {
        cat(
          "There",
          ifelse(sum(sample.info.record == "Error") > 1, "are", "is"),
          sum(sample.info.record == "Error"),
          ifelse(sum(sample.info.record == "Error") > 1, "Errors", "Error"),
          "in your sample.info. Please check it according to the information.\n"
        )
      }
      }
    stat <- as.data.frame(stat)
    stat[, 1] <- as.character(stat[, 1])
    stat[, 2] <- as.numeric(stat[, 2])
    stat[, 3] <- as.numeric(stat[, 3])
    invisible(stat)
  })


