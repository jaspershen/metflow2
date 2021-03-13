#' @title Show the duplicated samples of peak tables from Thermo QE instrument
#' @description Show the duplicated samples of peak tables from Thermo QE instrument.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param peak_table A peak_table from XCMS or other software.
#' @return A list of duplicated samples.
#' @export
#' @importFrom magrittr %>%

show_duplicated_samples = function(
  peak_table
){
  idx <-
    which(stringr::str_detect(colnames(peak_table), "_[0-9]{10,16}"))
  if (length(idx) == 0) {
    return(cat(crayon::red("No duplicated samples.\n")))
  }
  duplicated_name <- colnames(peak_table)[idx]
  duplicated_name <-
    duplicated_name %>%
    stringr::str_replace("_[0-9]{10,16}", "") %>%
    unique()
  
  result <-
    lapply(duplicated_name, function(x) {
      temp_idx1 <- match(x, colnames(peak_table))
      temp_idx2 <- stringr::str_detect(colnames(peak_table), 
                                       paste(x, "_", sep = "")) %>% 
        which()
      temp_idx <- c(temp_idx1, temp_idx2)
      temp_idx <- temp_idx[!is.na(temp_idx)]
      x <-
        colnames(peak_table)[temp_idx]
      x
    })
  names(result) <-
    duplicated_name
  return(result)
}



#' @title Remove the duplicated samples of peak tables from Thermo QE instrument
#' @description Remove the duplicated samples of peak tables from Thermo QE instrument.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param peak_table A peak_table from XCMS or other software.
#' @return A new peak table without duplicated samples.
#' @export
#' @importFrom magrittr %>%

remove_duplicated_name = function(peak_table){
  idx <-
    which(stringr::str_detect(colnames(peak_table), "_[0-9]{10,16}"))
  if (length(idx) == 0) {
    return(cat(crayon::red("No duplicated samples.\n")))
  }
  duplicated_name <- colnames(peak_table)[idx]
  duplicated_name <-
    duplicated_name %>%
    stringr::str_replace("_[0-9]{10,16}", "") %>%
    unique()
  remove_name <- NULL
  for(x in duplicated_name){
    temp_idx1 <- match(x, colnames(peak_table))
    temp_idx2 <- stringr::str_detect(colnames(peak_table), 
                                     paste(x, "_", sep = "")) %>% 
      which()
    temp_idx <- c(temp_idx1, temp_idx2)
    temp_idx <- temp_idx[!is.na(temp_idx)]
    
    cat(crayon::yellow(rep("-", 13)), "\n")
    cat(crayon::yellow("--->"), crayon::green(x), "\n")
    if (length(temp_idx) == 1) {
      cat(crayon::red("Only one sample"), "\n")
      # colnames(peak_table)[temp_idx] <-
      #   colnames(peak_table)[temp_idx] %>%
      #   stringr::str_replace("_[0-9]{10,16}", "")
    } else{
      na_number <-
        apply(peak_table[, temp_idx], 2, function(x)
          sum(is.na(x)))
      info <- data.frame(
        sample = names(na_number),
        NA_number = unname(na_number),
        stringsAsFactors = FALSE
      )
      print(info)
      cat("\n")
      remove_idx <- temp_idx[-which.min(na_number)]
      remove_name <- c(remove_name, colnames(peak_table)[remove_idx])
      # remain_idx <- temp_idx[which.min(na_number)]
      cat(crayon::green(paste(colnames(peak_table)[remove_idx], collapse = ";"),
                        "are removed.\n"))
    }
  }
  
  peak_table <-
    peak_table %>%
    dplyr::select(-c(remove_name))
  
  colnames(peak_table) <-
    colnames(peak_table) %>%
    stringr::str_replace("_[0-9]{10,16}", "")
  peak_table
}

