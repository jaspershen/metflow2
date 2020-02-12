#' @title Show the duplicated samples of peak tables from Thermo QE instrument
#' @description Show the duplicated samples of peak tables from Thermo QE instrument.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param peak_table A peak_table from XCMS or other software.
#' @return A list of duplicated samples.
#' @export
#' @importFrom magrittr %>% 

setGeneric(
  name = "show_duplicated_samples",
  def = function(peak_table) {
    idx <-
      which(stringr::str_detect(colnames(peak_table), "_20[0-9]{12}"))
    if (length(idx) == 0) {
      return(cat(crayon::red("No duplicated samples.\n")))
    }
    duplicated_name <- colnames(peak_table)[idx]
    duplicated_name <-
      duplicated_name %>%
      stringr::str_replace("_20[0-9]{12}", "") %>%
      unique()
    
    result <-
      lapply(duplicated_name, function(x) {
        x <-
          stringr::str_detect(colnames(peak_table), paste(x, "_", sep = "")) %>%
          which() %>%
          `[`(colnames(peak_table), .)
        x
      })
    names(result) <-
      duplicated_name
    return(result)
  }
)

#' @title Remove the duplicated samples of peak tables from Thermo QE instrument
#' @description Remove the duplicated samples of peak tables from Thermo QE instrument.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param peak_table A peak_table from XCMS or other software.
#' @return A new peak table without duplicated samples.
#' @export
#' @importFrom magrittr %>% 
setGeneric(
  name = "remove_duplicated_name",
  def = function(peak_table) {
    idx <-
      which(stringr::str_detect(colnames(peak_table), "_20[0-9]{12}"))
    if (length(idx) == 0) {
      return(cat(crayon::red("No duplicated samples.\n")))
    }
    duplicated_name <- colnames(peak_table)[idx]
    duplicated_name <-
      duplicated_name %>%
      stringr::str_replace("_20[0-9]{12}", "") %>%
      unique()
    
    
    lapply(duplicated_name, function(x) {
      temp_idx <-
        stringr::str_detect(colnames(peak_table), paste(x, "_", sep = "")) %>%
        which()
      cat(crayon::yellow(rep("-", 13)), "\n")
      cat(crayon::yellow("--->"), crayon::green(x), "\n")
      if (length(temp_idx) == 1) {
        cat(crayon::red("Only one sample"), "\n")
        colnames(peak_table)[temp_idx] <-
          colnames(peak_table)[temp_idx] %>%
          stringr::str_replace("_20[0-9]{12}", "")
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
        # remain_idx <- temp_idx[which.min(na_number)]
        cat(crayon::green(colnames(peak_table)[remove_idx], " are removed.\n"))
        peak_table <- 
          peak_table %>% 
          select(-colnames(peak_table)[remove_idx])
      }
      
    })
    
    colnames(peak_table) <- 
      colnames(peak_table) %>% 
      stringr::str_replace("_20[0-9]{12}", "")
    peak_table
  }
)


