#'@title transID
#'@description Translate metabolite ID.
#'@author Xiaotao Shen
#'\email{shenxt1990@@163.com}
#'@param query The ID of metabolite you want to translate.
#'@param from The databases of metabolites. Supported database can be shown using databaseName("from").
#'@param to The databases of metabolites. Supported database can be shown using databaseName("to").
#'@param top How many results should be returned?
#'@param server server.
#'@return A data frame.
#'@import xml2
#'@import rvest
#'@import tidyverse
#'@import stringr
#'@export

setGeneric(
  name = "transID",
  def = function(query = "C00001",
                 from = "KEGG",
                 to = "PubChem SID",
                 top = 1,
                 server = "http://cts.fiehnlab.ucdavis.edu/service/convert") {
    top <- as.numeric(top)
    if (is.na(top)) {
      top <- 1
    }
    
    url <- paste(server, from, to, query, sep = "/")
    url <- stringr::str_replace_all(url, " ", "%20")
    
    result <-
      try(expr = xml2::read_html(url, encoding = "UTF-8"),
          silent = TRUE)
    if (any(class(result) %in% "try-error")) {
      warning(
        'Please check you query, from and to again.
You can use databaseName() function to check the databases this package support.'
      )
      result <- NA
    } else{
      result <-
        try(result %>%
              html_nodes("p") %>%
              html_text(trim = TRUE) %>%
              stringr::str_split("\n") %>%
              `[[`(1) %>%
              sapply(function(x) {
                x <- stringr::str_trim(x, "both")
                x <-
                  stringr::str_replace_all(string = x,
                                           pattern = '\"',
                                           replacement = "")
                x <-
                  stringr::str_replace_all(string = x,
                                           pattern = ',',
                                           replacement = "")
                x
              }) %>%
              unname() %>%
              data.frame(name = ., stringsAsFactors = FALSE) %>%
              dplyr::filter(!name %in% c("[", "]", "{", "}", "result:")) %>%
              dplyr::filter(!stringr::str_detect(name, "fromIdentifier|searchTerm|toIdentifier")) %>%
              dplyr::pull(name))
      
      if (any(class(result) %in% "try-error")) {
        result <- NA
      }
    }
    
    if (top > length(result)) {
      top <- length(result)
    }
    result <- result[1:top]
    result <-
      data.frame(query, result, stringsAsFactors = FALSE)
    colnames(result) <- c(from, to)
    result
    
  }
)


#'@title databaseName
#'@description Whate databases are supported.
#'@author Xiaotao Shen
#'\email{shenxt1990@@163.com}
#'@param from.to From or to.
#'@return A vector..
#'@import xml2
#'@import rvest
#'@import tidyverse
#'@import stringr
#'@export


setGeneric(
  name = "databaseName",
  def = function(from.to = c("from", "to")) {
    from.to <- match.arg(from.to)
    if (from.to == "from") {
      chemical_database = xml2::read_html("http://cts.fiehnlab.ucdavis.edu/service/conversion/fromValues")
    } else{
      chemical_database = xml2::read_html("http://cts.fiehnlab.ucdavis.edu/service/conversion/toValues")
    }
    
    chemical_database <-
      chemical_database %>%
      rvest::html_nodes("p") %>%
      rvest::html_text(TRUE) %>%
      stringr::str_split(pattern = "\n") %>%
      `[[`(1) %>%
      sapply(function(x) {
        x <- stringr::str_trim(x, "both")
        x <-
          stringr::str_replace_all(string = x,
                                   pattern = '\"',
                                   replacement = "")
        x <-
          stringr::str_replace_all(string = x,
                                   pattern = ',',
                                   replacement = "")
      }) %>%
      unname() %>%
      data.frame(name = ., stringsAsFactors = FALSE) %>%
      dplyr::filter(!name %in% c("[", "]", "{", "}")) %>%
      dplyr::pull(name)
    chemical_database
    
    
  }
)
