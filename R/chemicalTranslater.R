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
#'@export

transID = function(query = "C00001",
                   from = "KEGG",
                   to = "PubChem SID",
                   top = 1,
                   server = c(
                     "http://cts.fiehnlab.ucdavis.edu/service/convert",
                     "https://www.chemspider.com/InChI.asmx"
                   )){
  cat(crayon::yellow("transID() is deprecated, please use the trans_ID() in tinyTools.\n"))
  server <- match.arg(server)
  top <- as.numeric(top)
  if (is.na(top)) {
    top <- 1
  }
  if (server == "http://cts.fiehnlab.ucdavis.edu/service/convert") {
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
    tibble::as_tibble(result)
  } else{
    from_to <- paste(from, to, sep = "_")
    result <- data.frame(
      name = c(
        "csid_mol",
        "inchikey_csid",
        "inchikey_inchi",
        "inchikey_mol",
        "inchi_csid",
        "inchi_inchikey",
        "inchi_mol",
        "inchi_smiles",
        "smiles_inchi"
      ),
      server = c(
        "https://www.chemspider.com/InChI.asmx/CSIDToMol?",
        "https://www.chemspider.com/InChI.asmx/InChIKeyToCSID?",
        "https://www.chemspider.com/InChI.asmx/InChIKeyToInChI?",
        "https://www.chemspider.com/InChI.asmx/InChIKeyToMol?",
        "https://www.chemspider.com/InChI.asmx/InChIToCSID",
        "https://www.chemspider.com/InChI.asmx/InChIToInChIKey",
        "https://www.chemspider.com/InChI.asmx/InChIToMol",
        "https://www.chemspider.com/InChI.asmx/InChIToSMILES",
        "https://www.chemspider.com/InChI.asmx/SMILESToInChI"
      ),
      stringsAsFactors = FALSE
    )
    if (!from_to %in% result$name) {
      cat(
        crayon::red(
          from_to,
          "is not supported for server",
          "https://www.chemspider.com/InChI.asmx\n"
        )
      )
      return(NA)
    }
    baseurl <- 
      result %>%  
      dplyr::filter(name == from_to) %>% 
      dplyr::pull(server)
    # baseurl <-
    #   "https://www.chemspider.com/InChI.asmx/InChIToSMILES"
    Sys.sleep(rgamma(1, shape = 15, scale = 1 / 45))
    
    body <- 
      switch(EXPR = from,
             csid = list(csid = query),
             inchikey = list(inchi_key = query),
             inchi = list(inchi = query),
             smiles = list(smiles = query),
             mol = list(mol = query)
      )
    
    res <-
      try(httr::POST(url = baseurl, 
                     body = body,
                     encode = "form"),
          silent = TRUE)
    
    if (inherits(res, "try-error")) {
      warning("Problem with service... Returning NA.")
      return(NA)
    }
    
    out <- try(xml2::read_xml(httr::content(res, "raw")), silent = TRUE)
    if (inherits(out, "try-error")) {
      warning("inchi not found... Returning NA.")
      return(NA)
    }
    
    out <- xml2::xml_text(out)
    out
    
  }
  
}


#' @title databaseName
#' @description Whate databases are supported.
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param server server.
#' @return A vector..
#' @export
#' @examples 
#' \dontrun{
#' databaseNames(server = "http://cts.fiehnlab.ucdavis.edu/service/convert")
#' databaseNames(server = "https://www.chemspider.com/InChI.asmx")
#' }

databaseName = function(server = c(
  "http://cts.fiehnlab.ucdavis.edu/service/convert",
  "https://www.chemspider.com/InChI.asmx"
)){
  cat(crayon::yellow("databaseName() is deprecated, please use the trans_id_database() in tinyTools.\n"))
  server <- match.arg(server)
  
  if (server == "http://cts.fiehnlab.ucdavis.edu/service/convert") {
    chemical_database_from = xml2::read_html("http://cts.fiehnlab.ucdavis.edu/service/conversion/fromValues")
    
    chemical_database_to = xml2::read_html("http://cts.fiehnlab.ucdavis.edu/service/conversion/toValues")
    
    chemical_database_from <-
      chemical_database_from %>%
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
    
    
    chemical_database_to <-
      chemical_database_to %>%
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
    
    
    chemical_database_from <- sort(chemical_database_from)
    chemical_database_to <- sort(chemical_database_to)
    cat(
      crayon::green(
        length(chemical_database_from),
        "databases are supported in server http://cts.fiehnlab.ucdavis.edu/service/convert for 'from'.\n"
      )
    )
    cat(
      crayon::green(
        length(chemical_database_to),
        "databases are supported in server http://cts.fiehnlab.ucdavis.edu/service/convert for 'to'.\n"
      )
    )
    result <-
      list(
        From =  tibble::as_tibble(
          data.frame(From = chemical_database_from,
                     stringsAsFactors = FALSE)
        ),
        To = tibble::as_tibble(
          data.frame(From = chemical_database_to,
                     stringsAsFactors = FALSE)
        )
      )
    
    result
  } else{
    result <- c(
      "csid_mol",
      "inchikey_csid",
      "inchikey_inchi",
      "inchikey_mol",
      "inchi_csid",
      "inchi_inchikey",
      "inchi_mol",
      "inchi_smiles",
      "smiles_inchi"
    )
    
    result <-
      stringr::str_split(result, "_") %>%
      dplyr::bind_cols() %>%
      t() %>%
      tibble::as_tibble(.name_repair = "minimal")
    colnames(result) <- c("From", "To")
    cat(
      crayon::green(
        nrow(result),
        "are supported in server https://www.chemspider.com/InChI.asmx.\n"
      )
    )
    result
  }
  
}

