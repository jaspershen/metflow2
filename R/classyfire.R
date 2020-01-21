#'@title get_metclass
#'@description Get the class information of a metabolite using classyfire.
#'@author Xiaotao Shen
#'\email{shenxt1990@@163.com}
#'@param inchikey The inchikey ID of a metabolite.
#'@param server server.
#'@param sleep Sleep time for system. Second.
#'@return A classyfire class object.
#'@import xml2
#'@import rvest
#'@import tidyverse
#'@import stringr
#'@import crayon
#'@import cli
#'@export

setGeneric(
  name = "get_metclass",
  def = function(inchikey = "QZDWODWEESGPLC-UHFFFAOYSA-N",
                 server = "http://classyfire.wishartlab.com/entities/",
                 sleep = 5) {
    url <- paste(server, inchikey, sep = "")
    Sys.sleep(time = sleep)
    result <- try(expr = xml2::read_html(url), silent = TRUE)
    if (class(result)[1] == "try-error") {
      message(crayon::red(
        clisymbols::symbol$cross,
        inchikey,
        "is not available in website.\nPlease check this link:\n",url,'\n'
      ))
      return(NA)
    }
    
    result <-
      try(result %>%
            rvest::html_nodes(css = ".main_card"),
          silent = TRUE)
    
    if (class(result)[1] == "try-error") {
      message(crayon::red(
        clisymbols::symbol$cross,
        inchikey,
        "is not available in website.\nPlease check this link:\n", url, '\n'
      ))
      return(NA)
    }
    
    result <-
      try(rvest::html_text(x = result, trim = TRUE))
    
    if (class(result)[1] == "try-error") {
      message(crayon::red(
        clisymbols::symbol$cross,
        inchikey,
        "is not available in website.\nPlease check this link:\n",url, '\n'
      ))
      return(NA)
    }
    
    compound_info <-
      try(result[[1]] %>%
            stringr::str_replace_all("\n", "{}") %>%
            stringr::str_split('\\{\\}') %>%
            `[[`(1) %>%
            stringr::str_trim(side = "both") %>%
            tibble::enframe(name = NULL) %>% 
            dplyr::filter(value != "") %>%
            dplyr::pull(value) %>%
            lapply(function(x) {
              if (x %in% c("SMILES", "InChIKey", "Formula", "Mass")) {
                tibble::tibble(name = x, value = .[which(x == .) + 1])
              }
            }) %>%
            do.call(rbind, .) %>%
            tibble::as_tibble() %>%
            dplyr::distinct(name, value),
          silent = TRUE)
    
    classification_info <-
      try(result[[2]] %>%
            stringr::str_replace_all("\n", "{}") %>%
            stringr::str_split('\\{\\}') %>%
            `[[`(1) %>%
            stringr::str_trim(side = "both") %>%
            tibble::enframe(name = NULL) %>% 
            dplyr::filter(value != "") %>%
            dplyr::filter(!value %in% c("Taxonomic Classification", "Taxonomy Tree")) %>%
            dplyr::pull(value),
          silent = TRUE)
    
    idx <- try(classification_info %>%
                 `==`("Kingdom") %>%
                 which())
    
    # if (length(idx) == 2) {
    taxonomy_tree <-
      try(classification_info[idx[1]:(idx[2] - 1)] %>%
            matrix(ncol = 2, byrow = TRUE) %>%
            tibble::as_tibble(.name_repair = "minimal"))
    try(colnames(taxonomy_tree) <- c("name", "value"))
    # }
    
    if (class(taxonomy_tree)[1] == "try-error") {
      taxonomy_tree <-
        tibble::tibble(
          name = c("Kingdom",
                   "Superclass",
                   "Class",
                   "Subclass"),
          value = rep(NA, 4)
        )
      message(crayon::red(
        clisymbols::symbol$cross,
        inchikey,
        "is not available in website.\n"
      ))
    }else{
      message(crayon::green(
        clisymbols::symbol$tick,
        inchikey,
        "is available in website.\n"
      ))
    }
    
    
    classification_info <-
      try(result[[2]] %>%
            stringr::str_replace_all("\n", "{}") %>%
            stringr::str_split('\\{\\}') %>%
            `[[`(1) %>%
            stringr::str_trim(side = "both") %>%
            tibble::enframe(name = NULL) %>% 
            dplyr::filter(value != "") %>%
            dplyr::pull(value) %>%
            lapply(function(x) {
              if (x %in% c(
                "Kingdom",
                "Superclass",
                "Class",
                "Subclass",
                "Intermediate Tree Nodes",
                "Direct Parent",
                "Alternative Parents",
                "Molecular Framework",
                "Substituents"
              )) {
                tibble::tibble(name = x, value = .[which(x == .) + 1])
              }
            }) %>%
            do.call(rbind, .) %>%
            tibble::as_tibble() %>%
            dplyr::distinct(name, value),
          silent = TRUE)
    
    description <-
      try(result[[3]] %>%
            stringr::str_replace_all("\n", "{}") %>%
            stringr::str_split('\\{\\}') %>%
            `[[`(1) %>%
            stringr::str_trim(side = "both") %>%
            tibble::enframe(name = NULL) %>% 
            dplyr::filter(value != "") %>%
            dplyr::pull(value) %>%
            lapply(function(x) {
              if (x %in% c("Description")) {
                tibble::tibble(name = x, value = .[which(x == .) + 1])
              }
            }) %>%
            do.call(rbind, .) %>%
            tibble::as_tibble() %>%
            dplyr::distinct(name, value),
          silent = TRUE)
    
    external_descriptors <-
      try(result[[4]] %>%
            stringr::str_replace_all("\n", "{}") %>%
            stringr::str_split('\\{\\}') %>%
            `[[`(1) %>%
            stringr::str_trim(side = "both") %>%
            tibble::enframe(name = NULL) %>% 
            dplyr::filter(value != "") %>%
            dplyr::pull(value) %>%
            lapply(function(x) {
              if (x %in% c("External Descriptors")) {
                tibble::tibble(name = x, value = .[which(x == .) + 1])
              }
            }) %>%
            do.call(rbind, .) %>%
            tibble::as_tibble() %>%
            dplyr::distinct(name, value),
          silent = TRUE)
    
    if (class(compound_info)[1] == "try-error") {
      compound_info <-
        tibble::tibble(
          name = c("SMILES", "InChIKey", "Formula", "Mass"),
          value = rep(NA, 4)
        )
    }
    
    if (class(classification_info)[1] == "try-error") {
      classification_info <-
        tibble::tibble(
          name = c(
            "Kingdom",
            "Superclass",
            "Class",
            "Subclass",
            "Intermediate Tree Nodes",
            "Direct Parent",
            "Alternative Parents",
            "Molecular Framework",
            "Substituents"
          ),
          value = rep(NA, 9)
        )
    }
    
    if (class(description)[1] == "try-error") {
      description <- tibble::tibble(name = "Description",
                                    value = NA)
    }
    
    if (class(external_descriptors)[1] == "try-error") {
      external_descriptors <-  tibble::tibble(name = "External Descriptors",
                                              value = NA)
    }
    
    result <- new(Class = "classyfire")
    result@compound_info <- compound_info
    result@taxonomy_tree <- taxonomy_tree
    result@classification_info <- classification_info
    result@description <- description
    result@external_descriptors <- external_descriptors
    
    return(result)
    
  }
)

setClass(
  Class = 'classyfire',
  representation = representation(
    compound_info = 'tbl_df',
    taxonomy_tree = 'tbl_df',
    classification_info = 'tbl_df',
    description = 'tbl_df',
    external_descriptors = 'tbl_df'
  )
)



setMethod('show',
          signature = 'classyfire',
          function(object) {
            cat(cli::rule(
              left = crayon::bold('classyfire Object'),
              right = paste0('metflow2 v', utils::packageVersion('metflow2'))
            ), '\n')
            
            cat(crayon::red(
              'Object Size:',
              format(utils::object.size(object), units = 'Kb'),
              '\n',
              '\n'
            ))
            
            cat(crayon::green('Information:'), '\n')
            
            cat('SMILES: ', dplyr::pull(object@compound_info, "value")[1], '\n')
            cat('InChIKey: ', dplyr::pull(object@compound_info, "value")[2], '\n')
            cat('Formula: ', dplyr::pull(object@compound_info, "value")[3], '\n')
            cat('Mass: ', dplyr::pull(object@compound_info, "value")[4], '\n')
            
            
            tree_list <-
              object@taxonomy_tree
            
            tree_df <- data.frame(
              stringsAsFactors = FALSE,
              id = tree_list$name,
              connections = I(c(
                as.list(tree_list$name)[-1], list(character(0))
              ))
            )
            
            tree_df$label <-
              paste0(crayon::bold(tree_df$id),
                     ' : ',
                     cli::col_cyan(tree_list$value))
            
            print(cli::tree(tree_df))
            
          })
