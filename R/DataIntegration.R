#' @title integrateData
#' @description Integrate different batch datasets together.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object object.
#' @param method Which method do you want to use?
#' subject.mean or qc.mean, default is qc.mean.
#' @return Return a MetFlowData which has been integrated.
#' @details Data integration is a necessary step for multiple batch dataset.
#' \href{https://www.readcube.com/library/fe13374b-5bc9-4c61-9b7f-6a354690947e:abe41368-d08d-4806-871f-3aa035d21743}{Dunn's}
#' method has been used in this function.
#' @export
#' @import tidyverse
#' @import plyr
#' @import purrr

## data integration
setGeneric(
  name = "integrateData",
  def = function(object,
                 method = c("qc.mean",
                            "qc.median",
                            "subject.mean",
                            "subject.median")) {
    
    method <- match.arg(method)
    if (class(object) != "metflowClass") {
      stop("Only the metflowClass is supported!\n")
    }
    
    ms1_data <- object@ms1.data
    if(length(ms1_data) > 1){
      stop("Please align batch first.\n")
    }
    
    ms1_data <- ms1_data[[1]]
    sample_info <- object@sample.info
    
    if(
      sample_info$batch %>% 
      unique() %>% 
      length() == 1
    ){
      return(object)
    }
    
    qc_data <- getData(object = object, slot = "QC")
    subject_data <- getData(object = object, slot = "Subject")
    
    if(sum(is.na(qc_data)) +  sum(is.na(subject_data)) > 0){
      stop("Please impute MV first.\n")
    }
    
    if(is.null(qc_data)){
      if(method %in% c("qc.mean", "qc.median")){
        stop("No QC samples in your data, please change other method.\n")
      }
    }
    
    sample_info <- object@sample.info
    sample_info <- 
      sample_info %>% 
      dplyr::filter(class %in% c('QC', 'Subject'))
    
    ms1_data <- object@ms1.data[[1]]
    ms1_data <- 
      ms1_data %>% 
      dplyr::select(one_of(sample_info$sample.name))
    
    ###split data according to batch
    ##sample_info is a list
    sample_info <- 
      plyr::dlply(sample_info, .variables = .(batch))
    
    subject_data <- 
      lapply(sample_info, function(x){
        temp_subject_data <- 
          ms1_data %>% 
          dplyr::select(one_of(x$sample.name[x$class == "Subject"]))
      })
    
    qc_data <- 
      lapply(sample_info, function(x){
        temp_subject_data <- 
          ms1_data %>% 
          dplyr::select(one_of(x$sample.name[x$class == "QC"]))
      })
    
    correct_factor <- 
      ifelse(stringr::str_detect(method, "subject"), list(subject_data), list(qc_data))[[1]] %>% 
      lapply(., 
             function(x){
               apply(x, 1, ifelse(stringr::str_detect(method, "mean"), mean, median))
             }
      )
    
    correct_factor <- 
      lapply(correct_factor[-1], function(x){
        correct_factor[[1]] / x
      })
    
    
    subject_data[-1] <-
      purrr::map2(subject_data[-1], correct_factor, function(x , y){
        x * y
      }) 
    
    subject_data <- 
      subject_data %>% 
      dplyr::bind_cols()
    
    qc_data <- 
      qc_data %>% 
      dplyr::bind_cols()
    
    object@process.info$integrateData <- list()
    object@process.info$integrateData$method <- method
    sample_info <- object@sample.info
    subject_qc_data <- cbind(qc_data, subject_data)
    
    subject_qc_name <- dplyr::filter(.data = sample_info, class %in% c("Subject", "QC")) %>% 
      dplyr::pull(., sample.name)
    
    subject_qc_data <- subject_qc_data[, match(subject_qc_name, colnames(subject_qc_data))]
    ms1_data <- 
      object@ms1.data[[1]]
    ms1_data[,match(subject_qc_name, colnames(ms1_data))] <- subject_qc_data
    ms1_data <- list(ms1_data)
    object@ms1.data <- ms1_data
    invisible(object)
    
  }
)
