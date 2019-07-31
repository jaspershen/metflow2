#' #' @title DataIntegration
#' #' @description Integrate different batch datasets together.
#' #' @author Xiaotao Shen
#' #' \email{shenxt@@sioc.ac.cn}
#' #' @param MetFlowData MetFlowData.
#' #' @param method Which method do you want to use?
#' #' subject.mean or qc.mean, default is qc.mean.
#' #' @return Return a MetFlowData which has been integrated.
#' #' @details Data integration is a necessary step for multiple batch dataset.
#' #' \href{https://www.readcube.com/library/fe13374b-5bc9-4c61-9b7f-6a354690947e:abe41368-d08d-4806-871f-3aa035d21743}{Dunn's}
#' #' method has been used in this function.
#' #' @examples
#' #' \donttest{
#' #'}
#' 
#' ## data integration
#' integrateData <- function(object,
#'                           method = c("qc.mean",
#'                                      "qc.median",
#'                                      "subject.mean",
#'                                      "subject.median")) {
#'   method <- match.arg(method)
#'   subject <- MetFlowData@"subject"
#'   qc <- MetFlowData@"qc"
#'   subject.info <- MetFlowData@"subject.info"
#'   qc.info <- MetFlowData@"qc.info"
#'   subject.name <- subject.info[, 1]
#'   qc.name <- qc.info[, 1]
#'   
#'   subject.batch <- subject.info[, 4]
#'   qc.batch <- qc.info[, 4]
#'   
#'   if (length(unique(subject.batch)) == 1) {
#'     return(MetFlowData)
#'   }
#'   
#'   subject1 <- list()
#'   qc1 <- list()
#'   
#'   data <- SplitBatch(MetFlowData = MetFlowData)
#'   subject1 <- data[[1]]
#'   qc1 <- data[[2]]
#'   
#'   subject1.mean <- lapply(subject1, function(x) {
#'     apply(x, 1, mean)
#'   })
#'   qc1.mean <- lapply(qc1, function(x) {
#'     apply(x, 1, mean)
#'   })
#'   
#'   ref.qc <- lapply(qc1.mean, function(x) {
#'     qc1.mean[[1]] / x
#'   })
#'   ref.subject <-
#'     lapply(subject1.mean, function(x) {
#'       qc1.mean[[1]] / x
#'     })
#'   
#'   subject2 <- as.list(rep(NA, length(subject1)))
#'   qc2 <- as.list(rep(NA, length(qc1)))
#'   
#'   if (method == "qc.mean") {
#'     for (i in seq_along(subject1)) {
#'       subject2[[i]] <- subject1[[i]] * ref.qc[[i]]
#'       qc2[[i]] <- qc1[[i]] * ref.qc[[i]]
#'     }
#'   }
#'   
#'   if (method == "subject.mean") {
#'     for (i in seq_along(subject1)) {
#'       subject2[[i]] <- subject1[[i]] * ref.subject[[i]]
#'       qc2[[i]] <- qc1[[i]] * ref.qc[[i]]
#'     }
#'   }
#'   
#'   subject3 <- subject2[[1]]
#'   qc3 <- qc2[[1]]
#'   
#'   for (i in 2:length(qc2)) {
#'     qc3 <- cbind(qc3, qc2[[i]])
#'     subject3 <- cbind(subject3, subject2[[i]])
#'   }
#'   
#'   subject3[is.na(subject3)] <- 0
#'   qc3[is.na(qc3)] <- 0
#'   
#'   subject3[is.infinite(subject3)] <- 0
#'   qc3[is.infinite(qc3)] <- 0
#'   
#'   MetFlowData@"subject" <- as.matrix(subject3)
#'   MetFlowData@"qc" <- as.matrix(qc3)
#'   MetFlowData@"data.integration" <- "yes"
#'   MetFlowData@"data.integration.method" <- method
#'   return(MetFlowData)
#' }
