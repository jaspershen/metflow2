#' @title dotsne
#' @description tSNE analysis for metflowClass object.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object A metflowClass object.
#' @param slot Group of data.
#' @param dims The number of dimensions.
#' @param perplexity perplexity.
#' @param colour.index colour.index
#' @param verbose verbose
#' @return A ggplot object.

setGeneric(name = "dotsne",
           function(object,
                    slot = c("QC", "Subject"),
                    dims = 2,
                    perplexity = 30,
                    colour.index = c("class", "group", "batch"),
                    verbose = TRUE) {
             # requireNamespace("tidyverse")
             # requireNamespace("ggplot2")
             colour.index <- match.arg(colour.index)
             if (class(object) != "metflowClass") {
               stop("Only the metflowClass is supported!\n")
             }
             
             if (length(object@ms1.data) != 1) {
               stop("Please align batches first.\n")
             }
             
             if (any(!slot %in% c("QC", "Subject"))) {
               stop("Slot can only be QC and Subject.\n")
             }
             
             tag <- getData(object = object, slot = "Tags")
             name <- dplyr::pull(.data = tag, name)
             
             data <- lapply(slot, function(x) {
               x <- getData(object = object, slot = x)
               if (is.null(x)) {
                 return(x)
               }
               rownames(x) <- name
               as.data.frame(t(x))
             })
             
             remain_idx <- which(!unlist(lapply(data, is.null)))
             slot <- slot[remain_idx]
             data <- data[remain_idx]
             
             class <- mapply(function(x, y) {
               rep(y, nrow(x))
             },
             x = data,
             y = slot)
             
             class <- unlist(class)
             
             data <- do.call(rbind, data)
             
             if (sum(is.na(data)) != 0) {
               stop("Please impute MV first.\n")
             }
             
             sample_info <- object@sample.info %>%
               filter(., sample.name %in% rownames(data))
             
             sample_info <- sample_info$sample.name %>%
               match(., rownames(data)) %>%
               `[`(sample_info, ., )
             colour.index <- pull(sample_info, var = colour.index)
             # colour.index <- factor(colour.index)
             tsne_object <- Rtsne::Rtsne(
               X = as.matrix(data),
               dims = dims,
               perplexity = perplexity,
               verbose = verbose
             )
             
             Y <- tsne_object$Y
             Y <-
               data.frame(Y, "Class" = colou.index, stringsAsFactors = FALSE)
             
             plot <- ggplot(Y, aes(X1, X2, colour = Class)) +
               geom_point() +
               labs(x = "Dimension 1",
                    y = "Dimension 2") +
               theme_bw() +
               scale_colour_manual(values = c("#E64B35FF", "#4DBBD5FF")) +
               theme(
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 12),
                 legend.title = element_text(size = 15),
                 legend.text = element_text(size = 12),
                 strip.background = element_rect(fill = "#0099B47F"),
                 strip.text = element_text(color = "white", size = 15)
               )
             plot
           })
