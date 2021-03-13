#' @title doPCA
#' @description PCA analysis for MetFlowData.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object A metflowClass object.
#' @param scale.method Scale method.
#' @param slot Class of data.
#' @export
#' @return A ggplot object.

doPCA = function(object,
                 scale.method = c("no", "auto", "pareto", "center"),
                 slot = c("QC", "Subject")) {
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
  
  data <- sxtScale(df = data, method = scale.method)
  
  data <- data.frame(data, class, stringsAsFactors = FALSE)
  data <- data.frame(data, name = rownames(data),
                     stringsAsFactors = FALSE)
  
  pca_object <- prcomp(select(data, -one_of(c('name', 'class'))))
  
  plot <-
    ggplot2::autoplot(
      pca_object,
      data = data,
      colour = 'class',
      frame = TRUE,
      frame.type = "norm"
    ) +
    # scale_y_continuous(limits = c(0, 0.015)) +
    # scale_x_continuous(limits = c(-0.1, 0.1)) +
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
}
