#' @title Get peak intensity distributation plot
#' @description Get peak intensity distributation plot.
#' @author Xiaotao Shen
#' \email{shenxt@@sioc.ac.cn}
#' @param object A metflowClass object.
#' @param peak_name Peak name.
#' @param interactive Interactive plot or not.
#' @export
#' @return A ggplot2 object.

get_peak_int_distribution = function(object, 
                                     peak_name,
                                     interactive = TRUE){
  # browser()
  if (class(object) != "metflowClass") {
    stop("Only the metflowClass is supported!\n")
  }
  
  if (length(object@ms1.data) > 1) {
    stop("Please align batches first.\n")
  }
  
  if (!peak_name %in% object@ms1.data[[1]]$name) {
    stop(peak_name, " is not in your data.\n")
  }
  
  int <-
    object@ms1.data[[1]] %>%
    dplyr::filter(name == peak_name) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("sample.name")
  
  colnames(int)[2] <- "int"
  
  int <-
    object@sample.info %>%
    dplyr::left_join(int, by = "sample.name") %>%
    dplyr::filter(class %in% c("QC", "Subject")) %>%
    dplyr::mutate(int = as.numeric(int))
  
  plot <-
    int %>%
    ggplot(aes(x = injection.order, y = int)) +
    geom_point(aes(color = class)) +
    # ggsci::scale_color_d3() +
    geom_smooth(aes(color = class)) +
    labs(x = "Injection order", y = 'Intensity') +
    theme_bw() +
    theme(axis.title = element_text(size = 13),
          axis.text = element_text(size = 12))
  
  if(interactive){
    plot <- plotly::ggplotly(plot)
  }
  
  return(plot)
}

