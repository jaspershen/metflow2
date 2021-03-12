#' @title metflow2
#' @description metflow2
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @export

metflow2 <- function(){
  cat(crayon::yellow(
    "`metflow2()` is deprecated, use `metflow2_logo()`."
  ))  
  cat(crayon::green(
    c("                 _    __ _              ___  ", "                | |  / _| |            |__ \\ ",
      "  _ __ ___   ___| |_| |_| | _____      __ ) |", " | '_ ` _ \\ / _ \\ __|  _| |/ _ \\ \\ /\\ / // / ",
      " | | | | | |  __/ |_| | | | (_) \\ V  V // /_ ", " |_| |_| |_|\\___|\\__|_| |_|\\___/ \\_/\\_/|____|",
      "                                             ", "                                             "
    )
    
  ), sep = "\n")
}

