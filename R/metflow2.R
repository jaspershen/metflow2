#' @title metflow2
#' @description metflow2
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @importFrom magrittr %>%
#' @export

metflow2 <- function(){
  cat("Thank you for using metflow2!\n")
}

.onAttach <- function(libname, pkgname){
  packageStartupMessage("metflow2,
More information can be found at https://jaspershen.github.io/metflow2/
Authors: Xiaotao Shen (shenxt@stanford.edu)
Maintainer: Xiaotao Shen.
Version 0.0.4 (20190805)")
}

packageStartupMessage("metflow2,
More information can be found at https://jaspershen.github.io/metflow2/
Authors: Xiaotao Shen (shenxt@stanford.edu)
Maintainer: Xiaotao Shen.
Version 0.0.4 (20190805)")
