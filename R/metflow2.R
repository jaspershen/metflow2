#' @title metflow2
#' @description metflow2
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @export

metflow2 <- function(){
  cat(crayon::green("Thank you for using metflow2!\n"))
  cat(crayon::green("Version 0.1.2 (20201230)\n"))
  cat(crayon::green("Bug fixing\n"))
  cat(crayon::green("More information can be found at https://jaspershen.github.io/metflow2/\n"))
  cat(crayon::green(
    c("                 _    __ _              ___  ", "                | |  / _| |            |__ \\ ",
      "  _ __ ___   ___| |_| |_| | _____      __ ) |", " | '_ ` _ \\ / _ \\ __|  _| |/ _ \\ \\ /\\ / // / ",
      " | | | | | |  __/ |_| | | | (_) \\ V  V // /_ ", " |_| |_| |_|\\___|\\__|_| |_|\\___/ \\_/\\_/|____|",
      "                                             ", "                                             "
    )
    
  ), sep = "\n")
}

.onAttach <- function(libname, pkgname){
  packageStartupMessage(crayon::green(
    "metflow2,
More information can be found at https://jaspershen.github.io/metflow2/
Authors: Xiaotao Shen (shenxt@stanford.edu)
Maintainer: Xiaotao Shen.
Version 0.1.2 (20201230)
Small bug fix."

  ),
  cat(crayon::green(
    c("                 _    __ _              ___  ", "                | |  / _| |            |__ \\ ",
      "  _ __ ___   ___| |_| |_| | _____      __ ) |", " | '_ ` _ \\ / _ \\ __|  _| |/ _ \\ \\ /\\ / // / ",
      " | | | | | |  __/ |_| | | | (_) \\ V  V // /_ ", " |_| |_| |_|\\___|\\__|_| |_|\\___/ \\_/\\_/|____|",
      "                                             ", "                                             "
    )

  ), sep = "\n")
  )
}

# packageStartupMessage(crayon::green(
# "metflow2,
# More information can be found at https://jaspershen.github.io/metflow2/
# Authors: Xiaotao Shen (shenxt@stanford.edu)
# Maintainer: Xiaotao Shen.
# Version 0.0.8 (20200113)"
# ),
# cat(crayon::green(
#   c("                 _    __ _              ___  ", "                | |  / _| |            |__ \\ ",
#     "  _ __ ___   ___| |_| |_| | _____      __ ) |", " | '_ ` _ \\ / _ \\ __|  _| |/ _ \\ \\ /\\ / // / ",
#     " | | | | | |  __/ |_| | | | (_) \\ V  V // /_ ", " |_| |_| |_|\\___|\\__|_| |_|\\___/ \\_/\\_/|____|",
#     "                                             ", "                                             "
#   )
# ), sep = "\n")
# )
# 

# library(cowsay)
#https://onlineasciitools.com/convert-text-to-ascii-art
# writeLines(capture.output(say("Hello"), type = "message"), con = "ascii_art.txt")
# art <- readLines("ascii_art.txt")
# dput(art)
# metflow_logo <-
#   c("                 _    __ _              ___  ", "                | |  / _| |            |__ \\ ",
#     "  _ __ ___   ___| |_| |_| | _____      __ ) |", " | '_ ` _ \\ / _ \\ __|  _| |/ _ \\ \\ /\\ / // / ",
#     " | | | | | |  __/ |_| | | | (_) \\ V  V // /_ ", " |_| |_| |_|\\___|\\__|_| |_|\\___/ \\_/\\_/|____|",
#     "                                             ", "                                             "
#   )
# cat(metflow_logo, sep = "\n")
# 
# 
# library(asciify)
# bayes_img <- ascii_data("bayes.png")      # path to the bayes image
# bayes_map <- ascii_map(file = bayes_img)  # construct ASCII map
# bayes_map
# ascii_plot(bayes_map, charsize = 2)



# cities <- c("S\u00e3o Paulo", "Reykjav\u00edk")
# print(cities)
# ASCIIfy(cities, 1)
# ASCIIfy(cities, 2)
# 
# athens <- "\u0391\u03b8\u03ae\u03bd\u03b1"
# print(athens)
# ASCIIfy(athens)





  
  
