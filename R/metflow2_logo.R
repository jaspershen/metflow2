#' @title metflow2_logo
#' @description metflow2_logo
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @importFrom xcms CentWaveParam findChromPeaks adjustRtime ObiwarpParam chromatogram PeakDensityParam groupChromPeaks 
#' @importFrom xcms featureChromatograms groupnames featureDefinitions featureValues fillChromPeaks
#' @importFrom Biobase featureData
#' @importFrom crayon yellow red green bold bgRed
#' @import ggplot2
#' @importFrom pbapply pblapply pboptions
#' @importFrom xml2 read_html read_xml
#' @importFrom rvest html_nodes html_text
#' @importFrom stringr str_split str_replace_all str_trim str_detect str_extract
#' @importFrom dplyr filter select pull everything distinct one_of left_join mutate bind_cols arrange
#' @importFrom tibble as_tibble enframe tibble rownames_to_column
#' @importFrom httr POST content
#' @importFrom clisymbols symbol
#' @importFrom cli rule col_cyan tree
#' @importFrom utils packageVersion object.size write.csv tail
#' @importFrom purrr map map2
#' @importFrom plyr dlply .
#' @importFrom RColorBrewer brewer.pal
#' @importFrom Biobase featureData
#' @importFrom MSnbase readMSData
#' @importFrom readr read_csv cols
#' @importFrom readxl read_excel
#' @importFrom ggrepel geom_text_repel
#' @importFrom tinyTools get_os
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
#' @importFrom e1071 svm
#' @importFrom magrittr %>%
#' @importFrom plotly ggplotly
#' @importFrom BiocGenerics basename
#' @importFrom impute impute.knn
#' @importFrom missForest missForest
#' @importFrom pcaMethods pca completeObs
#' @importFrom patchwork plot_layout
#' @import patchwork
#' @importFrom Rtsne Rtsne
#' @importFrom stats coefficients lm loess median predict 
#' @importFrom stats rgamma rt sd cor p.adjust prcomp t.test wilcox.test
#' @export

metflow2_logo <- function(){
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

  )
  )
}

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


#' # art <- readLines("metflo2_logo.txt")
#' # dput(art)
#' 
#' #https://www.text-image.com/convert/ascii.html
#' 
#' #' @title metflow2_logo
#' #' @description Show logo of metflow2.
#' #' @author Xiaotao Shen
#' #' \email{shenxt@@sioc.ac.cn}
#' #' @param logo.index Which logo you want to use.
#' #' @param colour Logo colour.
#' #' @return Ascii art logo.
#' #' @export
#' #' @import crayon
#' 
#' metflow2_logo <- 
#'   function(logo.index = 1, colour = c("yellow", "white", "green", "red")){
#'     colour <- match.arg(colour)
#'     if(logo.index == 1){
#'       logo <- 
#'         c("                                    `:smms:`                                    ", 
#'           "                                 .+hmy/``/ymh+.                                 ", 
#'           "                             `:smdo-        -odms:                              ", 
#'           "                          .+hmy/`     `://-    `/ymh+.                          ", 
#'           "                       :smdo-       `yMMMMMN+      -odms:                       ", 
#'           "                   .+hmy/`          hMMMMMMMM/        `/ymh+.                   ", 
#'           "                :smdo-              hMMMMMMMM/            -odms:                ", 
#'           "            .+hmy/`                 `mMMMMMm/                `/ymh/.            ", 
#'           "         :sddo-                     sMy-::.                      -oddo-         ", 
#'           "     `/ymy/`                      .dM+                              `/ymy/`     ", 
#'           "  -oddo-                         :Nm-                                   :smdo-  ", 
#'           "ymh+.                       `-:.oMh`                                       .+hmy", 
#'           "Mo                        :dMMMMMN.             .--`                          sM", 
#'           "Mo                       .MMMMMMMMd          `yNMMMMh-                        sM", 
#'           "Mo                       -MMMMMMMMMNNmdhyyso+dMMMMMMMM.                       sM", 
#'           "Mo                        /NMMMMMd-  `.-:/+osmMMMMMMMM-                       sM", 
#'           "Mo                          -/+/.            .dMMMMMN/                        sM", 
#'           "Mo                                            -Nm+/-                          sM", 
#'           "Mo                                           .NN.                             sM", 
#'           "Mo                                          `mN-                              sM", 
#'           "Mo                                         `dM:                               sM", 
#'           "Mo                                     `/oodM/                                sM", 
#'           "Mo                                    yMMMMMMd.                               sM", 
#'           "Mo                                   +MMMMMMMMd                               sM", 
#'           "Mo                                   :MMMMMMMMy                               sM", 
#'           "Mo                                    -hMMMMm+                                sM", 
#'           "Mo                                       ..`                                  sM", 
#'           "Mo                                    `::``-                                  sM", 
#'           "Mo                             -.    sm+/.dd                          syyydy. sM", 
#'           "Mo       ``    ``      `.`    :M:   :M:  -M/    `.`                  `.   `My sM", 
#'           "Mo   hdossmm:yssNh  `odyohm-`yNmyy:yNmyy`yN  `odysym+ .M:  `mM`  +N-     `sN- sM", 
#'           "Mo  .Mh`  +Ms`  yN `mh` .+M/ .M+   .M+  `Mo  dm`   oM. N+ .moM- /N-    :ydo`  sM", 
#'           "Mo  sN`   dd   `Ms /Mhyys+.  sN`   sN`  oM` :M:    yN` ds.m/ N//N.  .sdo.     sM", 
#'           "Mo `Ns   -M/   +M. -Mo   -`  Ny ` `Ns   mh  -Ms  `sN:  ydm:  dhm.  /M+`````   sM", 
#'           "ymysh.   /s    o+   -syyyo   +yy: +M.   +yy` -syyy+`   :s:   /s.   ssssssss-+hmy", 
#'           "  -oddo-                          mh                                    -oddo-  ", 
#'           "     .+hmy/`                      .`                                `/ymy/`     ", 
#'           "         :smdo-                                                  -odds:         ", 
#'           "            .+hmy/`                                          `/ymh+.            ", 
#'           "                :smdo-                                    -odms:                ", 
#'           "                   .+hmy/`                            `/ymh+.                   ", 
#'           "                       :smdo-                      -odms:                       ", 
#'           "                          .+hmy/`              `/ymh+.                          ", 
#'           "                             `:smdo-        -odms:                              ", 
#'           "                                 .+hmy/``/ymh+.                                 ", 
#'           "                                    `:smms:`                                    "
#'         )  
#'     }else{
#'       logo <- 
#'         c("                 _    __ _              ___  ", "                | |  / _| |            |__ \\ ", 
#'           "  _ __ ___   ___| |_| |_| | _____      __ ) |", " | '_ ` _ \\ / _ \\ __|  _| |/ _ \\ \\ /\\ / // / ", 
#'           " | | | | | |  __/ |_| | | | (_) \\ V  V // /_ ", " |_| |_| |_|\\___|\\__|_| |_|\\___/ \\_/\\_/|____|", 
#'           "                                             ", "                                             "
#'         )
#'     }
#'     
#'     
#' 
#'     switch(colour, 
#'            yellow = cat(crayon::yellow(logo), sep = "\n"),
#'            white = cat(crayon::white(logo), sep = "\n"),
#'            green = cat(crayon::green(logo), sep = "\n"),
#'            red = cat(crayon::red(logo), sep = "\n")
#'            )
#'       
#'     
#'     
#'   }




  
  
