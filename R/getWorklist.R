#' @title getWorklist
#' @description Generate sample and worklist..
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param table.name Sample name.
#' @param instrument instrument.
#' @param each.mode.number each.mode.number.
#' @param randommethod See xcms.
#' @param samplenumber See xcms.
#' @param QCstep See xcms.
#' @param conditionQCnumber See xcms.
#' @param qc.index.from Number of threads.
#' @param dir dir.
#' @param method.path method.path.
#' @param ms1.method.pos ms1.method.pos.
#' @param ms1.method.neg ms1.method.neg.
#' @param ms2.method.pos ms2.method.pos.
#' @param ms2.method.neg ms2.method.neg.
#' @param path The working directory.
#' @return Peak table.
#' @export

setGeneric(
  name = "getWorklist",
  def = function(table.name = "batch.xlsx",
                 instrument = c("Thermo", "Agilent", "AB"),
                 each.mode.number = 32,
                 randommethod = c("no", "position", "injection"),
                 samplenumber = NULL,
                 QCstep = 8,
                 conditionQCnumber = 8,
                 qc.index.from = 1,
                 dir = "D:\\Liang\\data\\PS4U\\HILIC\\batch3\\",
                 method.path = "D:\\Liang\\Method\\urine\\HILIC\\",
                 ms1.method.pos = "ZIC-HILIC_MS_pos",
                 ms1.method.neg = "ZIC-HILIC_MS_neg",
                 ms2.method.pos = c(
                   "ZIC-HILIC_MSMS_pos_NCE25",
                   "ZIC-HILIC_MSMS_pos_NCE25",
                   "ZIC-HILIC_MSMS_pos_NCE25",
                   "ZIC-HILIC_MSMS_pos_NCE25"
                 ),
                 ms2.method.neg = c(
                   "ZIC-HILIC_MSMS_neg_NCE25",
                   "ZIC-HILIC_MSMS_neg_NCE25",
                   "ZIC-HILIC_MSMS_neg_NCE25",
                   "ZIC-HILIC_MSMS_neg_NCE25"
                 ),
                 path = ".") {
    dir.create(path = path, showWarnings = FALSE)
    instrument <- match.arg(instrument)
    randommethod <- match.arg(randommethod)
    options(warn = -1)
    file <- dir(path)
    if(all(file != table.name)){
      stop("No ", table.name, " in ", path)
    }
    if (instrument == "Thermo") {
      batch <- readxl::read_excel(file.path(path, table.name))
      batch <- batch[, 1]
      if (randommethod == "no") {
        ###add position
        position <-
          unlist(lapply(c("B", "G", "R"), function(y) {
            paste(y,
                  unlist(lapply(LETTERS[1:5], function(x) {
                    paste(x, 1:8, sep = "")
                  })),
                  sep = "")
          }))
        
        position <- position[-c(1:8)]
        
        condition.qc <-
          matrix(rep(
            c("Condition_QC", 'Condition_QC', "BA1"),
            conditionQCnumber
          ),
          ncol = 3,
          byrow = TRUE)
        
        condition.qc <-
          as.data.frame(condition.qc, stringsAsFactors = FALSE)
        
        colnames(condition.qc) <-
          c("File.Name", "Sample.ID", 'Position')
        
        blank.qc <-
          data.frame(
            "File.Name" = c("Blank", "QC"),
            "Sample.ID" = c("Blank", "QC"),
            "Position" = c("BA2", "BA1"),
            stringsAsFactors = FALSE
          )
        
        blank <-
          data.frame(
            "File.Name" = c("Blank", "Blank", "Blank"),
            "Sample.ID" = c("Blank", "Blank", "Blank"),
            "Position" = c("BA2", "BA2", "BA2"),
            stringsAsFactors = FALSE
          )
        
        dilution.qc <-
          data.frame(
            "File.Name" = c(
              "DL_QC1_1",
              "DL_QC2_1",
              "DL_QC4_1",
              "DL_QC8_1",
              "DL_QC1_2",
              "DL_QC2_2",
              "DL_QC4_2",
              "DL_QC8_2"
            ),
            "Sample.ID" = c(
              "DL_QC1_1",
              "DL_QC2_1",
              "DL_QC4_1",
              "DL_QC8_1",
              "DL_QC1_2",
              "DL_QC2_2",
              "DL_QC4_2",
              "DL_QC8_2"
            ),
            "Position" = c("BA3", "BA4", "BA5", "BA6", "BA3", "BA4", "BA5", "BA6"),
            stringsAsFactors = FALSE
          )
        
        ms2.qc <-
          data.frame(
            "File.Name" = c(
              "QC_MS2_NCE25_1",
              "QC_MS2_NCE25_2",
              "QC_MS2_NCE25_3",
              "QC_MS2_NCE25_4",
              "QC_MS2_NCE50_1",
              "QC_MS2_NCE50_2",
              "QC_MS2_NCE50_3",
              "QC_MS2_NCE50_4"
            ),
            "Sample.ID" = c(
              "QC_MS2_NCE25_1",
              "QC_MS2_NCE25_2",
              "QC_MS2_NCE25_3",
              "QC_MS2_NCE25_4",
              "QC_MS2_NCE50_1",
              "QC_MS2_NCE50_2",
              "QC_MS2_NCE50_3",
              "QC_MS2_NCE50_4"
            ),
            "Position" = c("BA1", "BA1", "BA1", "BA1", "BA1", "BA1", "BA1", "BA1"),
            stringsAsFactors = FALSE
          )
        
        batch <- cbind(batch, batch)
        batch <- data.frame(batch,
                            "Position" = rep(position, nrow(batch))[1:nrow(batch)],
                            stringsAsFactors = FALSE)
        colnames(batch) <-
          c("File.Name", "Sample.ID", "Position")
        temp.class <- sort(rep(1:nrow(batch), QCstep))
        temp.class <- temp.class[1:nrow(batch)]
        batch <-
          data.frame(temp.class, batch, stringsAsFactors = FALSE)
        batch <-
          plyr::dlply(.data = batch,
                      .variables = plyr::.(temp.class))
        
        batch <- lapply(batch, function(x) {
          rbind(blank.qc, x[, -1, drop = FALSE])
        })
        
        batch <- do.call(rbind, batch)
        batch <- rbind(batch, blank.qc)
        rownames(batch) <- NULL
        
        batch.pos <- data.frame(
          batch,
          Path = paste(dir, "POS", sep = ""),
          Instrument.Method = paste(method.path, ms1.method.pos, sep = "\\")
        )
        
        batch.neg <- data.frame(
          batch,
          Path = paste(dir, "NEG", sep = ""),
          Instrument.Method = paste(method.path, ms1.method.neg, sep = "\\")
        )
        
        batch.pos <-
          batch.pos[, c("File.Name",
                        "Sample.ID",
                        "Path",
                        "Instrument.Method",
                        "Position")]
        batch.neg <-
          batch.neg[, c("File.Name",
                        "Sample.ID",
                        "Path",
                        "Instrument.Method",
                        "Position")]
        
        condition.qc.pos <- data.frame(
          condition.qc,
          "Path" = paste(dir, 'POS', sep = ""),
          Instrument.Method = paste(method.path, ms1.method.pos, sep = "\\"),
          stringsAsFactors = FALSE
        )
        
        condition.qc.neg <- data.frame(
          condition.qc,
          "Path" = paste(dir, 'NEG', sep = ""),
          Instrument.Method = paste(method.path, ms1.method.neg, sep = "\\"),
          stringsAsFactors = FALSE
        )
        
        blank.pos <- data.frame(
          blank,
          "Path" = paste(dir, 'POS', sep = ""),
          Instrument.Method = paste(method.path, ms1.method.pos, sep = "\\"),
          stringsAsFactors = FALSE
        )
        
        blank.neg <- data.frame(
          blank,
          "Path" = paste(dir, 'NEG', sep = ""),
          Instrument.Method = paste(method.path, ms1.method.neg, sep = "\\"),
          stringsAsFactors = FALSE
        )
        
        
        dilution.qc.pos <- data.frame(
          dilution.qc,
          "Path" = paste(dir, 'POS', sep = ""),
          Instrument.Method = paste(method.path, ms1.method.pos, sep = "\\"),
          stringsAsFactors = FALSE
        )
        
        dilution.qc.neg <- data.frame(
          dilution.qc,
          "Path" = paste(dir, 'NEG', sep = ""),
          Instrument.Method = paste(method.path, ms1.method.neg, sep = "\\"),
          stringsAsFactors = FALSE
        )
        
        
        ms2.qc.pos <- data.frame(
          ms2.qc,
          "Path" = paste(dir, 'POS', sep = ""),
          Instrument.Method = paste(method.path, ms2.method.pos, sep = "\\"),
          stringsAsFactors = FALSE
        )
        
        ms2.qc.neg <- data.frame(
          ms2.qc,
          "Path" = paste(dir, 'NEG', sep = ""),
          Instrument.Method = paste(method.path, ms2.method.neg, sep = "\\"),
          stringsAsFactors = FALSE
        )
        
        condition.qc.pos <- condition.qc.pos[, colnames(batch.pos)]
        condition.qc.neg <- condition.qc.neg[, colnames(batch.neg)]
        
        blank.pos <- blank.pos[, colnames(batch.pos)]
        blank.neg <- blank.neg[, colnames(batch.neg)]
        
        dilution.qc.pos <- dilution.qc.pos[, colnames(batch.pos)]
        dilution.qc.neg <- dilution.qc.neg[, colnames(batch.neg)]
        
        ms2.qc.pos <- ms2.qc.pos[, colnames(batch.pos)]
        ms2.qc.neg <- ms2.qc.neg[, colnames(batch.neg)]
        
        batch.pos <-
          rbind(condition.qc.pos,
                dilution.qc.pos,
                ms2.qc.pos,
                batch.pos,
                blank.pos)
        batch.neg <-
          rbind(condition.qc.neg,
                dilution.qc.neg,
                ms2.qc.neg,
                batch.neg,
                blank.neg)
        
        ###rename
        batch.pos[which(batch.pos[, 1] == "Condition_QC"), 1] <-
          paste(batch.pos[which(batch.pos[, 1] == "Condition_QC"), 1],
                1:length(which(batch.pos[, 1] == "Condition_QC")), sep = "_")
        
        batch.pos[which(batch.pos[, 1] == "QC"), 1] <-
          paste("QC", qc.index.from:(sum(batch.pos[, 1] == "QC") + qc.index.from -
                                       1),
                sep = "_")
        
        batch.pos[which(batch.pos[, 1] == "Blank"), 1] <-
          paste(batch.pos[which(batch.pos[, 1] == "Blank"), 1],
                1:length(which(batch.pos[, 1] == "Blank")), sep = "_")
        
        batch.pos[, 2] <- batch.pos[, 1]
        
        
        batch.neg[which(batch.neg[, 1] == "Condition_QC"), 1] <-
          paste(batch.neg[which(batch.neg[, 1] == "Condition_QC"), 1],
                1:length(which(batch.neg[, 1] == "Condition_QC")), sep = "_")
        
        batch.neg[which(batch.neg[, 1] == "QC"), 1] <-
          paste("QC", qc.index.from:(sum(batch.neg[, 1] == "QC") + qc.index.from -
                                       1),
                sep = "_")
        
        batch.neg[which(batch.neg[, 1] == "Blank"), 1] <-
          paste(batch.neg[which(batch.neg[, 1] == "Blank"), 1],
                1:length(which(batch.neg[, 1] == "Blank")), sep = "_")
        
        batch.neg[, 2] <- batch.neg[, 1]
        
        write.csv(batch.pos, file.path(path, "worklist.pos.csv"), row.names = FALSE)
        write.csv(batch.neg, file.path(path, "worklist.neg.csv"), row.names = FALSE)
        
        
        ###pos and neg
        idx_blank <-
          which(stringr::str_detect(batch.pos$File.Name, "^Blank_[0-9]{1,3}"))
        idx_qc <-
          which(stringr::str_detect(batch.pos$File.Name, "^QC_[0-9]{1,3}"))
        
        # batch.pos$File.Name[idx_blank]
        # batch.pos$File.Name[idx_qc]
        
        ##for positive mode
        batch.pos.head <-
          batch.pos[1:(idx_blank[1] - 1),]
        
        batch.pos.tail <-
          batch.pos[(tail(idx_qc, 1) + 1):nrow(batch.pos),]
        
        batch.pos.middle <-
          batch.pos[idx_blank[1]:tail(idx_qc, 1),]
        
        temp_class <- rep(1:(QCstep + 2))
        
        idx1 <- grep("^Blank", batch.pos.middle$File.Name)
        idx1 <- idx1[-length(idx1)]
        
        
        idx1 <-
          c(idx1[seq(1, length(idx1), by = round(each.mode.number / QCstep))], tail(idx1, 1)) %>%
          unique
        
        idx2 <- c(idx1[-1] - 1, nrow(batch.pos.middle))
        
        
        batch.pos.middle <-
          mapply(function(x, y) {
            list(batch.pos.middle[x:y,])
          },
          x = idx1,
          y = idx2)
        
        
        
        
        batch.pos.middle <-
          lapply(batch.pos.middle, function(x) {
            temp.idx <- grep("Blank", x$File.Name)
            temp.idx <- temp.idx[-1]
            x <- x[-temp.idx, , drop = FALSE]
            x
          })
        
        batch.pos.middle <-
          lapply(batch.pos.middle, function(x) {
            if (length(grep("QC", tail(x$File.Name))) == 0) {
              x <- rbind(x, x[grep("QC", x$File.Name)[1], ])
              x
            } else{
              x
            }
          })
        
        
        #rename
        
        
        
        batch.neg.head <-
          batch.neg[1:(idx_blank[1] - 1),]
        
        batch.neg.tail <-
          batch.neg[(tail(idx_qc, 1) + 1):nrow(batch.neg),]
        
        batch.neg.middle <-
          batch.neg[idx_blank[1]:tail(idx_qc, 1),]
        
        temp_class <- rep(1:(QCstep + 2))
        
        idx1 <- grep("^Blank", batch.neg.middle$File.Name)
        idx1 <- idx1[-length(idx1)]
        
        idx1 <-
          c(idx1[seq(1, length(idx1), by = round(each.mode.number / QCstep))], tail(idx1, 1)) %>%
          unique
        
        idx2 <- c(idx1[-1] - 1, nrow(batch.neg.middle))
        
        
        batch.neg.middle <-
          mapply(function(x, y) {
            list(batch.neg.middle[x:y,])
          },
          x = idx1,
          y = idx2)
        
        batch.neg.middle <-
          lapply(batch.neg.middle, function(x) {
            temp.idx <- grep("Blank", x$File.Name)
            temp.idx <- temp.idx[-1]
            x <- x[-temp.idx, , drop = FALSE]
            x
          })
        
        batch.neg.middle <-
          lapply(batch.neg.middle, function(x) {
            if (length(grep("QC", tail(x$File.Name))) == 0) {
              x <- rbind(x, x[grep("QC", x$File.Name)[1], ])
              x
            } else{
              x
            }
          })
        
        batch.neg.head <-
          batch.neg.head[-grep("Condition_QC", batch.neg.head$File.Name),]
        
        
        batch.pos.middle[[1]] <-
          rbind(batch.pos.head, batch.pos.middle[[1]])
        
        batch.neg.middle[[1]] <-
          rbind(batch.neg.head, batch.neg.middle[[1]])
        
        batch.neg.middle[[length(batch.neg.middle)]] <-
          rbind(batch.neg.middle[[length(batch.neg.middle)]],
                batch.neg.tail)
        
        batch <- mapply(function(x, y) {
          list(rbind(x, y))
        },
        x = batch.pos.middle,
        y = batch.neg.middle)
        
        batch <- do.call(rbind, batch)
        
        
        # batch$File.Name
        
        mode <- rep(NA, nrow(batch))
        mode[grep("POS", batch$Path)] <- "POS"
        mode[grep("NEG", batch$Path)] <- "NEG"
        
        blank_pos <-
          which(
            mode == "POS" &
              stringr::str_detect(string = batch$File.Name, pattern = "^Blank_")
          )
        
        blank_neg <-
          which(
            mode == "NEG" &
              stringr::str_detect(string = batch$File.Name, pattern = "^Blank_")
          )
        
        
        QC_pos <-
          which(
            mode == "POS" &
              stringr::str_detect(string = batch$File.Name, pattern = "^QC_[0-9]{1,3}")
          )
        
        QC_neg <-
          which(
            mode == "NEG" &
              stringr::str_detect(string = batch$File.Name, pattern = "^QC_[0-9]{1,3}")
          )
        
        batch$File.Name[blank_pos] <-
          batch$Sample.ID[blank_pos] <-
          paste("Blank", 1:length(blank_pos), sep = "_")
        
        batch$File.Name[blank_neg] <-
          batch$Sample.ID[blank_neg] <-
          paste("Blank", 1:length(blank_neg), sep = "_")
        
        batch$File.Name[QC_pos] <-
          batch$Sample.ID[QC_pos] <-
          paste("QC", qc.index.from:(length(QC_pos) + qc.index.from - 1), sep = "_")
        
        batch$File.Name[QC_neg] <-
          batch$Sample.ID[QC_neg] <-
          paste("QC", qc.index.from:(length(QC_pos) + qc.index.from - 1), sep = "_")
        
        
        
        write.csv(batch, file.path(path, "worklist.csv"), row.names = FALSE)
        
      }
    }
  }
)
