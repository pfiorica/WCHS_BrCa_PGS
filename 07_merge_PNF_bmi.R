#!/usr/bin/env Rscript
# args = commandArgs(trailingOnly=TRUE)
# idx <- as.numeric(args[1])
# paste0("idx=",idx)
# paste0("class of idx is ",class(idx))
library(data.table)
library(openxlsx)
library(tidyverse)

dir_PGS <- "/projects/rpci/wchs/pnfioric/pgs_wchs/WCHS_BrCa_pgs/pgs_weights/bmi_weights"
PGS.files = list.files(dir_PGS, full.names = T)
# PGS.scoreID <- substr(PGS.files, start = 59, stop = 67)
PGS.scoreID <- unlist(regmatches(PGS.files, regexpr("PGS.*(?=\\.txt)", PGS.files, perl = TRUE))) 
# PGS.scoreNM <- gsub(" ", ".", substr(PGS.files, start = 69, stop = nchar(PGS.files) - 4))
#PGS.scoreNM <- c("Breast_Cancer_313", "AABCGS_black", "Gao2022_BrCa", "Gao2022_ERNEG", "Gao2022_ERPOS", "Shieh2023")
PGS.scoreNM <- PGS.scoreID   
# df_PGS_Pathway<-as.data.frame(matrix(NA, ncol = length(PGS.scoreNM)+1, nrow = 1))
# names(df_PGS_Pathway)<-c("samples", paste0("PGS_",PGS.scoreNM))

df_ls <- list()
for (i in 1:length(PGS.scoreID)) {
  for (pop in c("main")) {
    mywd = paste("/projects/rpci/wchs/pnfioric/pgs_wchs/WCHS_BrCa_pgs/", PGS.scoreID[i], "/",  pop, "/", sep = '')
    setwd(mywd)
    cat(getwd(), sep = '\n')
  
    score.files <- list.files()
    score.files <- score.files[!grepl("n", score.files)]
    score.list <- lapply(score.files, function(x) {read.table(x, header = F, sep = " ")})
    for (j in 1:length(score.list)) {
      names(score.list[[j]]) <- c("samples", substr(score.files[j], 1, nchar(score.files[j]) - 10))
    }
    all.score <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "samples", all.x = TRUE),
                      score.list)
    all.score$total <- rowSums(all.score[, -which(names(all.score) == "samples")])
    assign(paste0("all.score.", pop), all.score)
  }
  df_ls[[i]] <- all.score.main[, c('samples', "total")]
  names(df_ls[[i]])[2] <- paste0("PGS_", PGS.scoreNM[i])
  cat(paste0("i=", i, "; ", PGS.scoreNM[i], " is done."), sep = '\n')
}

df_PGS_Pathway <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "samples", all.x = TRUE),
                         df_ls)

geno_to_remove <- c("A006788", "A018806" , "A038363",  "A038456")
for (sample in geno_to_remove) {
  df_PGS_Pathway <- df_PGS_Pathway[!grepl(sample, df_PGS_Pathway$samples), ]
}

library(openxlsx)
write.xlsx(df_PGS_Pathway, "/projects/rpci/wchs/pnfioric/pgs_wchs/WCHS_BrCa_pgs/WCHS_BrCa_BMI_5_PGS.xlsx", colWidths = "auto")
write.csv(df_PGS_Pathway, "/projects/rpci/wchs/pnfioric/pgs_wchs/WCHS_BrCa_pgs/WCHS_BrCa_BMI_5_PGS.csv", row.names = F)

