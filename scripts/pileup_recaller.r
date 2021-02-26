#Check installation of required packages
if (require("optparse") == FALSE){
  install.packages("optparse")
}

if (require("tidyverse") == FALSE){
  install.packages("tidyverse")
}

#library loading

library("optparse")
library("tidyverse")


option_list <- list(make_option(c("--pileup_file"), default = "NA", type = "character", help = "Pileup file")
                    , make_option(c("--metrics"), default = "NA", type = "character", help = "Metrics file")
                    , make_option(c("--FNs_table"), default = "NA", type = "character", help = "Table of FNs")
                    , make_option(c("--FPs_table"), default = "NA", type = "character", help = "Table of FPs")
                    , make_option(c("--TPs_table"), default = "NA", type = "character", help = "Table of TPs")
                    , make_option(c("--bases"), default = "NA", type = "character", help = "number of bases covered by the panel high conf")
                    , make_option(c("--out"), default = "NA", type = "character", help = "Output directory")                 
)

opt = parse_args(OptionParser(option_list=option_list))

#generate data frames for pileup, TP, FP, FN and metrics

pileup_FNs <- read.delim(opt$pileup_file, 
                         , sep = "\t"
                         , header = F)

TPs <- read.delim(opt$TPs_table
                  , sep = " "
                  , header = F)

FPs <- read.delim(opt$FPs_table
                  , sep = " "
                  , header = F)

FNs <- read.delim(opt$FNs_table
                  , sep = " "
                  , header = F)

metrics <- read.delim(opt$metrics
                  , sep = " "
                  , header = T)

#if provided takes into account the number of high confidence bases covered by the panel
if (opt$bases != "NA"){
    options(scipen=999)
    bases <- read.delim(opt$bases #high conf regions covered by the panel
                    , sep = " "
                    , header = F)
bases$V1 <- as.numeric(bases$V1) 
covered_bases <- bases$V1
}
 
if (is.null(pileup_FNs) == FALSE){
  dp = str_match_all(pileup_FNs$V13, "DP(.*?);")
  dp = sapply(dp, "[[", 1)
  dp = sapply(str_split(dp, "DP="), "[[",2)
  dp = as.numeric(sapply(str_split(dp, ";"), "[[",1))

  median_dp <- median(dp)
  print(paste0("Median Depth in FNs reviewed is ", median_dp))

  #Create the string for the comparison
  pileup_fns <- paste0(paste0(paste0(paste0(paste0(paste0(paste0(paste0(pileup_FNs$V1, "_"), pileup_FNs$V2), "_"), pileup_FNs$V3),"_"), pileup_FNs$V4), "_"), pileup_FNs$V5)
  fns <- paste0(paste0(paste0(paste0(paste0(paste0(paste0(paste0(FNs$V1, "_"), FNs$V2), "_"), FNs$V3), "_"), FNs$V4), "_"), FNs$V5)

  #remove possible duplicates
  string_pileup <- pileup_fns[!duplicated(pileup_fns)]
  string_fns <- fns[!duplicated(fns)]

  #check how many FNs from pileup are in the FNs list
  #length(string_pileup[string_pileup %in% string_fns])
  #length(string_fns) - length(string_pileup[string_pileup %in% string_fns])

  #list of FNs that are present in the bam file
  in_the_bam <- data.frame(gsub("_", " ", string_fns[string_fns %in% string_pileup]))


  #Create the not in function
  '%!in%' <- function(x,y)!('%in%'(x,y))
  #check how many FNs are not present in the bam
  #length(string_pileup[string_pileup %!in% string_fns])

  #list of FNs that are not present in the bam file
  not_in_bam <- data.frame(gsub("_", " ",string_fns[string_fns %!in% string_pileup]))

  #actual FNs are the ones that are not present in the pileup
  FNs_real <- string_fns[string_fns %!in% string_pileup]
  #TPs updated are FNs not present in actual FNs added to TPs already listed
  TPs_to_add <- string_fns[string_fns %!in% FNs_real] 

  FNs_real <- data.frame(gsub("_", " ", FNs_real))




  write.table(not_in_bam, paste0(opt$out, "Variants_not_in_bam.txt")
              , sep = " "
              , col.names = F
              , row.names = F 
              , eol = "\n"
              , quote = F)


  write.table(FNs_real, paste0(opt$out, "FNs_updated.txt")
              , sep = " "
              , col.names = F
              , row.names = F
              , eol = "\n"
              , quote = F)

  write.table(data.frame(gsub("_", " ", TPs_to_add)), paste0(opt$out, "TPs_to_add.txt")
              , sep = " "
              , col.names = F
              , row.names = F
              , eol = "\n"
              , quote = F)

  write.table(in_the_bam, paste0(opt$out, "Variants_in_bam_not_called.txt")
              , sep = " "
              , col.names = F
              , row.names = F
              , eol = "\n"
              , quote = F)

  TPs_to_add<- read.delim(paste0(opt$out, "TPs_to_add.txt")
                    , sep = " "
                    , header = F)

  FNs_real <- read.delim(paste0(opt$out, "FNs_updated.txt")
                    , sep = " "
                    , header = F)


  TPs_updated <- rbind(TPs, TPs_to_add)
  Recall_updated <- as.numeric(dim(TPs_updated)[1])/(as.numeric(dim(TPs_updated)[1])+as.numeric(dim(FNs_real)[1]))
  F1_score_updated = 2*(as.numeric(metrics$Precision)*as.numeric(Recall_updated))/(as.numeric(metrics$Precision)+as.numeric(Recall_updated))
  #set differences from end and start columns in order to detect SNV and INDELs
  FNs_real$V4 <- (as.numeric(FNs_real$V3) - as.numeric(FNs_real$V2))
  FPs$V4 <- (as.numeric(FPs$V3) - as.numeric(FPs$V2))
  TPs_updated$V4 <- (as.numeric(TPs_updated$V3) - as.numeric(TPs_updated$V2))

  #if difference is 0 it is a SNV so, convert it to 1
  TPs_updated$V4[which(TPs_updated$V4 == 0)] <-  1

  FPs$V4[which(FPs$V4 == 0)] <-  1

  FNs_real$V4[which(FNs_real$V4 == 0)] <-  1

  FDR_updated <- dim(FPs)[1]/(dim(TPs_updated)[1]+dim(FPs)[1])

  # total bases to be subtracted from covered bases
  total_bases <- sum(as.numeric(TPs_updated$V4)) + sum(as.numeric(FPs$V4)) + sum(as.numeric(FNs_real$V4))

  #if bases is provided compute the TNs 
  if (opt$bases != "NA"){
      TNs <- (as.numeric(covered_bases) - as.numeric(total_bases))
      
  }



  if (exists("TNs")== TRUE){
      Specificity <- TNs/(TNs + metrics$FP)
      write.table(data.frame(observed = (as.numeric(dim(TPs_updated)[1])+ as.numeric(dim(FPs)[1])), expected = metrics$expected, Recall_updated = round(Recall_updated,3), Precision = round(metrics$Precision,3)
      , Specificity = round(Specificity,3)
      , FDR = round(FDR_updated,3), F1_score = round(F1_score_updated,3), TP = as.numeric(dim(TPs)[1]), TP_bam = as.numeric(dim(TPs_to_add)[1]), TP_updated = as.numeric(dim(TPs_updated)[1]), FP = as.numeric(dim(FPs)[1]), FN_updated = as.numeric(dim(FNs_real)[1]), TN = TNs)
      , paste0(opt$out, "metrics_updated.txt")
      , sep = "\t"
      , col.names = T
      , row.names = F
      , quote = F)

  }else{
      write.table(data.frame(observed = metrics$observed, expected = metrics$expected, Recall = round(Recall_updated,3), Precision = round(metrics$Precision,3)
      , FDR = round(FDR_updated,3), F1_score = round(F1_score_updated,3), TP = as.numeric(dim(TPs)[1]), , TP_bam = as.numeric(dim(TPs_to_add)[1]), TP_updated = as.numeric(dim(TPs_updated)[1]), FP = as.numeric(dim(FPs)[1]), FN_updated = as.numeric(dim(FNs_real)[1]), TN = TNs)
      , paste0(opt$out, "metrics_updated.txt")
      , sep = "\t"
      , col.names = T
      , row.names = F
      , quote = F)

  }

}else{
  print("No FNs... no updated metrics")
}