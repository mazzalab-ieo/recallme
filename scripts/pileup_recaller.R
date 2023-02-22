#Check installation of required packages
if (require("optparse") == FALSE){
  install.packages("optparse")
}

if (require("tidyverse") == FALSE){
  install.packages("tidyverse")
}

#library loading

library("optparse")
suppressMessages(library("tidyverse"))


option_list <- list(make_option(c("--pileup_file_snv"), default = "NA", type = "character", help = "Pileup file for SNVs")
                    , make_option(c("--pileup_file_indel"), default = "NA", type = "character", help = "Pileup file for SNVs")
                    , make_option(c("--metrics_snv"), default = "NA", type = "character", help = "Metrics file")
                    , make_option(c("--metrics_indel"), default = "NA", type = "character", help = "Metrics file")
                    , make_option(c("--FNs_table_snv"), default = "NA", type = "character", help = "Table of FNs")
                    , make_option(c("--FNs_table_indel"), default = "NA", type = "character", help = "Table of FNs")
                    , make_option(c("--FPs_table_snv"), default = "NA", type = "character", help = "Table of FPs")
                    , make_option(c("--FPs_table_indel"), default = "NA", type = "character", help = "Table of FPs")
                    , make_option(c("--TPs_table_snv"), default = "NA", type = "character", help = "Table of TPs")
                    , make_option(c("--TPs_table_indel"), default = "NA", type = "character", help = "Table of TPs")
                    , make_option(c("--bases"), default = "NA", type = "character", help = "number of bases covered by the panel high conf")
                    , make_option(c("--out"), default = "NA", type = "character", help = "Output directory")
                    , make_option(c("--caller"), default = "NA", type = "character", help = "Caller which produced the query VCF")                 
)

opt = parse_args(OptionParser(option_list=option_list))
#initialize the rds object
metrics_rds = list()

#generate data frames for pileup, TP, FP, FN and metrics
empty_cols = c("V1", "V2", "V3", "V4", "V5")
#check the existence of the table files and read them 
info = file.info(opt$pileup_file_snv)
if(info$size != 0){
  pileup_FNs_snv <- read.delim(opt$pileup_file_snv, 
                         , sep = "\t"
                         , header = F)
}else{
  pileup_FNs_snv <- data.frame(matrix(nrow = 0, ncol = 5))
  colnames(pileup_FNs_snv) = empty_cols
}

info = file.info(opt$pileup_file_indel)
if(info$size != 0){
  pileup_FNs_indel <- read.delim(opt$pileup_file_indel, 
                         , sep = "\t"
                         , header = F)
}else{
  pileup_FNs_indel <- data.frame(matrix(nrow = 0, ncol = 5))
  colnames(pileup_FNs_indel) = empty_cols
}

info = file.info(opt$TPs_table_snv)
if(info$size != 0){
  TPs_snv <- read.delim(opt$TPs_table_snv
                  , sep = " "
                  , header = F)
}else{
  TPs_snv <- data.frame(matrix(nrow = 0, ncol = 5))
  colnames(TPs_snv) = empty_cols
}

info = file.info(opt$TPs_table_indel)
if(info$size != 0){
  TPs_indel <- read.delim(opt$TPs_table_indel
                  , sep = " "
                  , header = F)
}else{
  TPs_indel <- data.frame(matrix(nrow = 0, ncol = 5))
  colnames(TPs_indel) = empty_cols
}
info = file.info(opt$FPs_table_snv)
if(info$size != 0){
  FPs_snv <- read.delim(opt$FPs_table_snv
                  , sep = " "
                  , header = F)
}else{
  FPs_snv <- data.frame(matrix(nrow = 0, ncol = 5))
  colnames(FPs_snv) = empty_cols
}

info = file.info(opt$FPs_table_indel)
if(info$size != 0){
  FPs_indel <- read.delim(opt$FPs_table_indel
                  , sep = " "
                  , header = F)
}else{
  FPs_indel <- data.frame(matrix(nrow = 0, ncol = 5))
  colnames(FPs_indel) = empty_cols
}

info = file.info(opt$FNs_table_snv)
if(info$size != 0){
  FNs_snv <- read.delim(opt$FNs_table_snv
                  , sep = " "
                  , header = F)
}else{
  FNs_snv <- data.frame(matrix(nrow = 0, ncol = 5))
  colnames(FNs_snv) = empty_cols
}

info = file.info(opt$FNs_table_indel)
if(info$size != 0){
  FNs_indel <- read.delim(opt$FNs_table_indel
                  , sep = " "
                  , header = F)
}else{
  FNs_indel <- data.frame(matrix(nrow = 0, ncol = 5))
  colnames(FNs_indel) = empty_cols
}

metrics_snv <- read.delim(opt$metrics_snv
                  , sep = " "
                  , header = T)

metrics_indel <- read.delim(opt$metrics_indel
                  , sep = " "
                  , header = T)
#store variants into the rds list
metrics_rds[["TPs_snv"]] = TPs_snv
metrics_rds[["FPs_snv"]] = FPs_snv
metrics_rds[["FNs_snv"]] = FNs_snv
metrics_rds[["TPs_indel"]] = TPs_indel
metrics_rds[["FPs_indel"]] = FPs_indel
metrics_rds[["FNs_indel"]] = FNs_indel
#if provided takes into account the number of high confidence bases covered by the panel
if (opt$bases != "NA"){
    options(scipen=999)
    bases <- read.delim(opt$bases #high conf regions covered by the panel
                    , sep = " "
                    , header = F)
bases$V1 <- as.numeric(bases$V1) 
covered_bases <- bases$V1
}

####### SNV metrics ########### 
if (is.null(pileup_FNs_snv) == FALSE){
    
    #Create the string for the comparison
    pileup_fns <- paste0(paste0(paste0(paste0(paste0(paste0(paste0(paste0(pileup_FNs_snv$V1, "_"), pileup_FNs_snv$V2), "_"), pileup_FNs_snv$V3),"_"), pileup_FNs_snv$V4), "_"), pileup_FNs_snv$V5)
    fns <- paste0(paste0(paste0(paste0(paste0(paste0(paste0(paste0(FNs_snv$V1, "_"), FNs_snv$V2), "_"), FNs_snv$V3), "_"), FNs_snv$V4), "_"), FNs_snv$V5)
    
    #remove possible duplicates
    string_pileup <- pileup_fns[!duplicated(pileup_fns)]
    string_fns <- fns[!duplicated(fns)]

    #check how many FNs from pileup are in the FNs list
    #length(string_pileup[string_pileup %in% string_fns])
    #length(string_fns) - length(string_pileup[string_pileup %in% string_fns])

    #list of FNs that are present in the bam file
    in_the_bam <- data.frame(gsub("_", " ", string_fns[string_fns %in% string_pileup]))
    col_name = colnames(in_the_bam)
    in_the_bam = separate(in_the_bam, col = col_name, into = c("V1", "V2", "V3", "V4", "V5"))
    pileup_FNs_snv$V1 = as.character(pileup_FNs_snv$V1)
    pileup_FNs_snv$V2 = as.numeric(pileup_FNs_snv$V2)
    pileup_FNs_snv$V3 = as.numeric(pileup_FNs_snv$V3)
    in_the_bam$V1 = as.character(in_the_bam$V1)
    in_the_bam$V2 = as.numeric(in_the_bam$V2)
    in_the_bam$V3 = as.numeric(in_the_bam$V3)
    in_the_bam <- inner_join(in_the_bam, pileup_FNs_snv, by = c("V1" = "V1","V2" = "V2","V3" = "V3","V4" = "V4","V5" = "V5"))
    colnames(in_the_bam) = c("chr", "position", "position", "ref", "base", "depth", "count", "vaf","positive_strand", "negative_strand", "percent_bias")


    #Create the not in function
    '%!in%' <- function(x,y)!('%in%'(x,y))

    #list of FNs that are not present in the bam file
    not_in_bam <- data.frame(gsub("_", " ",string_fns[string_fns %!in% string_pileup]))

    #actual FNs are the ones that are not present in the pileup
    FNs_real <- string_fns[string_fns %!in% string_pileup]
    #TPs updated are FNs not present in actual FNs added to TPs already listed
    TPs_to_add <- string_fns[string_fns %!in% FNs_real] 

    FNs_real <- data.frame(gsub("_", " ", FNs_real))

    write.table(not_in_bam, paste0(opt$out, "Variants_not_in_bam_snv.txt")
            , sep = " "
            , col.names = F
            , row.names = F 
            , eol = "\n"
            , quote = F)


    write.table(FNs_real, paste0(opt$out, "FNs_updated_snv.txt")
            , sep = " "
            , col.names = F
            , row.names = F
            , eol = "\n"
            , quote = F)

    write.table(data.frame(gsub("_", " ", TPs_to_add)), paste0(opt$out, "TPs_to_add_snv.txt")
            , sep = " "
            , col.names = F
            , row.names = F
            , eol = "\n"
            , quote = F)

    write.table(in_the_bam, paste0(opt$out, "Variants_in_bam_not_called_snv.txt")
            , sep = " "
            , col.names = T
            , row.names = F
            , eol = "\n"
            , quote = F)

    TPs_to_add<- read.delim(paste0(opt$out, "TPs_to_add_snv.txt")
                , sep = " "
                , header = F)
    #check the existence of the table file
    info = file.info(paste0(opt$out, "FNs_updated_snv.txt"))
    if(info$size != 0){
      FNs_real <- read.delim(paste0(opt$out, "FNs_updated_snv.txt")
                , sep = " "
                , header = F)
    }
    #if RecallME found all of the fns in the bam file the FNs_real dataframe is empty, so a check on dataframe existence is needed
    if (!exists('FNs_real')){
      FNs_real = data.frame()
    }
    TPs_updated <- rbind(TPs_snv, TPs_to_add)
    Recall_updated <- as.numeric(dim(TPs_updated)[1])/(as.numeric(dim(TPs_updated)[1])+as.numeric(dim(FNs_real)[1]))
    F1_score_updated = 2*(as.numeric(metrics_snv$Precision)*as.numeric(Recall_updated))/(as.numeric(metrics_snv$Precision)+as.numeric(Recall_updated))
    #set differences from end and start columns in order to detect SNV and INDELs
    if (exists('FNs_real')){
    FNs_real$V4 <- (as.numeric(FNs_real$V3) - as.numeric(FNs_real$V2))
    }
    FPs_snv$V4 <- (as.numeric(FPs_snv$V3) - as.numeric(FPs_snv$V2))
    TPs_updated$V4 <- (as.numeric(TPs_updated$V3) - as.numeric(TPs_updated$V2))
    
    #if difference is 0 it is a SNV so, convert it to 1 (to compute specificity the program has to sum all of the bases covered by the variants so 1 for SNVs)
    TPs_updated$V4[which(TPs_updated$V4 == 0)] <-  1

    FPs_snv$V4[which(FPs_snv$V4 == 0)] <-  1
    if (exists('FNs_real')){
    FNs_real$V4[which(FNs_real$V4 == 0)] <-  1
    }
    FDR_updated <- dim(FPs_snv)[1]/(dim(TPs_updated)[1]+dim(FPs_snv)[1])
    if (exists('FNs_real')){
    # total bases to be subtracted from covered bases
    total_bases <- sum(as.numeric(TPs_updated$V4)) + sum(as.numeric(FPs_snv$V4)) + sum(as.numeric(FNs_real$V4))
    }else{
      total_bases <- sum(as.numeric(TPs_updated$V4)) + sum(as.numeric(FPs_snv$V4))
    }
    #if bases is provided compute the TNs 
    if (opt$bases != "NA"){
    TNs_snv <- (as.numeric(covered_bases) - as.numeric(total_bases))

    }
    if (exists("TNs_snv")){
        Specificity <- TNs_snv/(TNs_snv + metrics_snv$FP)
        metrics = data.frame(observed = metrics_snv$observed, expected = metrics_snv$expected, Recall_sequencing = round(Recall_updated,3), Recall = metrics_snv$Recall, Precision = round(metrics_snv$Precision,3)
        , Specificity = round(Specificity,3)
        , FDR = round(FDR_updated,3), F1_score = round(F1_score_updated,3), TPs = as.numeric(dim(TPs_snv)[1]), TPs_in_bam = as.numeric(dim(TPs_to_add)[1]), TPs_updated = as.numeric(dim(TPs_updated)[1]), FPs = as.numeric(dim(FPs_snv)[1]), FNs_not_in_bam = as.numeric(dim(FNs_real)[1]), TNs = TNs_snv)

        write.table(metrics
        , paste0(opt$out, "metrics_snv_seq_evaluation.txt")
        , sep = "\t"
        , col.names = T
        , row.names = F
        , quote = F)

    }else{
        metrics = data.frame(observed = metrics_snv$observed, expected = metrics_snv$expected, Recall_sequencing = round(Recall_updated,3), Recall = metrics_snv$Recall, Precision = round(metrics_snv$Precision,3)
        , FDR = round(FDR_updated,3), F1_score = round(F1_score_updated,3), TPs = as.numeric(dim(TPs_snv)[1]), , TPs_in_bam = as.numeric(dim(TPs_to_add)[1]), TPs_updated = as.numeric(dim(TPs_updated)[1]), FPs = as.numeric(dim(FPs_snv)[1]), FNs_not_in_bam = as.numeric(dim(FNs_real)[1]), TNs = TNs_snv)

        write.table(metrics
        , paste0(opt$out, "metrics_snv_seq_evaluation.txt")
        , sep = "\t"
        , col.names = T
        , row.names = F
        , quote = F)

    }

    }else{
    print("No FNs... no updated metrics")
    }

metrics_rds[["Variants_in_bam_not_called_snv"]] = in_the_bam
#check the existence of the table file
    info = file.info(paste0(opt$out, "Variants_not_in_bam_snv.txt"))
    if(info$size != 0){
metrics_rds[["Variants_not_in_bam_snv"]] = read.delim(paste0(opt$out, "Variants_not_in_bam_snv.txt"), sep = " ", header = F)
    }else{
      metrics_rds[["Variants_not_in_bam_snv"]] = data.frame()
    }
metrics_rds[["metrics_snv"]] = metrics

####### INDEL metrics ########### 
if (is.null(pileup_FNs_indel) == FALSE){
  
  #Create the string for the comparison
  pileup_fns <- paste0(paste0(paste0(paste0(paste0(paste0(paste0(paste0(pileup_FNs_indel$V1, "_"), pileup_FNs_indel$V2), "_"), pileup_FNs_indel$V3),"_"), pileup_FNs_indel$V4), "_"), pileup_FNs_indel$V5)
  fns <- paste0(paste0(paste0(paste0(paste0(paste0(paste0(paste0(FNs_indel$V1, "_"), FNs_indel$V2), "_"), FNs_indel$V3), "_"), FNs_indel$V4), "_"), FNs_indel$V5)
  
  #remove possible duplicates
  string_pileup <- pileup_fns[!duplicated(pileup_fns)]
  string_fns <- fns[!duplicated(fns)]

  #check how many FNs from pileup are in the FNs list
  #length(string_pileup[string_pileup %in% string_fns])
  #length(string_fns) - length(string_pileup[string_pileup %in% string_fns])

  #list of FNs that are present in the bam file
  in_the_bam <- data.frame(gsub("-", "00",gsub("_", " ", string_fns[string_fns %in% string_pileup]))) # "00" is necessary for the separate step because it is not able to read the "-"
  #in_the_bam <- data.frame(gsub("-", "00", )) 
  col_name = colnames(in_the_bam)
  in_the_bam = separate(in_the_bam, col = col_name, into = c("V1", "V2", "V3", "V4", "V5"))
  in_the_bam$V4 = gsub("00", "-", in_the_bam$V4)
  in_the_bam$V5 = gsub("00", "-", in_the_bam$V5)
  pileup_FNs_indel$V1 = as.character(pileup_FNs_indel$V1)
  pileup_FNs_indel$V2 = as.numeric(pileup_FNs_indel$V2)
  pileup_FNs_indel$V3 = as.numeric(pileup_FNs_indel$V3)
  in_the_bam$V1 = as.character(in_the_bam$V1)
  in_the_bam$V2 = as.numeric(in_the_bam$V2)
  in_the_bam$V3 = as.numeric(in_the_bam$V3)
  in_the_bam <- inner_join(in_the_bam, pileup_FNs_indel, by = c("V1" = "V1","V2" = "V2","V3" = "V3","V4" = "V4","V5" = "V5"))
  colnames(in_the_bam) = c("chr", "position", "position", "ref", "base", "depth", "count", "vaf","positive_strand", "negative_strand", "percent_bias")

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




  write.table(not_in_bam, paste0(opt$out, "Variants_not_in_bam_indel.txt")
              , sep = " "
              , col.names = F
              , row.names = F 
              , eol = "\n"
              , quote = F)


  write.table(FNs_real, paste0(opt$out, "FNs_updated_indel.txt")
              , sep = " "
              , col.names = F
              , row.names = F
              , eol = "\n"
              , quote = F)

  write.table(data.frame(gsub("_", " ", TPs_to_add)), paste0(opt$out, "TPs_to_add_indel.txt")
              , sep = " "
              , col.names = F
              , row.names = F
              , eol = "\n"
              , quote = F)

  write.table(in_the_bam, paste0(opt$out, "Variants_in_bam_not_called_indel.txt")
              , sep = " "
              , col.names = T
              , row.names = F
              , eol = "\n"
              , quote = F)
  info_tps_to_add <- file.info(paste0(opt$out, "TPs_to_add_indel.txt"))
  if (info_tps_to_add$size == 0){
    print("No TPs to add in INDELs")
    TPs_updated <- TPs_indel
    TPs_to_add <- data.frame()
  }else{
    TPs_to_add<- read.delim(paste0(opt$out, "TPs_to_add_indel.txt")
                    , sep = " "
                    , header = F)
    TPs_updated <- rbind(TPs_indel, TPs_to_add)
  }
  info = file.info(paste0(opt$out, "FNs_updated_indel.txt"))
  if(info$size != 0 ){
    FNs_real <- read.delim(paste0(opt$out, "FNs_updated_indel.txt")
                    , sep = " "
                    , header = F)
  }
  

  #if RecallME found all of the fns in the bam file the FNs_real dataframe is empty, so a check on dataframe existence is needed
    if (!exists('FNs_real')){
      FNs_real = data.frame()
    }
  Recall_updated <- as.numeric(dim(TPs_updated)[1])/(as.numeric(dim(TPs_updated)[1])+as.numeric(dim(FNs_real)[1]))
  F1_score_updated = 2*(as.numeric(metrics_indel$Precision)*as.numeric(Recall_updated))/(as.numeric(metrics_indel$Precision)+as.numeric(Recall_updated))
  #set differences from end and start columns in order to detect indel and INDELs
  if (exists('FNs_real')){
  FNs_real$V4 <- (as.numeric(FNs_real$V3) - as.numeric(FNs_real$V2))
  }
  FPs_indel$V4 <- (as.numeric(FPs_indel$V3) - as.numeric(FPs_indel$V2))
  TPs_updated$V4 <- (as.numeric(TPs_updated$V3) - as.numeric(TPs_updated$V2))

  #if difference is 0 it is a SNV so, convert it to 1
  TPs_updated$V4[which(TPs_updated$V4 == 0)] <-  1

  FPs_indel$V4[which(FPs_indel$V4 == 0)] <-  1
  if (exists('FNs_real')){
  FNs_real$V4[which(FNs_real$V4 == 0)] <-  1
  }
  FDR_updated <- dim(FPs_indel)[1]/(dim(TPs_updated)[1]+dim(FPs_indel)[1])

  # total bases to be subtracted from covered bases
  if (exists('FNs_real')){
  total_bases <- sum(as.numeric(TPs_updated$V4)) + sum(as.numeric(FPs_indel$V4)) + sum(as.numeric(FNs_real$V4))
  }else{
    total_bases <- sum(as.numeric(TPs_updated$V4)) + sum(as.numeric(FPs_indel$V4))
  }
  #if bases is provided compute the TNs 
  if (opt$bases != "NA"){
      TNs_indel <- (as.numeric(covered_bases) - as.numeric(total_bases))
      
  }

  if (exists("TNs_indel")== TRUE){
      Specificity <- TNs_indel/(TNs_indel + metrics_indel$FP)
      metrics = data.frame(observed = metrics_indel$observed, expected = metrics_indel$expected, Recall_sequencing = round(Recall_updated,3), Recall = metrics_indel$Recall,Precision = round(metrics_indel$Precision,3)
      , Specificity = round(Specificity,3)
      , FDR = round(FDR_updated,3), F1_score = round(F1_score_updated,3), TPs = as.numeric(dim(TPs_indel)[1]), TPs_in_bam = as.numeric(dim(TPs_to_add)[1]), TPs_updated = as.numeric(dim(TPs_updated)[1]), FPs = as.numeric(dim(FPs_indel)[1]), FNs_not_in_bam = as.numeric(dim(FNs_real)[1]), TNs = TNs_indel)
      write.table(metrics
      , paste0(opt$out, "metrics_indel_seq_evaluation.txt")
      , sep = "\t"
      , col.names = T
      , row.names = F
      , quote = F)

  }else{
      metrics = data.frame(observed = metrics_indel$observed, expected = metrics_indel$expected, Recall_sequencing = round(Recall_updated,3), Recall = metrics_indel$Recall, Precision = round(metrics_indel$Precision,3)
      , FDR = round(FDR_updated,3), F1_score = round(F1_score_updated,3), TPs = as.numeric(dim(TPs_indel)[1]), , TPs_in_bam = as.numeric(dim(TPs_to_add)[1]), TPs_updated = as.numeric(dim(TPs_updated)[1]), FPs = as.numeric(dim(FPs_indel)[1]), FNs_not_in_bam = as.numeric(dim(FNs_real)[1]), TNs = TNs_indel)
      write.table(metrics
      , paste0(opt$out, "metrics_indel_seq_evaluation.txt")
      , sep = "\t"
      , col.names = T
      , row.names = F
      , quote = F)

  }

}else{
  print("No FNs... no updated metrics")
}


metrics_rds[["Variants_in_bam_not_called_indel"]] = in_the_bam
info = file.info(paste0(opt$out, "Variants_not_in_bam_indel.txt"))
    if(info$size != 0){
metrics_rds[["Variants_not_in_bam_indel"]] = read.delim(paste0(opt$out, "Variants_not_in_bam_indel.txt"), sep = " ", header = F)
    }else{
      metrics_rds[["Variants_not_in_bam_indel"]] = data.frame()
    }
metrics_rds[["metrics_indel"]] = metrics
metrics_rds[['caller']] = as.character(opt$caller)
saveRDS(object = metrics_rds, file = paste0(opt$out, "metrics_to_upload.rds"))