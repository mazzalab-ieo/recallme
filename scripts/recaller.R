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

#parsing arguments
option_list <- list(make_option(c("-g", "--ground_truth"), default = "NA", type = "character", help = "Ground truth file")
    , make_option(c("-v", "--query_vcf"), default = "NA", type = "character", help = "Called VCF to compare")
    , make_option(c("--query_vaf"), default = "NA", type = "character", help = "Set VAF threshold for query VCF")
    , make_option(c("--gt_vaf"), default = "NA", type = "character", help = "Set VAF threshold for ground_truth VCF")
    , make_option(c("--out"), default = "NA", type = "character", help = "Output directory")
    , make_option(c("--caller"), default = "NA", type = "character", help = "Caller which produced the query VCF")                 
)

opt = parse_args(OptionParser(option_list=option_list))

print("Loading files...")

#loading files
ground_truth = read.delim(opt$ground_truth
    , sep = "\t"
    , header = F)

query_vcf = read.delim(opt$query_vcf
    , sep = "\t"
    , header = F)

print("Files loaded!")
print("Extracting VAF and DP from query VCF...")

#assign name to type column
colnames(query_vcf)[ncol(query_vcf)] <- "type"
colnames(ground_truth)[ncol(ground_truth)] <- "type"
query_vcf$type <- as.character(query_vcf$type)
ground_truth$type <- as.character(ground_truth$type)

if (opt$caller == "GATK" || opt$caller == "TVC" || opt$caller == "LoFreq" || opt$caller == "Freebayes"){
    vaf = str_match_all(query_vcf$V13, "AF(.*?);")
    vaf = sapply(vaf, "[[", 1)
    vaf = sapply(str_split(vaf, "AF="), "[[",2)
    vaf = as.numeric(sapply(str_split(vaf, ";"), "[[",1))

    dp = str_match_all(query_vcf$V13, "DP(.*?);")
    dp = sapply(dp, "[[", 1)
    dp = sapply(str_split(dp, "DP="), "[[",2)
    dp = as.numeric(sapply(str_split(dp, ";"), "[[",1))

    query_vcf$VAF = vaf
    query_vcf$DP = dp

}else if(opt$caller == "Deepvariant"){
    vaf = strsplit(as.character(query_vcf$V15), ":")
    vaf = sapply(vaf, "[[", 5)
    vaf = as.numeric(vaf)

    dp = strsplit(as.character(query_vcf$V15), ":")
    dp = sapply(dp, "[[", 3)
    dp = as.numeric(dp)

    query_vcf$VAF = vaf
    query_vcf$DP = dp
}else if (opt$caller == "VarScan2"){
    vaf = strsplit(as.character(query_vcf$V15), ":")
    vaf = sapply(vaf, "[[", 7)
    vaf = gsub("%", "", vaf)
    vaf = as.numeric(vaf)

    dp = strsplit(as.character(query_vcf$V15), ":")
    dp = sapply(dp, "[[", 4)
    dp = as.numeric(dp)

    query_vcf$VAF = vaf
    query_vcf$DP = dp

}else if (opt$caller == "VarDict"){
    vaf = strsplit(as.character(query_vcf$V15), ":")
    vaf = sapply(vaf, "[[", 7)
    vaf = as.numeric(vaf)

    dp = strsplit(as.character(query_vcf$V15), ":")
    dp = sapply(dp, "[[", 2)
    dp = as.numeric(dp)

    query_vcf$VAF = vaf
    query_vcf$DP = dp

}

print("VAF and DP extracted and added.")

if (opt$query_vaf == "NA"){
    call = query_vcf

}else{
    query_vcf = query_vcf %>% filter(query_vcf$VAF >= as.numeric(opt$query_vaf))
    call = query_vcf
}

if (opt$gt_vaf == "NA"){
    gt = ground_truth

}else{
    if (opt$caller == "GATK" || opt$caller == "TVC" || opt$caller == "LoFreq"){
    vaf = str_match_all(ground_truth$V13, "AF(.*?);")
    vaf = sapply(vaf, "[[", 1)
    vaf = sapply(str_split(vaf, "AF="), "[[",2)
    vaf = as.numeric(sapply(str_split(vaf, ";"), "[[",1))

    dp = str_match_all(ground_truth$V13, "DP(.*?);")
    dp = sapply(dp, "[[", 1)
    dp = sapply(str_split(dp, "DP="), "[[",2)
    dp = as.numeric(sapply(str_split(dp, ";"), "[[",1))

    ground_truth$VAF = vaf
    ground_truth$DP = dp

}else if(opt$caller == "Deepvariant"){
    vaf = strsplit(as.character(ground_truth$V15), ":")
    vaf = sapply(vaf, "[[", 5)
    vaf = as.numeric(vaf)

    dp = strsplit(as.character(ground_truth$V15), ":")
    dp = sapply(dp, "[[", 3)
    dp = as.numeric(dp)

    ground_truth$VAF = vaf
    ground_truth$DP = dp
}else if (opt$caller == "VarScan2"){
    vaf = strsplit(as.character(ground_truth$V15), ":")
    vaf = sapply(vaf, "[[", 7)
    vaf = gsub("%", "", vaf)
    vaf = as.numeric(vaf)

    dp = strsplit(as.character(ground_truth$V15), ":")
    dp = sapply(dp, "[[", 4)
    dp = as.numeric(dp)

    ground_truth$VAF = vaf
    ground_truth$DP = dp

}else if (opt$caller == "VarDict"){
    vaf = strsplit(as.character(ground_truth$V15), ":")
    vaf = sapply(vaf, "[[", 7)
    vaf = as.numeric(vaf)

    dp = strsplit(as.character(ground_truth$V15), ":")
    dp = sapply(dp, "[[", 2)
    dp = as.numeric(dp)

    ground_truth$VAF = vaf
    ground_truth$DP = dp

}
    ground_truth = ground_truth %>% filter(ground_truth$VAF >= as.numeric(opt$gt_vaf))
    gt = ground_truth
}
#create variables for SNVs and for INDELs
call_snv <- call %>% filter(type == "SNV")
gt_snv <- gt %>% filter(type == "SNV")
call_indel <- call %>% filter(type == "INDEL")
gt_indel <- gt %>% filter(type == "INDEL")
call_string <- call_snv
gt_string <- gt_snv
print("Preparing for comparison...")
call_string <- paste0(paste0(paste0(paste0(paste0(paste0(paste0(paste0(call_string$V1, "_"), call_string$V2), "_"), call_string$V3), "_"), call_string$V4), "_"), call_string$V5)
gt_string <- paste0(paste0(paste0(paste0(paste0(paste0(paste0(paste0(gt_string$V1, "_"), gt_string$V2), "_"), gt_string$V3), "_"), gt_string$V4), "_"), gt_string$V5)
print("Removing duplicates...")
call_string <- call_string[!duplicated(call_string)]
gt_string <- gt_string[!duplicated(gt_string)]
#create function "not in"
'%!in%' <- function(x,y)!('%in%'(x,y))
#compute metrics for SNVs
Recall <- length(call_string[call_string %in% gt_string])/length(gt_string)
Precision <- length(call_string[call_string %in% gt_string])/(length(call_string[call_string %in% gt_string])+length(call_string[call_string %!in% gt_string]))
FDR <- length(call_string[call_string %!in% gt_string])/length(call_string)
F1_score <- 2*(Precision * Recall)/(Precision + Recall)

TP <- call_string[call_string %in% gt_string]
FP <- call_string[call_string %!in% gt_string]
FN <- gt_string[gt_string %!in% call_string]

observed <- length(call_string)
expected <- length(gt_string)

write.table(data.frame(observed = observed, expected = expected, Recall = round(Recall,3), Precision = round(Precision, 3)
, FDR = round(FDR,3), F1_score = round(F1_score,3), TPs = length(TP), FPs = length(FP), FNs = length(FN))
, paste0(opt$out, "bioinfo_metrics_snv.txt")
, sep = " "
, col.names = T
, row.names = F
, quote = F)

FN <- gsub("_", " ", FN)
TP <- gsub("_", " ", TP)
FP <-gsub("_", " ", FP)
TP <- data.frame(TP)
write.table(TP
, paste0(opt$out, "TPs_snv.txt")
, sep = " "
, col.names = F
, row.names = F
, quote = F)

FP <- data.frame(FP)
write.table(FP
, paste0(opt$out, "FPs_snv.txt")
, sep = " "
, col.names = F
, row.names = F
, quote = F)

FN <- data.frame(FN)
write.table(FN
, paste0(opt$out, "FNs_snv.txt")
, sep = " "
, col.names = F
, row.names = F
, quote = F)

###### INDELS metrics

call_string <- call_indel
gt_string <- gt_indel

print("Preparing for comparison...")
call_string <- paste0(paste0(paste0(paste0(paste0(paste0(paste0(paste0(call_string$V1, "_"), call_string$V2), "_"), call_string$V3), "_"), call_string$V4), "_"), call_string$V5)
gt_string <- paste0(paste0(paste0(paste0(paste0(paste0(paste0(paste0(gt_string$V1, "_"), gt_string$V2), "_"), gt_string$V3), "_"), gt_string$V4), "_"), gt_string$V5)
print("Removing duplicates...")
call_string <- call_string[!duplicated(call_string)]
gt_string <- gt_string[!duplicated(gt_string)]

#create function "not in"
#'%!in%' <- function(x,y)!('%in%'(x,y))


Recall <- length(call_string[call_string %in% gt_string])/length(gt_string)
Precision <- length(call_string[call_string %in% gt_string])/(length(call_string[call_string %in% gt_string])+length(call_string[call_string %!in% gt_string]))
FDR <- length(call_string[call_string %!in% gt_string])/length(call_string)
F1_score <- 2*(Precision * Recall)/(Precision + Recall)

TP <- call_string[call_string %in% gt_string]
FP <- call_string[call_string %!in% gt_string]
FN <- gt_string[gt_string %!in% call_string]

observed <- length(call_string)
expected <- length(gt_string)


print("Generating table with metrics...")
write.table(data.frame(observed = observed, expected = expected, Recall = round(Recall,3), Precision = round(Precision, 3)
, FDR = round(FDR,3), F1_score = round(F1_score,3), TPs = length(TP), FPs = length(FP), FNs = length(FN))
, paste0(opt$out, "bioinfo_metrics_indel.txt")
, sep = " "
, col.names = T
, row.names = F
, quote = F)

print("Generating table for TPs, FPs and FNs...")

FN <- gsub("_", " ", FN)
TP <- gsub("_", " ", TP)
FP <-gsub("_", " ", FP)

TP <- data.frame(TP)
write.table(TP
, paste0(opt$out, "TPs_indel.txt")
, sep = " "
, col.names = F
, row.names = F
, quote = F)

FP <- data.frame(FP)
write.table(FP
, paste0(opt$out, "FPs_indel.txt")
, sep = " "
, col.names = F
, row.names = F
, quote = F)

FN <- data.frame(FN)
write.table(FN
, paste0(opt$out, "FNs_indel.txt")
, sep = " "
, col.names = F
, row.names = F
, quote = F)
