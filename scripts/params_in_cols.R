library(stringr)
library(stringi)

params_in_cols = function(set){
  vaf = str_match_all(set$V13, "AF(.*?);")
  vaf = sapply(vaf, "[[", 1)
  vaf = sapply(str_split(vaf, "AF="), "[[",2)
  vaf = as.numeric(sapply(str_split(vaf, ";"), "[[",1))
  
  dp = str_match_all(set$V13, "DP(.*?);")
  dp = sapply(dp, "[[", 1)
  dp = sapply(str_split(dp, "DP="), "[[",2)
  dp = as.numeric(sapply(str_split(dp, ";"), "[[",1))
  
  qd = str_match_all(set$V13, "QD(.*?);")
  qd = sapply(qd, "[", 1)
  qd = sapply(str_split(qd, "QD="), "[",2)
  qd = as.numeric(sapply(str_split(qd, ";"), "[",1))
  
  ao = str_match_all(set$V13, "AO(.*?);")
  ao = sapply(ao, "[", 1)
  ao = sapply(str_split(ao, "AO="), "[",2)
  ao = as.numeric(sapply(str_split(ao, ";"), "[",1))
  
  ro = str_match_all(set$V13, "RO(.*?);")
  ro = sapply(ro, "[", 1)
  ro = sapply(str_split(ro, "RO="), "[",2)
  ro = as.numeric(sapply(str_split(ro, ";"), "[",1))
  
  sssb = str_match_all(set$V13, "SSSB(.*?);")
  sssb = sapply(sssb, "[", 1)
  sssb = sapply(str_split(sssb, "SSSB="), "[",2)
  sssb = as.numeric(sapply(str_split(sssb, ";"), "[",1))
  
  stb = str_match_all(set$V13, "STB(.*?);")
  stb = sapply(stb, "[", 1)
  stb = sapply(str_split(stb, "STB="), "[",2)
  stb = as.numeric(sapply(str_split(stb, ";"), "[",1))
  
  stbp = str_match_all(set$V13, "STBP(.*?);")
  stbp = sapply(stbp, "[", 1)
  stbp = sapply(str_split(stbp, "STBP="), "[",2)
  stbp = as.numeric(sapply(str_split(stbp, ";"), "[",1))
  
  if(dim(set)[1] == length(vaf)) set$VAF = vaf else print("Not same length. Skipping...")
  if(dim(set)[1] == length(dp)) set$DP = dp else print("Not same length. Skipping...")
  if(dim(set)[1] == length(qd)) set$QD = qd else print("Not same length. Skipping...")
  if(dim(set)[1] == length(stb)) set$STB = stb else print("Not same length. Skipping...")
  #if(dim(set)[1] == length(ao)) set$AO = ao else print("Not same length. Skipping...")
  #if(dim(set)[1] == length(ro)) set$RO = ro else print("Not same length. Skipping...")
  #if(dim(set)[1] == length(stbp)) set$STBP = stbp else print("Not same length. Skipping...")
  #if(dim(set)[1] == length(sssb)) set$SSSB = sssb else print("Not same length. Skipping...")
  
  return(set)
}

#'%!in%' <- function(x,y){!('%in%'(x,y))}




