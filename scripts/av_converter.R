args = commandArgs(trailingOnly = TRUE)
#load the table
bam_rc_fn = read.delim(args[1],
                     header = T,
                     sep = ",")
library(dplyr)
#remove zero vaf variants
bam_rc_fn_nozero = bam_rc_fn %>% filter(count > 0)

#load the parsed dataframe for putative fns
#create the avinput like dataframe for FNs
tmp = data.frame()
df = data.frame()
for (line in 1:dim(bam_rc_fn_nozero)[1]){
  tmp = bam_rc_fn_nozero[line,c("chr", "position", "position","ref", "base", "depth", "count", "vaf","positive_strand", "negative_strand", "percent_bias")]
  tmp$base = gsub("\\)", "",gsub("\\(", "", tmp$base))
  df = rbind(df,tmp)
}

#get the indexes of insertions and deletions

idx_ins = rownames(df[grepl("\\+", df$base),])
idx_del = rownames(df[grepl("\\-", df$base),])



#check the length of indexes vectors to find the putative fn indels in the df
#df_safe = df
#df = df_safe
#convert the columns in characters from factors
df$position = as.character(df$position)
df$position.1 = as.character(df$position.1)
df$ref = as.character(df$ref)
df$base = as.character(df$base)

##if there are indels convert the notation as avinputs
#deletions conversion
if (length(idx_del) == 0){
  df = df
}else if(idx_del > 0){
  for (i in idx_del){
    #print(i)
    base = as.character(gsub("-", "",df[i, "base"]))
    df[i,"ref"] = base
    df[i,"base"] = "-"
  }
}


#insertions conversion
if (length(idx_ins) == 0){
  df = df
}else if(length(idx_ins > 0)) {
  for (i in idx_ins){
    #print(i)
    base = as.character(gsub("\\+", "",df[i, "base"]))
    df[i,"ref"] = "-"
    df[i,"base"] = base
  }
}

#print out the table
write.table(df, args[2],
            sep = "\t",
            col.names = F,
            row.names = F,
            quote = F,
            eol = "\n")




