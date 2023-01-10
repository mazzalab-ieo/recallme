#Output converter for SNVs retrieved by bam-readcont
args = commandArgs(trailingOnly = TRUE)
library(dplyr)
bam_rc_fn = read.delim(args[1],
                       header = T,
                       sep = ",")

#remove zero vaf variants
bam_rc_fn_nozero = bam_rc_fn %>% filter(vaf > 0)

#load the parsed dataframe for putative fns
#create the avinput like dataframe for FNs
tmp = data.frame()
df = data.frame()
for (line in 1:dim(bam_rc_fn_nozero)[1]){
  tmp = bam_rc_fn_nozero[line,c("chr", "position", "position", "ref", "base")]
  tmp$base = gsub("\\)", "",gsub("\\(", "", tmp$base))
  df = rbind(df,tmp)
}

#get only snv
df = df[!grepl("\\+", df$base),]
df = df[!grepl("\\-", df$base),]


#print out the table
write.table(df, args[2],
            sep = "\t",
            col.names = F,
            row.names = F,
            quote = F,
            eol = "\n")

