if (require("flexdashboard") == FALSE){
    install.packages("flexdashboard", repos='http://cran.us.r-project.org', dependencies = TRUE)
}
if (require("rmarkdown") == FALSE){
    install.packages("rmarkdown", repos='http://cran.us.r-project.org', dependencies = TRUE)
}

if (require("reactable") == FALSE){
    install.packages("reactable", repos='http://cran.us.r-project.org', dependencies = TRUE)
}

if (require("plotly") == FALSE){
    install.packages("plotly", repos='http://cran.us.r-project.org', dependencies = TRUE)
}
args = commandArgs(trailingOnly=TRUE)
rmarkdown::render(args[1], params = list(query = args[2],
TP_snv = args[3],
TP_indels = args[4],
FP_snv = args[5],
FP_indels = args[6],
FN_snv = args[7],
FN_indels = args[8],
metrics_snv = args[9],
metrics_indel = args[10],
caller = args[11]
), "flex_dashboard")

