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
rmarkdown::render(paste0(getwd(),"/scripts/ReportME.Rmd"), params = list(query = args[1],
TP_snv = args[2],
TP_indels = args[3],
FP_snv = args[4],
FP_indels = args[5],
FN_snv = args[6],
FN_indels = args[7],
metrics_snv = args[8],
metrics_indel = args[9],
caller = args[10]
), "flex_dashboard")

