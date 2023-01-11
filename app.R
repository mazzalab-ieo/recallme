#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shiny.semantic)
library(ggpubr)
library(dplyr)
library(stringr)
library(stringi)
library(semantic.dashboard)
library(reactable)
library(gridExtra)
source(paste0(getwd(),"/scripts/params_in_cols.R"))


#create boxes with not in bam variants
bam_out = tags$div(
    tags$div(style="display: inline-block;vertical-align:top; width: 50%;",box(shinycssloaders::withSpinner(reactableOutput("snv_nobam")),title = "SNV variants not in bam", color = "red", collapsible = T)),
    tags$div(style="display: inline-block;vertical-align:top; width: 50%;",box(shinycssloaders::withSpinner(reactableOutput("indel_nobam")), title = "INDEL variants not in bam", color = "red",collapsible = T)),
    tags$div(style="display: inline-block;vertical-align:top; width: 50%;",box(shinycssloaders::withSpinner(reactableOutput("snv_bam")),title = "SNV variants in bam", color = "green", collapsible = T)),
    tags$div(style="display: inline-block;vertical-align:top; width: 50%;",box(shinycssloaders::withSpinner(reactableOutput("indel_bam")), title = "INDEL variants in bam", color = "green",collapsible = T))
)

#diagnostic plots
diag_plot_snv = tags$div(
    tags$div(style="display: inline-block;vertical-align:top; width: 1000px;",box(shinycssloaders::withSpinner(plotOutput("diag_plot_snv"), type = 1),title = "SNV Parameters Diagnostics", color = "blue", collapsible = T))
)

diag_plot_indel = tags$div(
    tags$div(style="display: inline-block;vertical-align:top; width: 1000px;",box(shinycssloaders::withSpinner(plotOutput("diag_plot_indel")),title = "INDEL Parameters Diagnostics", color = "blue", collapsible = T))
)

#reactable wilcox test params
tab_wilcox_snv = tags$div(style="display: inline-block;vertical-align:top; width: 1000px;",box(shinycssloaders::withSpinner(reactableOutput("wilcox_tab_snv")),title = "SNV - Mann-Whitney Test", color = "blue", collapsible = T))

tab_wilcox_indel = tags$div(style="display: inline-block;vertical-align:top; width: 1000px;",box(shinycssloaders::withSpinner(reactableOutput("wilcox_tab_indel")),title = "INDEL - Mann-Whitney Test", color = "blue", collapsible = T))

##create div and boxes for tp, fp and fn in SNVs tabs
snv_tabs = tags$div(
    downloadLink("download_tp_snv", "Download"),
    tags$div(style="display: inline-block;vertical-align:top; width: 100%;",box(shinycssloaders::withSpinner(reactableOutput("tp_tab_snv")),title = "TP variants in SNV", color = "green", collapsible = T)),
    downloadLink("download_fp_snv", "Download"),
    tags$div(style="display: inline-block;vertical-align:top; width: 100%;",box(shinycssloaders::withSpinner(reactableOutput("fp_tab_snv")), title = "FP variants in SNV", color = "orange",collapsible = T)),
    downloadLink("download_fn_snv", "Download"),
    tags$div(style="display: inline-block;vertical-align:top; width: 100%;",box(shinycssloaders::withSpinner(reactableOutput("fn_tab_snv")), title = "FN variants in SNV", color = "red",collapsible = T))
    
)

##create div and boxes for tp, fp and fn in INDELs tabs
indel_tabs = tags$div(
    downloadLink("download_tp_indel", "Download"),
    tags$div(style="display: inline-block;vertical-align:top; width: 100%;",box(shinycssloaders::withSpinner(reactableOutput("tp_tab_indel")),title = "TP variants in INDEL", color = "green", collapsible = T)),
    downloadLink("download_fp_indel", "Download"),
    tags$div(style="display: inline-block;vertical-align:top; width: 100%;",box(shinycssloaders::withSpinner(reactableOutput("fp_tab_indel")), title = "FP variants in INDEL", color = "orange",collapsible = T)),
    downloadLink("download_fn_indel", "Download"),
    tags$div(style="display: inline-block;vertical-align:top; width: 100%;",box(shinycssloaders::withSpinner(reactableOutput("fn_tab_indel")), title = "FN variants in INDEL", color = "red",collapsible = T))
    
)

#create content for the first page (Upload)
cont_tab1 = tags$div(
    tags$div(style="display: inline-block;vertical-align:top; width: 30%;",file_input("query", label = "Query callset",button_label = "avinput", placeholder = "Select")),
    tags$div(style="display: inline-block;vertical-align:top; width: 30%;",file_input("gt", label = "GT callset", button_label = "avinput", placeholder = "Select")),
    tags$div(style="display: inline-block;vertical-align:top; width: 40%;",text_input("folder", label = "ResultsDir", value = "path/to/metrics/")),
    tags$div(style="display:inline-block",bam_out)
)
#create contents for parameters page
##slider inputs for parameters
cont_params = tags$div(
    sliderInput(
        inputId = "vaf_snv",
        label = "Set VAF threshold",
        min = 0,
        max = 1,
        value = 0.01,
        step = 0.001
        
    ), 
    tags$br(),
    sliderInput(
        inputId = "dp_snv",
        label = "Set DP threshold",
        min = 0,
        max = 1000,
        value = 50,
        step = 1
        
    ),
    
    tags$br(),
    sliderInput(
        inputId = "qd_snv",
        label = "Set QD threshold",
        min = 0,
        max = 100,
        value = 9,
        step = 0.5
        
    ),tags$br(),
    sliderInput(
        inputId = "stb_snv",
        label = "Set STB threshold",
        min = 0,
        max = 1,
        value = 0.95,
        step = 0.01
        
    )
)

cont_params_indel = tags$div(
    sliderInput(
        inputId = "vaf_indel",
        label = "Set VAF threshold",
        min = 0,
        max = 1,
        value = 0.01,
        step = 0.01
        
    ),
    tags$br(),
    
    sliderInput(
        inputId = "dp_indel",
        label = "Set DP threshold",
        min = 0,
        max = 1000,
        value = 50,
        step = 1
        
    ),
    tags$br(),
    
    sliderInput(
        inputId = "qd_indel",
        label = "Set QD threshold",
        min = 0,
        max = 100,
        value = 9,
        step = 0.5
        
    ),
    tags$br(),
    sliderInput(
        inputId = "stb_indel",
        label = "Set STB threshold",
        min = 0,
        max = 1,
        value = 0.95,
        step = 0.01
        
    )
)


params_snv = tags$div(
                  tags$div(style="display:inline-block", cont_params,reactableOutput("snv"))
)


params_indel = tags$div(
    tags$div(style="display:inline-block", cont_params_indel, reactableOutput("indel"))
)

footer = div(
    class = "footer",
    includeHTML("www/footer.html")
)

# Define UI for application that draws a histogram
ui <- dashboardPage(
    dashboardHeader(color = "blue", title = "RecallME v.0.1", inverted = TRUE,logo_path = "recallme.png", logo_align = "left"),
    dashboardSidebar(
        size = "thin", color = "teal",
        sidebarMenu(
            menuItem(tabName = "upload", "Upload"),
            menuItem(tabName = "diagnostics", "Diagnostics"),
            menuItem(tabName = "set_params-snv", "Set parameters - SNV"),
            menuItem(tabName = "set_params-indel", "Set parameters - INDEL")
        )
    ),
    dashboard_body(
        #create body items
        tabItems(
            tabItem(
                tabName = "upload", #upload and tables page
                cont_tab1,
                tags$br(),
                footer
                
            ),
            
            tabItem(
                tabName = "diagnostics", #diagnostic page
                diag_plot_snv,
                diag_plot_indel,
                tab_wilcox_snv,
                tab_wilcox_indel,
                tags$br(),
                footer
            ),
            
            tabItem(
                tabName = "set_params-snv", #SNV parameters page
                params_snv,
                tags$div(style="display: inline-block;vertical-align:top; width: 20%;",box(title = "TP calls", h1(verbatimTextOutput("tp_snv")), color = "green")),
                tags$div(style="display: inline-block;vertical-align:top; width: 20%;",box(title = "FP calls", h1(verbatimTextOutput("fp_snv")), color = "orange")),
                tags$div(style="display: inline-block;vertical-align:top; width: 20%;",box(title = "FN calls", h1(verbatimTextOutput("fn_snv")), color = "red")),
                tags$br(),
                snv_tabs,
                tags$br(),
                footer
                
            ),
            tabItem(
                tabName = "set_params-indel", #INDEL parameters page
                params_indel,
                tags$div(style="display: inline-block;vertical-align:top; width: 20%;",box(title = "TP calls", h1(verbatimTextOutput("tp_indel")), color = "green")),
                tags$div(style="display: inline-block;vertical-align:top; width: 20%;",box(title = "FP calls", h1(verbatimTextOutput("fp_indel")), color = "orange")),
                tags$div(style="display: inline-block;vertical-align:top; width: 20%;",box(title = "FN calls", h1(verbatimTextOutput("fn_indel")), color = "red")),
                indel_tabs,
                tags$br(),
                footer
        )
                
            
    )
        
    )
)


# Define server logic required to draw a histogram
server <- function(input, output, session) {
    ##increase maximum allowed size for input files to 30 MB
    options(shiny.maxRequestSize=30*1024^2) 
    
    ###get query and GT variants
    query_reactive_snv <- reactive({
        req(input$query)
        query_snv <- read.delim(input$query$datapath, sep="\t", header = F)
        query_snv
    })
    
    query_reactive_indel <- reactive({
        req(input$query)
        query_indel <- read.delim(input$query$datapath, sep="\t", header = F)
        query_indel
    })
    
    gt_reactive_snv <- reactive({
        req(input$gt)
        gt_snv <- read.delim(input$gt$datapath, sep="\t", header = F)
        gt_snv
    })
    
    gt_reactive_indel <- reactive({
        req(input$gt)
        gt_indel <- read.delim(input$gt$datapath, sep="\t", header = F)
        gt_indel
    })
    
    ##metrics snv
    metrics_snv_reac <- reactive({
        req(input$query)
        req(input$gt)
        query_snv <- query_reactive_snv()
        gt_snv <- gt_reactive_snv()
        colnames(query_snv)[ncol(query_snv)] <- "type"
        colnames(gt_snv)[ncol(gt_snv)] <- "type"
        query_snv$type <- as.character(query_snv$type)
        gt_snv$type <- as.character(gt_snv$type)
        
        query_snv = params_in_cols(query_snv)
        
        #query_snv = query_snv %>% filter(query_snv$VAF >= as.numeric(input$vaf_snv))
        #call_snv = query_snv
        
        #check that params columns have proper values
        if (!is.na(query_snv$VAF)){
            query_snv = query_snv %>% filter(query_snv$VAF >= as.numeric(input$vaf_snv))
            call_snv = query_snv
        }else{
            call_snv = query_snv
        }
        
        if (!is.na(query_snv$QD)){
            query_snv = query_snv %>% filter(query_snv$QD >= as.numeric(input$qd_snv))
            call_snv = query_snv
        }else{
            call_snv = query_snv
        }
        
        if (!is.na(query_snv$DP)){
            query_snv = query_snv %>% filter(query_snv$DP >= as.numeric(input$dp_snv))
            call_snv = query_snv
        }else{
            call_snv = query_snv
        }
        
        
        if (!is.na(query_snv$STB)){
            query_snv = query_snv %>% filter(query_snv$STB >= as.numeric(input$STB_snv))
            call_snv = query_snv
        }else{
            call_snv = query_snv
        }
        
        #query_snv = query_snv %>% filter(query_snv$STBP >= as.numeric(input$stbp_snv))
        #call_snv = query_snv
        
        call_snv <- call_snv %>% filter(type == "SNV")
        gt_snv <- gt_snv %>% filter(type == "SNV")
        
        '%!in%' <- function(x,y)!('%in%'(x,y))
        call_snv <- paste0(paste0(paste0(paste0(paste0(paste0(paste0(paste0(call_snv$V1, "_"), call_snv$V2), "_"), call_snv$V3), "_"), call_snv$V4), "_"), call_snv$V5)
        gt_snv <- paste0(paste0(paste0(paste0(paste0(paste0(paste0(paste0(gt_snv$V1, "_"), gt_snv$V2), "_"), gt_snv$V3), "_"), gt_snv$V4), "_"), gt_snv$V5)
        call_snv <- call_snv[!duplicated(call_snv)]
        gt_snv <- gt_snv[!duplicated(gt_snv)]
        #compute metrics for SNVs
        Recall_snv <- length(call_snv[call_snv %in% gt_snv])/length(gt_snv)
        Precision_snv <- length(call_snv[call_snv %in% gt_snv])/(length(call_snv[call_snv %in% gt_snv])+length(call_snv[call_snv %!in% gt_snv]))
        FDR_snv <- length(call_snv[call_snv %!in% gt_snv])/length(call_snv)
        F1_score_snv <- 2*(Precision_snv * Recall_snv)/(Precision_snv + Recall_snv)
        
        TP_snv <- call_snv[call_snv %in% gt_snv]
        FP_snv <- call_snv[call_snv %!in% gt_snv]
        FN_snv <- gt_snv[gt_snv %!in% call_snv]
        
        observed_snv <- length(call_snv)
        expected_snv <- length(gt_snv)
        metrics_snv <- data.frame(observed = observed_snv, expected = expected_snv, Recall = round(Recall_snv,3), Precision = round(Precision_snv, 3), FDR = round(FDR_snv,3), F1_score = round(F1_score_snv,3), TPs = length(TP_snv), FPs =length(FP_snv), FNs = length(FN_snv))
        metrics_snv
    })
 
    ###metrics in INDELs
    metrics_indel_reac <- reactive({
        req(input$query)
        req(input$gt)
        query_indel <- query_reactive_indel()
        gt_indel <- gt_reactive_indel()
        colnames(query_indel)[ncol(query_indel)] <- "type"
        colnames(gt_indel)[ncol(gt_indel)] <- "type"
        query_indel$type <- as.character(query_indel$type)
        gt_indel$type <- as.character(gt_indel$type)
        query_indel = params_in_cols(query_indel)
        
        
        #check that params columns have proper values
        
        if (!is.na(query_indel$VAF)){
            query_indel = query_indel %>% filter(query_indel$VAF >= as.numeric(input$vaf_indel))
            call_indel = query_indel
        }else{
            call_indel = query_indel
        }
        
        if (!is.na(query_indel$QD)){
            query_indel = query_indel %>% filter(query_indel$QD >= as.numeric(input$qd_indel))
            call_indel = query_indel
        }else{
            call_indel = query_indel
        }
        
        if (!is.na(query_indel$DP)){
            query_indel = query_indel %>% filter(query_indel$DP >= as.numeric(input$dp_indel))
            call_indel = query_indel
        }else{
            call_indel = query_indel
        }
        
        
        if (!is.na(query_indel$STB)){
            query_indel = query_indel %>% filter(query_indel$STB >= as.numeric(input$STB_indel))
            call_indel = query_indel
        }else{
            call_indel = query_indel
        }
        
        #query_indel = query_indel %>% filter(query_indel$STBP >= as.numeric(input$stbp_indel))
        #call_indel = query_indel
        
        call_indel <- call_indel %>% filter(type == "INDEL")
        gt_indel <- gt_indel %>% filter(type == "INDEL")
        
        print("Preparing for comparison...")
        call_indel <- paste0(paste0(paste0(paste0(paste0(paste0(paste0(paste0(call_indel$V1, "_"), call_indel$V2), "_"), call_indel$V3), "_"), call_indel$V4), "_"), call_indel$V5)
        gt_indel <- paste0(paste0(paste0(paste0(paste0(paste0(paste0(paste0(gt_indel$V1, "_"), gt_indel$V2), "_"), gt_indel$V3), "_"), gt_indel$V4), "_"), gt_indel$V5)
        print("Removing duplicates...")
        call_indel <- call_indel[!duplicated(call_indel)]
        gt_indel <- gt_indel[!duplicated(gt_indel)]
        #create function "not in"
        '%!in%' <- function(x,y)!('%in%'(x,y))
        #compute metrics for INDELs
        Recall_indel <- length(call_indel[call_indel %in% gt_indel])/length(gt_indel)
        Precision_indel <- length(call_indel[call_indel %in% gt_indel])/(length(call_indel[call_indel %in% gt_indel])+length(call_indel[call_indel %!in% gt_indel]))
        FDR_indel <- length(call_indel[call_indel %!in% gt_indel])/length(call_indel)
        F1_score_indel <- 2*(Precision_indel * Recall_indel)/(Precision_indel + Recall_indel)
        
        TP_indel <- call_indel[call_indel %in% gt_indel]
        FP_indel <- call_indel[call_indel %!in% gt_indel]
        FN_indel <- gt_indel[gt_indel %!in% call_indel]
        
        observed_indel <- length(call_indel)
        expected_indel <- length(gt_indel)
        metrics_indel <- data.frame(observed = observed_indel, expected = expected_indel, Recall = round(Recall_indel,3), Precision = round(Precision_indel, 3), FDR = round(FDR_indel,3), F1_score = round(F1_score_indel,3), TPs = length(TP_indel), FPs =length(FP_indel), FNs = length(FN_indel))
        metrics_indel
    })
    output$snv <- renderReactable(reactable(metrics_snv_reac()))
    output$indel <- renderReactable(reactable(metrics_indel_reac()))
    ####VALUE BOXES IN PARAMETERS PAGES
    #number of TPs per value box snv
    output$tp_snv = renderText({
        req(input$query)
        req(input$gt)
        metrics_snv_reac()[,"TPs"]
    })
    #number of FPs per value box snv
    output$fp_snv = renderText({
        req(input$query)
        req(input$gt)
        metrics_snv_reac()[,"FPs"]
    })
    #number of FNs per value box snv
    output$fn_snv = renderText({
        req(input$query)
        req(input$gt)
        metrics_snv_reac()[,"FNs"]
    })
    
    ###INDELs
    
    #number of TPs per value box snv
    output$tp_indel = renderText({
        req(input$query)
        req(input$gt)
        metrics_indel_reac()[,"TPs"]
    })
    #number of FPs per value box snv
    output$fp_indel = renderText({
        req(input$query)
        req(input$gt)
        metrics_indel_reac()[,"FPs"]
    })
    #number of FNs per value box snv
    output$fn_indel = renderText({
        req(input$query)
        req(input$gt)
        metrics_indel_reac()[,"FNs"]
    })
    
    #Variants not in bam
    var_not_in_bam_snv = reactive({
           print(input$metrics)
           nobam_snv = read.delim(paste0(input$folder, "/Variants_not_in_bam_snv.txt"), sep = " ", header = F)
           colnames(nobam_snv) = c("Chr", "Start", "End", "Ref", "Alt")
           nobam_snv
           })
    
    var_not_in_bam_indel = reactive({
        nobam_indel = read.delim(paste0(input$folder, "/Variants_not_in_bam_indel.txt"), sep = " ", header = F)
        colnames(nobam_indel) = c("Chr", "Start", "End", "Ref", "Alt")
        nobam_indel
    })
    
    #Variants in bam
    var_in_bam_snv = reactive({
        print(input$metrics)
        bam_snv = read.delim(paste0(input$folder, "/TPs_to_add_snv.txt"), sep = " ", header = F)
        colnames(bam_snv) = c("Chr", "Start", "End", "Ref", "Alt")
        bam_snv
    })
    
    var_in_bam_indel = reactive({
        print(input$metrics)
        bam_indel = read.delim(paste0(input$folder, "/TPs_to_add_indel.txt"), sep = " ", header = F)
        colnames(bam_indel) = c("Chr", "Start", "End", "Ref", "Alt")
        bam_indel
    })
    
    output$snv_nobam  = renderReactable({
        if (input$folder == "path/to/metrics/"){
            return(NULL)
        }else{
            reactable(var_not_in_bam_snv()) 
        }
        
        })
        
    output$indel_nobam = renderReactable({
        if (input$folder == "path/to/metrics/"){
            return(NULL)
        }else{
            reactable(var_not_in_bam_indel())
        }
        })
    
    #output variants in bam file 
    output$snv_bam  = renderReactable({
        if (input$folder == "path/to/metrics/"){
            return(NULL)
        }else{
            reactable(var_in_bam_snv()) 
        }
        
    })
    
    output$indel_bam = renderReactable({
        if (input$folder == "path/to/metrics/"){
            return(NULL)
        }else{
            reactable(var_in_bam_indel())
        }
    })
    
##############create tables with variants classes and plots for paramaters
    ##### tables
    query_reac = reactive({
        req(input$query)
        if (is.null(input$query)){
            return(NULL)
        }else{
            query = read.delim(input$query$datapath, sep="\t", header = F)#read whole table
            tmp = as.data.frame(query[,ncol(query)])
            colnames(tmp) = "Type"
            query = query[,c("V1", "V2", "V3", "V4", "V5", "V13")]
            query = cbind(query, tmp)
            #query = query[,c("V1", "V2", "V3", "V4", "V5", "V13", "V16")]#retain VCF columns that are needed
            query = params_in_cols(query)
            query
        }
        
    })
    ##SNVs
    tp_tab_snv = reactive({
        req(input$query)
        req(input$folder)
        tp_snv = read.delim(paste0(input$folder,"/TPs_snv.txt"), header = F, sep = " ")
        tp_snv = inner_join(query_reac(), tp_snv, by = c("V1", "V2", "V3", "V4", "V5"))
        colnames(tp_snv) = c("Chr", "Start", "End", "Ref", "Alt","Info", "Type", "VAF", "DP", "QD", "STB")
        tp_snv$Info = NULL
        tp_snv$Type = NULL
        #colnames(tp_snv) = c("Chr", "Start", "End", "Ref", "Alt","Chr_VCF", "Pos_VCF", "ID_VCF", "Ref_VCF", "Alt_VCF", "Qual_VCF", "Filter_VCF", "Info", "Info2", "Info3","Type")
        #tp_snv = params_in_cols(tp_snv)
        #tp_snv = tp_snv[,c("Chr", "Start", "End", "Ref", "Alt","VAF", "DP", "QD", "STB","Info", "Type")]
        #tp_snv$Info = NULL
        
        if(!is.na(tp_snv$VAF)){
            tp_snv =  tp_snv %>% filter(tp_snv$VAF >= as.numeric(input$vaf_snv))
        }else{
            tp_snv =  tp_snv
        }
        
        if(!is.na(tp_snv$DP)){
            tp_snv =  tp_snv %>% filter(tp_snv$DP >= as.numeric(input$dp_snv))
        }else{
            tp_snv =  tp_snv
        }
        
        if(!is.na(tp_snv$QD)){
            tp_snv =  tp_snv %>% filter(tp_snv$QD >= as.numeric(input$qd_snv))
        }else{
            tp_snv =  tp_snv
        }
        
        if(!is.na(tp_snv$STB)){
            tp_snv =  tp_snv %>% filter(tp_snv$STB >= as.numeric(input$stb_snv))
        }else{
            tp_snv =  tp_snv
        }
        
        tp_snv
    }) 
        
    fp_tab_snv = reactive({
        req(input$query)
        req(input$folder)
        fp_snv = read.delim(paste0(input$folder,"/FPs_snv.txt"), header = F, sep = " ")
        fp_snv = inner_join(query_reac(), fp_snv, by = c("V1", "V2", "V3", "V4", "V5"))
        colnames(fp_snv) = c("Chr", "Start", "End", "Ref", "Alt","Info", "Type", "VAF", "DP", "QD", "STB")
        fp_snv$Info = NULL
        fp_snv$Type = NULL
        if(!is.na(fp_snv$VAF)){
            fp_snv =  fp_snv %>% filter(fp_snv$VAF >= as.numeric(input$vaf_snv))
        }else{
            fp_snv =  fp_snv
        }
        
        if(!is.na(fp_snv$DP)){
            fp_snv =  fp_snv %>% filter(fp_snv$DP >= as.numeric(input$dp_snv))
        }else{
            fp_snv =  fp_snv
        }
        
        if(!is.na(fp_snv$QD)){
            fp_snv =  fp_snv %>% filter(fp_snv$QD >= as.numeric(input$qd_snv))
        }else{
            fp_snv =  fp_snv
        }
        
        if(!is.na(fp_snv$STB)){
            fp_snv =  fp_snv %>% filter(fp_snv$STB >= as.numeric(input$stb_snv))
        }else{
            fp_snv =  fp_snv
        }
        
        fp_snv
    })
    fn_tab_snv = reactive({
        req(input$query)
        req(input$folder)
        #fn_snv = read.delim(paste0(input$folder,"/FNs_snv.txt"), header = F, sep = " ")
        tmp = query_reac() %>% filter(Type == "SNV")
        gt_tmp = gt_reactive_snv() %>% filter(V16 == "SNV")
        if(!is.na(tmp$VAF)){
            tmp =  tmp %>% filter(tmp$VAF >= as.numeric(input$vaf_snv))
        }else{
            tmp =  tmp
        }
        
        if(!is.na(tmp$DP)){
            tmp =  tmp %>% filter(tmp$DP >= as.numeric(input$dp_snv))
        }else{
            tmp =  tmp
        }
        
        if(!is.na(tmp$QD)){
            tmp =  tmp %>% filter(tmp$QD >= as.numeric(input$qd_snv))
        }else{
            tmp =  tmp
        }
        
        if(!is.na(tmp$STB)){
            tmp =  tmp %>% filter(tmp$STB >= as.numeric(input$stb_snv))
        }else{
            tmp =  tmp
        }
        
        fn_snv = anti_join(gt_tmp, tmp, by = c("V1", "V2", "V3", "V4", "V5"))
        fn_snv = fn_snv[,c(1:5)]
        colnames(fn_snv) = c("Chr", "Start", "End", "Ref", "Alt")
        #colnames(fn_snv) = c("Chr", "Start", "End", "Ref", "Alt","Info", "Type", "VAF", "DP", "QD", "STB")
        fn_snv
    }) 
    
    output$tp_tab_snv = renderReactable({
        if(input$folder == "path/to/metrics/"){
            return(NULL)
        }else{
            reactable(tp_tab_snv())
        }
        
    })
    
    output$fp_tab_snv = renderReactable({
        if(input$folder == "path/to/metrics/"){
            return(NULL)
        }else{
            reactable(fp_tab_snv())
        }
    })
    
    output$fn_tab_snv = renderReactable({
        if(input$folder == "path/to/metrics/"){
            return(NULL)
        }else{
            reactable(fn_tab_snv())
        }
    })
    
##INDELs
    tp_tab_indel = reactive({
        req(input$query)
        req(input$folder)
        tp_indel = read.delim(paste0(input$folder,"/TPs_indel.txt"), header = F, sep = " ")
        tp_indel = inner_join(query_reac(), tp_indel, by = c("V1", "V2", "V3", "V4", "V5"))
        colnames(tp_indel) = c("Chr", "Start", "End", "Ref", "Alt","Info", "Type", "VAF", "DP", "QD", "STB")
        tp_indel$Info = NULL
        tp_indel$Type = NULL
        if(!is.na(tp_indel$VAF)){
            tp_indel =  tp_indel %>% filter(tp_indel$VAF >= as.numeric(input$vaf_indel))
        }else{
            tp_indel =  tp_indel
        }
        
        if(!is.na(tp_indel$DP)){
            tp_indel =  tp_indel %>% filter(tp_indel$DP >= as.numeric(input$dp_indel))
        }else{
            tp_indel =  tp_indel
        }
        
        if(!is.na(tp_indel$QD)){
            tp_indel =  tp_indel %>% filter(tp_indel$QD >= as.numeric(input$qd_indel))
        }else{
            tp_indel =  tp_indel
        }
        
        if(!is.na(tp_indel$STB)){
            tp_indel =  tp_indel %>% filter(tp_indel$STB >= as.numeric(input$stb_indel))
        }else{
            tp_indel =  tp_indel
        }
        
        tp_indel
    }) 
    
    fp_tab_indel = reactive({
        req(input$query)
        req(input$folder)
        fp_indel = read.delim(paste0(input$folder,"/FPs_indel.txt"), header = F, sep = " ")
        fp_indel = inner_join(query_reac(), fp_indel, by = c("V1", "V2", "V3", "V4", "V5"))
        colnames(fp_indel) = c("Chr", "Start", "End", "Ref", "Alt","Info", "Type", "VAF", "DP", "QD", "STB")
        fp_indel$Info = NULL
        fp_indel$Type = NULL
        if(!is.na(fp_indel$VAF)){
            fp_indel =  fp_indel %>% filter(fp_indel$VAF >= as.numeric(input$vaf_indel))
        }else{
            fp_indel =  fp_indel
        }
        
        if(!is.na(fp_indel$DP)){
            fp_indel =  fp_indel %>% filter(fp_indel$DP >= as.numeric(input$dp_indel))
        }else{
            fp_indel =  fp_indel
        }
        
        if(!is.na(fp_indel$QD)){
            fp_indel =  fp_indel %>% filter(fp_indel$QD >= as.numeric(input$qd_indel))
        }else{
            fp_indel =  fp_indel
        }
        
        if(!is.na(fp_indel$STB)){
            fp_indel =  fp_indel %>% filter(fp_indel$STB >= as.numeric(input$stb_indel))
        }else{
            fp_indel =  fp_indel
        }
        fp_indel
    })
    
    fn_tab_indel = reactive({
        req(input$query)
        req(input$folder)
        #fn_indel = read.delim(paste0(input$folder,"/FNs_indel.txt"), header = F, sep = " ")
        #colnames(fn_indel) = c("Chr", "Start", "End", "Ref", "Alt")
        tmp = query_reac() %>% filter(Type == "INDEL")
        gt_tmp = gt_reactive_indel() %>% filter(V16 == "INDEL")
        if(!is.na(tmp$VAF)){
            tmp =  tmp %>% filter(tmp$VAF >= as.numeric(input$vaf_indel))
        }else{
            tmp =  tmp
        }
        
        if(!is.na(tmp$DP)){
            tmp =  tmp %>% filter(tmp$DP >= as.numeric(input$dp_indel))
        }else{
            tmp =  tmp
        }
        
        if(!is.na(tmp$QD)){
            tmp =  tmp %>% filter(tmp$QD >= as.numeric(input$qd_indel))
        }else{
            tmp =  tmp
        }
        
        if(!is.na(tmp$STB)){
            tmp =  tmp %>% filter(tmp$STB >= as.numeric(input$stb_indel))
        }else{
            tmp =  tmp
        }
        fn_indel = anti_join(gt_tmp, tmp, by = c("V1", "V2", "V3", "V4", "V5"))
        fn_indel = fn_indel[,c(1:5)]
        colnames(fn_indel) = c("Chr", "Start", "End", "Ref", "Alt")
        fn_indel
    }) 
    
    output$tp_tab_indel = renderReactable({
        if(input$folder == "path/to/metrics/"){
            return(NULL)
        }else{
            reactable(tp_tab_indel())
        }
        
    })
    
    output$fp_tab_indel = renderReactable({
        if(input$folder == "path/to/metrics/"){
            return(NULL)
        }else{
            reactable(fp_tab_indel())
        }
    })
    
    output$fn_tab_indel = renderReactable({
        if(input$folder == "path/to/metrics/"){
            return(NULL)
        }else{
            reactable(fn_tab_indel())
        }
    })
    #diagnostic plots
    output$diag_plot_snv = renderPlot({
        tmp_tp = tp_tab_snv()
        tmp_fp = fp_tab_snv()
        tmp_tp$Type = rep("TP", dim(tmp_tp)[1])
        tmp_fp$Type = rep("FP", dim(tmp_fp)[1])
        tmp = rbind(tmp_tp, tmp_fp)
       p = grid.arrange(
            ggdensity(tmp, x = "VAF",
                      add = "mean", rug = TRUE,
                      color = "Type", fill= "Type",
                      palette = c("#FF9AA2", "#B5EAD7"))+theme(axis.text = element_text(size = 20), axis.title.y = element_text(size = 22), axis.title.x = element_text(size = 22), axis.text.y = element_text(size = 22),plot.title = element_text(size = 24), legend.text = element_text(size = 20), legend.title = element_text(size = 22)),
            
            ggdensity(tmp, x = "DP",
                      add = "mean", rug = TRUE,
                      color = "Type", fill= "Type",
                      palette = c("#FF9AA2", "#B5EAD7"))+theme(axis.text = element_text(size = 20), axis.title.y = element_text(size = 22), axis.title.x = element_text(size = 22), axis.text.y = element_text(size = 22),plot.title = element_text(size = 24), legend.text = element_text(size = 20), legend.title = element_text(size = 22)),
            
            
            ggdensity(tmp, x = "QD",
                      add = "mean", rug = TRUE,
                      color = "Type", fill= "Type",
                      palette = c("#FF9AA2", "#B5EAD7"))+theme(axis.text = element_text(size = 20), axis.title.y = element_text(size = 22), axis.title.x = element_text(size = 22), axis.text.y = element_text(size = 22),plot.title = element_text(size = 24), legend.text = element_text(size = 20), legend.title = element_text(size = 22)),
            
            ggdensity(tmp, x = "STB",
                      add = "mean", rug = TRUE,
                      color = "Type", fill= "Type",
                      palette = c("#FF9AA2", "#B5EAD7"))+theme(axis.text = element_text(size = 20), axis.title.y = element_text(size = 22), axis.title.x = element_text(size = 22), axis.text.y = element_text(size = 22),plot.title = element_text(size = 24), legend.text = element_text(size = 20), legend.title = element_text(size = 22)),
            
            
            ncol=2, nrow=2
            
        ) 
        print(p)
    })
    
    output$diag_plot_indel = renderPlot({
        tmp_tp = tp_tab_indel()
        tmp_fp = fp_tab_indel()
        tmp_tp$Type = rep("TP", dim(tmp_tp)[1])
        tmp_fp$Type = rep("FP", dim(tmp_fp)[1])
        tmp = rbind(tmp_tp, tmp_fp)
       p = grid.arrange(
            ggdensity(tmp, x = "VAF",
                      add = "mean", rug = TRUE,
                      color = "Type", fill= "Type",
                      palette = c("#FF9AA2", "#B5EAD7"))+theme(axis.text = element_text(size = 20), axis.title.y = element_text(size = 22), axis.title.x = element_text(size = 22), axis.text.y = element_text(size = 22),plot.title = element_text(size = 24), legend.text = element_text(size = 20), legend.title = element_text(size = 22)),
            
            ggdensity(tmp, x = "DP",
                      add = "mean", rug = TRUE,
                      color = "Type", fill= "Type",
                      palette = c("#FF9AA2", "#B5EAD7"))+theme(axis.text = element_text(size = 20), axis.title.y = element_text(size = 22), axis.title.x = element_text(size = 22), axis.text.y = element_text(size = 22),plot.title = element_text(size = 24), legend.text = element_text(size = 20), legend.title = element_text(size = 22)),
            
            
            ggdensity(tmp, x = "QD",
                      add = "mean", rug = TRUE,
                      color = "Type", fill= "Type",
                      palette = c("#FF9AA2", "#B5EAD7"))+theme(axis.text = element_text(size = 20), axis.title.y = element_text(size = 22), axis.title.x = element_text(size = 22), axis.text.y = element_text(size = 22),plot.title = element_text(size = 24), legend.text = element_text(size = 20), legend.title = element_text(size = 22)),
            
            ggdensity(tmp, x = "STB",
                      add = "mean", rug = TRUE,
                      color = "Type", fill= "Type",
                      palette = c("#FF9AA2", "#B5EAD7"))+theme(axis.text = element_text(size = 20), axis.title.y = element_text(size = 22), axis.title.x = element_text(size = 22), axis.text.y = element_text(size = 22),plot.title = element_text(size = 24), legend.text = element_text(size = 20), legend.title = element_text(size = 22)),
            
            
            ncol=2, nrow=2
            #top=grid::textGrob("Parameters distributions in TPs and FPs - INDELs",gp=grid::gpar(fontsize=20,font=3))
            
        ) 
       print(p)
    })
    
    ##Wilcoxon test on params
    output$wilcox_tab_snv = renderReactable({
        params = c("VAF", "DP", "QD", "STB")
        tmp = data.frame()
        param_tab = data.frame()
        for (param in params){
            tmp_tp = tp_tab_snv()
            tmp_fp = fp_tab_snv()
            if(is.na(tmp_tp[,param]) || is.na(tmp_fp[,param])){
                next()
            }else{
                w.t = wilcox.test(as.numeric(tmp_tp[,param]), as.numeric(tmp_fp[,param]), paired = F, alternative = "two.side")
                tmp = data.frame(Param=param, p_value = round(w.t$p.value,3))
                param_tab = rbind(param_tab, tmp) 
            }
            
        }
        reactable(param_tab)
    })
    
    output$wilcox_tab_indel = renderReactable({
        params = c("VAF", "DP", "QD", "STB")
        tmp = data.frame()
        param_tab = data.frame()
        for (param in params){
            tmp_tp = tp_tab_indel()
            tmp_fp = fp_tab_indel()
            if(is.na(tmp_tp[,param]) || is.na(tmp_fp[,param])){
                next()
            }else{
                w.t = wilcox.test(as.numeric(tmp_tp[,param]), as.numeric(tmp_fp[,param]), paired = F, alternative = "two.side")
                tmp = data.frame(Param=param, p_value = round(w.t$p.value,3))
                param_tab = rbind(param_tab, tmp) 
            }
        }
        reactable(param_tab)
    })
    
    #Download metric tabs (TP, FP, FN)
    ##SNV
    output$download_tp_snv<-downloadHandler(
        
        filename=function(){
            
            paste(paste(paste(paste(paste0("SNV_TP_vaf_dp_qd_stb", input$vaf_snv), input$dp_snv, sep = "_"),input$qd_snv, sep = "_" ),input$stb_snv, sep = "_"),".csv",sep="")
            
        },
        
        content=function(file){
            write.csv(tp_tab_snv(),file=file)
        }
        
    )
    
    output$download_fp_snv<-downloadHandler(
        
        filename=function(){
            
            paste(paste(paste(paste(paste0("SNV_FP__vaf_dp_qd_stb", input$vaf_snv), input$dp_snv, sep = "_"),input$qd_snv, sep = "_" ),input$stb_snv, sep = "_"),".csv",sep="")
            
        },
        
        content=function(file){
            write.csv(fp_tab_snv(),file=file)
        }
        
    )
    
    output$download_fn_snv<-downloadHandler(
        
        filename=function(){
            
            paste(paste(paste(paste(paste0("SNV_FN_vaf_dp_qd_stb", input$vaf_snv), input$dp_snv, sep = "_"),input$qd_snv, sep = "_" ),input$stb_snv, sep = "_"),".csv",sep="")
            
        },
        
        content=function(file){
            write.csv(fn_tab_snv(),file=file)
        }
        
    )
    
    ##INDELs
    output$download_tp_indel<-downloadHandler(
        
        filename=function(){
            
            paste(paste(paste(paste(paste0("INDEL_TP_vaf_dp_qd_stb", input$vaf_indel), input$dp_indel, sep = "_"),input$qd_indel, sep = "_" ),input$stb_indel, sep = "_"),".csv",sep="")
            
        },
        
        content=function(file){
            write.csv(tp_tab_indel(),file=file)
        }
        
    )
    
    output$download_fp_indel<-downloadHandler(
        
        filename=function(){
            
            paste(paste(paste(paste(paste0("INDEL_FP_vaf_dp_qd_stb", input$vaf_indel), input$dp_indel, sep = "_"),input$qd_indel, sep = "_" ),input$stb_indel, sep = "_"),".csv",sep="")
            
        },
        
        content=function(file){
            write.csv(fp_tab_indel(),file=file)
        }
        
    )
    
    output$download_fn_indel<-downloadHandler(
        
        filename=function(){
            
            paste(paste(paste(paste(paste0("INDEL_FN_vaf_dp_qd_stb", input$vaf_indel), input$dp_indel, sep = "_"),input$qd_indel, sep = "_" ),input$stb_indel, sep = "_"),".csv",sep="")
            
        },
        
        content=function(file){
            write.csv(fn_tab_indel(),file=file)
        }
        
    )
    
}

# Run the application 
shinyApp(ui = ui, server = server)
