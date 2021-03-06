list.of.packages <- c("shiny","ggplot2","DT","readr","shinythemes","readxl","plotly",
                      "forcats","crosstalk","heatmaply","knitr","shinydashboard","forcats")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos = "http://cran.us.r-project.org" )

lapply(list.of.packages, require, character.only = TRUE, warn.conflicts = FALSE,quietly=TRUE)

options(shiny.maxRequestSize = 100*1024^2,bitmapType='cairo') # increase the size of the file can be uploaded

header <- dashboardHeader(
  titleWidth = 300,
  title = "shiny-clone"
)

sidebar <- ## Sidebar content
  dashboardSidebar(
    sidebarMenu(id = "tabs",
                menuItem("Upload data", icon = icon("upload"), tabName = "upload"),
                menuItem("Lineplots by patient", icon = icon("chart-line"), tabName = "lineplots"),
                menuItem("Heatmaps by patient", icon = icon("th"), tabName = "heatmaps"),
                menuItem("Lineplots by gene", icon = icon("chart-line"), tabName = "lineplots_gene"),
                menuItem("Heatmaps by gene", icon = icon("th"), tabName = "heatmaps_gene"),
                menuItem("Cohort summary", icon = icon("th"), tabName = "c_summary")
    ),
    width = 200
  )

body <- dashboardBody(
  
  #########
  ## Upload
  #########
  
  tabItems(
    
    tabItem(
      tabName = "upload",
      box(title = "Explore the mutations in your cohort", width = 6, solidHeader = TRUE, status = "primary",
          collapsible = TRUE,
          h4("Purpose"),
          p("shiny-clone was developed to explore the mutations in a cohort. Even though it was developed to visualise multiple samples for a patient over time, this is not a requirement."),
          
          h4("Input requirement"),  
            p("It requires a ",code('.csv')," file as input which contains specific column names: ",
                     strong('PID')," for the patient ID; ",
                     strong('SampleName')," to identify a sample within a patient (set this the same as PID if only one sample per time point is availavle); ",
                     strong('Time')," defines the different times available (e.g. Time1, Time2 etc..); ",
                     strong('SYMBOL')," is the gene identifier; ",
                     strong('mutation_key')," uniquely identifies a mutation and can be defined as gene:chr:pos:ref_allele:alt_allele; ",
                     strong('mutation_det') ," usually set as the consequence of the mutations on the protein but can be any detail to annotate the mutation; ",
                     strong('VAF'), " is the variant allele frequency; ",
                     strong('variant_type')," to specify if the variant is a SNV, a CNV, an INDEL; ",
                     strong('Blast')," is thought for cancer samples, it is an estimate of the tumour content;",
              strong('Outcome') ," defines the overall outcome of a patient (e.g. Responder, Refractory etc..)"),
          
          p(strong('Blast')," and ", strong('Outcome')," shouldn't be input requirements and I will work to update the app."),
    
    h4("Get started"),  
    p("Once a ",code('.csv')," file with the columns defined above has been created, follow the steps on the rigth to get started.")),
      
      fluidRow(
        column(2, 
               
               shiny::wellPanel(shiny::fileInput("file", "1. Choose CSV File",
                                                 multiple = FALSE,
                                                 accept = c(".csv")),
                                
                                shiny::textOutput("needs")
               )
        ), 
        
        column(4,   shiny::wellPanel( 
          
          shiny::uiOutput("time_order"),
          shiny::textOutput("time_order_text")
        )
        ),
        
        column(2,   shiny::wellPanel( 
          
          #conditionalPanel(condition = "input.Run",
                  shiny::textOutput('ready')
          #)
        #)
        )
      ),
      
      #fluidRow(
      #  column(2,actionButton(inputId = "Run", label = h4("Create plots!")))
      #)
      #),
      
      fluidRow(
        column(12, DT::dataTableOutput("contents"))
      )
    )
    ),
    
    ####################
    ## Tab lineplots
    ####################
    
    tabItem(
      tabName = "lineplots",
      
      box(title = "Produce line plots by patient", width = 12, 
          solidHeader = TRUE, status = "primary",
          
          fluidRow(
            column(2, wellPanel(
              #conditionalPanel(condition = "input.Run",
                               shiny::uiOutput('patient_select')
              #)
            )
            ),
            
            column(8, 
                   #conditionalPanel(condition = "input.Run",
                                    plotlyOutput("scatter1",height = "600px",width = "1000px")
                   #)
            )
          ), 
          
          fluidRow(
            
            column(8, 
                   downloadButton("downloadDataPanel1", "Download Table"),
                   DT::dataTableOutput("table_by_pid")
                   )
            
          )
          
      )
    ),
    
    
    ####################
    ## Tab heatmaps
    ####################
    
    tabItem(
      tabName = "heatmaps",
      box(title = "Produce heatmaps by patient", width = 12, solidHeader = TRUE, status = "primary",
          
          
          tabPanel("Heatmaps by patient",
                   
                   fluidRow(
                     column(2, 
                            # conditionalPanel(condition = "input.Run",
                                             shiny::uiOutput('patient_select1')
                           # )
                     )
                   ),
                   
                   fluidRow(
                     
                     column(10,
                            #shiny::conditionalPanel(condition = "input.Run",
                                                       plotly::plotlyOutput("scatter3",height = "800px",width = "800px")
                     #)
                     ),
                     
                     column(10, DT::dataTableOutput("table_heat")
                     )
                   )
          )
      )
    ),
    
    
    ####################
    ## Tab lineplots by gene
    ####################
    
    tabItem(
      tabName = "lineplots_gene",
      
      box(title = "Produce line plots by gene", width = 12, 
          solidHeader = TRUE, status = "primary",
          
          fluidRow(
            column(2, wellPanel(
              shiny::uiOutput('gene_select_line')
            )
            ),
            
            column(8, 
                   #conditionalPanel(condition = "input.Run",
                   plotlyOutput("scatter_genes",height = "600px",width = "1000px")
                   #)
            )
          ), 
          
          fluidRow(
            
            column(8, 
                   downloadButton("downloadDataPanel2", "Download Table"),
                   DT::dataTableOutput("table_by_gene")
            )
            
          )
          
      )
    ),
    
    
    ####################
    ## Tab heatmaps by gene
    ####################
    
    tabItem(
      tabName = "heatmaps_gene",
      box(title = "Produce heatmaps by gene", width = 12, solidHeader = TRUE, status = "primary",
          
          
          tabPanel("Heatmaps by gene",
                   
                   #conditionalPanel(condition = "input.Run",
                                    
                        plotlyOutput("scatter4",height = "900px",width = "900px"),
                   #)
                   
                   
                   fluidRow(
                     column(4,
                            #conditionalPanel(condition = "input.Run",
                                               
                                               shiny::wellPanel(
                                                 shiny::uiOutput('gene_select'),
                                                 shiny::textOutput("explain_heat")
                                               )
                     #)
                     ),
                     
                     column(10,
                            #conditionalPanel(condition = "input.Run",
                                                DT::dataTableOutput("contents_gene")
                     #)
                     )
                   )
                   
          )
      )
    ),
    
    ####################
    ## summary cohort
    ####################
    
    tabItem(
      tabName = "c_summary",
      box(title = "Sumamry of mutations in the cohort", width = 12, solidHeader = TRUE, status = "primary",
          
          
          tabPanel("Overall mutations",
                   
                   #conditionalPanel(condition = "input.Run",
                                    
                                    #shiny::wellPanel(
                                    #  shiny::uiOutput('variable_select')
                                    #), 
                                    
                                    sidebarPanel(
                                      h4("Plot of the number of patients with a mutation on a gene. 
                                         Click on one dot to explore the gene behaviour across patients."),
                                      
                                      plotly::plotlyOutput("overallscatter",height = "900px",width = "600px"),
                                      width = 6),
                                    
                                    mainPanel(
                                      tabsetPanel(
                                        tabPanel("Mutation details", plotOutput("summarygene"),height = "600px",width = "500px"),
                                        tabPanel("Table", DT::dataTableOutput("summarytable"))
                                      ),
                                      width = 5)
                   #)
                   
          )
      )
    )
    
    
    ## end tabs
    
    
  )
)



server <- function(input, output, session) {
  
  # Load data
  # Until a file hasn't been given in input display message Upload data
  readInput <- shiny::reactive({
    
        validate(
          need(!is.null(input$file),"Upload data")
        )
        
        # validated
        
        inFile <- input$file
        
        inFile <- read.csv(inFile$datapath)
        
        # column requirements
        need_columns <- c("PID","SampleName","SYMBOL","mutation_key","mutation_det",
                          "VAF","variant_type","Blast","Outcome")
        
        missing <- need_columns[!(need_columns %in% colnames(inFile))]
        
        validate(
          need(sum(colnames(inFile) %in% need_columns) == length(need_columns), 
               paste0('Check Column names! Make sure the following are present: ',
                      paste0(missing,collapse=","))
          ))    
    
        
        return(inFile)
        
  })
  
  # If no data is uploaded then use the one provided whose existence I have already checked
  output$needs <- shiny::renderText({
  
    validate(
      need(readInput() != "Upload data",NULL)
    )
      
      need_columns <- c("PID","SampleName","SYMBOL","mutation_key","mutation_det",
                        "VAF","variant_type","Blast","Outcome")
      d <- data.frame(readInput())
      missing <- need_columns[!(need_columns %in% colnames(d))]
      
      validate(
        need(sum(colnames(d) %in% need_columns) == length(need_columns), 
             paste0('Check Column names! Make sure the following are present: ',
                    paste0(missing,collapse=","))
        ))
      
      "Input data looks good!\n\nNow go to step 2 and select order of time \n\nor import a new variant CSV file.\n\n"  
    
  })
  
  output$ready <- shiny::renderText({
    
    # If a file is provided do not provide error message
    validate(
      need(!is.null(input$file) & length(input$time) == length(unique(readInput()$Time)),"You forgot some steps, upload data or check that you selected all time points available.")
    )
    
    "You're ready to plot.\n
    Select lineplots or heatmaps from the side meanu."
    
  })
  
  
  # Choose columns to display in the first panel

  
  # Display table with all data
  output$contents <- DT::renderDataTable({
    
    # If a file is provided do not provide error message
    validate(
      need(!is.null(input$file),NULL)
    )
    
    need_columns <- c("PID","SampleName","SYMBOL","mutation_key","mutation_det",
                      "VAF","variant_type","Blast","Outcome")
    d <- data.frame(readInput())
    missing <- need_columns[!(need_columns %in% colnames(d))]
    
    validate(
      need(sum(colnames(d) %in% need_columns) == length(need_columns), NULL))
    
    # Validated
    
    d <- d %>%
      dplyr::mutate(VAF = round(VAF,2)) %>%
      dplyr::select(PID,SampleName,PID,Time,Blast,Outcome) %>%
      unique()
    
    DT::datatable(d,filter = "top",rownames = FALSE,width = 20,
                  options = list(scrollX = TRUE,pageLength = 30))
  },server = TRUE)
  
  
  ##############################
  # Panel1: Lineplots by patient
  ##############################
  
  # Choose PID reactively
  output$patient_select <- renderUI({
    
    validate(
      need(!is.null(input$file) & length(input$time) == length(unique(readInput()$Time)),"You forgot some steps! Go back to `Upload data`")
    )
    
    d <- readInput()
    selectInput('PID','Patients',choices =  unique(d$PID),selected = unique(d$PID)[1])
  })
  
  # Choose time order reactively - in upload page
  output$time_order <- renderUI({
    
    d <- readInput()
    selectInput('time','2. Select order of time to use for plotting',
                choices =  unique(d$Time),
                selected = NULL,
                multiple = TRUE)
  })
  
  # Display selected time
  output$time_order_text <- renderText({
    
    validate(
      need(!is.null(input$time),NULL)
    )
    
    
    time_order <- paste0(input$time,collapse=",")
    # d <- data.frame(readInput()) 
    
    return(time_order)
    
  })
  
  # Subset data for one patient and order time points
  selectInputPID <- reactive({
    
    # If a file is provided do not provide error message
    validate(
      need(!is.null(input$file),NULL)
    )
    
    need_columns <- c("PID","SampleName","Time","SYMBOL","mutation_key","mutation_det","VAF","variant_type","Blast")
    d <- data.frame(readInput())
    missing <- need_columns[!(need_columns %in% colnames(d))]
    
    validate(
      need(sum(colnames(d) %in% need_columns) == length(need_columns), 
           paste0('Check Column names! Make sure the following are present: ',
                  paste0(missing,collapse=","))
      ))
    
    # Validate
    
    d <- data.frame(readInput())  %>% 
      dplyr::mutate(Time = factor(Time,levels=input$time)) %>%
      dplyr::mutate(SampleName = forcats::fct_reorder(SampleName,as.numeric(Time)))
    
    #d <- data.frame(readInput())
    
    pid_mut <- subset(d , PID %in% input$PID) %>%
      dplyr::mutate(VAF = round(VAF,2)) %>%
      dplyr::mutate(SYMBOL = ifelse(variant_type %in% "CNA", mutation_det,SYMBOL))
    
    pid_mut <- pid_mut[!duplicated(pid_mut),]
    
    pid_mut
    
  })
  
  # Create shared data
  shared_ti <- SharedData$new(selectInputPID,~mutation_det)
  
  output$scatter1 <- renderPlotly({
    
    validate(
      need(!is.null(input$PID),"Loading data...")
    )
    
    g1=ggplot(shared_ti,aes(x = SampleName,y = VAF, group = mutation_key,
                            text = mutation_det,colour=variant_type)) + 
      geom_point() + geom_line() + theme_bw() + 
      ggtitle(paste0("Mutations for patient ",input$PID)) +
      theme(legend.position = "bottom") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    
    gg1 = ggplotly(g1,tooltip = "text")
    highlight(gg1,selected = "text",on = "plotly_click",off = "plotly_doubleclick")
  })
  
  
  output$table_by_pid <- DT::renderDataTable({
    
    DT::datatable(shared_ti,rownames = FALSE,width = 20,
                  filter = "top",options = list(scrollX = TRUE,pageLength = 30))
  },server = FALSE)
  
  output$downloadDataPanel1 <- downloadHandler(
    filename = function() {
      paste0("mutations_patient_",input$PID, ".csv")
    },
    content = function(file) {
      write.csv(selectInputPID(), file, row.names = FALSE)
    }
  )
  
  
  ##############################
  # Panel2: Heatmap by patient
  ##############################
  
  # Choose PID reactively
  output$patient_select1 <- renderUI({
    
    validate(
      need(!is.null(input$file) & length(input$time) == length(unique(readInput()$Time)),
           "You forgot some steps! Go back to `Upload data`")
    )
    
    d <- readInput()
    selectInput('PID1','Patients',choices =  unique(d$PID),selected = unique(d$PID)[1])
  })
  
  # Subset data for one patient
  heatmapdata <- reactive({
    
    # Re-order times
    d <- data.frame(readInput())  %>% 
      dplyr::mutate(Time = factor(Time,levels=input$time)) %>%
      dplyr::mutate(SampleName = forcats::fct_reorder(SampleName,as.numeric(Time)))
    
    pid_mut <- d %>%
      dplyr::filter(PID %in% input$PID1)
    
    pid_mut
    
  })
  
  
  output$scatter3 <- renderPlotly({
    
    validate(
      need(!is.null(input$PID1),"Loading data...")
    )
    
    # Already filtered by patient but contains many CNA depending on how many genes are in that CNA
    data <- heatmapdata()
    data <- data %>% 
      dplyr::mutate(SYMBOL = ifelse(variant_type %in% "CNA", mutation_det,SYMBOL)) %>%
      dplyr::mutate(mutation_det = ifelse(variant_type %in% "CNA", paste(mutation_det, chrom, sep = " "),as.character(mutation_det)))
    
    # Keep only one CNA per patient
    data <- data[!duplicated(data),] %>%
      dplyr::select(mutation_key,mutation_det,variant_type,SampleName,VAF) %>%
      dplyr::mutate(VAF = round(VAF,2)) %>%
      tidyr::spread(SampleName,value = VAF) 
    
    heatmaply(data[,4:ncol(data)],row_dend_left=FALSE,limits = c(0,1),colors = BuPu,
              margins = c(60,400,40,200),Rowv=TRUE,Colv = FALSE,
              labRow = data$mutation_det,
              row_side_colors = data$variant_type,
              main = paste0("Mutations for patient ",input$PID1))
    
  })
  
  output$table_heat <- DT::renderDataTable({
    data <- heatmapdata() %>%
      dplyr::mutate(VAF = round(VAF,2)) %>%
      dplyr::select(SampleName, Time,chrom,pos,mutation_key,mutation_det,SYMBOL,VAF,
                    variant_type,tot_depth,SYMBOL)
    datatable(data,
              filter = "top",options = list(pageLength = 30,scrollX = TRUE))
  },server = TRUE)
  
  
  ##########################
  ## Lineplot by gene
  ##########################
  
  # Subset data for one gene and order time points
  
  # Choose symbol/CNA reactively
  output$gene_select_line <- shiny::renderUI({
    
    d <- readInput()
    
    d <- d %>% 
      dplyr::filter(!(variant_type %in% "clones"))

    # Create gene selection
    shiny::selectInput('gene_line','Select gene',
                       choices =  unique(d$SYMBOL),
                       selected = unique(d$SYMBOL)[1],multiple=TRUE)
    
  })
  
  
  
  lineplots_gene <- reactive({
    
    # If a file is provided do not provide error message
    validate(
      need(!is.null(input$file),NULL)
    )
    
    need_columns <- c("PID","SampleName","SYMBOL","mutation_key","mutation_det","VAF","variant_type","Blast")
    d <- data.frame(readInput())
    missing <- need_columns[!(need_columns %in% colnames(d))]
    
    validate(
      need(sum(colnames(d) %in% need_columns) == length(need_columns), 
           paste0('Check Column names! Make sure the following are present: ',
                  paste0(missing,collapse=","))
      ))
    
    # Validated
    
    d <- data.frame(readInput())  %>% 
      dplyr::mutate(Time = factor(Time,levels=input$time)) %>%
      dplyr::mutate(SampleName = forcats::fct_reorder(SampleName,as.numeric(Time)))
    
    gene_mut_line <- d %>%
      dplyr::filter(SYMBOL %in% input$gene_line) %>%
      dplyr::group_by(PID,Time,SYMBOL,mutation_det,mutation_key,variant_type) %>%
      dplyr::mutate(VAF = round(mean(VAF),2),
                    Outcome = as.character(Outcome)) %>%
      ungroup() %>%
      dplyr::mutate(SYMBOL = ifelse(variant_type %in% "CNA", mutation_det,SYMBOL),
                    mutation_det_id = paste(mutation_det,PID,sep="-"))
    
    return(gene_mut_line)
    
  })

  
  # Create shared data
  shared_ti_gene <- SharedData$new(lineplots_gene,~mutation_det_id)
  
  output$scatter_genes <- renderPlotly({
    
    validate(
      need(!is.null(input$file) & length(input$time) == length(unique(readInput()$Time)),"You forgot some steps! Go back to `Upload data`")
    )
    
    validate(
      need(!is.null(input$gene_line) & !is.null(input$file) & length(input$time) == length(unique(readInput()$Time)),"Loading data...")
    )
    
  
    g1 = ggplot(shared_ti_gene,aes(x = Time,y = VAF, group = interaction(PID,mutation_key),
                            text = mutation_det_id,colour = Outcome)) + 
      geom_point() + geom_line() + theme_bw() + 
      ggtitle(paste0("Mutations found on gene ",input$gene_line)) +
      theme(legend.position = "bottom") + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    
    gg1 = ggplotly(g1,tooltip = "text")
    highlight(gg1,selected = "text",on = "plotly_click",off = "plotly_doubleclick")
  })
  
  
  output$table_by_gene <- DT::renderDataTable({
    
    DT::datatable(shared_ti_gene,rownames = FALSE,width = 20,
                  filter = "top",options = list(scrollX = TRUE,pageLength = 30))
  },server = FALSE)
  
  output$downloadDataPanel2 <- downloadHandler(
    filename = function() {
      paste0("mutations_genes_",input$gene_line, ".csv")
    },
    content = function(file) {
      write.csv(lineplots_gene(), file, row.names = FALSE)
    }
  )
  
  
  
  
  ##############################
  # Panel3: Heatmap by gene
  ##############################
  # I am only plotting CNA and genes not clones
  # I am computing the mean VAF for a patient at a time point, including tissue replicates
  # Because of how CNas are reported, I attach chromosome name to the CNA
  
  output$explain_heat <- shiny::renderText({
    
    "For every variant, the mean VAF is computed for a patient at a specific time point. This means that if multiple samples are available for a patient at a time point (either tissue or technical replicates), then their mean VAF will be computed."
    
  })
  
  
  # Choose symbol/CNA reactively
  output$gene_select <- shiny::renderUI({
    
    d <- readInput()
    
    d <- d %>% 
      dplyr::filter(!(variant_type %in% "clones"))
    #  dplyr::mutate(mutation_det = ifelse(variant_type %in% "CNA",
    #                                      paste(mutation_det,chrom),
    #                                      mutation_det)) %>%
    #  dplyr::filter(!(variant_type %in% "clones")) %>%
    #  dplyr::mutate(SYMBOL = ifelse(variant_type %in% "CNA", mutation_det,SYMBOL))
    
    # Create gene selection
    shiny::selectInput('gene','Select mutation',
                       choices =  unique(d$SYMBOL),
                       selected = unique(d$SYMBOL)[1],multiple=TRUE)
    
  })
  
  
  output$show_genes <- shiny::renderText({
    
    return(paste0(input$gene,collapse=","))
    
  })
  
  
  # Subset data for one patient
  # Paste chromosome to CNA spec; filter out clones
  heatmapdata_gene <- reactive({
    
    d <- readInput() %>% 
      dplyr::filter(!(variant_type %in% "clones"))
    
    gene_mut <- d %>%
      dplyr::filter(SYMBOL %in% input$gene)
    
    return(gene_mut)
    
  })
  
  # tables with selected genes
  output$contents_gene <- DT::renderDataTable({
    
    data <- heatmapdata_gene()
    
    DT::datatable(data,filter = "top",options = list(scrollX = TRUE,pageLength = 30))
  },server = TRUE)
  
  
  # Need to improve - labeling of clones/labeling of mutations
  # Probably all CNV on chr 9 should go together
  
  output$scatter4 <- renderPlotly({
    
    validate(
      need(!is.null(input$file) & length(input$time) == length(unique(readInput()$Time)),"You forgot some steps! Go back to `Upload data`")
    )
    
    validate(
      need(!is.null(input$gene) & !is.null(input$file) & length(input$time) == length(unique(readInput()$Time)),"Loading data...")
    )
    
    data <- heatmapdata_gene()
    
    data <- data %>%
      dplyr::group_by(PID,Time,SYMBOL,mutation_key,mutation_det,variant_type,Outcome) %>%
      dplyr::summarise(VAF = round(mean(VAF),2)) %>%
      dplyr::select(mutation_key,SYMBOL,mutation_det,PID,variant_type,Time,VAF,Outcome) %>%
      ungroup() %>%
      dplyr::mutate(Time = factor(Time,levels = input$time),
                    Outcome = factor(Outcome)) %>%
      dplyr::mutate(mutation_det = paste(PID, mutation_key)) %>%
      tidyr::spread(Time,value = VAF) 
    
    row_den <- ifelse(nrow(data) > 1,TRUE,FALSE)
    
    heatmaply(data[,5:ncol(data)],row_dend_left=FALSE,limits = c(0,1),colors = BuPu,
              margins = c(60,400,40,200),Rowv=FALSE,Colv = FALSE,
              labRow = data$mutation_det,
              main = paste0("Mutations"))
    
  })
  
  
  #################################
  # Panel3: Plots/Table by mutation by variable
  #################################
  
  # Choose PID reactively
  output$variable_select <- renderUI({
    
    validate(
      need(!is.null(input$file) & length(input$time) == length(unique(readInput()$Time)),"You forgot some steps! Go back to `Upload data`")
    )
    
    
    d <- readInput() %>%
      dplyr::select(-mutation_key,-mutation_det,-SYMBOL,-PID,-SampleName,-VAF)
   
    varSelectInput("var", "Select variable:", d,selected = NULL)
  })
  
  
  variants_sum <- shiny::reactive({
    
    validate(
      need(!is.null(input$file) & length(input$time) == length(unique(readInput()$Time)),"You forgot some steps! Go back to `Upload data`")
    )
    
    
    d <- readInput() %>%
      dplyr::mutate(Outcome = as.character(Outcome))

    # Number of patients with that mutation
    d %>%
      dplyr::group_by(SYMBOL,Outcome) %>%
      dplyr::summarise(n.patients = length(unique(PID))) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(SYMBOL = forcats::fct_reorder(SYMBOL,n.patients))
    
  })
  
  
  output$overallscatter <- renderPlotly({
    
    d <- variants_sum()
    
    # Plot number of patients with a specific gene mutated
    g1 <- ggplot(d,aes(x=SYMBOL,y=n.patients,key=SYMBOL))  +
      geom_jitter(position = position_jitter(w = 0, h = 0.1),size=1.8,alpha=0.8,aes(colour = Outcome)) + 
      theme_bw() + coord_flip() + 
      ggtitle("Number of overall mutations found on genes") + scale_color_brewer(palette="Set2") +
      theme(legend.position = "bottom")
    
    ggplotly(g1,source = "allmutplot")
    
    
  })
  
  
  output$selection <- renderPrint({
    s <- event_data("plotly_click", source = "allmutplot")
    if (length(s) == 0) {
      "Click on a dot in the scatter plot to visualise the mutations on this gene across patients"
    } else {
      cat("You selected: \n\n")
      as.list(s)
    }
  })
  
  output$click <- renderPrint({
    d <- event_data("plotly_click", source = "allmutplot")
    d
  })
  
  
  output$summarygene <- renderPlot({
    
    s <- event_data("plotly_click", source = "allmutplot")
    
    if (length(s)) {
      
      symbol <- s$key
      
      d <- readInput()
      
      sel <- d[d$SYMBOL %in% symbol,] 
      
      sel_data <- sel %>% 
        dplyr::group_by(mutation_det,Outcome) %>%
        dplyr::summarise(n.patients = length(unique(PID))) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(mutation_det = forcats::fct_reorder(mutation_det,n.patients)) %>%
        dplyr::mutate(Outcome = as.character(Outcome)) 
      
      
        ggplot(sel_data, aes(x = mutation_det,y = n.patients,colour=Outcome)) + 
        geom_jitter(position = position_jitter(w = 0, h = 0.1),alpha=0.8,aes(colour=Outcome),size=3) + theme_bw() + 
        coord_flip() + ggtitle(paste0("Mutation types on ",symbol)) + scale_color_brewer(palette="Set2") + 
        scale_y_continuous(breaks = seq(1,max(sel_data$n.patients),1),limits = c(0,max(sel_data$n.patients))) +
        labs(x="Mutation details") 
      
      
    } else {
      NULL
    }
  })
  
  output$summarytable <- DT::renderDataTable({
    
    d <- readInput()
    
    s <- event_data("plotly_click", source = "allmutplot")
    if (length(s)) {
      symbol <- s$key
      
      sel <- d[d$SYMBOL %in% symbol,] %>%
        dplyr::select(SampleName,mutation_det,VAF,mutation_key,variant_type) %>%
        dplyr::mutate(VAF = round(VAF,2))
      
      datatable(sel,filter = "top")
      
    } else {
      NULL
    }
  })
  

  
  
  
  
}

ui <- dashboardPage(header, sidebar, body, skin = "black",
                    title = "Explore mutations")

shinyApp(ui, server)
