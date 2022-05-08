library(shiny)
library(shinythemes)
library(shinyjs)
library(tidyverse)
library(DT)

# Define UI for application that draws a histogram
ui <- navbarPage(theme = shinytheme("flatly"), title = "Kidney C1", fluid=TRUE, id = "tabs",
    useShinyjs(),
    
    tabPanel("Home",

    ),

    tabPanel("Details on Risk Calculation", 
             
        tabsetPanel(id = "tabset",

            tabPanel(
              title = "Model Weaknesses",
            ),
            
            tabPanel(
              title = "Consideration of Fairness",
            ),
            
            tabPanel(
              title = "Advanced Details",
            )
        )
    ),

    # file upload page
    tabPanel("File Upload",
        
        # accepts a csv file
        fileInput('target_upload', 'Choose file to upload',
                  accept = c(
                    'text/csv',
                    'text/comma-separated-values',
                    '.csv'
                  )),
        
        # IF a csv file uploaded, printing out file contents for initial inspect
        DTOutput('tbl'),
        
        # IF a csv file uploaded, creating button that brings user to prediction page
        actionButton("showTab", "Show Prediction"),
        
    )
    
)
                 
                 
# Define server logic required to draw a histogram
server <- function(input, output) {
  
   # if action button clicked, insert a new prediction tab and hide action button
   # so duplicate tabs aren't added
   observeEvent(input$showTab, {
     
       # select = TRUE immediately brings user to prediction page
       insertTab(inputId = "tabs",
                 tabPanel("Prediction", position ="after"), select = TRUE)
        shinyjs::hide("showTab")
   })
  
  # hide the action button initially
  observe({
    shinyjs::hide("showTab")
  })
  
  # if csv added, show action button
  observeEvent(input$target_upload, {
    shinyjs::show("showTab")
  })
  
  # reading csv file and returning data frame
  # DO CHECKS AND IMPUTATION HERE
  mydata <- reactive({
    
    # string file name
    inFile <- input$target_upload
    
    if (is.null(inFile))
      return(NULL)
    
    # reading file
    tbl <- read.csv(inFile$datapath)
      
    if(ncol(tbl)!=2){
      shinyalert("Oops!", "Column is missing", type = "error")
    }else{
    for(i in tbl[1]){
    if(!is.character(i))
      shinyalert("Oops!", "The first column should be name of genes", type = "error")
      break
    }
    for(i in tbl[2]){
    if(!is.numeric(i) || !is.double(i))
      shinyalert("Oops!", "The second column should be numeric expression set value", type = "error")
      break
    }
    }
      
    
    return (tbl)
  })
  
  # rendering table in File Upload page
  output$tbl = DT::renderDT({
    mydata()
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
                 
