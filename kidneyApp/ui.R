library(shiny)
library(shinythemes)
library(shinyjs)
library(DT)
library(shinyalert)

# Define UI for application that draws a histogram
shinyUI(navbarPage(tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "mystyle.css")
),
                   theme = shinytheme("flatly"),
                   title = img(src = "logo1.jpg",inline=TRUE, id="logo"),
                   fluid=TRUE, 

     useShinyjs(),
     
     tabPanel("Home",
              
     ),
     
     tabPanel("ABMR Analysis", 
              
              tabsetPanel(id = "tabset",
                          
                          tabPanel(
                            title = "Overview",
                          ),
                          
                          tabPanel(
                            title = "Annotation Analysis",
                            
                            tabsetPanel(
                                        tabPanel(
                                          title = "Network Plot",
                                          
                                          plotOutput("networkPlotABMR")
                                          
                                        ),
                                        
                                        tabPanel(
                                          title = "Dot Plot",
                                          
                                          plotOutput("dotPlotABMR")
                                          
                                        ),
                                        
                                        tabPanel(
                                          title = "Tree Plot",
                                          
                                          plotOutput("treePlotABMR")
                                        ),
                                        
                                        tabPanel(
                                          title = "KEGG Enrichment Analysis",
                                          tags$head(
                                                  tags$link(rel = "stylesheet", type = "text/css", href = "a.css")
                                          ),
                                          mainPanel(
                                                  img(id="one",src = "ABMR1.png", height = 300, width = 300),
                                                  img(id="two",src = "ABMR2.png", height = 300, width = 300),
                                                  img(id="three",src = "ABMR3.png", height = 300, width = 300),
                                                  img(id="four",src = "ABMR4.png", height = 300, width = 300)
                                          )
                                        )
                            )
                          ),
                          
                          tabPanel(
                                  title = "Prediction",
                                  
                                  tabsetPanel(id = "tabset4",
                                              
                                  tabPanel(
                                  title= "File Input",
                                  sidebarLayout(
                                          
                                          sidebarPanel(
                                                  
                                                  titlePanel("Uploading Files"),
                                                  
                                                  fileInput('target_upload1', 'Choose file to upload',
                                                            accept = c(
                                                                    'text/csv',
                                                                    'text/comma-separated-values',
                                                                    '.csv'
                                                            )),
                                                  tags$hr(),
                                                  radioButtons("dis1", "Display",
                                                               choices = c(Head = "head",
                                                                           All = "all"),
                                                               selected = "head"),
                                                  tags$hr(),
                                                  actionButton("showTab", "Show Prediction"),
                                          ),
                                          mainPanel(
                                                  DTOutput("ab"),
                                                  
                                          )
                                  )
                          ),
                          tabPanel(
                                  title="Manual Input",
                          )
                          
              )
     )
     )
     ),

     tabPanel("TCMR Analysis", 
              
              tabsetPanel(id = "tabset2",
                          
                          tabPanel(
                            title = "Overview",
                          ),
                          
                          tabPanel(
                            title = "Annotation Analysis",
                            
                            tabsetPanel(
                              tabPanel(
                                title = "Network Plot",
                                
                                plotOutput("networkPlotTCMR")
                              ),
                              
                              tabPanel(
                                title = "Dot Plot",
                                
                                plotOutput("dotPlotTCMR")
                              ),
                              
                              tabPanel(
                                title = "Tree Plot",
                                
                                plotOutput("treePlotTCMR")
                              ),
                              
                              tabPanel(
                                title = "KEGG Enrichment Analysis",
                                tags$head(
                                        tags$link(rel = "stylesheet", type = "text/css", href = "a.css")
                                ),
                                        mainPanel(
                                                img(id="1",src = "TCMR1.png", height = 300, width = 300),
                                                img(id="2",src = "TCMR2.png", height = 300, width = 300),
                                                img(id="3",src = "TCMR3.png", height = 300, width = 300),
                                                img(id="4",src = "TCMR4.png", height = 300, width = 300),
                                                img(id="5",src = "TCMR5.png", height = 300, width = 300)
                                        )
                              )
                            )
                          ),
                          
                          tabPanel(
                            title = "Prediction",
                            tabsetPanel(id = "tabset3",
                            
                            tabPanel(
                            title= "File Input",
                            sidebarLayout(
                                    
                                    sidebarPanel(
                                            
                            titlePanel("Uploading Files"),
                            
                            fileInput('target_upload', 'Choose file to upload',
                                      accept = c(
                                              'text/csv',
                                              'text/comma-separated-values',
                                              '.csv'
                                      )),
                            tags$hr(),
                            radioButtons("dis", "Display",
                                         choices = c(Head = "head",
                                                     All = "all"),
                                         selected = "head"),
                            tags$hr(),
                            actionButton("showTab", "Show Prediction"),
                                    ),
                            mainPanel(
                                    DTOutput("abc"),
                                    
                            )
                            )
                          ),
                          tabPanel(
                          title="Manual Input",
                          sidebarLayout(
                                  
                                  sidebarPanel(
                                          uiOutput("sliders"),
                                  ),
                                  
                                  # Main panel for displaying outputs ----
                                  mainPanel(
                                          
                                          # Output: Table summarizing the values entered ----
                                          tableOutput("values"),
                                          
                                  )
                          )
              )
     )
     
)))))
