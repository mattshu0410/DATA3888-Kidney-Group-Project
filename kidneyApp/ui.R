library(shiny)
library(shinythemes)
library(shinyjs)
library(DT)
library(shinyalert)

# Define UI for application that draws a histogram
shinyUI(navbarPage(theme = shinytheme("flatly"), 
                   title = "Kidney C1", 
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
                                title = "Enrichment Analysis",
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
                          )
              )
     )
     
))))
