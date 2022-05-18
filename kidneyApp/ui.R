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
              h3("Introduction"),
              p("Welcome to Kidney C1’s Shiny App! This app is a Kidney rejection risk calculation tool for clinical researchers that aims to increase our understanding of the influence of different genes on Kidney graft acceptance or rejection. The tool demonstrates the influence that different genes have on T-Cell mediated rejection and Antibody-Mediated rejection by providing a prediction model that quantifies the risk of an individual's kidney rejection given their gene profile. This model will provide clinical researchers with the tools necessary to improve kidney acceptance rates, and maximise the efficient allocation of Kidneys."),
              br(),
              h3("Why"),
              p("Kidney rejection occurs within 10-15 patients per 100 transplants. Since more than 90,000 patients are in need of a kidney and only 20,000 kidney transplants occur each year, the utility of a kidney rejection risk calculator for clinical research is considerable ('Rejection of a transplanted kidney', 2022). It will provide advantageous information that will facilitate the efficient allocation of kidneys to patients and minimise the risk of rejection. With the goal of maximising the greater wellbeing of individuals in need of kidney transplants, this tool will provide advancements in knowledge for clinical researchers."),
              br(),
              h3("How to Use"),
              p("ABMR and TCMR Analysis"),
              p("Both the ‘Antibody-Mediated Rejection’ Analysis and ‘T-cell Mediated Rejection Analysis’ tabs have three different sections to navigate through."),
              
              p("Overview"),
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
                                                  actionButton("showTab", "Show Prediction"),
                                          ),
                                          mainPanel(
                                                  plotlyOutput("pcaplot1"),
                                                  
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
                            actionButton("showTab", "Show Prediction"),
                                    ),
                            mainPanel(
                                    #DTOutput("abc"),
                                    
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
