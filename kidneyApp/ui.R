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
              p("Welcome to Kidney C1's Shiny App! This app is a Kidney rejection risk calculation tool for clinical researchers that aims to increase our understanding of the influence of different genes on Kidney graft acceptance or rejection. The tool demonstrates the influence that different genes have on T-Cell mediated rejection and Antibody-Mediated rejection by providing a prediction model that quantifies the risk of an individual's kidney rejection given their gene profile. This model will provide clinical researchers with the tools necessary to improve kidney acceptance rates, and maximise the efficient allocation of Kidneys."),
              h3("Why"),
              p("Kidney rejection occurs within 10-15 patients per 100 transplants. Since more than 90,000 patients are in need of a kidney and only 20,000 kidney transplants occur each year, the utility of a kidney rejection risk calculator for clinical research is considerable ('Rejection of a transplanted kidney', 2022). It will provide advantageous information that will facilitate the efficient allocation of kidneys to patients and minimise the risk of rejection. With the goal of maximising the greater wellbeing of individuals in need of kidney transplants, this tool will provide advancements in knowledge for clinical researchers."),
              h3("How to Use"),
              p("ABMR and TCMR Analysis"),
              p("Both the Antibody-Mediated Rejection Analysis and T-cell Mediated Rejection Analysis tabs have three different sections to navigate through."),
              tags$ul(
                      tags$li("Overview"),
                      tags$ul(
                              tags$li("This section provides a brief introduction to ABMR or TCMR and how to navigate through the tab"),
                              ),
                      tags$li("Annotation Analysis"),
                      tags$ul(
                              tags$li("This component includes a Network plot, Dot Plot, Tree plot and a KEGG enrichment analysis on each type of rejection"),
                      ),
                      tags$li("Prediction"),
                      tags$ul(
                              tags$li("Provides researchers a prediction of kidney acceptance/rejection based upon uploaded individual gene profiles or custom gene profiles which researchers can create using the provided sliders."),
                      ),
                      
              ),
              h3("Disclaimer"),
              p("We have utilised the numerous datasets in the creation of our model. They were all collected from different sources and compiled in tandem to create our tool. Within this report the data collected was from academic sources and utilised for educational purposes. These sources have been properly cited, and the details of each are listed below."),
              h4("Dataset"),
              p("Data Sets
The core utility of our datasets is the similarity between the two. Uniformity between these datasets is maximised because the data was collected by the same people using the same affymetrix platform.
"),
              
              h5("GSE36059"),
              p("The first dataset selected was from a comparative study conducted by the American Society of Nephrology that entails the ‘Disappearance of T Cell-Mediated Rejection Despite Continued Antibody-Mediated Rejection in Late Kidney Transplant Recipients’. The dataset provides a vast amount of data, with 54675 features from 411 samples and clearly defines rejection type.
"),
              h5("GSE48581"),
              p("The GSE48581 data shows great similarity to our other dataset and hence is highly valid in our model. This dataset was selected from a clinical trial conducted by The American Society of Transplantation and the American Society of Transplant Surgeons. This is the  INTERCOM study involving ‘Potential impact of microarray diagnosis of T cell-mediated rejection in kidney transplants’. INTERCOM 300 included biopsies from 264 patients."),
     ),
     
     tabPanel("ABMR Analysis", 
              
              tabsetPanel(id = "tabset",
                          
                          tabPanel(
                            title = "Overview",
                            h3("What is ABMR?"),
                            p("Antibody mediated rejection (ABMR)  is a more common form of rejection and presents a significant challenge for long-term graft survival. (Kojc & Haler, 2022) ABMR occurs when anti donor specific antibodies such as anti-HLA antibodies work against the transplanted kidney, and incite its rejection. Unlike TCMR, the risk of kidney rejection from ABMR is a long term concern that increases overtime. Due to advancements in medical technology,  TCMR is no longer the primary cause of graft failure, and hence ABMR now poses the largest threat to kidney rejection, accounting for graft failure in 64% of rejected transplants. ABMR peaks 5 years after transplantation, and therefore is more difficult to manage.(Kwon, H. et al, 2021)"),
                            h3("Annotation Analysis"),
                            tags$ul(
                                    tags$li("Network Plot - This plot depicts the biological states and processes of our gene profiles through the aggregation of gene data sets."),
                                    tags$li("Dot plot - Depicts the proportion of significant genes present out of the total count of genesis associated with the annotated molecular process"),
                                    tags$li("Tree Plot - The grouping of genes  portrays different families of the molecular process. The branches represent the similarity of gene thar are markers of multiple processes"),
                                    tags$li("KEGG Enrichment Analysis - Shows the physiological process of different gene pathways to reveal which molecules are being upregulated and down regulated. Furthermore, it generates a visual understanding of how rejection occurs in different pathways to aid clinical researchers understanding."),
                                    
                            ),
                            
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
                                          includeCSS("www/a.css"),
                                          mainPanel(
                                                  tags$div(class = "myclass",
                                                  img(id="one",src = "ABMR1.png", height = 300, width = 300),
                                                  img(id="two",src = "ABMR2.png", height = 300, width = 300),
                                                  img(id="three",src = "ABMR3.png", height = 300, width = 300),
                                                  img(id="four",src = "ABMR4.png", height = 300, width = 300)
                                                  )
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
                                                  tabsetPanel(type = "tabs",
                                                              tabPanel("knn", plotlyOutput("knn")),
                                                              tabPanel("log",plotlyOutput("log")),
                                                              tabPanel("rf",plotlyOutput("rf")),
                                                              tabPanel("svm",plotlyOutput("svm")),
                                                              tabPanel("decisiontree",plotlyOutput("tree"))
                                                              
                                                  )
                                                  
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
                            h3("What is TCMR?"),
                            p("The next rejection type explored through our app is Acute  T-Cell mediated rejection (TCMR). TCMR occurs when the instersitium is infiltrated by T cells and macrophages and is characterised by decreased urine output and proteinuria. (Halloran P. F., 2010). TCMR affects about 10-12% of patients and  is a less common cause of rejection in comparison to other conditions. This is because TCMR can be easily controlled as numerous medicinal methods of management have been introduced. TCMR develop as soon as a week after transplantation and its risk of rejection reduces over time (Halloran P. F., 2010). As such the risk of rejection due to TCMR peaks in the short term and slowly reduces over time."),
                            h3("Why TCMR?"),
                            p("Despite the minimal risk of acute kidney rejection arising from pure TCMR, it is still a vital component to analyse within this report. It is basic biological knowledge that when the immune system is fighting forreign cells, its first response is to mobilise the utility of T-cells. Subsequently, since T-cells are necessary for the activation of B-cells, ABMR cannot occur without the preexistence of TCMR. Accordingly, to account for this ideology, this app will also target the identification of genes that attribute TCMR, because it is a contributing factor to the inducement of AMBR. Furthermore, this can also assist clinical researchers in determining new methods to optimise the level of immunosuppressive therapy used  for each patient."),
                            h3("Annotation Analysis"),
                            tags$ul(
                                    tags$li("Network Plot - This plot depicts the biological states and processes of our gene profiles through the aggregation of gene data sets."),
                                    tags$li("Dot plot - Depicts the proportion of significant genes present out of the total count of genesis associated with the annotated molecular process"),
                                    tags$li("Tree Plot - The grouping of genes  portrays different families of the molecular process. The branches represent the similarity of gene thar are markers of multiple processes"),
                                    tags$li("KEGG Enrichment Analysis - Shows the physiological process of different gene pathways to reveal which molecules are being upregulated and down regulated. Furthermore, it generates a visual understanding of how rejection occurs in different pathways to aid clinical researchers understanding."),
                                    
                            ),
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
                                        mainPanel(
                                                tags$div(class = "myclass",
                                                img(id="1",src = "TCMR1.png", height = 300, width = 300),
                                                img(id="2",src = "TCMR2.png", height = 300, width = 300),
                                                img(id="3",src = "TCMR3.png", height = 300, width = 300),
                                                img(id="4",src = "TCMR4.png", height = 300, width = 300),
                                                img(id="5",src = "TCMR5.png", height = 300, width = 300)
                                                )
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
