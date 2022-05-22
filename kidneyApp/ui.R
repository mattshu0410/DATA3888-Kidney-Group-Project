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
              h2("Introduction"),
              p("Welcome to Kidney C1's Shiny app! This tool has been built for clinical researchers to demonstrate the influence that different genes have on T-Cell Mediated Rejection and Antibody-Mediated Rejection, by providing a predictive model that quantifies the risk of a kidney rejection given an individual's gene profile. It also allows researchers to explore significant pathways and gene clusters. Inevitably, this model will provide clinical researchers with the tools necessary to improve kidney acceptance rates and post-transplant quality of life."),
              h3("How to Use"),
              h4("ABMR and TCMR Analysis"),
              p("Both the Antibody-Mediated Rejection Analysis and T-Cell Mediated Rejection Analysis tabs have three different sections to navigate through."),
              tags$ul(
                      tags$li("Overview"),
                      tags$ul(
                              tags$li("This section provides a brief introduction to ABMR or TCMR and how to navigate through the tab."),
                              ),
                      tags$li("Annotation Analysis"),
                      tags$ul(
                              tags$li("This component includes a Network plot, Dot plot, Tree plot and KEGG enrichment analysis for each type of rejection."),
                      ),
                      tags$li("Prediction"),
                      tags$ul(
                              tags$li("Provides researchers a prediction of kidney acceptance/rejection based upon uploaded individual gene profiles, or custom gene profiles created using the provided sliders."),
                      ),
                      
              ),
              h3("Datasets Used"),
              h5("GSE36059"),
              p("The first dataset selected was from a comparative study conducted by the American Society of Nephrology that entails the 'Disappearance of T Cell-Mediated Rejection Despite Continued Antibody-Mediated Rejection in Late Kidney Transplant Recipients'. The dataset provides a vast amount of data, with 54675 features from 411 samples, and clearly defines rejection types."),
              h5("GSE48581"),
              p("The GSE48581 data was created as part of a clinical trial conducted by The American Society of Transplantation and the American Society of Transplant Surgeons. This is an INTERCOM study involving 'Potential impact of microarray diagnosis of T cell-mediated rejection in kidney transplants'. It includes biopsies from 264 patients."),
              h3("Disclaimer"),
              p("We have utilised numerous datasets in the creation of our model that were collected from different sources, and compiled in tandem to create our tool. The data was utilised for educational purposes, and the app should only be used by clinical researchers. It accepts no liability for the quality or accuracy of its predictions/information.")
     ),
     
     tabPanel("Antibody-Mediated Rejection Analysis", 
              
              tabsetPanel(id = "tabset",
                          
                          tabPanel(
                            title = "Overview",
                            h3("What is ABMR?"),
                            p("Antibody-Mediated Rejection (ABMR) is a more common form of rejection and presents a significant challenge for long-term graft survival (Kojc & Haler, 2022). ABMR occurs when anti-donor specific antibodies such as anti-HLA antibodies work against the transplanted kidney, and incite its rejection. Unlike TCMR, the risk of kidney rejection from ABMR is a long term concern that increases overtime. Due to advancements in medical technology, ABMR now poses the largest threat to kidney rejection, accounting for graft failure in 64% of rejected transplants. ABMR peaks 5 years after transplantation, and therefore is more difficult to manage.(Kwon, H. et al, 2021)"),
                            h3("Annotation Analysis"),
                            tags$ul(
                                    tags$li("Network Plot - This plot depicts the biological states and processes of gene profiles through the aggregation of gene data sets."),
                                    tags$li("Dot plot - Depicts the proportion of significant genes present out of the total count of genes that are associated with the annotated molecular process."),
                                    tags$li("Tree Plot - The grouping of genes portrays different families of the molecular process. The branches represent the similarity of genes that are markers of multiple processes."),
                                    tags$li("KEGG Enrichment Analysis - Shows the physiological process of different gene pathways to reveal which molecules are being up-regulated and down-regulated. Furthermore, it generates a visual understanding of how rejection occurs in different pathways."),
                                    
                            ),
                            
                          ),
                          
                          tabPanel(
                            title = "Annotation Analysis",
                            
                            tabsetPanel(
                                        tabPanel(
                                          title = "Network Plot",
                                          includeCSS("www/d.css"),
                                          h2("Significant Genes selected by Over-representation Analysis (ORA)"),
                                          tags$div(class="d1",
                                          tags$div(class="m",
                                          h4("Overlapping genes between human biological states and processes."),
                                          p("The network plot shows the intersection of over-represented genes (ORA) between both GSE36059 & GSE48581 after regressions on antibody-mediated rejection versus healthy. Central categories represent the top six (by p.value) significant biological states where size shows the number of representative ORA genes which are present. Notice the large fold changes of genes which participate in Interferon gamma, Interferon alpha, Inflammatory responses and IL6 JAK STAT3 signalling."),
                                          ),
                                          tags$div(class = "a",
                                          plotOutput("networkPlotABMR")
                                          ),
                                        ),
                                        ),
                                        
                                        tabPanel(
                                          title = "Dot Plot",
                                          includeCSS("www/d.css"),
                                          h2("Gene Ratio of Significant Biological States/Processes"),
                                          tags$div(class="d1",
                                                   tags$div(class="m",
                                                            h4("Proportional representation of biological states and processes."),
                                                            p("This dot plot shows all the biological states and process that were found to be significant given the intersection of over-represented genes (ORA) between both GSE36059 & GSE48581. Gene ratio denotes the ratio of relevant ORA genes to all possible relevant genes for a given biological process. Notice that Gamma response and Allograft reject are very well represented."),
                                                   ),
                                                   tags$div(class = "a",
                                                            plotOutput("dotPlotABMR")
                                                   ),
                                          ),
                                          
                                          
                                        ),
                                        
                                        tabPanel(
                                          title = "Tree Plot",
                                          includeCSS("www/d.css"),
                                          h2("Gene Ratio of Significant Biological States/Processes"),
                                          tags$div(class="d1",
                                                   tags$div(class="m",
                                                            h4("Similarity of genes clustered by Jaccard similarity index."),
                                                            p("This tree plot shows the clustering of ORA genes that represent for significant biological states and processes. The Jaccard index compares sets by taking the ratio of the intersection with the union. The highlight shows a family of process that are closely related."),
                                                   ),
                                                   tags$div(class = "a",
                                                            plotOutput("treePlotABMR")
                                                   ),
                                          ),
                                          
                                        ),
                                        
                                        tabPanel(
                                          title = "KEGG Enrichment Analysis",
                                          h2("KEGG Pathway Explanation"),
                                          p("KEGG is a database of biological molecular pathways that allow researchers to ascribe meaningful physiological annotations to genomic information. These images show the significant pathways based on the intersection of over-represented genes (ORA) between both GSE36059 & GSE48581. Where there is a significantly high fold ratio of a gene responsible/related to the transcribing of a biological signalling molecule or messenger, the molecule is marked red. e.g. Notice that there is a consistent representation of MHC-I, MHC-II present which primarily commonly involved in autoimmune dysfunction."),
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
                                                  
                                                  titlePanel("Upload Files"),
                                                  p("Note: Only Affymetrix data will produce valid results."),
                                                  fileInput('target_upload1', 'Choose file to upload',
                                                            accept = c(
                                                                    'text/csv',
                                                                    'text/comma-separated-values',
                                                                    '.csv'
                                                            )),
                                                  tags$hr(),
                                                  uiOutput('mybutton')
                                          ),
                                          mainPanel(
                                                  tabsetPanel(type = "tabs",
                                                              tabPanel("KNN", plotlyOutput("knn"),plotlyOutput("knnr")),
                                                              tabPanel("Logistic Regression",plotlyOutput("log"),plotlyOutput("logr")),
                                                              tabPanel("Random Forest",plotlyOutput("rf"),plotlyOutput("rfr")),
                                                              tabPanel("Support Vector Machine",plotlyOutput("svm"),plotlyOutput("svmr")),
                                                              tabPanel("Decision Tree",plotlyOutput("tree"),plotlyOutput("treer"))
                                                              
                                                  )
  
                                          )
                                  )
                          ),
                          tabPanel(
                                  title="Manual Input",
                                  sidebarLayout(
                                          
                                          sidebarPanel(
                                                  uiOutput("mysliders1"),
                                          ),
                                          
                                          mainPanel(
                                                  tabsetPanel(type = "tabs",
                                                              tabPanel("KNN", plotlyOutput("knn3"),plotlyOutput("knnr2")),
                                                              tabPanel("Logistic Regression",plotlyOutput("log3"),plotlyOutput("logr2")),
                                                              tabPanel("Random Forest",plotlyOutput("rf3"),plotlyOutput("rfr2")),
                                                              tabPanel("Support Vector Machine",plotlyOutput("svm3"),plotlyOutput("svmr2")),
                                                              tabPanel("Decision Tree",plotlyOutput("tree3"),plotlyOutput("treer2"))
                                                              
                                                  )
                                                  
                                          )
                                  )
                          )
                          
              )
     )
     )
     ),

     tabPanel("T-Cell Mediated Rejection Analysis", 
              
              tabsetPanel(id = "tabset2",
                          
                          tabPanel(
                            title = "Overview",
                            h3("What is TCMR?"),
                            p("Acute T-Cell Mediated Rejection (TCMR) occurs when the instersitium is infiltrated by T cells and macrophages, and is characterised by decreased urine output and proteinuria. (Halloran P. F., 2010). TCMR affects about 10-12% of patients and is a less common cause of rejection in comparison to other conditions. This is because TCMR can be more easily controlled with immunosupression. TCMR develops as soon as a week after transplantation, and its risk of rejection reduces over time (Halloran P. F., 2010). As such the risk of rejection due to TCMR peaks in the short term and slowly reduces over time."),
                            h3("Why Explore TCMR?"),
                            p("Despite the minimal risk of acute kidney rejection arising from pure TCMR, it is still a vital component to analyse in research. Given that T-cells are often involved in the activation of B-cells, it is assumed that ABMR cannot occur without the pre-existence of TCMR. To account for this, this app will also target the identification of genes that are attributed to TCMR. Furthermore, this can also assist clinical researchers in determining new methods to optimise the level of immunosuppressive therapy used for each patient."),
                            h3("Annotation Analysis"),
                            tags$ul(
                              tags$li("Network Plot - This plot depicts the biological states and processes of gene profiles through the aggregation of gene data sets."),
                              tags$li("Dot plot - Depicts the proportion of significant genes present out of the total count of genes that are associated with the annotated molecular process."),
                              tags$li("Tree Plot - The grouping of genes portrays different families of the molecular process. The branches represent the similarity of genes that are markers of multiple processes."),
                              tags$li("KEGG Enrichment Analysis - Shows the physiological process of different gene pathways to reveal which molecules are being up-regulated and down-regulated. Furthermore, it generates a visual understanding of how rejection occurs in different pathways."),
                                    
                            ),
                          ),
                          
                          tabPanel(
                            title = "Annotation Analysis",
                            
                            tabsetPanel(
                              tabPanel(
                                title = "Network Plot",
                                includeCSS("www/d.css"),
                                h2("Significant Genes selected by Over-representation Analysis (ORA)"),
                                tags$div(class="d1",
                                         tags$div(class="m",
                                                  h4("Overlapping genes between human biological states and processes."),
                                                  p("The network plot shows the intersection of over-represented genes (ORA) between both GSE36059 & GSE48581 after regressions on T-cell mediated rejection versus healthy. Central categories represent the top six (by p.value) significant biological states where size shows the number of representative ORA genes which are present. Notice the large fold changes of genes which participate in Interferon gamma, Interferon alpha and Inflammatory responses."),
                                         ),
                                         tags$div(class = "a",
                                                  plotOutput("networkPlotTCMR")
                                         ),
                                ),
                                
            
                              ),
                              
                              tabPanel(
                                title = "Dot Plot",
                                includeCSS("www/d.css"),
                                h2("Gene Ratio of Significant Biological States/Processes"),
                                tags$div(class="d1",
                                         tags$div(class="m",
                                                  h4("Proportional representation of biological states and processes."),
                                                  p("This dot plot shows all the biological states and process that were found to be significant given the intersection of over-represented genes (ORA) between both GSE36059 & GSE48581. Gene ratio denotes the ratio of relevant ORA genes to all possible relevant genes for a given biological process. Notice that Gamma response and Allograft reject are very well represented."),
                                         ),
                                         tags$div(class = "a",
                                                  plotOutput("dotPlotTCMR")
                                         ),
                                ),
                              ),
                              
                              tabPanel(
                                title = "Tree Plot",
                                includeCSS("www/d.css"),
                                h2("Gene Ratio of Significant Biological States/Processes"),
                                tags$div(class="d1",
                                         tags$div(class="m",
                                                  h4("Similarity of genes clustered by Jaccard similarity index."),
                                                  p("This tree plot shows the clustering of ORA genes that represent for significant biological states and processes. The Jaccard index compares sets by taking the ratio of the intersection with the union. The highlight shows a family of process that are closely related."),
                                         ),
                                         tags$div(class = "a",
                                                  plotOutput("treePlotTCMR")
                                         ),
                                ),
                                
                              ),
                              
                              tabPanel(
                                title = "KEGG Enrichment Analysis",
                                h2("KEGG Pathway Explanation"),
                                p("KEGG is a database of biological molecular pathways that allow researchers to ascribe meaningful physiological annotations to genomic information. These images show the significant pathways based on the intersection of over-represented genes (ORA) between both GSE36059 & GSE48581. Where there is a significantly high fold ratio of a gene responsible/related to the transcribing of a biological signalling molecule or messenger, the molecule is marked red. e.g. Notice that there is a consistent representation of MHC-I, MHC-II present which primarily commonly involved in autoimmune dysfunction."),
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
                            p("Note: Only Affymetrix data will produce valid results."),
                            fileInput('target_upload', 'Choose file to upload',
                                      accept = c(
                                              'text/csv',
                                              'text/comma-separated-values',
                                              '.csv'
                                      )),
                            tags$hr(),
                            uiOutput('mybutton1')
                                    ),
                            mainPanel(
                                    tabsetPanel(type = "tabs",
                                                tabPanel("KNN", plotlyOutput("knn1"),plotlyOutput("knnr1")),
                                                tabPanel("Logistic Regression",plotlyOutput("log1"),plotlyOutput("logr1")),
                                                tabPanel("Random Forest",plotlyOutput("rf1"),plotlyOutput("rfr1")),
                                                tabPanel("Support Vector Machine",plotlyOutput("svm1"),plotlyOutput("svmr1")),
                                                tabPanel("Decision Tree",plotlyOutput("tree1"),plotlyOutput("treer1"))
                                                
                                    )
                                    
                            )
                            )
                          ),
                          tabPanel(
                          title="Manual Input",
                          sidebarLayout(
                                  
                                  sidebarPanel(
                                    uiOutput("mysliders"),
                                  ),
                                  
                                  mainPanel(
                                          tabsetPanel(type = "tabs",
                                                      tabPanel("KNN", plotlyOutput("knn2"),plotlyOutput("knnr3")),
                                                      tabPanel("Logistic Regression",plotlyOutput("log2"),plotlyOutput("logr3")),
                                                      tabPanel("Random Forest",plotlyOutput("rf2"),plotlyOutput("rfr3")),
                                                      tabPanel("Support Vector Machine",plotlyOutput("svm2"),plotlyOutput("svmr3")),
                                                      tabPanel("Decision Tree",plotlyOutput("tree2"),plotlyOutput("treer3"))
                                                      
                                          )
                                          
                                  )
                          )
              )
     )
     
)))))
