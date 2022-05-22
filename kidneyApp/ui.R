library(shiny)
library(shinythemes)
library(shinyjs)
library(DT)
library(shinyalert)
library(prompter)

# Define UI for application that draws a histogram
shinyUI(navbarPage(tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "mystyle.css")
),
                   theme = shinytheme("flatly"),
                   title = img(src = "logo1.jpg",inline=TRUE, id="logo"),
                   fluid=TRUE, 

     useShinyjs(),
     use_prompt(),

     tabPanel("Home",
              h2("Introduction"),
              p("Welcome to Kidney C1’s Shiny App! This app is a kidney rejection analysis tool for pharmaceutical researchers that aims to increase our understanding of the influence of different genes on kidney graft acceptance or rejection. The tool demonstrates the influence that different genes have on T-Cell Mediated and Antibody-Mediated Rejection by predicting the rejection of different gene profiles and comparing gene pathways on data analysis plots. By revealing the genes that attribute to ABMR/TCMR and acceptance, this app will provide pharmaceutical researchers with the information necessary to create medication that targets the  gene pathways responsible for rejection, imcluding a less destructive immunosuppressant, or new medication that can reduce the effects of rejection ."),
              h2("Why?"),
              p("Kidney rejection occurs within 10-15 patients per 100 transplants. Since more than 90,000 patients are in need of a kidney and only 20,000 kidney transplants occur each year, the utility of a kidney rejection analysis tool for pharmaceutical researchers is considerable ('Rejection of a transplanted kidney', 2022).  The creation of a tool that can be used in the drug discovery pipeline to aid the identification of new transcription-specific targets will aid the mitigation of kidney rejection risk by providing information that will improve the treatment of ABMR and TCMR. With the goal of maximising the greater wellbeing of individuals in need of kidney transplants, this tool will provide advancements in knowledge for pharmaceuticals researchers."),
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
              h3("Disclaimer"),
              p("We have utilised numerous datasets in the creation of our model that were collected from different sources, and compiled in tandem to create our tool. The data was utilised for educational purposes, and the app should only be used by pharmaceutical researchers. It accepts no liability for the quality or accuracy of its predictions/information."),
              h3("Datasets Used"),
              h5("GSE36059"),
              p("Halloran et al. collected 703 unselected biopsies 3 days to 35 years after the transplant. Using microarrays and conventional methods they sought to discern the adaptive changes in the alloimmune response over time.  Through conventional methods it was found that 228 biopsy specimens were rejected, where 67 had pure TCMR, 110 AMBR and 28 were mixed between the two. With microarrays it was observed that 228 kidneys were rejected, with 76 pure TCMR, 124 pure ABMR and 28 being mixed. Consequently it was noted that ABMR was strongly associated with increased kidney loss, whilst TCMR was not, and TCMR appeared early but disappeared overtime, whilst ABMR appeared usually a year after transplant."),
              h5("GSE48581"),
              p("Halloran et al. collected international biopsies samples from from 264 patients. These were analysed using microarrays to test for TCMR and then assigning scores based on an algorithm that was established using a previously collected reference set of 403 biopsies. These scores were also compared to the local histological evaluation. The accuracy between both the reference dataset and the international clinical trial were similar (89) and (87) respectively, and discrepancies arose due to known histology limitations due to biopsies with scarring or inflammation caused by other diseases. However the TMCR score and the histologic TCMR were not associated with kidney loss, yet the TMCR score can provide greater insight in situations where histology is enigmatic or possibly deceptive."),
        
     ),
     
     tabPanel("Antibody-Mediated Rejection Analysis", 
              
              tabsetPanel(id = "tabset",
                          
                          tabPanel(
                            title = "Overview",
                            h3("What is ABMR?"),
                            p("Antibody-Mediated Rejection (ABMR) is a more common form of rejection and presents a significant challenge for long-term graft survival (Kojc & Haler, 2022). ABMR occurs when anti-donor specific antibodies such as anti-HLA antibodies work against the transplanted kidney, and incite its rejection. Unlike TCMR, the risk of kidney rejection from ABMR is a long term concern that increases overtime. Due to advancements in medical technology, ABMR now poses the largest threat to kidney rejection, accounting for graft failure in 64% of rejected transplants. ABMR peaks 5 years after transplantation, and therefore is more difficult to manage.(Kwon, H. et al, 2021)"),
                            h3("Annotation Analysis"),
                            tags$ul(
                                    tags$li("Network Plot - This plot depicts the biological states and processes of our gene profiles through the aggregation of gene data sets."),
                                    tags$li("Dot plot - Depicts the proportion of significant genes present out of the total count of genes associated with the annotated molecular process."),
                                    tags$li("Tree Plot - The grouping of genes  portrays different families of the molecular process. The branches represent the similarity of genes across multiple molecular processes."),
                                    tags$li("KEGG Enrichment Analysis - Shows the physiological process of different gene pathways to reveal which functional-level proteins/protein complexes are being upregulated and down regulated. Furthermore, it generates a visual understanding of how rejection occurs in different pathways to aid pharmaceutical researchers understanding."),
                                    
                            ),
                            
                          ),
                          
                          tabPanel(
                            title = "Annotation Analysis",
                            
                            tabsetPanel(
                                        tabPanel(
                                          title = "Network Plot",
                                          includeCSS("www/d.css"),
                                          h2("Significant Genes Selected by Over-Representation Analysis (ORA)"),
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
                                  title = "Classification",
                                  
                                  tabsetPanel(id = "tabset4",
                                              
                                  tabPanel(
                                  title= "File Input",
                                  sidebarLayout(
                                          
                                          sidebarPanel(
                                                  
                                                  titlePanel("Upload Files"),
                                                  p("Note: Only Affymetrix data will produce valid results."),
                                                  add_prompt(
                                                  fileInput('target_upload1', 'Choose file to upload',
                                                            accept = c(
                                                                    'text/csv',
                                                                    'text/comma-separated-values',
                                                                    '.csv'
                                                            )),
                                                  position = "bottom", 
                                                  size = "medium",
                                                  message = "Please only upload a csv file or a text file (comma separated values) with Affymetrix data that uses the control probe AFFX-TrpnX-M_at.",
                                                  ),
                                                  tags$hr(),
                                                  uiOutput('mybutton'),
                                                  
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
                            p("The next rejection type explored through our app is Acute T-Cell mediated rejection (TCMR). TCMR occurs when the instersitium is infiltrated by T cells and macrophages and is characterised by decreased urine output and proteinuria. (Halloran P. F., 2010). TCMR affects about 10-12% of patients and  is a less common cause of rejection in comparison to other conditions. This is because TCMR can be easily controlled as numerous medicinal methods of management have been introduced. TCMR develops as soon as a week after transplantation and its risk of rejection reduces over time (Halloran P. F., 2010). As such the risk of rejection due to TCMR peaks in the short term and slowly reduces over time."),
                            h3("Why Explore TCMR?"),
                            p("Despite the minimal risk of acute kidney rejection arising from pure TCMR, it is still a vital component to analyse within this report. It is basic biological knowledge that when the immune system is fighting forreign cells, its first response is to mobilise the utility of T-cells. Subsequently, since T-cells are necessary for the activation of B-cells, ABMR cannot occur without the preexistence of TCMR. Accordingly, to account for this ideology, this app will also target the identification of genes that attribute to TCMR, because it is a contributing factor to the inducement of AMBR. Furthermore, this can also assist pharmaceutical  researchers in determining new methods to optimise the level and method of immunosuppressive therapy used  for each patient."),
                            h3("Annotation Analysis"),
                            tags$ul(
                              tags$li("Network Plot - This plot depicts the biological states and processes of our gene profiles through the aggregation of gene data sets."),
                              tags$li("Dot plot - Depicts the proportion of significant genes present out of the total count of genes associated with the annotated molecular process."),
                              tags$li("Tree Plot - The grouping of genes  portrays different families of the molecular process. The branches represent the similarity of genes across multiple molecular processes."),
                              tags$li("KEGG Enrichment Analysis - Shows the physiological process of different gene pathways to reveal which functional-level proteins/protein complexes are being upregulated and down regulated. Furthermore, it generates a visual understanding of how rejection occurs in different pathways to aid pharmaceutical researchers understanding."),
                              
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
                            title = "Classification",
                            tabsetPanel(id = "tabset3",
                            
                            tabPanel(
                            title= "File Input",
                            
                            sidebarLayout(
                                    
                                    sidebarPanel(
                                            
                            titlePanel("Uploading Files"),
                            p("Note: Only Affymetrix data will produce valid results."),
                            add_prompt(
                            fileInput('target_upload', 'Choose file to upload',
                                      accept = c(
                                              'text/csv',
                                              'text/comma-separated-values',
                                              '.csv'
                                      )),
                            position = "bottom",
                            size = "medium",
                            message = "Please only upload a csv file or a text file (comma separated values) with Affymetrix data that uses the control probe AFFX-TrpnX-M_at."
                            ),
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
