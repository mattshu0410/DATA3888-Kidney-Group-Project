library(shiny)
library(shinythemes)
library(shinyjs)
library(DT)
library(shinyalert)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

  output$dotPlotABMR <- renderPlot({
    make_dotplot(FC_ABMR_ENTREZID)
  })  
  
  output$networkPlotABMR <- renderPlot({
    make_network_plot(FC_ABMR_ENTREZID)
  })  
  
  output$treePlotABMR <- renderPlot({
    make_treeplot(FC_ABMR_ENTREZID)
  })  
  
  output$dotPlotTCMR <- renderPlot({
    make_dotplot(FC_TCMR_ENTREZID)
  })  
  
  output$networkPlotTCMR <- renderPlot({
    make_network_plot(FC_TCMR_ENTREZID)
  })  
  
  output$treePlotTCMR <- renderPlot({
    make_treeplot(FC_TCMR_ENTREZID)
  }) 
  
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
    current_data = mydata()
    if (!is.null(current_data)) {
      shinyjs::show("showTab")
    }
  })
  
  mydata <- reactive({
    File <- input$target_upload1
    ab <- read.csv(File$datapath)
    a<-get_pairwise_differences_probe_id(abmr_nonrej_features,ab)
    return(a)
  })
  
# class_model = {"log", "svm", "tree", "rf", "knn"}
  output$knn <-renderPlotly({
    get_PCA_plot(abmr_nonrej_features,abmr_nonrej_outcome,mydata(),"ABMR","knn")
  })
  output$log<-renderPlotly({
    get_PCA_plot(abmr_nonrej_features,abmr_nonrej_outcome,mydata(),"ABMR","log")
  }) 
  output$svm<-renderPlotly({
    get_PCA_plot(abmr_nonrej_features,abmr_nonrej_outcome,mydata(),"ABMR","svm")
  })
  output$tree<-renderPlotly({
    get_PCA_plot(abmr_nonrej_features,abmr_nonrej_outcome,mydata(),"ABMR","tree")
  })
  output$rf<-renderPlotly({
    get_PCA_plot(abmr_nonrej_features,abmr_nonrej_outcome,mydata(),"ABMR","rf")
  })
  
  
  output$sliders <- renderUI({
    a=get_genes_for_sliders(tcmr_nonrej_features)
    sliders <- lapply(1:length(a), function(i) {
      inputName <- a[i]
      sliderInput(inputName, inputName, min=5, max=20, value=10)
    })
    do.call(tagList, sliders)
  })
  
  
})
