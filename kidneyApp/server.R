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
  
  
  mydata <- reactive({
    File <- input$target_upload1
    if (is.null(File))
      return(NULL)
    ab <- read.csv(File$datapath)
    a<-get_pairwise_differences_probe_id(abmr_nonrej_features,ab)
    return(a)
  })
  
  
  output$mybutton <- renderUI({
    if (is.null(mydata())) return(NULL)
    actionButton("action", "Show Prediction")
  })
  
# class_model = {"log", "svm", "tree", "rf", "knn"}
  observeEvent(input$action, {
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
    
    output$knnr<-renderPlotly({
      get_cross_val_plot(nrow(abmr_nonrej_features), 5, 10, abmr_nonrej_features, abmr_nonrej_outcome,"k-Nearest-Neighbours")
    })
    
    output$logr<-renderPlotly({
      get_cross_val_plot(nrow(abmr_nonrej_features), 5, 10, abmr_nonrej_features, abmr_nonrej_outcome,"Logistic Regression")
    })
    
    output$svmr<-renderPlotly({
      get_cross_val_plot(nrow(abmr_nonrej_features), 5, 10, abmr_nonrej_features, abmr_nonrej_outcome,"Simple Vector Machine")
    })
    
    output$treer<-renderPlotly({
      get_cross_val_plot(nrow(abmr_nonrej_features), 5, 10, abmr_nonrej_features, abmr_nonrej_outcome,"Tree")
    })
    
    output$rfr<-renderPlotly({
      get_cross_val_plot(nrow(abmr_nonrej_features), 5, 10, abmr_nonrej_features, abmr_nonrej_outcome,"Random Forest")
    })
    
  })
  
  mydata1<- reactive({
    File <- input$target_upload
    if (is.null(File))
      return(NULL)
    abc <- read.csv(File$datapath)
    d<-get_pairwise_differences_probe_id(tcmr_nonrej_features,abc)
    return(d)
  })
  
  output$mybutton1 <- renderUI({
    if (is.null(mydata1())) return(NULL)
    actionButton("action1", "Show Prediction")
  })
  
  
  observeEvent(input$action1, {
  # class_model = {"log", "svm", "tree", "rf", "knn"}
  output$knn1 <-renderPlotly({
    get_PCA_plot(tcmr_nonrej_features,tcmr_nonrej_outcome,mydata1(),"TCMR","knn")
  })
  output$log1<-renderPlotly({
    get_PCA_plot(tcmr_nonrej_features,tcmr_nonrej_outcome,mydata1(),"TCMR","log")
  }) 
  output$svm1<-renderPlotly({
    get_PCA_plot(tcmr_nonrej_features,tcmr_nonrej_outcome,mydata1(),"TCMR","svm")
  })
  output$tree1<-renderPlotly({
    get_PCA_plot(tcmr_nonrej_features,tcmr_nonrej_outcome,mydata1(),"TCMR","tree")
  })
  output$rf1<-renderPlotly({
    get_PCA_plot(tcmr_nonrej_features,tcmr_nonrej_outcome,mydata1(),"TCMR","rf")
  })
  
  
  output$knnr1<-renderPlotly({
    get_cross_val_plot(nrow(tcmr_nonrej_features), 5, 10, tcmr_nonrej_features, tcmr_nonrej_outcome,"k-Nearest-Neighbours")
  })
  
  output$logr1<-renderPlotly({
    get_cross_val_plot(nrow(tcmr_nonrej_features), 5, 10, tcmr_nonrej_features, tcmr_nonrej_outcome,"Logistic Regression")
  })
  
  output$svmr1<-renderPlotly({
    get_cross_val_plot(nrow(tcmr_nonrej_features), 5, 10, tcmr_nonrej_features, tcmr_nonrej_outcome,"Simple Vector Machine")
  })
  
  output$treer1<-renderPlotly({
    get_cross_val_plot(nrow(tcmr_nonrej_features), 5, 10, tcmr_nonrej_features, tcmr_nonrej_outcome,"Tree")
  })
  
  output$rfr1<-renderPlotly({
    get_cross_val_plot(nrow(tcmr_nonrej_features), 5, 10, tcmr_nonrej_features, tcmr_nonrej_outcome,"Random Forest")
  })
  
  })
  
  # mydata1<- reactive({
  #   File <- input$target_upload
  #   if (is.null(File))
  #     return(NULL)
  #   abc <- read.csv(File$datapath)
  #   d<-get_pairwise_differences_probe_id(tcmr_nonrej_features,abc)
  #   return(d)
  # })
  
  output$mysliders <- renderUI({
    a=get_genes_for_sliders(tcmr_nonrej_features)
    sliders <- lapply(1:length(a), function(i) {
      inputName <- a[i]
      sliderInput(inputName, inputName, min=5, max=20, value=10)
    })
    do.call(tagList, sliders)
  })
  
  
  mydata3<- reactive({
    d=c()
    e=c()
    a=get_genes_for_sliders(tcmr_nonrej_features)
    for(i in 1:length(a)){
      d=c(d,a[i])
      e=c(e,input[[a[i]]])
    }
    data3 <-data.frame(Column1 = d,Column2=e)
    m<-get_pairwise_differences_gene_symbol(tcmr_nonrej_features,data3)
    return(m)
  })
  

  output$knn2 <-renderPlotly({
     get_PCA_plot(tcmr_nonrej_features,tcmr_nonrej_outcome,mydata3(),"TCMR","knn")
   })
  output$log2 <-renderPlotly({
    get_PCA_plot(tcmr_nonrej_features,tcmr_nonrej_outcome,mydata3(),"TCMR","log")
  })
  output$svm2 <-renderPlotly({
    get_PCA_plot(tcmr_nonrej_features,tcmr_nonrej_outcome,mydata3(),"TCMR","svm")
  })
  output$tree2 <-renderPlotly({
    get_PCA_plot(tcmr_nonrej_features,tcmr_nonrej_outcome,mydata3(),"TCMR","tree")
  })
  output$rf2 <-renderPlotly({
    get_PCA_plot(tcmr_nonrej_features,tcmr_nonrej_outcome,mydata3(),"TCMR","rf")
  })
  
  output$mysliders1 <- renderUI({
    a=get_genes_for_sliders(abmr_nonrej_features)
    sliders <- lapply(1:length(a), function(i) {
      inputName <- a[i]
      sliderInput(inputName, inputName, min=5, max=20, value=10)
    })
    do.call(tagList, sliders)
  })
  
  
  mydata4<- reactive({
    d=c()
    e=c()
    a=get_genes_for_sliders(abmr_nonrej_features)
    for(i in 1:length(a)){
      d=c(d,a[i])
      e=c(e,input[[a[i]]])
    }
    data3 <-data.frame(Column1 = d,Column2=e)
    m<-get_pairwise_differences_gene_symbol(abmr_nonrej_features,data3)
    return(m)
  })
  
  
  output$knn3 <-renderPlotly({
    get_PCA_plot(abmr_nonrej_features,abmr_nonrej_outcome,mydata4(),"TCMR","knn")
  })
  output$log3 <-renderPlotly({
    get_PCA_plot(abmr_nonrej_features,abmr_nonrej_outcome,mydata4(),"TCMR","log")
  })
  output$svm3 <-renderPlotly({
    get_PCA_plot(abmr_nonrej_features,abmr_nonrej_outcome,mydata4(),"TCMR","svm")
  })
  output$tree3 <-renderPlotly({
    get_PCA_plot(abmr_nonrej_features,abmr_nonrej_outcome,mydata4(),"TCMR","tree")
  })
  output$rf3 <-renderPlotly({
    get_PCA_plot(abmr_nonrej_features,abmr_nonrej_outcome,mydata4(),"TCMR","rf")
  })
})
