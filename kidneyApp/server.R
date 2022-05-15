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
  
})