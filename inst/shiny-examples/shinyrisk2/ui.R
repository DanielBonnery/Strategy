  library(shiny)
  library(shinythemes)
  library(leaflet)
  library(ggplot2)
  library(Strategy)
  library(RLeafletTools)
  
  shinyUI(
    fluidPage(
      # Application title
      #setBackgroundColor("black"),
      navbarPage(
        "Type 2 Risk at the field level",id="main_navbar",theme = shinytheme("slate"),
        tabPanel(
          ###################### Tabpanel 1
          "Formula",
          sidebarLayout(
            # Sidebar with a slider and selection inputs
            sidebarPanel(
              fluidRow(
                column(
                  6,style=list("padding-right: 5px;"),
                  h4("Population parameters"),
                  sliderInput("N","N, number of trees in the field",min = 1,  max = 100000,  value = 32369)),
                column(
                  6,style=list("padding-right: 5px;"),
                  h4("Epidemic parameters"),
                  sliderInput("m","infected trees",min = 1,  max = 300,  value = 100),
                  sliderInput("foyers","infection sources",min = 1,  max = 10,  value = 1))),
              h4("Single leaf level risks"),
              fluidRow(
                column(
                  4,style=list("padding-right: 5px;"),
                  sliderInput("betac","False negative",min = 0,  max = 1,  value = .1,round=F)),
                column(
                  4,style=list("padding-right: 5px;"),
                  sliderInput("alphac","False positive",min = 0,  max = 1,  value = 0,round=F))),
              h4("Sampling parameters"),
              h5("Number of selected ..."),
              fluidRow(
                column(
                  4,style=list("padding-right: 5px;"),
                  sliderInput("n0",h5("...bulks,"),min = 1,  max = 300,  value = 100)),
                column(
                  4,style=list("padding-left: 5px;"),
                  sliderInput("n1",h5("...trees / bulk,"),min = 1,  max = 100,  value = 8)),
                column(
                  4,style=list("padding-left: 5px;"),
                  sliderInput("n2",h5("...leaves / tree."),min = 1,  max = 300,  value = 4))),
              fluidRow(
                column(
                  6,style=list("padding-right: 5px;"),
                  radioButtons("samplemethod", 
                               h5("Selection method:"),
                               choices = list("Simple random sampling" = "srs", "Stratified/systematic sampling" = "systematic"),
                               selected = "srs")),
                column(
                  6,style=list("padding-right: 5px;"),
                  actionButton("refresh", "Refresh map"),
                  actionButton("rerun", "Re-draw 1000 samples",color="red")))),
            mainPanel(
              tabsetPanel(
                type="tabs",
                tabPanel(
                  "Formula",
                  plotOutput("plot.1",height = 500),
                  h3(textOutput("risk")),
                  tableOutput("counts")),
                tabPanel("Visualisation",
                         leafletOutput("bbmap1.4", height=600))
              ))))#,
        ###################### Tabpanel 2
        #tabPanel(
        #  "Real field simulations",
        #  id="map1",
        #  sidebarLayout(
        #    # Sidebar with a slider and selection inputs
        #    sidebarPanel(
        #      fluidRow(h3("Simulation parameters:")),
        #      fluidRow("Simple Random Sampling"),
        #      fluidRow("m=100"),
        #      fluidRow("N=82000"),
        #      fluidRow("n x n'' / n'=100")),
        #      mainPanel(leafletOutput("bbmap2", height=600))))
        ###################### Tabpanel 3
        #tabPanel(
        #  "Square field simulations",id="map2",
        #  mainPanel(
        #    leafletOutput("bbmap2", height=300)))
      )))