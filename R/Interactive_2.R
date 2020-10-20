Interactive2<-function(){
  
  #' runCompare
  #' 
  #' @description Shiny App to 
  #' @export
  #' @examples
  #' package1<-NULL
  #' package2<-NULL
  #' runCompare()
    library(dplyr)
    library(gridExtra)
    library(ggplot2)
    library(shinythemes)
    library(Strategy)
    library(shiny)
    
    library(dplyr)
    
    library(leaflet)
    
    library(DT)
    library(Strategy)
    data(Avo_fields,package="Strategy")
    data(U2E2,package="Strategy")
    
    U2E2<-U2E2[c(2:4,12:112)]
    names(U2E2)[1:2]<-c("long","lat")
    #Plot the trees
    Avo_fields$Source_yr<-addNA(as.factor(Avo_fields$Source_yr))
    
    QLD<-Avo_fields$State=="Qld"
    
    Avo_ids<-unique(Avo_fields[QLD,]$Avo_id)[1:5]
    QLD<-QLD&is.element(Avo_fields$Avo_id,Avo_ids)
    QLDt<-is.element(U2$id,Avo_ids)
    
    
    yearpal <- colorFactor(heat.colors(5),domain = levels(Avo_fields$Source_yr),na.color = "#aaff56")
    
    AllTrees<-merge(U2E2,Avo_fields[c("Avo_id","State")],by.x="id",by.y="Avo_id")
    if(FALSE){
      leaflet(Avo_fields[QLD,]) %>%
        addProviderTiles('Esri.WorldImagery',options = providerTileOptions(minZoom = 1, maxZoom = 21,maxNativeZoom=19)) %>% 
        addProviderTiles("CartoDB.PositronOnlyLabels")%>% 
        addPolylines(fillOpacity = 1, weight = 3, smoothFactor = 0.5,opacity = 1.0,
                     color=~yearpal(Avo_fields[QLD,]$Source_yr),
                     fillColor=~yearpal(Avo_fields[QLD,]$Source_yr))%>%
        addpiechartclustermarkers(.data=Trees(),
                                  ,.colors=c("green","red","orange","purple","black"),
                                  group="I081")}
    
    ui <- navbarPage(theme = shinytheme("slate"),
                     title="Detect epidemic",
                     id="main",
                     #                   tabPanel(title= "Interactive Map",
                     #                            leafletOutput("bbmap", height=1000)),
                     tabPanel("Interactive Map",
                              sidebarLayout(
                                sidebarPanel(width = 3, 
                                             sliderInput("time",
                                                         "Time (in days)",
                                                         min = 0,
                                                         max = 100,
                                                         value = 0),
                                             checkboxInput("piechart", "Pie cluster markers", value = TRUE, width = NULL),
                                             selectInput(inputId = "Status", label = "Disease status", selected = "infected",multiple = T,
                                                         choices = c("susceptible","exposed","cryptic","infected","removed")),
                                             selectInput(inputId = "State", label = "State", selected="Qld",
                                                         choices = unique(Avo_fields$State)),
                                ),
                                
                                mainPanel(plotOutput('toto'),
                                          textOutput('totox'),
                                          leafletOutput("bbmap", height=1000)
                                          
                                )
                              )
                     ),
                     tabPanel(title="Data", DT::dataTableOutput("data")))
    
    ################################# Server
    server <- function(input, output) {
      
      #library(datagovuk)
      library("sp")
      library("rgdal")
      library(dplyr)
      library(leaflet)
      
      
      infectionattime <- reactive({
        y<-paste0("I",formatC(input$time, width =3, format = "d", flag = "0"))
      })
      
      
      StateAvoFieldsdata <- reactive({
        State_Avo_fields<-Avo_fields[Avo_fields$State==input$State,]
      })
      
      yearpal <- colorFactor(heat.colors(5),domain = levels(Avo_fields$Source_yr),na.color = "#aaff56")
      Trees<-reactive({AllTrees[AllTrees$State==input$State,]})
      
      subsetData <- reactive({
        A<-Trees()[c("long","lat",paste0("I",formatC(input$time, width =3, format = "d", flag = "0")))]
        names(A)[3]<-"Imap"
        A[is.element(A$Imap,input$Status),]})
      
      #y<-paste0("I",formatC(input$time, width =3, format = "d", flag = "0"))
      #D<-if(input$Status=="Pie-chart"){Trees}else{Trees[is.element(Trees[[y]],input$Status),]}
      #D$II<-D[[y]]
      #D
      
      output$toto <- renderPlot({
        ggplot(subsetData(), aes(x=Imap)) + geom_bar()
      })
      output$totox<-reactive({capture.output(summary(subsetData()["Imap"]))})
      
      data(UE,package="Strategy")
      #bb_data<-bb_data[bb_data$ReportedDay<input$time,]
      # create a color paletter for category type in the data file
      
      #  popbins<-quantile(shapeData$population,(seq_len(11)-1)/10) 
      # poppal <- colorBin(heat.colors(5), bins=popbins, na.color = "#aaff56",reverse = T)
      #  popbins<-quantile(shapeData$population,(seq_len(11)-1)/10) 
      #  poppal <- colorBin(heat.colors(5), bins=popbins, na.color = "#aaff56",reverse = T)
      library(leaflet)
      library(Strategy)
      
      
      output$bbmap <- renderLeaflet({
        
        leaflet(Avo_fields[Avo_fields$State=="Qld",]) %>%
          addProviderTiles('Esri.WorldImagery',options = providerTileOptions(minZoom = 1, maxZoom = 21,maxNativeZoom=19)) %>% 
          addProviderTiles("CartoDB.PositronOnlyLabels")%>% 
          addPolylines(fillOpacity = 1, weight = 3, smoothFactor = 0.5,opacity = 1.0,
                       color=~yearpal(StateAvoFieldsdata()$Source_yr),
                       fillColor=~yearpal(StateAvoFieldsdata()$Source_yr))
      })
      
      observe({
        proxy<-leafletProxy("bbmap",data=StateAvoFieldsdata())
      })
      
      
      
      observe({
        proxy<-leafletProxy("bbmap")%>%clearMarkerClusters()%>%clearMarkers()
        #proxy%>%addAwesomeMarkers(data=subsetData(),group="Imap")
        if(input$piechart){proxy%>%Strategy::addpiechartclustermarkers(.data=subsetData(),.colors=c("green","red","orange","purple","black"),group="Imap")}else{
          proxy%>%addMarkers(data = subsetData(),clusterOptions = markerClusterOptions())}
        #
      })
      
      #create a data object to display data
      
      output$data <-DT::renderDataTable(subsetData())
      
      
    }
    shinyApp(ui = ui, server = server)
  }
  

