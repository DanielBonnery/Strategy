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