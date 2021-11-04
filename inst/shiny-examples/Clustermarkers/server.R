server <- function(input, output) {
  
  #library(datagovuk)
  library("sp")
  library("rgdal")
  library(dplyr)
  library(leaflet)
  
  
  #y<-paste0("I",formatC(input$time, width =3, format = "d", flag = "0"))
  #D<-if(input$Status=="Pie-chart"){Trees}else{Trees[is.element(Trees[[y]],input$Status),]}
  #D$II<-D[[y]]
  #D
  
  #bb_data<-bb_data[bb_data$ReportedDay<input$time,]
  # create a color paletter for category type in the data file
  
  #  popbins<-quantile(shapeData$population,(seq_len(11)-1)/10) 
  # poppal <- colorBin(heat.colors(5), bins=popbins, na.color = "#aaff56",reverse = T)
  #  popbins<-quantile(shapeData$population,(seq_len(11)-1)/10) 
  #  poppal <- colorBin(heat.colors(5), bins=popbins, na.color = "#aaff56",reverse = T)
  library(leaflet)
  library(Strategy)
  data("breweries91",package="leaflet")
  breweries91.2<-breweries91
  breweries91.2$goodbeer<-sample(as.factor(c("terrific","marvelous","culparterretaping")),nrow(breweries91),replace=T)
  beer=do.call(data.frame,plyr::rlply(3,sample(as.factor(c("terrific","marvelous","culparterretaping")),nrow(breweries91),replace=T)))
  names(beer)<-paste0("beer",1:3);
  breweries91.3<-data.frame(breweries91, 
                            beer)
  
  variableattime <- reactive({
    y<-paste0("beer",input$time)
  })
  
  output$bbmap2 <- renderLeaflet({
    leaflet(breweries91.3) %>%
      addTiles() %>%
      addpiechartclustermarkers(
        .data=breweries91.3,
        .colors=c("red","green","blue"),
        group=paste0("beer",input$time))
  })

  
  output$bbmap3 <- renderLeaflet({
    leaflet(breweries91.3) %>%
      addTiles() %>%
      addpiechartclustermarkers(
        .data=breweries91.3,
        .colors=c("red","green","blue"),
        group=paste0("beer1"))
  })
  
 # observeEvent(input$time,{
#        leafletProxy("bbmap2")%>%
#          clearMarkerClusters()%>%
#          clearMarkers()},priority=2)

  observeEvent(input$time,{
    vv=paste0("beer",input$time)
    leafletProxy("bbmap3")%>%
      #clearMarkerClusters()%>%
      #clearMarkers()%>%
      addpiechartclustermarkers(
        .data=breweries91.3,
        .colors=c("orange","purple","black","green","red"),
        group=vv)%>%
      addpiechartclustermarkers(
        .data=breweries91.3,
        .colors=c("yellow","gray","pink","green","red"),
        group=vv)})

  
  
  output$bbmap <- renderLeaflet({
    leaflet(breweries91.2) %>%
      addTiles() %>%
      addpiechartclustermarkers(.data=breweries91.2,
                                .colors=c("red","green","blue"),
                                group="goodbeer")
  })
  
}