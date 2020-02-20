library(Strategy)
data(U)
data(UE)

data(parish110217popest,package="dataONS")
data("mtcty150217population",package="dataONS")
shapeData2<-dataONS::dataParishes_December_2011_Boundaries_EW_BFC()
yy<-unique(get(data(Output_Area_to_Parish_to_Local_Authority_District_December_2011_Lookup_in_England_and_Wales,package="dataONS"))[c("PAR11CD","LAD11NM")])
names(yy)<-tolower(names(yy))
shapeData<-sp::merge(shapeData2,yy,by="par11cd",duplicateGeoms = TRUE)
parish110217popest2<-parish110217popest[
  is.element(parish110217popest$PAR11CD,
             shapeData$par11cd)&
    parish110217popest$year=="mid_2006",
  c("PAR11CD","Population")]
names(parish110217popest2)<-tolower(names(parish110217popest2))
shapeData=sp::merge(shapeData,parish110217popest2,by="par11cd",duplicateGeoms = TRUE)
shapeData$population[is.na(shapeData$population)]<-mean(shapeData$population,na.rm=TRUE)
shapeData<-subset(shapeData,is.element(lad11nm ,c("Allerdale", "Barrow-in-Furness", "Carlisle", "Copeland", "Eden","South Lakeland")))



#' runCompare
#' 
#' @description Shiny App to 
#' @export
#' @examples
#' package1<-NULL
#' package2<-NULL
#' runCompare()
test<-function(){
  library(dplyr)
  library(gridExtra)
  library(ggplot2)
  library(shinythemes)
  
  library(shiny)
  
  library(dplyr)
  
  library(leaflet)
  
  library(DT)
  
  
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
                                                       min = 1,
                                                       max = ncol(UE)-4,
                                                       value = 30),
                                           selectInput("testapproach", 
                                                       "Test", 
                                                       c("Frequentist","Bayesian")) 
                              ),
                              
                              mainPanel(
                                leafletOutput("bbmap", height=1000)
                                
                              )
                            )
                   ),
                   tabPanel(title="Data", DT::dataTableOutput("data")))
  
  
  server <- function(input, output) {
    
    #library(datagovuk)
    library("sp")
    library("rgdal")
    library(dplyr)
    library(leaflet)
    
    
    subsetData <- reactive({
      return(UE[UE[,input$time+4]=="sick",])
    })
    
    data(UE,package="Strategy")
    #bb_data<-bb_data[bb_data$ReportedDay<input$time,]
    # create a color paletter for category type in the data file
    
    popbins<-quantile(shapeData$population,(seq_len(11)-1)/10) 
    poppal <- colorBin(heat.colors(5), bins=popbins, na.color = "#aaff56",reverse = T)
    popbins<-quantile(shapeData$population,(seq_len(11)-1)/10) 
    poppal <- colorBin(heat.colors(5), bins=popbins, na.color = "#aaff56",reverse = T)
    library(leaflet)
    
    output$bbmap <- renderLeaflet({
    leaflet(UE) %>% 
      addPolygons(data=shapeData,
                  stroke=TRUE,
                  weight=1,
                  color="black",
                  fillOpacity=5,
                  fillColor=~poppal(shapeData$population)) %>% 
      addTiles() %>% 
      addLegend(title = "Population count", pal=poppal, 
                values=shapeData$population,
                opacity=1, 
                na.label = "Not Available") %>% 
      addMarkers(data = subsetData(),lng = ~x, lat = ~y,
                 clusterOptions = markerClusterOptions())})
    
    
    #create a data object to display data
    
    output$data <-DT::renderDataTable(datatable(UE))
    
    
  }
  shinyApp(ui = ui, server = server)
}

test()