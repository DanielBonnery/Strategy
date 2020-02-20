

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
                                                       max = 200,
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
      new_data <- bb_data[bb_data$ReportedDay<input$time,]
      return(new_data)
    })
    
    data(UE,package="Strategy")
    #bb_data<-bb_data[bb_data$ReportedDay<input$time,]
    # create a color paletter for category type in the data file
    
    pal <-
      colorNumeric(
        palette = "plasma",
        domain = bb_data$ReportedDay
      )
    popbins<-quantile(shapeData$population,(seq_len(11)-1)/10) 
    poppal <- colorBin(heat.colors(5), bins=popbins, na.color = "#aaff56",reverse = T)
    
    # create the leaflet map  
    output$bbmap <- renderLeaflet({
      leaflet(bb_data) %>% 
        addPolygons(data=shapeData,
                    stroke=TRUE,
                    weight=1,
                    color="black",
                    fillOpacity=.2,
                    fillColor=~poppal(shapeData$population)) %>% 
        addCircles(lng = Farms$Longitude, lat = Farms$Latitude, color = "black",fill = "white",opacity=5) %>% 
        addCircles(data = subsetData(),lng = ~Longitude, lat = ~Latitude, 
                   color = ~pal(ReportedDay)) %>% 
        addTiles() %>%
        addCircleMarkers(data = subsetData(), lat =  ~Latitude, lng =~Longitude, 
                         radius = 3, popup = ~as.character(cntnt), 
                         color = ~pal(ReportedDay),
                         stroke = FALSE, 
                         fillOpacity = 5)%>%
        addLegend(title="Reported Day",pal=pal, 
                  values=bb_data$ReportedDay,
                  opacity=1, 
                  na.label = "Not Available")%>%
        addLegend(title = "Population count", pal=poppal, 
                  values=shapeData$population,
                  opacity=1, 
                  na.label = "Not Available")%>%
        addEasyButton(easyButton(
          icon="fa-crosshairs", title="ME",
          onClick=JS("function(btn, map){ map.locate({setView: true}); }")))
    })
    
    
    #create a data object to display data
    
    output$data <-DT::renderDataTable(datatable(
      bb_data[-(4:5)],filter = 'top',
      colnames = c("X", "Y", "ReportedDay")
    ))
    
    
  }
  shinyApp(ui = ui, server = server)
}

test()