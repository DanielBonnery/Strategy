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
library(shiny)
library(dplyr)
library(leaflet)
library(DT)
library(Strategy)
 library(leaflet)
 library(dplyr)
 
ui <- navbarPage(theme = shinytheme("slate"),
                 title="Nice pie chart cluster markers",
                 id="main",
                 #                   tabPanel(title= "Interactive Map",
                 #                            leafletOutput("bbmap", height=1000)),
                 tabPanel("Interactive Map 1",
                          mainPanel(
                            leafletOutput("bbmap", height=1000))),
                 tabPanel("Interactive Map 2",
                          sidebarLayout(
                            sidebarPanel(width = 3, 
                                         sliderInput("time",
                                                     "Time (in days)",
                                                     min = 1,
                                                     max = 3,
                                                     value = 1)),
                            mainPanel(
                                      leafletOutput("bbmap2", height=1000)
                            )
                          )))
