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


ui <- navbarPage(theme = shinytheme("slate"),
                 title="Detect epidemic",
                 id="main",
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
                                                     choices = unique(Avo_fields$State))),
                            mainPanel(plotOutput('toto'),
                                      textOutput('totox'),
                                      leafletOutput("bbmap", height=1000)
                            )
                          )
                 ),
                 tabPanel(title="Data", DT::dataTableOutput("data")))