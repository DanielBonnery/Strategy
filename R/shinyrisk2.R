#' Shiny app - Simulation of infection and detection process
#' 
#' @param Fields a Spatial Data frame with polygons
#' @example 
#' 
library(Strategy)
library(avocado)
Shiny.FieldLevelType2risk(Donovan)

Shiny.FieldLevelType2risk <- function (Fields0=get(data(Fake,package="Strategy"))) {
require(shiny)
require(shinythemes)
require(leaflet)
require(ggplot2)
require(Strategy)
require(RLeafletTools)
ui <-function () {
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
      ))}
  server <- function(input , output,session) {
    library("sp")
    library("sf")
    library("rgdal")
    library(dplyr)
    library(leaflet)
    library(RLeafletTools)
    theme_dark2 = function() {
      theme_grey() %+replace%
        theme(
          # Specify axis options
          axis.line = element_blank(),  
          axis.text.x = element_text(color = "white",angle=90),  
          axis.text.y = element_text(color = "white"),  
          axis.ticks = element_line(color = "white"),  
          axis.title.x = element_text(color = "white"),  
          axis.title.y = element_text(color = "white"),
          # Specify legend options
          legend.background = element_rect(color = NA, fill = " gray10"),  
          legend.key = element_rect(color = "white",  fill = " gray10"),  
          legend.text = element_text(color = "white"),  
          legend.title = element_text(color = "white"),
          # Specify panel options
          panel.background = element_rect(fill = " gray10", color  =  NA),  
          panel.border = element_rect(fill = NA, color = "white"),  
          panel.grid.major = element_line(color = "grey35"),  
          panel.grid.minor = element_line(color = "grey20"),
          # Specify facetting options
          strip.background = element_rect(fill = "grey30", color = "grey10"),  
          strip.text.x = element_text(color = "white"),  
          strip.text.y = element_text(color = "white"),  
          # Specify plot options
          plot.background = element_rect(color = " gray10", fill = " gray10"),  
          plot.title = element_text(color = "white"),  
          plot.subtitle = element_text(color = "white"),  
          plot.caption = element_text( color = "white")
        )}
    #ggplot(data=data.frame(x=1:15,z=1,y=factor(1:3)),aes(x=x,y=z,color=y))+geom_point()
    th<-theme_dark2()
    theme_set(th)
    my_palette <- c('lightblue', 'red', 'white')
    names(my_palette)<-c('Original','Transformed',NA)
    assign("scale_colour_discrete", function(..., values = my_palette) scale_colour_manual(..., values = values), globalenv())
    assign("scale_fill_discrete", function(..., values = my_palette) scale_fill_manual(..., values = values), globalenv())
    #assign("scale_fill_ordinal", function(..., values = my_palette) scale_fill_manual(..., values = values), globalenv())
    #assign("scale_colour_ordinal", function(..., values = my_palette) scale_fill_manual(..., values = values), globalenv())
    colScale <- scale_colour_manual(name = "Origin",values = c('Original'='lightblue', 'Transformed'='red','white'))
    colScale2 <- scale_fill_manual(name = "Origin",values = c('Original'='lightblue', 'Transformed'='red','white'))
    
    #Risk with sampling with replacement approximation.
    approximatedrisk1<-function(N,m,n0,n1,n2,betac){(1-(m/N)*(1-betac^(n2)))^(n1*n0)}
    approximatedrisk2<-function(N,m,n0,n1,n2,betac){.5*(1-(m/N)*(1-betac^(n2)))^(n1*n0)}
    #for testing: 
    if(FALSE){
      input<-list(N=nrow(Trees),m=nrow(Trees),n0=2,n1=1,n2=1,betac=.5,alphac=0,samplemethod="srs",alpha=.1,beta=.1)
      attach(input)
      }
    drawLeaf<-function(Trees2,n2){
      Trees2$LeavesStatusSample<-"Healthy"
      Trees2$LeavesStatusSample[Trees2$HealthStatus=="Infected"]<-
        c("False negative","True positive")[2-rbinom(sum(Trees2$HealthStatus=="Infected"),
                                                     size=1,
                                                     prob=Trees2$beta[Trees2$HealthStatus=="Infected"]^n2)]
      Trees2$LeavesStatusSample[Trees2$HealthStatus=="Healthy"]<-
        c("True negative","False positive")[2-rbinom(sum(Trees2$HealthStatus=="Healthy"),
                                                     size=1,
                                                     prob=(1-Trees2$alpha[Trees2$HealthStatus=="Healthy"])^n2)]
      Trees2$LeavesStatusSample<-factor(Trees2$LeavesStatusSample,
                                        levels = c("True negative","False negative","False positive","True positive"),
                                        ordered=T)
      #levels(Trees2$LeavesStatusSample)
    Trees2$LeavesStatusSample
    }
    
    
    drawSystematic<-function(Trees2,n0,n1){
    Trees2$index<-1:nrow(Trees2)
    Trees2$Stratanames<-paste0(Trees2$Field,Trees2$Name)
    totalsamplesize=n0*n1
    ratios<-c(totalsamplesize*table(Trees2$Stratanames))
    allocation<-ratios%/%nrow(Trees2)
    prioritaires<-sort(ratios%%nrow(Trees2),decreasing = T)
    restants<-totalsamplesize-sum(allocation)
    if(restants>0){allocation[names(prioritaires)[1:restants]]<-allocation[names(prioritaires)[1:restants]]+1}
    sum(allocation)
    Sample<-plyr::adply(names(allocation)[allocation>0],1,function(stratum){
      NN<-nrow(Trees2[Trees2$Stratanames==stratum,])
      samplesize<-allocation[stratum]
      nbsamples<-NN%/%samplesize
      starter<-sample(1:nbsamples,1,T)
      Sample<-Trees2[Trees2$Stratanames==stratum,]$index[unique((((starter+(1:samplesize)*(nbsamples))-1)%%NN)+1)]
      data.frame(Sample=Sample)
    })$Sample
    Sample}
    
    drawHealthStatus<-function(Trees,m,foyers,.distriskhalf=10^(-4),delta=0.001){
      if(m<300){
        Trees$x=Trees$long
        Trees$y=Trees$lat
        infec<-which(Generate_Constrained_Epidemic(
          U=Trees,.distriskhalf=.distriskhalf,
                   jumprisk=0,
                   delta=delta,
                   numberinfected=m,
                   foyersaleatoires=foyers)$Iinf=="sick")
      }else{infec<-sample(nrow(Trees),min(m,nrow(Trees)),F)}
      status<-rep("Healthy",nrow(Trees))
      status[infec]<-"Infected"
      status
    }
    
    
    if(FALSE){
      Trees$x<-Trees$long
      Trees$y<-Trees$lat
      U<-Trees
      numberinfected<-m<-100;
      .distriskhalf=5*10^(-4)
      delta=.05
      foyersaleatoires<-foyers<-10;
      z<-drawHealthStatus(Trees2,m,foyers,.distriskhalf = 10^(-3),delta=.002)
      x=which(z=="Infected")
      ggplot2::ggplot(data.frame(x=x),aes(x=x))+geom_histogram(color="white")
      
    }

        drawSample<-function(Trees2, n0,n1,samplemethod){
    if(samplemethod=="srs"){
      Sample<-sample(1:nrow(Trees2),min(n0*n1,nrow(Trees2)),F)            
    }else{
      Sample<-drawSystematic(Trees2,n0,n1)
    }
    is.element(1:nrow(Trees2),Sample)}
    

        

    
    
    Trees.1000f<-function(TreesSim,n0=100,n1=1,n2=1,samplemethod="srs"){
      plyr::adply(
        plyr::raply(1000,
                  (function(){
                    TreesSim$Sample<-drawSample(TreesSim,n0,n1,samplemethod)
                    TreesSim$LeavesStatusSample<-drawLeaf(TreesSim,n2)
                    c("Number of samples that contain no leaf from an infected tree"      =!any(is.element(TreesSim[TreesSim$Sample,"LeavesStatusSample"],c("False negative","True positive"))),
                      "Number of samples that contain no positives bulks"      =!any(is.element(TreesSim[TreesSim$Sample,"LeavesStatusSample"],c("False positive","True positive"))),
                      "Number of samples that contain no true positives bulks" =!any(is.element(TreesSim[TreesSim$Sample,"LeavesStatusSample"],c("True positive"))))
                  })()),2,sum)}  
    
    Trees.f<-function(Trees2,m=100,n0=100,n1=1,n2=1,alpha=.1,beta=.1,foyers=1,samplemethod="srs"){
      
      Trees2$HealthStatus<-drawHealthStatus(Trees2,m,foyers)
      
      
     Trees2$beta<-beta
      Trees2$alpha<-alpha
      Trees2$LeavesStatusSample<-drawLeaf(Trees2,n2)
      
      Trees2$HealthStatus<-factor(Trees2$HealthStatus)
      Trees2["HealthStatusSample"]<-Trees2$HealthStatus
      Trees2$Sample<-drawSample(Trees2, n0,n1,samplemethod)
      Trees2}
    
    
    
    
    leafletmap<-function(TreesSim){
      leaflet(data=Fields0)%>%
        addPolygons(col=Fields0$color,fillColor = Fields0$color)%>% 
        addProviderTiles('Esri.WorldImagery',
                         options = providerTileOptions(minZoom = 1, maxZoom = 22,maxNativeZoom=18))%>%
        addProviderTiles("CartoDB.PositronOnlyLabels")%>%
        RLeafletTools::addpiechartclustermarkers2(.data=TreesSim,
                                   .colors=c("green","red"),
                                   group="HealthStatus",
                                   clusterId="toto2")%>%
        RLeafletTools::addpiechartclustermarkers2(.data=TreesSim[TreesSim$Sample,],
                                   .colors=c("green","red"),
                                   group="HealthStatusSample",
                                   clusterId="toto4")%>%
        RLeafletTools::addpiechartclustermarkers2(.data=TreesSim[TreesSim$Sample,],
                                   .colors=c("green","blue","orange","red"),
                                   group="LeavesStatusSample",
                                   clusterId="toto5")%>%
        addLayersControl(
          baseGroups = c("Clear",
                         "HealthStatus",
                         "HealthStatusSample",
                         "LeavesStatusSample"),
          options = layersControlOptions(collapsed = FALSE))%>%
        addLegend("bottomright", 
                  pal = colorFactor(c("blue","orange","green","red"), 
                                    domain = NULL), 
                  values = levels(TreesSim$LeavesStatusSample),
                  title = "Leaf level result",
                  opacity = .8
        )
    }
        
    
    # Panel 1
    observeEvent(input$N,  {
      updateSliderInput(session = session, inputId = "m", max = input$N)
      updateSliderInput(session = session, inputId = "n0", max = input$N%/%input$n1)
      #updateActionButton(session=session,inputId="rerun",label="Rerun")
    })
    
    observeEvent(input$n1,  {
      updateSliderInput(session = session, inputId = "n0", max = input$N%/%input$n1)
    })
    
    
    betaseq=seq(0,1,length.out=1000);
    riskbetac<-reactive({
      rbind(data.frame(Approximation=TRUE,
                       beta=input$betac,
                       risk=approximatedrisk1(N=input$N,m=input$m,n0=input$n0,n1=input$n1,n2=input$n2,betac=input$betac)),
            data.frame(Approximation=FALSE,
                       beta=input$betac,
                       risk=approximatedrisk2(N=input$N,m=input$m,n0=input$n0,n1=input$n1,n2=input$n2,betac=input$betac)))})
    dat.1 <- reactive({
      rbind(
        data.frame(
          Approximation=TRUE,
          beta=betaseq,
          risk=approximatedrisk1(N=input$N,m=input$m,n0=input$n0,n1=input$n1,n2=input$n2,betac=betaseq)),
        data.frame(
          Approximation=FALSE,
          beta=betaseq,
          risk=approximatedrisk2(N=input$N,m=input$m,n0=input$n0,n1=input$n1,n2=input$n2,betac=betaseq)))})
    output$plot.1<-renderPlot({
      ggplot(dat.1(),aes(x=beta,y=risk,color=Approximation,xintercept=beta,yintercept=risk))+
        geom_line()+
        geom_point(data=riskbetac())+
        geom_vline(data=riskbetac(),aes(xintercept=beta),color="white")+
        geom_hline(data=riskbetac(),aes(yintercept=risk,color=Approximation))+
        ylim(0,1)+colScale+colScale2},height = 400,width = 600)
    output$risk<-renderText({paste0("The risk is lower than ",signif(riskbetac()$risk,3))})

    
    Trees0<-reactive({
      Fields0$newcount<-round(input$N/sum(Fields0$Trees)*Fields0$Trees)
      Trees0<-plyr::adply(1:nrow(Fields0), 1, function(i) {
        X<-as.data.frame(sp::spsample(as(Fields0[i, ],"Spatial"), 
                                      as.data.frame(Fields0)[i, "newcount"], "regular"))
        names(X)<-c("long","lat")
        #SpatialPointsDataFrame(X,Fields0[i,setdiff(names(Fields0),"geometry")])
        merge(X,Fields0[i,setdiff(names(Fields0),"geometry")],all.x=T)})
      Trees0$Variety<-as.factor(Trees0$Variety)
      Trees0["VarietySample"]<-Trees0$Variety
      Trees0
    })
        
    Trees1.4<-eventReactive(input$refresh,{
      Trees.f(Trees0(),
              m=input$m,
              n0=input$n0,
              n1=input$n1,
              n2=input$n2,
              beta=input$betac,
              alpha=input$alphac,
              foyers=input$foyers,
              samplemethod = input$samplemethod)})
    
    
    Trees1000.1.4<-eventReactive(input$rerun,{
      Trees3<-Trees.f(Trees0(),
                      m=input$m,
              n0=input$n0,
              n1=input$n1,
              n2=input$n2,
              beta=input$betac,
              foyers=input$foyers,
              alpha=input$alphac,
              samplemethod = input$samplemethod)
      Trees.1000f(Trees3, n0=input$n0,n1=input$n1,n2=input$n2)
    })
    
        
    output$bbmap1.4<- renderLeaflet({leafletmap(Trees1.4())})
    
    output$counts<-renderTable(Trees1000.1.4())
    # Panel 2
  }
  runApp(list(ui = ui ,server = server))
}


