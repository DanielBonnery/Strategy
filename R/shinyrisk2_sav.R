#' Shiny app - Simulation of infection and detection process
#' 
#' @param Fields a Spatial Data frame with polygons
#' @example 
#' DonovanFarm<-sf::read_sf(file.path(find.package("avocado"),'extdata','DonovanFarm.shp'))
#' DonovanFarm$color=as.factor(DonovanFarm$Variety)
#' levels(DonovanFarm$color)=c("green","blue","black")
#' DonovanFarm$color<-as.character(DonovanFarm$color)


#' Fields0<-DonovanFarm
#' source(file.path(Mydirectories::googledrive.teamdrive.directory(),"Avocado Sunblotch/Programs_and_papers/Strategy/R/shinyrisk2.R")); Shiny.FieldLevelType2risk(Fields0)

Shiny.FieldLevelType2risk <- function (Fields0=get(data(Fake,package="Strategy"))) {
require(shiny)
require(shinythemes)
require(leaflet)
require(ggplot2)
require(Strategy)
require(RLeafletTools)
require(plyr)
require(dplyr)
  
  
  
  ui <-function () {
    fluidPage(
      # Application title
      #setBackgroundColor("black"),
      navbarPage(
        "Type 2 Risk at the field level",id="main_navbar",theme = shinytheme("slate"),
        tabPanel(
          ###################### Tabpanel 1
          "Sample size computation and simulations",
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
                  sliderInput("r","r, infection rate (in %)",step=0.1,min = 0,  max = 10,  value = .5),
                  sliderInput("foyers","Infection sources",min = 1,  max = 10,  value = 1))),
              h4("Single leaf level risks"),
              fluidRow(
                column(
                  4,style=list("padding-right: 5px;"),
                  withMathJax(),
                  sliderInput("betacp","\\(\\beta\\), false negative (%)",min = 0,  max = 100,  value = 67,round=F)),
                column(
                  4,style=list("padding-right: 5px;"),
                  sliderInput("alphacp","False positive (%)",min = 0,  max = 100,  value = 0,round=F))),
              h4("Sampling parameters"),
              h5("Number of selected ..."),
              fluidRow(
                column(
                  4,style=list("padding-right: 5px;"),
                  sliderInput("n0",h5("...bulks (\\(n_0\\)),"),min = 1,  max = 300,  value = 120)),
                column(
                  4,style=list("padding-left: 5px;"),
                  sliderInput("n1",h5("...trees / bulk (\\(n_1\\)),"),min = 1,  max = 100,  value = 4)),
                column(
                  4,style=list("padding-left: 5px;"),
                  sliderInput("n2",h5("...leaves / tree (\\(n_2\\))."),min = 1,  max = 300,  value = 10))),
              h4("Risk requirements"),
              sliderInput("B",h5("B, maximal acceptable risk at the field level (%)"),min = 0,  max = 20,  value = 5),
              fluidRow(
                column(
                  6,style=list("padding-right: 5px;"),
                  radioButtons("samplemethod", 
                               h5("Selection method:"),
                               choices = list("Simple random sampling" = "srs", 
                                              "Stratified/systematic sampling" = "systematic",
                                              "Cluster sampling"="cluster"),
                               selected = "srs")),
                column(
                  6,style=list("padding-right: 5px;")))),
            mainPanel(
              tabsetPanel(
                type="tabs",
                tabPanel(
                  "Computation of theoretical risk and simulations",
                  h3("Theoretical risk of not detecting the viroid in the field"),
                  textOutput("risk"),
                  h3("Simulations"),
                  p("To display for the first time or refresh after changing control parameters, click 'Draw 1000 samples' button and wait for 1 minute."),
                  actionButton("rerun", "Draw 1000 samples",color="red"),
                  htmlOutput("tableparameters"),
                  htmlOutput("The lowest the percentages, the best the design"),
                  tableOutput("counts"),
                  h3("Risk as a function of the percentage of non infected leaves per tree"),
                  plotOutput("plot.1",height = 400)),
                tabPanel(
                  "Minimal sample size",
                  h3("Sample size"),
                  htmlOutput("samplesize3"),
                  h3("More leaves per tree options"),
                  htmlOutput("requirements1"),
                  tableOutput("samplesize4")),
                tabPanel("Visualisation",
                         p("To display for the first time or resimulate (with new or same parameters) the map, please click the 'Refresh map' button. "),
                         actionButton("refresh", "Refresh map"),
                         leafletOutput("bbmap1.4", height=600)),
                tabPanel(
                  "Biosecurity NZ",
                  h3("Sample size and sampling rate as a function of field size"),
                  p("Figure below shows that sampling rate and samping sizes as a function of the field size.
                    Those plotted sampling rate and samping sizes allow to detect with a probability of 95% the presence of the viroid when 0.5% of the trees of the field have a proportion of leaves without detectable positive RNA not greater than 67%.
                    The sampling rate decreases to 0 when the field size increases to infinity, whereas the sample size converges to a bounded number.
                    This means that the sample size to detect 0.5% of infected trees in a large field is the same than the sample size to detected 0.5% of infected trees in the whole Australia.
                    Biosecurity NZ does not impone an explicit limit on the field size in the document we consulted. However, the sampling rates are only given for fields up to around 5000 trees.
                    If Biosecurity NZ imposed that fields be divided in to 5000 trees size portions and that each portion be tested independently, this would mean that the rate of 14% would need to be applied to each portion."),
                  plotOutput("plot.2",height = 400),
                  h3("Official sampling rates"),
                  tableOutput("samplesize5"),
                  fluidRow(
                    column(
                      4,style=list("padding-right: 5px;"),
                      tableOutput("samplesize5.1")),
                    column(
                      4,style=list("padding-left: 5px;"),
                      tableOutput("samplesize5.2")),
                    column(
                      4,style=list("padding-left: 5px;"),
                      tableOutput("samplesize5.3"))),
                  h3("Parameters"),
                  withMathJax(),
                  htmlOutput("requirements2")
                  ),
                tabPanel(
                  "Formulae",
                  h3("Risk"),
                  withMathJax(),
                  p('The formula used to majorate the risk of not detecting when sampling with replacement or with systematic sampling when 
                           $$B(N,r,n_0,n_1,n_2,\\beta)=\\sum_{k=0}^{n_0\\times n_1} \\beta^{n_2\\times k} h(N,\\max(\\lfloor r\\times N\\rfloor,1),n_0\\times n_1,k),$$
                           where \\(\\lfloor.\\rfloor\\) is the ceiling function and \\(h(N,m,n,k)=\\left({{N \\choose n}}\\right)^{-1}{m \\choose k}{{N-m} \\choose {{n}-k}}\\) if \\( n+m-N\\leq k\\leq m\\), \\(0\\) otherwise.'),
                  h3("Sample size"),
                  withMathJax(),
                  p('The formula used to compute the minimum sample size to ensure that the risk is smaller than \\(B\\) is:
                           $$n(N,r,n_2,B)=\\min\\left\\{n\\leq N; B(N,r,n,1,n_2,\\beta)\\leq B\\right\\}.$$')))
        )))))}
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
  my_palette <- c('lightblue', 'red', 'white','green')
  names(my_palette)<-c('Original','Transformed',NA)
  assign("scale_colour_discrete", function(..., values = my_palette) scale_colour_manual(..., values = values), globalenv())
  assign("scale_fill_discrete", function(..., values = my_palette) scale_fill_manual(..., values = values), globalenv())
  #assign("scale_fill_ordinal", function(..., values = my_palette) scale_fill_manual(..., values = values), globalenv())
  #assign("scale_colour_ordinal", function(..., values = my_palette) scale_fill_manual(..., values = values), globalenv())
  #colScale <- scale_colour_manual(name = "Origin",values = c('Original'='lightblue', 'Transformed'='red','white'))
  #colScale2 <- scale_fill_manual(name = "Origin",values = c('Original'='lightblue', 'Transformed'='red','white'))
  
  #Risk with sampling with replacement approximation.
  approximatedrisk1<-function(N,m,n0,n1,n2,betac){(1-(m/N)*(1-betac^(n2)))^(n1*n0)}
  approximatedrisk2<-function(N,m,n0,n1,n2,betac){.5*(1-(m/N)*(1-betac^(n2)))^(n1*n0)}
  riskwithR        <-function(N,m,n0,n1,n2,betac){(1-(m/N)*(1-(betac)^(n2)))^(n1*n0)}
  riskwithoutR     <-function(N,m,n0,n1,n2,betac){
    possiblexs<-0:min(n0*n1,m)
    sum(betac^(n2*possiblexs)*dhyper(x=possiblexs,m=m,n=N-m,k=n0*n1))}
  
  requirednumberoftreeswithR<-function(N,m,n1,n2,beta,risk){
    ceiling(log(risk)/(log(1-(m/N)*(1-beta^n2))))
  }
  
  #'@param N population size
  #'@param m number of infected
  #'@param n0 number of bulks in sample
  #'@param n1 number of trees in a bulk
  #'@param n2 number of colected leaves per sampled tree
  #'@param beta risk at the leaf level, or 1-proportion of leaves with detectable RNA in an infected tree 
  #'@return an integer, number of required leaves
  #'@examples
  #' requirednumberoftreeswithR(1237,6,1,1,0.63,0.05)
  #' requirednumberoftreeswithoutR(1237,6,1,1,0.63,0.05)
  #' riskwithoutR(1237,6,1236,1,0.63)
  #' requirednumberoftreeswithoutR(1237,6,1,10,0.63,0.05)
  #' requirednumberoftreeswithoutR(1237,6,1,10,0,0.05)
  requirednumberoftreeswithoutR<-function(N,m,n1,n2,beta,risk){
    n<-0
    risk2<-1
    while(n<N&risk2>risk){
      n<-n+1
      risk2<-riskwithoutR(N,m,n0=1,n,n2,beta)}
    n}
  
  
  
  
  
      #for testing: 
    if(FALSE){
      TreesSim<-Trees3
      input<-list(N=nrow(Trees),
                  r=.5,
                  foyers=1,
                  betacp=67,
                  alphacp=0,
                  m=nrow(Trees),
                  n0=2,n1=1,n2=1,
                  betac=.5,alphac=0,
                  samplemethod="srs",
                  alpha=.1,
                  beta=.1)
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

    drawSystematic2<-function(Trees2,n0,n1){
      NN<-nrow(Trees2)
      (sample(NN,1)+(1:(n0*n1))*(NN%/%(n0*n1)))%%NN+1
      }
    
    
    drawcluster2<-function(Trees2,n0,n1){
      NN<-nrow(Trees2)
      (sample(NN,1)+(1:(n0*n1)))%%NN+1}
    
        
    drawcluster<-function(Trees2,n0,n1){
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
        starter<-sample(NN,1,T)
        Sample<-Trees2[Trees2$Stratanames==stratum,]$index[unique((((starter+(1:samplesize))-1)%%NN)+1)]
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
    }
    if(samplemethod=="systematic"){
            Sample<-drawSystematic2(Trees2,n0,n1)            
          }
    if(samplemethod=="cluster"){
          Sample<-drawcluster2(Trees2,n0,n1)
    }
    is.element(1:nrow(Trees2),Sample)}
    

        

    
    
    Trees.1000f<-function(TreesSim,n0=100,n1=1,n2=1,samplemethod="srs",nsim=1000){
      X<-plyr::aaply(
        plyr::raply(nsim,
                  (function(){
                    TreesSim$Samplesrs<-drawSample(TreesSim,n0,n1,"srs")
                    TreesSim$Samplesys<-drawSample(TreesSim,n0,n1,"systematic")
                    TreesSim$Samplecluster  <-drawSample(TreesSim,n0,n1,"cluster")
                    TreesSim$LeavesStatusSample<-drawLeaf(TreesSim,n2)
                    X<-cbind(
                      srs=c(!any(is.element(TreesSim[TreesSim$Samplesrs,"LeavesStatusSample"],c("False positive","True positive"))),
                      !any(is.element(TreesSim[TreesSim$Samplesrs,"LeavesStatusSample"],c("False negative","True positive"))),
                      !any(is.element(TreesSim[TreesSim$Samplesrs,"LeavesStatusSample"],c("True positive")))),
                      systematic=c(
                        !any(is.element(TreesSim[TreesSim$Samplesys,"LeavesStatusSample"],c("False positive","True positive"))),
                        !any(is.element(TreesSim[TreesSim$Samplesys,"LeavesStatusSample"],c("False negative","True positive"))),
                        !any(is.element(TreesSim[TreesSim$Samplesys,"LeavesStatusSample"],c("True positive")))),
                      cluster=c(
                        !any(is.element(TreesSim[TreesSim$Samplecluster,"LeavesStatusSample"],c("False positive","True positive"))),
                        !any(is.element(TreesSim[TreesSim$Samplecluster,"LeavesStatusSample"],c("False negative","True positive"))),
                        !any(is.element(TreesSim[TreesSim$Samplecluster,"LeavesStatusSample"],c("True positive")))))
                    rownames(X)<-c("Samples for which all bulks tested negative (Type 2 error)",
                                                        "Samples not containing any leaf from an infected tree",
                                                        "Samples not containing any true positives bulks")
                    X})()),2:3,sum)
      X<-plyr::aaply(100*X/nsim,1:2,paste,"%")
      X<-data.frame(rownames(X),X)
      rownames(X)<-NULL
      names(X)<-c("","Simple random Sampling","Systematic","Cluster sampling")
      X
      }  
    
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
        RLeafletTools::addpiechartclustermarkers2(
                                   .data=TreesSim[TreesSim$Sample,],
                                   .colors=c("green","blue","orange","red"),
                                   group="LeavesStatusSample",
                                   clusterId="toto5")%>%
        addLayersControl(
          baseGroups = c("Clear",
                         "HealthStatus",
                         "HealthStatusSample",
                         "LeavesStatusSample"),
          options = layersControlOptions(collapsed = FALSE))%>%
        leaflet::addLegend("bottomright", 
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
    
    samplesize<-reactive({
      samplesize2=  requirednumberoftreeswithoutR(N=input$N,
                                                 m=ceiling(input$N*input$r/100),
                                                 n1=input$n1,
                                                 n2=input$n2,
                                                 beta=(input$betacp/100),
                                                 risk=(input$B/100))
      data.frame(samplesize=samplesize2,
                 bulks=ceiling(samplesize2/input$n1))})

    
    
    tablesamplesize<-reactive({
      n2<-1:10
      tabsamplesize=  sapply(n2,function(nx2){requirednumberoftreeswithoutR(N=input$N,
                                                  m=ceiling(input$N*input$r/100),
                                                  n1=input$n1,
                                                  n2=nx2,
                                                  beta=(input$betacp/100),
                                                  risk=(input$B/100))})
      n1<-round(50/(1:10))
      X<-data.frame(x=1:10,
                    n1=as.integer(n1),
                    n1n2=as.integer(n1*n2),
                    samplesize=as.integer(tabsamplesize),
                 bulks=as.integer(ceiling(tabsamplesize/(n1))),
                 leaves=as.integer(tabsamplesize*n2))
      names(X)<-c("Number of leaves per tree","Trees per bulk","Leaves per bulk","Required number of trees","Required number of bulks","Required total number of leaves")
      X
      })
    
    
        
    riskbetac<-reactive({
      data.frame(     beta=(input$betacp/100),
                       risk=riskwithoutR(N=input$N,
                                              m=ceiling(input$N*input$r/100),
                                              n0=input$n0,
                                              n1=input$n1,
                                              n2=input$n2,
                                              betac=(input$betacp/100)))})
    dat.1 <- reactive({
        data.frame(
          beta=betaseq,
          risk=sapply(betaseq,function(beta){riskwithoutR(N=input$N,m=ceiling(input$N*input$r/100),n0=input$n0,n1=input$n1,n2=input$n2,betac=beta)}))})
    output$plot.1<-renderPlot({
      ggplot(dat.1(),aes(x=beta,y=risk,xintercept=beta,yintercept=risk))+
        geom_line(color="red")+
        geom_point(data=riskbetac(),color="white")+
        geom_vline(data=riskbetac(),aes(xintercept=beta),color="white")+
        geom_hline(data=riskbetac(),aes(yintercept=risk),color="white")+
        ylim(0,1)+xlab("beta: maximal proportion of leaves WITHOUT detectable NRA in infected trees")},height = 400,width = 600)
        output$toto<-renderText({paste0("The risk of not detecting the desease is lower than ",signif(riskbetac()$risk,3)," (approximately ",round(riskbetac()$risk,3)*100,"%)")})
        output$risk<-renderText({paste0("The risk of not detecting the desease is lower than ",signif(riskbetac()$risk,3)," (approximately ",round(riskbetac()$risk,3)*100,"%)")})
    output$samplesize3<-renderUI({HTML(paste0("For the risk of no detection at the field level to be inferior to <b>",input$B,"%</b>, when<ul>
                                              <li>",input$n2," leaves are collected from each sampled tree,</li>
                                              <li> each bulk is made of ",input$n1*input$n2," leaves coming from ",input$n1, " distinct trees,</li>
                                              <li> material from the same tree does not go into different bulks,</li>
                                              <li> At least ",input$r,"% of the trees have a proportion of leaves without detactable positive RNA NOT EXCEEDING <b>",input$betacp,"%</b></li>
                                              </ul>
                                              <p style='font-size:2ex;'>the required number of trees is <b>",
                                           samplesize()$samplesize,".</b></p><br/>
                                           This corresponds to <b>",ceiling(samplesize()$samplesize/input$n1),"</b> bulks of ",input$n1*input$n2," leaves each, and <b>",samplesize()$samplesize*input$n2,"</b> leaves in total."))})
    output$samplesize4<-renderTable({tablesamplesize()})
    
    
    output$requirements1<-renderUI({HTML(paste0("For the risk of no detection at the field level to be inferior to ",input$B,"%, when
                                               at least ",input$r,"% of the trees have a proportion of leaves without detactable positive RNA NOT EXCEEDING ",input$betacp,"%, the table below gives different options in function of the number of leaves sampled from each sampled tree."))})
    
    output$requirements2<-renderUI({
      withMathJax(HTML(paste0("Biosecurity New Zealand standards requires that <ul>
                                              <li>10 leaves are collected from each sampled tree</li>
                                              <li> each bulk is made of 40 leaves coming from 4 distinct trees</li>
                                              <li> material from the same tree does not go into different bulks</li>
                                              </ul> 
The number of trees to sample is such that there is less than 5% chance of not detecting the viroid when at least 0.5% of the trees have a proportion of leaves without detactable positive RNA NOT EXCEEDING 67%.</br></br>
                                         Where does the number 0.67 come from ?</br>
There is an ambiguity in the Biosecurity NZ document. Biosecurity NZ explains that 0.98% chance of detection at the tree level with 10 leaves  corresponds to \\(\\beta=0.63\\). 
However the relationship between the confidence of 0.98 and \\(\\beta\\) is the following: $$\\mathrm{confidence}=\\left(1-\\beta^{10}\\right).$$ 
If we solve the equation \\(0.98=\\left(1-\\beta^{10}\\right)\\), the solution we obtain is  \\(\\beta=\\exp(\\log(0.02)/10)=",round(exp(log(0.02)/10),5),"\\), so approximately ", round(100*exp(log(0.02)/10)),"%. 
If we compute \\(\\left(1-\\beta^{10}\\right)\\) with \\(\\beta=0.63\\), we get \\(\\mathrm{confidence}=",round(1-(.63^10),5),"\\). 
To get a confidence of exactly 99%, the coefficient \\(\\beta\\) must be equal to  \\(\\exp(\\log(0.01)/10)=",round(exp(log(0.01)/10),5),"\\), 
which is close to \\(\\beta=0.63\\).
                                                To be conservative, we use the threshold ",100*round(exp(log(0.02)/10),5),"%, or even 80%.")))})
    
    
    Trees0.f<-reactive({
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
      Trees.f(Trees0.f(),
              m=ceiling(input$N*input$r/100),
              n0=input$n0,
              n1=input$n1,
              n2=input$n2,
              beta=(input$betacp/100),
              alpha=(input$alphacp/100),
              foyers=input$foyers,
              samplemethod = input$samplemethod)})
    
    
    Trees1000.1.4<-eventReactive(input$rerun,{
      Trees3<-Trees.f(Trees0.f(),
                      m=ceiling(input$N*input$r/100),
              n0=input$n0,
              n1=input$n1,
              n2=input$n2,
              beta=(input$betacp/100),
              foyers=input$foyers,
              alpha=(input$alphacp/100),
              samplemethod = input$samplemethod)
      Trees.1000f(Trees3, n0=input$n0,n1=input$n1,n2=input$n2)
    })
    
    tableparamslag<-eventReactive(input$rerun,HTML(
    paste0("
    Parameters used for the table<br/><ul> 
               <li>Trees in the field: ",input$N,"</li>
               <li>Trees in the sample :",input$n0*input$n1," (",input$n0," bulks x ",input$n1," trees/bulk)</li>
               <li>Number and percentage of infected trees: ",floor(input$r*input$N/100)," (",input$r,"% of the total)</li>
               <li>Infected trees contain exactly : ",input$betacp,"% of  leaves WITHOUT detectable RNA.</li>")))
    
    output$bbmap1.4<- renderLeaflet({leafletmap(Trees1.4())})
    
    output$counts<-renderTable(Trees1000.1.4())

      output$tableparameters<-renderUI(tableparamslag())

      bionz<-data.frame(
        countmax=c(206,412,618,824,1030,1236,1443,1649,1855,2061,2267,2473,2679,2885,3092,3298,3504,3710,3916,4122,4535),
        samplingrate=c(1,.96,.78,.64,.53,.45,.40,.35,
                       .32,.29,.26,.24,.23,.21,.20,
                       .19,.18,.17,.16,.15,.14))%>%
          mutate(      
            countmin=c(1,countmax[-length(countmax)]+1),
            infectedtrees=pmax(1,floor(.005*countmin)),
            noninfectedtrees=countmin-floor(.005*countmin),
            "Sampling rate"=paste0(100*samplingrate,"%"),
            "Field size interval"=paste0("[",countmin,",",countmax,"]"))%>%
          select("Field size interval","Sampling rate")
      
      output$samplesize5.1<-renderTable({bionz[1:7,]})
      output$samplesize5.2<-renderTable({bionz[8:14,]})
      output$samplesize5.3<-renderTable({bionz[15:21,]})
      plot2<-function(){
        
        table2<-data.frame(
          countmax=c(206,412,618,824,1030,1236,1443,1649,1855,2061,2267,2473,2679,2885,3092,3298,3504,3710,3916,4122,4535,Inf),
          samplingrate=c(1,.96,.78,.64,.53,.45,.40,.35,
                         .32,.29,.26,.24,.23,.21,.20,
                         .19,.18,.17,.16,.15,.14,NA))%>%
          dplyr::mutate(      
            countmin=c(1,countmax[-length(countmax)]+1),
            infectedtrees=pmax(1,floor(.005*countmin)),
            noninfectedtrees=countmin-floor(.005*countmin),
            samplesize=ceiling(samplingrate*countmin))
        
        
        table2.2<-plyr::adply(table2[1:(nrow(table2)-1),],
                              1,
                              function(d){
                                data.frame(count=d$countmin:d$countmax)
                              })%>%
          dplyr::mutate(
            recommendedsamplesize=ceiling(samplingrate*count),
            infectedtrees=pmax(1,floor(.005*count)),
            noninfectedtrees=count-infectedtrees)
        
        
        table2.2<-plyr::adply(table2.2,
                              1,
                              function(d){
                                data.frame(optimalsamplesize=requirednumberoftreeswithoutR(d$count,d$infectedtrees,n1=NA,n2=10,beta=exp(log(0.02)/10),risk=0.05))
                              })
        table2.2<-table2.2%>%
          mutate(optimalsamplingrate=optimalsamplesize/count)
        
        table2.3<-table2.2%>%
          reshape2::melt(id.vars="count",
                         measure.vars=c("optimalsamplesize",
                                        "recommendedsamplesize",
                                        "samplingrate",
                                        "infectedtrees",
                                        "optimalsamplingrate"))%>%
          mutate(type1=
                   recode(variable,"optimalsamplesize"="Sample sizes",
                          "recommendedsamplesize"="Sample sizes",
                          "samplingrate"="Sampling rates",
                          "infectedtrees"="Detectable infected trees",
                          "optimalsamplingrate"="Sampling rates"),
                 Type=recode(variable,
                             "optimalsamplesize"="Optimal",
                             "recommendedsamplesize"="Recommended by MAF",
                             "samplingrate"="Recommended by MAF",
                             "infectedtrees"="Infected trees",
                             "optimalsamplingrate"="Optimal"))
        ggplot(data=table2.3,
               aes(x=count,y=value,
                   group=variable,colour=Type))+
          geom_line()+
          theme(legend.position = "bottom")+
          xlab("Number of trees in the field")+
          ylab("Sample size")+
          facet_wrap(~type1, scales = "free")+
          ggtitle("Biosecurity New Zealand Surveillance Standard ",
                  "Sample size vs number of trees in the field, to achieve a risk of 5% with 10 leaves per tree")+ 
          labs(fill="")+scale_colour_discrete(name="")
      }
      
      
      plot2.0<-plot2()  
      

      output$plot.2<-renderPlot({plot2.0})
    # Panel 2
  }
  runApp(list(ui = ui ,server = server))
}


