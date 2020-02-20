#' Generate spatial data that matches population counts
#' 
#' @param SpatialData : an object of class that includes 
#' @param type : argument to be passed to sp::spsample
#' @examples 
#' data(parish110217popest,package="dataONS")
#' data("mtcty150217population",package="dataONS")
#' shapeData2<-dataONS::dataParishes_December_2011_Boundaries_EW_BFC()
#' yy<-unique(get(data(Output_Area_to_Parish_to_Local_Authority_District_December_2011_Lookup_in_England_and_Wales,package="dataONS"))[c("PAR11CD","LAD11NM")])
#' names(yy)<-tolower(names(yy))
#' shapeData<-sp::merge(shapeData2,yy,by="par11cd",duplicateGeoms = TRUE)
#' parish110217popest2<-parish110217popest[
#'  is.element(parish110217popest$PAR11CD,
#'               shapeData$par11cd)&
#'                   parish110217popest$year=="mid_2006",
#'                     c("PAR11CD","Population")]
#'                     names(parish110217popest2)<-tolower(names(parish110217popest2))
#'                     shapeData=sp::merge(shapeData,parish110217popest2,by="par11cd",duplicateGeoms = TRUE)
#'                     shapeData$population[is.na(shapeData$population)]<-mean(shapeData$population,na.rm=TRUE)
#'                     shapeData<-subset(shapeData,is.element(lad11nm ,c("Allerdale", "Barrow-in-Furness", "Carlisle", "Copeland", "Eden","South Lakeland")))
#'
#' U<-Generate_U(shapeData,.id="par11cd",.spatialobject="st_areasha",type="random")
#' popbins<-quantile(shapeData$population,(seq_len(11)-1)/10) 
#' poppal <- colorBin(heat.colors(5), bins=popbins, na.color = "#aaff56",reverse = T)
#' library(leaflet)
#' 
#' leaflet(U) %>% 
#'   addPolygons(data=shapeData,
#'               stroke=TRUE,
#'               weight=1,
#'               color="black",
#'               fillOpacity=5,
#'               fillColor=~poppal(shapeData$population)) %>% 
#'   addTiles() %>% 
#'  addLegend(title = "Population count", pal=poppal, 
#'              values=shapeData$population,
#'               opacity=1, 
#'               na.label = "Not Available")

Generate_U<-function(SpatialData,.id=NULL,.spatialobject,type="random"){
  plyr::adply(1:nrow(shapeData),1,function(i){
    xx<-as.data.frame(sp::spsample(shapeData[i,.spatialobject],as.data.frame(shapeData)[i,"population"],type))
    if(!is.null(.id)){xx["id"]=as.data.frame(shapeData)[i,.id]}
    xx})} 



risktobeinfectedbydistancetooneinfectedunit<-function(.dist,.distriskhalf=5*10^(-4)){
  exp(-.dist/(log(2)*.distriskhalf))
}

risktobeinfectedbydistancetoallinfectedunit<-function(.dist,nI,.distriskhalf=5*10^(-4),jumprisk=10^-6){
  1-(prod(1-risktobeinfectedbydistancetooneinfectedunit(.dist,.distriskhalf))*(1-jumprisk)^nI)
}


library(fields)
#' @examples 
#' delta<-.005
#' sicks<-(1:nrow(U))[y=="sick"]
#' closedistances=updatedist(NULL,U,sicks)

updatedist<-function(closedistances=NULL,U,sicks,new.sicks=NULL,delta=0.005){
  if(is.null(closedistances)){
    closedistances=list(ind=matrix(NA,0,2),ra=vector())}
  
  if(is.null(new.sicks)){new.sicks<-setdiff(sicks,unique(closedistances$ind[,2]))}
  
  stillfine<-setdiff(1:nrow(U),sicks)
   if(length(new.sicks)>0){
    yy<-fields::fields.rdist.near(x1=U[stillfine,c("x","y"),drop=FALSE],
                                  x2=U[new.sicks,c("x","y"),drop=FALSE], delta=delta,max.points=length(stillfine)*length(new.sicks))
    closedistances$ind<-rbind(closedistances$ind,cbind(stillfine[yy$ind[,1]],new.sicks[yy$ind[,2]]))
    closedistances$ra=c(closedistances$ra,yy$ra)}
  keep=is.element(closedistances$ind[,1],stillfine)
  closedistances$ind<-closedistances$ind[keep,]
  closedistances$ra<-closedistances$ra[keep]
  closedistances}


#' @examples 
#' y=rep("Sane",nrow(U));y[sample(length(y),10)]<-"sick"
#' jumprisk=10^-6 
#' .distriskhalf=10^-6
risktobeinfected<-function(U,closedistances=NULL,sicks,new.sicks=NULL,.distriskhalf=5*10^(-4),jumprisk=10^-6,delta=0.01){
  stillfine<-setdiff(1:nrow(U),sicks)
  nI=length(sicks)
  closedistances=updatedist(closedistances=closedistances,U=U,sicks=sicks,new.sicks=new.sicks,delta=0.05)
  #risk<-plyr::aaply(Matrix::sparseMatrix(i=closedistances$ind[,1],j=closedistances$ind[,2],x=closedistances$ra),1,risktobeinfectedbydistancetoallinfectedunit,nI=length(sicks))
  exposed=intersect(stillfine,unique(closedistances$ind[,1]))
  risk=plyr::aaply(exposed,1,function(i){risktobeinfectedbydistancetoallinfectedunit(closedistances$ra[closedistances$ind[,1]==i],nI =nI )},.progress="text")
  return(list(exposed=exposed,risk=risk,closedistances=closedistances))
}

r<-function(risk){
  rbinom(nrow(y),size = 1,prob =risk )
}

#' @examples 
#' .distriskhalf=5*10^(-4);jumprisk=10^-6;delta=0.05 TT=10
#' UE<-Generate_Discrete_Time_Epidemic(U,3)

Generate_Discrete_Time_Epidemic<-function(U,TT,.distriskhalf=5*10^(-4),jumprisk=10^-6,delta=0.05){
  y0<-paste0("I",formatC(0, width = 1+floor(log(TT)/log(10)), format = "d", flag = "0"))
  y1<-paste0("I",formatC(1, width = 1+floor(log(TT)/log(10)), format = "d", flag = "0"))
  U[[y0]]<-factor(c("sane","sick"))[1]
  U[[y1]]<-U[[y0]]
  U[[y1]][sample(nrow(U),10)]<-"sick"
  closedistances=NULL
  for (tt in 2:TT){
    y<-paste0("I",formatC(tt, width = 1+floor(log(TT)/log(10)), format = "d", flag = "0"))
    y_1<-paste0("I",formatC(tt-1, width = 1+floor(log(TT)/log(10)), format = "d", flag = "0"))
    y_2<-paste0("I",formatC(tt-2, width = 1+floor(log(TT)/log(10)), format = "d", flag = "0"))
    sicks_2=(1:nrow(U))[U[[y_2]]=="sick"]
    sicks_1=(1:nrow(U))[U[[y_1]]=="sick"]
    new.sicks<-setdiff(sicks_1,sicks_2)
    U[[y]]<-U[[y_1]]
    R<-risktobeinfected(U,closedistances=closedistances,sicks=sicks,new.sicks=new.sicks,.distriskhalf=.distriskhalf,jumprisk=jumprisk,delta=delta)
    contamination<-rbinom(length(R$exposed),size = 1,prob=R$risk)==1
    U[[y]][R$exposed][contamination]<-"sick"
    save(U,file="U.rda")
  }
U
  }






Generate_O<-function(U,E,S){}


