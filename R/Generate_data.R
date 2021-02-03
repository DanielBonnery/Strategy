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
  plyr::adply(1:nrow(SpatialData),1,function(i){
    xx<-try(as.data.frame(sp::spsample(SpatialData[i,],as.data.frame(SpatialData)[i,"population"],type)))
    if(inherits(xx,"try-error")){xx<-NULL;print(SpatialData[i,"population"])}
    if(!is.null(.id)){xx["id"]=as.data.frame(SpatialData)[i,.id]}
    xx})} 


#' Risk to be infected by a neighbour at a distance x
#' 
#' @param .dist : a distance
#' @param .disthalfrisk : distance for which the risk is one half
#' @return  exp(-.dist/(log(2)*.distriskhalf))
#' @examples 
#' #Risk to be ingfected 2 m from the victim when the 50%risk distance is 1 m:
#' risktobeinfectedbydistancetooneinfectedunit(2,1)
risktobeinfectedbydistancetooneinfectedunit<-function(.dist,.distriskhalf=5*10^(-4)){
  exp(-.dist/(log(2)*.distriskhalf))
}

#' Risk to be infected by many neighbours a neighbour at a certain distance 
#' 
#' @param .dist : a vector (distances)
#' @param .disthalfrisk : distance for which the risk is one half
#' @param nI: total number of infected
#' @param jumprisk: probability to be infected by one person, no matter how far he(she) is
#' @return  1-(prod(1-risktobeinfectedbydistancetooneinfectedunit(.dist,.distriskhalf))*(1-jumprisk)^nI)
#' @examples 
#' #Risk to be ingfected 2 m from the victim when the 50%risk distance is 1 m:
#' risktobeinfectedbydistancetooneinfectedunit(2,1)

risktobeinfectedbydistancetoallinfectedunit<-function(.dist,nI,.distriskhalf=5*10^(-4),jumprisk=10^-6){
  1-(prod(1-risktobeinfectedbydistancetooneinfectedunit(.dist,.distriskhalf))*(1-jumprisk)^nI)
}





#' Computes the risk to be infected
#' 
#' @examples 
#' y=rep("Sane",nrow(U));y[sample(length(y),10)]<-"sick"
#' jumprisk=10^-6 
#' .distriskhalf=10^-6
risktobeinfected<-function(U,closedistances=NULL,sicks,new.sicks=NULL,.distriskhalf=5*10^(-4),jumprisk=10^-6,delta=0.01,
                           previouslyexposed=c(),previousrisk=NULL){
  stillfine<-setdiff(1:nrow(U),sicks)
  nI=length(sicks)
  #print(paste0(nI," sick."))
  newdistances=updatedist(closedistances=closedistances,U=U,sicks=sicks,new.sicks=new.sicks,delta=delta)
  #risk<-plyr::aaply(Matrix::sparseMatrix(i=closedistances$ind[,1],j=closedistances$ind[,2],x=closedistances$ra),1,risktobeinfectedbydistancetoallinfectedunit,nI=length(sicks))
  
  exposed=intersect(stillfine,unique(newdistances$ind[,1]))
  #print(paste0(length(exposed), " exposed."))
  newadditionalrisk=unlist(lapply(exposed,function(i){
    risktobeinfectedbydistancetoallinfectedunit(.dist = newdistances$ra[newdistances$ind[,1]==i],.distriskhalf =.distriskhalf, nI =nI )}))
  risk=newadditionalrisk
  if(!is.null(previousrisk)){risk[is.element(exposed,previouslyexposed)]=1-(1-previousrisk[is.element(previouslyexposed,exposed)])*(1-risk[is.element(exposed,previouslyexposed)])}
  return(list(exposed=exposed,risk=risk,newadditionalrisk=newadditionalrisk,closedistances=closedistances))
}

r<-function(risk){
  rbinom(nrow(y),size = 1,prob =risk )
}





#' Generate epidemic
#'
#' @param U a data.frame
#' @param TT an integer
#' @param .distriskhalf a positive number(default 5*10^(-4))
#' @param jumprisk =10^-6 a positive number
#' @param delta =0.05 a positive number
#' @examples 
#' .distriskhalf=5*10^(-4);jumprisk=10^-6;delta=0.05; TT=10
#' UE<-Generate_Discrete_Time_Epidemic(U,3)

Generate_Discrete_Time_Epidemic<-function(U,TT,.distriskhalf=5*10^(-4),jumprisk=10^-6,delta=0.05){
  y0<-paste0("I",formatC(0, width = 1+floor(log(TT)/log(10)), format = "d", flag = "0"))
  y1<-paste0("I",formatC(1, width = 1+floor(log(TT)/log(10)), format = "d", flag = "0"))
  U[[y0]]<-factor(c("sane","sick"))[1]
  U[[y1]]<-U[[y0]]
  U[[y1]][sample(nrow(U),10,replace=TRUE)]<-"sick"
  closedistances=NULL
  h<-neighbourhoods(U,delta)
  U$hexagon<-paste0(h@cID)
  dist_areas<-dist_areas_f(U,delta,h)
  for (tt in 2:TT){
    y<-paste0("I",formatC(tt, width = 1+floor(log(TT)/log(10)), format = "d", flag = "0"))
    y_1<-paste0("I",formatC(tt-1, width = 1+floor(log(TT)/log(10)), format = "d", flag = "0"))
    y_2<-paste0("I",formatC(tt-2, width = 1+floor(log(TT)/log(10)), format = "d", flag = "0"))
    sicks_2=(1:nrow(U))[U[[y_2]]=="sick"]
    sicks=sicks_1=(1:nrow(U))[U[[y_1]]=="sick"]
    new.sicks<-setdiff(sicks_1,sicks_2)
    U[[y]]<-U[[y_1]]
    R<-risktobeinfected(U,closedistances=closedistances,sicks=sicks_1,new.sicks=new.sicks,.distriskhalf=.distriskhalf,jumprisk=jumprisk,delta=delta)
    if(length(R$exposed)>0){
      closedistances<-R$closedistances
    contamination<-rbinom(length(R$exposed),size = 1,prob=R$risk)==1
    U[[y]][R$exposed][contamination]<-"sick"}
  }
  U
}




#' Generate epidemic under size constraint
#'
#' @param U a data.frame with x and y
#' @param .distriskhalf a positive number(default 5*10^(-4))
#' @param jumprisk =10^-6 a positive number
#' @param delta =0.05 a positive number
#' @param numberinfected target number of  infected elements of the population,
#' @param foyersaleatoires number of sources at random at the start
#' @examples 
#' .distriskhalf=5*10^(-4);jumprisk=10^-6;delta=0.05; TT=10
#' UE<-Generate_Discrete_Time_Epidemic(U,3)

Generate_Constrained_Epidemic<-
  function(U,.distriskhalf=5*10^(-4),
           jumprisk=10^-6,
           delta=0.05,
           numberinfected=10,
           foyersaleatoires=2){
  TT<-numberinfected
  y0<-paste0("I",formatC(0, width = 1+floor(log(TT)/log(10)), format = "d", flag = "0"))
  y1<-paste0("I",formatC(1, width = 1+floor(log(TT)/log(10)), format = "d", flag = "0"))
  U[[y0]]<-factor(c("sane","sick"))[1]
  U[[y1]]<-U[[y0]]
  U[[y1]][sample(nrow(U),foyersaleatoires,replace=F)]<-"sick"
  closedistances=NULL
  h<-neighbourhoods(U,delta)
  U$hexagon<-paste0(h@cID)
  dist_areas<-dist_areas_f(U,delta,h)
  conta=foyersaleatoires
  tt=1
  while (max(conta,tt)<numberinfected){
    tt<-tt+1
    y<-paste0("I",formatC(tt, width = 1+floor(log(TT)/log(10)), format = "d", flag = "0"))
    y_1<-paste0("I",formatC(tt-1, width = 1+floor(log(TT)/log(10)), format = "d", flag = "0"))
    y_2<-paste0("I",formatC(tt-2, width = 1+floor(log(TT)/log(10)), format = "d", flag = "0"))
    sicks_2=(1:nrow(U))[U[[y_2]]=="sick"]
    sicks=sicks_1=(1:nrow(U))[U[[y_1]]=="sick"]
    new.sicks<-setdiff(sicks_1,sicks_2)
    U[[y]]<-U[[y_1]]
    R<-risktobeinfected(U,closedistances=closedistances,sicks=sicks_1,new.sicks=new.sicks,.distriskhalf=.distriskhalf,jumprisk=jumprisk,delta=delta)
    if(length(R$exposed)>0){
      closedistances<-R$closedistances
      contamination<-rbinom(length(R$exposed),size = 1,prob=R$risk)==1
      if(sum(contamination)+sum(U[[y]]=="sick")>numberinfected){
        contamination=sample(length(R$exposed),size=numberinfected-sum(U[[y]]=="sick"),
                                                     prob=R$risk,replace=F)
      }
      if(sum(contamination)>0){U[[y]][R$exposed][contamination]<-"sick"}else{tt<-TT}
      }else{tt<-TT}
    conta<-sum(U[[y]]=="sick")
  }
  U[["Iinf"]]<-U[[y]]
  if(conta<numberinfected){U$Iinf[U$Iinf=="sane"][sample(sum(U$Iinf=="sane"),numberinfected-conta,F)]<-"sick"}
  U
}
