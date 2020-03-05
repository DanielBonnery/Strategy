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

#' hexagonal bins
#' 
#' @param U : a dataframe containing the numerical variables x and y 
#' @param delta: controls the bin diameter
#' @return a hexbin object  hexagonal bins
#' @examples 
#' # plot the hex bins of cumbria
#' data(U)
#' plot(neighbourhoods(U,.002))
neighbourhoods<-function(U,delta=(range(U$x,na.rm=TRUE)[2]-range(U$x,na.rm=TRUE)[1])/100){
  hexbin::hexbin(U$x, U$y, xbins =ceiling((range(U$x,na.rm=TRUE)[2]-range(U$x,na.rm=TRUE)[1])/2*delta), xbnds=range(U$x,na.rm=TRUE),ybnds= range(U$y,na.rm=TRUE), IDs = TRUE)}

#' Distances between hexagonal bins 
#' 
#' @param U : a dataframe containing the numerical variables x and y and preferable hexagon
#' @param delta : needed if hexagon is not a variable of U: bins will be recomputed
#' @return a named matrix
#' @examples 
#' data(U) 
#' dist_areas_f(U)
#' 
#' delta<-0.01
#' h<-neighbourhoods(U,delta)
#' U$hexagon<-paste0(h@cID)
#' hD<-dist_areas_f(U,h)
#'
#'sss1=sample(nrow(U),1000)
#'sss2=sample(nrow(U),1000)
#'x=sapply(1:1000,function(i){dist(U[c(sss1[i],sss2[i]),c("x","y")])})
#'y<-sapply(1:1000,function(i){hD[U$hexagon[sss1[i]],U$hexagon[sss2[i]]]})
#'plot(x,y,pch=".")

dist_areas_f<-function(U,delta=(range(U$x)[2]-range(U$x)[1])/100,h=neighbourhoods(U,delta)){
  if(is.null(U$hexagon)){U$hexagon<-h@cID}
  dd<-plyr::ddply(U,~hexagon, function(d){data.frame(mx=median(d$x),my=median(d$y))})
  centers<-hexbin::hcell2xy(h)
  dd$x<-sapply(dd$mx,function(y,x){x[which.min(abs(x - y))]},x=unique(centers$x))
  dd$y<-sapply(dd$my,function(y,x){x[which.min(abs(x - y))]},x=unique(centers$y))
  hD<-as.matrix(dist(dd[c("x","y")]))
  #ggplot(dd, aes(x, y)) +geom_segment(aes(x = x, y = y, xend = mx, yend = my, colour = "segment"), data = dd)
  dimnames(hD)<-list(dd$hexagon,dd$hexagon)
  hD
}




#' @examples 
#' delta<-.005
#' sicks<-(1:nrow(U))[y=="sick"]
#' closedistances=updatedist(NULL,U,sicks)

updatedist<-function(closedistances=NULL,U,sicks,new.sicks=NULL,delta=0.005,dist_areas=dist_areas_f(U,delta)){
  if(is.null(closedistances)){
    closedistances=list(ind=matrix(NA,0,2),ra=vector())}
  
  if(is.null(new.sicks)){new.sicks<-setdiff(sicks,unique(closedistances$ind[,2]))}
  
  stillfine<-setdiff(1:nrow(U),sicks)
  if(length(new.sicks)>0){
    new.sickhexagons<-unique(U$hexagon[new.sicks])
    
    L<- parallel::mclapply(new.sickhexagons,function(hexagon){
      exposedhexagons<-dimnames(dist_areas)[[1]][dist_areas[hexagon,]<delta]
      new.sickinhexagon<-new.sicks[is.element(U$hexagon[new.sicks],hexagon)]
      stillfineexposedhexagon<-stillfine[is.element(U$hexagon[stillfine],exposedhexagons)]
      if(length(stillfineexposedhexagon)*length(new.sickinhexagon)>0){
        yy<-fields::fields.rdist.near(x1=U[stillfineexposedhexagon,c("x","y"),drop=FALSE],
                                      x2=U[new.sickinhexagon,c("x","y"),drop=FALSE], delta=delta,max.points=length(stillfineexposedhexagon)*length(new.sickinhexagon))
        return(if(!identical(yy$ra,-1)){
          list(ind=cbind(stillfine[yy$ind[,1]],new.sicks[yy$ind[,2]]),ra=yy$ra)}else{list(ind=NULL,ra=NULL)})
      }})
    closedistances$ra=c(closedistances$ra,do.call(c,lapply(L,function(x){x$ra})))
    closedistances$ind<-rbind(closedistances$ind,do.call(rbind,lapply(L,function(x){x$ind})))
    keep=is.element(closedistances$ind[,1],stillfine)
    closedistances$ind<-closedistances$ind[keep,]
    closedistances$ra<-closedistances$ra[keep]}
  closedistances}


#' @examples 
#' y=rep("Sane",nrow(U));y[sample(length(y),10)]<-"sick"
#' jumprisk=10^-6 
#' .distriskhalf=10^-6
risktobeinfected<-function(U,closedistances=NULL,sicks,new.sicks=NULL,.distriskhalf=5*10^(-4),jumprisk=10^-6,delta=0.01){
  stillfine<-setdiff(1:nrow(U),sicks)
  nI=length(sicks)
  print(paste0(nI," sick."))
  closedistances=updatedist(closedistances=closedistances,U=U,sicks=sicks,new.sicks=new.sicks,delta=delta)
  #risk<-plyr::aaply(Matrix::sparseMatrix(i=closedistances$ind[,1],j=closedistances$ind[,2],x=closedistances$ra),1,risktobeinfectedbydistancetoallinfectedunit,nI=length(sicks))
  exposed=intersect(stillfine,unique(closedistances$ind[,1]))
  print(paste0(length(exposed), " exposed."))
  risk=unlist(parallel::mclapply(exposed,function(i){risktobeinfectedbydistancetoallinfectedunit(closedistances$ra[closedistances$ind[,1]==i],nI =nI )},.progress="text"))
  return(list(exposed=exposed,risk=risk,closedistances=closedistances))
}

r<-function(risk){
  rbinom(nrow(y),size = 1,prob =risk )
}


#' @examples 
#' .distriskhalf=5*10^(-4);jumprisk=10^-6;delta=0.05; TT=10
#' UE<-Generate_Discrete_Time_Epidemic(U,3)

Generate_Discrete_Time_Epidemic<-function(U,TT,.distriskhalf=5*10^(-4),jumprisk=10^-6,delta=0.05){
  y0<-paste0("I",formatC(0, width = 1+floor(log(TT)/log(10)), format = "d", flag = "0"))
  y1<-paste0("I",formatC(1, width = 1+floor(log(TT)/log(10)), format = "d", flag = "0"))
  U[[y0]]<-factor(c("sane","sick"))[1]
  U[[y1]]<-U[[y0]]
  U[[y1]][sample(nrow(U),10)]<-"sick"
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
    contamination<-rbinom(length(R$exposed),size = 1,prob=R$risk)==1
    U[[y]][R$exposed][contamination]<-"sick"
    if(tt%%10==0){save(U,file="UE.rda")}
  }
  U
}








