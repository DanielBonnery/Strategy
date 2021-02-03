#' Generate a population of avocado trees
#' 
#' @examples
#' data(U2,package="Strategy")
#' data(Avo_fields,package="Strategy")
#' data(U2,package="Strategy")
#' U2<-U2[2:4]
#' #Plot the trees
#' Avo_fields$Source_yr<-addNA(as.factor(Avo_fields$Source_yr))
#' QLD<-Avo_fields$State=="Qld"
#' Avo_ids<-unique(Avo_fields[QLD,]$Avo_id)[c(1,3,4)]
#' QLD<-QLD&is.element(Avo_fields$Avo_id,Avo_ids)
#' QLDt<-is.element(U2$id,Avo_ids)
#' yearpal <- colorFactor(heat.colors(5),domain = levels(Avo_fields$Source_yr),na.color = "#aaff56")
#' leaflet(Avo_fields[QLD,]) %>%
#' addProviderTiles('Esri.WorldImagery',options = providerTileOptions(minZoom = 1, maxZoom = 21,maxNativeZoom=19)) %>% 
#' addProviderTiles("CartoDB.PositronOnlyLabels")%>% 
#' addPolylines(fillOpacity = 1, weight = 3, smoothFactor = 0.5,opacity = 1.0,
#'             color=~yearpal(Avo_fields[QLD,]$Source_yr),
#'                         fillColor=~yearpal(Avo_fields[QLD,]$Source_yr))%>% 
#'                         addMarkers(lng = U2[QLDt,]$x1,lat = U2[QLDt,]$x2,clusterOptions = markerClusterOptions())


generate_U2<-function(){
  data(CLUM_Commodities_2018_v2,package="dataACLUMP")
  sh<-CLUM_Commodities_2018_v2[CLUM_Commodities_2018_v2$Commod_dsc=="avocados",]
  set.seed(1)
  p=.1+.9*exp(-sh$Area_ha)
  #plot(sh$Area_ha,p)
  sh$population<-2+rbinom(n=nrow(sh),size=round(200*sh$Area_ha/p),p)
  #summary(sh$population)
  sum(sh$population)
  #quantile(sh$population,.9+(0:10)/100)
  sh$Avo_id<-1:nrow(sh)
  U2<-Generate_U(sh,.id="Avo_id",.spatialobject="Area_ha",type="regular")
  return(list(U2=U2,sh=sh))}