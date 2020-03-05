library(Strategy)
data(Avo_fields,package="Strategy")
data(U2,package="Strategy")

U2<-U2[2:4]
#Plot the trees
Avo_fields$Source_yr<-addNA(as.factor(Avo_fields$Source_yr))

QLD<-Avo_fields$State=="Qld"
QLDt<-is.element(U2$id,Avo_fields$Avo_id[QLD])

yearpal <- colorFactor(heat.colors(5),domain = levels(Avo_fields$Source_yr),na.color = "#aaff56")

leaflet(Avo_fields[QLD,]) %>%
  addProviderTiles('Esri.WorldImagery',options = providerTileOptions(minZoom = 1, maxZoom = 21,maxNativeZoom=19)) %>% 
  addProviderTiles("CartoDB.PositronOnlyLabels")%>% 
  addPolylines(fillOpacity = 1, weight = 3, smoothFactor = 0.5,opacity = 1.0,
               color=~yearpal(Avo_fields[QLD,]$Source_yr),
               fillColor=~yearpal(Avo_fields[QLD,]$Source_yr))%>% 
  addMarkers(lng = U2[QLDt,]$x1,lat = U2[QLDt,]$x2,clusterOptions = markerClusterOptions())




if(FALSE){
data(UE)
library(ggplot2)
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

popbins<-quantile(shapeData$population,(seq_len(11)-1)/10) 
poppal <- colorBin(heat.colors(5), bins=popbins, na.color = "#aaff56",reverse = T)
library(leaflet)

whox<-UE$I002=="sick"

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
  addMarkers(lng = UE[whox,]$x,lat = UE[whox,]$y,
             clusterOptions = markerClusterOptions()
  )

}