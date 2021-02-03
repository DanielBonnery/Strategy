if(FALSE){library(rgdal)
library(raster)
library(mapview)
library(leaflet)
library(leafem)
library(Strategy)
library(plyr)
library("sp")
library(widgetframe) #for frameWidget
##############################
#1. Raw Data
##1.1. CLUM data
data(CLUM_Commodities_2018_v2,package="dataACLUMP")

##1.2. QLUMP data
QLUMP.filepath<-extract_DP_QLD_LANDUSE_June_2019()
QLUMP<-get(load(QLUMP.filepath))
QLUMP.sf<-QLUMP[!is.na(QLUMP$Commodity)&QLUMP$Commodity=="avocados",]
QLUMP.sf$Id<-1:nrow(QLUMP.sf)
#Donovan corresponds to #feature ID: 64554,63806

variety<-sf::read_sf(file.path(find.package("avocado"),'extdata','DonovanFarm.shp'))

#2. Harmonisation of CRS
QLUMP.sp<-as(QLUMP.sf,"Spatial")
QLUMP.sp$Id<-lapply(QLUMP.sp@polygons,function(l){l@ID})
QLUMP.sp4<-spTransform(QLUMP.sp,CRS("+proj=longlat +datum=WGS84"))
QLUMP.sp.Donovan<-QLUMP.sp4[is.element(QLUMP.sp$Id,c("64554","63806")),]
QLUMP.sf.Donovan<-st_as_sf(QLUMP.sp.Donovan)

raster::crs(spdf)
raster::crs(QLUMP.sf4)
raster::crs(QLUMP.sp4)

#st_intersection(QLUMP.sp4,spdf2)


#over(spdf2,QLUMP.sp4,returnList = T)
#over(spdf2,geometry(QLUMP.sp4))
#over(SpatialPoints(spdf2),QLUMP.sp4,returnList = T)
#what does not work:
#over(SpatialPoints(coords = rawpoints,
#                   proj4string  = proj4string(spdf2)),
#QLUMP.sp4,returnList = T)



Trees<-plyr::adply(1:nrow(variety), 1, function(i) {
  X<-as.data.frame(sp::spsample(as(variety[i, ],"Spatial"), 
                   as.data.frame(variety)[i, "Trees"], "regular"))
  names(X)<-c("long","lat")
  X$Name=variety$Name[i]
  X$Variety=variety$Variety[i]
  X})

#3. plots
variety$color=as.factor(variety$Variety)
levels(variety$color)=c("red","blue","black")
variety$color<-as.character(variety$color)
leaflet(data=variety)%>%addPolygons(col=variety$color,fillColor = variety$color)%>% 
  addProviderTiles('Esri.WorldImagery',
                   options = providerTileOptions(minZoom = 1, maxZoom = 21,maxNativeZoom=19))%>%
  addProviderTiles("CartoDB.PositronOnlyLabels")
pal <- colorFactor(c("navy", "red","black"), domain = c("Hass", "Shepard","Unknown"))
leaflet(data=variety)%>%addPolygons()%>% 
  addProviderTiles('Esri.WorldImagery',
                   options = providerTileOptions(minZoom = 1, maxZoom = 22,maxNativeZoom=18))%>%
  addProviderTiles("CartoDB.PositronOnlyLabels")%>%
  #addCircleMarkers(lng = Trees$x1,lat = Trees$x2,
  #addMarkers(lng = Trees$x1,lat = Trees$x2,
  #           clusterOptions = markerClusterOptions(col = ~pal(Trees$Variety)))
  addpiechartclustermarkers(.data=Trees,
                            .colors=c("green","blue"),
                            group="Variety")

}