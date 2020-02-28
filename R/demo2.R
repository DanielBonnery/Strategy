if(FALSE){# Direct application of https://nceas.github.io/oss-lessons/spatial-data-gis-law/3-mon-intro-gis-in-r.html

#Load all the libraries
library(rgdal)
library(raster)
library(ggplot2)
library(rgeos)
library(mapview)
library(leaflet)
library(broom) # if you plot with ggplot and need to turn sp data into dataframes
options(stringsAsFactors = FALSE)
#Declare tempfile and temp directories where the data will be saved
tmpf  <-tempfile()
tmpd  <-tempdir()
#Download and unzip
unzip(system.file("inst/extdata","shapefile_currency_clum_50m1218m.zip",package="dataACLUMP"),exdir=tmpd)
x=file.path("~/daniel.bonnery@gmail.com/Telecargaments/Spatial_data_zip/CLUM_Commodities_2018_v2.shp")


aus.sf<-sf::read_sf(x)
aus.df<-data.frame(aus.sf)
# load the data 
Rshp <- readOGR(x,verbose=TRUE)
unique(Rshp$Commod_dsc)
Rshpavocado<-Rshp[Rshp$Commod_dsc=="avocados",]

Rshptree<-Rshp[is.element(Rshp$Tertiary,c("5.1.1 Production nurseries","3.4.1 Tree fruits","4.4.1 Irrigated tree fruits")),]
nrow(Rshptree)
class(Rshp)
extent(Rshp)
crs(Rshp)
#plot(Rshp, main = "ACLUMP")
library(ggmap)

leaflet(Rshpavocado) %>%
 addProviderTiles('Esri.WorldImagery') %>% 
  addProviderTiles("CartoDB.PositronOnlyLabels")%>% 
  addPolylines(fillOpacity = 1,color = "red", weight = 3, smoothFactor = 0.5,
               opacity = 1.0)

sum(area(Rshpavocado))
sum(Rshpavocado$Area_ha)
nrow(Rshpavocado)

nrow(Rshptree)

leaflet(Rshptree) %>%
  addProviderTiles('Esri.WorldImagery') %>% 
  addProviderTiles("CartoDB.PositronOnlyLabels")%>% 
  addPolylines(fillOpacity = 1,color = ~Commod_dsc, weight = 3, smoothFactor = 0.5,
               opacity = 1.0)

joliepalette<-scales::hue_pal()(nlevels(as.factor(Rshptree$Commod_dsc)))

map <-  leaflet(Rshptree) %>%
  addProviderTiles('Esri.WorldImagery') %>% 
  addProviderTiles("CartoDB.PositronOnlyLabels") 
  
Rshptree$col=joliepalette[as.factor(Rshptree$Commod_dsc)]

for( group in levels(as.factor(Rshptree$Commod_dsc))){
  map <- addPolylines(map,data=Rshptree[as.factor(Rshptree$Commod_dsc)==group,], color=~col)
}
map

#Now loading Queensland only.
y=file.path("~/daniel.bonnery@gmail.com/Telecargaments/QLD_LANDUSE_June_2019/QLD_LANDUSE_June_2019.gdb")
ql <- readOGR(dsn = y,layer = "QLD_LANDUSE_CURRENT_X" ,require_geomType="wkbPolygon",drop_unsupported_fields = T,dropNULLGeometries = T)

ogrListLayers(y)

rgdal::ogrInfo(y, layer =  "QLD_LANDUSE_CURRENT_X",require_geomType="wkbPolygon")
}