---
title: "Generic functions for plant epidemiology"
output: pdf_document
editor_options: 
  chunk_output_type: console
always_allow_html: yes
---

# Simulations for desease detection


```
## Error in library(Strategy): there is no package called 'Strategy'
```

```
## Error in library(RLeafletTools): there is no package called 'RLeafletTools'
```


## Installation
To install  the package, execute:

```r
cred <- git2r::cred_user_pass("username", "password")
devtools::install_git(
  "http://pine--is.grid.private.cam.ac.uk:8888/gitlab/dbb31/Strategy.git", 
  credentials = cred)
```
where '"username"' and '"password"' need to be replaced by correct strings.





## Generate a tree population that matches the ACLUMP data:

The following script creates a tree population based on the ACLUMP data. 

```r
data(CLUM_Commodities_2018_v2, package = "dataACLUMP")
    sh <- CLUM_Commodities_2018_v2[CLUM_Commodities_2018_v2$Commod_dsc == 
        "avocados", ]
    set.seed(1)
    p = 0.1 + 0.9 * exp(-sh$Area_ha)
    sh$population <- 2 + rbinom(n = nrow(sh), size = round(200 * 
        sh$Area_ha/p), p)
    sum(sh$population)
    sh$Avo_id <- 1:nrow(sh)
    U2 <- Generate_U(sh, .id = "Avo_id", .spatialobject = "Area_ha", 
        type = "regular")
```

To visualize the data, one can use code like the following (only three fields out of 797 are selected:


```r
library(Strategy)
```

```
## Error in library(Strategy): there is no package called 'Strategy'
```

```r
data(Avo_fields,package="Strategy")
```

```
## Error in find.package(package, lib.loc, verbose = verbose): there is no package called 'Strategy'
```

```r
data(U2,package="Strategy")
```

```
## Error in find.package(package, lib.loc, verbose = verbose): there is no package called 'Strategy'
```

```r
U2<-U2[2:4]
```

```
## Error in eval(expr, envir, enclos): object 'U2' not found
```

```r
#Plot the trees
Avo_fields$Source_yr<-addNA(as.factor(Avo_fields$Source_yr))
```

```
## Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'as.factor': object 'Avo_fields' not found
```

```r
QLD<-Avo_fields$State=="Qld"
```

```
## Error in eval(expr, envir, enclos): object 'Avo_fields' not found
```

```r
Avo_ids<-unique(Avo_fields[QLD,]$Avo_id)[c(1,3,4)]
```

```
## Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'unique': object 'Avo_fields' not found
```

```r
QLD<-QLD&is.element(Avo_fields$Avo_id,Avo_ids)
```

```
## Error in eval(expr, envir, enclos): object 'QLD' not found
```

```r
QLDt<-is.element(U2$id,Avo_ids)
```

```
## Error in is.element(U2$id, Avo_ids): object 'U2' not found
```

```r
yearpal <- colorFactor(heat.colors(5),
                       domain = levels(Avo_fields$Source_yr),
                       na.color = "#aaff56")
```

```
## Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'levels': object 'Avo_fields' not found
```

```r
leaflet(Avo_fields[QLD,]) %>%
  addProviderTiles('Esri.WorldImagery',
                   options = providerTileOptions(minZoom = 1, 
                                                 maxZoom = 21,
                                                 maxNativeZoom=19)) %>% 
  addProviderTiles("CartoDB.PositronOnlyLabels")%>% 
  addPolylines(fillOpacity = 1, weight = 3, smoothFactor = 0.5,opacity = 1.0,
               color=~yearpal(Avo_fields[QLD,]$Source_yr),
               fillColor=~yearpal(Avo_fields[QLD,]$Source_yr))%>% 
  addMarkers(lng = U2[QLDt,]$x1,
             lat = U2[QLDt,]$x2,
             clusterOptions = markerClusterOptions())
```

```
## Error in structure(list(options = options), leafletData = data): object 'Avo_fields' not found
```

We have 797 avocado fields in Queensland from ACLUMP.
The desease is spreading by root, so we can consider that root to root transmission is not possible for trees separated by more than the maximum avocado root reach. I read that the avocado should be planted at least 30 feet from houses, so I think it is safe to say that fields separated by more than 200 m will not contaminate each other via their roots.
So we can simulate the in-field spread independantly for each of these fields, which means we can parallelize computations.
We then need to compute the distances between each fields.
I wrote some functions that compute distances between polygons.


The ACLUMp data uses the following projection system:

```r
data(Avo_fields)
sp::proj4string(Avo_fields)
```

```
## [1] "+proj=longlat +ellps=GRS80 +no_defs"
```


The following code converts the longitude latitude values using the projection EPSG:3107 (GDA94/SA Lambert).

```r
projAvo_fields<-sp::spTransform(Avo_fields,CRS("+init=epsg:3107 +units=m"))
```



The following code extracts all the projected polygons.


```r
list.poly<-Strategy::extractpolygonsaslist(projAvo_fields)
```


The following code gets the list of couple of fields and their  distances.


```r
Australia_avo_fields_small_dist_matrix<-
  Australia_avo_fields_small_dist_matrix_f(delta=Inf)
```


From this list, one can establish a clustering of the different fields, in which simulation of root-to-root infection can be done independently.


```r
Australia_connected_fields200<-Australia_connected_fields_f(delta=Inf)
```

We print the number of fields and the number of bins:

```r
data(Australia_connected_fields200)
lapply(Australia_connected_fields200,function(x){length(unique(x))})
```

```
## $polygon
## [1] 797
## 
## $bin
## [1] 388
```
This means we can parallelize the simulation over 388 clusters.

This is what we do when we simulate an epidemy with


```r
U2E2<-generate_Avo_Epi_1()
```
The result on a selection of fields can be seen by executing



```r
library(Strategy)
```

```
## Error in library(Strategy): there is no package called 'Strategy'
```

```r
data(Avo_fields,package="Strategy")
```

```
## Error in find.package(package, lib.loc, verbose = verbose): there is no package called 'Strategy'
```

```r
data(U2E2,package="Strategy")
```

```
## Error in find.package(package, lib.loc, verbose = verbose): there is no package called 'Strategy'
```

```r
U2E2<-U2E2[c(2:4,12:112)]
```

```
## Error in eval(expr, envir, enclos): object 'U2E2' not found
```

```r
names(U2E2)[1:2]<-c("long","lat")
```

```
## Error in names(U2E2)[1:2] <- c("long", "lat"): object 'U2E2' not found
```

```r
#Plot the trees
Avo_fields$Source_yr<-addNA(as.factor(Avo_fields$Source_yr))

QLD<-Avo_fields$State=="Qld"

Avo_ids<-unique(Avo_fields[QLD,]$Avo_id)[c(1:10)]
QLD<-QLD&is.element(Avo_fields$Avo_id,Avo_ids)
QLDt<-is.element(U2$id,Avo_ids)
```

```
## Error in is.element(U2$id, Avo_ids): object 'U2' not found
```

```r
yearpal <- colorFactor(heat.colors(5),
                       domain = levels(Avo_fields$Source_yr),
                       na.color = "#aaff56")

leaflet(Avo_fields[QLD,]) %>%
  addProviderTiles('Esri.WorldImagery',
                   options = providerTileOptions(minZoom = 1, 
                                                 maxZoom = 21,maxNativeZoom=19)) %>% 
  addProviderTiles("CartoDB.PositronOnlyLabels")%>% 
  addPolylines(fillOpacity = 1, weight = 3, smoothFactor = 0.5,opacity = 1.0,
               color=~yearpal(Avo_fields[QLD,]$Source_yr),
               fillColor=~yearpal(Avo_fields[QLD,]$Source_yr))%>%
 addpiechartclustermarkers(.data=U2E2[is.element(U2E2$id,Avo_ids),]
                           ,.colors=c("green","red","orange","purple","black"),
                           group="I012")
```

```
## Error in addpiechartclustermarkers(., .data = U2E2[is.element(U2E2$id, : could not find function "addpiechartclustermarkers"
```



```r
library(Strategy)
```

```
## Error in library(Strategy): there is no package called 'Strategy'
```

```r
data(Avo_fields,package="Strategy")
```

```
## Error in find.package(package, lib.loc, verbose = verbose): there is no package called 'Strategy'
```

```r
data(U2E2,package="Strategy")
```

```
## Error in find.package(package, lib.loc, verbose = verbose): there is no package called 'Strategy'
```

```r
U2E2<-U2E2[c(2:4,12:112)]
```

```
## Error in eval(expr, envir, enclos): object 'U2E2' not found
```

```r
names(U2E2)[1:2]<-c("long","lat")
```

```
## Error in names(U2E2)[1:2] <- c("long", "lat"): object 'U2E2' not found
```

```r
#Plot the trees
Avo_fields$Source_yr<-addNA(as.factor(Avo_fields$Source_yr))

QLD<-Avo_fields$State=="Qld"

Avo_ids<-unique(Avo_fields[QLD,]$Avo_id)[c(1:10)]
QLD<-QLD&is.element(Avo_fields$Avo_id,Avo_ids)
QLDt<-is.element(U2$id,Avo_ids)
```

```
## Error in is.element(U2$id, Avo_ids): object 'U2' not found
```

```r
yearpal <- colorFactor(heat.colors(5),domain = levels(Avo_fields$Source_yr),na.color = "#aaff56")

leaflet(Avo_fields[QLD,]) %>%
  addProviderTiles('Esri.WorldImagery',
                   options = providerTileOptions(minZoom = 1, 
                                                 maxZoom = 21,maxNativeZoom=19)) %>% 
  addProviderTiles("CartoDB.PositronOnlyLabels")%>% 
  addPolylines(fillOpacity = 1, weight = 3, smoothFactor = 0.5,opacity = 1.0,
               color=~yearpal(Avo_fields[QLD,]$Source_yr),
               fillColor=~yearpal(Avo_fields[QLD,]$Source_yr))%>%
 addpiechartclustermarkers(.data=U2E2[is.element(U2E2$id,Avo_ids),]
                           ,.colors=c("green","red","orange","purple","black"),
                           group="I085")
```

```
## Error in addpiechartclustermarkers(., .data = U2E2[is.element(U2E2$id, : could not find function "addpiechartclustermarkers"
```


## Details and code demo



### Data

I wrote a package `dataACLUMP` a contain functions (`dataACLUMP::get_data_from_web`) to download the ACLUMP data and create data files.

After installing the package dataACLUMP from gitlab, the last commodities data files can be obtained by executing

```r
data(CLUM_Commodities_2018_v2, package = "dataACLUMP")
```




### Plotting the epidemic
This function 'addpiechartclustermarkers' allows to replace cluster markers by cluster pie charts on leaflet.

```r
library(Strategy)
```

```
## Error in library(Strategy): there is no package called 'Strategy'
```

```r
data("breweries91",package="leaflet")
breweries91$goodbear<-sample(as.factor(c("terrific","marvelous","culparterretaping")),
                             nrow(breweries91),
                             replace=T)
leaflet(breweries91) %>%
 addTiles() %>%
 addpiechartclustermarkers(.data=breweries91,.colors=c("red","green","blue"),group="goodbear")
```

```
## Error in addpiechartclustermarkers(., .data = breweries91, .colors = c("red", : could not find function "addpiechartclustermarkers"
```


### Computing distances between polygons

The data from ACLUMP contains longitudes and latitudes of vertices of fields, given as polygons.
To compute distances between two points identified with their longitudes or latitudes, one can 
use native R functions from the packages rgdal or .
To compute the minimum distance between two fields, I had to write down my own functions.
The distance between two polygons can be computed by the minimum distance separating the different edges of the two polygons.
The different functions described in this section allow to compute this distance.


A function of two triangles have the same orientation.


```r
example(triangleorientation,echo=F)
```

A function to tell if two segments intersect.


```r
example(segment.intersect,echo=F)
```




```r
#example(Generate_Discrete_Time_Epidemics,echo=F)
```



```r
example(ranges.gap,echo=F)
```



```r
example(rangesoverlap,echo=F)
```


```r
example(Generate_U,echo=F)
```




```r
example(risktobeinfectedbydistancetooneinfectedunit)
```



`projpointonseg_a` is a function that gives the position
of the closest point of a segment to a specific point on the plane. The position is the ratio of the distance of the closest point to one extremity of the segment divided by the length of the segment.

```r
example(projpointonseg_a,echo=F)
```


```r
example(closestpointonpolygon,echo=F)
```


```r
example(closestpointsontwosegments,echo=F)
```


```r
example(closestpointsontwosegments_n,echo=F)
```


```r
example(closestpointsontwopolygons,echo=F)
```


```r
example(closestpointsontwopolygons_n,echo=F)
```


```r
example(distpointtoseg,echo=F)
```




```r
example(distpointtopoly,echo=F)
```



`distpointtoseg` is a function that returns the distance between a point and a segment.

```r
example(distpointtoseg,echo=F)
```


```r
example(segment.intersect,echo=F)
```

```r
example(distsegmenttosegment,echo=F)
```


```r
example(distsegmenttopoly)
```


```r
example(distpolytopoly,ask=F)
```


```r
example(polydistmat,echo=F)
```


```r
example(polysmalldistmat,echo=F)
```

```r
example(connectedpop)
```


```r
example(extractpolygonsaslist)
```



```r
example(neighbourhoods,echo=F)
```


```r
example(dist_areas_f)
```


```r
example(newdist)
```


```r
example(updatedist)
```
