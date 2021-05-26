#' Create hexagons to cover an area
#'
#' @param U a spatial object  of class sf.
#' @param delta = a positive number:length in meters of the side of smallest hexagon.
#' @examples 
#' library(sf)
#' .spatialdata<-st_as_sf(LocustAnalysis::Swarms)
#' hexagonize(.spatialdata) 
hexagonize<-function(.spatialdata,hexagondiameterinangle=5){
  centroid<-st_coordinates(st_centroid(st_combine(.spatialdata)))
  A=matrix(c(3/2,-sqrt(3)/4,3/2,sqrt(3)/4),2,2)
  projmat<-solve(t(A)%*%A)%*%t(A)
  alpha<-(function(A,.scale=hexagondiameterinangle){
    round(t(projmat%*%(t(A)-c(centroid))/.scale))})(st_coordinates(.spatialdata)[,1:2])
  colnames(alpha)<-c("alpha1","alpha2")
  cbind(.spatialdata,alpha)}


triangularise<-function(.spatialdata,trianglebasis=5){
  centroid<-st_coordinates(st_centroid(st_combine(.spatialdata)))
  A=matrix(c(3/2,-sqrt(3)/4,3/2,sqrt(3)/4),2,2)
  projmat<-solve(t(A)%*%A)%*%t(A)
  alpha<-(function(A,.scale=trianglebasis){
    round(t(projmat%*%(t(A)-c(centroid))/.scale))})(st_coordinates(.spatialdata)[,1:2])
  colnames(alpha)<-c("alpha1","alpha2")
  cbind(.spatialdata,alpha)}


#'@example 
#'squarise(LHSwarms)
squarise<-function(.spatialdata,squarebasis=5){
  alpha<-squarebasis*round(st_coordinates(.spatialdata)/squarebasis)
  squareid<-plyr::aaply(alpha,1,paste,collapse="-")
  squares<-data.frame(squareid=squareid,alpha,stringsAsFactors = F)
  uniquesquares<-unique(squares)
  squares<-st_sf(cbind(data.frame(id=1:nrow(.spatialdata)),squares),geom=st_geometry(.spatialdata))
  uniquesquares<-st_as_sf(as.data.frame(uniquesquares),coords=intersect(c("X","Y","Z"),names(uniquesquares)), crs = sf::st_crs(.spatialdata))
  list(squares=squares,
       uniquesquares=uniquesquares)}
#'@examples
#'sparse_st_distances(.spatialdata1,.spatialdata2)
sparse_st_distances<-
  function(.spatialdata1,
           .spatialdata2=NULL,
           maxdistfilter=100000,
           angleforsquares=maxdistfilter/111000,
           id1=NULL,
           id2=NULL){
    if(is.null(.spatialdata2)){
      squares1<-squarise(.spatialdata1,squarebasis=angleforsquares)
      squares1$squaresid<-1:nrow(.spatialdata1)
      distm<-st_distance(squares1$uniquesquares)
      distm<-units::drop_units(distm)
      dimnames(distm)<-list(1:length(squares1$uniquesquares$squareid),
                            1:length(squares1$uniquesquares$squareid))
      distm<-reshape2::melt(distm)
      distm<-distm[distm$value<2*maxdistfilter&distm$Var1<=distm$Var2,]
      distm$Var1<-squares1$uniquesquares$squareid[distm$Var1]
      distm$Var2<-squares1$uniquesquares$squareid[distm$Var2]
      pointdistm<-plyr::ddply(distm,~Var1+Var2,function(d){
        dd<-st_distance(x = squares1$squares[squares1$squares$squareid==d$Var1,],
                        y=squares1$squares[squares1$squares$squareid==d$Var2,])
        dd<-units::drop_units(dd)
        dd<-reshape2::melt(dd)
        dd$Var1<-squares1$squares$id[squares1$squares$squareid==d$Var1]
        if(d$Var1!=d$Var2){dd}else{dd[dd$Var1<dd$Var2,]}
      },.progress="text")
      if(!is.null(id1)){
        pointdistm$Var1<-.spatialdata[[id1]][pointdistm$Var1];
        pointdistm$Var2<-.spatialdata[[id1]][pointdistm$Var2]}
      }else{
        squares1<-squarise(.spatialdata1,squarebasis=angleforsquares)
        squares2<-squarise(.spatialdata2,squarebasis=angleforsquares)
        
        squares1$squaresid<-1:nrow(.spatialdata1)
        squares2$squaresid<-1:nrow(.spatialdata2)
        distm<-st_distance(squares1$uniquesquares,squares2$uniquesquares)
        distm<-units::drop_units(distm)
        dimnames(distm)<-list(1:length(squares1$uniquesquares$squareid),
                              1:length(squares2$uniquesquares$squareid))
        distm<-reshape2::melt(distm)
        distm<-distm[distm$value<2*maxdistfilter,]
        distm$Var1<-squares1$uniquesquares$squareid[distm$Var1]
        distm$Var2<-squares2$uniquesquares$squareid[distm$Var2]
        pointdistm<-plyr::ddply(distm,~Var1+Var2,function(d){
          dd<-st_distance(x = squares1$squares[squares1$squares$squareid==d$Var1,],
                          y=squares2$squares[squares2$squares$squareid==d$Var2,])
          dd<-units::drop_units(dd)
          dd<-reshape2::melt(dd)
          if(nrow(dd)>0){
            dd$Var1<-squares1$squares$id[squares1$squares$squareid==d$Var1]
            dd$Var2<-squares2$squares$id[squares2$squares$squareid==d$Var2]}
          dd},.progress="text")
        if(!is.null(id1)){pointdistm$Var1<-.spatialdata[[id1]][pointdistm$Var1]}
        if(!is.null(id2)){pointdistm$Var2<-.spatialdata[[id2]][pointdistm$Var2]}}
    pointdistm[pointdistm$value<maxdistfilter,]}





  
hexagonizecheck<-function(.spatialdata,hexagondiameterinangle=5){
  if(!is.element("sf",class(.spatialdata))){.spatialdata<-st_as_sf(.spatialdata)}
  #creates envolope of x
  envelope<-sf::st_convex_hull(st_combine(.spatialdata))
  #leaflet(envelope)%>%addTiles()%>%leafem::addFeatures()
  rangex<-range(st_coordinates(envelope)[,"X"])
  rangey<-range(st_coordinates(envelope)[,"Y"])
  distancebetweencentroidsonsamerow<-hexagondiameterinangle*3/2
  distancebetweencentroidsonsamecolumn<-hexagondiameterinangle*sqrt(3)/2
  ncentroidperrow<-(rangex[2]-rangex[1])%/%distancebetweencentroidsonsamerow+2
  ncentroidpercolumn<-((rangey[2]-rangey[1])%/%distancebetweencentroidsonsamecolumn)+1
  hexagoncentroidevenline.x<-rangex[1]+((-1):(ncentroidperrow-1))*distancebetweencentroidsonsamerow
  hexagoncentroidoddline.x<-hexagoncentroidevenline.x+distancebetweencentroidsonsamerow/2
  hexagoncentroidevenline.y<-rangey[1]+(-1:ncentroidpercolumn-1)*distancebetweencentroidsonsamecolumn
  hexagoncentroidoddline.y<-hexagoncentroidevenline.y+(distancebetweencentroidsonsamecolumn/2)
  hexcentroids.mat<-rbind(expand.grid(
    X=hexagoncentroidevenline.x[-1],
    Y=hexagoncentroidevenline.y[-1],Z=0),
    expand.grid(
      X=hexagoncentroidoddline.x,
      Y=hexagoncentroidoddline.y,Z=0))%>%as.matrix
  get_centroid<-function(x,y){
    (x-rangex[1])%/%distancebetweencentroidsonsamecolumn
    (y-rangey[1])%/%distancebetweencentroidsonsamecolumn
  }
  A=matrix(c(3/2,-sqrt(3)/4,3/2,sqrt(3)/4),2,2)
  projmat<-solve(t(A)%*%A)%*%t(A)
  coordinatesinbasis<-function(A){
    round(t(projmat%*%(t(A)-hexcentroids.mat[1,1:2])/hexagondiameterinangle))
  }
  centroidscoordinatesinbasis<-coordinatesinbasis(hexcentroids.mat[,1:2])
  
  
  hexagonvertices<-hexagondiameterinangle/2*rbind(cos(((0:6)/3)*pi),round(sin((((0:6)/3)*pi)),12),rep(0,7))
if(FALSE){  hexagons<-
    plyr::alply(hexcentroids,1,function(point){t(point+hexagonvertices)})%>%
    st_polygon()%>%
    st_sfc()  %>%
  st_set_crs(st_crs(.spatialdata))
  hexagons<-
    plyr::alply(hexcentroids,1,function(point){list(t(point+hexagonvertices))})%>%
    st_multipolygon()%>%
    st_sfc()  %>%
    st_set_crs(st_crs(.spatialdata))
  hexagons%>%st_sf()}

    hexagons<-
    do.call(rbind,
            plyr::llply(1:nrow(hexcentroids.mat),function(i){
              point=hexcentroids.mat[i,]
              st_sf(hexagon=list(t(point+hexagonvertices))%>%
                st_polygon()%>%
                st_sfc()  %>%
                st_set_crs(st_crs(.spatialdata)),
                HexId=i)}))
  hexcentroids<-
    do.call(rbind,
            plyr::llply(1:nrow(hexcentroids.mat),function(i){
              point=hexcentroids.mat[i,]
              st_sf(centroid=st_point(point)%>%
                      st_sfc()  %>%
                      st_set_crs(st_crs(.spatialdata)),
                    HexId=i)}))
  plot(hexagons,col="white")
  plot(hexcentroids,add=T,pch=19,col="black")
  mapview::mapview(hexagons)
  leaflet(envelope)%>%
    addTiles()%>%
    leafem::addFeatures()%>%
    leafem::addFeatures(data = hexagons)%>%
    leafem::addFeatures(data = hexcentroids,color='black',lwd=.3)
  
  
  hexcentroids<-st_multipoint(hexcentroids,dim="XY")%>% st_sfc()%>%st_set_crs(st_crs(.spatialdata))
  leaflet(hexcentroids)%>%addTiles()%>%addMarkers(hexcentroids)
  
  triangles<-sf::st_triangulate(st_combine(envelope))
  leaflet(.spatialdata)%>%addTiles()%>%leafem::addFeatures()
  hexcentroids<-
  
  plot(hexcentroids)
  plot(hexagons)
}