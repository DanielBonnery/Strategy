#distances on the sphere


degtorad<-function(x){pi*x/180}

latlongtoe<-function(latlong){
  LatLong<-degtorad(latlong)
  ca<-cos(LatLong[1])
  sa<-sin(LatLong[1])
  co<-cos(LatLong[2])
  so<-sin(LatLong[2])
  c(co*ca,so*ca,sa)
}

Raddist<-function(s){
  normalvectors<-apply(s,1,laglongtoe)
  acos(crossprod(normalvectors[1,],normalvectors[2,]))}


Greatcircledist<-function(s,R=6371){
  R*Raddist(s)}
