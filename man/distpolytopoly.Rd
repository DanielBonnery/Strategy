% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Distances.R
\name{distpolytopoly}
\alias{distpolytopoly}
\title{Computes distance between two polygons}
\usage{
distpolytopoly(poly1, poly2)
}
\arguments{
\item{poly1}{a polygon (a n x 2 numerical matrix)}

\item{poly2}{a polygon (a n x 2 numerical matrix)}
}
\value{
a positive number, the distance between the two polygons
}
\description{
Computes distance between two polygons
}
\examples{
zz<-function(){
poly1=matrix(sample(0:6,6,rep=T),3,2)[c(1:3,1),]
poly2=matrix(sample(0:6,6,rep=T),3,2)[c(1:3,1),]
s<-rbind(poly1,poly2)
dd<-distpolytopoly(poly1,poly2)
plot(s,cex=.5,main=paste0("Distance: ", signif(dd,3)),asp=1,xlim=range(s),ylim=range(s),xaxt='n',yaxt='n',xlab='',ylab='')
points(poly1,type="l",lwd=2)
points(poly2,type="l",lwd=2)
for(cc in closestpointsontwopolygons_n(poly1,poly2)){
points(cc,type="l",col="red",lty=3)}}

par(mfrow=c(2,2),oma=c(0,0,1,0),mar=c(0.1,0.1,1,0.1))
set.seed(2);replicate(4,zz())
}
\seealso{
distsegmenttopoly
}
