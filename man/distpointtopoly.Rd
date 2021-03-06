% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Distances.R
\name{distpointtopoly}
\alias{distpointtopoly}
\title{computes the distance between a point and a polygon}
\usage{
distpointtopoly(p, .poly)
}
\arguments{
\item{p}{a numeric vector of length 2}

\item{.poly}{a n x2 numeric matrix, representing a polygon. Each row of the matrix are the coordinates of a verticeof the polygon.}
}
\description{
computes the distance between a point and a polygon
}
\examples{
zz<-function(){
p<-sample(0:6,2,rep=T)
.poly<-matrix(sample(0:6,6,rep=T),3,2)[c(1:3,1),]
plot(rbind(.poly,p),
xlab="",ylab="",
cex=.2,main=paste0("Distance: ",signif(distpointtopoly(p,.poly),3)))
points(.poly,type='l')
points(x=p[1],y=p[2] ,col="red",cex=1)
points(closestpointonpolygon(p,.poly)[1],closestpointonpolygon(p,.poly)[2],col="red",cex=1)
points(rbind(p,closestpointonpolygon(p,.poly)),col="red",lty=3,type='l')}
par(mfrow=c(3,3),oma=c(0,0,1,0),mar=c(2,2.1,1,0.1))
set.seed(1);replicate(9,zz())
}
