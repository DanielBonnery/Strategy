#' hexagonal bins
#' 
#' @param U : a dataframe containing the numerical variables x and y 
#' @param delta: controls the bin diameter
#' @return a hexbin object  hexagonal bins
#' @examples 
#' # plot the hex bins of cumbria
#' data(U,package="Strategy")
#' library("hexbin")
#' plot(neighbourhoods(U,.1))
#' plot(neighbourhoods(U,.01))
#' plot(neighbourhoods(U,.001))
neighbourhoods<-function(U,delta=(range(U$x,na.rm=TRUE)[2]-range(U$x,na.rm=TRUE)[1])/100){
  hexbin::hexbin(U$x, U$y, 
                 xbins =ceiling((range(U$x,na.rm=TRUE)[2]-range(U$x,na.rm=TRUE)[1])/(2*delta)), 
                 xbnds=range(U$x,na.rm=TRUE),
                 ybnds= range(U$y,na.rm=TRUE), 
                 IDs = TRUE)}

#' Distances between hexagonal bins 
#' 
#' @param U : a dataframe containing the numerical variables x and y and preferable hexagon
#' @param delta : needed if hexagon is not a variable of U: bins will be recomputed
#' @return a named matrix
#' @examples 
#' data(U) 
#' dist_areas_f(U)[1:3,1:3]
#' dist_areas_f(U,0.03)[1:3,1:3]

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



#' compute distances between new infected and exposed
#' 
#' @param closedistances NULL, or a named list with 2 named elements: closedistances$ra, closedistances$id 
#' @param U a data.frame with the variables hexagon (can be any bin identifier), x, y : coordinates, 
#' @param sicks a vector of integers
#' @param new.sicks a vector of integers
#' @param delta a positive number : a threshold
#' @param dist_areas: a function between 
#' @return NULL, or a named list with 2 named elements: closedistances$ra, closedistances$id 
#' @examples 
#' data(UE,package="Strategy")
#' delta<-.005
#' sicks<-(1:nrow(UE))[UE$I001=="sick"]
#' closedistances=newdist(NULL,UE,sicks)
#' do.call(cbind,closedistances)[1:3,]
newdist<-function(closedistances=NULL,U,sicks,new.sicks=NULL,delta=0.005,dist_areas=dist_areas_f(U,delta)){
  
  if(is.null(new.sicks)){new.sicks<-if(!is.null(closedistances)){setdiff(sicks,unique(closedistances$ind[,2]))}else{sicks}}
  
  stillfine<-setdiff(1:nrow(U),sicks)
  if(length(new.sicks)>0){
    new.sickhexagons<-unique(U$hexagon[new.sicks])
    
    L<- lapply(new.sickhexagons,function(hexagon){
      exposedhexagons<-dimnames(dist_areas)[[1]][dist_areas[hexagon,]<delta]
      new.sickinhexagon<-new.sicks[is.element(U$hexagon[new.sicks],hexagon)]
      stillfineexposedhexagon<-stillfine[is.element(U$hexagon[stillfine],exposedhexagons)]
      if(length(stillfineexposedhexagon)*length(new.sickinhexagon)>0){
        yy<-fields::fields.rdist.near(x1=U[stillfineexposedhexagon,c("x","y"),drop=FALSE],
                                      x2=U[new.sickinhexagon,c("x","y"),drop=FALSE], delta=delta,max.points=length(stillfineexposedhexagon)*length(new.sickinhexagon))
        return(if(!identical(yy$ra,-1)){
          if((length(stillfineexposedhexagon)*length(new.sickinhexagon)==1)|!is.matrix(yy$ind)){yy$ind<-matrix(yy$ind,1,2)}
          list(ind=cbind(stillfine[yy$ind[,1]],new.sicks[yy$ind[,2]]),ra=yy$ra)}else{list(ind=NULL,ra=NULL)})
      }})
    ra=do.call(c,lapply(L,function(x){x$ra}))
    ind<-do.call(rbind,lapply(L,function(x){x$ind}))
    if(!is.matrix(ind)&!is.null(ind)){ind<-matrix(ind,1,2)}
    if(!is.null(ind)){
      keep=is.element(ind[,1],stillfine)
      ind=ind[keep,];ra=ra[keep]
    }
    return(list(ind=ind,ra=ra))}
  else{return(list(ind=NULL,ra=NULL))}}
#' update the list of the already nown distances between subjects with distances between new infected and exposed
#' 
#' @param closedistances NULL, or a named list with 2 named elements: closedistances$ra, closedistances$id 
#' @param U a data.frame with the variables hexagon (can be any bin identifier), x, y : coordinates, 
#' @param sicks a vector of integers indicating the row numbers in U for sicks
#' @param new.sicks a vector of integers indicating the row numbers in U for new sicks
#' @param delta a positive number : a threshold
#' @param dist_areas: a function between 
#' @return NULL, or a named list with 2 named elements: closedistances$ra, closedistances$id 
#' @examples 
#' data(UE,package="Strategy")
#' delta<-.005
#' sicks<-(1:nrow(UE))[UE$I001=="sick"]
#' closedistances=updatedist(NULL,UE,sicks)
#' do.call(cbind,closedistances)[1:3,]
updatedist<-function(closedistances=NULL,U,sicks,new.sicks=NULL,delta=0.005,dist_areas=dist_areas_f(U,delta)){
  if(is.null(closedistances)){closedistances=list(ind=matrix(NA,0,2),ra=vector())}
  L<-newdist(closedistances,U,sicks,new.sicks,delta,dist_areas)
  list(ind=rbind(closedistances$ind,L$ind),ra=c(closedistances$ra,L$ra))}










#' computes the position of a projected point in the basis formed by a segment
#' @param p a numeric vector of length 2
#' @param s a 2x2 numeric matrix, representing a segment. Each row of the matrix are the coordinates of a extreme point of the segment.
#' @examples
#' zz<-function(p){
#' s=matrix(c(0,1,0,0),2,2)
#' plot(s,type="l",xlim=c(-.5,1.5),
#' xlab="",ylab="",
#' xaxt='n',yaxt='n')
#' points(x=p[1],y=p[2] ,col="red",cex=2)
#' points(closestpointonseg(p,s)[1],closestpointonseg(p,s)[2],col="red",cex=2)
#' segments(x0 = p[1],y0=p[2],x1=closestpointonseg(p,s)[1],y1=closestpointonseg(p,s)[2],col="red")
#' text(projpointonseg_a(p,s),-.8,paste0("a=",projpointonseg_a(p,s)))}
#' par(oma=c(0,0,0,0),mfrow=c(1,4))
#' zz(c(-.5,1))
#' zz(0:1)
#' zz(c(.5,1))
#' zz(c(1.5,1))
projpointonseg_a<-function(p,s,method="euclidean"){
  x1<-s[1,]
  X1<-p-x1
  X2<-s[2,]-x1
  den<-sum(X2^2)
  a0=crossprod(X1,X2)/(den+(den==0))
  min(1,max(0,a0))
}


#' computes the coordinates of the closest point on a segment to a point in the plane
#' @param p a numeric vector of length 2
#' @param s a 2x2 numeric matrix, representing a segment. Each row of the matrix are the coordinates of a extreme point of the segment.
#' @examples
#' zz<-function(p){
#' s=matrix(c(0,1,0,0),2,2)
#' plot(s,type="l",xlim=c(-.5,1.5),
#' xlab="",ylab="",
#' xaxt='n',yaxt='n')
#' points(x=p[1],y=p[2] ,col="red",cex=.5)
#' points(closestpointonseg(p,s)[1],closestpointonseg(p,s)[2],col="red",cex=.5)
#' segments(x0 = p[1],y0=p[2],x1=closestpointonseg(p,s)[1],y1=closestpointonseg(p,s)[2],col="red")}
#' par(mfrow=c(3,3))
#' set.seed(1);replicate(9,zz(c(runif(1,-.5,1.5),runif(1,-1,1))))

closestpointonseg<-function(p,s){
  s[1,]+projpointonseg_a(p,s)*(s[2,]-s[1,])
}


#' computes the coordinates of the closest point on the border of a polygon to a point in the plane
#' 
#' @param p a numeric vector of length 2
#' @param .poly a nx2 numeric matrix, representing a polygon. Each row of the matrix are the coordinates of a vertice of the polygon.
#' @return the coordinates of the closest point on the polygon
#' @examples  
#' zz<-function(){
#' .poly=matrix(sample(0:4,6,rep=T),3,2)[c(1:3,1),]
#' p<-sample(0:4,2,rep=T)
#' dd<-distpointtopoly(p,.poly)
#' plot(rbind(p,.poly),cex=.5,main=paste0("Distance: ", signif(dd,3)),
#' asp=1,xlim=range(cbind(p,.poly)),ylim=range(cbind(p,.poly)),xaxt='n',yaxt='n',xlab='',ylab='')
#' points(.poly,type="l",lwd=2)
#' cc<-rbind(p,closestpointonpolygon(p,.poly))
#' points(cc,col="red",cex=2)
#' points(cc,type="l",lty=3,col="red")
#' }
#' 
#' par(oma=c(0,0,0,0),mfrow=c(2,2))
#' set.seed(3);replicate(4,zz())
closestpointonpolygon<-function(p,.poly){
  dd<-plyr::aaply(1:(nrow(.poly)-1),1,function(i){distpointtoseg(p,.poly[i:(i+1),])})
  i=which(dd==min(dd))[1]
  closestpointonseg(p,.poly[i:(i+1),])}
#' computes the coordinates of the closest point on a segment to a point in the plane
#' @param s1 a 2x2 numeric matrix, representing a segment. Each row of the matrix are the coordinates of a extreme point of the segment.
#' @param s2 a 2x2 numeric matrix, representing a segment. Each row of the matrix are the coordinates of a extreme point of the segment.
#' @example 
#'zz<-function(){
#' s1=matrix(sample(0:4,4,rep=T),2,2)
#' s2=matrix(sample(0:4,4,rep=T),2,2)
#' s<-rbind(s1,s2)
#' dd<-distsegmenttosegment(s1,s2)
#' plot(s,cex=.5,main=paste0("Distance: ", signif(dd,3)),asp=1,xlim=range(s),ylim=range(s),xaxt='n',yaxt='n',xlab='',ylab='')
#' points(s1,type="l",lwd=2)
#' points(s2,type="l",lwd=2)
#' cc<-closestpointsontwosegments(s1,s2)
#' points(cc ,col="red",cex=2)
#' points(cc,type="l",col="red",lty=3)
#' }
#' 
#' par(oma=c(0,0,0,0),mfrow=c(3,3))
#' set.seed(3);replicate(9,zz())


closestpointsontwosegments<-function(s1,s2){
  dd<-distsegmenttosegment(s1,s2)
  if(dd>0){
    l<-which(c(distpointtoseg(s1[1,],s2),distpointtoseg(s1[2,],s2),distpointtoseg(s2[1,],s1),distpointtoseg(s2[2,],s1))==dd)[1]
    s<-if(l<=2){s2}else{s1}
    p1<-rbind(s1,s2)[l,]
    p2<-closestpointonseg(p1,s)}else{
      A=cbind(s1[2,]-s1[1,],s2[2,]-s2[1,])
      if(det(A)!=0){p1<-p2<-s1[1,]+(solve(A)%*%(s2[2,]-s1[1,]))[1]*A[,1]}else{
        pp<-rbind(s1[2,]-s1[1,],s2[1,]-s1[1,],s2[2,]-s1[1,])
        nn<-apply(pp,1,function(x){sum(x^2)})
        zz<-which(nn==max(nn))[1]
        l<-order(pp%*%pp[zz,])[2]
        p1<-p2<-s1[1,]+pp[l,]
      }}
  rbind(p1,p2)}

#'  returns a list of one or more  pair of points, one on each of two segments, with minimal distance
#' @param s1 a 2x2 numeric matrix, representing a segment. Each row of the matrix are the coordinates of a extreme point of the segment.
#' @param s2 a 2x2 numeric matrix, representing a segment. Each row of the matrix are the coordinates of a extreme point of the segment.
#' @return a list of  2x2 numeric matrices, representing a segment.
#'@examples
#'zz<-function(){
#' s1=matrix(sample(0:4,4,rep=T),2,2)
#' s2=matrix(sample(0:4,4,rep=T),2,2)
#' s<-rbind(s1,s2)
#' dd<-distsegmenttosegment(s1,s2)
#' plot(s,cex=.5,main=paste0("Distance: ", signif(dd,3)),asp=1,xlim=range(s),ylim=range(s),xaxt='n',yaxt='n',xlab='',ylab='')
#' points(s1,type="l",lwd=2)
#' points(s2,type="l",lwd=2)
#' cc<-closestpointsontwosegments_n(s1,s2)
#' for(ccc in cc){points(ccc ,col="red",cex=2)
#' points(ccc,type="l",col="red",lty=3)}
#' }
#' 
#' par(oma=c(0,0,0,0),mfrow=c(3,3))
#' set.seed(3);replicate(9,zz())
#'closestpointsontwosegments_n(matrix(c(0,3,0,0),2,2),matrix(c(1,2,0,0),2,2))
closestpointsontwosegments_n<-function(s1,s2){
  dd<-distsegmenttosegment(s1,s2)
  if(dd>0){
    ll<-which(c(distpointtoseg(s1[1,],s2),distpointtoseg(s1[2,],s2),distpointtoseg(s2[1,],s1),distpointtoseg(s2[2,],s1))==dd)
    L<-lapply(ll,function(l){
      s<-if(l<=2){s2}else{s1}
      p1<-rbind(s1,s2)[l,]
      p2<-closestpointonseg(p1,s)
      rbind(p1,p2)})}else{
        A=cbind(s1[2,]-s1[1,],s2[2,]-s2[1,])
        if(det(A)!=0){p1<-p2<-s1[1,]+(solve(A)%*%(s2[2,]-s1[1,]))[1]*A[,1]
        L<-list(rbind(p1,p2))}else{
          pp<-rbind(c(0,0),s1[2,]-s1[1,],s2[1,]-s1[1,],s2[2,]-s1[1,])
          nn<-apply(pp,1,function(x){sum(x^2)})
          zz<-which(nn==max(nn))[1]
          ll<-order(pp%*%pp[zz,])[2:3]
          L<-lapply(ll,function(l){
            p1<-p2<-s1[1,]+pp[l,]
            rbind(p1,p2)})
        }}
  L}


#' computes the coordinates of the closest point on the border of a polygon to a point in the plane
#' 
#' @param poly1 a nx2 numeric matrix, representing a polygon. Each row of the matrix are the coordinates of a vertice of the polygon.
#' @param poly2 a nx2 numeric matrix, representing a polygon. Each row of the matrix are the coordinates of a vertice of the polygon.
#' @return the coordinates of the closest point on the polygon
#' @examples 
#' zz<-function(){
#' poly1=matrix(sample(0:6,6,rep=T),3,2)[c(1:3,1),]
#' poly2=matrix(sample(0:6,6,rep=T),3,2)[c(1:3,1),]
#' s<-rbind(poly1,poly2)
#' dd<-distpolytopoly(poly1,poly2)
#' plot(s,cex=.5,main=paste0("Distance: ", signif(dd,3)),asp=1,xlim=range(s),ylim=range(s),xaxt='n',yaxt='n',xlab='',ylab='')
#' points(poly1,type="l",lwd=2)
#' points(poly2,type="l",lwd=2)
#' cc<-closestpointsontwopolygons(poly1,poly2)
#' points(cc ,col="red",cex=2)
#' points(cc,type="l",col="red",lty=3)
#' }
#' 
#' par(oma=c(0,0,0,0),mfrow=c(2,2))
#' set.seed(2);replicate(4,zz())


closestpointsontwopolygons<-function(poly1,poly2){
  x=Inf;
  for (i in 1:(nrow(poly1)-1)){
    dd<-plyr::aaply(1:(nrow(poly2)-1),1,function(j){distsegmenttosegment(poly1[i:(i+1),],poly2[j:(j+1),])})
    if(min(dd)<x){
      x<-min(dd)
      j=which(dd==x)[1]
      s1=poly1[i:(i+1),];s2=poly2[j:(j+1),]
    }
  }
  closestpointsontwosegments(s1,s2)
}

#' Compute minimum distance between two polygons
#' @param poly1 a nx2 numeric matrix, representing a polygon. Each row of the matrix are the coordinates of a vertice of the polygon.
#' @param poly2 a nx2 numeric matrix, representing a polygon. Each row of the matrix are the coordinates of a vertice of the polygon.
#' @examples 
#' zz<-function(){
#' poly1=matrix(sample(0:6,6,rep=T),3,2)[c(1:3,1),]
#' poly2=matrix(sample(0:6,6,rep=T),3,2)[c(1:3,1),]
#' s<-rbind(poly1,poly2)
#' dd<-distpolytopoly(poly1,poly2)
#' plot(s,cex=.5,main=paste0("Distance: ", signif(dd,3)),asp=1,xlim=range(s),ylim=range(s),xaxt='n',yaxt='n',xlab='',ylab='')
#' points(poly1,type="l",lwd=2)
#' points(poly2,type="l",lwd=2)
#' for(cc in closestpointsontwopolygons_n(poly1,poly2)){
#' points(cc ,col="red",cex=2)
#' points(cc,type="l",col="red",lty=3)}
#' }
#' 
#' par(oma=c(0,0,0,0),mfrow=c(2,2))
#' set.seed(2);replicate(4,zz())

closestpointsontwopolygons_n<-function(poly1,poly2){
  x=Inf;
  L<-list()
  for (i in 1:(nrow(poly1)-1)){
    dd<-plyr::aaply(1:(nrow(poly2)-1),1,function(j){distsegmenttosegment(poly1[i:(i+1),],poly2[j:(j+1),])})
    if(min(dd)==x){
      jj=which(dd==x)
      L<-c(L,lapply(jj,function(j){
        s1=poly1[i:(i+1),];
        s2=poly2[j:(j+1),]
        rbind(s1,s2)}))
    }
    if(min(dd)<x){
      x<-min(dd)
      jj=which(dd==x)
      L<-lapply(jj,function(j){
        s1=poly1[i:(i+1),];
        s2=poly2[j:(j+1),]
        rbind(s1,s2)})
    }
  }
  do.call(c,lapply(L,function(l){
    closestpointsontwosegments_n(l[1:2,],l[3:4,])}))
}

#' computes the distance between a point and a segment
#' 
#' @param p a numeric vector of length 2
#' @param s a 2x2 numeric matrix, representing a segment. Each row of the matrix are the coordinates of a extreme point of the segment.
#' @examples
#' zz<-function(p){
#' s=matrix(c(0,3,0,0),2,2)
#' plot(s,type="l",xlim=c(-2,5),ylim=c(-2,2),
#' xlab="",ylab="",
#' xaxt='n',yaxt='n')
#' points(x=p[1],y=p[2] ,col="red",cex=.5)
#' points(closestpointonseg(p,s)[1],closestpointonseg(p,s)[2],col="red",cex=.5)
#' segments(x0 = p[1],y0=p[2],x1=closestpointonseg(p,s)[1],y1=closestpointonseg(p,s)[2],col="red")
#' text((closestpointonseg(p,s)[1]+p[1])/2,(closestpointonseg(p,s)[2]+p[2])/2,round(distpointtoseg(p,s),2))}
#' par(mfrow=c(3,3))
#' set.seed(1);replicate(9,zz(c(sample(-2:3,1),sample(-2:2,1))))

distpointtoseg<-function(p,s){
  X1<-p-s[1,]
  X2<-s[2,]-s[1,]
  sqrt(sum((X1-(projpointonseg_a(p,s)*X2))^2))
  #dist(rbind(X1,projpointonseg_a(p,s)*X2))
}

#' computes the orientation of a triangle
#' 
#' @param s a 3x2 numeric matrix, representing a triangle. Each row of the matrix are the coordinates of an extreme point of the triangle.
#' @examples
#' zz<-function(p){
#' s=matrix(sample(0:4,6,rep=T),3,2)
#' plot(s[c(1:3,1),],type="l",main=if(triangleorientation(s)==1){"+"}else{if(triangleorientation(s)==-1){"-"}else{"aligned"}},xlim=c(-1,5),ylim=c(-1,5),xlab="",ylab="",xaxt='n', yaxt='n')
#' text(s[,1],s[,2],toupper(letters[1:3]),cex=1,col="red")
#' }
#' par(mfrow=c(2,2),mar=c(5,3,2,2))
#' set.seed(7);replicate(4,zz(c(sample(-2:3,1),sample(-2:2,1))))
triangleorientation<-function(s){sign(det(s[2:3,]-s[c(1,1),]))}


#' test if two segments intersect
#' 
#' @param s1 a 2x2 numeric matrix, representing a segment. Each row of the matrix are the coordinates of an extreme point of the segment.
#' @param s2 a 2x2 numeric matrix, representing a segment. Each row of the matrix are the coordinates of an extreme point of the segment.
#' @examples
#' zz<-function(s1=matrix(sample(0:3,4,rep=T),2,2),s2=matrix(sample(0:3,4,rep=T),2,2)){
#' si<-segment.intersect(s1,s2)
#' s<-rbind(s1,s2)
#' plot(s,cex=.5,main=if(si){"Intersect"}else{"Disjoint"},xlab="",ylab="",xaxt='n', yaxt='n',xlim=range(s)+c(-1,1),ylim=range(s)+c(-1,1))
#' points(s1,type="l")
#' points(s2,type="l")
#' text(s[,1],s[,2],toupper(letters[1:4]),cex=1,col="red")
#' 
#' }
#' par(mfrow=c(3,4),mar=c(5,3,2,2))
#' set.seed(12);replicate(4,zz())
#' zz(matrix(0,2,2),matrix(0,2,2))
#' zz(matrix(0,2,2),matrix(1,2,2))
#' zz(matrix(c(1,1,1,1),2,2),matrix(c(2,3,2,3),2,2))
#' zz(matrix(c(1,1,0,0),2,2),matrix(c(0,1,0,0),2,2))
#' zz(matrix(c(1,3,1,3),2,2),matrix(c(2,4,2,4),2,2))
#' zz(matrix(c(0,1,0,1),2,2),matrix(c(2,3,2,3),2,2))
#' zz(matrix(c(0,4,0,4),2,2),matrix(c(1,3,1,3),2,2))
#' zz(matrix(c(0,1,0,1),2,2),matrix(c(0,3,0,3),2,2))
segment.intersect<-function(s1,s2){
  x=c(triangleorientation(rbind(s1[1,],s2)),
      triangleorientation(rbind(s1[2,],s2)),
      triangleorientation(rbind(s2[1,],s1)),
      triangleorientation(rbind(s2[2,],s1)))
  y<-max(prod(x[1:2]),prod(x[3:4]))
  if(y==1){FALSE}else{if(y==-1){TRUE}else{
    w<-s1[2,]-s1[1,]
    if(any(w!=0)){
      z=c(sum(w^2),sort(c(crossprod(s2[1,]-s1[1,],w),crossprod(s2[2,]-s1[1,],w))))
      (z[3]>z[1]&z[2]<z[1])|(z[3]<z[1]&z[3]>0)}else{
        w<-s2[2,]-s2[1,]
        if(any(w!=0)){
          z=c(sum(w^2),sort(c(crossprod(s1[1,]-s2[1,],w),crossprod(s1[2,]-s2[1,],w))))
          (z[3]>=z[1]&z[2]<=z[1])|(z[3]<=z[1]&z[3]>=0)
        }else{all(s1[1,]==s2[1,])}
      }
  }}}


#' tells if two ranges overlap
#' @param r1 a range a length 2 numerical vector
#' @param r2 a range a length 2 numerical vector
#' @retrun TRUE if two ranges overlap
#' @examples
#' par(mfrow=c(1,4),oma=c(0,0,0,0))
#' set.seed(8);replicate(4,(function(){
#' r1<-sample(1:8,2);r2<-sample(1:5,2)
#' plot(cbind(c(r1,r2),0),yaxt='n',xlab='',ylab='',xaxt='n',
#' main=paste0(if(rangesoverlap(r1,r2)){"O"}else{"Does not o"},"verlap"))
#' points(cbind(r1,0),type="l",lwd=3)
#' points(cbind(r2,0),type="l",col="red")
#' })())
rangesoverlap<-function(r1,r2){
  !((min(r2)>max(r1))|(min(r1)>max(r2)))
}

#' Interval length between two ranges (0 if they overlap)
#' @param r1 a range a length 2 numerical vector
#' @param r2 a range a length 2 numerical vector
#' @retrun TRUE if two ranges overlap
#' 
#' @examples
#' par(mfrow=c(1,4),oma=c(0,0,0,0))
#' set.seed(10);replicate(4,(function(){
#' r1<-sample(1:8,2);r2<-sample(1:5,2)
#' plot(cbind(c(r1,r2),0),yaxt='n',xlab='',ylab='',
#' main=paste0("Gap: ",ranges.gap(r1,r2)))
#' points(cbind(r1,0),type="l",lwd=3)
#' points(cbind(r2,0),type="l",col="red")
#' })())
ranges.gap<-function(r1,r2){
  max(0,min(r2)-max(r1),min(r1)-max(r2))
}



#' distance segment to segment

#' @param s1 a 2x2 numeric matrix, representing a segment. Each row of the matrix are the coordinates of a extreme point of the segment.
#' @param s2 a 2x2 numeric matrix, representing a segment. Each row of the matrix are the coordinates of a extreme point of the segment.
#' @return a number
#' @examples
#' zz<-function(){
#' s1=matrix(sample(0:4,4,rep=T),2,2)
#' s2=matrix(sample(0:4,4,rep=T),2,2)
#' s<-rbind(s1,s2)
#' dd<-distsegmenttosegment(s1,s2)
#' plot(s,cex=.5,main=paste0("Distance: ", signif(dd,3)),asp=1,xlim=range(s),ylim=range(s),xaxt='n',yaxt='n',xlab='',ylab='')
#' points(s1,type="l",lwd=2)
#' points(s2,type="l",lwd=2)
#' cc<-closestpointsontwosegments(s1,s2)
#' points(cc ,col="red",cex=2)
#' points(cc,type="l",col="red",lty=3)
#' }
#' 
#' par(mfrow=c(3,3),oma=c(0,0,0,0),mar=c(1.1,1.2,1,1.1))
#' set.seed(3);replicate(9,zz())

distsegmenttosegment<-function(s1,s2){
  if(segment.intersect(s1,s2)){0}else{
    min(c(plyr::aaply(s1,1,distpointtoseg,s=s2),plyr::aaply(s2,1,distpointtoseg,s=s1)),na.rm=T)}
}

#' computes the distance between a point and a polygon
#' 
#' @param p a numeric vector of length 2
#' @param .poly a n x2 numeric matrix, representing a polygon. Each row of the matrix are the coordinates of a verticeof the polygon.
#' @examples
#' zz<-function(){
#' p<-sample(0:6,2,rep=T)
#' .poly<-matrix(sample(0:6,6,rep=T),3,2)[c(1:3,1),]
#' plot(rbind(.poly,p),
#' xlab="",ylab="",
#' cex=.2,main=paste0("Distance: ",signif(distpointtopoly(p,.poly),3)))
#' points(.poly,type='l')
#' points(x=p[1],y=p[2] ,col="red",cex=1)
#' points(closestpointonpolygon(p,.poly)[1],closestpointonpolygon(p,.poly)[2],col="red",cex=1)
#' points(rbind(p,closestpointonpolygon(p,.poly)),col="red",lty=3,type='l')}
#' par(mfrow=c(3,3),oma=c(0,0,1,0),mar=c(2,2.1,1,0.1))
#' set.seed(1);replicate(9,zz())

distpointtopoly<-function(p,.poly){
  min(plyr::aaply(1:(nrow(.poly)-1),1,function(i){distpointtoseg(p,.poly[i:(i+1),])}))}


#' Distance between a segment and a polygon
#' 
#' @param .poly a polygon (a nx2 matrix, each line is a point)
#' @param s a 2x2 numeric matrix, representing a segment. Each row of the matrix are the coordinates of a extreme point of the segment.
#' @examples
#' zz<-function(){
#' .poly<-matrix(sample(0:6,6,T),3,2)[c(1:3,1),]
#' s<-matrix(sample(0:6,6,T),2,2)
#' plot(rbind(.poly,s),xlab="",yaxt="n")
#' points(.poly,type="l");points(s,type="l")
#' points(closestpointsontwopolygons(s,.poly),lty=3,col="red")}
#' 


distsegmenttopoly<-function(s,.poly){
  min(plyr::aaply(1:(nrow(.poly)-1),1,function(i){distsegmenttosegment(s,.poly[i:(i+1),])}))}


#' Computes distance between two polygons
#' @param poly1 a polygon (a n x 2 numerical matrix)
#' @param poly2 a polygon (a n x 2 numerical matrix)
#' @return a positive number, the distance between the two polygons
#' @seealso distsegmenttopoly
#' @examples
#' zz<-function(){
#' poly1=matrix(sample(0:6,6,rep=T),3,2)[c(1:3,1),]
#' poly2=matrix(sample(0:6,6,rep=T),3,2)[c(1:3,1),]
#' s<-rbind(poly1,poly2)
#' dd<-distpolytopoly(poly1,poly2)
#' plot(s,cex=.5,main=paste0("Distance: ", signif(dd,3)),asp=1,xlim=range(s),ylim=range(s),xaxt='n',yaxt='n',xlab='',ylab='')
#' points(poly1,type="l",lwd=2)
#' points(poly2,type="l",lwd=2)
#' for(cc in closestpointsontwopolygons_n(poly1,poly2)){
#' points(cc,type="l",col="red",lty=3)}}
#' 
#' par(mfrow=c(2,2),oma=c(0,0,1,0),mar=c(0.1,0.1,1,0.1))
#' set.seed(2);replicate(4,zz())
distpolytopoly<-function(poly1,poly2){
  min(plyr::aaply(1:(nrow(poly1)-1),1,function(i){distsegmenttopoly(poly1[i:(i+1),],poly2)}))
}

#' Computes distance between two polygons, only works for polygons that do not intersect
#' @param poly1 a polygon (a n x 2 numerical matrix)
#' @param poly2 a polygon (a n x 2 numerical matrix)
#' @return a positive number, the distance between the two polygons
distpolytopoly2<-function(poly1,poly2){
  min(plyr::aaply(poly1,1,distpointtopoly,.poly=poly2),
      plyr::aaply(poly2,1,distpointtopoly,.poly=poly1))}

#' converts a shapefile to list of  polygons (nx2 matrices)
extractpolygonsaslist<-function(shp){
  lapply(1:nrow(shp),function(i){shp[i,]@polygons[[1]]@Polygons[[1]]@coords})}

#' Compute distance matrix for a list of polygons
#' @param list.poly a list of nx2 numeric matrices
#' @return a  (n*(n-1)/2)x 3 matrix  
#' @examples
#' zz<-function(){
#' list.poly=plyr::alply(cbind(rep(0:8,9),rep(0:8,each=9))[sample(81,4),],1,function(x){
#' cbind(x[1]+c(0,0,.5,.5,0),x[2]+c(0,.5,.5,0,0))})
#' gradients=cbind(c(0,1),c(1,0),c(1,1),c(1,-1))
#' par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(0.1,0.1,.1,0.1))
#' plot(do.call(rbind,list.poly),xlab="",yaxt="n",ylab="",cex=.1)
#' for(i in 1:length(list.poly)){.poly=list.poly[[i]]
#' points(.poly,type="l")
#' text(mean(.poly[,1]),mean(.poly[,2]),as.roman(i))
#' }
#' X=polydistmat(list.poly)
#' X<-cbind(X,floor(rank(X[,3])))
#' colorlink=topo.colors(2*max(X[,4]))[X[,4]]
#' for(i in 1:nrow(X)){
#' cc<-closestpointsontwopolygons(list.poly[[X[i,1]]],list.poly[[X[i,2]]])
#' points(cc,col=colorlink[i],type="l",lty=3)
#' text(mean(cc[,1]),mean(cc[,2]),signif(X[i,3],3))
#' }
#' colnames(X)<-c("polygon 1","polygon 2", "distance","col")
#' X[,1:3]}
#' set.seed(1);zz()
polydistmat<-function(list.poly){
  L<-plyr::alply(1:(length(list.poly)-1),1,function(i){
    do.call(rbind,
            parallel::mclapply((i+1):length(list.poly),function(j){
              c(i,j,distpolytopoly(list.poly[[i]],list.poly[[j]]))}))})
  do.call(rbind,L)}



#' Compute range of a polygon along certain gradients
#' 
#' @param list.poly a list of nx2 numeric matrices
#' @param delta  a positive number
#' @param gradients a 2x n matrix each column representing a vector. 
#' @examples
#' a=c(1-1/sqrt(2),1/sqrt(2));.poly=cbind(1+c(0,0,a,1,1,a[2:1]),1+c(a,1,1,a[2:1],0,0))[c(1:8,1),]
#' plot(.poly,type='l',xlab='',ylab='')
#' rangesbygradients_f(.poly,gradients=cbind(c(0,1),c(1,0),c(1,-1)))


rangesbygradients_f<-function(.poly,gradients=-apply(cbind(c(0,1),c(1,0),c(1,1),c(1,-1)),2,function(x){x/(sqrt(sum(x^2)))})){
apply(.poly%*%(gradients),2,range)
}

#' Compute distance matrix for a list of polygons
#' @param list.poly a list of nx2 numeric matrices
#' @param delta  a positive number
#' @param gradients a 2x n matrix each column representing a vector. 
#' @return a  (n*(n-1)/2)x 3 matrix  
#' @examples
#' zz<-function(delta){
#' list.poly=plyr::alply(cbind(rep(0:8,9),rep(0:8,each=9))[sample(81,20),],1,function(x){
#' cbind(x[1]+c(0,0,.5,.5,0),x[2]+c(0,.5,.5,0,0))})
#' gradients=cbind(c(0,1),c(1,0),c(1,1),c(1,-1))
#' par(mfrow=c(1,1),oma=c(0,0,1,0),mar=c(0.1,0.1,1,0.1))
#' plot(do.call(rbind,list.poly),xlab="",yaxt="n",ylab="",cex=.1,main=paste0("Match polygons distant less than ",delta))
#' for(.poly in list.poly){points(.poly,type="l")}
#' X=polysmalldistmat(list.poly,delta)
#' for(i in 1:nrow(X)){
#' points(closestpointsontwopolygons(list.poly[[X[i,1]]],list.poly[[X[i,2]]]),col="red",type="l",lty=3)
#' }}
#' set.seed(1);zz(.5)
#' set.seed(1);zz(1)
#' set.seed(1);zz(2)
polysmalldistmat<-function(list.poly,delta,gradients=apply(cbind(c(0,1),c(1,0),c(1,1),c(1,-1)),2,function(x){x/(sqrt(sum(x^2)))})){
  if(delta==Inf){polydistmat(list.poly)}else{
  n<-ncol(gradients)
  gradients=apply(gradients,2,function(x){x/(sqrt(sum(x^2)))})
  rangesbygradients<-lapply(list.poly,rangesbygradients_f,gradients=gradients)
  print("ranges by gradient computed")
  print(paste0(length(list.poly)," polygons."))
  L<-parallel::mclapply(1:(length(list.poly)-1),function(i,.rangesbygradients,.delta,.list.poly,.n){
    smalls<-sapply((i+1):length(.list.poly),function(j){
      all(sapply(1:.n,function(k){
        ranges.gap(.rangesbygradients[[i]][,k],.rangesbygradients[[j]][,k])<=.delta}))})
    if(any(smalls)){
      do.call(rbind,
              lapply(((i+1):length(.list.poly))[smalls],function(j){
                c(i,j,Strategy::distpolytopoly(.list.poly[[i]],.list.poly[[j]]))}))}else{matrix(NA,0,3)}},
    .rangesbygradients=rangesbygradients,
    .delta=delta,
    .list.poly=list.poly,
    .n=n,
    mc.cores=parallel::detectCores())
  D<-do.call(rbind,L)
  D[D[,3]<=delta,]}}


#'@example
#'data(Avo_fields)
#'list.poly<-extractpolygonsaslist(Avo_fields)
#'dd<-polysmalldistmat(list.poly,.1)
#'save(dd,file="~/daniel.bonnery@gmail.com/dd.ra")





#'Connected populations
#' @param MM a distance matrix (column one is an integer id, column 2 is an integer id, column 3 is a distance) in a certain unit u
#' @param delta a distance in the same unit u
#' @param n the number of initial bins
#'@examples
#' zz<-function(delta){
#' 
#' list.poly=plyr::alply(cbind(rep(0:8,9),rep(0:8,each=9))[sample(81,9),],1,function(x){
#' cbind(x[1]+c(0,0,.5,.5,0),x[2]+c(0,.5,.5,0,0))})
#' MM=polysmalldistmat(list.poly,delta+1)
#' bins<-connectedpop(MM,delta,n=9)
#' 
#' par(mfrow=c(1,1),oma=c(0,0,0,0),mar=c(0.1,0.1,1.1,0.1))
#' plot(do.call(rbind,list.poly),xlab="",yaxt="n",ylab="",cex=.1,main=paste0("regroup when dist<",delta))
#' for(i in 1:length(list.poly)){.poly=list.poly[[i]]
#' points(.poly,type="l",col=as.factor(bins$bin)[i])
#' text(mean(.poly[,1]),mean(.poly[,2]),as.roman(i))
#' }
#' colorlink=sapply(MM[,3],function(x){if(x<=delta){"green"}else{"red"}})
#' for(i in 1:nrow(MM)){
#' cc<-closestpointsontwopolygons(list.poly[[MM[i,1]]],list.poly[[MM[i,2]]])
#' points(cc,col=colorlink[i],type="l",lty=3)
#' text(mean(cc[,1]),mean(cc[,2]),signif(MM[i,3],3))
#' 
#' }
#' colnames(MM)<-c("polygon 1","polygon 2", "distance")
#' MM[,1:3]
#' }
#' set.seed(4);delta=2;zz(delta)

connectedpop<-function(MM,delta,n=max(MM[,1:2])){
  x=MM[MM[,3]<delta,1:2]
  bins<-data.frame(polygon=1:n,
                   bin=1:n)
  someremain=TRUE
  nextbins<-sort(unique(x[,1]))
  while(someremain){
    nextbin=nextbins[1]
    islinked<-(x[,1]==nextbin)
    while(any(islinked)){
      newtobin<-x[islinked,2]
      bins$bin[newtobin]<-nextbin
      if(any(!islinked)){
        x<-x[!islinked,,drop=FALSE]
        x[is.element(x,newtobin)]<-nextbin
        x<-plyr::aaply(x,1,sort,.drop = FALSE)
        islinked<-x[,1]==nextbin
        nextbins<-setdiff(nextbins,c(nextbin,newtobin))}else{islinked=FALSE;nextbins=c()}
    }
    someremain=length(nextbins)>0
  }
  bins
}



