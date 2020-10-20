#' All ordered n-sized subsets of 1...N 
#' n 
#'@examples 
#'subsets(1,10)
#'subsets(2,10)
#'subsets(3,10)
#'subsets(8,10)
#'subsets(9,10)
#'subsets(10,10)

subsets<-function(n,N,subsets_1=if(n>=2){subsets(n-1,N)}else{NULL}){
  if(n<=0|n>N){M<-matrix()}
  if(n==1){M<-matrix(1:N,N,1)}
  if(n>=2&n<=N/2){M<-do.call(rbind,plyr::alply(subsets_1[subsets_1[,n-1]!=N,],1,function(x){do.call(cbind,c(as.list(x),list((max(x)+1):N)))}))}
  if(n>N/2&n<N){M<-plyr::aaply(subsets(N-n,N),1,function(x){setdiff(1:N,x)})}
  if(n==N){M<-matrix(1:N,1,N)}
  M}


figures1..3<-function(){
  
  library(ggplot2)
  N=15
  nrep=6
  set.seed(1)
  LL<-plyr::alply(1:nrep,1,function(Rep){
    #N<-rpois(1,N)
    X=matrix(runif(N*2),N,2)
    colnames(X)<-c("x","y")
    XX<-exp(-4*as.matrix(dist(X,upper=T,diag = T)))
    objective<-apply(XX,1,mean)
    optimalsamples<-lapply(1:5,function(n){
      subsets_n<-subsets(n,N)
      prob<-plyr::aaply(subsets_n,1,function(x){
        mean(1-plyr::aaply(1-XX[,x],1,prod))
      })
      argmax=subsets_n[which(prob==max(prob)),]
    })
    argmaxes<-do.call(cbind,lapply(optimalsamples,function(x){is.element(1:N,x)}))
    argmax=(objective==max(objective))
    
    infected=plyr::raply(10000,
                         .expr = (function(){infected0=  infected0=is.element(1:N,sample(N,1));
                         cbind(infected0,infected1=rbinom(1:N,1,prob=XX[infected0,]))})())
    list(.data=data.frame(id=1:N,X,
                          objective=objective,
                          argmaxes=argmaxes,
                          rep=paste0("Population ",Rep),
                          meaninfected0=plyr::aaply(infected[,,1],2,mean),
                          meaninfected1=plyr::aaply(infected[,,2],2,mean)))
  })
  .data=do.call(rbind,lapply(LL, function(L){L$.data}))
  figure1<-ggplot(data=.data,aes(x=x,y=y,color=objective,fill=objective))+
    geom_point(colour="black",shape=21,size=4)+ 
    scale_fill_gradient(low = "#FFFFFF",high = "#000000")+
    geom_point(data=.data[.data$argmaxes.1==1,],mapping=aes(x=x,y=y),shape=10,size=10,color="black")+
    facet_wrap(~rep)+ theme(legend.position="bottom")
  
  print(figure1)
  
  figure2<-ggplot(data=.data,aes(x=x,y=y,color=objective,fill=objective))+
    geom_point(colour="black",shape=21,size=4)+ 
    scale_fill_gradient(low = "#FFFFFF",high = "#000000")+
    geom_point(data=.data[.data$argmaxes.2==1,],mapping=aes(x=x,y=y),shape=10,size=10,color="black")+
    facet_wrap(~rep)+ theme(legend.position="bottom")
  
  
  .data2<-reshape2::melt(.data[c("id","rep","x","y","objective",paste0("argmaxes.",1:5))],id=c("id","rep","x","y","objective"))
  .data2$variable<-gsub(pattern = "argmaxes.",replacement="n=",.data2$variable)
  
  figure3<-ggplot(data=.data2,aes(x=x,y=y,color=objective,fill=objective))+
    geom_point(colour="black",shape=21,size=2)+ 
    scale_fill_gradient(low = "#FFFFFF",high = "#000000")+
    geom_point(data=.data2[.data2$value,],mapping=aes(x=x,y=y),shape=10,size=4,color="black")+
    facet_grid(vars(rep),vars(variable))+ theme(legend.position="none",axis.title.x=element_blank(),
                                                axis.text.x=element_blank(),
                                                axis.ticks.x=element_blank(),,axis.title.y=element_blank(),
                                                axis.text.y=element_blank(),
                                                axis.ticks.y=element_blank())
  return(list(figure1=figure1,figure2=figure2,figure3=figure3))
}




figures4..6<-function(){
  
  library(ggplot2)
  Nh=4
  H=4
  N=Nh*H
  nrep=2
  set.seed(1)
  LL<-plyr::alply(1:nrep,1,function(Rep){
    #N<-rpois(1,N)
    X=matrix(runif(H*2),H,2)
    X=do.call(rbind,plyr::alply(X,1,function(x){matrix(x[rep(1:2,each=Nh)]+.1*runif(Nh*2),Nh,2)}))
    colnames(X)<-c("x","y")
    XX<-exp(-4*as.matrix(dist(X,upper=T,diag = T)))
    objective<-apply(XX,1,mean)
    optimalsamples<-lapply(1:5,function(n){
      subsets_n<-subsets(n,N)
      prob<-plyr::aaply(subsets_n,1,function(x){
        mean(1-plyr::aaply(1-XX[,x],1,prod))
      })
      argmax=subsets_n[which(prob==max(prob)),]
    })
    argmaxes<-do.call(cbind,lapply(optimalsamples,function(x){is.element(1:N,x)}))
    argmax=(objective==max(objective))
    
    #infected=plyr::raply(10000,
    #                     .expr = (function(){infected0=  infected0=is.element(1:N,sample(N,1));
    #                     cbind(infected0,infected1=rbinom(1:N,1,prob=XX[infected0,]))})())
    list(.data=data.frame(id=1:N,X,
                          objective=objective,
                          argmaxes=argmaxes,
                          rep=paste0("Population ",Rep)#,
                          #meaninfected0=plyr::aaply(infected[,,1],2,mean),
                          #meaninfected1=plyr::aaply(infected[,,2],2,mean)
    ))
  })
  .data=do.call(rbind,lapply(LL, function(L){L$.data}))
  figure4<-ggplot(data=.data,aes(x=x,y=y,color=objective,fill=objective))+
    geom_point(colour="black",shape=21,size=4)+ 
    scale_fill_gradient(low = "#FFFFFF",high = "#000000")+
    geom_point(data=.data[.data$argmaxes.1==1,],mapping=aes(x=x,y=y),shape=10,size=10,color="black")+
    facet_wrap(~rep)+ theme(legend.position="bottom")
  
  print(figure4)
  
  figure5<-ggplot(data=.data,aes(x=x,y=y,color=objective,fill=objective))+
    geom_point(colour="black",shape=21,size=4)+ 
    scale_fill_gradient(low = "#FFFFFF",high = "#000000")+
    geom_point(data=.data[.data$argmaxes.2==1,],mapping=aes(x=x,y=y),shape=10,size=10,color="black")+
    facet_wrap(~rep)+ theme(legend.position="bottom")
  
  
  .data2<-reshape2::melt(.data[c("id","rep","x","y","objective",paste0("argmaxes.",1:5))],id=c("id","rep","x","y","objective"))
  .data2$variable<-gsub(pattern = "argmaxes.",replacement="n=",.data2$variable)
  
  figure6<-ggplot(data=.data2,aes(x=x,y=y,color=objective,fill=objective))+
    geom_point(colour="black",shape=21,size=2)+ 
    scale_fill_gradient(low = "#FFFFFF",high = "#000000")+
    geom_point(data=.data2[.data2$value,],mapping=aes(x=x,y=y),shape=10,size=4,color="black")+
    facet_grid(vars(rep),vars(variable))+ theme(legend.position="none",axis.title.x=element_blank(),
                                                axis.text.x=element_blank(),
                                                axis.ticks.x=element_blank(),,axis.title.y=element_blank(),
                                                axis.text.y=element_blank(),
                                                axis.ticks.y=element_blank())
  return(list(figure4=figure4,figure5=figure5,figure6=figure6))
}






figures7..9<-function(){
  
  library(ggplot2)
  Nh=4
  H=4
  N=Nh*H
  nrep=2
  set.seed(1)
  LL<-plyr::alply(1:nrep,1,function(Rep){
    #N<-rpois(1,N)
    X=matrix(runif(H*2),H,2)
    X=do.call(rbind,plyr::alply(X,1,function(x){matrix(x[rep(1:2,each=Nh)]+.1*runif(Nh*2),Nh,2)}))
    colnames(X)<-c("x","y")
    XX<-exp(-4*as.matrix(dist(X,upper=T,diag = T)))
    XX[(row(XX)-1)%/%(N/2)!=(col(XX)-1)%/%(N/2)]<-0
    objective<-apply(XX,1,mean)
    optimalsamples<-lapply(1:5,function(n){
      subsets_n<-subsets(n,N)
      prob<-plyr::aaply(subsets_n,1,function(x){
        mean(1-plyr::aaply(1-XX[,x],1,prod))
      })
      argmax=subsets_n[which(prob==max(prob)),]
    })
    argmaxes<-do.call(cbind,lapply(optimalsamples,function(x){is.element(1:N,x)}))
    argmax=(objective==max(objective))
    
    # infected=plyr::raply(10000,
    #                       .expr = (function(){infected0=  infected0=is.element(1:N,sample(N,1));
    #                       cbind(infected0,infected1=rbinom(1:N,1,prob=XX[infected0,]))})())
    list(.data=data.frame(id=1:N,owner=as.factor((0:(N-1))%/%(N/2)),X,
                          objective=objective,
                          argmaxes=argmaxes,
                          rep=paste0("Population ",Rep)#,
                          #meaninfected0=plyr::aaply(infected[,,1],2,mean),
                          #meaninfected1=plyr::aaply(infected[,,2],2,mean)
    ))
  })
  .data=do.call(rbind,lapply(LL, function(L){L$.data}))
  figure7<-ggplot(data=.data,aes(x=x,y=y,color=objective,fill=objective))+
    geom_point(colour="black",shape=21,size=4)+ 
    scale_fill_gradient(low = "#FFFFFF",high = "#000000")+
    geom_point(data=.data[.data$argmaxes.1==1,],mapping=aes(x=x,y=y),shape=10,size=10,color="black")+
    facet_wrap(~rep)+ theme(legend.position="bottom")
  
  print(figure4)
  
  figure8<-ggplot(data=.data,aes(x=x,y=y,color=objective,fill=objective))+
    geom_point(colour="black",shape=21,size=4)+ 
    scale_fill_gradient(low = "#FFFFFF",high = "#000000")+
    geom_point(data=.data[.data$argmaxes.2==1,],mapping=aes(x=x,y=y),shape=10,size=10,color="black")+
    facet_wrap(~rep)+ theme(legend.position="bottom")
  
  
  .data2<-reshape2::melt(.data[c("id","rep","x","y","owner","objective",paste0("argmaxes.",1:5))],id=c("id","rep","owner","x","y","objective"))
  .data2$variable<-gsub(pattern = "argmaxes.",replacement="n=",.data2$variable)
  
  figure9<-ggplot(data=.data2,aes(x=x,y=y,color=objective,fill=objective,shape=owner))+
    geom_point(colour="black",size=2)+ 
    scale_fill_gradient(low = "#FFFFFF",high = "#000000")+
    geom_point(data=.data2[.data2$value,],mapping=aes(x=x,y=y),shape=10,size=4,color="black")+
    facet_grid(vars(rep),vars(variable))+ theme(legend.position="none",axis.title.x=element_blank(),
                                                axis.text.x=element_blank(),
                                                axis.ticks.x=element_blank(),,axis.title.y=element_blank(),
                                                axis.text.y=element_blank(),
                                                axis.ticks.y=element_blank())
  return(list(figure7=figure7,figure8=figure8,figure9=figure9))
}










figures10<-function(){
  
  library(ggplot2)
  
  N=10
  nrep=2
  set.seed(1)
  LLL<-plyr::alply(1:nrep,1,function(Rep){
    #N<-rpois(1,N)
    X=matrix(runif(N*2),N,2)
    colnames(X)<-c("x","y")
    DD<-as.matrix(dist(X,upper=T,diag = T))
    XX<-exp(-4*DD)
    
    
    XXy<-function(y,.N=N,.XX=XX){
      plyr::aaply(subsets(y,.N),1,function(x){
        1-plyr::aaply(1-.XX[x,,drop=FALSE],2,prod)
      })}
    #any(abs(XXy(1)-XX)>1e-10*max(XX))
    #XXy(2)
    
    LL<-plyr::alply(1:3,1,function(y){
      XXY<-XXy(y)       
      objective<-apply(XXY,2,mean)
      optimalsamples<-lapply(1:5,function(n){
        subsets_n<-subsets(n,N)
        prob<-plyr::aaply(subsets_n,1,function(x){
          mean(1-plyr::aaply(1-XXY[,x],1,prod))
        })
        maxprob=max(prob)
        list(maxprob=maxprob,
             prob=prob,
             optimalsample=subsets_n[which(prob==maxprob),])
      })
      
      
      argmaxes<-do.call(cbind,lapply(optimalsamples,function(x){is.element(1:N,x$optimalsample)}))
      
      #infected=plyr::raply(10000,
      #                     .expr = (function(){infected0=  infected0=is.element(1:N,sample(N,1));
      #                     cbind(infected0,infected1=rbinom(1:N,1,prob=XX[infected0,]))})())
      list(points=data.frame(id=1:N,X,
                             initiallyinfected=paste0("Init: ",y),
                             objective=objective,
                             argmaxes=argmaxes,
                             rep=paste0("Pop. ",Rep)),
           probdist=do.call(rbind,lapply(1:5,function(n){
             data.frame(
               y=y,
               rep=Rep,
               n=n,
               prob=optimalsamples[[n]]$prob,
               srs=mean(optimalsamples[[n]]$prob),
               maxprob=max(optimalsamples[[n]]$prob),
               minprob=min(optimalsamples[[n]]$prob))}))
      )
    })
    list(points=do.call(rbind,lapply(LL,function(ll){ll$points})),
         probs=do.call(rbind,lapply(LL,function(ll){ll$probdist})))
  })
  .data=do.call(rbind,lapply(LLL,function(lll){lll$points}))
  
  .data2<-reshape2::melt(.data[c("id","rep","x","y","objective","initiallyinfected",paste0("argmaxes.",1:5))],
                         id=c("id","rep","x","y","objective","initiallyinfected"))
  .data2$variable<-gsub(pattern = "argmaxes.",replacement="n=",.data2$variable)
  
  figure10<-ggplot(data=.data2,aes(x=x,y=y,color=objective,fill=objective))+
    geom_point(colour="black",shape=21,size=2)+ 
    scale_fill_gradient(low = "#FFFFFF",high = "#000000")+
    geom_point(data=.data2[.data2$value,],mapping=aes(x=x,y=y),shape=10,size=4,color="black")+
    facet_grid(vars(rep,initiallyinfected),vars(variable))+ theme(legend.position="none",axis.title.x=element_blank(),
                                                                  axis.text.x=element_blank(),
                                                                  axis.ticks.x=element_blank(),,axis.title.y=element_blank(),
                                                                  axis.text.y=element_blank(),
                                                                  axis.ticks.y=element_blank())
  
  .dataprob=do.call(rbind,lapply(LLL,function(lll){lll$probs}))
  
  .dataprob$improvement=with(.dataprob,prob/srs)
  .dataprob$maximprovement=with(.dataprob,maxprob/srs)
  .dataprob$minimprovement=with(.dataprob,minprob/srs)
  .dataprob$initiallyinfected=with(.dataprob,paste0("Init: ",y))
  .dataprob$samplesize=with(.dataprob,paste0("$n=",n,"$"))
  .dataprob$srs=1
  figure11<-
    ggplot(data=.dataprob,aes(x=improvement))+
    geom_histogram(binwidth = function(x){(max(x)-1)/5.5},center=1)+ 
    geom_vline(xintercept = 1,color="red")+
    geom_vline(aes(xintercept = maximprovement),colour="green")+
    #geom_vline(aes(xintercept = minimprovement),colour="blue")+
    facet_grid(vars(rep,initiallyinfected),vars(samplesize),scales = "free_y")+
    theme(legend.position="none",
          axis.title.x=element_blank(),
          #axis.text.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    scale_x_continuous(breaks = c(.8,1,1.2))
  
  
  if(FALSE){
    figure11+stat_density(trim=TRUE)
    figure11+geom_histogram(bins=10)
    figure11+geom_histogram(bins=10,breaks=function(x){min(x)+(0:10)*(max(x)-min(x))/10})
    figure11+geom_histogram(binwidth = function(x){(max(x)-min(x))/10})
    figure11+geom_histogram(binwidth = function(x){(max(x)-min(x))/10},breaks=function(x){min(x)+(0:10)*(max(x)-min(x))/10})
    figure11+geom_histogram(binwidth = function(x){(max(x)-1)/5},center=1)
    figure11+geom_histogram(center=(min(prob)+(max(prob)-min(prob))/20),binwidth = function(x){(max(x)-min(x))/10})
    figure11+geom_density(trim=TRUE)
    figure11+geom_histogram(stat = "density",adjust=0.05)
    figure11+geom_histogram(stat = "density")
    figure11+geom_density(bw=.003)}
  
  tab1<-unique(.dataprob[c("rep","n","y","srs","minprob","maxprob","minimprovement","maximprovement")])
  return(list(figure10=figure10,figure11=figure11,tab1=tab1))
}





figures12<-function(){
  
  library(ggplot2)
  
  N=30
  set.seed(1)
  #N<-rpois(1,N)
  X=matrix(runif(N*2),N,2)
  colnames(X)<-c("x","y")
  DD<-as.matrix(dist(X,upper=T,diag = T))
  XX<-exp(-4*DD)
  
  
  XXy<-function(y,.N=N,.XX=XX){
    plyr::aaply(subsets(y,.N),1,function(x){
      1-plyr::aaply(1-.XX[x,,drop=FALSE],2,prod)
    })}
  #any(abs(XXy(1)-XX)>1e-10*max(XX))
  #XXy(2)
  
  y=1
  XXY<-XXy(y)       
  objective<-apply(XXY,2,mean)
  n=4
  
  subsets_n<-subsets(n,N)
  prob<-plyr::aaply(subsets_n,1,function(x){
    mean(1-plyr::aaply(1-XXY[,x],1,prod))})
  
  maxprob=max(prob)
  x=list(maxprob=maxprob,
         prob=prob,
         optimalsample=subsets_n[which(prob==maxprob),])
  argmaxes<-is.element(1:N,x$optimalsample)
  
  points=data.frame(id=1:N,X,
                    initiallyinfected=paste0("Init: ",y),
                    objective=objective,
                    argmaxes=argmaxes)
  probdist=
    data.frame(
      y=y,
      n=n,
      prob=prob,
      srs=mean(prob),
      maxprob=max(prob),
      minprob=min(prob))
  
  
  Pairs<-function(Samp,U){do.call(rbind,lapply(Samp,function(x){cbind(x,U)}))}
  
  pairswithoptimal<-Pairs((1:N)[argmaxes],1:N)
  prob2<-as.data.frame(cbind(pairswithoptimal,plyr::aaply(pairswithoptimal,1,function(xx){
    mean(plyr::aaply(unname(XXY[,xx]),1,prod))/
      mean(XXY[,xx[1]])
  })))
  names(prob2)<-c("cluster","k","prob")
  groups<-plyr::ddply(.data = prob2,"k",
                      function(d){
                        data.frame(cluster=d$cluster[d$prob==max(d$prob)])})
  points$cluster=factor(groups$cluster,levels=unique(groups$cluster)[c(2:1,4:3)],ordered = TRUE)
  points$x<-with(points,(x-mean(x))/sd(x))
  points$y<-with(points,(y-mean(y))/sd(y))
  
  
  library(tidyverse)  # data manipulation
  library(cluster)    # clustering algorithms
  library(factoextra)
  set.seed(2)
  k2<-kmeans(X,centers=4,nstart=25)
  
  figure12.b<-fviz_cluster(k2, data = X)+
    ggtitle("K-means clustering")+ 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    geom_point(data=points[points$argmaxes,],mapping=aes(x=x,y=y),shape=10,size=4,color="black")
  
  
  
  
  figure12.a<-
    ggpubr::ggscatter(points, "x", "y", color = "cluster",point =TRUE,ellipse = TRUE, ellipse.type = "convex", 
                      ellipse.level = 0.95, ellipse.alpha = 0.2)+
    theme_gray()+ 
    geom_point(data=points[points$argmaxes,],mapping=aes(x=x,y=y),shape=10,size=4,color="black")+
    ggtitle("Clustering around sampled units")+ 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  
  grid.arrange(figure12.a, figure12.b, nrow = 1)
  
  return(list(figure12.a=figure12.a,figure12.b=figure12.b))
}




figures13<-function(){
  
  library(ggplot2)
  
  N=36
  X=expand.grid(x=(0:5)/5,y=(0:5)/5)
  DD<-as.matrix(dist(X,upper=T,diag = T))
  XX<-exp(-4*DD)
  
  
  XXy<-function(y,.N=N,.XX=XX){
    plyr::aaply(subsets(y,.N),1,function(x){
      1-plyr::aaply(1-.XX[x,,drop=FALSE],2,prod)
    })}
  #any(abs(XXy(1)-XX)>1e-10*max(XX))
  #XXy(2)
  
  y=1
  XXY<-XXy(y)       
  objective<-apply(XXY,2,mean)
  n=4
  
  subsets_n<-subsets(n,N)
  prob<-plyr::aaply(subsets_n,1,function(x){
    mean(1-plyr::aaply(1-XXY[,x],1,prod))})
  
  maxprob=max(prob)
  x=list(maxprob=maxprob,
         prob=prob,
         optimalsample=subsets_n[which(prob==maxprob)[1],])
  argmaxes<-is.element(1:N,x$optimalsample)
  
  points=data.frame(id=1:N,X,
                    initiallyinfected=paste0("Init: ",y),
                    objective=objective,
                    argmaxes=argmaxes)
  probdist=
    data.frame(
      y=y,
      n=n,
      prob=prob,
      srs=mean(prob),
      maxprob=max(prob),
      minprob=min(prob))
  
  
  Pairs<-function(Samp,U){as.matrix(expand.grid(Samp,U))}
  
  pairswithoptimal<-Pairs((1:N)[argmaxes],1:N)
  prob2<-as.data.frame(cbind(pairswithoptimal,plyr::aaply(pairswithoptimal,1,function(xx){
    mean(plyr::aaply(unname(XXY[,xx]),1,prod))/
      mean(XXY[,xx[1]])
  })))
  names(prob2)<-c("cluster","k","prob")
  groups<-plyr::ddply(.data = prob2,"k",
                      function(d){data.frame(cluster=d$cluster[d$prob==max(d$prob)])})
  points$cluster=factor(groups$cluster,levels=unique(groups$cluster)[c(2:1,4:3)],ordered = TRUE)
  points$x<-with(points,(x-mean(x))/sd(x))
  points$y<-with(points,(y-mean(y))/sd(y))
  
  
  library(tidyverse)  # data manipulation
  library(cluster)    # clustering algorithms
  library(factoextra)
  set.seed(2)
  k2<-kmeans(X,centers=4,nstart=25)
  
  figure13.b<-fviz_cluster(k2, data = X)+
    ggtitle("K-means clustering")+ 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())+
    geom_point(data=points[points$argmaxes,],mapping=aes(x=x,y=y),shape=10,size=4,color="black")
  
  
  
  
  figure13.a<-
    ggpubr::ggscatter(points, "x", "y", color = "cluster",point =TRUE,ellipse = TRUE, ellipse.type = "convex", 
                      ellipse.level = 0.95, ellipse.alpha = 0.2)+
    theme_gray()+ 
    geom_point(data=points[points$argmaxes,],mapping=aes(x=x,y=y),shape=10,size=4,color="black")+
    ggtitle("Clustering around sampled units")+ 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  
  grid.arrange(figure13.a, figure13.b, nrow = 1)
  
  return(list(figure13.a=figure13.a,figure13.b=figure13.b))
}





approxoptimalsample<-function(X,n,M=1000){
  N=nrow(X)
  DD<-as.matrix(dist(X,upper=T,diag = T))
  XX<-exp(-4*DD)
  y=1
  XXy<-function(y,.N=N,.XX=XX){
    plyr::aaply(subsets(y,.N),1,function(x){
      1-plyr::aaply(1-.XX[x,,drop=FALSE],2,prod)
    })}
  XXY<-XXy(y)       
  k=vector()
  for(i in 1:n){
    prob<-plyr::aaply(setdiff(1:N,k),1,function(x){
      mean(1-plyr::aaply(1-XXY[,c(x,k),drop=FALSE],1,prod))})    
    k<-c(k,which(prob==max(prob))[1])}
  k
  alpha=0
  beta=0
  while(alpha<n & beta<M){
    prob<-plyr::aaply(setdiff(1:N,unname(k)[2:n]),1,function(x){
      mean(1-plyr::aaply(1-XXY[,c(x,unname(k)[2:n]),drop=FALSE],1,prod))})
    kprime<-unname((setdiff(1:N,unname(k)[2:n]))[which(prob==max(prob))])
    alpha<-(alpha+1)*(is.element(k[1],kprime))
    beta<-beta+1
    k<-unname(c(k[2:n],kprime[1]))
  }
  k
}

if(FALSE){
subsets_n<-subsets(n,N)
prob<-plyr::aaply(subsets_n,1,function(x){
  mean(1-plyr::aaply(1-XXY[,x],1,prod))})


#any(abs(XXy(1)-XX)>1e-10*max(XX))
#XXy(2)

objective<-apply(XXY,2,mean)
n=4


}




figures14<-function(){
  
  library(ggplot2)
  N=36
  X=expand.grid(x=(0:5)/5,y=(0:5)/5)
  DD<-as.matrix(dist(X,upper=T,diag = T))
  XX<-exp(-4*DD)
  
  
  XXy<-function(y,.N=N,.XX=XX){
    plyr::aaply(subsets(y,.N),1,function(x){
      1-plyr::aaply(1-.XX[x,,drop=FALSE],2,prod)
    })}
  #any(abs(XXy(1)-XX)>1e-10*max(XX))
  #XXy(2)
  
  y=1
  XXY<-XXy(y)       
  objective<-apply(XXY,2,mean)
  n=4
  
  subsets_n<-subsets(n,N)
  prob<-plyr::aaply(subsets_n,1,function(x){
    mean(1-plyr::aaply(1-XXY[,x],1,prod))})
  
  maxprob=max(prob)
  x=list(maxprob=maxprob,
         prob=prob,
         optimalsample=subsets_n[which(prob==maxprob)[1],])
  argmaxes<-is.element(1:N,x$optimalsample)
  approxoptimal<-approxoptimalsample(X,4)
  
  .points=data.frame(id=1:N,X,
                     initiallyinfected=paste0("Init: ",y),
                     objective=objective,
                     argmaxes=argmaxes,
                     argmaxes2=is.element(1:N,approxoptimal))
  probdist=
    data.frame(
      y=y,
      n=n,
      prob=prob,
      srs=mean(prob),
      maxprob=max(prob),
      minprob=min(prob))
  
  
  
  Pairs<-function(Samp,U){as.matrix(expand.grid(Samp,U))}
  
  pairswithoptimal<-Pairs((1:N)[argmaxes],1:N)
  pairswithoptimal2<-Pairs(approxoptimal,1:N)
  prob2<-as.data.frame(cbind(pairswithoptimal,plyr::aaply(pairswithoptimal,1,function(xx){
    mean(plyr::aaply(unname(XXY[,xx]),1,prod))/
      mean(XXY[,xx[1]])
  })))
  prob3<-as.data.frame(cbind(pairswithoptimal2,plyr::aaply(pairswithoptimal2,1,function(xx){
    mean(plyr::aaply(unname(XXY[,xx]),1,prod))/
      mean(XXY[,xx[1]])
  })))
  names(prob2)<-c("cluster","k","prob")
  names(prob3)<-c("cluster","k","prob")
  groups<-plyr::ddply(.data = prob2,"k",
                      function(d){data.frame(cluster=d$cluster[d$prob==max(d$prob)])})
  .points$cluster=factor(groups$cluster,levels=unique(groups$cluster)[c(2:1,4:3)],ordered = TRUE)
  .points$x<-with(.points,(x-mean(x))/sd(x))
  .points$y<-with(.points,(y-mean(y))/sd(y))
  
  
  groups2<-plyr::ddply(.data = prob3,"k",
                       function(d){data.frame(cluster=d$cluster[d$prob==max(d$prob)])})
  .points$cluster2=factor(groups2$cluster,levels=unique(groups2$cluster),ordered = TRUE)
  
  
  
  library(tidyverse)  # data manipulation
  library(cluster)    # clustering algorithms
  library(factoextra)
  set.seed(2)
  k2<-kmeans(X,centers=4,nstart=25)
  
  
  figure14.b<-
    ggpubr::ggscatter(.points, "x", "y", color = "cluster2",point =TRUE,ellipse = TRUE, ellipse.type = "convex", 
                      ellipse.level = 0.95, ellipse.alpha = 0.2)+
    theme_gray()+ 
    geom_point(data=.points[.points$argmaxes2,],mapping=aes(x=x,y=y),shape=10,size=4,color="black")+
    ggtitle("Clustering around approximate optimal sample")+ 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  
  
  
  figure14.a<-
    ggpubr::ggscatter(.points, "x", "y", color = "cluster",point =TRUE,ellipse = TRUE, ellipse.type = "convex", 
                      ellipse.level = 0.95, ellipse.alpha = 0.2)+
    theme_gray()+ 
    geom_point(data=.points[.points$argmaxes,],mapping=aes(x=x,y=y),shape=10,size=4,color="black")+
    ggtitle("Clustering around sampled units")+ 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  
  grid.arrange(figure14.a, figure14.b, nrow = 1)
  
  return(list(figure14.a=figure14.a,
              figure14.b=figure14.b))
}

figures15<-function(){
  
  library(ggplot2)
  N=30
  set.seed(1)
  #N<-rpois(1,N)
  X=matrix(runif(N*2),N,2)
  colnames(X)<-c("x","y")
  
  DD<-as.matrix(dist(X,upper=T,diag = T))
  XX<-exp(-4*DD)
  
  
  XXy<-function(y,.N=N,.XX=XX){
    plyr::aaply(subsets(y,.N),1,function(x){
      1-plyr::aaply(1-.XX[x,,drop=FALSE],2,prod)
    })}
  #any(abs(XXy(1)-XX)>1e-10*max(XX))
  #XXy(2)
  
  y=1
  XXY<-XXy(y)       
  objective<-apply(XXY,2,mean)
  n=4
  
  subsets_n<-subsets(n,N)
  prob<-plyr::aaply(subsets_n,1,function(x){
    mean(1-plyr::aaply(1-XXY[,x],1,prod))})
  
  maxprob=max(prob)
  x=list(maxprob=maxprob,
         prob=prob,
         optimalsample=subsets_n[which(prob==maxprob)[1],])
  argmaxes<-is.element(1:N,x$optimalsample)
  approxoptimal<-approxoptimalsample(X,4)
  
  .points=data.frame(id=1:N,X,
                     initiallyinfected=paste0("Init: ",y),
                     objective=objective,
                     argmaxes=argmaxes,
                     argmaxes2=is.element(1:N,approxoptimal))
  probdist=
    data.frame(
      y=y,
      n=n,
      prob=prob,
      srs=mean(prob),
      maxprob=max(prob),
      minprob=min(prob))
  
  
  
  Pairs<-function(Samp,U){as.matrix(expand.grid(Samp,U))}
  
  pairswithoptimal<-Pairs((1:N)[argmaxes],1:N)
  pairswithoptimal2<-Pairs(approxoptimal,1:N)
  prob2<-as.data.frame(cbind(pairswithoptimal,plyr::aaply(pairswithoptimal,1,function(xx){
    mean(plyr::aaply(unname(XXY[,xx]),1,prod))/
      mean(XXY[,xx[1]])
  })))
  prob3<-as.data.frame(cbind(pairswithoptimal2,plyr::aaply(pairswithoptimal2,1,function(xx){
    mean(plyr::aaply(unname(XXY[,xx]),1,prod))/
      mean(XXY[,xx[1]])
  })))
  names(prob2)<-c("cluster","k","prob")
  names(prob3)<-c("cluster","k","prob")
  groups<-plyr::ddply(.data = prob2,"k",
                      function(d){data.frame(cluster=d$cluster[d$prob==max(d$prob)])})
  .points$cluster=factor(groups$cluster)
  .points$x<-with(.points,(x-mean(x))/sd(x))
  .points$y<-with(.points,(y-mean(y))/sd(y))
  
  
  groups2<-plyr::ddply(.data = prob3,"k",
                       function(d){data.frame(cluster=d$cluster[d$prob==max(d$prob)])})
  .points$cluster2=factor(groups2$cluster)
  
  
  
  library(tidyverse)  # data manipulation
  library(cluster)    # clustering algorithms
  library(factoextra)
  set.seed(2)
  
  
  p1<-ggplot(as.data.frame(X), aes(x=x, y=y) ) +xlim(c(-.25,1.25))+ylim(c(-.25,1.25))+
    stat_density_2d(aes(fill = ..level..), geom = "polygon",h=.25)+
    geom_point(col="black")
  
  p2<-ggplot(as.data.frame(X), aes(x=x, y=y) ) +xlim(c(-.25,1.25))+ylim(c(-.25,1.25))+
    stat_density_2d(data=as.data.frame(X[argmaxes,]),aes(x=x,y=y,fill = ..level..), geom = "polygon",h=.25)+
    geom_point(col="black")
  
  grid.arrange(p1,p2,nrow=1)
  
  figure15.b<-
    ggpubr::ggscatter(.points, "x", "y", color = "cluster2",point =TRUE,ellipse = TRUE, ellipse.type = "convex", 
                      ellipse.level = 0.95, ellipse.alpha = 0.2)+
    theme_gray()+ 
    geom_point(data=.points[.points$argmaxes2,],mapping=aes(x=x,y=y),shape=10,size=4,color="black")+
    ggtitle("Clustering around approximate optimal sample")+ 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  
  
  
  figure15.a<-
    ggpubr::ggscatter(.points, "x", "y", color = "cluster",point =TRUE,ellipse = TRUE, ellipse.type = "convex", 
                      ellipse.level = 0.95, ellipse.alpha = 0.2)+
    theme_gray()+ 
    geom_point(data=.points[.points$argmaxes,],mapping=aes(x=x,y=y),shape=10,size=4,color="black")+
    ggtitle("Clustering around sampled units")+ 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  
  grid.arrange(figure15.a, figure15.b, nrow = 1)
  
  return(list(figure15.a=figure15.a,
              figure15.b=figure15.b))
}






figure15.a<-function(){
  ggpubr::ggscatter(.points, "x", "y", color = "cluster",point =TRUE,ellipse = TRUE, ellipse.type = "convex", 
                    ellipse.level = 0.95, ellipse.alpha = 0.2)+
  theme_gray()+ 
  geom_point(data=.points[.points$argmaxes,],mapping=aes(x=x,y=y),shape=10,size=4,color="black")+
  ggtitle("Clustering around sampled units")+ 
  theme(legend.position="none",
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())


grid.arrange(figure15.a, figure15.b, nrow = 1)

return(list(figure15.a=figure15.a,
            figure15.b=figure15.b))
}







#quizz
figures16<-function(){
  
  library(ggplot2)
  N=12
  i=10
  i=i+1
  set.seed(i)
  #N<-rpois(1,N)
  X=matrix(runif(N*2),N,2)
  colnames(X)<-c("x","y")
  
  DD<-as.matrix(dist(X,upper=T,diag = T))
  XX<-exp(-4*DD)
  
  
  XXy<-function(y,.N=N,.XX=XX){
    plyr::aaply(subsets(y,.N),1,function(x){
      1-plyr::aaply(1-.XX[x,,drop=FALSE],2,prod)
    })}
  #any(abs(XXy(1)-XX)>1e-10*max(XX))
  #XXy(2)
  
  y=1
  XXY<-XXy(y)       
  objective<-apply(XXY,2,mean)
  n=2
  
  subsets_n<-subsets(n,N)
  prob<-plyr::aaply(subsets_n,1,function(x){
    mean(1-plyr::aaply(1-XXY[,x],1,prod))})
  
  maxprob=max(prob)
  x=list(maxprob=maxprob,
         prob=prob,
         optimalsample=subsets_n[which(prob==maxprob)[1],])
  argmaxes<-is.element(1:N,x$optimalsample)
  approxoptimal<-approxoptimalsample(X,2)
  
  .points=data.frame(id=1:N,X,
                     initiallyinfected=paste0("Init: ",y),
                     objective=objective,
                     argmaxes=argmaxes,
                     argmaxes2=is.element(1:N,approxoptimal))
  probdist=
    data.frame(
      y=y,
      n=n,
      prob=prob,
      srs=mean(prob),
      maxprob=max(prob),
      minprob=min(prob))
  
  
  
  Pairs<-function(Samp,U){as.matrix(expand.grid(Samp,U))}
  
  pairswithoptimal<-Pairs((1:N)[argmaxes],1:N)
  pairswithoptimal2<-Pairs(approxoptimal,1:N)
  prob2<-as.data.frame(cbind(pairswithoptimal,plyr::aaply(pairswithoptimal,1,function(xx){
    mean(plyr::aaply(unname(XXY[,xx]),1,prod))/
      mean(XXY[,xx[1]])
  })))
  prob3<-as.data.frame(cbind(pairswithoptimal2,plyr::aaply(pairswithoptimal2,1,function(xx){
    mean(plyr::aaply(unname(XXY[,xx]),1,prod))/
      mean(XXY[,xx[1]])
  })))
  names(prob2)<-c("cluster","k","prob")
  names(prob3)<-c("cluster","k","prob")
  groups<-plyr::ddply(.data = prob2,"k",
                      function(d){data.frame(cluster=d$cluster[d$prob==max(d$prob)])})
  .points$cluster=factor(groups$cluster)
  .points$x<-with(.points,(x-mean(x))/sd(x))
  .points$y<-with(.points,(y-mean(y))/sd(y))
  
  
  groups2<-plyr::ddply(.data = prob3,"k",
                       function(d){data.frame(cluster=d$cluster[d$prob==max(d$prob)])})
  .points$cluster2=factor(groups2$cluster)
  
  
  
  library(tidyverse)  # data manipulation
  library(cluster)    # clustering algorithms
  library(factoextra)
  set.seed(2)
  
  
  
  
  results=data.frame(subsets_n,prob=prob)
  results<-merge(results,.points[,1:3],by.x="X1",by.y="id")
  results<-merge(results,.points[,1:3],by.x="X2",by.y="id")
  
  
  
  alpha=10
  z<-strsplit(
    "(5,8)
(5,12)
(10,11)
(10,11)
(10,12)
(10,12)
(1,8)
(1,5)",split="
")[[1]]
  z<-gsub(")","",z)
  z<-gsub("(","",z,fixed = T)
  answers<-data.frame(t(sapply(strsplit(z,","),strtoi)))
  answers<-data.frame(table(answers))
  answers<-answers[answers$Freq!=0,]
  answers$X1<-levels(answers$X1)[answers$X1]
  answers$X2<-levels(answers$X2)[answers$X2]
  answers<-merge(answers,results)
  
  results$prob2<-with(results,(1+(is.element(X1,c(6,9))*!is.element(X2,c(6,9))+is.element(X2,c(6,9))*!is.element(X1,c(6,9))))/2)
  
  
  figure16.a<-
    ggplot(data=.points,aes(x=x,y=y,label=id))+geom_point()+
    theme_gray()+
    geom_text(hjust=-1,vjust=.5)+ 
    ggtitle("Trees")+ 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  
  
  figure16.b<-
    ggpubr::ggscatter(.points, "x", "y", color = "cluster",label = "id",point =TRUE,ellipse = TRUE, ellipse.type = "convex", 
                      ellipse.level = 0.95, ellipse.alpha = 0.2)+
    theme_gray()+ 
    geom_point(data=.points[.points$argmaxes,],mapping=aes(x=x,y=y),shape=10,size=4,color="black")+
    ggtitle("One solution")+ 
    theme(legend.position="none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  
  figure16.c<-
    ggplot(data=answers,aes(x=x.x,y=y.x))+
    geom_line(data=.points[.points$argmaxes,],mapping=aes(x=x,y=y),size=4,color="red")+
    geom_segment(aes(x = x.x, y = y.x, xend = x.y, yend = y.y,size =Freq),lineend = "round",color="black",alpha=0.2)+
    geom_point(data=.points,aes(x=x,y=y))+
    theme_gray()+
    geom_text(data=.points,aes(x=x,y=y,label=id),hjust=-1,vjust=.5)+ 
    geom_point(data=.points[.points$argmaxes,],
               mapping=aes(x=x,y=y),shape=10,size=4,color="black")+
    ggtitle("Answers")+ theme(axis.title.x=element_blank(),
                              axis.text.x=element_blank(),
                              axis.ticks.x=element_blank(),
                              axis.title.y=element_blank(),
                              axis.text.y=element_blank(),
                              axis.ticks.y=element_blank())
  
  
  
  figure16.d<-
    ggplot(rbind(data.frame(prob=results$prob,Design="At random"),data.frame(prob=answers$prob,Design="Team selection")),aes(x=prob,color=Design))+stat_ecdf()+
    theme_gray()+ 
    ggtitle("CDF of prior probability for 2-sized samples")+ylab("")+xlab("")+
    theme(legend.position="bottom")
  
  
  figure16.e<-
    ggplot(data=results,aes(x=x.x,y=y.x))+
    geom_point()+
    geom_segment(aes(x = x.x, y = y.x, xend = x.y, yend = y.y,color =prob,alpha=prob))+
    scale_colour_gradient2(low = "gray",
                           mid = "yellow",
                           high = "red",midpoint = mean(range(results$prob)))+
    geom_text(data=.points,aes(x=x,y=y,label=id),hjust=-1,vjust=.5)+ 
    ggtitle("Answers")+ theme(axis.title.x=element_blank(),
                              axis.text.x=element_blank(),
                              axis.ticks.x=element_blank(),
                              axis.title.y=element_blank(),
                              axis.text.y=element_blank(),
                              axis.ticks.y=element_blank())
  
  
  figure16.f<-
    ggplot(data=results[order(results$prob2),],aes(x=x.x,y=y.x))+
    geom_point()+
    geom_segment(aes(x = x.x, y = y.x, xend = x.y, yend = y.y,color =prob2))+
    scale_colour_gradient2(low = "gray",
                           mid = "yellow",
                           high = "red",midpoint = .75)+
    geom_text(data=.points,aes(x=x,y=y,label=id),hjust=-1,vjust=.5)+ 
    ggtitle("Optimal samples")+ theme(legend.position="none",
                                      axis.title.x=element_blank(),
                                      axis.text.x=element_blank(),
                                      axis.ticks.x=element_blank(),
                                      axis.title.y=element_blank(),
                                      axis.text.y=element_blank(),
                                      axis.ticks.y=element_blank())
  
  
  figure16<-grid.arrange(figure16.c,figure16.a, figure16.b, nrow = 1)
  
  
  
  return(list(figure16=figure16,
              results=results,
              answers=answers,
              figure16.a=figure16.a,
              figure16.b=figure16.b,
              figure16.c=figure16.c,
              figure16.d=figure16.d,
              figure16.e=figure16.e,
              figure16.f=figure16.f))
}

if(FALSE){

library(plot3D)
scatter3D(x=results$x.x, y=results$y.x, z=results$prob)
text3D(x, y, z, labels, colvar = NULL, add = FALSE)
points3D(x, y, z, ...)
lines3D(x=results$x.x, y=results$y.x, z=prob)+
  lines3D(x=results$x.y, y=results$y.y, z=prob)
segments3D(x0=results$x.x, y0=results$y.x, z0=prob,
           x1=results$x.y, y1=results$y.y, z1=prob)
}
