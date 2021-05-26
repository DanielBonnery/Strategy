approximatedrisk1<-function(N,m,n0,n1,n2,betac){(1-(m/N)*(1-betac^(n2)))^(n1*n0)}
requirednumberofbulks<-function(N,m,n1,n2,beta,risk){
  ceiling(log(risk)/(n1*log((1-(m/N)*(1-beta^n2)))))
}


parameters<-expand.grid(
  N=32903,
  m=c(10,30,100,300),
  beta=c(.5,.8),
  n1=c(25,50),
  n2=1:2,
  acceptablerisk=c(.01,.05)
)
Results<-plyr::adply(parameters,1,function(x){
  n0=requirednumberofbulks(x$N,x$m,x$n1,x$n2,x$beta,x$acceptablerisk)
})%>%dplyr::filter(n1*n2==50)
