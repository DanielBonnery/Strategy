Australia_avo_fields_small_dist_matrix_f<-function(delta=Inf){
  data("Avo_fields",package="Strategy")
newprojAvo_fields<-spTransform(Avo_fields,CRS("+init=epsg:3107 +units=m"))
list.poly<-Strategy::extractpolygonsaslist(newprojAvo_fields)
gradients=cbind(c(0,1),c(1,0),c(1,1),c(1,-1))
polysmalldistmat(list.poly,delta,gradients=gradients)
}


Australia_connected_fields_f<-function(delta=200){
  data(Avo_fields,package="Strategy")
  data(Australia_avo_fields_small_dist_matrix,package="Strategy")
  connectedpop(Australia_avo_fields_small_dist_matrix,delta,nrow(Avo_fields))  
}

