# Simulations for desease detection

## Installation
To install  the package, execute:

```r
cred <- git2r::cred_user_pass("username", "password")
devtools::install_git("http://pine--is.grid.private.cam.ac.uk:8888/gitlab/dbb31/Strategy.git", credentials = cred)
```
where '"username"' and '"password"' need to be replaced by correct strings.


## Generate a population that matches the ONS 2006 population estimates for the parishes of Cumbria

## Details and demo


```r
library(Strategy)
example(addpiechartclustermarkers)
```

```
FALSE Loading required package: ggplot2
```

```
FALSE Loading required package: leaflet
```

```
FALSE Loading required package: spatstat
```

```
FALSE Loading required package: spatstat.data
```

```
FALSE Loading required package: nlme
```

```
FALSE Loading required package: rpart
```

```
FALSE 
FALSE spatstat 1.63-2       (nickname: 'I'm sorry Dave, I'm afraid I can't do that') 
FALSE For an introduction to spatstat, type 'beginner'
```

```
FALSE Loading required package: sf
```

```
FALSE Linking to GEOS 3.6.2, GDAL 2.2.3, PROJ 4.9.3
```

```
FALSE Loading required package: sp
```

```
FALSE 
FALSE Attaching package: 'dplyr'
```

```
FALSE The following object is masked from 'package:nlme':
FALSE 
FALSE     collapse
```

```
FALSE The following objects are masked from 'package:stats':
FALSE 
FALSE     filter, lag
```

```
FALSE The following objects are masked from 'package:base':
FALSE 
FALSE     intersect, setdiff, setequal, union
```


```r
example(segment.intersect)
```

```
## 
## sgmnt.> zz<-function(){
## sgmnt.+ s1=matrix(sample(0:2,4,rep=T),2,2)
## sgmnt.+ s2=matrix(sample(0:2,4,rep=T),2,2)
## sgmnt.+ si<-segment.intersect(s1,s2)
## sgmnt.+ s<-rbind(s1,s2)
## sgmnt.+ plot(s,cex=.5,main=if(si){"Intersect"}else{"Disjoint"})
## sgmnt.+ points(s1,type="l")
## sgmnt.+ points(s2,type="l")
## sgmnt.+ text(s[,1],s[,2],toupper(letters[1:4]),cex=2,col="red")
## sgmnt.+ 
## sgmnt.+ }
## 
## sgmnt.> par(mfrow=c(3,3))
## 
## sgmnt.> set.seed(1);replicate(9,zz())
```

<img src="figure/unnamed-chunk-3-1.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" width="100%" />

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
## 
## [[3]]
## NULL
## 
## [[4]]
## NULL
## 
## [[5]]
## NULL
## 
## [[6]]
## NULL
## 
## [[7]]
## NULL
## 
## [[8]]
## NULL
## 
## [[9]]
## NULL
```



```r
example(triangleorientation)
```

```
## 
## trnglr> zz<-function(p){
## trnglr+ s=matrix(sample(0:2,6,rep=T),3,2)
## trnglr+ plot(s[c(1:3,1),],type="l",main=(if(triangleorientation(s)==1){"+"}else{if(triangleorientation(s)==-1){"-"}else{"aligned"}}))
## trnglr+ text(s[,1],s[,2],toupper(letters[1:3]),cex=2,col="red")
## trnglr+ 
## trnglr+ }
## 
## trnglr> par(mfrow=c(3,3))
## 
## trnglr> replicate(9,zz(c(sample(-2:3,1),sample(-2:2,1))))
```

<img src="figure/unnamed-chunk-4-1.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" width="100%" />

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
## 
## [[3]]
## NULL
## 
## [[4]]
## NULL
## 
## [[5]]
## NULL
## 
## [[6]]
## NULL
## 
## [[7]]
## NULL
## 
## [[8]]
## NULL
## 
## [[9]]
## NULL
```


