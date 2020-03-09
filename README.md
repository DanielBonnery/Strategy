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

This function allows to replace cluster markers by cluster pie charts on leaflet.

```r
library(Strategy)
example(addpiechartclustermarkers)
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


A function two triangles have the same orientation.


```r
example(triangleorientation,echo=F)
```

<img src="figure/unnamed-chunk-4-1.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" width="100%" />

A function to tell if two segments intersect.


```r
example(segment.intersect,echo=F)
```

<img src="figure/unnamed-chunk-5-1.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" width="100%" />




```r
#example(Generate_Discrete_Time_Epidemics,echo=F)
```



```r
example(Generate_U,echo=F)
```

```
## rgdal: version: 1.4-8, (SVN revision 845)
##  Geospatial Data Abstraction Library extensions to R successfully loaded
##  Loaded GDAL runtime: GDAL 2.2.3, released 2017/11/20
##  Path to GDAL shared files: /usr/share/gdal/2.2
##  GDAL binary built with GEOS: TRUE 
##  Loaded PROJ.4 runtime: Rel. 4.9.3, 15 August 2016, [PJ_VERSION: 493]
##  Path to PROJ.4 shared files: (autodetected)
##  Linking to sp version: 1.3-2
```

```
## OGR data source with driver: ESRI Shapefile 
## Source: "/tmp/RtmpBGszg8/Parishes_December_2011_Boundaries_EW_BFC.shp", layer: "Parishes_December_2011_Boundaries_EW_BFC"
## with 11358 features
## It has 7 fields
## Integer64 fields read as strings:  objectid
```




```r
example(risktobeinfectedbydistancetooneinfectedunit)
```

```
## 
## rsktbn> #Risk to be ingfected 2 m from the victim when the 50%risk distance is 1 m:
## rsktbn> risktobeinfectedbydistancetooneinfectedunit(2,1)
## [1] 0.05583301
```




```r
example(dist_areas_f)
```

```
## 
## dst_r_> data(U) 
## 
## dst_r_> dist_areas_f(U)
##              18         67         68         69         70        118
##             119        120        121        122        124        125
##             167        168        169        170        171        172
##             173        175        176        177        186        218
##             219        220        221        222        223        224
##             225        226        227        228        229        231
##             232        233        236        237        238        266
##             267        268        269        270        271        272
##             273        274        275        276        277        278
##             279        280        281        282        283        284
##             285        286        287        288        289        290
##             291        292        317        318        319        320
##             321        322        323        324        325        326
##             327        328        329        330        331        332
##             333        334        335        336        337        338
##             339        340        341        342        343        344
##             345        367        368        369        370        371
##             372        373        374        375        376        377
##             378        379        380        381        382        383
##             384        385        386        387        388        389
##             390        391        392        393        394        395
##             396        397        398        399        417        418
##             419        420        421        422        423        424
##             425        426        427        428        429        430
##             431        432        433        434        435        436
##             437        438        439        440        441        442
##             443        444        445        446        447        448
##             449        450        451        452        453        454
##             467        468        469        470        471        472
##             473        474        475        476        477        478
##             479        480        481        482        483        484
##             485        486        487        488        489        490
##             491        492        493        494        495        496
##             497        498        499        500        501        502
##             503        504        518        519        520        521
##             522        523        524        525        526        527
##             528        529        530        531        532        533
##             534        535        536        537        538        539
##             540        541        542        543        544        545
##             546        547        548        549        550        551
##             552        553        554        555        556        569
##             570        571        572        573        574        575
##             576        577        578        579        580        581
##             582        583        584        585        586        587
##             588        589        590        591        592        593
##             594        595        596        597        598        599
##             600        601        602        603        604        605
##             606        620        621        622        623        624
##             625        626        627        628        629        630
##             631        632        633        634        635        636
##             637        638        639        640        641        642
##             643        644        645        646        647        648
##             649        650        651        652        653        654
##             655        656        657        658        670        671
##             672        673        674        675        676        677
##             678        679        680        681        682        683
##             684        685        686        687        688        689
##             690        691        692        693        694        695
##             696        697        698        699        700        701
##             702        703        704        705        706        707
##             708        721        722        723        724        725
##             726        727        728        729        730        731
##             732        733        734        735        736        737
##             738        739        740        741        742        743
##             744        745        746        747        748        749
##             750        751        752        753        754        755
##             756        757        758        759        771        772
##             773        774        775        776        777        778
##             779        780        781        782        783        784
##             785        786        787        788        789        790
##             791        792        793        794        795        796
##             797        798        799        800        801        802
##             803        804        805        806        807        808
##             809        810        811        822        823        824
##             825        826        827        828        829        830
##             831        832        833        834        835        836
##             837        838        839        840        841        842
##             843        844        845        846        847        848
##             849        850        851        852        853        854
##             855        856        857        858        859        860
##             861        862        872        873        874        875
##             876        877        878        879        880        881
##             882        883        884        885        886        887
##             888        889        890        891        892        893
##             894        895        896        897        898        899
##             900        901        902        903        904        905
##             906        907        908        909        910        911
##             912        913        922        923        924        925
##             926        927        928        929        930        931
##             932        933        934        935        936        937
##             938        939        940        941        942        943
##             944        945        946        947        948        949
##             950        951        952        953        954        955
##             956        957        958        959        960        961
##             962        963        964        965        972        973
##             974        975        976        977        978        979
##             980        981        982        983        984        985
##             986        987        988        989        990        991
##             992        993        994        995        996        997
##             998        999       1000       1001       1002       1003
##            1004       1005       1006       1007       1008       1009
##            1010       1011       1012       1013       1014       1015
##            1016       1017       1018       1019       1023       1024
##            1025       1026       1027       1028       1029       1030
##            1031       1032       1033       1034       1035       1036
##            1037       1038       1039       1040       1041       1042
##            1043       1044       1045       1046       1047       1048
##            1049       1050       1051       1052       1053       1054
##            1055       1056       1057       1058       1059       1060
##            1061       1062       1063       1064       1065       1066
##            1067       1068       1069       1070       1071       1072
##            1073       1074       1075       1076       1077       1078
##            1079       1080       1081       1082       1083       1084
##            1085       1086       1087       1088       1089       1090
##            1091       1092       1093       1094       1095       1096
##            1097       1098       1099       1100       1101       1102
##            1103       1104       1105       1106       1107       1108
##            1109       1110       1111       1112       1113       1114
##            1115       1116       1117       1118       1119       1120
##            1121       1123       1124       1125       1126       1127
##            1128       1129       1130       1131       1132       1133
##            1134       1135       1136       1137       1138       1139
##            1140       1141       1142       1143       1144       1145
##            1146       1147       1148       1149       1150       1151
##            1152       1153       1154       1155       1156       1157
##            1158       1159       1160       1161       1162       1163
##            1164       1165       1166       1167       1168       1169
##            1170       1171       1172       1173       1174       1176
##            1177       1178       1179       1180       1181       1182
##            1183       1184       1185       1186       1187       1188
##            1189       1190       1191       1192       1193       1194
##            1195       1196       1197       1198       1199       1200
##            1201       1202       1203       1204       1205       1206
##            1207       1208       1209       1210       1211       1212
##            1213       1214       1215       1216       1217       1218
##            1219       1220       1221       1222       1223       1228
##            1229       1230       1231       1232       1233       1234
##            1235       1236       1237       1238       1239       1240
##            1241       1242       1243       1244       1245       1246
##            1247       1248       1249       1250       1251       1252
##            1253       1254       1255       1256       1257       1258
##            1259       1260       1261       1262       1263       1264
##            1265       1266       1267       1268       1269       1270
##            1271       1272       1273       1277       1278       1279
##            1280       1281       1282       1283       1284       1285
##            1286       1287       1288       1289       1290       1291
##            1292       1293       1294       1295       1296       1297
##            1298       1299       1300       1301       1302       1303
##            1304       1305       1306       1307       1308       1309
##            1310       1311       1312       1313       1314       1315
##            1316       1317       1318       1319       1320       1321
##            1322       1323       1329       1330       1331       1332
##            1333       1334       1335       1336       1337       1338
##            1339       1340       1341       1342       1343       1344
##            1345       1346       1347       1348       1349       1350
##            1351       1352       1353       1354       1355       1356
##            1357       1358       1359       1360       1361       1362
##            1363       1364       1365       1366       1367       1368
##            1369       1370       1371       1372       1373       1380
##            1381       1382       1383       1384       1385       1386
##            1387       1388       1389       1390       1391       1392
##            1393       1394       1395       1396       1397       1398
##            1399       1400       1401       1402       1403       1404
##            1405       1406       1407       1408       1409       1410
##            1411       1412       1413       1414       1415       1416
##            1417       1418       1419       1420       1421       1422
##            1423       1431       1432       1433       1434       1435
##            1436       1437       1438       1439       1440       1441
##            1442       1443       1444       1445       1446       1447
##            1448       1449       1450       1451       1452       1453
##            1454       1455       1456       1457       1458       1459
##            1460       1461       1462       1463       1464       1465
##            1466       1467       1468       1469       1470       1471
##            1472       1473       1474       1482       1483       1484
##            1485       1486       1487       1488       1489       1490
##            1491       1492       1493       1494       1495       1496
##            1497       1498       1499       1500       1501       1502
##            1503       1504       1505       1506       1507       1508
##            1509       1510       1511       1512       1513       1514
##            1515       1516       1517       1518       1519       1520
##            1521       1522       1523       1524       1525       1533
##            1534       1535       1536       1537       1538       1539
##            1540       1541       1542       1543       1544       1545
##            1546       1547       1548       1549       1550       1551
##            1552       1553       1554       1555       1556       1557
##            1558       1559       1560       1562       1563       1564
##            1565       1566       1567       1568       1569       1570
##            1571       1572       1573       1574       1575       1576
##            1577       1585       1586       1587       1588       1589
##            1590       1591       1592       1593       1594       1595
##            1596       1597       1598       1599       1600       1601
##            1602       1603       1604       1605       1606       1607
##            1608       1609       1610       1613       1614       1615
##            1616       1617       1618       1619       1620       1621
##            1622       1623       1624       1625       1637       1638
##            1639       1640       1641       1642       1643       1644
##            1645       1646       1647       1648       1649       1650
##            1651       1652       1653       1654       1655       1656
##            1657       1658       1659       1660       1661       1662
##            1664       1665       1666       1667       1668       1669
##            1670       1671       1672       1673       1674       1675
##            1676       1677       1688       1689       1690       1691
##            1692       1693       1694       1695       1696       1697
##            1698       1699       1700       1701       1702       1703
##            1704       1705       1706       1707       1708       1709
##            1710       1711       1712       1713       1714       1715
##            1716       1717       1718       1719       1720       1721
##            1722       1723       1724       1725       1726       1727
##            1728       1740       1741       1742       1743       1744
##            1745       1746       1747       1748       1749       1750
##            1751       1752       1753       1754       1755       1756
##            1757       1758       1759       1760       1761       1762
##            1763       1764       1765       1766       1767       1768
##            1769       1770       1771       1772       1773       1774
##            1775       1776       1777       1778       1779       1792
##            1793       1794       1795       1796       1797       1798
##            1799       1800       1801       1802       1803       1804
##            1805       1806       1807       1808       1809       1810
##            1811       1812       1813       1814       1815       1816
##            1817       1818       1819       1820       1821       1822
##            1823       1824       1825       1826       1827       1828
##            1829       1830       1844       1845       1846       1847
##            1848       1849       1850       1851       1852       1853
##            1854       1855       1856       1857       1858       1859
##            1860       1861       1862       1863       1864       1865
##            1866       1867       1868       1869       1870       1871
##            1872       1873       1874       1875       1876       1877
##            1878       1879       1880       1881       1882       1894
##            1895       1896       1897       1898       1899       1900
##            1901       1902       1903       1904       1905       1906
##            1907       1908       1909       1910       1911       1912
##            1913       1914       1915       1916       1917       1918
##            1919       1920       1921       1922       1923       1924
##            1925       1926       1927       1928       1929       1930
##            1931       1932       1946       1947       1948       1949
##            1950       1951       1952       1953       1954       1955
##            1956       1957       1958       1959       1960       1961
##            1962       1963       1964       1965       1966       1967
##            1968       1969       1970       1971       1972       1973
##            1974       1975       1976       1977       1978       1979
##            1980       1981       1982       1983       1997       1998
##            1999       2000       2001       2002       2003       2004
##            2005       2006       2007       2008       2009       2010
##            2011       2012       2013       2014       2015       2016
##            2017       2018       2019       2020       2021       2022
##            2023       2024       2025       2026       2029       2030
##            2031       2032       2049       2050       2051       2052
##            2053       2054       2055       2056       2057       2058
##            2059       2060       2061       2062       2063       2064
##            2065       2066       2067       2068       2069       2070
##            2071       2072       2073       2074       2075       2076
##            2077       2082       2083       2100       2101       2102
##            2103       2104       2105       2106       2107       2108
##            2109       2110       2111       2112       2113       2114
##            2115       2116       2117       2118       2119       2120
##            2121       2122       2123       2124       2125       2126
##            2127       2152       2153       2154       2155       2156
##            2157       2158       2159       2160       2161       2162
##            2163       2164       2165       2168       2169       2170
##            2171       2172       2173       2174       2175       2176
##            2177       2178       2179       2203       2204       2205
##            2206       2207       2208       2209       2210       2211
##            2212       2213       2214       2215       2216       2217
##            2218       2219       2220       2221       2222       2223
##            2224       2225       2226       2227       2228       2229
##            2230       2256       2257       2258       2259       2260
##            2261       2262       2263       2264       2265       2266
##            2267       2268       2269       2270       2271       2272
##            2273       2274       2275       2276       2277       2278
##            2279       2280       2281       2282       2307       2308
##            2309       2310       2311       2312       2313       2314
##            2315       2316       2317       2318       2319       2320
##            2321       2322       2323       2324       2325       2326
##            2327       2328       2329       2330       2331       2332
##            2362       2365       2366       2367       2368       2369
##            2370       2371       2372       2373       2374       2375
##            2376       2377       2378       2379       2380       2381
##            2382       2383       2417       2418       2419       2420
##            2421       2422       2423       2424       2425       2426
##            2427       2428       2429       2430       2431       2432
##            2433       2434       2469       2470       2471       2472
##            2473       2474       2475       2476       2477       2478
##            2479       2480       2481       2482       2483       2484
##            2485       2520       2521       2522       2523       2524
##            2525       2526       2527       2528       2529       2530
##            2531       2532       2533       2534       2535       2536
##            2571       2572       2573       2574       2575       2576
##            2577       2578       2579       2580       2581       2582
##            2583       2584       2585       2586       2587       2588
##            2589       2590       2622       2623       2624       2625
##            2626       2627       2628       2629       2630       2631
##            2632       2633       2634       2635       2636       2637
##            2638       2639       2640       2677       2678       2679
##            2680       2681       2682       2683       2684       2685
##            2686       2687       2688       2689       2690       2691
##            2692       2729       2730       2731       2732       2733
##            2734       2735       2736       2737       2738       2739
##            2740       2742       2781       2782       2783       2784
##            2785       2786       2787       2788       2789       2790
##            2833       2834       2835       2836       2837       2838
##            2839       2840       2886       2887       2888       2889
##            2890       2938       2939       2940       2991
##  [ reached getOption("max.print") -- omitted 1811 rows ]
## 
## dst_r_> delta<-0.01
## 
## dst_r_> h<-neighbourhoods(U,delta)
## 
## dst_r_> U$hexagon<-paste0(h@cID)
## 
## dst_r_> hD<-dist_areas_f(U,h)
```

```
## Error in 2 * delta: non-numeric argument to binary operator
```


```r
example(newdist)
```

```
## 
## newdst> delta<-.005
## 
## newdst> sicks<-(1:nrow(U))[y=="sick"]
```

```
## Error in eval(ei, envir): object 'y' not found
```


```r
example(updatedist)
```

```
## 
## updtds> delta<-.005
## 
## updtds> sicks<-(1:nrow(U))[y=="sick"]
```

```
## Error in eval(ei, envir): object 'y' not found
```
`projpointonseg_a` is a function that gives the position
of the closest point of a segment to a specific point on the plane. The position is the ratio of the distance of the closest point to one extremity of the segment divided by the length of the segment.

```r
example(projpointonseg_a,echo=F)
```

<img src="figure/unnamed-chunk-12-1.png" title="plot of chunk unnamed-chunk-12" alt="plot of chunk unnamed-chunk-12" width="100%" /><img src="figure/unnamed-chunk-12-2.png" title="plot of chunk unnamed-chunk-12" alt="plot of chunk unnamed-chunk-12" width="100%" />
`projpointonseg` is a function of a point and a segment that returns the closest point to the input point in the segment.

```r
example(projpointonseg,echo=F)
```

<img src="figure/unnamed-chunk-13-1.png" title="plot of chunk unnamed-chunk-13" alt="plot of chunk unnamed-chunk-13" width="100%" />
`distpointtoseg` is a function that returns the distance between a point and a segment.

```r
example(distpointtoseg,echo=F)
```

<img src="figure/unnamed-chunk-14-1.png" title="plot of chunk unnamed-chunk-14" alt="plot of chunk unnamed-chunk-14" width="100%" />


```r
example(segment.intersect,echo=F)
```

<img src="figure/unnamed-chunk-15-1.png" title="plot of chunk unnamed-chunk-15" alt="plot of chunk unnamed-chunk-15" width="100%" />

```r
example(distsegmenttosegment,echo=F)
```

<img src="figure/unnamed-chunk-16-1.png" title="plot of chunk unnamed-chunk-16" alt="plot of chunk unnamed-chunk-16" width="100%" />


```r
example(distsegmenttopoly)
```

```
## 
## dstsgm> data(Avo_fields,package="Strategy")
## 
## dstsgm> polygon1<-Avo_fields[1,]
## 
## dstsgm> A<-polygon1@polygons
## 
## dstsgm> B<-A[[1]]@Polygons[[1]]@coords
## 
## dstsgm> s<-cbind(runif(2,min = min(B[,1]),max=max(B[,1])),runif(2,min=min(B[,2]),max=max(B[,2])))
## 
## dstsgm> plot(B,type='l')
```

```
## 
## dstsgm> s1<-s
## 
## dstsgm> points(s,type="l",lwd=4,col="green")
## 
## dstsgm> x<-vector()
## 
## dstsgm> for(i in 1:(nrow(B)-1)){
## dstsgm+ s2<-B[(i:(i+1)),]
## dstsgm+ dd<-distsegmenttosegment(s1,s2)
## dstsgm+ l<-which(c(distpointtoseg(s1[1,],s2),distpointtoseg(s1[2,],s2),distpointtoseg(s2[1,],s1),distpointtoseg(s2[2,],s1))==dd)[1]
## dstsgm+ min(as.matrix(dist(rbind(s1,s2),diag=T,upper = T))[3:4,1:2]);dd
## dstsgm+ sk<-if(l<=2){s2}else{s1}
## dstsgm+ p=rbind(s1,s2)[l,]
## dstsgm+ points(x=p[1],y=p[2] ,col="red",cex=2)
## dstsgm+ points(projpointonseg(p,sk)[1],projpointonseg(p,sk)[2],col="red",cex=2)
## dstsgm+ segments(x0 = p[1],y0=p[2],x1=projpointonseg(p,sk)[1],y1=projpointonseg(p,sk)[2])
## dstsgm+ x=c(x,dd)
## dstsgm+ }
## 
## dstsgm> min(x)
## [1] 5.214712e-05
## 
## dstsgm> distsegmenttopoly(s,B)
## [1] 5.214712e-05
## 
## dstsgm> distpolytopoly()
```

```
## Error in nrow(poly1): argument "poly1" is missing, with no default
```

<img src="figure/unnamed-chunk-17-1.png" title="plot of chunk unnamed-chunk-17" alt="plot of chunk unnamed-chunk-17" width="100%" />


```r
example(distpolytopoly,ask=F)
```


```r
example(polydistmat)
```


```r
example(connectedpop)
```



