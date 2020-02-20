# Data: Annual Survey of Public Employment & Payroll

`dataONS` is an R package that contains data from the Annual Survey of Public Employment & Payroll.
The data was pulled from the [ONS website](https://www.ons.gov.uk).



## Installation
To install  the package, execute:

```r
devtools::install_github("DanielBonnery/dataONS")
```
Note that installation is slow, because part of the installation process is the downloading of data from the Census Bureau website.


##Documentation
To see the documentation, execute:


##How was the data pulled ?
To pull the data again, execute:

```r
alldata<-dataONS::get_data_from_web()
```

##How to use the data ?
To see an example of data use, execute:


```r
demo(Popbyparish,package = "dataONS")
```
