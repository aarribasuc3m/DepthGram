#### Examples

library(fda.usc)
library(fda)
library(roahd)
source("depthGram.R")
source("depthGramPlot.R")

#### 1. Canadian weather data set (from fda):
# Bi-variate functional data set with daily temperature and precipitation at 
# 35 different locations in Canada averaged over 1960 to 1994
Data<-list()
Data[[1]]<-t(CanadianWeather$dailyAv[,,1])
Data[[2]]<-t(CanadianWeather$dailyAv[,,2])
names=row.names(Data[[1]])

##DepthGram
DG<-depthGram(Data,marg.out=T,ids=names)
DGplot<-depthGramPlot(DG,limits=T,ids=names,plotly=T,plot.title="Canadian weather data set")
colDG<-DGplot$color

##Plotting the data to visualize outliers

#Sort by average temperature to define gray scale
temp.av<- apply(Data[[1]],1,mean)
sorted.loc<-setdiff(sort(temp.av, index.return=T)$ix,DGplot$out)
grays=gray.colors(length(sorted.loc), start = 0, end = 0.9, gamma = 1.5, alpha = NULL)
colDG[sorted.loc]<-grays
ids=c(sorted.loc,DGplot$out)
Data=lapply(Data,function(x){x[ids,]})
colDG=colDG[ids]

plot( mfData( 1:365, Data ),
      xlab = 'Date', ylab = list( 'C', 'mm' ),
      main = list( 'Temperature', 'Precipitation' ) , col=colDG)


#### 2. Aemet data set (from fda.usc):
# 3-variate functional data set with daily temperature, log-precipitation and 
# wind speed at 73 Spanish weather stations averaged over 1980 to 2009
data(aemet)

Data<-list()
Data[[1]]<-aemet[[2]]$data
Data[[2]]<-aemet[[3]]$data
Data[[3]]<-aemet[[4]]$data
names=row.names(Data[[1]])

DG<-depthGram(Data,marg.out=T,ids=names)
DGplot<-depthGramPlot(DG,limits=T,ids=names,plotly=F,plot.title="Aemet data set")
colDG<-DGplot$color

#Sort by average temperature to define gray scale
temp.av<- apply(Data[[1]],1,mean)
sorted.loc<-setdiff(sort(temp.av, index.return=T)$ix,DGplot$out)
grays=gray.colors(length(sorted.loc), start = 0, end = 0.9, gamma = 1.5, alpha = NULL)
colDG[sorted.loc]<-grays
ids=c(sorted.loc,DGplot$out)
Data=lapply(Data,function(x){x[ids,]})
colDG=colDG[ids]

plot( mfData( 1:365, Data),
      xlab = 'Date', ylab = list( 'C', 'mm', 'm/s^2' ),
      main = list( 'Temperature', 'Log - Precipitation', 'Wind speed' ) , col=colDG)


#### 3. 8-Lead ECG trace of healthy subjects (from roahd):
# 8-variate functional data set containing the 8-Lead ECG traces of 50 healthy subjects 

Data<-lapply(mfD_healthy$fDList, function(x){x$values})
names=row.names(Data[[1]])

DG<-depthGram(Data,marg.out=T,ids=names)
DGplot<-depthGramPlot(DG,limits=T,ids=names,plotly=T,plot.title="EGG healthy data set")
colDG<-DGplot$color

#Sort by average 1st lead to define gray scale
first.av<- apply(Data[[1]],1,mean)
sorted.loc<-setdiff(sort(first.av, index.return=T)$ix,DGplot$out)
grays=gray.colors(length(sorted.loc), start = 0, end = 0.9, gamma = 1.5, alpha = NULL)
colDG[sorted.loc]<-grays
ids=c(sorted.loc,DGplot$out)
Data=lapply(Data,function(x){x[ids,]})
colDG=colDG[ids]

plot( mfData( 1:1024, Data),
      xlab = 'Time (ms)', ylab='Potential (mV)',
      main = list( 'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'D1', 'D2' ) , col=colDG)

#### 4. 8-Lead ECG trace of subjects suffering from Left-Bundle-Branch-Block (from roahd):
# 8-variate functional data set containing the 8-Lead ECG traces of LBBB patients

Data<-lapply(mfD_LBBB$fDList, function(x){x$values})
names=row.names(Data[[1]])

DG<-depthGram(Data,marg.out=T,ids=names)
DGplot<-depthGramPlot(DG,limits=T,ids=names,plotly=T,plot.title="EGG LBBB data set")
colDG<-DGplot$color

#Sort by average 1st lead to define gray scale
first.av<- apply(Data[[1]],1,mean)
sorted.loc<-setdiff(sort(first.av, index.return=T)$ix,DGplot$out)
grays=gray.colors(length(sorted.loc), start = 0, end = 0.9, gamma = 1.5, alpha = NULL)
colDG[sorted.loc]<-grays
ids=c(sorted.loc,DGplot$out)
Data=lapply(Data,function(x){x[ids,]})
colDG=colDG[ids]

plot( mfData( 1:1024, Data),
      xlab = 'Time (ms)', ylab='Potential (mV)',
      main = list( 'V1', 'V2', 'V3', 'V4', 'V5', 'V6', 'D1', 'D2' ) , col=colDG)

