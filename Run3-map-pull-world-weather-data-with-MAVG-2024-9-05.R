###############################################################################
rm(list=ls()); graphics.off(); gc()
###############################################################################
plotmonth=F
###############################################################################
source('my_weather_setup.R')
keep=ls()
###############################################################################
my_require(c("doBy","chron","sp"))
my_require(c("readr","data.table")); setDTthreads(1)
###############################################################################

###############################################################################
# load Maps.Rdata and keep maps 'WMaps','USMapsS48','USMapsC48'
(fname=paste(datadir,'Maps.RData',sep='')); load(fname)
keep=c(keep,'WMaps','USMapsS48','USMapsC48','USMapsS50','USCAMaps','USCAMXMaps')
###############################################################################
# Shorter map names
Smaps = USMapsS48; Cmaps = USMapsC48; Smaps50=USMapsS50
###############################################################################
# Set color palette
palette('R3')
###############################################################################


###############################################################################
# load long-term TEMP world weather station metadata
###############################################################################
(fname=paste(stationsdir,'WORLD-TEMP-Stations-1900-2023.RData',sep=''))
load(fname)
head(WORLD_stations)
###############################################################################
# Subset station fields, eliminate duplicate records and store results in
# data frame with shorter name
df0=unique(WORLD_stations[,c('id','long','lat','elevation','name')])
head(df0)
dim(df0)
###############################################################################




###############################################################################
# download some data
###############################################################################
# check which stations' .Rdata files we have already downloaded
###############################################################################
# Move to downloads sub-directory and Find files in downdir which have '.Rdata'
# in file name
###############################################################################
files = list.files(path=tempdir,pattern='.RData')
#files = list.files(path=RDatadir,pattern='.RData')
# Examine name of first file found
files[1]
# Remove text '.RData' from file list and keep station ids. Examine first id
# found
st_ids = gsub('.RData','',files)
st_ids[1]
# Assign colors for plotting station locations with default color black (1)
# and a location id in the downloaded station list st_ids and return to
# working directory
df0$color = ifelse(is.na(match(df0$id,st_ids)),4,2)
###############################################################################



###############################################################################
my_x11()
plot(WMaps)
points(df0$long,df0$lat,pch=20,cex=1.5,col=4)
plot(USCAMXMaps,add=T,lwd=1)

# select  bounding points and plot new map
mtext('click on two points for zoomed map',cex=2)
tmp=locator(2)

plot(WMaps,xlim=sort(tmp$x),ylim=sort(tmp$y))
points(df0$long,df0$lat,pch=20,cex=2,col=df0$color)
plot(WMaps,add=T,lwd=2)
plot(USCAMXMaps,add=T,lwd=3)
mtext('red points have already been downloaded',line=-1,cex=1.5,font=2)
# select a new point and then pull data for point(s) near that point
mtext('click on one point',cex=2)
tmp=locator(1)

points(tmp$x,tmp$y,pch=20,cex=4,col=3)
df0$dist=distcalc(tmp$x,tmp$y,df0$long,df0$lat)
df0=orderBy(~dist,data=df0)
(df00=df0[1,])
###############################################################################
# # Note for homework  9/29/2021
# # closest ten stations
# (df1=df0[1:10,])
# text(df1$long,df1$lat-0.15,paste(1:10),font=2)
# 
# #stations within 150 miles
# (df2=subset(df0,dist<=150))
# text(df2$long,df2$lat-0.15,paste(1:nrow(df2)),font=2,col=2)
# nrow(df2)
###############################################################################






###############################################################################
# pull stationid from station closest to selected point
(stationid=df0$id[1])
(station_name=df0$name[1])
# check to see if we have downloaded station data and stored in downloads folder
(fname = paste(tempdir,stationid,'.RData',sep=''))
#(fname = paste(RDatadir,stationid,'.RData',sep=''))

# example of using the file.exists function
file.exists(fname)
!file.exists(fname)

# if file exists load it
if(file.exists(fname)) load(fname)

# if file does not exist download data from ghcn website and store data
if(!file.exists(fname)){
 (dwndate = Sys.Date())
 names1=c('id','year','month','element')
 names2=paste(rep(c('VALUE','MFLAG','QFLAG','SFLAG'),31),rep(1:31,each=4),sep='')
 cnames=c(names1,names2)
 widths=c(11,4,2,4,rep(c(5,1,1,1),31))
 ctypes=c('c','i','i','c',rep(c('i','c','c','c'),31))

 wdata=my_ghcnd(stationid=stationid)
 str(wdata)
 wdata[wdata<(-9998)] = NA
} # end else option
################################################################################

################################################################################


################################################################################
# subset  and process the TEMP data
dfTMP = subset(wdata,element=='TMIN'|element=="TMAX")
dim(dfTMP)
head(dfTMP,2)
###############################################################################
# Stack up data
# NOTE: The package data.table has nicer features for "stacking" this data
# I am demonstrating how we can do this using base R functionality
###############################################################################
# Create empty or NULL list object to store results
TMINMAX = list(NULL)
# Subset the dfTMP data frame into tmin and tmax data frames and keep only
# data for years 1900 and later
df_tmin=subset(dfTMP,element=='TMIN'&year>=1900)
df_tmax=subset(dfTMP,element=='TMAX'&year>=1900)
# Examine first rows of tmin and tmax data set
head(df_tmin,1)
head(df_tmax,1)

##############################################################################
# Stack data - pull and stack data for days 1 - 31
# We will discuss details in class
j=1
for (j in 1:31) {
 (txt_j = paste(c('VALUE','MFLAG','QFLAG','SFLAG'),j,sep=''))
 (cnames=c('id','year','month',txt_j))

 tmin=df_tmin[,cnames]; tmin$day = j
 head(tmin)
 names(tmin)=c('id','yr','mo','tmin','tmin_mflag','tmin_qflag','tmin_sflag','day')
 tmin$jday=julian(tmin$mo,tmin$day,tmin$yr)
 tmin=tmin[,c('id','jday','yr','mo','day','tmin','tmin_mflag','tmin_qflag','tmin_sflag')]
 head(tmin)

 tmax=df_tmax[,cnames]; tmax$day = j
 head(tmax)
 names(tmax)=c('id','yr','mo','tmax','tmax_mflag','tmax_qflag','tmax_sflag','day')
 tmax$jday=julian(tmax$mo,tmax$day,tmax$yr)
 tmax=tmax[,c('id','jday','yr','mo','day','tmax','tmax_mflag','tmax_qflag','tmax_sflag')]
 head(tmax)

 tminmax=merge(tmin,tmax,all.x=T,all.y=T)
 head(tminmax)
 # store results from day j in dynamically created j'th cell of list object
 # TMINMAX
 TMINMAX[[j]] = tminmax
} # end loop for (j in 1:31)
###############################################################################
# Use data.table function to "rbind" the data in the TMINMAX list
tminmax = rbindlist(TMINMAX)
# Be sure "stacked" data is in jday order
tminmax=orderBy(~id+jday,data=tminmax)
# Examine resulting data frame
head(tminmax)
tail(tminmax)

dim(tminmax)
###############################################################################
# drop observations with flagged quality issues with either tmin or tmax
dim(tminmax)
tminmax=subset(tminmax,tmin_qflag==' '&tmax_qflag==' ')
dim(tminmax)
# drop days without both tmin and tmax observations
tminmax=subset(tminmax,tmin!='NA'&tmax!='NA')
dim(tminmax)
###############################################################################

###############################################################################
# Save processed TMIN-TMAX results in the stationdata sub-directory 9/29/2021
(fname=paste(tempdir,stationid,'.RData',sep=''))
#wdata = tminmax
save(stationid,station_name,wdata,dwndate,df00,file=fname)
###############################################################################


###############################################################################
# Process and display results
# convert  tenths of degrees C to degrees fahrenheit
tminmax$tmin=(((tminmax$tmin/10) * (9/5))+32)
tminmax$tmax=(((tminmax$tmax/10) * (9/5))+32)


my_x11()
plot(tminmax$tmin,tminmax$tmax)
abline(0,1,lwd=3,col=2)
###############################################################################
tminmax$tavg=(tminmax$tmin+tminmax$tmax)/2
###############################################################################

my_x11()
#par(mfrow=c(2,1))
boxplot(tmax~yr,data=tminmax,col=2, main=paste('TMAX',station_name))
abline(h=median(tminmax$tmax,na.rm=T),lwd=3)
summary(lm(tmax~yr,data=tminmax))


my_x11()
boxplot(tmin~yr,data=tminmax,col=4, main=paste('TMIN',station_name))
abline(h=median(tminmax$tmin,na.rm=T),lwd=3)
summary(lm(tmin~yr,data=tminmax))


my_x11()
boxplot(tavg~yr,data=tminmax,col=3, main=paste('TAVG',station_name))
abline(h=median(tminmax$tavg,na.rm=T),lwd=3)
summary(lm(tavg~yr,data=tminmax))


df0[1,]
###############################################################################
# compute 10 year moving averages of annual average temps
tdata=summaryBy(tmin+tavg+tmax~yr,data=tminmax,FUN=mean,na.rm=T)
tdata2=summaryBy(tmin+tavg+tmax~yr,data=tminmax,FUN=length.na)

tdata=merge(tdata,tdata2)

head(tdata)

tdata$MA_tmin=apply(lagMat(tdata$tmin.mean,10),1,mean)
tdata$MA_tavg=apply(lagMat(tdata$tavg.mean,10),1,mean)
tdata$MA_tmax=apply(lagMat(tdata$tmax.mean,10),1,mean)

tdata=na.omit(tdata)

head(tdata)
tail(tdata)

(vlines=seq(1900,2040,10))


my_x11()
par(mfrow=c(4,1),mar=c(2,3,1,1))
plot(tdata$yr,tdata$tmin.length.na,type='l',lwd=2,
    main=paste('Data Observations by Year',station_name))
abline(v=vlines,lty=3)

   
plot(tdata$yr,tdata$MA_tmax,type='l',lwd=2,col=2,main='TMAX - 10 Year Moving AVG')
abline(v=vlines,lty=3)

plot(tdata$yr,tdata$MA_tavg,type='l',lwd=2,col=3,main='TAVG - 10 Year Moving AVG')
abline(v=vlines,lty=3)

plot(tdata$yr,tdata$MA_tmin,type='l',lwd=2,col=4,main='TMIN - 10 Year Moving AVG')
abline(v=vlines,lty=3)
###############################################################################






###############################################################################
if(plotmonth==T){
###############################################################################  
askum=function(data=c(1:12,0),
     txt="Type number of month to plot OR type 13  to stop"){
 select.list(data, title=txt)
} # end function askum
##############################################################################
moplot=7
loop = 0
while(moplot>=1&loop<=12){
# my_x11()
# plot(-1:13,-1:13,type='n')
# text(0:12,0:12,paste(0:12),cex=2)
# (txt=c('stop',rep('mo',12)))
# text(1:13,0:12,txt,cex=1)
# mtext('click on month to plot or zero to stop',cex=2)
# (tmp=locator(1))
# moplot=round(tmp$x,0)
 moplot=askum()

 if(moplot>0) {

  tmp=subset(tminmax,mo==moplot)
   my_x11()
   boxplot(tmax~yr,data=tmp,col=2, main=paste('TMAX',station_name))
   mtext(paste('Month = ',month.name[moplot]),cex=1.5)
   abline(h=median(tmp$tmax,na.rm=T),lwd=3)
   summary(lm(tmax~yr,data=tmp))


   my_x11()
   boxplot(tmin~yr,data=tmp,col=4, main=paste('TMIN',station_name))
   mtext(paste('Month = ',month.name[moplot]),cex=1.5)
   abline(h=median(tmp$tmin,na.rm=T),lwd=3)
   summary(lm(tmin~yr,data=tmp))

 loop=loop+1
 }# end loop if(moplot>0)
} #end while
###############################################################################
}# end if(plotmonth==T){
###############################################################################
df00
###############################################################################
# summary(lm(tmin~yr,data=wdata))
# summary(MASS::rlm(tmin~yr,data=wdata))
# 
# summary(lm(tmax~yr,data=wdata))
# summary(MASS::rlm(tmax~yr,data=wdata))
# ###############################################################################
