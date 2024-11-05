###############################################################################
# Load the required library
library(digest)
rm(list = ls())
load('WorkSpace.RData'); cleanup()
###############################################################################
my_require(c('sp'))
my_require('readr')
### get all stations from hcndir in directory
Each_Station <- list.files(path = hcndir, full.names = F)
head(Each_Station)
length(Each_Station)
## Pull all files in stationdir saved to harddrive
Station_Info_Files <- list.files(path = stationsdir,full.names = T)
### Separate the RData station info for the Temp and Precip Data
US_PRCP_Station_Info <- Station_Info_Files[4:8]
US_Temp_Station_Info <- Station_Info_Files[8:11]

### show all the paths to each saved RData info on Temp Station
US_Temp_Station_Info


### Load only the ones that has data since 1900
load(US_Temp_Station_Info[1])

dim(US_stations)
length(unique(US_stations$id))

##############################################################################
# create datestamp
(dwndate=Sys.Date())
###############################################################################
# Using information in readme.txt file (line 100) download .dly weather
# data from the  "https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/all/" directory
###############################################################################
(names1=c('id','year','month','element'))
(names2=paste(rep(c('VALUE','MFLAG','QFLAG','SFLAG'),31),rep(1:31,each=4),sep=''))
(cnames=c(names1,names2))
# specify column widths rather than beginning and ending columns for data fields
(widths=c(11,4,2,4,rep(c(5,1,1,1),31)))
#compare to last column in columns field in readme.txt
cumsum(widths)
# specify variable types
(ctypes=c('c','i','i','c',rep(c('i','c','c','c'),31)))


## Look at what you have already downloaded into the RData file
Rdata_downloaded <- list.files(paste0(getwd(),"/data/RData"),full.name = T)

### Removing the files that have been saved so do not redownload them
unlink(paste0(getwd(),"/data/RData/*"))
unlink(paste0(getwd(),"/data/downloads/*"))
## Pull 100 random ID's to keep so I do not have to download 800+ datasets
random_sample <- sample(1:length(US_stations$id),100, replace = F)

for(i in 1:100){
  # download data
  (fname1=paste0('https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/all/',US_stations$id[i] ,'.dly'))
  (fname2=paste(downdir,paste0(US_stations$id[i],'.dly'),sep=''))
  (fname3=paste(RDatadir,paste0(US_stations$id[i],'.RData'),sep=''))
  
  # time to download file and place in downloads directory
  system.time({
    download.file(fname1,fname2)
  })
  
  # time to download the file directly into R and save .RData file
  system.time({
    wdata=as.data.frame(read_fwf(fname1,fwf_widths(widths,cnames),trim_ws=F))
    wdata=my_convert(wdata,ctypes)
    wdata[wdata<(-9998)] = NA
    save(wdata,dwndate,file=fname3)
  })
}

