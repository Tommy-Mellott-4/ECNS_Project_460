###############################################################################
rm(list = ls())
load('WorkSpace.RData'); cleanup()
###############################################################################
my_require('readr')
###############################################################################
ls()

###############################################################################
# Get ready to read weather data directly from website
###############################################################################
# Function for converting data frame fields
# (Already loaded but printed here for students to see function)
###############################################################################
# my_convert=function(df,ctypes){
#   for(j in 1:length(ctypes)) {
#    if(ctypes[j]=='c') df[,j] = as.character(df[,j])
#    if(ctypes[j]=='n') df[,j] = as.numeric(as.character(df[,j]))
#    if(ctypes[j]=='i') df[,j] = as.integer(as.character(df[,j]))
#   }
#   df
# } # end function my_convert
###############################################################################

###############################################################################
# download readme and other information files from the website
# https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/
###############################################################################
filelist=c('readme.txt','ghcnd-countries.txt','ghcnd-inventory.txt',
     'ghcnd-states.txt','ghcnd-stations.txt')
###############################################################################
newfiles=T
###############################################################################
(file=filelist[1])
for(file in filelist){
 (fname1=paste('https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/',file,sep=''))
 (fname2=paste(datadir,file,sep=''))
 if(!file.exists(fname2)|newfiles==T) download.file(fname1,fname2)
} # end loop for(file in filelist)

# ghcnd_hcn data
file='ghcnd_hcn.tar.gz'
(fname1=paste('https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/',file,sep=''))
(fname2=paste(datadir,file,sep=''))
file.exists(fname2)
 
 if(!file.exists(fname2)) {
 	download.file(fname1,fname2)
   setwd(datadir)
   untar(fname2)
   setwd(workdir)
 } # end if(!file.exists(fname2))	
 

# ghcn_all data
#file='ghcnd_all.tar.gz'
#(fname1=paste('https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/',file,sep=''))
#(fname2=paste(datadir,file,sep=''))
#file.exists(fname2)

# if(!file.exists(fname2)) {
# 	download.file(fname1,fname2)
#   setwd(datadir)
#   (txt=paste('tar xzvf',file))
#   shell(txt)
# #  shell('tar xzvf ghcnd_all.tar.gz')
#   #untar(fname2)
#   setwd(workdir)
# } # end if(!file.exists(fname2))	

#getwd()	

###############################################################################
## Move to datadir, open readme.txt in default .txt editor and return to workdir
# setwd(datadir)
# shell('readme.txt',wait=F)
# setwd(workdir)
###############################################################################


###############################################################################
# Read and merge stations and stations inventory files to replicate output
# from rnoaa function ghcnd_stations()
###############################################################################
# ghcnd-stations data using file format in readme.txt file (line 406)
###############################################################################
(cnames1=c('id','latitude','longitude','elevation','state','name','gsn',
        'ncn_crn','wmo') )
cstart1 = c(1,13,22,32,39,42,73,77,81)
cend1 = c(11,20,30,37,40,71,75,79,85)
ctypes1 = c('c','n','n','n','c','c','c','c','c')
# ?read_fwf n=real, i=number, c = character field

(fname1=paste(datadir,'ghcnd-stations.txt',sep=''))

stations1=as.data.frame(read_fwf(file=fname1,
          fwf_positions(cstart1,cend1,cnames1),trim_ws = F))
##
str(stations1)
# Note numeric fields have been read as character
# Convert to numeric
stations1=my_convert(stations1,ctypes1)
str(stations1)
(fname2=paste(stationsdir,'ghcnd-stations.csv',sep=''))
write.csv(stations1,file=fname2,row.names=F)
(fname2=paste(stationsdir,'ghcnd-stations.RData',sep=''))
stations=stations1
save(stations,file=fname2)
rm(stations)
##
###############################################################################
# Station inventory data using file format in readme.txt file (line 535)
###############################################################################
wnames2 = c('id','lat','long','element','first_year','last_year')
cstart2 = c(1,13,22,32,37,42)
cend2 = c(11,20,30,35,40,45)
ctypes2 = c('c','n','n','c','i','i')

(fname2=paste(datadir,'ghcnd-inventory.txt',sep=''))

stations2=as.data.frame(read_fwf(file=fname2,
          fwf_positions(cstart2,cend2,wnames2),trim_ws = F))

str(stations2)
stations2=my_convert(stations2,ctypes2)
str(stations2)
##
# Merge stations1 and stations2 data
# check for unintended matches in object names
names(stations1)
names(stations2)
##
stations = merge(stations1,stations2)
head(stations)
str(stations)
dim(stations)
summary(stations)
(fname3=paste(stationsdir,'stations_inventory.RData',sep=''))
save(stations,file=fname3)
###############################################################################



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


# test time to load file from web
fname='https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/all/AE000041196.dly'

system.time({
df1=as.data.frame(read_fwf(fname,fwf_widths(widths,cnames),trim_ws=F))
})

# examine
str(df1[,1:10])
df2=my_convert(df1,ctypes)
str(df2[,1:10])

# -9999 indicates missing data
df2[df2<(-9998)] = NA
str(df2[,1:10])


df1[1:6,1:10]
df2[1:6,1:10]
###############################################################################



##############################################################################
# pull Bozeman weather station
# USC00241044
##############################################################################
# create datestamp
(dwndate=Sys.Date())

# download data
(fname1='https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/all/USC00241044.dly')
(fname2=paste(downdir,'USC00241044.dly',sep=''))
(fname3=paste(RDatadir,'USC00241044.RData',sep=''))

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
###############################################################################
# compare time to download file using rnoaa package
# create place to temporarily download data
#(dumpdir = paste(datadir,'dumpdata/',sep='')); dir.create(dumpdir)
#setwd(dumpdir)
# shell('del *.dly')
#setwd(workdir)

#my_require('rnoaa')
#ghcnd_cache$cache_path_set(full_path=dumpdir)

#time0 = my_seconds()
#  tmp = as.data.frame(ghcnd(stationid = 'USC00241044'))
#time1 = my_seconds()
# time to download the file directly into R and save .RData file
 wdata=as.data.frame(read_fwf(fname1,fwf_widths(widths,cnames),trim_ws=F))
 wdata=my_convert(wdata,ctypes)
 wdata[wdata<(-9998)] = NA
 save(wdata,dwndate,file=fname3)
time2 = my_seconds()

# rnoaa time
(time_rnoaa=time1-time0)
# our time
(time_our=time2-time1)
# compare
time_rnoaa/time_our
##############################################################################


##############################################################################
# Note: in the downloads subdirectory we see space savings of using .Rdata
# We now test the time savings of loading .RData rather than reading  and
# cleaning up wdata files
##############################################################################
time0=my_seconds()
 wdata1=as.data.frame(read_fwf(fname2,fwf_widths(widths,cnames),trim_ws=F))
 wdata1=my_convert(wdata1,ctypes)
 wdata1[wdata1<(-9998)] = NA
time1 = my_seconds()
 wdata2=load(fname3)
time2=my_seconds()

(time_read = time1-time0)
(time_load = time2-time1)

time_read/time_load
##############################################################################

