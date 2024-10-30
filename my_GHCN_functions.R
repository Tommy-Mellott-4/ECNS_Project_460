###############################################################################
# rm(list=ls()); graphics.off(); gc()
# load('WorkSpace.RData')
###############################################################################
# my_require(c("doBy","chron","readr","data.table")); setDTthreads(1)
###############################################################################

###############################################################################
# for debugging only
# datadir=NULL;stationsdir=NULL;newfile=T
###############################################################################
my_ghcnd_stations=function(datadir=NULL,stationsdir=NULL,newfile=F){
###############################################################################
 require(readr)
# check if datadir has been passed, if the directory exists, and create
# datadir directory if it does not exist
 if(is.null(datadir)){
  (datadir = paste(getwd(),'/data/',sep=''))  # create name for data directory
 } # end if(is.null(datadir))

 if(is.null(stationsdir)){
  (stationsdir = paste(datadir,'station_info/',sep=''))  # create name for data directory
 } # end if(is.null(datadir))

 if(!dir.exists(datadir)) dir.create(datadir)
 if(!dir.exists(stationsdir)) dir.create(stationsdir)
###############################################################################


###############################################################################
# check if station data exists and download/process it
###############################################################################
# create name of desired stations data file
 (fname=paste(stationsdir,'ghcnd-stations.RData',sep=''))
# if file does not exist or user has instructed to download new files
# download and process the station files
 if(!file.exists(fname)|newfile==T){
  # create list of files to download
  filelist=c('ghcnd-inventory.txt','ghcnd-stations.txt')
#loop through list, download files in data directory
  (file=filelist[1])
  for(file in filelist){
   (fname1=paste('https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/',file,sep=''))
   (fname2=paste(datadir,file,sep=''))
   download.file(fname1,fname2)
  } # end loop for(file in filelist)
###############################################################################


###############################################################################
# Function for converting data frame fields
###############################################################################
  my_convert=function(df,ctypes){
   for(j in 1:length(ctypes)) {
    if(ctypes[j]=='c') df[,j] = as.character(df[,j])
    if(ctypes[j]=='n') df[,j] = as.numeric(as.character(df[,j]))
    if(ctypes[j]=='i') df[,j] = as.integer(as.character(df[,j]))
   } # end for(j in 1:length(ctypes))
   df
  } # end function my_convert
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

  (fname1=paste(datadir,'ghcnd-stations.txt',sep=''))

  stations1=as.data.frame(read_fwf(file=fname1,
          fwf_positions(cstart1,cend1,cnames1),trim_ws = F))
##
  str(stations1)
# Note numeric fields have been read as character
# Convert to numeric
  stations=my_convert(stations1,ctypes1)
  str(stations)
  
  (file=paste(stationsdir,'ghcnd-stations.csv',sep=''))
  write.csv(stations,file=file,row.names=F)
  (file=paste(stationsdir,'ghcnd-stations.RData',sep=''))
  save(stations,file=file)
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
# Merge stations and stations2 data
# check for unintended matches in object names
  names(stations)
  names(stations2)
##
  stations = merge(stations,stations2)
  head(stations)
  str(stations)
  dim(stations)
  summary(stations)
  (file=paste(stationsdir,'stations_inventory.RData',sep=''))
  save(stations,file=file)
###############################################################################
 } # end if(!file.exists(fname)|newfile==T)
###############################################################################
} # end function my_ghcnd_stations
###############################################################################


###############################################################################
# for debugging my_ghcnd
#  stationid='USC00241044';datadir=NULL; downdir=NULL;RDatadir=NULL
#  newfile=F;METADATA=NULL
#  bulk=T;bulkdir=NULL
###############################################################################
my_ghcnd=function(stationid='USC00241044',datadir=NULL,downdir=NULL,
	RDatadir=NULL,bulkdir=NULL,bulk=T,newfile=F,METADATA=NULL){
#		'D:/my_weather/data/ghcnd_all/'){
###############################################################################
 require(chron); require(readr)
###############################################################################
# Function for converting data frame fields
###############################################################################
  my_convert=function(df,ctypes){
   for(j in 1:length(ctypes)) {
    if(ctypes[j]=='c') df[,j] = as.character(df[,j])
    if(ctypes[j]=='n') df[,j] = as.numeric(as.character(df[,j]))
    if(ctypes[j]=='i') df[,j] = as.integer(as.character(df[,j]))
   } # end for(j in 1:length(ctypes))
   df
  } # end function my_convert
###############################################################################

###############################################################################
if(is.null(datadir)){
  (datadir = paste(getwd(),'/data/',sep=''))  # create name for data directory
 } # end if(is.null(datadir))

 if(is.null(downdir)){
  (downdir = paste(datadir,'downloads/',sep=''))  
 } # end if(is.null(downdir))
  
 if(is.null(RDatadir)){
  (RDatadir = paste(datadir,'Rdata/',sep=''))  
 } # end if(is.null(RDatadir))

 if(is.null(bulkdir)){
  (bulkdir = paste(datadir,'ghcnd_all/',sep=''))  
 } # end if(is.null(bulkdir))
   
 if(!dir.exists(datadir)) dir.create(datadir)
 if(!dir.exists(downdir)) dir.create(downdir)
 if(!dir.exists(RDatadir)) dir.create(Rdatadir)
 if(!dir.exists(bulkdir)) dir.create(bulkdir)

# check to see if we have downloaded station data and stored in RDatadir folder
 (fname = paste(RDatadir,stationid,'.RData',sep=''))

 # if file exists and newfile!=F  load it
 if(file.exists(fname)&newfile!=T) load(fname)
 	

 # if file does not exist or newfile==T download data from ghcn website and store data
 if(!file.exists(fname)|newfile==T){
  # create metadata for download time and other information
  (date_dwn = Sys.Date())
  (jday_dwn = as.numeric(julian(date_dwn)))
  names1=c('id','year','month','element')
  names2=paste(rep(c('VALUE','MFLAG','QFLAG','SFLAG'),31),rep(1:31,each=4),sep='')
  cnames=c(names1,names2)
  widths=c(11,4,2,4,rep(c(5,1,1,1),31))
  ctypes=c('c','i','i','c',rep(c('i','c','c','c'),31))

  filetest=F
  if(bulk==T){
   if(dir.exists(bulkdir)) {  
    (fnameD=paste(bulkdir,stationid,'.dly',sep=''))
	if(file.exists(fnameD)) filetest=T
   } # end if on bulkdrive
  } # end if(bulk==T)

  if(bulk!=T|filetest==F|newfile==T){ 
   (fnameD=paste(
   'https://www1.ncdc.noaa.gov/pub/data/ghcn/daily/all/',stationid,'.dly',sep=''))
   print(paste('download',stationid))  
   }

   
   WDATA=as.data.frame(
     read_fwf(fnameD,fwf_widths(widths,cnames),trim_ws=F,show_col_types = F))

   WDATA=my_convert(WDATA,ctypes)
   WDATA[WDATA<(-9998)] = NA
   save(WDATA,date_dwn,jday_dwn,METADATA,file=fname)
  } # end if(!file.exists(fname)|newfile==T))
###############################################################################
 WDATA
} #end function my_ghcnd
###############################################################################




###############################################################################
# debugon
#  stationid='USC00241044';elements=c('TMIN','TMAX')
#  datadir=NULL;downdir=NULL;RDatadir=NULL;bulkdir=NULL;newfile=F;METADATA=NULL
#  bulk=T;bulkdir=NULL
###############################################################################

###############################################################################
ghcnd_elements=function(stationid='USC00241044',elements=c('TMIN','TMAX'),
  datadir=NULL,downdir=NULL,RDatadir=NULL,bulkdir=NULL,bulk=T,newfile=F,
	METADATA=NULL,elementdir=tempdir){
###############################################################################
 pkglist=c('chron','doBy','data.table')
 for (p in pkglist) {
  if(!require(p,character.only = TRUE)) {
    install.packages(p)
    require(p,character.only = TRUE)}
 } # end loop for (p in pkglist)
 setDTthreads(1)
#######################################
if(is.null(datadir)){
  (datadir = paste(getwd(),'/data/',sep=''))  # create name for data directory
 } # end if(is.null(datadir))

 if(is.null(downdir)){
  (downdir = paste(datadir,'downloads/',sep=''))  
 } # end if(is.null(downdir))
  
 if(is.null(RDatadir)){
  (RDatadir = paste(datadir,'Rdata/',sep=''))  
 } # end if(is.null(RDatadir))

 if(is.null(bulkdir)){
  (bulkdir = paste(datadir,'ghcnd_all/',sep=''))  
 } # end if(is.null(bulkdir))
   
 if(is.null(elementdir)){
  (elementdir = paste(datadir,'temp_data/',sep=''))  
 } # end if(is.null(elementdir))

 
 if(!dir.exists(datadir)) dir.create(datadir)
 if(!dir.exists(downdir)) dir.create(downdir)
 if(!dir.exists(RDatadir)) dir.create(Rdatadir)
 if(!dir.exists(bulkdir)) dir.create(bulkdir)
 if(!dir.exists(elementdir)) dir.create(elementdir)

 
 wdata=my_ghcnd(stationid=stationid,datadir=datadir,downdir=downdir,
	RDatadir=RDatadir,bulkdir=bulkdir,bulk=bulk,newfile=newfile,
 	METADATA=METADATA)

 ELEM_LIST=list(NULL)

 ie=1
 for (ie in 1:length(elements)){
  (element_ie=elements[ie])
  (ename=tolower(element_ie))
  ################################################################################
   # subset  and process the TEMP data
   dftmp = subset(wdata,element==element_ie)
   ###############################################################################
   # Stack up data
   # NOTE: The package data.table has nicer features for "stacking" this data
   # I am demonstrating how we can do this using base R functionality
   ###############################################################################
   # Create empty or NULL list object to store results
   ELEMJ = list(NULL)
   ##############################################################################
   # Stack data - pull and stack data for days 1 - 31
   # We will discuss details in class
   j=1
   for (j in 1:31) {
    (txt_j = paste(c('VALUE','MFLAG','QFLAG','SFLAG'),j,sep=''))
    (cnames=c('id','year','month',txt_j))

    tmp=dftmp[,cnames]; tmp$day = j
    head(tmp)
    (names(tmp)=c('id','yr','mo',paste(rep(ename,4),
                c('','_mflag','_qflag','_sflag'),sep=''),'day'))

    tmp$jday=julian(tmp$mo,tmp$day,tmp$yr)
    head(tmp);tail(tmp)
    # eliminate missing data
    tmp=tmp[!is.na(tmp[,4]),]
    head(tmp);tail(tmp)

    # store results from day j in dynamically created j'th cell of list object
    # TMINMAX
    ELEMJ[[j]] = tmp
   } # end loop for (j in 1:31)
   ############################################################################
   # Use data.table function to "rbind" the data in the TMINMAX list
   elemj = rbindlist(ELEMJ)
   # Examine resulting data frame
   head(elemj)
   tail(elemj)

   tmp=elemj[,c('id','jday')]
   elemj=merge(tmp,elemj)

   head(elemj)
   tail(elemj)

   dim(elemj)
  #############################################################################
  ELEM_LIST[[ie]]=elemj
 } # end loop #for (ie in 1:length(elements))
###############################################################################
  head(ELEM_LIST[[1]])
  if(length(ELEM_LIST)>1) head(ELEM_LIST[[2]])
  wdata = ELEM_LIST[[1]]
  if(length(ELEM_LIST)>1) {
    for(i in 2:length(ELEM_LIST)) {
     wdata = merge(wdata,ELEM_LIST[[i]],by=c('id','jday','yr','mo','day'),all=T)
    } # end loop for(i in 2:length(ELEM_LIST))
  } # end if(length(ELEM_LIST>1))

  head(wdata)
  tail(wdata)

  (date_dwn = Sys.Date())
  (jday_dwn = as.numeric(julian(date_dwn)))

  (fname=paste(elementdir,stationid,'.RData',sep=''))
  
  save(wdata,date_dwn,jday_dwn,METADATA,file=fname)
  wdata
###############################################################################
} #end function ghcnd_elements
###############################################################################






###############################################################################
#tests
###############################################################################
#my_ghcnd_stations()
###############################################################################

###############################################################################
# (stationid='USC00250070')
#
# system.time({
#  wdata=my_ghcnd(stationid,dirname=NULL,newfile=F)
# })
#
# system.time({
#  wdata=my_ghcnd(stationid,dirname=NULL,newfile=T)
# })
#
# system.time({
#  wdata=my_ghcnd(stationid,dirname='downloads',newfile=T)
# })
#
#
# system.time({
#  wdata=my_ghcnd(stationid,dirname=paste(getwd(),'/dumpdata/',sep=''),newfile=T)
# })
################################################################################
# wdata=ghcnd_elements(stationid='USC00241044',newfile=F)
################################################################################
#(stationdir = paste(getwd(),'/prcp_data/',sep=''))
#
#wdata=ghcnd_elements(stationid='USC00241044',elements=c('PRCP'),
#      downdir=NULL,stationdir=stationdir,newfile=F)
###############################################################################
#(stationid='USC00250070'); elements=c('PRCP')
# downdir=NULL; stationdir=stationdir; newfile=F
#bulk=T;bulkdrive='D:/Atwood/GHCN/daily/ghcnd_all/USCADATA/DATA0/'
  
