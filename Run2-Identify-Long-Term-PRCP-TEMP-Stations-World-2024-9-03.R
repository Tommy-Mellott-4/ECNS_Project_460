###############################################################################
rm(list = ls())
load('WorkSpace.RData'); cleanup()
###############################################################################
my_require(c('sp'))
load(paste(datadir,'Maps.RData',sep=''))
###############################################################################

###############################################################################
# Load ghcn weather station metadata
###############################################################################
(fname=paste(stationsdir,'stations_inventory.RData',sep=''))
load(fname)
stations[1,]
###############################################################################



###############################################################################
(ELIST=rev(c('PRCP','TMIN')))
(YRMINS=rev(c(1900,1925,1950,1980)))
YRMAX=2024
###############################################################################
(ELEMENT=ELIST[1])
for(ELEMENT in ELIST){
 (YRMIN=YRMINS[1])
 for(YRMIN in YRMINS){	
  
# pull world long-term element stations and inventory
  df1=subset(stations,(first_year<=YRMIN)&last_year>=YRMAX
 	   &(element==ELEMENT))

  dim(df1)
  summary(df1)
# construct table with data counts of first year of data
  table(df1$first_year)
# construct table with data counts of last year of data
  table(df1$last_year)
###############################################################################
# How many stations in the WORLD have data since YRMIN or before
  (nWD=length(unique(df1$id)))
###############################################################################
# Make named copy of world long-term stations with element inventory
  WORLD_stations=df1
# Save WORLD_stations in compressed .RData file
  if(ELEMENT=='PRCP'){
   (fname3=paste(stationsdir,'WORLD-PRCP-Stations-',YRMIN,'-',YRMAX,'.RData',sep=''))
  }	
  if(ELEMENT=='TMIN'){
   (fname3=paste(stationsdir,'WORLD-TEMP-Stations-',YRMIN,'-',YRMAX,'.RData',sep=''))
  }	
  
 save(WORLD_stations,file=fname3)
###############################################################################
 } # end loop  for(YRMIN in YRMINS){	
}  # end loop for(ELEMENT in ELIST){
###############################################################################
