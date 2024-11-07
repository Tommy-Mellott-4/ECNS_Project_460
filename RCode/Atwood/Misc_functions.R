###############################################################################
# debugoff
PolyfindplotM_sf=function(XY,map,polynum=10,npoints=20,plotmap=F,plotmod=10,
plotnew=T,nsecs=0.0){
############################################################################### 
time11=my_seconds()
P_CLOSE=list(NULL)
# common for all points to be identified
# pull polygon boundary points for all polygons in map
Cpts=as.data.frame(st_coordinates(map))
Cpts$index=1:nrow(Cpts)
dim(Cpts)

# L3 = major polygon
# L2 = minor polygon with major polygon L3

# Specify endpoints on each polygon segment
Cpts$LX=my_lag(Cpts$X)
Cpts$LY=my_lag(Cpts$Y)
# Delete "replicated close the loop point in each polygon"  
Cpts$LL3=my_lag(Cpts$L3)
Cpts=subset(Cpts,!is.na(LX)&L3==LL3)
# compute the changes between x and y when we put extrapolated npoints between 
# the ends of each polygon segment
Cpts$dX=(Cpts$LX-Cpts$X)/npoints
Cpts$dY=(Cpts$LY-Cpts$Y)/npoints
# construct a matrix with xy coordinates for each point between end points 
# of a given polygon segment.  Row = polygon segment. Col = extrapolated point 
X=matrix(Cpts$X,nrow(Cpts),npoints)
Y=matrix(Cpts$Y,nrow(Cpts),npoints)
# construct extrapolated points
j=0
for(j in 0:(npoints-1)){
 X[,j+1] = Cpts$X +j*Cpts$dX
 Y[,j+1] = Cpts$Y +j*Cpts$dY 
} # end loop on j

longrad=X/57.2958
latrad=Y/57.2958
sin.latrad=sin(latrad)
cos.latrad=cos(latrad)

my_seconds()-time11
#########################################
# pull centroid of each map polygon
LL=as.data.frame(st_coordinates(suppressWarnings(st_centroid(map))))
names(LL)=c('XC','YC')
# index the centroids to the major polygons numbers in boundary Cpts
LL$L3=1:nrow(LL)
########################################
time12=my_seconds()
########################################
# find if points are in polygons
tmp1=st_as_sf(as.data.frame(XY),coords=c('x','y'))
tmp2=st_join(tmp1,map,join=st_within)
######################################




######################################
i=1
# debug off
for(i in 1:nrow(XY)){

# loop through for each point in original data set
# Note: Explore whether we could do this with parallel loop in the future       

(x0=XY[i,1]); (y0=XY[i,2])
 longrad0=x0/57.2958
 latrad0=y0/57.2958

 (polyin=tmp2$index[i])
#######################################
if(!is.na(polyin)){
 (pick0=which(LL$L3==polyin))
 LL2=LL[pick0,]
 dim(LL2)
} # end if(!is.na(polyin))

if(is.na(polyin)){ 
 # compute distance from given point to centroid of each map major polygon      
 LL$dist=distcalc(x0,y0,LL$XC,LL$YC)
 # sort in ascending distance
 LL=orderBy(~dist,data=LL)
 # pull identify closest polynum polygons
 LL2=LL[1:polynum,]
 dim(LL2)
} # end if(is.na(polyin))

 pick=which(Cpts$L3 %in% LL2$L3)
 Cpts2=Cpts[pick,]
 X2=X[pick,]
 Y2=Y[pick,]

longrad2=longrad[pick,]
latrad2=latrad[pick,]
sin.latrad2=sin.latrad[pick,]
cos.latrad2=cos.latrad[pick,]

cosD=(sin(latrad0)*sin.latrad2+cos(latrad0)*cos.latrad2*cos(longrad0-longrad2))
D=acos(cosD)
DIST=3963*D

#DIST2=matrix(distcalc(x0,y0,X2,Y2),nrow(X2),ncol(X2))
#
#summary(as.vector(DIST-DIST2))

dist=apply(DIST,1,min)
(dist2=min(dist)[1])
(pick1=which(dist==dist2)[1])
DIST[pick1,]
(pick2=which(DIST[pick1,]==dist2)[1])
(x2=X2[pick1,pick2])
(y2=Y2[pick1,pick2])
(polyclose=Cpts2$L3[pick1])
#(tmp=as.data.frame(cbind(i,x0,y0,polyclose,x2,y2,dist2)))
(tmp3=as.data.frame(cbind(x0,y0,polyin,polyclose,x2,y2,dist2)))
P_CLOSE[[i]]=tmp3
####################################### 


if(plotmap==T &(i==1|i%%plotmod==0|i==nrow(XY))){
 if(plotnew==T) my_x11()        
 map2=map[LL2$L3,]
 map2$color=3
 (pick3=which(LL2$L3==polyclose))
 map2$color[pick3]=2
 Cpts2=as.data.frame(st_coordinates(map2))
 dim(Cpts2)
 head(Cpts2)
 (xlim=range(c(x0,Cpts2$X)))
 (ylim=range(c(y0,Cpts2$Y)))

#points(x0,y0,pch=20,cex=4)
 plot(st_geometry(map2),lwd=1,col=map2$color,xlim=xlim,ylim=ylim,
         main=paste('point',i,'of',nrow(XY)))
 points(x0,y0,pch=20,cex=2)
 points(c(tmp3$x0,tmp3$x2),c(tmp3$y0,tmp3$y2),type='l')
 wait(nsecs)
 } # end if if(plotmap==T &(i==1|i%%plotmod==0|i==nrow(XY)))


} # end loop on i
#######################################
time13=my_seconds()

P_CLOSE=bind_rows(P_CLOSE)
###############################################################################
} # end function PolyfindplotM_sf
###############################################################################





################################################################################
# require(jcode64)
# X0=pts
# shapes=map0
# polynum=5
# npoints=10
# plotum=T
# pointbuffer=5
# nsecs=0.0
# dfunction=distcalc
# mytext=''
# newmap=F
################################################################################
# debugumoff
PolyfindplotM=function(X0,shapes,polynum=5,npoints=10,plotum=F,pointbuffer=5,
   nsecs=0.0,dfunction=distcalc,mytext='',newmap=F,plotmod=1){
################################################################################
 my_require('dplyr')
 shapes@data$index=1:nrow(shapes@data)
 spoly=shapes@polygons
 sdata=shapes@data

 # The number of polygons to check can't exceed number available
 (polynum=min(polynum,nrow(sdata)))

 # get centroid locations for each polygon and bind them to sdata
 coords=data.frame(coordinates(shapes))
 names(coords)=c('x','y')
 sdata=cbind(sdata,coords)

 # pull major centroid coordinates for shapefile slots
 longlat2=data.frame(coordinates(shapes))
 longlat2$slot1=1:nrow(longlat2)
 names(longlat2)=c('long2','lat2','slot1')
 longlat2=longlat2[,c('slot1','long2','lat2')]


 # Find the number of subpolygons in each major polygon
 # Count number of minor slots in  each major slot (polygon) in the shapefile
 sdata$npolys=sapply(spoly, function(xs) length(slot(xs,"Polygons")))

 #######################################################
 # create data-frame labels and pull extreme points
 slot1=rep(sdata$index,times=sdata$npolys)
 tmp1=data.frame(slot1)
 tmp1$index=1:nrow(tmp1)
 tmp2=summaryBy(index~slot1,data=tmp1,FUN=min)
 tmp1=merge(tmp1,tmp2)
 tmp1$slot2=tmp1$index-tmp1$index.min+1
 #tmp1
 #####################
 tmp1=tmp1[,c('slot1','slot2')]
 #tmp1
 tmp1=as.matrix(tmp1)
 (x=tmp1[1,])

 # function to pull data from polygon i  subslot j
 myfun1=function(x,mypoly=NA){
  i=as.numeric(x[1])
  j=as.numeric(x[2])
  tmp2=slot(slot((slot(mypoly, 'polygons'))[[i]],'Polygons')[[j]], 'coords')
  x=tmp2[,1]
  y=tmp2[,2]
  slot1=rep(i,nrow(tmp2))
  slot2=rep(j,nrow(tmp2))
  list(slot1=slot1,slot2=slot2,x=x,y=y)
 }

 x=tmp1[1,]
 myfun1(x,mypoly=shapes)
 # pull the points - the results will be in a list
 tmp2=apply(tmp1,1,myfun1,mypoly=shapes)

 #
 # compile/compress list items into a data frame
 df1=do.call(rbind, lapply(tmp2, function(z) cbind(z$slot1,z$slot2,z$x,z$y)))
 df1=as.data.frame(df1)
 names(df1)=c('slot1','slot2','x1','y1')
 ######################################################
 # Generate matrix with extrapolated points between subpolygon extreme points
 df1$slot11=my_lag(df1$slot1)
 df1$slot22=my_lag(df1$slot2)
 df1$x2=my_lag(df1$x1)
 df1$y2=my_lag(df1$y1)

 df2=subset(df1,slot1==slot11&slot2==slot22)
 df2$xdif=df2$x2-df2$x1
 df2$ydif=df2$y2-df2$y1
 ##########################
 # extrapolated extreme points
 theta=seq(0,1,length.out=npoints)
 X2=matrix(NA,nrow(df2),npoints)
 Y2=matrix(NA,nrow(df2),npoints)
 for(j in 1:npoints){
  X2[,j]=df2$x1+theta[j]*df2$xdif
  Y2[,j]=df2$y1+theta[j]*df2$ydif
 }
 # compute sin and cos of common variables
 # to avoid need to recompute again for every point in X0
 longrad2=X2/57.2958
 latrad2=Y2/57.2958
 sin.latrad2=sin(latrad2)
 cos.latrad2=cos(latrad2)
 ###############################################################################
 # put data in list
 mdata=list(longlat2=longlat2,df2=df2,X2=X2,Y2=Y2,longrad2=longrad2,
 sin.latrad2=sin.latrad2,cos.latrad2=cos.latrad2)
 ###############################################################################
 Polyclose=function(x,mdata=mdata,polynum=polynum){
  (long0=as.numeric(x[1]))
  (lat0=as.numeric(x[2]))

  longrad0=long0/57.2958
  latrad0=lat0/57.2958


  longlat2=mdata$longlat2
  longlat2$dist=dfunction(long0,lat0,longlat2$long,longlat2$lat)
  longlat2=orderBy(~dist,data=longlat2)
  # pull closets polynum slots
  longlat3=longlat2[1:polynum,]


  keep=!is.na(match(mdata$df2$slot1,longlat3$slot1))

  df3=mdata$df2[keep,]
  sin.latrad3=mdata$sin.latrad2[keep,]
  cos.latrad3=mdata$cos.latrad2[keep,]
  longrad3=mdata$longrad2[keep,]
  X3=mdata$X2[keep,]
  Y3=mdata$Y2[keep,]

  cosD=(sin(latrad0)*sin.latrad3+cos(latrad0)*cos.latrad3*cos(longrad0-longrad3))
  D=acos(cosD)
  distance=3963*D
  (dmin=min(distance))
  if(dmin=="NaN") dmin=0
  #distance[distance=="NaN"]=0
  (tmp=which(distance==dmin))
  (jmin=as.integer((tmp[1]-.1)/nrow(distance))+1)
  (imin=tmp[1]-(jmin-1)*nrow(distance))
  distance[imin,jmin]
  dmin
  x3=X3[imin,jmin]
  y3=Y3[imin,jmin]
  (slot1=df3$slot1[imin])

  c(slot1,x3,y3,dmin)
 }                                  # end function Polyclose
########################################################################
 #Find which points are inPolys
 X1=X0
 # find number of points
 nobpoints=nrow(X1)
 X1=data.frame(X1)
 names(X1)=c('x0','y0')
 # Find which points are in polygons and which are not
 dftmp=(X1)
 coordinates(dftmp) = ~x0+y0
 proj4string(dftmp) = proj4string(shapes)
 X1$polyin=over(dftmp,shapes)$index

 (x=as.matrix(X1[1,1:2]))
 Polyclose(x,mdata=mdata,polynum=polynum)

 if(plotum!=T){
  (X12=as.matrix(X1[,c('x0','y0')]))
 (tmp=t(apply(X12,1,Polyclose,mdata=mdata,polynum=polynum)))
  x2=tmp[,2]
  y2=tmp[,3]
  dist2=tmp[,4]
  polyclose=as.integer(tmp[,1])

 X2=data.frame(cbind(X1,polyclose,x2,y2,dist2))

 }                          # end if on plotum !=T



 if(plotum==T){
  X1block=list(NULL)
  nobsX=nrow(X1)
  js=1
  for(js in 1:nobsX){
   x1=X1[js,]
   x0=x1$x0[1]
   y0=x1$y0[1]

   (polyin=x1$polyin[1])

   tmp=Polyclose(c(x0,y0),mdata=mdata,polynum=polynum)
   x2=tmp[2]; y2=tmp[3]; dist2=tmp[4]; polyclose=as.integer(tmp[1])
   blockup=cbind(x1,polyclose,x2,y2,dist2)
#   X1block=rbind(X1block,blockup)
   X1block[[js]]=blockup 
   if(plotum!=F) {
    xlim=c(x0-pointbuffer,x0+pointbuffer)
    ylim=c(y0-pointbuffer,y0+pointbuffer)
    shapes@data$color=3
    polypick=ifelse(is.na(polyin),polyclose,polyin)
    shapes@data$color[polypick]=2
    if(js==1|js%%plotmod==0|js==nobsX){
     if(newmap==T) x11()
     sp::plot(shapes,col=shapes@data$color,xlim=xlim,ylim=ylim)
     mtext(paste(mytext,'POINT NPOINTS',js,nobsX),side=3,font=2,cex=2)
     points(x0,y0,pch=3,cex=3,lwd=2,font=2)
     wait(nsecs=nsecs)
    } # end  if(js==1|js%%plotmod==0|js==nobsX)
   }                 # end if on plotum

#  wait(nsecs=nsecs)

  }                  # end loop on js locations

#  X2=data.frame(X1block)
  X2=as.data.frame(bind_rows(X1block))
  
 }                                  # end if on plotum == T


 X2$polyclose=ifelse(!is.na(X2$polyin),X2$polyin,X2$polyclose)


return(X2)
}
# end function PolyfindM
######################################################################################

################################################################################
#experiments
#require(maptools)
#require(jcode64)
#require(doBy)
#load('C:/maps/mapdata/Maps.RData')
#require(maptools)
################################################################################
## Not Run
## long-lat of Bozeman Montana     = -111.047,45.68
## long-lat of Seattle Washington  = -122.35,47.62
## long-lat of Atlanta, Georgia    =  -84.39,33.75
#
## long-lat of Bozeman Montana     = -111.047,45.68
## long-lat of Seattle Washington  = -122.35,47.62
## long-lat of Atlanta, Georgia    =  -84.39,33.75
## location off Georgia Coast      =  -81.09,31.35
#
#data(USMapsS48)
#x0=c(-111.047,-122.35,-84.39,-81.09)
#y0=c(45.68,47.62,33.75,31.35)
#X=data.frame(cbind(x0,y0))
#graphics.off()
#(tmp=PolyfindplotM(X,USMapsS48,plotum=T,nsecs=0,newmap=F))
#        x0    y0 polyin polyclose         x2       y2     dist2
#1 -111.047 45.68      2         2 -111.05343 44.99570 47.332577
#2 -122.350 47.62      1         1 -122.38908 47.60989  1.951914
#3  -84.390 33.75     44        44  -85.33210 33.63763 54.770122
#4  -81.090 31.35     NA        44  -81.23607 31.43287 10.355521
################################################################################
#X0=X
#shapes=USMapsS48
#polynum=5
#npoints=10
#plotum=T
#pointbuffer=5
#nsecs=2
#dfunction=distcalc
#mytext=''
##################
#require(jcode64)
#data(USMapsS48)
#plot(USMapsS48)
#x0=c(-86.70)
#y0=c(41.046)
#X0=data.frame(cbind(x0,y0))
#################
#shapes=USMapsS48
#polynum=3
#npoints=10
#plotum=T
#pointbuffer=10
#dfunction=distcalc
#tmp=PolyfindplotM(X0,shapes,polynum=5,npoints=10,plotum=T,pointbuffer=5)
#my_seconds()
#X1=PolyfindM(X0,CAMaps,polynum=2,npoints=10)
#my_seconds()
#X1=PolyfindM(X0,CAMaps,polynum=5,npoints=10)
#my_seconds()
#X1=PolyfindM(X0,CAMaps,polynum=10,npoints=10)
#my_seconds()
################################################################################



################################################################################
distcalc=function(long0,lat0,long1,lat1,rads=F){
  if(rads!=T) {
  longrad0=long0/57.2958;latrad0=lat0/57.2958
  longrad1=long1/57.2958;latrad1=lat1/57.2958}
  if(rads==T) {
  longrad0=long0;latrad0=lat0
  longrad1=long1;latrad1=lat1}

  cosD=(sin(latrad0)*sin(latrad1)+cos(latrad0)*cos(latrad1)*cos(longrad0-longrad1))
  D=acos(cosD)
  distance=3963*D
  distance=ifelse(distance=="NaN",0,distance)
  return(distance)
} # end function distcalc
################################################################################



################################################################################
SP2SPDF = function(SP){
  # convert SpatialPolygons object to SpatialPolygonsDataFrame object
  (ncells=length(SP))
  tmp=as.data.frame(coordinates(SP))
  names(tmp)=c('long','lat')
  index=1:ncells
  
  df0=as.data.frame(cbind(tmp,index))
  head(df0)
  
  
  SPDF=SpatialPolygonsDataFrame(SP,df0)
} # end  function SP2SPDF
################################################################################


################################################################################
gradn =function(obj,b,epsmin=.00000001,...){
# print('gradum line 3')
 nobs=length(b)
 g=b*0
 for(i in 1:nobs){
  b1=b
  eps=max(abs(epsmin*b[i]),epsmin)
  b1[i]=b[i]-eps
  g[i]=(obj(b1,...)-obj(b,...))/eps
 }
 g
}
########################################
gradp = function(obj,b,epsmin=.00000001,...){
# print('gradum line 3')
 nobs=length(b)
 g=b*0
 for(i in 1:nobs){
  b1=b
  eps=max(abs(epsmin*b[i]),epsmin)
  b1[i]=b[i]+eps
  g[i]=(obj(b1,...)-obj(b,...))/eps
 }
 g
}
########################################
gradum=function(obj,b,epsmin=.00000001,...){
# print('gradum line 3')
 nobs=length(b)
 g=b*0
 for(i in 1:nobs){
  b1=b
  b2=b
  eps=max(abs(epsmin*b[i]),epsmin)
  b1[i]=b[i]+eps
  b2[i]=b[i]-eps
  g[i]=(obj(b1,...)-obj(b2,...))/(2*eps)
 }
 g
}
########################################
hessum = function(obj,b,epsmin=.00000001,...){
# print('hessum line 16')
 eps=sqrt(epsmin)
 nobs=length(b)
 H=matrix(0,nobs,nobs)
 ID=diag(1,nobs)
 i=1
 j=1
 for(i in 1:nobs){
  for(j in 1:i){
   eps1=max(abs(eps*b[i]),eps)
   eps2=max(abs(eps*b[j]),eps)
   hij=(obj(b+eps1*ID[,i]+eps2*ID[,j],...)-obj(b+eps1*ID[,i]-eps2*ID[,j],...)-
   obj(b-eps1*ID[,i]+eps2*ID[,j],...)
      + obj(b-eps1*ID[,i]-eps2*ID[,j],...))/(4*eps1*eps2)
   (H[i,j]=hij)
   (H[j,i]=hij)
  }
 }
 H
} # end function hessum
###############################################################################
###############################################################################
maxgrad=function(f0,b,itlim=500,epstop=0.0000001,
        b.low=rep(-10^10,length(b)),b.up=rep(10^10,length(b)),...)   {
###############################################################################
 require(MASS)
# source("grad-hess.r")
 npars=length(b)
 stopum=0
 iter=1
 (fval0=f0(b,...))
 fhist=fval0
 while(stopum==0&iter<=itlim) {         # first while
  (fval0=f0(b,...))
  (grad=gradum(f0,b,...))
  (gradtest=as.numeric(sum(abs(grad))))
  if(gradtest<=epstop){stopum=1}
  # start loops on steps
  loopgo=1
  stepi=1.0
  while(loopgo==1 & stepi>=epstop){     # second while
  (b1=b+stepi*grad)
  b1=ifelse(b1<b.low,b.low,b1)
  b1=ifelse(b1>b.up,b.up,b1)
  b1
  (fval1=f0(b1,...))

  if(fval1>(fval0+epstop)) {b=b1; loopgo=0}
   (stepi=0.1*stepi)
   }                                    # end second while
##########################################################

# check in individual gradient directions

   loopgo
  (fval0=f0(b,...))
  (grad=gradum(f0,b,...))
  (gradtest=as.numeric(sum(abs(grad))))

if(loopgo==1 & gradtest>0) {
  stepi=1.0
  while(loopgo==1 & stepi>=epstop){     # third while
   j=0
   while(j<npars){
    (j=j+1)
    b1=b
    (b1[j]=b[j]+stepi*grad[j])
    b1=ifelse(b1<b.low,b.low,b1)
    b1=ifelse(b1>b.up,b.up,b1)
    b
    b1
    (fval1=f0(b1,...))
    fval0
    if(fval1>(fval0+epstop)) {fval0=fval1; b=b1; loopgo=0; j=npars}
    loopgo
   }                                    # end loop on npars

     (stepi=0.1*stepi)
  }                                     # end third while

}

loopgo

##########################################


  if(loopgo==1) stopum=1
  iter=iter+1
  (fval0=f0(b,...))
  fhist[iter]=fval0
 } # end first while


  fval1
  (fval0=f0(b,...))
  (grad=gradum(f0,b,...))
  (gradtest=as.numeric(sum(abs(grad))))

  b
  (fval0=f0(b,...))
  (grad=gradum(f0,b,...))
  (gradtest=as.numeric(sum(abs(grad))))
  (hess0=hessum(f0,b,...))
  (sigma=-ginv(hess0))
  status=1
  if(gradtest<=epstop) {status=0}
  status

list(b=b,obj=fval0,grad=grad,hess=hess0,sigma=sigma,status=status,iter=iter)
} # end function maxgrad
###############################################################################

################################################################################
maxhess=function(f0,b,itlim=500,epstop=0.0000001,
   b.low=rep(-10^10,length(b)),b.up=rep(10^10,length(b)),...)   {
################################################################################
 require(MASS)
# source("grad-hess.r")
 npars=length(b)
 stopum=0
 iter=1
 (fval0=f0(b,...))
 fhist=fval0
 while(stopum==0&iter<=itlim) {         # first while
  (fval0=f0(b,...))
  (grad=gradum(f0,b,...))
  (gradtest=as.numeric(sum(abs(grad))))
  (grad=gradum(f0,b,...))
  (hess=hessum(f0,b,...))
  if(npars==1) {hessnd=-1*abs(hess);hesstest=hess}
  if(npars>=2) {
   (ehess=eigen(hess))
   (v=ehess$vectors)
   (d=diag(ehess$values))
   (hesstest=max(d))
   (d=-1*abs(d)  )
   (hessnd=v %*% d %*% t(v))
   }
  (if(gradtest<=epstop&hesstest<0){stopum=1})
  d=-1* ginv(hessnd) %*% grad
  # start loops on steps
  loopgo=1
  stepi=1.0
  while(loopgo==1 & stepi>=epstop){     # second while
  (b1=b+stepi*d)
  b1=ifelse(b1<b.low,b.low,b1)
  b1=ifelse(b1>b.up,b.up,b1)
  b1
  (fval1=f0(b1,...))

  if(fval1>(fval0+epstop)) {b=b1; loopgo=0}
   stepi=0.1*stepi
   }                                    # end second while

   stepi
   epstop
   fval0
   fval1


  if(loopgo==1) stopum=1
  iter=iter+1
  (fval0=f0(b,...))
  fhist[iter]=fval0
 }                                      # end first while

  b
  (fval0=f0(b,...))
  (grad=gradum(f0,b,...))
  (gradtest=as.numeric(sum(abs(grad))))
  (hess=hessum(f0,b,...))
  (sigma=-ginv(hess))
   ehess=eigen(hess)
   v=ehess$vectors
   (hesstest=max(ehess$values))
  status=1
  if(gradtest<=epstop) {status=0}
  status

list(b=b,obj=fval0,grad=grad,hess=hess,sigma=sigma,status=status,iter=iter,
     hesstest=hesstest)
} # end function maxhess
###############################################################################



###############################################################################
maxgdhess = function(obj,b,itlim=100,epstop=0.0000001,plotum=F,
ptitle='',ndirs=1,useDIRg=F,...)   {
###############################################################################
 #library(MASS)
 #library(numDeriv)
 dirgen=function(w,d1,d2){(1-w)*d1+w*d2}
  w=as.matrix(seq(0,1,1/ndirs))
 npars=length(b)
 one.w=matrix(1,(length(w)),1)
 if(useDIRg==T) one.w=matrix(1,(length(w)+npars),1)
 stopum=0
 iter=1
 (fval0=obj(b,...))
# (fval0=obj(b))
 fhist=fval0

 while(stopum==0&iter<=itlim) {         # first while
  (fval0=obj(b,...))
#  (fval0=obj(b))
  (gradb=grad(obj,b,...))
#  (gradb=grad(obj,b))
  (gradtest=as.numeric(sum(abs(gradb))))

  (hess=hessian(obj,b,...))
#  (hess=hessian(obj,b))

  if(npars==1) {hessnd=-1*abs(hess);hesstest=hess}
  if(npars>=2) {
   (ehess=eigen(hess))
   (v=ehess$vectors)
   (d=diag(ehess$values))
   (hesstest=max(d))
   (d=-1*abs(d)  )
   (hessnd=v %*% d %*% t(v))
   }
  (if(gradtest<=epstop&hesstest<0){stopum=1})
  d1=-1* ginv(hessnd) %*% gradb
  d2=gradb
  (DIR=apply(w,1,'dirgen',d1=d1,d2=d2))  # generates an npar by ndirs+1 matrix of directions
  (fval0=obj(b,...))
#  (fval0=obj(b))
  (gradb=grad(obj,b,...))
#  (gradb=grad(obj,b))
  (gradtest=as.numeric(sum(abs(gradb))))

  (hess=hessian(obj,b,...))
#  (hess=hessian(obj,b))

  if(npars==1) {hessnd=-1*abs(hess);hesstest=hess}
  if(npars>=2) {
   (ehess=eigen(hess))
   (v=ehess$vectors)
   (d=diag(ehess$values))
   (hesstest=max(d))
   (d=-1*abs(d)  )
   (hessnd=v %*% d %*% t(v))
   }
  (if(gradtest<=epstop&hesstest<0){stopum=1})
  d1=-1* ginv(hessnd) %*% gradb
  d2=gradb
  (DIR=apply(w,1,'dirgen',d1=d1,d2=d2))  # generates matrix of directions

   DIRg=diag(gradb)
   if(useDIRg==T)  DIR=cbind(DIR,DIRg)

  B=b %*% t(one.w)

  # start loops on steps
  loopgo=1
  stepi=1.0

  while(loopgo==1 & stepi>=epstop){     # second while

  (B1=B+stepi*DIR)
  FVAL=apply(B1,2,'obj',...)
#  (FVAL=apply(B1,2,'obj'))

  (fval1=max(FVAL))

  if(fval1>(fval0+epstop)) {b=B1[,min(which(FVAL==fval1))]; loopgo=0}
   stepi=0.1*stepi
   }                                    # end second while

   stepi
   epstop
   fval0
   fval1

  if(loopgo==1) stopum=1
  iter=iter+1
  (fval0=obj(b,...))
#  (fval0=obj(b))

  fhist[iter]=fval0
if(plotum==T&(iter%%10==0|iter>=itlim)) plot(1:iter,fhist[1:iter],
     main=paste(ptitle,'iter',iter,round(fval0,2),round(gradtest,2)))

 }                                      # end first while

  (fval0=obj(b,...))
#  (fval0=obj(b))

  (gradb=grad(obj,b,...))
#  (gradb=grad(obj,b))

  (gradpos=gradp(obj,b,...))
  (gradneg=gradn(obj,b,...))


  (gradtest=as.numeric(sum(abs(gradb))))
  (hess=hessian(obj,b,...))
#  (hess=hessian(obj,b))

  (sigma=-ginv(hess))
   ehess=eigen(hess)
   v=ehess$vectors
   (hesstest=max(ehess$values))
  status=1
  if(gradtest<=epstop) {status=0}
  status


#  dftmp=obj(b,returnpars=T,...)
#  b.restrict=dftmp$b
#print(b.restrict)


returnum=list(b=b,obj=fval0,gradb=gradb,gradpos=gradpos,gradneg=gradneg,
hess=hess,sigma=sigma,status=status,iter=iter,hesstest=hesstest)
returnum
} # end function maxgdhess
###############################################################################

###############################################################################
maxgdhess2=function(obj,b,itlim=100,epstop=0.0000001,plotum=F,ptitle='',
           ndirs=1,useDIRg=F,...)   {
###############################################################################
# library(MASS)
# library(numDeriv)
 dirgen=function(w,d1,d2){(1-w)*d1+w*d2}
  w=as.matrix(seq(0,1,1/ndirs))
 npars=length(b)
 one.w=matrix(1,(length(w)),1)
 if(useDIRg==T) one.w=matrix(1,(length(w)+npars),1)
 stopum=0
 iter=1
 (fval0=obj(b,...))
# (fval0=obj(b))
 fhist=fval0

 while(stopum==0&iter<=itlim) {         # first while


#  dftmp=obj(b,returnpars=T,...)
#  b.restrict=dftmp$b
#print(b.restrict)
#  b=b.restrict


  (fval0=obj(b,...))
#  (fval0=obj(b))
  (gradb=gradum(obj,b,...))
#  (gradb=gradum(obj,b))
  (gradtest=as.numeric(sum(abs(gradb))))

  (hess=hessum(obj,b,...))
#  (hess=hessum(obj,b))

  if(npars==1) {hessnd=-1*abs(hess);hesstest=hess}
  if(npars>=2) {
   (ehess=eigen(hess))
   (v=ehess$vectors)
   (d=diag(ehess$values))
   (hesstest=max(d))
   (d=-1*abs(d)  )
   (hessnd=v %*% d %*% t(v))
   }
  (if(gradtest<=epstop&hesstest<0){stopum=1})
  d1=-1* ginv(hessnd) %*% gradb
  d2=gradb
  (DIR=apply(w,1,'dirgen',d1=d1,d2=d2))  # generates an npar by ndirs+1 matrix of directions
  (fval0=obj(b,...))
#  (fval0=obj(b))
  (gradb=gradum(obj,b,...))
#  (gradb=gradum(obj,b))
  (gradtest=as.numeric(sum(abs(gradb))))

  (hess=hessum(obj,b,...))
#  (hess=hessum(obj,b))

  if(npars==1) {hessnd=-1*abs(hess);hesstest=hess}
  if(npars>=2) {
   (ehess=eigen(hess))
   (v=ehess$vectors)
   (d=diag(ehess$values))
   (hesstest=max(d))
   (d=-1*abs(d)  )
   (hessnd=v %*% d %*% t(v))
   }
  (if(gradtest<=epstop&hesstest<0){stopum=1})
  d1=-1* ginv(hessnd) %*% gradb
  d2=gradb
  (DIR=apply(w,1,'dirgen',d1=d1,d2=d2))  # generates an npar by ndirs+1 matrix of directions

   DIRg=diag(gradb)
   if(useDIRg==T)  DIR=cbind(DIR,DIRg)

  B=b %*% t(one.w)

  # start loops on steps
  loopgo=1
  stepi=1.0

  while(loopgo==1 & stepi>=epstop){     # second while

  (B1=B+stepi*DIR)
  FVAL=apply(B1,2,'obj',...)
#  (FVAL=apply(B1,2,'obj'))


  (fval1=max(FVAL))

  if(fval1>(fval0+epstop)) {b=B1[,min(which(FVAL==fval1))]; loopgo=0}
   stepi=0.1*stepi
   }                                    # end second while

   stepi
   epstop
   fval0
   fval1


  if(loopgo==1) stopum=1
  iter=iter+1
  (fval0=obj(b,...))
#  (fval0=obj(b))

  fhist[iter]=fval0
if(plotum==T&(iter%%10==0|iter>=itlim)) plot(1:iter,fhist[1:iter],
  main=paste(ptitle,'iter',iter,round(fval0,2),round(gradtest,2)))

 }                                      # end first while

  (fval0=obj(b,...))
#  (fval0=obj(b))

  (gradb=gradum(obj,b,...))
#  (gradb=gradum(obj,b))

  (gradpos=gradp(obj,b,...))
  (gradneg=gradn(obj,b,...))


  (gradtest=as.numeric(sum(abs(gradb))))
  (hess=hessum(obj,b,...))
#  (hess=hessum(obj,b))

  (sigma=-ginv(hess))
   ehess=eigen(hess)
   v=ehess$vectors
   (hesstest=max(ehess$values))
  status=1
  if(gradtest<=epstop) {status=0}
  status


#  dftmp=obj(b,returnpars=T,...)
#  b.restrict=dftmp$b
#print(b.restrict)


returnum=list(b=b,obj=fval0,gradb=gradb,gradpos=gradpos,gradneg=gradneg,
hess=hess,sigma=sigma,status=status,iter=iter,hesstest=hesstest)
returnum
} # end function maxgdhess2
################################################################################

################################################################################
maxgdhessr = function(obj,b,itlim=100,epstop=0.0000001,plotum=F,ptitle='',
              ndirs=1,useDIRg=F,...)   {
################################################################################
# library(MASS)
# library(numDeriv)
 dirgen=function(w,d1,d2){(1-w)*d1+w*d2}
  w=as.matrix(seq(0,1,1/ndirs))
 npars=length(b)
 one.w=matrix(1,(length(w)),1)
 if(useDIRg==T) one.w=matrix(1,(length(w)+npars),1)
 stopum=0
 iter=1
 (fval0=obj(b,...))
# (fval0=obj(b))
 fhist=fval0

 while(stopum==0&iter<=itlim) {         # first while

#debugum
  dftmp=obj(b,returnpars=T,...)
  b.restrict=dftmp$b
  b=b.restrict


  (fval0=obj(b,...))
#  (fval0=obj(b))
  (gradb=grad(obj,b,...))
#  (gradb=grad(obj,b))


  (gradtest=as.numeric(sum(abs(gradb))))

  (hess=hessian(obj,b,...))
#  (hess=hessian(obj,b))

  if(npars==1) {hessnd=-1*abs(hess);hesstest=hess}
  if(npars>=2) {
   (ehess=eigen(hess))
   (v=ehess$vectors)
   (d=diag(ehess$values))
   (hesstest=max(d))
   (d=-1*abs(d)  )
   (hessnd=v %*% d %*% t(v))
   }
  (if(gradtest<=epstop&hesstest<0){stopum=1})
  d1=-1* ginv(hessnd) %*% gradb
  d2=gradb
  (DIR=apply(w,1,'dirgen',d1=d1,d2=d2))  # generates matrix of directions
  (fval0=obj(b,...))
#  (fval0=obj(b))
  (gradb=grad(obj,b,...))
#  (gradb=grad(obj,b))
  (gradtest=as.numeric(sum(abs(gradb))))

  (hess=hessian(obj,b,...))
#  (hess=hessian(obj,b))

  if(npars==1) {hessnd=-1*abs(hess);hesstest=hess}
  if(npars>=2) {
   (ehess=eigen(hess))
   (v=ehess$vectors)
   (d=diag(ehess$values))
   (hesstest=max(d))
   (d=-1*abs(d)  )
   (hessnd=v %*% d %*% t(v))
   }
  (if(gradtest<=epstop&hesstest<0){stopum=1})
  d1=-1* ginv(hessnd) %*% gradb
  d2=gradb
  (DIR=apply(w,1,'dirgen',d1=d1,d2=d2))  # generates matrix of directions

   DIRg=diag(gradb)
   if(useDIRg==T)  DIR=cbind(DIR,DIRg)


  B=b %*% t(one.w)
  # start loops on steps
  loopgo=1
  stepi=1.0

  while(loopgo==1 & stepi>=epstop){     # second while

  (B1=B+stepi*DIR)
  FVAL=apply(B1,2,'obj',...)
#  (FVAL=apply(B1,2,'obj'))


  (fval1=max(FVAL))

  if(fval1>(fval0+epstop)) {b=B1[,min(which(FVAL==fval1))]; loopgo=0}
   stepi=0.1*stepi
   }                                    # end second while

   stepi
   epstop
   fval0
   fval1


  if(loopgo==1) stopum=1
  iter=iter+1
  (fval0=obj(b,...))
#  (fval0=obj(b))

  fhist[iter]=fval0
if(plotum==T&(iter%%10==0|iter>=itlim))  plot(1:iter,fhist[1:iter],
 main=paste(ptitle,'iter',iter,round(fval0,2),round(gradtest,6)))

 }                                      # end first while

  (fval0=obj(b,...))
#  (fval0=obj(b))

  (gradb=grad(obj,b,...))
#  (gradb=grad(obj,b))

  (gradpos=gradp(obj,b,...))
#  (gradpos=gradp(obj,b))

  (gradneg=gradn(obj,b,...))
#  (gradneg=gradn(obj,b))


  (gradtest=as.numeric(sum(abs(gradb))))
  (hess=hessian(obj,b,...))
#  (hess=hessian(obj,b))

  (sigma=-ginv(hess))
   ehess=eigen(hess)
   v=ehess$vectors
   (hesstest=max(ehess$values))
  status=1
  if(gradtest<=epstop) {status=0}
  status


  dftmp=obj(b,returnpars=T,...)
  b.restrict=dftmp$b
#print(b.restrict)


list(b=b,obj=fval0,gradb=gradb,gradpos=gradpos,gradneg=gradneg,hess=hess,
    sigma=sigma,status=status,iter=iter,hesstest=hesstest,b.restrict=b.restrict)

} # end function maxgdhessr
################################################################################

################################################################################
maxgdhessr2=function(obj,b,itlim=100,epstop=0.0000001,plotum=F,
            ptitle='',ndirs=1,useDIRg=F,...)   {
################################################################################
# library(MASS)
# library(numDeriv)
 dirgen=function(w,d1,d2){(1-w)*d1+w*d2}
  w=as.matrix(seq(0,1,1/ndirs))
 npars=length(b)
 one.w=matrix(1,(length(w)),1)
 if(useDIRg==T) one.w=matrix(1,(length(w)+npars),1)
 stopum=0
 iter=1
 (fval0=obj(b,...))
# (fval0=obj(b))
 fhist=fval0

 while(stopum==0&iter<=itlim) {         # first while

#debugum
  dftmp=obj(b,returnpars=T,...)
  b.restrict=dftmp$b
  b=b.restrict


  (fval0=obj(b,...))
#  (fval0=obj(b))
  (gradb=gradum(obj,b,...))
#  (gradb=gradum(obj,b))

  (gradtest=as.numeric(sum(abs(gradb))))

  (hess=hessum(obj,b,...))
#  (hess=hessum(obj,b))

  if(npars==1) {hessnd=-1*abs(hess);hesstest=hess}
  if(npars>=2) {
   (ehess=eigen(hess))
   (v=ehess$vectors)
   (d=diag(ehess$values))
   (hesstest=max(d))
   (d=-1*abs(d)  )
   (hessnd=v %*% d %*% t(v))
   }
  (if(gradtest<=epstop&hesstest<0){stopum=1})
  d1=-1* ginv(hessnd) %*% gradb
  d2=gradb
  (DIR=apply(w,1,'dirgen',d1=d1,d2=d2))  # generates an matrix of directions
  (fval0=obj(b,...))
#  (fval0=obj(b))
  (gradb=gradum(obj,b,...))
#  (gradb=gradum(obj,b))
  (gradtest=as.numeric(sum(abs(gradb))))

  (hess=hessum(obj,b,...))
#  (hess=hessum(obj,b))

  if(npars==1) {hessnd=-1*abs(hess);hesstest=hess}
  if(npars>=2) {
   (ehess=eigen(hess))
   (v=ehess$vectors)
   (d=diag(ehess$values))
   (hesstest=max(d))
   (d=-1*abs(d)  )
   (hessnd=v %*% d %*% t(v))
   }
  (if(gradtest<=epstop&hesstest<0){stopum=1})
  d1=-1* ginv(hessnd) %*% gradb
  d2=gradb
  (DIR=apply(w,1,'dirgen',d1=d1,d2=d2))  # generates matrix of directions

   DIRg=diag(gradb)
   if(useDIRg==T)  DIR=cbind(DIR,DIRg)


  B=b %*% t(one.w)
  # start loops on steps
  loopgo=1
  stepi=1.0

  while(loopgo==1 & stepi>=epstop){     # second while

  (B1=B+stepi*DIR)
  FVAL=apply(B1,2,'obj',...)
#  (FVAL=apply(B1,2,'obj'))


  (fval1=max(FVAL))

  if(fval1>(fval0+epstop)) {b=B1[,min(which(FVAL==fval1))]; loopgo=0}
   stepi=0.1*stepi
   }                                    # end second while

   stepi
   epstop
   fval0
   fval1


  if(loopgo==1) stopum=1
  iter=iter+1
  (fval0=obj(b,...))
#  (fval0=obj(b))

  fhist[iter]=fval0
if(plotum==T&(iter%%10==0|iter>=itlim))  plot(1:iter,fhist[1:iter],
  main=paste(ptitle,'iter',iter,round(fval0,2),round(gradtest,6)))

 }                                      # end first while

  (fval0=obj(b,...))
#  (fval0=obj(b))

  (gradb=gradum(obj,b,...))
#  (gradb=gradum(obj,b))

  (gradpos=gradp(obj,b,...))
#  (gradpos=gradp(obj,b))

  (gradneg=gradn(obj,b,...))
#  (gradneg=gradn(obj,b))


  (gradtest=as.numeric(sum(abs(gradb))))
  (hess=hessum(obj,b,...))
#  (hess=hessum(obj,b))

  (sigma=-ginv(hess))
   ehess=eigen(hess)
   v=ehess$vectors
   (hesstest=max(ehess$values))
  status=1
  if(gradtest<=epstop) {status=0}
  status


  dftmp=obj(b,returnpars=T,...)
  b.restrict=dftmp$b
#print(b.restrict)

returnum=list(b=b,obj=fval0,gradb=gradb,gradpos=gradpos,gradneg=gradneg,
hess=hess,sigma=sigma,status=status,iter=iter,hesstest=hesstest,
b.restrict=b.restrict)
returnum

} # end function maxgdhessr2
###############################################################################


###############################################################################
my_over=function(XY,map0) {
 XY = as.data.frame(XY)
 names(XY)=c("x0","y0")
 coordinates(XY) = ~x0+y0
class(XY)
class(map0)

proj4string(XY) = proj4string(map0)
return(over(XY,map0))
} # end function my_over
################################################################################
 wait = function(nsecs=5){
 sec0=my_seconds()
 secdif=0
 while(secdif<nsecs){
  secdif=my_seconds()-sec0
 }
} # end function wait
########################################
my_seconds = function(){
 time1=proc.time()
 as.numeric(time1[3])
}
################################################################################


################################################################################
#http://conjugateprior.org/2015/06/identifying-the-os-from-r/
################################################################################
get_os = function(){
  sysinf = Sys.info()
  if (!is.null(sysinf)){
    os = sysinf['sysname']
    if (os == 'Darwin')
      os = "osx"
  } else { ## mystery machine
    os = .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os = "osx"
    if (grepl("linux-gnu", R.version$os))
      os = "linux"
  }
  tolower(os)
}  # end function get_os
################################################################################
#(opsys = as.character(get_os()))
################################################################################

################################################################################
my_x11=function(display="",width=7,height=7,pointsize=12,
 gamma=1,bg="transparent",canvas="white",xpos=NA,ypos=NA,
 title=""){

  (opsys = as.character(get_os()))
  if(opsys!="osx"){
   X11(display=display,width=width,height=height,
     pointsize=pointsize,gamma=gamma,bg=bg,canvas=canvas,
     xpos=xpos,ypos=ypos,title=title)
   }
  else {
   quartz(title=title,width=width,height=height,
     pointsize=pointsize,bg=bg,canvas=canvas)
  }
} # end function my_x11
################################################################################
#my_x11()
################################################################################





################################################################################
colorum.gplots=function(scores,ncoLs=25,coL1='blue',coL2='yellow',coL3='red',
Nlgd=11,Lround=0,scoremin=-(10^10),scoremax=10^10){
################################################################################
colorpanel.local = function(n,low,mid,high){
    if(missing(mid) || missing(high) )
      {
        ## convert to rgb
        low = col2rgb(low)
        if(missing(high))
          high = col2rgb(mid)
        else
          high = col2rgb(high)

        red    = seq(low[1,1], high[1,1], length=n)/255
        green  = seq(low[3,1], high[3,1], length=n)/255
        blue   = seq(low[2,1], high[2,1], length=n)/255
      }
    else # use a center color
      {
        isodd = odd.local(n)
        if(isodd)
          {
            n = n+1
          }

        ## convert to rgb
        low = col2rgb(low)
        mid = col2rgb(mid)
        high = col2rgb(high)

        ## determine length of each component
        lower = floor(n/2)
        upper = n - lower

        red  = c(
                  seq(low[1,1], mid [1,1], length=lower),
                  seq(mid[1,1], high[1,1], length=upper)
                  )/255

        green = c(
                   seq(low[3,1], mid [3,1], length=lower),
                   seq(mid[3,1], high[3,1], length=upper)
                   )/255

        blue = c(
                  seq(low[2,1], mid [2,1], length=lower),
                  seq(mid[2,1], high[2,1], length=upper)
                  )/255

        if(isodd)
          {
            red   = red  [-(lower+1)]
            green = green[-(lower+1)]
            blue  = blue [-(lower+1)]
          }
      }

    rgb(red,blue,green)
  }

############################
odd.local=function(x){
 y=x %% 2
 y!=0
}
############################
even.local=function(x){
 y=x %% 2
 y==0
}
#############################

# Generate red-to-green colorscale
redgreen.local = function(n) colorpanel.local(n, 'red', 'black', 'green')
greenred.local = function(n) colorpanel.local(n, 'green', 'black', 'red' )
bluered.local  = function(n) colorpanel.local(n, 'blue','white','red')
redblue.local  = function(n) colorpanel.local(n, 'red','white','blue')
################################################################################
if(is.numeric(coL2)){
        if(abs(coL2)<0.00001) coL2=0.00001
} 

if(!is.null(coL2)) coLs=colorpanel.local(ncoLs,low=coL1,mid=coL2,high=coL3)
if(is.null(coL2))  coLs=colorpanel.local(ncoLs,low=coL1,high=coL3)

if(scoremin<=(-(10^10))){ if(scoremin<(min.na(scores))) scoremin=min.na(scores)}
if(scoremax>=10^10){ if(scoremax>max.na(scores)) scoremax=max.na(scores)}


 scores=ifelse(scores<scoremin,scoremin,scores)
 scores=ifelse(scores>scoremax,scoremax,scores)

 breaks=seq(scoremin,scoremax,length.out=ncoLs)
 fr=cut(scores,breaks=breaks,include=T)
 datacolors=coLs[fr]

 Lbreaks=seq(scoremin,scoremax,length.out=Nlgd)
 Lbreaks=round(Lbreaks,Lround)
 Lbreaks[1]=min.na(breaks)
 Lbreaks[Nlgd]=max.na(breaks)
 Lfr=cut(Lbreaks,breaks=breaks,include=T)

 Lcolors=coLs[Lfr]
 Ltext=paste(round(Lbreaks,Lround))

 Ltext[1]=paste('<=',Ltext[1],sep='')
 Ltext[Nlgd]=paste('>=',Ltext[Nlgd],sep='')

return(list(datacolors=datacolors,Lcolors=Lcolors,Ltext=Ltext))
} # end function colorum.gplots
################################################################################

################################################################################
##debugum
#N=50

#(ncoLs=25)
#(coL1='blue')
#(coL2='yellow')
#(coL3='red')
#(Nlgd=11)
#(Lround=2)
#(scoremin=-(10^10))
#(scoremax=10^10)
#scores=1:N
#
#
#joe=colorum.gplots(scores=1:N,Lround=0)
#
#plot(1:N,1:N,col=joe$datacolors,pch=20,cex=2)
#
#legend('bottomright',legend=joe$Ltext,fill=joe$Lcolors,cex=0.5)
################################################################################

################################################################################
min.na=function(x,skipzero=F){
 if(skipzero!=F) x=x[x!=0]
 x=x[!is.na(x)]
 nobsx=length.na(x)
 if(nobsx==0) minx=NA
 if(nobsx>=1) minx=min(x)
 minx
}
###########################
mean.na = function(x,skipzero=F){
 if(skipzero!=F) x=x[x!=0]
 x=x[!is.na(x)]
 nobsx=length.na(x)
 if(nobsx==0) meanx=NA
 if(nobsx>=1) meanx=mean(x)
 meanx
}
###########################
median.na = function(x,skipzero=F){
 if(skipzero!=F) x=x[x!=0]
 x=x[!is.na(x)]
 nobsx=length.na(x)
 if(nobsx==0) medx=NA
 if(nobsx>=1) medx=median(x)
 medx
}
###########################
max.na = function(x,skipzero=F){
 if(skipzero!=F) x=x[x!=0]
 x=x[!is.na(x)]
 nobsx=length.na(x)
 if(nobsx==0) maxx=NA
 if(nobsx>=1) maxx=max(x)
 maxx
}
###########################
sd.na = function(x,skipzero=F){
 if(skipzero!=F) x=x[x!=0]
 x=x[!is.na(x)]
 nobsx=length.na(x)
 if(nobsx==0) sdx=NA
 if(nobsx>=1) sdx=sd(x)
 sdx
}
###########################
sum.na=function(x){
 sumx=NA
 x=x[!is.na(x)]
 if(length(x)>=1) sumx=sum(x)
 sumx
}
############################
length.na=function(x,skipzero=F) {
 if(skipzero!=F) x=x[x!=0]
 x=x[!is.na(x)]
 length(x)
}
################################################################################


################################################################################
meanOA=function(x,ncut=1){
 xbar=NA
 x=x[!is.na(x)]
 n=length(x)
 if(n>=(2*ncut+1)){
  x=x[(ncut+1):(n-ncut)]
  xbar=mean(x)
 }
 xbar
}
################################################################################
#x=1:10
#x[2:3]=NA
#meanOA(x,1)
#meanOA(x,2)
#meanOA(x,4)
################################################################################

################################################################################
my_lag=function(x,nlag=1) {
  nobx=length(x)
  if(nlag<0) x=rev(x)
  x0=rep(NA,abs(nlag))
  x1=x[1:(nobx-abs(nlag))]
  lagx=c(x0,x1)
  if(nlag<0){
    x=rev(x)
    lagx=rev(lagx)
  }
  lagx
} # end function my_lag
###############################################################################

###############################################################################
lagMat=function(x,nlag=2) {
  (nobx=length(x))
  if(nlag<0) x=rev(x)
  
  X=matrix(0,nobx,(abs(nlag)+1))
  X[,1]=x
  j=1
  for(j in 1:abs(nlag)){
    (jon=j+1)
    (x0=rep(NA,j))
    (x1=x[1:(nobx-j)])
    (lagx=c(x0,x1))
    (X[,jon]=lagx)
  }
  
  if(nlag<0){
    x=rev(x)
    X=apply(X,2,rev)
  }
  X
} # end function lagMat
###############################################################################


###############################################################################
moveAVG=function(x,nlag=5,na.rm=T,include=T){
 tmp=lagMat(x,nlag)
 if(include==T) mavg=apply(tmp[,1:nlag],1,mean,na.rm=na.rm)
 if(include!=T) mavg=apply(tmp[,2:(nlag+1)],1,mean,na.rm=na.rm)
 mavg
}
################################################################################
#Note that the above code computes a moving
#average of the most recent nlag observations INCLUDING
#the last observation....Be sure that this function is
#doing what you want.  If you want the most "recent" value
#to be 'excluded' from the moving average set include=F
# example:
################################################################################
#x=c(1,2,3,4,5)
#moveAVG(x,3)
#moveAVG(x,3,include=F)
################################################################################




################################################################################
ulong=function(x,na.rm=F){
 if(na.rm==T) x=x[!is.na(x)]
 length(unique(x))
}

#x1=rep(1:5,each=3)
#x2=rep(1:5,3)
#x2[5]=NA
#ulong(x2)
#ulong(x2,na.rm=T)
#tmp=data.frame(x1,x2)
#summaryBy(x2~x1,data=tmp,FUN=ulong)
#summaryBy(x2~x1,data=tmp,FUN=ulong,na.rm=T)
################################################################################

################################################################################
pvalx=function(x,q=mean(x)){
 x=x[!is.na(x)]
 length(x[x<q])/length(x)
}
############################
#set.seed(2012)
#x=rnorm(100000,100,25)
#pvalx(x)
#pvalx(x,90)
#pnorm(90,100,25)
############################
pvalx_V=function(x,q){
 q=as.matrix(q)
 apply(q,1,pvalx,x=x)
} # end function Pvalx_V
############################
#x=rnorm(1000000,0,10)
#pvalx_V(x=x,q=c(-30,-20,-10))
#pnorm(c(-30,-20,-10),0,10)
################################################################################


################################################################################
cor2cov=function(R,V){
    if(is.matrix(V)) {
     nrc=dim(V)
     if(nrc[1]!=nrc[2]) V=as.vector(V)
     if(nrc[1]==nrc[2]) V = diag(V)
    } # end if
    s=sqrt(V)
    V2 = R
    V2[] = s * R * rep(s, each = nrow(R))
    V2
} # end function cor2cov
########################################
cor2rho=function(R){
  (m=ncol(R))
  (rho=R[1,2:m])
  (rnames=paste('rho_',1,'_',2:m,sep=''))
  if(m>2){
   i=2
   for(i in 2:(m-1)) {
    (rho=c(rho,R[i,(i+1):m]))
    (rnames=c(rnames,paste('rho_',i,'_',(i+1):m,sep='')))
   } # end for(for(i in 2:(m-1))
  } # end if(m>2)
  tmp=data.frame(t(rho))
  names(tmp)=rnames
  tmp
} # end function cor2rho
########################################
rho2cor=function(rho){
 (k=length(rho))
 (n=as.integer((1+sqrt(1+4*2*k))/2))
 R=matrix(0,n,n)
 diag(R)=1
 ron=0
 for(j in 1:(n-1)){
  for(i in (j+1):n){
   ron=ron+1
   R[i,j]=rho[ron]
   R[j,i]=rho[ron]
  }  # end loop on i
 }  # end loop on j
R
} # end function rho2cor
################################################################################

###############################################################################
my_upperT=function(x=1:10){
  (k=length(x))
#  (n=as.integer((1+sqrt(1+4*2*k))/2-1))
  (n=(sqrt(1+8*k)-1)/2)
  if((n-as.integer(n))>0) stop('Length of x not consistent with UT square matrix')
  R=matrix(0,n,n)
  ron=0
  for(i in 1:n){
    for(j in i:n){
      ron=ron+1
      R[i,j]=x[ron]
    }  # end loop on i
  }  # end loop on j
  R
} # end function my_upper
###############################################################################
my_lowerT=function(x){ 
t(my_upperT(x))
} # end function my_lowerT
###############################################################################
#x=1:6
#my_upperT(x);my_lowerT(x)
###############################################################################
#my_upperT(1:6)
#my_upperT(1:10)
#my_upperT(1:15)
#my_upperT(1:21)
###############################################################################
#my_upperT(1:11)
###############################################################################

###############################################################################
my_rank=function(x,ties.method='first',na.last='keep'){
 rank(x,ties.method=ties.method,na.last=na.last)  
} # end function my_rank
###############################################################################

###############################################################################
rankorder=function(Y,RM){
###############################################################################
 if(class(Y)[1]!='matrix')  stop('Y is not a matrix')
 if(class(RM)[1]!='matrix')  stop('RM is not a matrix')
 if(nrow(Y)!=nrow(RM)) stop('Y and RM must have an equal number of rows')
 if(ncol(Y)!=ncol(RM)) stop('Y and RM must have an equal number of columns')
 # put columns of Y in same rank order as the columns in RM  
 for(j in 1:ncol(Y)){
  ys=sort(Y[,j])
  Y[,j] = ys[rank(RM[,j],ties.method='first')]
 } # end loop for(j in 1:ncol(YC)) 
 Y  
###############################################################################  
}# end function rankorder
###############################################################################

###############################################################################
rankMAT = function(X){
############################################################################### 
 if(class(X)[1]!='matrix') stop('X is not a matrix')
 n=nrow(X); m = ncol(X)
 for(j in 1:m) X[,j] = rank(X[,j],ties.method='first')/(n+1)
 X
############################################################################### 
} # end function rankMAT
###############################################################################


###############################################################################
ImanConover = function(Y,sigma,Icor=F,seedval=NA,whiten=T) {
###############################################################################  
  if(ncol(Y)!=ncol(sigma)) stop('Dimensions of Y and sigma do not match')
  if(!is.na(seedval)) set.seed(seedval)
  (nc=ncol(Y))
  (nr=nrow(Y))
  Zi=matrix(rnorm(nr*nc),nr,nc)
  #############################################
  # spectral decompostion to guarantee positive definite sigma matrix
  tmp=eigen(sigma)
  v=tmp$vectors
  (d=tmp$values)
  d[d<=0.0001]=0.0001
  #############################################
  if(whiten!=T){
    VT=t(v)
    for(i in 1:nrow(VT)) VT[i,]=VT[i,]*sqrt(d[i])
    Zc=Zi %*% VT
  } # end if  if(whiten!=T)
  #############################################
  if(whiten==T){
   (d=diag(d))
   (sigma = v %*% d %*% t(v))  
   if(Icor==T) (sigma= cov2cor(sigma))
    Zi=Zi%*%solve(chol(cov(Zi)))
    for(j in 1:ncol(Zi)) Zi[,j] = Zi[,j] - mean(Zi[,j])
    Zc=Zi%*%(chol(sigma))
  } # end if(whiten==T)
  ############################################
  for(j in 1:nc) {
    Zc[,j]=Zc[,j]-mean(Zc[,j])
    ys=sort(Y[,j])
    Y[,j] = ys[rank(Zc[,j],ties.method="first",na.last='keep')]
    
  }  
  Y
} # end function ImanConover
################################################################################
ImanConoverYcZc=function(Y,sigma,Icor=F,seedval=NA,whiten=T) {
  if(ncol(Y)!=ncol(sigma)) stop('Dimensions of Y and sigma do not match')
  if(!is.na(seedval)) set.seed(seedval)
  (nc=ncol(Y))
  (nr=nrow(Y))
  Zi=matrix(rnorm(nr*nc),nr,nc)
  ##############################################
  # spectral decompostion to guarantee positive definite sigma matrix
  tmp=eigen(sigma)
  v=tmp$vectors
  (d=tmp$values)
  d[d<=0.0001]=0.0001
#  (d=diag(d))
  #############################################
  if(whiten!=T){
    VT=t(v)
    for(i in 1:nrow(VT)) VT[i,] = VT[i,]*sqrt(d[i])
    Zc=Zi %*% VT
  } # end if  if(whiten!=T)
  #############################################
  if(whiten==T){
    (d=diag(d))
    (sigma = v %*% d %*% t(v))
    if(Icor==T) (sigma= cov2cor(sigma))
    Zi=Zi%*%solve(chol(cov(Zi)))
    for(j in 1:ncol(Zi)) Zi[,j] = Zi[,j] - mean(Zi[,j])
    Zc=Zi%*%(chol(sigma))
  } # end if(whiten==T)
  ############################################
  for(j in 1:nc) {
    Zc[,j]=Zc[,j]-mean(Zc[,j])
    ys=sort(Y[,j])
    Y[,j] = ys[rank(Zc[,j],ties.method="first",na.last='keep')]
  }
 list(Yc=Y,Zc=Zc,Zi=Zi,sigma=sigma)
} # end function ImanConoverYcZc
################################################################################


################################################################################
IC_FactorYcZc = function(Y,f_rhos=0.5,seedval=NA) {
  nr=nrow(Y); nc=ncol(Y); nf_rhos=length(f_rhos)
  if(nf_rhos>1 & nf_rhos!=nc) stop('Dimensions of Y and f_rhos do not match') 
  if(nf_rhos==1) f_rhos = rep(f_rhos,nc)
  if(!is.na(seedval)) set.seed(seedval)
  Yc=Y
  zF = rnorm(nr)
  Zc = matrix(0,nr,nc)
  for(j in 1:nc){
   Zc[,j] = f_rhos[j] * zF  + sqrt(1-f_rhos[j]^2) * rnorm(nr)  
   ys=sort(Y[,j])
   Yc[,j] = ys[rank(Zc[,j],ties.method="first",na.last='keep')]
  } # end loop  for(j in 1:nc)
  list(Yc=Yc,Zc=Zc,zF=zF)
} # end function FactorYcZc
###############################################################################


###############################################################################
PDEF=function(sigma,Icor=F,eps=1E-5){
 tmp=eigen(sigma)
 v=tmp$vectors
 (d=tmp$values)
 d[d<=eps]=eps
 (d=diag(d))
 sigma = v %*% d %*% t(v)
 if(Icor==F) sigma=cov2cor(sigma)
 sigma
} # end function PDEF
################################################################################



################################################################################
my_require=function(pkglist){
 for (p in pkglist) {
  if(!require(p,character.only = TRUE)) {
    install.packages(p)
    require(p,character.only = TRUE)}
 }
} # end function my_require
################################################################################

###############################################################################
# my_hist function allows user to plot prob = counts/sum(counts) ot y axis
###############################################################################
my_hist=function(x,probs=F,pticks=5,plab='Prob',box=T,...){
 if(probs==F) hist(x,...)
 if(probs==T) {
  #tmp=hist(x,freq=T,yaxt='n',ylab='Prob')  # for debugging i.e. no ... dots
  tmp=hist(x,freq=T,yaxt='n',ylab=plab,...)
#   names(tmp)
#   [1] "breaks"   "counts"   "density"  "mids"     "xname"    "equidist"  
   (n=sum(tmp$counts))
   probs=round(tmp$counts/n,3)
   (at=seq(0,max(tmp$counts),length.out=pticks))
   (probs=seq(0,max(probs),length.out=pticks))
   axis(2,at=at,labels= probs)
 } # end if(probs==T) 
 if(box==T) box()  
} # end function my_hist
###############################################################################
# Examples with my_hist
# set.seed(1001)
# x=rnorm(10000)
# x11()
# my_hist(x,xlab='X',breaks=50)
# x11()
# my_hist(x,probs=T,xlab='X',breaks=50)
# x11()
# my_hist(x,probs=T,xlab='X',plab='PROBS',breaks=50)
# #############################################################################

###############################################################################
ddnorm=function(x,N=100000,text1='',main='',seedval=NA,Lplot='topright',
              Lcex=0.75,lwdx=6,lwdN=3,colx=1,colN=2){
############################################################################### 
 if(!is.na(seedval)) set.seed(seedval)
 n=max(length.na(x),N)  
 (mu=mean(x,na.rm=T))
 (sdev=sd(x,na.rm=T))
 x=(x-mu)/sdev
 xnorm=rnorm(N,0,1)
 tmp1=density(x); tmp2=density(xnorm)
 (xlim=range(c(tmp1$x,tmp2$x)))
 (ylim=range(c(tmp1$y,tmp2$y)))
 plot(tmp1,type='l',lwd=lwdx,main=main,xlim=xlim,ylim=ylim,xlab='',col=colx)
 points(tmp2,type='l',lwd=lwdN,col=colN)
 legend(Lplot,legend=c('normal'),fill=2,cex=Lcex)
###############################################################################
} # end function ddnorm
###############################################################################




################################################################################
#' merge_lists: Merges two list objects.
#' This function appends any objects in list2 and not in list 1
#' to list1 with priority given to list 1 components.
#'
#' @param  list1    A list object
#' @param  list2    A list object
#' @return list3    A merged object
merge_lists = function(list1, list2){
#####################################################################
# This function appends any objects in list2 and not in list 1
# to list1 with priority given to list 1 components.
#####################################################################
   (names1 = names(list1))
   (names2 = names(list2))

   (list3 = list1)
   (n1=length(list1))

   (pick2 = which(is.na(match(names2, names1))))

   if (length(pick2) > 0) {
     for(j in 1:length(pick2)) list3 = c(list3,list2[pick2[j]])
   }

   sort_list(list3)
} # end function merge_lists
#####################################################################
# # examples (not run)
# L1=list(a=1:3,b=4:6,c=7:9,e=10:12)
# L2=list(b=13:15,d=16:18,e=19:21,f=22:24)
# (L3=merge_lists(L1,L2))
# (L4=merge_lists(L2,L1))
###############################################################################

###############################################################################
#' sort_list: Sorts objects in list by object name
#'
#' @param  list1    A list object
#' @return list2    A resorted list object
sort_list=function(list1){
   (n=length(list1))
   df1=data.frame(Lnames=names(list1))
   df1$cell=1:n;  df1$cell=as.numeric(df1$cell)
   df1=df1[order(df1$Lnames),]
   list2=list1[df1$cell[1]]
   for(j in 2:n) list2=c(list2,list1[df1$cell[j]])
   list2
} # end function sort_list
################################################################################

################################################################################
betaN_calc=function(mu,sdv,k1=0,k2=1){
 ktest = min.na(k2 - k1)
 if (ktest <= 0)  stop("k2 must exceed k1")
 O1=(mu-k1)/sdv
 O2=(k2-mu)/sdv
 ahat=(O1/(O1+O2))*(O1*O2-1)
 bhat=(O2/(O1+O2))*(O1*O2-1)
 list(ahat=ahat,bhat=bhat)
}
########################################
betaN_est=function(x,k1=0,k2=1){
 ktest = min.na(k2 - k1)
 if (ktest <= 0)
     stop("k2 must exceed k1")
 if (min.na(x)<k1 | max(x)>k2)
     stop("The vector values must lie between k1 and k2")

 (mu=mean(x))
 (sdv=sd(x))
 O1=(mu-k1)/sdv
 O2=(k2-mu)/sdv
 ahat=(O1/(O1+O2))*(O1*O2-1)
 bhat=(O2/(O1+O2))*(O1*O2-1)
 list(ahat=ahat,bhat=bhat)
}
########################################
dbetaN=function(x,al=5,bt=5,k1=0,k2=1){
    ktest=min.na(k2-k1)
    if (ktest<=0)
      stop("k2 must exceed k1")
    atest=min.na(al)
    if (atest < 0)
      stop("The alpha parameter must exceed zero")
    btest=min.na(bt)
    if (btest < 0)
      stop("The beta parameter must exceed zero")

    xtest1=min.na(x-k1)
    xtest2=min.na(k2-x)

    if(xtest1 < 0 | xtest2 < 0)
      stop("The vector values must lie between k1 and k2")

 z=(x-k1)/(k2-k1)
 (1/(k2-k1))*dbeta(z,al,bt)
}
########################################
LLbetaN = function(x,al=5,bt=5,k1=0,k2=1){
 if(k1>=k2) stop('k2 must exceed k1')
 if(al<0) stop('The alpha parameter must exceed zero')
 if(bt<0) stop('The beta parameter must exceed zero')
 if(min(x)<k1|max(x)>k2) stop('The vector values must lie between k1 and k2')
 z=(x-k1)/(k2-k1)
 dbetaN=(1/(k2-k1))*dbeta(z,al,bt)
 dbetaN=dbetaN[dbetaN>0]
 sum(log(dbetaN))
}
########################################
pbetaN=function(x,al=5,bt=5,k1=0,k2=1){
    ktest=min.na(k2-k1)
    if (ktest<=0)
      stop("k2 must exceed k1")
    atest=min.na(al)
    if (atest < 0)
      stop("The alpha parameter must exceed zero")
    btest=min.na(bt)
    if (btest < 0)
      stop("The beta parameter must exceed zero")

    xtest1=min.na(x-k1)
    xtest2=min.na(k2-x)

    if(xtest1 < 0 | xtest2 < 0)
      stop("The vector values must lie between k1 and k2")

 z=(x-k1)/(k2-k1)
 pbeta(z,al,bt)
}
########################################
qbetaN=function(p,al=5,bt=5,k1=0,k2=1){
    ktest=min.na(k2-k1)
    if (ktest<=0)
      stop("k2 must exceed k1")
    atest=min.na(al)
    if (atest < 0)
      stop("The alpha parameter must exceed zero")
    btest=min.na(bt)
    if (btest < 0)
      stop("The beta parameter must exceed zero")

    ptest1=min.na(p)
    ptest2=max.na(p)

    if(ptest1 < 0 | ptest2 > 1)
      stop("The p vector values must satisfy 0<=p<=1")

  z=qbeta(p,al,bt)
 k1+z*(k2-k1)
}
########################################
rbetaN=function(n=100,al=5,bt=5,k1=0,k2=1){
    ktest=min.na(k2-k1)
    if (ktest<=0)
      stop("k2 must exceed k1")
    atest=min.na(al)
    if (atest < 0)
      stop("The alpha parameter must exceed zero")
    btest=min.na(bt)
    if (btest < 0)
      stop("The beta parameter must exceed zero")


 z=rbeta(n,al,bt)
 k1+z*(k2-k1)
}
########################################
# NOT Run
#x=sort(rbetaN(10000,2,2,0,10))
#graphics.off()
#tmp1=betaN_est(x=x,k1=0,k2=10)
#tmp2=betaN_est(x=x,k1=0,k2=20)
#plot(density(x),ylim=c(0,0.35))
#xsim1=rbetaN(10000,tmp1$ahat,tmp1$bhat,0,10)
#xsim2=rbetaN(10000,tmp2$ahat,tmp2$bhat,0,20)
#points(density(xsim1),col=2,type='l')
#points(density(xsim2),col=3,type='l')
#mean(x);sd(x)
#mean(xsim1);sd(xsim1)
#mean(xsim2);sd(xsim2)
################################################################################



################################################################################
bignine.na=function(df1,q1=0.10,q2=0.90,skipzero=F,Tpose=F){ # start function
  if(!is.data.frame(df1)) df1=as.data.frame(df1)
  Bigblock=NULL
  if(q1>1) q1=q1/100
  if(q2>1) q2=q2/100
  dnames=names(df1)
  (rnames=c('Min',paste('Q',q1*100,sep=''),'Quartile1','Median','Mean',
            'Stdev','Quartile2',
            paste('Q',q2*100,sep=''),'Max'))
  cnames=NULL
  (nc=ncol(df1))
  j=1
  for(j in 1:nc){                                # start loop j on columns
    x=df1[,j]
    if(is.numeric(x)){                            # start if on is.numeric
      (cnames=c(cnames,dnames[j]))
      x=x[!is.na(x)]
      if(skipzero==T) x=x[x!=0]
      (y=c(min(x),as.numeric(quantile(x,q1)),as.numeric(quantile(x,0.25)),
           median(x),mean(x),sd(x),as.numeric(quantile(x,0.75)),as.numeric(quantile(x,q2)),max(x)))
      Bigblock=cbind(Bigblock,y)
    }                                             # end loop on is.numeric
  }                                              # end loop j on columns
  
  df2=data.frame(Bigblock)
  rownames(df2)=rnames
  colnames(df2)=cnames
  
  
  if(Tpose==T){
    (rnames=colnames(df2))
    (cnames=rownames(df2))
    tmp=t(as.matrix(df2))
    df2=as.data.frame(tmp)
    colnames(df2)=cnames
    rownames(df2)=rnames
  }
  
  df2
}  # end function bignine.na
################################
bigeight.na=function(df1,q1=0.10,q2=0.90,skipzero=F,Tpose=F){ # start function
 if(!is.data.frame(df1)) df1=as.data.frame(df1)
 Bigblock=NULL
 if(q1>1) q1=q1/100
 if(q2>1) q2=q2/100
 dnames=names(df1)
 (rnames=c('Min',paste('Q',q1*100,sep=''),'Quartile1','Median','Mean','Quartile2',
   paste('Q',q2*100,sep=''),'Max'))
 cnames=NULL
 (nc=ncol(df1))
 j=1
 for(j in 1:nc){                                # start loop j on columns
  x=df1[,j]
  if(is.numeric(x)){                            # start if on is.numeric
   (cnames=c(cnames,dnames[j]))
   x=x[!is.na(x)]
   if(skipzero==T) x=x[x!=0]
   (y=c(min(x),as.numeric(quantile(x,q1)),as.numeric(quantile(x,0.25)),
   median(x),mean(x),as.numeric(quantile(x,0.75)),as.numeric(quantile(x,q2)),max(x)))
   Bigblock=cbind(Bigblock,y)
  }                                             # end loop on is.numeric
 }                                              # end loop j on columns

 df2=data.frame(Bigblock)
 rownames(df2)=rnames
 colnames(df2)=cnames


 if(Tpose==T){
  (rnames=colnames(df2))
  (cnames=rownames(df2))
   tmp=t(as.matrix(df2))
  df2=as.data.frame(tmp)
  colnames(df2)=cnames
  rownames(df2)=rnames
 }

 df2
}  # end function bigeight.na
########################################
big4.na = function(x,skipzero=F){
 xmin=min.na(x,skipzero=skipzero)
 xmean=mean.na(x,skipzero=skipzero)
 xmax=max.na(x,skipzero=skipzero)
 xmed=median.na(x,skipzero=skipzero)
 c(min=xmin,mean=xmean,median=xmed,max=xmax)
} # end function big4.na
#######################################
bigseven.na=function(df1,q1=0.25,q2=0.75,skipzero=F,Tpose=F){ # start function
 if(!is.data.frame(df1)) df1=as.data.frame(df1)
 Bigblock=NULL
 if(q1>1) q1=q1/100
 if(q2>1) q2=q2/100
 dnames=names(df1)
 (rnames=c('Min',paste('Q',q1*100,sep=''),'Median','Mean','Stdev',
 paste('Q',q2*100,sep=''),'Max'))
 cnames=NULL
 (nc=ncol(df1))
 j=1
 for(j in 1:nc){                                # start loop j on columns
  x=df1[,j]
  if(is.numeric(x)){                            # start if on is.numeric
   (cnames=c(cnames,dnames[j]))
   x=x[!is.na(x)]
   if(skipzero==T) x=x[x!=0]
   (y=c(min(x),as.numeric(quantile(x,q1)),median(x),mean(x),sd(x),
      as.numeric(quantile(x,q2)),max(x)))
   Bigblock=cbind(Bigblock,y)
  }                                             # end loop on is.numeric
 }                                              # end loop j on columns

 df2=data.frame(Bigblock)
 rownames(df2)=rnames
 colnames(df2)=cnames


 if(Tpose==T){
  (rnames=colnames(df2))
  (cnames=rownames(df2))
   tmp=t(as.matrix(df2))
  df2=as.data.frame(tmp)
  colnames(df2)=cnames
  rownames(df2)=rnames
 }

 df2

}   # end function bigseven.na
#######################################
bigsix.na=function(df1,q1=0.25,q2=0.75,skipzero=F,Tpose=F){ # start function
 if(!is.data.frame(df1)) df1=as.data.frame(df1)
 Bigblock=NULL
 if(q1>1) q1=q1/100
 if(q2>1) q2=q2/100
 dnames=names(df1)
 (rnames=c('Min',paste('Q',q1*100,sep=''),'Median','Mean',paste('Q',q2*100,sep=''),'Max'))
 cnames=NULL
 (nc=ncol(df1))
 j=1
 for(j in 1:nc){                                # start loop j on columns
  x=df1[,j]
  if(is.numeric(x)){                            # start if on is.numeric
   (cnames=c(cnames,dnames[j]))
   x=x[!is.na(x)]
   if(skipzero==T) x=x[x!=0]
   (y=c(min(x),as.numeric(quantile(x,q1)),median(x),mean(x),
       as.numeric(quantile(x,q2)),max(x)))
   Bigblock=cbind(Bigblock,y)
  }                                             # end loop on is.numeric
 }                                              # end loop j on columns

 df2=data.frame(Bigblock)
 rownames(df2)=rnames
 colnames(df2)=cnames


 if(Tpose==T){
  (rnames=colnames(df2))
  (cnames=rownames(df2))
   tmp=t(as.matrix(df2))
  df2=as.data.frame(tmp)
  colnames(df2)=cnames
  rownames(df2)=rnames
 }

 df2
################################################################################
}  # end function  bigsix.na
################################################################################
big3.na=function(x,skipzero=F,Round=3){
  xmin=round(min.na(x,skipzero=skipzero),Round)
  xmean=round(mean.na(x,skipzero=skipzero),Round)
  xmax=round(max.na(x,skipzero=skipzero),Round)
  c(mean=xmean,min=xmin,max=xmax)  
}# end function  big3.na 
################################################################################
class_df=function(df0){
 (ncols=ncol(df0))
 bigblock=NULL
 nc=1
 for(nc in 1:ncols){
  (cname=as.character(names(df0))[nc])
  (cclass=class(df0[,nc]))
  blockup=cbind(cname,cclass)
  bigblock=rbind(bigblock,blockup)
 }
 t(bigblock)
}                       # end function class_df

################################################################################

###############################################################################
#require(dplyr)
###############################################################################
#' A2SM: Convert a matrix A to sparse matrix form
#'
#' @param A        The matrix A.
#' @param SMM      Sparse matrix method with SMM='A',CRI','RCI','CMO', or 'RMO'
#' @param ZINDEX   T = Use zero indexing or F = Use one indexing.
#' @param eps      Value to use in non-zero test.
#' @param NAflag   Number to uses as NA flag
#'
#' @return A list object containing the matrix components:
#' @return nnz     = the number of non-zero elements in ra.
#' @return nr      = the number of rows in matrix A.
#' @return nc      = the number of columns in matrix A.
#' @return ia      = the row index.
#' @return ja      = the column index.
#' @return ra      = the non-zero coefficients in A
#' @return rnames  = matrix row names -- may be  ""
#' @return cnames  = matrix column names -- may be ""
#' @return SMM     = the sparse matrix type.
#' @return ZINDEX with T = Use zero indexing or F = Use one indexing.
#' @examples
#' \dontrun{
#' (A = matrix(c(1,0,0,2,0,3,0,0,0,5,0,6),3,4,byrow=T))
#' (SM1 = A2SM(A,SMM='CMO',ZINDEX=T))
#' (SM2=SM2SM(SM1,SMM='CRI',ZINDEX=F))
#' SM2A(SM1)
#' SM2A(SM2)
#'
#' (A = matrix(c(1,0,0,2,0,3,0,0,0,5,0,6),3,4,byrow=T)); A[2,3]=NA; A
#' (ASM = A2SM(A,SMM='CMO',ZINDEX=T)); SM2A(ASM)
#' (ASM2 = SM2SM(ASM,SMM2='CRI',ZINDEX2=T)); A ; (A2 = SM2A(SM2))
#'
#' #clpAPI  documentation example
#' nr=5
#' nc=8
#' ra=c(3.0,5.6,1.0,2.0,1.1,1.0,-2.0,2.8,-1.0,1.0,1.0,-1.2,-1.0,1.9)
#' ia=c(0,4,0,1,1,2,0,3,0,4,2,3,0,4)
#' ja=c(0,2,4,6,8,10,11,12,14)
#' SMM='CMO'
#' ZINDEX=T
#' ASM=list(nr=nr,nc=nc,ra=ra,ia=ia,ja=ja,SMM=SMM,ZINDEX=ZINDEX)
#' (A=SM2A(ASM))
#' }
###############################################################################
A2SM = function(A, SMM = 'CRI',ZINDEX = F,
                eps = .Machine$double.eps,NA_flag=(-1E6)){
##############################################################################
 my_require('dplyr')
 A[is.na(A)] = NA_flag
##############################################################################
  (nr = nrow(A))
  (nc = ncol(A))
  if(is.null(rownames(A))) rownames(A)=rep('',nr)
  if(is.null(colnames(A))) colnames(A)=rep('',nc)
  (rnames=rownames(A))
  (cnames=colnames(A))
##############################################################################
# jatwood algorithm 10/29/2022
# initially uses SMM RCI
  tmp=list(NULL); eps=1E-9
  j=1
  for(j in 1:nc){
   a=as.vector(A[,j])
   ipick=which(abs(a)>eps)
   ra=a[ipick]
   ia=ipick
   ja=rep(j,length(ipick))
   as.data.frame(cbind(ia,ja,ra))
   tmp[[j]]=as.data.frame(cbind(ia,ja,ra))
  } # end loop for(j in 1:nc)
  
  tmp=bind_rows(tmp)
  
  (nnz=nrow(tmp))

  ASM=list(nnz=nnz,nr=nr,nc=nc,ia=tmp$ia,ja=tmp$ja,ra=tmp$ra,
         rnames=rnames,cnames=cnames,SMM='CRI',ZINDEX=F)
  
  rm(tmp)
  
# convert to user specified sparse matrix format     
  ASM=SM2SM(ASM,SMM2=SMM,ZINDEX2=ZINDEX)

  return(ASM)
} # end function A2SM
###############################################################################


###############################################################################
#' SM2A: Convert a sparse matrix object ASM into a matrix A
#'
#' ASM is a sparse matrix object containing:
#' nr = number of rows, nc = number of columns, ra = nonzero coeff,
#' ia = row indices, ja = col indices,
#' rnames = row names,cnames = column names # Note these may be missing or ""
#' SMM = space matrix method, and ZINDEX = use 0-indexing
#' @param ASM  A sparse matrix object
#' @param NAflag   Number to uses as NA flag
#'
#' @return The matrix A.
#' @examples
#' \dontrun{
#' (A = matrix(c(1,0,0,2,0,3,0,0,0,5,0,6),3,4,byrow=T))
#' (SM1 = A2SM(A,SMM='CMO',ZINDEX=T))
#' (SM2=SM2SM(SM1,SMM='CRI',ZINDEX=F))
#' SM2A(SM1)
#' SM2A(SM2)
#'
#' (A = matrix(c(1,0,0,2,0,3,0,0,0,5,0,6),3,4,byrow=T)); A[2,3]=NA; A
#' (ASM = A2SM(A,SMM='CMO',ZINDEX=T)); SM2A(ASM)
#' (ASM2 = SM2SM(ASM,SMM2='CRI',ZINDEX2=T)); A ; (A2 = SM2A(SM2))
#'
#' #clpAPI  documentation example
#' nr=5
#' nc=8
#' ra=c(3.0,5.6,1.0,2.0,1.1,1.0,-2.0,2.8,-1.0,1.0,1.0,-1.2,-1.0,1.9)
#' ia=c(0,4,0,1,1,2,0,3,0,4,2,3,0,4)
#' ja=c(0,2,4,6,8,10,11,12,14)
#' SMM='CMO'
#' ZINDEX=T
#' ASM=list(nr=nr,nc=nc,ra=ra,ia=ia,ja=ja,SMM=SMM,ZINDEX=ZINDEX)
#' (A=SM2A(ASM))
#' }
###############################################################################
SM2A = function(ASM,NA_flag=(-1E6)){
###############################################################################
 ASM$ra[is.na(ASM$ra)] = NA_flag
##
  nr=ASM$nr; nc=ASM$nc
  if(!exists('rnames',where=ASM)) ASM$rnames=rep('',nr)
  if(!exists('cnames',where=ASM)) ASM$cnames=rep('',nc)
##
 (SMM0=ASM$SMM)

 if(SMM0=='A') A = ASM$A

 if(SMM0!='A'){

  ASM2=SM2SM(ASM,SMM2='CRI',ZINDEX2=F)
  nr=ASM2$nr;nc=ASM2$nc;
  ia=ASM2$ia;ja=ASM2$ja;ra=ASM2$ra


  (ra=as.numeric(ra))
  (ia=as.integer(ia))
  (ja=as.integer(ja))

  (nr=as.integer(nr))
  (nc=as.integer(nc))

  A=matrix(0,nr,nc)
  A[cbind(ia,ja)]=ra

 } # end if(SMM0!='A')
##
  A[A <= NA_flag] = NA
##
  if(!is.null(ASM$rnames)&length(ASM$rnames)==nrow(A)) rownames(A)=ASM$rnames
  if(!is.null(ASM$cnames)&length(ASM$cnames)==ncol(A)) colnames(A)=ASM$cnames
##
  return(A)
} # end function SM2A
###############################################################################


###############################################################################
#' SM2SM: Convert a sparse matrix form into a different sparse matrix form.
#'
#' ASM is a sparse matrix object containing:
#' nr = number of rows, nc = number of columns, ra = nonzero coeff,
#'  ia = row indices, ja = col indices, SMM = space matrix method,
#'  ZINDEX=use 0-indexing
#' @param ASM   A sparse matrix object
#' @param SMM2  sparse matrix method with SMM2 ='A','CRI','RCI','CMO',or'RMO'
#' @param ZINDEX2   T = Use zero indexing or F = Use one indexing.
#' @param NAflag   Number to uses as NA flag
#'
#' @return A list object containing the following sparse matrix components:
#' @return nnz    = the number of non-zero elements in ra.
#' @return nr     = the number of rows in matrix.
#' @return nc     = the number of columns in matrix.
#' @return ia     = the revised row index.
#' @return ja     = the revised column index.
#' @return rnames = row names  -- may be missing or ""
#' @return cnames = column names -- may be missing or ""
#' @return ra     = the revised non-zero coefficients in A
#' @return SMM    = the revised sparse matrix type.
#' @return ZINDEX = T or F for the revised sparse form.
#' @examples
#' \dontrun{
#' (A = matrix(c(1,0,0,2,0,3,0,0,0,5,0,6),3,4,byrow=T))
#' (SM1 = A2SM(A,SMM='CMO',ZINDEX=T))
#' (SM2=SM2SM(SM1,SMM='CRI',ZINDEX=F))
#' SM2A(SM1)
#' SM2A(SM2)
#'
#' (A = matrix(c(1,0,0,2,0,3,0,0,0,5,0,6),3,4,byrow=T)); A[2,3]=NA; A
#' (ASM = A2SM(A,SMM='CMO',ZINDEX=T)); SM2A(ASM)
#' (ASM2 = SM2SM(ASM,SMM2='CRI',ZINDEX2=T)); A ; (A2 = SM2A(SM2))
#'
#' #clpAPI  documentation example
#' nr=5
#' nc=8
#' ra=c(3.0,5.6,1.0,2.0,1.1,1.0,-2.0,2.8,-1.0,1.0,1.0,-1.2,-1.0,1.9)
#' ia=c(0,4,0,1,1,2,0,3,0,4,2,3,0,4)
#' ja=c(0,2,4,6,8,10,11,12,14)
#' SMM='CMO'
#' ZINDEX=T
#' ASM=list(nr=nr,nc=nc,ra=ra,ia=ia,ja=ja,SMM=SMM,ZINDEX=ZINDEX)
#' (A=SM2A(ASM))
#' }
###############################################################################
SM2SM = function(ASM, SMM2='CRI', ZINDEX2=F, NA_flag=(-1E6)){
###############################################################################
  ASM$ra[is.na(ASM$ra)] = NA_flag
##
  nnz=ASM$nnz;nr=ASM$nr; nc=ASM$nc
  if(!exists('rnames',where=ASM)) ASM$rnames=rep('',nr)
  if(!exists('cnames',where=ASM)) ASM$cnames=rep('',nc)
##
  (rnames=ASM$rnames)
  (cnames=ASM$cnames)
## original sparse form
  (SMM0=ASM$SMM)
  
## if SMM0 indicates a matrix convert ASM to SMM='CRI' sparse form
  if(SMM0=='A') {
    (tmp=A2SM(ASM$A, SMM = 'CRI', ZINDEX = F))
    ASM=merge_lists(tmp,ASM)
  } # end if(SMM0=='A')

## Pull sparse matrix elements from ASM 
  (nr=ASM$nr);(nc=ASM$nc);(ia=ASM$ia);(ja=ASM$ja);(ra=ASM$ra)
  (SMM=ASM$SMM);(ZINDEX=ASM$ZINDEX)

## convert "zero" to "one" indexing if needed
  if(ZINDEX==T){
    ia=ia+1
    ja=ja+1
  }# end if (ZINDEX == T)

  ia;ja

## convert CMO ja vector or RMO ia vector to complete index vectors     
  if(SMM=='CMO') ja=rep(1:nc,diff(ja))
  if(SMM=='RMO') ia=rep(1:nr,diff(ia))
  
  tmp=as.data.frame(cbind(ia,ja,ra))
   
## convert ASM to ASM2 form 
  if(SMM2=='CRI'|SMM2=='CMO'){
    # resort to be sure in "complete and contiguous CRI form
    (tmp=tmp[with(tmp,order(ja,ia)),])
    # pull resorted data
    ra=tmp$ra; ia=tmp$ia; ja=tmp$ja
     
    if(SMM2=='CMO'){ 
     (ja=c(match(1:nc,ja),nnz+1))
    } # end if(SMM2=='CMO')     
         
  } # end if(SMM2=='CRI'|SMM2=='CMO')
  
  if(SMM2=='RCI'|SMM2=='RMO'){
    # resort to be sure in "complete and contiguous RCI form
    (tmp=tmp[with(tmp,order(ia,ja)),])
    # pull resorted data
    ra=tmp$ra; ia=tmp$ia; ja=tmp$ja
     
    if(SMM2=='RMO'){ 
     ia=c(match(1:nr,ia),nnz+1)   
    } # end if(SMM2=='RMO')     
    
  } # end if(SMM2=='RCI'|SMM2=='RMO)

  (ra=as.numeric(ra))
  (ia=as.integer(ia))
  (ja=as.integer(ja))
  (nr=as.integer(nr))
  (nc=as.integer(nc))

  if(ZINDEX2 == T) {ia = ia-1; ja = ja-1}
##
##
ASM2=ASM # inherits all objects in ASM
# put revised data into ASM2
ASM2$nnz=nnz;ASM2$nr=nr;ASM2$nc=nc;ASM2$ia=ia;ASM2$ja=ja;ASM2$ra=ra
ASM2$SMM=SMM2;ASM2$ZINDEX=ZINDEX2


if(SMM2=='A'){
 ASM2=ASM  # inherits all objects in  original ASM
 ASM2$A = SM2A(ASM)
 ASM2$SMM = SMM2
}
#(ASM2 = list(nnz=nnz,nr = nr,nc = nc,ia = ia2,ja = ja2,ra = ra2,
#            rnames=rnames,cnames=cnames,SMM = SMM2,ZINDEX = ZINDEX2))
###############################################################################
return(ASM2)
###############################################################################
} # end function SM2SM
###############################################################################


# ###############################################################################
# # Examples:
# ###############################################################################
# # create matrix with missing data 
# (A = matrix(c(1,0,0,2,0,3,0,0,0,5,0,6),3,4,byrow=T)); A[2,3]=NA; A
# # no names for in matrix A
# c(rownames(A),colnames(A))
# # Convert to sparse matrix with NA's replaced by NA_flag
# (SM1=A2SM(A,NA_flag=-1000))
# # add names to A matrix
# rownames(A)=paste('r',1:nrow(A),sep='');colnames(A)=paste('c',1:ncol(A),sep='')
# A; SM2=A2SM(A)
# # Note below: SM2A puts " NA's"  back in matrix if original matrix had them
# # of if user indicates levels in matrix that designate missing data
# SM2A(SM1);SM2A(SM2)
# #####################
# ## test heritability of SM objects
# SM4$M = matrix(1:6,2,3)
# SM4
# (SM5 = SM2SM(SM4,SMM2='CMO'))
# #####################
# (ASM = A2SM(A,SMM='CMO',ZINDEX=T)); SM2A(ASM)
# (ASM2 = SM2SM(ASM,SMM2='CRI',ZINDEX2=T)); A ; (A2 = SM2A(SM2))
# ##############################################################################
# 
# ##############################################################################
# A=t(matrix(1:15,5,3))
# A[cbind(c(3,1,2,3,1,3),c(1,2,3,3,4,5))]=0
# A
# #
# # (A=matrix(c(1,0,2,0,0,3,0,4,0,5,6,0),3,4,byrow=T))
# # (A=matrix(c(1,0,0,2,0,3,4,0,0,5,0,6),3,4,byrow=T))
# # (A=matrix(c(1,0,0,2,0,3,0,0,0,5,0,6),3,4,byrow=T))
# # (A=matrix(c(2,3,1,4,1,2,3,4,2),3,3,byrow=T))
# #
# A2SM(A,SMM='CMO',ZINDEX=T)
# 
# SMMlist=c('RCI','CRI','RMO','CMO')
# ZINDEXlist=c(F,T)
# 
# SMM1 = 'RCI'     # 'RCI', 'CRI', 'RMO', 'CMO'
# ZINDEX1 = F      # T or F
# SMM2='CRI'       # 'RCI', 'CRI', 'RMO', 'CMO'
# ZINDEX2 = F      # T or F
# 
# for(SMM1 in SMMlist){
#  for(ZINDEX1 in ZINDEXlist) {
#   for(SMM2 in SMMlist) {
#    for(ZINDEX2 in ZINDEXlist) {
#     print(paste(SMM1,ZINDEX1,SMM2,ZINDEX2))
#     (ASM=A2SM(A,SMM=SMM1,ZINDEX=ZINDEX1))
#     (ASM1=SM2SM(ASM,SMM2=SMM2,ZINDEX2=ZINDEX2))
#     (ASM2=A2SM(A,SMM=SMM2,ZINDEX=ZINDEX2))
#     print(summary(unlist(ASM1)==unlist(ASM2)))
#    } # end for(ZINDEX2 in ZINDEXlist)
#   } # end for(SMM2 in SMMlist)
#  } # end for(ZINDEX1 in ZINDEXlist)
# } # end for(SMM1 in SMMlist)
# 
# 
# summary(SM2A(ASM1)==SM2A(ASM2))
# #############################################################
# (A=matrix(c(1,0,2,0,0,3,0,4,0,5,6,0),3,4,byrow=T))
# 
# (ASM=A2SM(A,SMM='A'))
# ASM2=SM2SM(ASM,SMM2=SMM2,ZINDEX2=ZINDEX2)
# (A2=SM2A(ASM2))
# summary(as.vector(A)==as.vector(A2))
# 
#  for(SMM2 in SMMlist) {
#   for(ZINDEX2 in ZINDEXlist) {
#    print(paste(SMM2,ZINDEX2))
#    (ASM=A2SM(A=A,SMM='A'))
#    (ASM2=SM2SM(ASM,SMM2=SMM2,ZINDEX2=ZINDEX2))
#    (A2=SM2A(ASM2))
#    print(summary(as.vector(A)==as.vector(A2)))
#   } # end for(ZINDEX2 in ZINDEXlist)
#  } # end for(SMM2 in SMMlist)
# ###########################################################
# (A=matrix(c(1,0,2,0,0,3,0,4,0,5,6,0),3,4,byrow=T))
# SMM2='CRI'      #  'RCI', 'CRI', 'RMO', 'CMO'
# ZINDEX2 = F     # T or F
# 
# (ASM=A2SM(A,SMM='A'))
# ASM2=SM2SM(ASM,SMM2=SMM2,ZINDEX2=ZINDEX2)
# (A2=SM2A(ASM2))
# summary(as.vector(A)==as.vector(A2))
# (A3=SM2A(ASM))
# summary(as.vector(A)==as.vector(A3))
# ##########################################################
# # clpAPI  documentation example
# nr=5
# nc=8
# ra=c(3.0,5.6,1.0,2.0,1.1,1.0,-2.0,2.8,-1.0,1.0,1.0,-1.2,-1.0,1.9)
# ia=c(0,4,0,1,1,2,0,3,0,4,2,3,0,4)
# ja=c(0,2,4,6,8,10,11,12,14)
# SMM='CMO'
# ZINDEX=T
# ASM=list(nr=nr,nc=nc,ra=ra,ia=ia,ja=ja,SMM=SMM,ZINDEX=ZINDEX)
# (A=SM2A(ASM))
# #############################################################################





###############################################################################
#' cbindSM: "column bind" two sparse matrices
#'
#' @param   SM1     First sparse matrix object
#' @param   SM2     Second sparse matrix object
#' @param   SMM     Sparse matrix method with 'CRI','RCI','CMO',or'RMO'
#' @param   ZINDEX  T = Use zero indexing or F = Use one indexing.
#'
#' @return A 'column bound" sparse matrix object with components"
#' @return nnz     = the number of non-zero elements in ra.
#' @return nr      = the number of rows in matrix A.
#' @return nc      = the number of columns in matrix A.
#' @return ia      = the row index.
#' @return ja      = the column index.
#' @return ra      = the non-zero coefficients in A
#' @return rnames  = matrix row names -- may be ''
#' @return cnames  = matrix column names -- may be ''
#' @return SMM     = the sparse matrix type.
#' @return ZINDEX with T = Use zero indexing or F = Use one indexing.
###############################################################################
cbindSM=function(SM1,SM2,SMM='CRI',ZINDEX=F){
 nr1=SM1$nr; nc1=SM1$nc; nr2=SM2$nr; nc2=SM2$nc
 if(nr1!=nr2) stop('STOP: row numbers not compatible')
 tmp1=SM2SM(SM1,SMM='CRI',ZINDEX=F)
 tmp2=SM2SM(SM2,SMM='CRI',ZINDEX=F)
 nr=nr1; nc = nc1+nc2
 ia=c(tmp1$ia,tmp2$ia)
 ja=c(tmp1$ja,tmp2$ja+nc1)
 ra=c(tmp1$ra,tmp2$ra)
 tmp1$nnz=length(ra);tmp1$nr=nr;tmp1$nc=nc
 tmp1$ia=ia;tmp1$ja=ja;tmp1$ra=ra;tmp1$SMM='CRI';tmp1$ZINDEX=F

 if(length(tmp1$cnames)>1) {
   if(length(tmp2$cnames)>1)  tmp1$cnames=c(tmp1$cnames,tmp2$cnames)
   if(length(tmp2$cnames)<=1) tmp1$cnames=c(tmp1$cnames,rep('',nc2))
 } # end

 SM2SM(tmp1,SMM=SMM,ZINDEX=ZINDEX)
} # end function cbindSPM
###############################################################################

###############################################################################
#' rbindSM: "row bind" two sparse matrices
#'
#' @param   SM1     First sparse matrix object
#' @param   SM2     Second sparse matrix object
#' @param   SMM     Sparse matrix method with 'CRI','RCI','CMO',or'RMO'
#' @param   ZINDEX  T = Use zero indexing or F = Use one indexing.
#'
#' @return  A 'row bound" sparse matrix object with components"
#' @return nnz     = the number of non-zero elements in ra.
#' @return nr      = the number of rows in matrix A.
#' @return nc      = the number of columns in matrix A.
#' @return ia      = the row index.
#' @return ja      = the column index.
#' @return ra      = the non-zero coefficients in A
#' @return rnames  = matrix row names -- may be ''
#' @return cnames  = matrix column names -- may be ''
#' @return SMM     = the sparse matrix type.
#' @return ZINDEX with T = Use zero indexing or F = Use one indexing.
###############################################################################
rbindSM=function(SM1,SM2,SMM='CRI',ZINDEX=F){
 nr1=SM1$nr; nc1=SM1$nc; nr2=SM2$nr; nc2=SM2$nc
 if(nc1!=nc2) stop('STOP: column number not compatible')
 tmp1=SM2SM(SM1,SMM='RCI',ZINDEX=F)
 tmp2=SM2SM(SM2,SMM='RCI',ZINDEX=F)
 nr=nr1+nr2; nc = nc1
 ia=c(tmp1$ia,tmp2$ia+nr1)
 ja=c(tmp1$ja,tmp2$ja)
 ra=c(tmp1$ra,tmp2$ra)
 tmp1$nnz=length(ra);tmp1$nr=nr;tmp1$nc=nc
 tmp1$ia=ia;tmp1$ja=ja;tmp1$ra=ra;tmp1$SMM='RCI';tmp1$ZINDEX=F

 if(length(tmp1$rnames)>1) {
   if(length(tmp2$rnames)>1)  tmp1$rnames=c(tmp1$rnames,tmp2$rnames)
   if(length(tmp2$rnames)<=1) tmp1$rnames=c(tmp1$rnames,rep('',nr2))
 } # end

 SM2SM(tmp1,SMM=SMM,ZINDEX=ZINDEX)
} # end function rbindSPM
###############################################################################
#
#
#
################################################################################
# # Examples of cbindSM and rbindSM
#  set.seed(1001)
#  nr=7;nc=7
#  (A=matrix(round(runif(nr*nc),3),nr,nc))
#  A[A<0.5]=0
#  colnames(A)=paste('c',1:nc,sep='')
#  rownames(A)=paste('r',1:nr,sep='')
#  A
#
#  (A11=A[1:5,1:5]); (A12=A[1:5,6:7])
#  (A21=A[6:7,1:5]); (A22=A[6:7,6:7])
#
#
#  (SM11=A2SM(A11))
#  (SM12=A2SM(A12))
#  (SM21=A2SM(A21))
#  (SM22=A2SM(A22))
#
#  cbind(A11,A12)
#  SM2A(cbindSM(SM11,SM12))
#
#
#  rbind(A11,A21)
#  SM2A(rbindSM(SM11,SM21))
#
#  A
#  SM2A(cbindSM(rbindSM(SM11,SM21),rbindSM(SM12,SM22)))
###############################################################################
#  A2SM(A,SMM='CMO')
#  cbindSM(rbindSM(SM11,SM21),rbindSM(SM12,SM22),SMM='CMO')
###############################################################################











###############################################################################
#' merge_lists: Merges two list objects.
#' This function appends any objects in list2 and not in list 1
#' to list1 with priority given to list 1 components.
#'
#' @param  list1    A list object
#' @param  list2    A list object
#' @return list3    A merged object
#' @examples
#' \dontrun{
#' L1=list(a=1:3,b=4:6,c=7:9,e=10:12)
#' L2=list(b=13:15,d=16:18,e=19:21,f=22:24)
#' (L3=merge_lists(L1,L2))
#' (L4=merge_lists(L2,L1))
#' }
###############################################################################
merge_lists = function(list1, list2){
#####################################################################
# This function appends any objects in list2 and not in list 1
# to list1 with priority given to list 1 components.
#####################################################################
   (names1 = names(list1))
   (names2 = names(list2))

   (list3 = list1)
   (n1=length(list1))

   (pick2 = which(is.na(match(names2, names1))))

   if (length(pick2) > 0) {
     for(j in 1:length(pick2)) list3 = c(list3,list2[pick2[j]])
   }

   sort_list(list3)
} # end function merge_lists
###############################################################################

###############################################################################
#' sort_list: Sorts objects in list by object name
#'
#' @param  list1    A list object
#' @return list2    A resorted list object
###############################################################################
sort_list=function(list1){
   (n=length(list1))
   df1=data.frame(Lnames=names(list1))
   df1$cell=1:n;  df1$cell=as.numeric(df1$cell)
   df1=df1[order(df1$Lnames),]
   list2=list1[df1$cell[1]]
   for(j in 2:n) list2=c(list2,list1[df1$cell[j]])
   list2
} # end function sort_list
################################################################################






################################################################################
jday_shift=function(jdays,baseyr,newyr){
 #require(chron)
 tmp=data.frame(month.day.year(jdays))
 tmp$year=tmp$year+(newyr-baseyr)
 julian(tmp$month,tmp$day,tmp$year)
}
############################################################
jmonth=function(mo=1,yr=2013){
 12*(yr-1)+mo
}
############################################################
mo.yr=function(jmo){
 yr=as.integer(jmo/12)+1
 mo=jmo-(yr-1)*12
 mo=ifelse(mo==0,12,mo)
 yr=ifelse(mo==12,yr-1,yr)
 data.frame(cbind(mo,yr))
}
############################################################
#set.seed(1001)
#nobs=10
#(mo=sample(1:12,nobs,replace=T))
#(yr=sample(1900:2013,nobs,replace=T))
#
#(jmo=jmonth(mo,yr))
#mo.yr(jmo)
#cbind(mo,yr)
############################################################

############################################################
jbimonth=function(mo=1,day=1,yr=2013,bimocut=15){
 jbimo=jmonth(mo,yr)*2-1
 jbimo=ifelse(day>=bimocut,jbimo+1,jbimo)
 jbimo
}
###########################################################
bimo.yr=function(jbimo,bimocut=15){
 (jmo=as.integer((jbimo+1)/2))
 (moyr=mo.yr(jmo))
 (mo=moyr$mo)
 (yr=moyr$yr)
 (daytest=((jbimo+1)/2)%%1)
 (day=ifelse(daytest<0.001,1,bimocut))
 list(mo=mo,day=day,yr=yr)
}
###########################################################
#(jbimo=jbimonth(2,14,1))
#bimo.yr(jbimo)
#data.frame(bimo.yr(1:12))

#jbimo=jbimonth(1,1,2013)
#(jbimo=jbimo:(jbimo+23))
#data.frame(bimo.yr(jbimo))
################################################################################


################################################################################
# Johnson distribution
################################################################################
JBpars=function(x,k1=0,k2=1){
  jtest=min.na(k2-k1)
  if (jtest<=0)  stop("k2 must exceed k1")
    xtest1=min.na(x-k1)
    xtest2=min.na(k2-x)
    if(xtest1 < 0 | xtest2 < 0) stop("The vector values must lie between k1 and k2")
 z=(x-k1)/(k2-x)
 u=log(z)
list(Jmn=mean(u),Jsd=sd(u))
}
###############################################
rJB=function(n,Jmn,Jsd,k1=0,k2=1){
  jtest=min.na(k2-k1)
  if (jtest<=0)  stop("k2 must exceed k1")
 u=rnorm(n,Jmn,Jsd)
 z=exp(u)
 (x=(k1+k2*z)/(1+z))
}
###############################################
pJB=function(q,Jmn,Jsd,k1=0,k2=1){
  jtest=min.na(k2-k1)
  if (jtest<=0)  stop("k2 must exceed k1")
  qtest1=min.na(q-k1)
  qtest2=min.na(k2-q)
  if(qtest1 < 0 | qtest2 < 0) stop("The vector values must lie between k1 and k2")
 z=(q-k1)/(k2-q)
 u=log(z)
 pnorm(u,Jmn,Jsd)
}
###############################################
qJB=function(p,Jmn,Jsd,k1=0,k2=1) {
  jtest=min.na(k2-k1)
  if (jtest<=0)  stop("k2 must exceed k1")
  ptest1=min.na(p)
  ptest2=min.na(1-p)
  if(ptest1<0|ptest2<0) stop("The vector values must lie between k1 and k2")

 u=qnorm(p,Jmn,Jsd)
 z=exp(u)
 (x=(k1+k2*z)/(1+z))
}
#################################################
dJB=function(x,Jmn,Jsd,k1=0,k2=1){
  jtest=min.na(k2-k1)
  if (jtest<=0)  stop("k2 must exceed k1")
  xtest1=min.na(x-k1)
  xtest2=min.na(k2-x)
  if(xtest1 < 0 | xtest2 < 0) stop("The vector values must lie between k1 and k2")

 xi=k1
 lambda=k2-k1
 delta=1/Jsd
 gamma=-delta*Jmn
# parms=list(gamma=gamma,delta=delta,xi=xi,lambda=lambda,type='SB')
# dJohnson(x,parms=joe1)
 z=(x-xi)/lambda
(delta/(lambda*sqrt(2*pi)*z*(1-z))) *exp(-0.5*(gamma+delta*log(z/(1-z)))^2)
}
################################################################################


################################################################################
longlatdistbear = function(long0,lat0,distance,bear,rads=F){

  if(rads!=T) {
  longrad0=long0/57.2958;latrad0=lat0/57.2958
  bearR=bear/57.2958}

  if(rads==T) {
  longrad0=long0;latrad0=lat0
  bearR=bear}

  (D=distance/3963)
  sinlatrad1=cos(bearR)*cos(latrad0)*sin(D)+sin(latrad0)*cos(D)
  latrad1=asin(sinlatrad1)
  calc4=(cos(D)-sin(latrad0)*sin(latrad1))/(cos(latrad0)*cos(latrad1))
  calc5=acos(calc4)
  longrad11=longrad0-calc5
  calc6=(-1*acos(calc4))
  longrad12=longrad0-calc6

  bearR11=bearcalc(longrad0,latrad0,longrad11,latrad1,rads=T)
  bearR12=bearcalc(longrad0,latrad0,longrad12,latrad1,rads=T)
  bearErr11=abs(bearR11-bearR)
  longrad1=ifelse(bearErr11<0.1,longrad11,longrad12)

  if(rads!=T){
  long1=longrad1*(360/(2*pi))
  long1=ifelse(long1<(-180),180+(long1+180),long1)
  lat1=latrad1*(360/(2*pi))
  }

  if(rads==T){
  long1=long1
  lat1=lat1
  }


 dftmp=data.frame(long1,lat1)
 dftmp
} # end function longlatdistbear
########################################
distbearcalc = function(long0,lat0,long1,lat1,rads=F){
  if(rads!=T) {
  longrad0=long0/57.2958;latrad0=lat0/57.2958
  longrad1=long1/57.2958;latrad1=lat1/57.2958}
  if(rads==T) {
  longrad0=long0;latrad0=lat0
  longrad1=long1;latrad1=lat1}

  cosD=(sin(latrad0)*sin(latrad1)+cos(latrad0)*cos(latrad1)*cos(longrad0-longrad1))
  D=acos(cosD)
  distance=3963*D
  distance=ifelse(distance=="NaN",0,distance)
  distance
 # calc bearing
   calc3=(sin(latrad1)-(sin(latrad0)*cos(D)))/(cos(latrad0)*sin(D))
   (signtest=sin(longrad0-longrad1))
   bearR=acos(calc3)
   bearR=ifelse(signtest=="NaN",0,bearR)
   bearR=ifelse(signtest>0,bearR,2*pi-bearR)
   bearR

if(rads!=T){bear=bearR*(360/(2*pi))}
if(rads==T){bear=bearR}
bear
return(list(distance=distance,bear=bear))
}
#######################################
distcalc =function(long0,lat0,long1,lat1,rads=F){
  if(rads!=T) {
  longrad0=long0/57.2958;latrad0=lat0/57.2958
  longrad1=long1/57.2958;latrad1=lat1/57.2958}
  if(rads==T) {
  longrad0=long0;latrad0=lat0
  longrad1=long1;latrad1=lat1}

  cosD=(sin(latrad0)*sin(latrad1)+cos(latrad0)*cos(latrad1)*cos(longrad0-longrad1))
  D=acos(cosD)
  distance=3963*D
  distance=ifelse(distance=="NaN",0,distance)
  distance
}

#######################################
## Example only
#load("CtyMaps50.RData")
#us=CtyData[,c('sfips','fips','long','lat')]
#mt=subset(us,sfips==30)
#mt=mt[order(mt$fips),]
#x=mt$long
#y=mt$lat
#zerodist=0.1
#ids=mt$fips
#
#dftmp=distmat(x,y,ids=mt$fips)
#dstack=dftmp$dstack
#dstack$fips=dstack$locto
#
#mt2=merge(mt,dstack)
#mt3=subset(mt2,locfrom==30001)
#library(maps)
#map('county','mont')
#text(mt3$long,mt3$lat,round(mt3$dist),font=2,cex=0.75)
################################################################################
distmat=function(x,y,ids=1:length(x),zerodist=0.1){

 if(length(x)!=length(y)) stop('Error - lengths of x and y not equal')
 nobs=length(x)
 locto=ids

 dstack=NULL
 dmat=NULL

 i=1
 for(i in 1:nobs){                              # start loop on i
  locfrom=ids[i]
  dftmp3=data.frame(locfrom,locto)
  dftmp3$dist=distcalc(x[i],y[i],x,y)
  dftmp3$dist=ifelse(dftmp3$dist<zerodist,0,dftmp3$dist)

  dstack=rbind(dstack,dftmp3)
  dmat=rbind(dmat,t(dftmp3$dist))
 }                                              # end loop on i

 dstack=data.frame(dstack)
 dmat=as.matrix(dmat)

 return(list(dstack=dstack,dmat=dmat))

}  # end function distmat
################################################################################


###########################################################################
lookup=function(x,getvals,xmax=10000000){
 x1=c(getvals[,1],xmax)
 getvals=as.data.frame(getvals)
 nobs=length(x1)
 pickum=as.numeric(cut(x,x1,include.lowest = TRUE,right=F))
 getvals[pickum,]
}
############################################################################
# Not Run
#cvg=seq(0.5,0.85,0.05)
#sub=c(0.67,0.64,0.64,0.59,0.59,0.55,0.48,0.38)
#df0=data.frame(cbind(cvg,sub))
#x=seq(0.49,0.86,0.01)
#cbind(x,lookup(x,getvals=df0))
##############################################################################


##############################################################################
maxgdhess = function(obj,b,itlim=100,epstop=0.0000001,plotum=F,
ptitle='',ndirs=1,useDIRg=F,...)   {
 #library(MASS)
 #library(numDeriv)
 dirgen=function(w,d1,d2){(1-w)*d1+w*d2}
  w=as.matrix(seq(0,1,1/ndirs))
 npars=length(b)
 one.w=matrix(1,(length(w)),1)
 if(useDIRg==T) one.w=matrix(1,(length(w)+npars),1)
 stopum=0
 iter=1
 (fval0=obj(b,...))
# (fval0=obj(b))
 fhist=fval0

 while(stopum==0&iter<=itlim) {         # first while
  (fval0=obj(b,...))
#  (fval0=obj(b))
  (gradb=grad(obj,b,...))
#  (gradb=grad(obj,b))
  (gradtest=as.numeric(sum(abs(gradb))))

  (hess=hessian(obj,b,...))
#  (hess=hessian(obj,b))

  if(npars==1) {hessnd=-1*abs(hess);hesstest=hess}
  if(npars>=2) {
   (ehess=eigen(hess))
   (v=ehess$vectors)
   (d=diag(ehess$values))
   (hesstest=max(d))
   (d=-1*abs(d)  )
   (hessnd=v %*% d %*% t(v))
   }
  (if(gradtest<=epstop&hesstest<0){stopum=1})
  d1=-1* ginv(hessnd) %*% gradb
  d2=gradb
  (DIR=apply(w,1,'dirgen',d1=d1,d2=d2))  # generates an npar by ndirs+1 matrix of directions
  (fval0=obj(b,...))
#  (fval0=obj(b))
  (gradb=grad(obj,b,...))
#  (gradb=grad(obj,b))
  (gradtest=as.numeric(sum(abs(gradb))))

  (hess=hessian(obj,b,...))
#  (hess=hessian(obj,b))

  if(npars==1) {hessnd=-1*abs(hess);hesstest=hess}
  if(npars>=2) {
   (ehess=eigen(hess))
   (v=ehess$vectors)
   (d=diag(ehess$values))
   (hesstest=max(d))
   (d=-1*abs(d)  )
   (hessnd=v %*% d %*% t(v))
   }
  (if(gradtest<=epstop&hesstest<0){stopum=1})
  d1=-1* ginv(hessnd) %*% gradb
  d2=gradb
  (DIR=apply(w,1,'dirgen',d1=d1,d2=d2))  # generates  matrix of directions

   DIRg=diag(gradb)
   if(useDIRg==T)  DIR=cbind(DIR,DIRg)

  B=b %*% t(one.w)

  # start loops on steps
  loopgo=1
  stepi=1.0

  while(loopgo==1 & stepi>=epstop){     # second while

  (B1=B+stepi*DIR)
  FVAL=apply(B1,2,'obj',...)
#  (FVAL=apply(B1,2,'obj'))

  (fval1=max(FVAL))

  if(fval1>(fval0+epstop)) {b=B1[,min(which(FVAL==fval1))]; loopgo=0}
   stepi=0.1*stepi
   }                                    # end second while

   stepi
   epstop
   fval0
   fval1

  if(loopgo==1) stopum=1
  iter=iter+1
  (fval0=obj(b,...))
#  (fval0=obj(b))

  fhist[iter]=fval0
if(plotum==T&(iter%%10==0|iter>=itlim)) plot(1:iter,fhist[1:iter],main=paste(ptitle,'iter',iter,round(fval0,2),round(gradtest,2)))

 }                                      # end first while

  (fval0=obj(b,...))
#  (fval0=obj(b))

  (gradb=grad(obj,b,...))
#  (gradb=grad(obj,b))

  (gradpos=gradp(obj,b,...))
  (gradneg=gradn(obj,b,...))


  (gradtest=as.numeric(sum(abs(gradb))))
  (hess=hessian(obj,b,...))
#  (hess=hessian(obj,b))

  (sigma=-ginv(hess))
   ehess=eigen(hess)
   v=ehess$vectors
   (hesstest=max(ehess$values))
  status=1
  if(gradtest<=epstop) {status=0}
  status


#  dftmp=obj(b,returnpars=T,...)
#  b.restrict=dftmp$b
#print(b.restrict)


returnum=list(b=b,obj=fval0,gradb=gradb,gradpos=gradpos,gradneg=gradneg,
hess=hess,sigma=sigma,status=status,iter=iter,hesstest=hesstest)
returnum
} # end function maxgdhess
#############################################

#############################################
maxgdhess2=function(obj,b,itlim=100,epstop=0.0000001,plotum=F,ptitle='',ndirs=1,useDIRg=F,...)   {
# library(MASS)
# library(numDeriv)
 dirgen=function(w,d1,d2){(1-w)*d1+w*d2}
  w=as.matrix(seq(0,1,1/ndirs))
 npars=length(b)
 one.w=matrix(1,(length(w)),1)
 if(useDIRg==T) one.w=matrix(1,(length(w)+npars),1)
 stopum=0
 iter=1
 (fval0=obj(b,...))
# (fval0=obj(b))
 fhist=fval0

 while(stopum==0&iter<=itlim) {         # first while


#  dftmp=obj(b,returnpars=T,...)
#  b.restrict=dftmp$b
#print(b.restrict)
#  b=b.restrict


  (fval0=obj(b,...))
#  (fval0=obj(b))
  (gradb=gradum(obj,b,...))
#  (gradb=gradum(obj,b))
  (gradtest=as.numeric(sum(abs(gradb))))

  (hess=hessum(obj,b,...))
#  (hess=hessum(obj,b))

  if(npars==1) {hessnd=-1*abs(hess);hesstest=hess}
  if(npars>=2) {
   (ehess=eigen(hess))
   (v=ehess$vectors)
   (d=diag(ehess$values))
   (hesstest=max(d))
   (d=-1*abs(d)  )
   (hessnd=v %*% d %*% t(v))
   }
  (if(gradtest<=epstop&hesstest<0){stopum=1})
  d1=-1* ginv(hessnd) %*% gradb
  d2=gradb
  (DIR=apply(w,1,'dirgen',d1=d1,d2=d2))  # generates an npar by ndirs+1 matrix of directions
  (fval0=obj(b,...))
#  (fval0=obj(b))
  (gradb=gradum(obj,b,...))
#  (gradb=gradum(obj,b))
  (gradtest=as.numeric(sum(abs(gradb))))

  (hess=hessum(obj,b,...))
#  (hess=hessum(obj,b))

  if(npars==1) {hessnd=-1*abs(hess);hesstest=hess}
  if(npars>=2) {
   (ehess=eigen(hess))
   (v=ehess$vectors)
   (d=diag(ehess$values))
   (hesstest=max(d))
   (d=-1*abs(d)  )
   (hessnd=v %*% d %*% t(v))
   }
  (if(gradtest<=epstop&hesstest<0){stopum=1})
  d1=-1* ginv(hessnd) %*% gradb
  d2=gradb
  (DIR=apply(w,1,'dirgen',d1=d1,d2=d2))  # generates an npar by ndirs+1 matrix of directions

   DIRg=diag(gradb)
   if(useDIRg==T)  DIR=cbind(DIR,DIRg)

  B=b %*% t(one.w)

  # start loops on steps
  loopgo=1
  stepi=1.0

  while(loopgo==1 & stepi>=epstop){     # second while

  (B1=B+stepi*DIR)
  FVAL=apply(B1,2,'obj',...)
#  (FVAL=apply(B1,2,'obj'))


  (fval1=max(FVAL))

  if(fval1>(fval0+epstop)) {b=B1[,min(which(FVAL==fval1))]; loopgo=0}
   stepi=0.1*stepi
   }                                    # end second while

   stepi
   epstop
   fval0
   fval1


  if(loopgo==1) stopum=1
  iter=iter+1
  (fval0=obj(b,...))
#  (fval0=obj(b))

  fhist[iter]=fval0
if(plotum==T&(iter%%10==0|iter>=itlim)) plot(1:iter,fhist[1:iter],main=paste(ptitle,'iter',iter,round(fval0,2),round(gradtest,2)))

 }                                      # end first while

  (fval0=obj(b,...))
#  (fval0=obj(b))

  (gradb=gradum(obj,b,...))
#  (gradb=gradum(obj,b))

  (gradpos=gradp(obj,b,...))
  (gradneg=gradn(obj,b,...))


  (gradtest=as.numeric(sum(abs(gradb))))
  (hess=hessum(obj,b,...))
#  (hess=hessum(obj,b))

  (sigma=-ginv(hess))
   ehess=eigen(hess)
   v=ehess$vectors
   (hesstest=max(ehess$values))
  status=1
  if(gradtest<=epstop) {status=0}
  status


#  dftmp=obj(b,returnpars=T,...)
#  b.restrict=dftmp$b
#print(b.restrict)


returnum=list(b=b,obj=fval0,gradb=gradb,gradpos=gradpos,gradneg=gradneg,
hess=hess,sigma=sigma,status=status,iter=iter,hesstest=hesstest)
returnum
}
#######################################
maxgdhessr = function(obj,b,itlim=100,epstop=0.0000001,plotum=F,ptitle='',
              ndirs=1,useDIRg=F,...)   {

# library(MASS)
# library(numDeriv)
 dirgen=function(w,d1,d2){(1-w)*d1+w*d2}
  w=as.matrix(seq(0,1,1/ndirs))
 npars=length(b)
 one.w=matrix(1,(length(w)),1)
 if(useDIRg==T) one.w=matrix(1,(length(w)+npars),1)
 stopum=0
 iter=1
 (fval0=obj(b,...))
# (fval0=obj(b))
 fhist=fval0

 while(stopum==0&iter<=itlim) {         # first while

#debugum
  dftmp=obj(b,returnpars=T,...)
  b.restrict=dftmp$b
  b=b.restrict


  (fval0=obj(b,...))
#  (fval0=obj(b))
  (gradb=grad(obj,b,...))
#  (gradb=grad(obj,b))


  (gradtest=as.numeric(sum(abs(gradb))))

  (hess=hessian(obj,b,...))
#  (hess=hessian(obj,b))

  if(npars==1) {hessnd=-1*abs(hess);hesstest=hess}
  if(npars>=2) {
   (ehess=eigen(hess))
   (v=ehess$vectors)
   (d=diag(ehess$values))
   (hesstest=max(d))
   (d=-1*abs(d)  )
   (hessnd=v %*% d %*% t(v))
   }
  (if(gradtest<=epstop&hesstest<0){stopum=1})
  d1=-1* ginv(hessnd) %*% gradb
  d2=gradb
  (DIR=apply(w,1,'dirgen',d1=d1,d2=d2))  # generates an npar by 11 matrix of directions
  (fval0=obj(b,...))
#  (fval0=obj(b))
  (gradb=grad(obj,b,...))
#  (gradb=grad(obj,b))
  (gradtest=as.numeric(sum(abs(gradb))))

  (hess=hessian(obj,b,...))
#  (hess=hessian(obj,b))

  if(npars==1) {hessnd=-1*abs(hess);hesstest=hess}
  if(npars>=2) {
   (ehess=eigen(hess))
   (v=ehess$vectors)
   (d=diag(ehess$values))
   (hesstest=max(d))
   (d=-1*abs(d)  )
   (hessnd=v %*% d %*% t(v))
   }
  (if(gradtest<=epstop&hesstest<0){stopum=1})
  d1=-1* ginv(hessnd) %*% gradb
  d2=gradb
  (DIR=apply(w,1,'dirgen',d1=d1,d2=d2))  # generates an npar by 11 matrix of directions

   DIRg=diag(gradb)
   if(useDIRg==T)  DIR=cbind(DIR,DIRg)


  B=b %*% t(one.w)
  # start loops on steps
  loopgo=1
  stepi=1.0

  while(loopgo==1 & stepi>=epstop){     # second while

  (B1=B+stepi*DIR)
  FVAL=apply(B1,2,'obj',...)
#  (FVAL=apply(B1,2,'obj'))


  (fval1=max(FVAL))

  if(fval1>(fval0+epstop)) {b=B1[,min(which(FVAL==fval1))]; loopgo=0}
   stepi=0.1*stepi
   }                                    # end second while

   stepi
   epstop
   fval0
   fval1


  if(loopgo==1) stopum=1
  iter=iter+1
  (fval0=obj(b,...))
#  (fval0=obj(b))

  fhist[iter]=fval0
if(plotum==T&(iter%%10==0|iter>=itlim)) {
plot(1:iter,fhist[1:iter],main=paste(ptitle,'iter',iter,round(fval0,2),round(gradtest,6)))}

 }                                      # end first while

  (fval0=obj(b,...))
#  (fval0=obj(b))

  (gradb=grad(obj,b,...))
#  (gradb=grad(obj,b))

  (gradpos=gradp(obj,b,...))
#  (gradpos=gradp(obj,b))

  (gradneg=gradn(obj,b,...))
#  (gradneg=gradn(obj,b))


  (gradtest=as.numeric(sum(abs(gradb))))
  (hess=hessian(obj,b,...))
#  (hess=hessian(obj,b))

  (sigma=-ginv(hess))
   ehess=eigen(hess)
   v=ehess$vectors
   (hesstest=max(ehess$values))
  status=1
  if(gradtest<=epstop) {status=0}
  status


  dftmp=obj(b,returnpars=T,...)
  b.restrict=dftmp$b
#print(b.restrict)


returnum=list(b=b,obj=fval0,gradb=gradb,gradpos=gradpos,gradneg=gradneg,hess=hess,sigma=sigma,status=status,iter=iter,hesstest=hesstest,b.restrict=b.restrict)
returnum

}
#######################################
maxgdhessr2=function(obj,b,itlim=100,epstop=0.0000001,plotum=F,
            ptitle='',ndirs=1,useDIRg=F,...)   {

# library(MASS)
# library(numDeriv)
 dirgen=function(w,d1,d2){(1-w)*d1+w*d2}
  w=as.matrix(seq(0,1,1/ndirs))
 npars=length(b)
 one.w=matrix(1,(length(w)),1)
 if(useDIRg==T) one.w=matrix(1,(length(w)+npars),1)
 stopum=0
 iter=1
 (fval0=obj(b,...))
# (fval0=obj(b))
 fhist=fval0

 while(stopum==0&iter<=itlim) {         # first while

#debugum
  dftmp=obj(b,returnpars=T,...)
  b.restrict=dftmp$b
  b=b.restrict


  (fval0=obj(b,...))
#  (fval0=obj(b))
  (gradb=gradum(obj,b,...))
#  (gradb=gradum(obj,b))

  (gradtest=as.numeric(sum(abs(gradb))))

  (hess=hessum(obj,b,...))
#  (hess=hessum(obj,b))

  if(npars==1) {hessnd=-1*abs(hess);hesstest=hess}
  if(npars>=2) {
   (ehess=eigen(hess))
   (v=ehess$vectors)
   (d=diag(ehess$values))
   (hesstest=max(d))
   (d=-1*abs(d)  )
   (hessnd=v %*% d %*% t(v))
   }
  (if(gradtest<=epstop&hesstest<0){stopum=1})
  d1=-1* ginv(hessnd) %*% gradb
  d2=gradb
  (DIR=apply(w,1,'dirgen',d1=d1,d2=d2))  # generates an npar by 11 matrix of directions
  (fval0=obj(b,...))
#  (fval0=obj(b))
  (gradb=gradum(obj,b,...))
#  (gradb=gradum(obj,b))
  (gradtest=as.numeric(sum(abs(gradb))))

  (hess=hessum(obj,b,...))
#  (hess=hessum(obj,b))

  if(npars==1) {hessnd=-1*abs(hess);hesstest=hess}
  if(npars>=2) {
   (ehess=eigen(hess))
   (v=ehess$vectors)
   (d=diag(ehess$values))
   (hesstest=max(d))
   (d=-1*abs(d)  )
   (hessnd=v %*% d %*% t(v))
   }
  (if(gradtest<=epstop&hesstest<0){stopum=1})
  d1=-1* ginv(hessnd) %*% gradb
  d2=gradb
  (DIR=apply(w,1,'dirgen',d1=d1,d2=d2))  # generates an npar by 11 matrix of directions

   DIRg=diag(gradb)
   if(useDIRg==T)  DIR=cbind(DIR,DIRg)


  B=b %*% t(one.w)
  # start loops on steps
  loopgo=1
  stepi=1.0

  while(loopgo==1 & stepi>=epstop){     # second while

  (B1=B+stepi*DIR)
  FVAL=apply(B1,2,'obj',...)
#  (FVAL=apply(B1,2,'obj'))


  (fval1=max(FVAL))

  if(fval1>(fval0+epstop)) {b=B1[,min(which(FVAL==fval1))]; loopgo=0}
   stepi=0.1*stepi
   }                                    # end second while

   stepi
   epstop
   fval0
   fval1


  if(loopgo==1) stopum=1
  iter=iter+1
  (fval0=obj(b,...))
#  (fval0=obj(b))

  fhist[iter]=fval0
if(plotum==T&(iter%%10==0|iter>=itlim))  plot(1:iter,fhist[1:iter],main=paste(ptitle,'iter',iter,round(fval0,2),round(gradtest,6)))

 }                                      # end first while

  (fval0=obj(b,...))
#  (fval0=obj(b))

  (gradb=gradum(obj,b,...))
#  (gradb=gradum(obj,b))

  (gradpos=gradp(obj,b,...))
#  (gradpos=gradp(obj,b))

  (gradneg=gradn(obj,b,...))
#  (gradneg=gradn(obj,b))


  (gradtest=as.numeric(sum(abs(gradb))))
  (hess=hessum(obj,b,...))
#  (hess=hessum(obj,b))

  (sigma=-ginv(hess))
   ehess=eigen(hess)
   v=ehess$vectors
   (hesstest=max(ehess$values))
  status=1
  if(gradtest<=epstop) {status=0}
  status


  dftmp=obj(b,returnpars=T,...)
  b.restrict=dftmp$b
#print(b.restrict)

returnum=list(b=b,obj=fval0,gradb=gradb,gradpos=gradpos,gradneg=gradneg,
hess=hess,sigma=sigma,status=status,iter=iter,hesstest=hesstest,
b.restrict=b.restrict)
returnum
}
################################################################################



################################################################################
#' pointpoint: Plots line through two points
#'
#' @param x1         First point.
#' @param x2         Second point.
#' @param lty        Line type. Default = 1.
#' @param lwd          Line width. Default = 1.
#' @param col        Line color. Default = 1.
#' @param font       Font. Default = 1.
#' @param xlim       x limits for line. No intial limits
#' @return nothing   Plots a line.
pointpoint=function(x1,x2,lty=1,lwd=1,col=1,font=1,
    xlim=c(-1e+06,1e+05),offset=0,angle=30) {
  (m=(x2[2]-x1[2])/(x2[1]-x1[1]))
  pointslope(x=x1,m=m,lty=lty,lwd=lwd,col=col,font=font,
      xlim=xlim)
}  # end function pointpoint
########################################

########################################
#' pointslope: Plots line and/or arrow through point given slope.
#'
#' @param x1         First point.
#' @param x2         Second point.
#' @param lty        Line type. Default = 1.
#' @param lwd          Line width. Default = 1.
#' @param col        Line color. Default = 1.
#' @param font       Font. Default = 1.
#' @param xlim       x limits for line. No intial limits
#' @param arrow      Plot arrows. Default = F
#' @param R          Radius (Length of arrow along main shaft).
#' @param length     Argument to R's arrows function. Default 0.15
#' @param offset     Argument to R's arrows function. Default 0.15
#' @param angle      Argument to R's arrows function. Default 30.
#' @return nothing   Plots a line.
pointslope=function(x,m,lty=1,lwd=1,col=1,font=1,xlim=c(-1e+06,1e+05),
  arrow=FALSE,R=1,length=0.15,offset=0,angle=30)  {
    a = x[2] - m * x[1]
    xp = xlim
    yp = a + m * xp
    if(arrow==F)  {
     if(abs(m)!=Inf) points(xp, yp, type = "l", lty = lty, lwd = lwd, col = col)
     if(abs(m)==Inf) abline(v=x[1],lty = lty, lwd = lwd, col = col)
     }
    if(arrow==T) {
     theta=pi/2
     if(abs(m)!=Inf) theta=atan(m)+offset
     x2=x+c(R*cos(theta),R*sin(theta))
     arrows(x[1],x[2],x2[1],x2[2],length=length,angle=angle,col=col,lwd=lwd)
    }
} # end function pointslope
########################################
## Examples for pointslope and pointpoint
# plot(10,10,xlim=c(-10,20),ylim=c(-10,20),type='n')
# abline(h=0,lwd=2,font=2)
# abline(v=0,lwd=2,font=2)
# (point1=c(10,10))
# points(point1[1],point1[2],pch=20,cex=2,col=4)
# pointslope(point1,m=2,lwd=3)
# pointslope(point1,m=1,lty=3,lwd=3,col=2,xlim=c(6,14))
# point2=c(15,5)
# pointpoint(point1,point2,col=3,lwd=2,lty=2,xlim=c(6,14))
# pointpoint(point1,point1+c(-1,3),col=1,lwd=2,lty=2,xlim=c(8,12))
# points(point1[1],point1[2],pch=20,cex=3,col=4)
###############################################################################


###############################################################################
radcalc=function(x1,y1,x0=0,y0=0,degrees=F){
 x=x1-x0;y=y1-y0
 theta=ifelse(x>0&y>=0,atan(y/x),NA)
 theta=ifelse(x<0&y>=0,pi-atan(y/-x),theta)
 theta=ifelse(x<0&y<0,(3/2)*pi-atan(-y/-x),theta)
 theta=ifelse(x<0&y<0,(3/2)*pi-atan(-y/-x),theta)
 theta=ifelse(x>0&y<0,(4/2)*pi-atan(-y/x),theta)
 theta=ifelse(x==0&y>0,(1/2)*pi,theta)
 theta=ifelse(x==0&y<0,(3/2)*pi,theta)
 if(degrees==T|degrees=='TRUE') theta=theta*360/(2*pi)
 theta
} # end function radcalc
########################################
 cos_d=function(dg){
  cos(dg/(360/(2*pi)))
 }
########################################
sin_d=function(dg){
  sin(dg/(360/(2*pi)))
 }
################################################################################


################################################################################
summaryDF=function(data, splits, operate.on.columns, functions = "mean",...){
   ## function to create a summary data frame
   ##  data is the original data to use
   ##  splits = column names which describe the rows to combine
   ##  operate.on.columns = column names which give the columns to summarize
   ##  functions = the function to apply to the elements of one column
   ##     having the same splits.
#########################################################################
###########################################################################
# example front end code
#df1=read.table(file.choose(),header=T,sep=",")
#df1$fips=df1$nst*1000+df1$ncty
#names(df1)
#[1] "yr"         "nst"        "crd"        "ncty"       "pacre"
#[6] "hacre"      "hyld"       "production" "fips"
#######################################
#source('c:\\rcode\\summaryDF.r')
#######################################
#df.tmp=summaryDF(df1,c("nst","yr"), c("hacre","hyld"), c("mean"))
#df.tmp=summaryDF(df1,c("nst","yr"), c("hacre","hyld"), c("mean"),na.rm=TRUE)
#or write your own function
#my.sum = function(x)  {sum(x, na.rm=TRUE)}
#df.tmp1=summaryDF(df7,c("mcty","myr"), c("pcpc","pcpc"),c("sum","my.sum"))
# rename columns
#dimnames(df.tmp1)[[2]] = c("county","yr","pcpc1","pcpc2")
###########################################################################
########################################################################
   if(!is.numeric(splits)){
     if( !all(splits %in% names(data)))
       stop("Second argument must be columns of the first argument.")
   }
   else {
     if(any(splits < 0) | any(splits > length(data))){


      stop("Second argument is < 0 or > dim(data)")
    }
   }
   split.names =  names(data[,splits, drop=FALSE ])


   if(!is.numeric(operate.on.columns)){


     if( !all(operate.on.columns %in% names(data)))
       stop("Third argument must be columns of the first argument.")
   }
   else if(any(operate.on.columns < 0) | any(operate.on.columns > length(data))){


      stop("Third argument is < 0 or > dim(data)")
    }

   if((op.len = length(operate.on.columns)) > length(functions))
       functions = rep(functions, op.len)[1:op.len]
       ##  Use the same functions repeatedly


   ##  start building up combinations of levels of split columns

   combos = as.character(unlist(data[,splits[1]]))

   if((s.len = length(splits)) == 1){
     ## only one split column
     id.columns = list(levels(as.factor(data[,splits])))
   }
   else {                         ## more than one split
      for(i in splits[-1] )
        combos = paste(combos,as.character(unlist(data[,i])),sep="@")
 #     uniqCombosplit =  unlist(lapply(list(unique(combos)), strsplit,"@")[[1]])
 #     id.columns = data.frame(
 #               uniqCombosplit[seq(1,length(uniqCombosplit), s.len)])
 #     if(is.numeric(data[,splits[1]]))
 #        id.columns[,1] = as.numeric(as.character(id.columns[,1]))
 #     for(i in 2:s.len ){
 #       id.columns = cbind(id.columns,
 #                uniqCombosplit[seq(i,length(uniqCombosplit), s.len)])
 #       if(is.numeric(data[,splits[i]]))
 #         id.columns[,i] = as.numeric(as.character(id.columns[,i]))
 #     }
    }


 #  names(id.columns) = splits
   ##  build first column to be output:
   col1 = as.numeric(
               tapply(data[,operate.on.columns[1]], combos, eval(functions[1]),...))



   df = data.frame( col1)
   ## grab names to match farms with proper rows
   if(s.len == 1){
     id.columns =  data.frame(c1=names(tapply(combos,combos,length)))
   }
   else{
       # split apart the combo strings used in the split
      rowlabels = strsplit(names(tapply(combos,combos,length)),"@")
        ## a list of split strings

      id.columns = data.frame( c1 =
                      unlist(lapply(rowlabels, function(x) x[1])))
       ## first split labels in column 1


      ## check to see if it's factor coded or numeric:
       if(is.numeric(data[,splits[1]]))
          id.columns[,1] = as.numeric(as.character(id.columns[,1]))



      for(i in 2:s.len ){
        id.columns = cbind(id.columns,
                          unlist( lapply(rowlabels,function(x) x[i])))
        ## additional split labels in more columns



        ## check to see if it's factor coded or numeric:
        if(is.numeric(data[,splits[i]]))
          id.columns[,i] = as.numeric(as.character(id.columns[,i]))
      }

    }
   ##  give a name to each column
   names(id.columns) = split.names


   if(op.len > 1){
     for(i in 2:op.len){
       df = cbind(df, as.numeric(
         tapply(data[,operate.on.columns[i]], combos,eval(functions[i]),...)))
     }
   }


   names(df) = paste(names(data[,operate.on.columns,drop=FALSE]),functions,sep=".")

  df2=data.frame(id.columns, df)

# Jim  I added the following.  Joe

  for(i in 1:length(splits))  {
   if(is.numeric(data[,splits[i]]))  {
     df2[,splits[i]] = as.numeric(as.character(df2[,splits[i]]))
      }
  }


   df2

 } # end function summaryDF
################################################################################



################################################################################
waittalk=function(nsecs=5,talk=1){
 sec0 = my_seconds()
 secdif = 0
 t1=0
 while(secdif<nsecs) {
  secdif=my_seconds()-sec0
  if((secdif-t1)>talk){
   print(paste('timeleft = ',round(nsecs-secdif,1)))
   t1=secdif
  }
 }
} # end function waittalk
###############################################################################


###############################################################################
# PVAF with discrete or continuous compounding
###############################################################################
PVAF=function(r,N,compound='D',FPP=1){
 if(r==0) pvaf=N        
 if(r>0){ 
        if(compound=='D'|compound=='d') {
    pvaf=((1+r)^(1-FPP)) * ((1-(1+r)^-N)/r)
  }
  if(compound=='C'|compound=='c') {
    pvaf = (exp(r)^(1-FPP)) * ((1-(exp(r))^-N)/(exp(r)-1))
  }
 } # end         if(r>0)
pvaf
###############################################################################
} # end function PVAF
###############################################################################

###############################################################################
# FVAF with discrete or continuous compounding
###############################################################################
FVAF=function(r,N,compound='C',FPP=1){
 if(r==0) fvaf=N        
 if(r>0){       
  pvaf=PVAF(r,N,compound,FPP)
  if(compound=='D'|compound=='d') {
    fvaf=((1+r)^N) * pvaf
  }
  if(compound=='C'|compound=='c') {
    fvaf=(exp(r)^N) * pvaf
  }
 } # end         if(r>0)
fvaf
} # end function FVAF
###############################################################################

###############################################################################
# debug my_NPV
#CF=c(-432.95,rep(100,5)); CFTimes=0:5; r=0.05
###############################################################################
my_NPV = function(CF=c(-432.95,rep(100,5)), CFTimes=0:5, r=0.05){
 if(length(CF)!=length(CFTimes)) stop('length CF != length(CFTimes')    
 sum(CF/((1+r)^CFTimes))        
} # end function my_NPV
###############################################################################
# test my_NPV
#my_NPV(CF=c(-432.95,rep(100,5)), CFTimes=0:5, r=0.05)
###############################################################################

###############################################################################
# debug my_NPV
#CF=c(-432.95,rep(100,5)); CFTimes=0:5; rLB=0; rUB=5; ngrid=1000
###############################################################################
my_IRR = function(CF=c(-432.95,rep(100,5)), CFTimes=0:5,rLB=0,rUB=5,ngrid=1000){
############################################################################### 
 if(length(CF)!=length(CFTimes)) stop('length CF != length(CFTimes')    
 rlist=seq(rLB,rUB,length.out=ngrid)
 ABSNPV=0
 for(j in 1:ngrid){
        ABSNPV[j]=abs(my_NPV(CF,CFTimes,rlist[j]))
        }

#plot(rlist,ABSNPV,type='l') 
(ABSmin=min(ABSNPV))    

(pick=which(ABSNPV==min(ABSNPV)))
(IRR=rlist[pick])
list(IRR=IRR,pick=pick,ABSNPV)
} # end function my_IRR
###############################################################################
# test my_IRR
# my_IRR(CF=c(-432.95,rep(100,5)), CFTimes=0:5,rLB=0,rUB=5,ngrid=1000)$IRR
# my_IRR(CF=c(-432.95,rep(100,5)), CFTimes=0:5,rLB=0,rUB=5,ngrid=10000)$IRR
###############################################################################



###############################################################################
pull_qmod = function(SYM='QQQ',date0="2000-01-01") {
###############################################################################
  my_require('quantmod')
  out = tryCatch(
    {
      # Just to highlight: if you want to use more than one
      # R expression in the "try" part then you'll have to
      # use curly brackets.
      # 'tryCatch()' will return the last evaluated expression
      # in case the "try" part was completed successfully

      #message("This is the 'try' part")
      as.data.frame(getSymbols(SYM,auto.assign = FALSE, from = date0))
    },
    error=function(cond) {
      message(paste("error", cond))
      #message("Here's the original error message:")
      #message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    warning=function(cond) {
      message(paste("warning:", cond))
      #message("Here's the original warning message:")
      #message(cond)
      ## Choose a return value in case of warning
      return(NA)
    },
    finally={
      ## NOTE:
      ## Here goes everything that should be executed at the end,
      ## regardless of success or error.
      ## If you want more than one expression to be executed, then you
      ## need to wrap them in curly brackets ({...}); otherwise you could
      ## just have written 'finally=<expression>'
      #message(paste("Processed URL:", url))
      #message("Some other message at the end")
    }
  )
  return(out)
################################################################################
}# end function pull_qmod
################################################################################

################################################################################
my_substrng = function(x,keepS=NULL,dropS=NULL){
 tmp=x
 if(!is.null(keepS)) {
   nK = length(keepS)
   for(i in 1:nK){
    pick = grep(keepS[i],tmp)
    if(length(pick)>=1) tmp=tmp[pick]
   }
 } # end   if(!is.null(keepS))

 if(!is.null(dropS)) {
   nD = length(dropS)
   for(i in 1:nD){
    pick = as.integer(grep(dropS[i],tmp))
    if(length(pick)>=1) tmp=tmp[-pick]
   }
 } # end   if(!is.null(dropS))
 tmp
} # end function my_substrng
################################################################################

################################################################################
my_cut=function(x,breaks=c(min(x),mean(x),max(x)),...){
  cut(x,breaks=breaks,include.lowest=T,...)
} # end function my_cut
################################################################################

###############################################################################
my_ecdf=function(x,xL=min(x,na.rm=T),xU=max(x,na.rm=T)){
 x=sort(x[!is.na(x)])
 n=length(x)
 x=c(x)
 F=c((1:n)/(n))
 xF = cbind(x,F)
 list(x=x,F=F,xF=xF)
 } # end function my_ecdf
###############################################################################


###############################################################################
plot_ecdf=function(x,xL=min(x,na.rm=T),xU=max(x,na.rm=T),main='ecdf(x)',
          pch=19,cex=2,lty=1,col=1,lwd=2,add=F){
###############################################################################
      x=sort(x[!is.na(x)])
      n=length(x)
      x=c(x)
      F=c((1:n)/(n))
      xF = cbind(x,F)
      tmp=list(x=x,F=F,xF=xF)
#     tmp=my_ecdf(x,xL,xU)
     (x=tmp$x)
     (F=tmp$F)
     (n=length(x))

     (xL=min(x)-1)
     (xU=max(x)+1)


     xpick=rep(1:n,each=2)
     Fpick=xpick

     (xp=c(xL,x[xpick],xU))
     (Fp=c(0,0,F[Fpick]))

     if(add!=T){
      plot(x,F,xlim=range(xp),ylim=range(Fp),main=main,
          pch=pch,cex=cex,col=col)
     } # end if
     points(xp,Fp,type='l',lwd=lwd,lty=lty,col=col)

} # end function plot_ecdf
###############################################################################

###############################################################################
#graphics.off()
#x=rnorm(50)
#y=rnorm(50)
#x11()
#plot(ecdf(x),lwd=3)
#plot(ecdf(y),lwd=3,add=T,col=2)
#
#x11()
#plot_ecdf(x=x,xL=min(x)-1,xU=max(x)+1,cex=1.5,pch=20,lwd=3)
#x11()
#plot_ecdf(x,xL=min(x)-1,xU=max(x)+1,cex=0.5,pch=20,lwd=3,main='')
#plot_ecdf(y,xL=min(y)-1,xU=max(y)+1,cex=0.5,pch=20,lwd=3,col=2,add=T)
###############################################################################


###############################################################################
# Function for converting data frame fields
###############################################################################
my_convert=function(df,ctypes){
  for(j in 1:length(ctypes)) {
   if(ctypes[j]=='c') df[,j] = as.character(df[,j])
   if(ctypes[j]=='n') df[,j] = as.numeric(as.character(df[,j]))
   if(ctypes[j]=='i') df[,j] = as.integer(as.character(df[,j]))
  }
  df
} # end function my_convert
###############################################################################


###############################################################################
# Functions to check objects contained in a .RData file. The file name 'fname'
# will need to include the path if the file being checked is not in the getwd()
# directory
###############################################################################
# queryRData=function(fname){
#   myfun2=function(fname){
#    load(fname)
#    rm(fname)
#    list=ls()
#   } # end functionmyfun2
#  list=myfun2(fname); gc()
#  list
# } # end function queryRData
###############################################################################
queryRData=function(fname,summary=F){
  myfun2=function(fname,summary=F){
    SUMMARY=NULL
    (list1=c(ls(),'list1'))
    load(fname)
    (list2=ls())
    (list3=setdiff(list2,list1))
    obj3=mget(list3)
#    if(summary==T){
     (obj_names = sapply(obj3,names))
     (obj_class =  sapply(obj3,class))
     (obj_summary=NULL) 
     if(summary==T)(obj_summary =  sapply(obj3,summary))
    SUMMARY=list(obj_names=obj_names,obj_class=obj_class,obj_summary=obj_summary)
 #   } # end if if(summary==T)
    list(fname=fname,objects=list3,SUMMARY=SUMMARY)
  } # end function myfun2
  tmp=myfun2(fname,summary); gc()
  tmp
} # end function queryRData
################################################################################
checkRData = function(file) {
 #' Function for listing objects object in an .RData file created by R's 
 #' save() command
 #' Inputs: RData file
 E = new.env()
 load(file=file, envir=E)
 objlist=ls(envir=E)
 rm(E); gc()
 return(objlist)
}  # end function checkRData
###############################################################################
getRData = function(file, object) {
 #' Function for extracting an object from a .RData file created by R's save() command
 #' Inputs: RData file, object name
 #' object names should be in quotes i.e. 'object'
 E = new.env()
 load(file=file, envir=E)
 return(get(object, envir=E, inherits=F))
}  # end function getRData
###############################################################################
loadRData=function(fname='test.Rdata',objects=NULL){
 #' Function for extracting an vector of named objects from an .RData file 
 #' created by R's save() command. Places set of objects on global environment
 #' Inputs: RData file, object names(character vector)
 if(is.null(objects)) objects=checkRData(fname)  
 (filelist=checkRData(fname))  
 load(fname)
 for(j in 1:length(objects)){
  object=pulldata[j]   
  tmp=list(get(object,inherits=F))
  names(tmp)=object
  list2env(tmp,envir=.GlobalEnv)
 }
 rm(list=ls())
} # end function loadRData
###############################################################################



##############################################################################
my_substrng = function(x,keepS=NULL,dropS=NULL){
 tmp=x
 if(!is.null(keepS)) {
   nK = length(keepS)
   for(i in 1:nK){
    pick = grep(keepS[i],tmp)
    if(length(pick)>=1) tmp=tmp[pick]
   }
 } # end   if(!is.null(keepS))

 if(!is.null(dropS)) {
   nD = length(dropS)
   for(i in 1:nD){
    pick = as.integer(grep(dropS[i],tmp))
    if(length(pick)>=1) tmp=tmp[-pick]
   }
 } # end   if(!is.null(dropS))
 tmp
} # end function my_substrng
###############################################################################

##################################################################################
try_catch_nassqs = function(params) {
  out = tryCatch(
    {
      # Just to highlight: if you want to use more than one
      # R expression in the "try" part then you'll have to
      # use curly brackets.
      # 'tryCatch()' will return the last evaluated expression
      # in case the "try" part was completed successfully

      #message("This is the 'try' part")
      nassqs(params)
    },
    error=function(cond) {
      message(paste("error", cond))
      #message("Here's the original error message:")
      #message(cond)
      # Choose a return value in case of error
      return(NA)
    },
    # warning=function(cond) {
    #   message(paste("warning:", cond))
    #   #message("Here's the original warning message:")
    #   #message(cond)
    #   ## Choose a return value in case of warning
    #   return(NA)
    # },
    finally={
      ## NOTE:
      ## Here goes everything that should be executed at the end,
      ## regardless of success or error.
      ## If you want more than one expression to be executed, then you
      ## need to wrap them in curly brackets ({...}); otherwise you could
      ## just have written 'finally=<expression>'
      #message(paste("Processed URL:", url))
      #message("Some other message at the end")
    }
  )
  return(out)
} # end function try_catch_nassqs
###############################################################################
my_nassqs=function(params,ntrys=5){
  trys=0; GO=T 
  while(trys<=ntrys&GO==T){
    trys=trys+1 
    print(paste('trying data download',trys))
    df0=try_catch_nassqs(params)
    if(class(df0)=='data.frame') GO=F
    if(class(df0)=='logical') df0=NULL
  } # end while
df0 
} # end function my_nassqs
###############################################################################

###############################################################################
# debug on
#year=df00$yr; values=df00$yld;niter=2000
###############################################################################
my_detrend=function(year,values,niter=2000) {
############################################################################### 
if(length(year)!=length(values)) stop('length(year)!=length(values)')
###############################################################################
# Dataframe for F tests
M=1:5; nM = 5
RM=rep(M,nM)
UM=rep(M,each=nM)
dF=data.frame(cbind(RM,UM))
dF=subset(dF,RM<UM)
dF=orderBy(~RM+UM,data=dF)
dF$KR=M[dF$RM]
dF$KU=M[dF$UM]
##############################################################################
BPAR_A=matrix(NA,1,1)
BPAR_B=matrix(NA,1,2)
BPAR_C=matrix(NA,1,3)
BPAR_D=matrix(NA,1,4)
BPAR_E=matrix(NA,1,5)
SSE=matrix(NA,1,5)
HETPAR=SSE
HETPVAL=SSE


DF0=vector('list',1)
MODELS=vector('list',5)
###############################################################################
graphics.off()
###############################################################################
my_x11()
###############################################################################
df0=as.data.frame(cbind(year=year,values=values))
(minYr = min(df0$year) - 1)
df0$yr = (df0$year - minYr) / 10
###############################################################################
# model A
mA = function(x, dataIn=df0) {
 ehat=dataIn$values - x[1]
 sum(ehat * ehat)
} # end function mB

# start value
(b = mean(df0$values))

(modelA=nlminb(start=b,objective=mA,control = list(eval.max=niter,iter.max=niter)))
modelA$b=modelA$par

(Ab=mean(df0$values))
df0$modelA = Ab
df0$EhatA=df0$values-df0$modelA

(regA1=summary(lm(abs(EhatA)~df0$yr,data=df0)))
(HETPAR[1,1]=regA1$coeff[2,1])
(HETPVAL[1,1]=regA1$coeff[2,4])


var(df0$values)*(nrow(df0)-1)
(SSE[1,1]=mA(Ab))

BPAR_A[1,]=Ab

################################################################
# Model B
mB = function(x, dataIn=df0) {
 ehat=dataIn$values - (x[1]+x[2]*dataIn$yr)
 sum(ehat * ehat)
} # end function mB

mA(Ab)
b=c(Ab[1], 0)
mB(b)

(modelB=nlminb(start=b,objective=mB,control = list(eval.max=niter,iter.max=niter)))
modelB$b=modelB$par


#check model1 versus ols
model1_ols = lm(values~yr, data=df0)
model1_ols$coef

(Bb=modelB$b)

x=Bb
df0$modelB = (x[1]+x[2]*df0$yr)

df0$EhatB=df0$values-df0$modelB
(regB1=summary(lm(abs(EhatB)~df0$yr,data=df0)))
(HETPAR[1,2]=regB1$coeff[2,1])
(HETPVAL[1,2]=regB1$coeff[2,4])

mB(Bb)
(SSE[1,2]=mB(Bb))


BPAR_B[1,]=modelB$b
if( abs(modelB$b[1] - model1_ols$coef[1]) > 0.0001 ) print("Coeficient differences should be zero - problem with intercept")
if( abs(modelB$b[1] - model1_ols$coef[1]) > 0.0001 ) print(cfips0)      # added to identify CRD if problem
if( abs(modelB$b[2] - model1_ols$coef[2]) > 0.0001 ) print("Coeficient differences should be zero - problem with slope")
if( abs(modelB$b[2] - model1_ols$coef[2]) > 0.0001 ) print(cfips0)      # added to identify CRD if problem

##############################################################
# Model C
mC = function(x, dataIn=df0) {
   ehat = dataIn$values - ( x[1] + x[2] * dataIn$yr^x[3] )
   sum(ehat * ehat)
}

# 4) Run Model C with starting values from A

mA(Ab)
b=c(Ab[1],0 ,1)
mC(b)
#

(modelC=nlminb(start=b,objective=mC,control = list(eval.max=niter,iter.max=niter)))
modelC$b=modelC$par


(ssetest = mC(modelC$b))
MODELS[[1]]=modelC




# 4) Run Model C with starting values from B
mB(Bb)
b = c(Bb[1],Bb[2] ,1)
mC(b)

(modelC=nlminb(start=b,objective=mC,control = list(eval.max=niter,iter.max=niter)))
modelC$b=modelC$par


(ssetest[2] = mC(modelC$b))
MODELS[[2]] = modelC


(pickum = which(ssetest==min(ssetest)))
(pickum = pickum[1])
modelC = MODELS[[pickum]]


(Cb=modelC$b)
x=Cb
df0$modelC = ( x[1] + x[2] * df0$yr^x[3] )

df0$EhatC = df0$values-df0$modelC

(regC1 = summary(lm(abs(EhatC)~df0$yr,data=df0)))
(HETPAR[1,3] = regC1$coeff[2,1])
(HETPVAL[1,3] = regC1$coeff[2,4])


mC(Cb)
(SSE[1,3] = mC(Cb))
BPAR_C[1,] = modelC$b
################################################################





###############################################################################
# Model D
###############################################################################
mD = function(x, dataIn=df0) {
  ehat = dataIn$values - (x[1] + (x[2]  / (x[3] + dataIn$yr^x[4])) )
  sum(ehat * ehat)
} # end function mD


#starting values from model A
mA(Ab)
b=c(Ab[1],0,1,0)
mD(b)  # should be same value

(modelD = nlminb(start=b,objective=mD,control = list(eval.max=niter,iter.max=niter)))
modelD$b = modelD$par

(ssetest = mD(modelD$b))
MODELS[[1]] = modelD


# starting values from model B
mB(Bb)
b = c(Bb[1],Bb[2],0,-1)
mD(b)  # should be same value

(modelD=nlminb(start=b,objective=mD,control = list(eval.max=niter,iter.max=niter)))
modelD$b=modelD$par

(ssetest[2] = mD(modelD$b))
MODELS[[2]] = modelD


# starting values from model C
mC(Cb)
b = c(Cb[1], Cb[2],0,(-1)* Cb[3])
mD(b)  # should be same value

(modelD = nlminb(start=b,objective=mD,control = list(eval.max=niter,iter.max=niter)))
modelD$b = modelD$par

(ssetest[3] = mD(modelD$b))
MODELS[[3]] = modelD


pickum=which(ssetest==min(ssetest))
(pickum=pickum[1])
modelD=MODELS[[pickum]]


(Db=modelD$b)
x=Db
df0$modelD = (x[1] + (x[2]  / (x[3] + df0$yr^x[4])) )

df0$EhatD = df0$values-df0$modelD
(regD1 = summary(lm(abs(EhatD)~df0$yr,data=df0)))
(HETPAR[1,4] = regD1$coeff[2,1])
(HETPVAL[1,4] = regD1$coeff[2,4])

summary(lm(abs(EhatD)~year,data=df0))


mD(Db)
(SSE[1,4] = mD(Db))
BPAR_D[1,] = modelD$b
###############################################################################






###############################################################################
# Model E

mE = function(x, dataIn=df0) {
ehat = dataIn$values - (x[1] + (x[2]*dataIn$yr^x[3]) / (1 + x[4]*dataIn$yr^x[5]) )
sum(ehat * ehat)
} # end function mE

#############################
# starting values model A
 mA(Ab)
 (b=c(Ab[1],0,1,0,1))
 mE(b)  # should be the same

modelE = NULL

tryCatch ({
  (modelE=nlminb(start=b,objective=mE,control = list(eval.max=niter,iter.max=niter)))
},         error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

if(!is.null(modelE)) modelE$b = modelE$par


(ssetest = ifelse(!is.null(modelE),mE(modelE$b),NA))
MODELS[[1]] = modelE

if(!is.null(modelE)) mE(modelE$b)

if(is.null(modelE)) {
 print(paste('Model EA error on loop',i))
}

####################################
# starting values model B
 mB(Bb)
 (b=c(Bb[1],Bb[2],1,0,1))
 mE(b)  # should be the same

modelE=NULL
#tryCatch ({
#          modelE = maxgdhess(obj=mE,b=b,itlim=niter,ndirs=1,plotum='N',plotmod=plotmod,useDIRg='Y',
#          stepcut=stepcut,fhistum=T,bhistum=F,ptitle = paste('Model EB cfips =',cfips0))
#},         error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

tryCatch ({
  (modelE=nlminb(start=b,objective=mE,control = list(eval.max=niter,iter.max=niter)))
},         error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

if(!is.null(modelE)) modelE$b = modelE$par

(ssetest[2] = ifelse(!is.null(modelE),mE(modelE$b),NA))
MODELS[[2]] = modelE

if(!is.null(modelE)) mE(modelE$b)

if(is.null(modelE)) {
 print(paste('Model EB error on loop',i))
}

#####################################
# starting values model C
 mC(Cb)
 (b = c(Cb[1],Cb[2],Cb[3],0,1))
 mE(b)  # should be the same

modelE=NULL

tryCatch ({
  (modelE=nlminb(start=b,objective=mE,control = list(eval.max=niter,iter.max=niter)))
},         error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

if(!is.null(modelE)) modelE$b = modelE$par

(ssetest[3] = ifelse(!is.null(modelE),mE(modelE$b),NA))
MODELS[[3]] = modelE

if(!is.null(modelE)) mE(modelE$b)

if(is.null(modelE)) {
 print(paste('Model EC error on loop',i))
}

#######################################
#starting values model D
#######################################
 mD(Db)
 (b=c(Db[1],Db[2],-Db[4],Db[3],-Db[4]))
 mE(b)  # should be the same

modelE=NULL

tryCatch ({
  (modelE=nlminb(start=b,objective=mE,control = list(eval.max=niter,iter.max=niter)))
},         error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

if(!is.null(modelE)) modelE$b = modelE$par

(ssetest[4] = ifelse(!is.null(modelE),mE(modelE$b),NA))
MODELS[[4]] = modelE

if(!is.null(modelE)) mE(modelE$b)
if(is.null(modelE)) {
 print(paste('Model ED error on loop',i))
}

ssetest

(pickum=which(ssetest==min.na(ssetest)))
(pickum=pickum[1])
if(!is.na(pickum)) modelE=MODELS[[pickum]]


if (!is.null(modelE)) {
 (Eb = modelE$b)
 mE(Eb)
 x = Eb
 df0$modelE = (x[1] + (x[2]*df0$yr^x[3]) / (1 + x[4]*df0$yr^x[5]) )

 df0$EhatE = df0$values-df0$modelE
 (regE1 = summary(lm(abs(EhatE)~df0$yr,data=df0)))
 (HETPAR[1,5] = regE1$coeff[2,1])
 (HETPVAL[1,5] = regE1$coeff[2,4])


} else{
 modelE$b=b
 modelE$obj=modelD$obj
 modelE$iter=0
}

(SSE[1,5] = modelE$obj)
BPAR_E[1,] = modelE$b
##################################################################


##################################################################
(nobs = nrow(df0))
(sse=SSE[1,])
dF$SSER=sse[dF$RM]
dF$SSEU=sse[dF$UM]
dF$FSTAT=((dF$SSER-dF$SSEU)/(dF$KU-dF$KR)) / (dF$SSEU/(nobs-dF$KU))
dF$PVAL = round(1 - pf(dF$FSTAT,df1 = dF$KU-dF$KR, df2 = nobs-dF$KU),3)
dF


(PVALS=summaryBy(PVAL~RM,data=dF,FUN=min,keep.names=T))
PVALS

 (pvals=PVALS$PVAL)

 (pvals=round(pvals,3))
 (txt=paste('pvals A=',pvals[1],',B=',pvals[2],',C=',pvals[3],',D=',pvals[4],sep=''))




 MODEL='A'
 if(PVALS$PVAL[1]<=0.05){
  MODEL='B'
   if(PVALS$PVAL[2]<=0.05){
    MODEL='C'
    if(PVALS$PVAL[3]<=0.05){
      MODEL='D'
      if(PVALS$PVAL[4]<=0.05){
        MODEL='E'
      }
    }
  }
 }

 MODEL


 plot(df0$year,df0$values,type='l',lwd=5)
 points(df0$year,df0$modelA,type='l',lty=2,lwd=3)
 points(df0$year,df0$modelB,type='l',lty=2,lwd=3,col=2)
 points(df0$year,df0$modelC,type='l',lty=2,lwd=3,col=3)
 points(df0$year,df0$modelD,type='l',lty=1,lwd=3,col=4)
 points(df0$year,df0$modelE,type='l',lty=2,lwd=3,col=5)

 legend('topleft',
        legend=c('Model A','Model B','Model C','Model D','Model E'),
        fill=1:5,cex=0.75)

 mtext(paste('Model Chosen = ',MODEL),side=3,line=-1,font=2)
 mtext(paste(txt),side=3,line=-2,cex=0.65,font=2)

##############################################################################

###############################################################################

df0$MODEL=MODEL
colname1=paste('model',MODEL,sep='')
colname2=paste('Ehat',MODEL,sep='')
df0$yhat=df0[,c(colname1)]
df0$ehat=df0[,c(colname2)]

head(df0)
##############################################################################  
list(df0=df0,modelA=modelA,modelB=modelB,modelC=modelC,
     modelD=modelD,modelE=modelE)  
##############################################################################
} # end function my_detrend
##############################################################################




###############################################################################
# for debugging
#n=10000;k=5;rho=-1/3;ZCOR=NULL;seedval=1001;strat=F;pmin=1E-6;whiten=T
#n=10000;k=5;rho=0;ZCOR=NULL;seedval=1001
###############################################################################
rMVN_Z=function(n,k=2,rho=0,ZCOR=NULL,seedval=NULL,strat=F,
                pmin=min(1E-6,1/(n+1)),whiten=T) {
###############################################################################                
 if(!is.null(seedval)) set.seed(seedval)
 if(is.null(ZCOR)){
   rho
  (rho=max(rho,-(1/(k-1))))
  ZCOR=matrix(rho,k,k)
  diag(ZCOR)=1
  ZCOR
 } # end if(!is.na(R)) 
  
 tmp=eigen(ZCOR)
 v=tmp$vectors
 (d=tmp$values)
 d[d<=1E-8]=1E-8
 (d=diag(d))
 (ZCOR = v %*% d %*% t(v))
 
 if(strat==T){
  p=seq(pmin,(1-pmin),length.out=n)
  y=qnorm(p)
  Z=matrix(0,n,k)
  for(j in 1:k) Z[,j]=sample(y,n)
 } # end if if(strat==T)
 
 if(strat!=T){
  Z=matrix(rnorm(n*k),n,k)
 } # end if if(strat==T)

 if(whiten==T){
  for(j in 1:k) Z[,j] = (Z[,j]-mean(Z[,j]))/sd(Z[,j])
  Z=(Z%*%solve(chol(cov(Z))))%*%(chol(ZCOR))
 } # end if(whiten=T)
 
 Z
############################################################################### 
} # end function rMVN_Z
###############################################################################

###############################################################################
Enorm=function(x) {sqrt(sum(x^2))}  # computes Euclidian norm of vector x
Escale=function(x,R=1){(R/Enorm(x))*x} # resales x to have Enorm of R
##############################################################################


#############################################################################
conplot=function(A,b,cols=rep(1,nrow(A)),xlim=c(0,10),ylim=c(0,10),
        lwd=3,xlab='',ylab='',main=''){
#############################################################################
 if(ncol(A)>2) stop('This function only plots 2 dimensional constraints')

 if(length(cols==1)) cols=rep(cols,(nrow(A)))
 plot(0,0,type='n',xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main) 
 abline(v=0,lwd=lwd)
 abline(h=0,lwd=lwd)
 i=1
 for(i in 1:nrow(A)){
  if(A[i,2]!=0){
   m=-A[i,1]/A[i,2]
   pt=c(0,b[i]/A[i,2])
   pointslope(pt,m,lwd=3)
  }
  if(A[i,2]==0&A[i,1]!=0){
   xcept=b[i]/A[i,1]
   abline(v=xcept,lwd=3)
  }
 } # end loop on i
} #end function conplot
###############################################################################

###############################################################################
conplot2=function(A2,b2,cols=rep(2,nrow(A)),lwd=3,lty=1){
###############################################################################
         if(nrow(A2)!=length(b2)){
        (nr2=nrow(A2)); (nb2=length(b2))
        if(nr2<length(b2)){
                tmp=matrix(0,length(b2),ncol(A2))
                for(i in 1:length(b2)) tmp[i,]=A2[1,]
                A2=tmp
                if(length(cols)!=length(b2)) cols=rep(cols[1],nb2)
        }
        if(nr2>length(b2)) b2=rep(b2[1],nr2)
 } # end  if(nrow(A2)!=length(b2))

i=1
  for(i in 1:nrow(A2)){
    if(A2[i,2]!=0){
      (m=-A2[i,1]/A2[i,2])
      (pt=c(0,b2[i]/A2[i,2]))
      pointslope(pt,m,col=cols[i],lwd=lwd,lty=lty)
    }
    if(A2[i,2]==0&A2[i,1]!=0){
      xcept=b2[i]/A2[i,1]
      abline(v=xcept,col=cols[i],lwd=lwd,lty=lty)
    }
  } # end loop on i
###############################################################################
} #end function conplot2
###############################################################################

###############################################################################
objplot=function(obj,objval=NA,pt=c(0,0),cols=3,lwd=3,lty=2,arrow=T,
              length=0.15,R=1){
###############################################################################
if(!is.na(objval)){
 i=1
  if(obj[2]!=0){
   m=-obj[1]/obj[2]
   pt0=c(0,objval/obj[2])
   pointslope(pt0,m,col=cols,lwd=lwd,lty=lty,R=R)
  }
  
  if(obj[2]==0&obj[1]!=0){
   xcept=objval/obj[1]
   abline(v=xcept,col=cols,lwd=lwd,lty=lty)
  }
} # end if(!is.na(objval))
# 
# if(is.na(objval)){
#       
#  if(obj[2]!=0){
#   (m=-obj[1]/obj[2])
#   pointslope(pt,m=m,col=cols,lwd=lwd,lty=lty)
#   if(arrow==T&m!=0) {
#       pointslope(pt,m=-1/m,arrow=T,length=0.15,col=3,lwd=lwd,lty=2,R=R)
#   }   
#   if(arrow==T&m==0) {
#       arrows(pt[1],pt[2],y1=pt[2]+R,length=length,col=cols,lty=1,lwd=lwd)
#   } 
#  } # end if(obj[2]!=0)
#  
#  if(obj[2]==0){
#   abline(v=pt[1],col=cols,lwd=lwd,lty=lty)
#   if(arrow==T) arrows(pt[1],pt[2],x1=pt[1]+R,length=length,col=cols,lty=1,lwd=lwd)      
#  } # end if(obj[2]!=0)
# } # end if(is.na(objval))

if(is.na(objval)){
 if(obj[2]!=0){
  (m=-obj[1]/obj[2])
  pointslope(pt,m=m,col=cols,lwd=lwd,lty=lty)
 } # end if(obj[2]!=0)
 
 if(obj[2]==0){
  abline(v=pt[1],col=cols,lwd=lwd,lty=lty)
 } # end if(obj[2]!=0)
        
 if(arrow==T) {
          objR=Escale(obj,R=R)
          pt2=pt+objR
        arrows(pt[1],pt[2],x1=pt2[1],y1=pt2[2],length=length,col=cols,lty=1,lwd=lwd)
  } 
} # end if(is.na(objval))
###############################################################################
} #end function objplot
###############################################################################



###############################################################################
newseed=function(seedval){
set.seed(seedval)
(u=runif(1))
(newseed=floor(u[1]*2147483645))
if(newseed==0) newseed=1
newseed
} # end function newseed
###############################################################################


###############################################################################
# functions for specialized pairs plots
###############################################################################
pairsrank = function(x) {
  rx = rank(x, na.last = 'keep',ties.method='first')
  rx / max(rx, na.rm = T)
} # end function pairsrank
##################################
ecop=function(x,y,col = 1,pch = 1,cex = 1) {
  usr = par("usr")
  on.exit(par(usr))
  par(usr = c(-0.1, 1.1,-0.1, 1.1))
    points(pairsrank(x),pairsrank(y),pch = pch,cex = cex,col = col)
} # end function ecop
##################################
ecop_y=function(x,y,col = 1,pch = 1,cex = 1) {
  usr = par("usr")
  on.exit(par(usr))
  par(usr = c(-0.1, 1.1,-0.1, 1.1))
  points(pairsrank(y),pairsrank(x),pch = pch,cex = cex,col = col)
} # end function ecop_y
##################################
ecop2=function(x,y,col = 1,pch = 1,cex = 1) {
  usr = par("usr")
  on.exit(par(usr))
  par(usr = c(-0.1, 1.1,-0.1, 1.1))
  points(pairsrank(x),pairsrank(y),pch = pch,cex = cex,col = col)
  abline(0,1,col=1,lwd=3)

  scor=round(cor(x,y,method='spearman',use='pairwise.complete.obs'),3)
  mtext(paste('SCOR =',scor),cex=0.75,font=2)

  x2=x[pairsrank(y)<0.5];y2=y[pairsrank(y)<0.5]
  scor2=round(cor(x2,y2,method='spearman',use='pairwise.complete.obs'),2)
  #text(0.90,0.1,paste('',scor2),cex=1.5,font=2)
  #text(0.10,0.975,paste('',scor2),cex=1.5,font=2)

  x3=x[pairsrank(y)>0.5];y3=y[pairsrank(y)>0.5]
  scor3=round(cor(x3,y3,method='spearman',use='pairwise.complete.obs'),2)
  #text(0.10,0.95,paste('',scor3),cex=1.5,font=2)
  #text(0.9,0.10,paste('',scor3),cex=1.5,font=2)

} # end function ecop2
##################################
ecop2_y=function(x,y,col = 1,pch = 1,cex = 1) {
  usr = par("usr")
  on.exit(par(usr))
  par(usr = c(-0.1, 1.1,-0.1, 1.1))
  points(pairsrank(y),pairsrank(x),pch = pch,cex = cex,col = col)
  abline(0,1,col=1,lwd=3)
  
  scor=round(cor(x,y,method='spearman',use='pairwise.complete.obs'),3)
  mtext(paste('SCOR =',scor),cex=0.75,font=2)
  
  x2=x[pairsrank(y)<0.5];y2=y[pairsrank(y)<0.5]
  scor2=round(cor(x2,y2,method='spearman',use='pairwise.complete.obs'),2)
  #text(0.90,0.1,paste('',scor2),cex=1.5,font=2)
  #text(0.10,0.975,paste('',scor2),cex=1.5,font=2)
  
  x3=x[pairsrank(y)>0.5];y3=y[pairsrank(y)>0.5]
  scor3=round(cor(x3,y3,method='spearman',use='pairwise.complete.obs'),2)
  #text(0.10,0.95,paste('',scor3),cex=1.5,font=2)
  #text(0.9,0.10,paste('',scor3),cex=1.5,font=2)
  
} # end function ecop2
##################################
ecop3=function(x,y,col = 1,pch = 1,cex = 1) {
  usr = par("usr")
  on.exit(par(usr))
  par(usr = c(-0.1, 1.1,-0.1, 1.1))
  points(pairsrank(x),pairsrank(y),pch = pch,cex = cex,col = col)

  scor=round(cor(x,y,method='spearman',use='pairwise.complete.obs'),2)
  mtext(paste('SCOR =',scor),cex=0.5,font=2)

  x2=x[pairsrank(y)<0.5];y2=y[pairsrank(y)<0.5]
  scor2=round(cor(x2,y2,method='spearman',use='pairwise.complete.obs'),2)
  #text(0.90,0.1,paste('',scor2),cex=1.5,font=2)
  #text(0.10,0.975,paste('',scor2),cex=1.5,font=2)

  x3=x[pairsrank(y)>0.5];y3=y[pairsrank(y)>0.5]
  scor3=round(cor(x3,y3,method='spearman',use='pairwise.complete.obs'),2)
  #text(0.10,0.95,paste('',scor3),cex=1.5,font=2)
  #text(0.9,0.10,paste('',scor3),cex=1.5,font=2)

} # end function ecop3
###############################################################################
ecop3_y=function(x,y,col = 1,pch = 1,cex = 1) {
  usr = par("usr")
  on.exit(par(usr))
  par(usr = c(-0.1, 1.1,-0.1, 1.1))
  points(pairsrank(y),pairsrank(x),pch = pch,cex = cex,col = col)
  
  scor=round(cor(x,y,method='spearman',use='pairwise.complete.obs'),2)
  mtext(paste('SCOR =',scor),cex=0.5,font=2)
  
  x2=x[pairsrank(y)<0.5];y2=y[pairsrank(y)<0.5]
  scor2=round(cor(x2,y2,method='spearman',use='pairwise.complete.obs'),2)
  #text(0.90,0.1,paste('',scor2),cex=1.5,font=2)
  #text(0.10,0.975,paste('',scor2),cex=1.5,font=2)
  
  x3=x[pairsrank(y)>0.5];y3=y[pairsrank(y)>0.5]
  scor3=round(cor(x3,y3,method='spearman',use='pairwise.complete.obs'),2)
  #text(0.10,0.95,paste('',scor3),cex=1.5,font=2)
  #text(0.9,0.10,paste('',scor3),cex=1.5,font=2)
  
} # end function ecop3
###############################################################################

###############################################################################
panel.hist = function(x, ...)
{
  usr = par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  h = hist(x, plot = FALSE)
  breaks = h$breaks
  nB = length(breaks)
  y = h$counts
  y = y / max(y)
  rect(breaks[-nB], 0, breaks[-1], y,  ...)
} # end function panel.hist
###############################################################################

###############################################################################
panel.hist1 = function(x, ...)
{
  usr = par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5))
  h = hist(x, plot = FALSE)
  breaks = h$breaks
  nB = length(breaks)
  y = h$counts
  y = y / max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = 1, ...)
} # end function panel.hist
###############################################################################


###############################################################################
# demonstrate the ecop pairs plots
###############################################################################
# palette('R3')
# #############################################################################
# set.seed(1001)
# n=100; K=5
# X=matrix(rnorm(n*K),n,K)
# 
# (rholist=seq(0.9,0,-0.1))
# rholist[10]=-0.9
# (R=rho2cor(rholist))
# eigen(R)$values
# tmp=ImanConoverYcZc(X,R)
# cov2cor(tmp$sigma)
# Y=tmp$Yc
# Y[1,1]=10; Y[1,2]=-10
# Y[2,2]=10; Y[2,3]=-10
# Y[3,3]=10; Y[3,4]=-10
# head(Y)
# for(j in 1:K) Y[,j] = 10*j*Y[,j]
# head(Y)
# 
# colors=rep(1,n)
# colors[1]=2
# colors[2]=3
# colors[3]=4
# 
# #############################################################################
# # Create the plots
# #############################################################################
# graphics.off()
# my_x11()
# pairs(Y,pch=20,cex=2,diag.panel=panel.hist,col=colors)
# #############################################################################
# my_x11()
# pairs(Y,lower.panel=ecop,pch=20,cex=2,main='lower diag = ecop',col=colors)
# #############################################################################
# my_x11()
# pairs(Y,lower.panel=ecop_y,pch=20,cex=2,main='ecop_y',col=colors)
# #############################################################################
# my_x11()
# pairs(Y, lower.panel = ecop2,pch=20,cex=2,main='ecop2',col=colors)
# #############################################################################
# my_x11()
# pairs(Y, lower.panel = ecop2_y,pch=20,cex=2,main='ecop2_y',col=colors)
# #############################################################################
# my_x11()
# pairs(Y, lower.panel = ecop3,pch=20,cex=2,main='ecop3',col=colors)
# #############################################################################
# my_x11()
# pairs(Y, lower.panel = ecop3_y,pch=20,cex=2,main='ecop3_y',col=colors)
# #############################################################################
# my_x11()
# pairs(Y, lower.panel = ecop2,pch=20,cex=2,diag.panel=panel.hist,main='ecop2')
# #############################################################################
# my_x11()
# pairs(Y,lower.panel = ecop2_y,pch=20,cex=2,diag.panel=panel.hist,main='ecop2_y')
# #############################################################################
# my_x11()
# pairs(Y, lower.panel = ecop3,pch=20,cex=2,diag.panel=panel.hist1,main='ecop3')
# #############################################################################
# my_x11()
# pairs(Y,lower.panel=ecop3_y,pch=20,cex=2,diag.panel=panel.hist1,main='ecop3_y')
# #############################################################################
#
# end domos for ecops pairs plots
#
###############################################################################


###############################################################################
my_dim=function(x){
  if(is.vector(x)) d=length(x)
  if(is.array(x)) d=dim(x)
  d
} # end function my_dim
###############################################################################

###############################################################################
my_grep=function(greplist,x,droplist=NULL){
  matches=NULL
  j=1
  for(j in 1:length(greplist)){
    matches=c(matches,grep(greplist[j],x)) 
  } # end loop on j
  
  if(!is.null(droplist)){
   j=1
   for(j in 1:length(droplist)){
   matches=matches[-grep(droplist[j],x[matches])] 
   } # end loop for(j in 1:length(droplist))    
  } # end  if(!is.null(droplist))
  matches  
} # end function my_grep
###############################################################################
# # Example use of my_grep
#  (cropnames=c("Corn","Sweet Corn","Pop Corn","Dbl Crop WinWht/Corn",
#     "Dbl Crop Oats/Corn","Dbl Crop Barley/Corn","Dbl Crop Corn/Soybeans",
#     "Dbl Crop Triticale/Corn","Soybeans","Dbl Crop WinWht/Soybeans",
#     "Dbl Crop Soybeans/Oats","Dbl Crop Barley/Soybeans",
#     "Dbl Crop Corn/Soybeans","Dbl Crop Soybeans/Cotton"))
# #
#  (keep="Corn")
# #
#  (drop=c("Sweet Corn","Pop Corn"))
# #
#  my_grep(keep,cropnames,drop)
# #
#  cropnames[my_grep(keep,cropnames,drop)]
###############################################################################

###############################################################################
my_source=function(files){
###############################################################################  
 i=1
 for(i in 1:length(files)){ 
  print(paste(files[i]))   
  source(files[i])
 }   
###############################################################################  
} # end function my_source
###############################################################################


###############################################################################
# for debugging
#sscale=0.20;YPROB=NULL;seedval=NA
#sscale=0.20;YPROB=1:nrow(Y);seedval=1001
###############################################################################
# historical copula binding approach
###############################################################################
my_HCOP=function(YSIM,Y,sscale=0.0,YPROB=NULL,seedval=NA){
###############################################################################
if(ncol(YSIM)!=ncol(Y)) stop('Error: ncol(YSIM)!=ncol(Y) ')
###############################################################################
#YSIM=as.data.frame(YSIM)
Y=as.data.frame(Y)
###############################################################################
if(!is.na(seedval)) set.seed(seedval)
############################################################################### 
(nc=ncol(Y))
(n=nrow(Y))
(N=nrow(YSIM))
###############################################################################
Y$index=1:nrow(Y)
###############################################################################  
  

###############################################################################
# generate historical copula
###############################################################################
if(is.null(YPROB)){
 (nreps=floor(N/n)+1)
 pick=rep(1:n,nreps)
 HCOP=Y[pick,]
 HCOP=HCOP[1:N,]
}
###################
if(!is.null(YPROB)){    
 if(length(YPROB)!=nrow(Y)) stop('ERROR: length(YPROB)!=nrow(Y)')
 picky=sample(x=1:nrow(Y),size=N,prob=YPROB,replace=T)
 HCOP=Y[picky,]
}               
###################
# add a little white noise for tie breaking
for(j in 1:nc) HCOP[,j]=HCOP[,j]+rnorm(N,0,1E-6)
###############################################################################
# "smooth" the historical copula
###############################################################################
(sdevs=apply(HCOP,2,sd))
for(j in 1:nc) HCOP[,j] = HCOP[,j] + sscale*rnorm(N,0,sdevs[j])
###############################################################################
# put simulated marginals in HCOP rank order
###############################################################################
j=1
for(j in 1:nc){
 y=sort(YSIM[,j])
 YSIM[,j]=y[my_rank(HCOP[,j])] 
}
index=HCOP$index
###############################################################################
list(YSIM=YSIM,HCOP=as.matrix(HCOP[,1:nc]),index-index)
###############################################################################
} # end function HCOP
###############################################################################


##############################################################################
## debug on for my_heat function 
# xlab='LogSI';ylab='LogCII';ncells=75;pch='.';cex=0.75
# ncoLs=25;coL1='blue';coL2='yellow';coL3='red'
# Nlgd=15;Lround=0;scoremin=-(10^10);scoremax=10^10
# pdfon=F;pdfname='pdfplot.pdf';main=NULL; newplot=F
###############################################################################
my_heat=function(x1,x2,xlab=NULL,ylab=NULL,ncells=30,pch='.',cex=0.75,
ncoLs=25,coL1='blue',coL2='yellow',coL3='red',
Nlgd=11,Lround=0,scoremin=-(10^10),scoremax=10^10,plot_legend=T,
pdfon=F,pdfname='pdfplot.pdf',main=NULL,newplot=F){
###############################################################################
if(newplot==T){
  if(pdfon!=T) x11()

  if(pdfon==T) pdf(pname)
#if(pdfon==T) png(pname)
} # end if if(newplot==T)
###############################################################################
if(plot_legend==T) par(mar=c(5,5,5,5))  
if(plot_legend==T) par(mar=c(5,5,5,6),xpd=T)
###############################################################################
# Heatmap
obs=1:length(x1)
###############################################################################
df00=as.data.frame(cbind(obs,x1,x2))
#str(df00)
#head(df00)
###############################################################################
#ncells=75
(w1=(max(x1)-min(x1))/ncells)
(w2=(max(x2)-min(x2))/ncells)
df00$cell1=floor((df00$x1-min(df00$x1))/w1)+1
df00$cell2=floor((df00$x2-min(df00$x1))/w2)+1
head(df00)


(mcell1=floor((mean(df00$x1)-min(df00$x1))/w1)+1)
(mcell2=floor((mean(df00$x2)-min(df00$x2))/w2)+1)

###############################################################################
df2=doBy::summaryBy(obs~cell1+cell2,data=df00,FUN=length)
head(df2)
summary(df2)
df2$obs.length=df2$obs.length+rnorm(nrow(df2),0,0.10)
summary(df2)


tmp=colorum.gplots(df2$obs.length,
   ncoLs=ncoLs,coL1=coL1,coL2=coL2,coL3=coL3,
   Nlgd=Nlgd,Lround=Lround,scoremin=scoremin,scoremax=scoremax)

df2$color=tmp$datacolors
head(df2)

dim(df00)
dim(df2)
df2=merge(df00,df2)
dim(df2)

plot(df2$x1,df2$x2,col=df2$color,pch=pch,cex=cex,main=main,
     xlab=xlab,ylab=ylab)

if(plot_legend==T){
legend('right',inset=c(-0.2,0.1),legend=tmp$Ltext,fill=tmp$Lcolors,cex=0.75)
}  
###############################################################################
} # end function my_heat
###############################################################################

###############################################################################
# Examples my_heat function
###############################################################################
###############################################################################
# rho = 0.75
# rho = 0
# rho = -0.75
# 
# set.seed(1001)
# n=100000
# X=matrix(rnorm(2*n),n,2)
# x11()
# #plot(X)
# R=rho2cor(rho)
# X=ImanConover(X,R)
# #plot(X)
# x=X[,1]
# y=X[,2]
# 
# x11()
# my_heat(x,y,pch='.',ncells=50,xlab='X',ylab='Y',main='NORM XY HEAT PLOT')
# 
# x11()
# my_heat(x,y,pch='.',ncells=50,xlab='X',ylab='Y',main='NORM XY HEAT PLOT',
#         plot_legend = F)
# 
# 
# x11()
# my_heat(x,y,pch='.',cex=5,ncells=50,xlab='X',ylab='Y',main='NORM XY HEAT PLOT')
# 
# 
# x11()
# my_heat(x,y,pch='.',cex=10,ncells=50,xlab='X',ylab='Y',main='NORM XY HEAT PLOT')
###############################################################################
# End Examples my_heat function
###############################################################################

###############################################################################
stop_quietly <- function() {
  opt <- options(show.error.messages = FALSE)
  on.exit(options(opt))
  stop()
}
###############################################################################

###############################################################################
my_sd=function(x,pop=T,...){
###############################################################################  
 tmp=sd(x,...)
 if(pop==T){
  n=length.na(x)
  tmp=sqrt((n-1)/n)*tmp
 }# end if(pop==T)
 tmp
############################################################################### 
} # end function my_sd
###############################################################################

###############################################################################
my_var=function(x,pop=T,...){
###############################################################################  
  tmp=var(x,...)
  if(pop==T){
    n=nrow(X)
    tmp=((n-1)/n)*tmp
  }# end if(pop==T)
  tmp
############################################################################### 
} # end function my_var
###############################################################################


###############################################################################
my_cov=function(x,pop=T,...){
###############################################################################  
  tmp=cov(x,...)
  if(pop==T){
    n=nrow(X)
    tmp=((n-1)/n)*tmp
  }# end if(pop==T)
  tmp
############################################################################### 
} # end function my_cov
###############################################################################

###############################################################################
my_cVaR=function(x,p,risk='U'){
###############################################################################  
  x=na.omit(x)
  CDF=pvalx_V(x,x)
  (qL=max(x[CDF<=p]))
  (qU=min(x[CDF>=p]))
  (cVaR_L = mean(x[CDF<p]))
  (cVaR_U = mean(x[CDF>p]))
  if(risk=='D') tmp=as.data.frame(list(p=p,VaR_L=qL,cVaR_L=cVaR_L))
  if(risk=='U') tmp=as.data.frame(list(p=p,VaR_U=qU,cVaR_U=cVaR_U))
  tmp
###############################################################################  
} # end function VaR_cVaR 
###############################################################################

###############################################################################
my_quantile=function(x,p){
###############################################################################  
  x=na.omit(x)  
  CDF=pvalx_V(x,x)
  (qL=max(x[CDF<=p])) 
}# end function my_quantile
###############################################################################



################################################################################
# for debugging
#plots=T;KEEP=NULL
################################################################################
cleanup=function(plots=T,KEEP=NULL){
# if(!is.null(KEEP)){
  (tmp=ls(envir = .GlobalEnv))
  if(!is.na(match('keep',tmp))) KEEP=c(KEEP,keep)
  if(length(KEEP)>0){
   (tmp=tmp[is.na(match(tmp,KEEP))])
   rm(list=tmp,envir = .GlobalEnv)
  } # end if(length(keep)>0)
# } #  if(!is.null(KEEP))
 if(plots==T) graphics.off()
} # end function cleanup
#####################################################
LoadFunction = function(file,...) {
  dots = match.call(expand.dots = FALSE)$...
  dots = sapply(dots, as.character)

  output = lapply(dots, function(x,file)
   {eval(parse(text=paste(x," = function(x) {0}",sep="")),envir = .GlobalEnv)
    suppressMessages(insertSource(file, functions=x))
    eval(parse(text=paste(x," = ",x,"@.Data",sep="")),envir = .GlobalEnv) },file=file)

} # end function LoadFunction
###############################################################################








###############################################################################
# Linear Programming (LP) Functions
###############################################################################


###############################################################################
# quantile regression function  
###############################################################################
# Note:  rq is the function call from the "quantreg" R package
###############################################################################
# debugon
#q=0.5;R=NULL;dirsR=NULL;rhsR=NULL;solver='clp'
###############################################################################
my_rq=function(y,X,q=0.5,R=NULL,dirsR=NULL,rhsR=NULL,solver='clp'){
###############################################################################
 (n=length(y)); (nrX=nrow(X)); (ncX=ncol(X))
 if(n!=nrX) stop('length(y) != nrow(X)')
 if(q<=0|q>=1) stop('q must satisfy 0<q<1')     
 if(is.null(colnames(X))) colnames(X)=paste('X',1:ncX,sep='')
 (xnames=colnames(X))
 
 A=cbind(X,-diag(n),diag(n))
 colnames(A)=c(xnames,paste('ep',1:n,sep=''),paste('en',1:n,sep=''))
 cnames=colnames(A)
 
 obj=c(rep(0,ncX),rep(q,n),rep(1-q,n))
 
 
  
 nrA=nrow(A); ncA=ncol(A)
 
 dirs = rep('=',n)
 rhs = y
 
 xlower = c(rep(-Inf,ncX),rep(0,n),rep(0,n))
 xupper = rep(Inf,ncA)
 
 if(!is.null(R)){
        nrR=nrow(R);ncR=ncol(R)
        if(ncR!=ncX) stop('ncol(R)!=ncol(X)')
        if(nrR!=length(dirsR)) stop('nrow(R)!=length(dirsR)')
        if(nrR!=length(rhsR)) stop('nrow(R)!=length(rhsR)')
        
        AR=cbind(R,matrix(0,nrR,n),matrix(0,nrR,n))
        
        A=rbind(A,AR)
        dirs=c(dirs,dirsR)
        rhs=c(rhs,rhsR)
 } # end  if(!is.null(R))
 
 
 LP=createLP()

  LP$nr=nrow(A)
  LP$nc=ncol(A)
  LP$obj=obj
  LP$A=A
  LP$dirs=dirs
  LP$rhs=rhs
  LP$xlower=xlower
  LP$xupper=xupper

  LP$cnames=cnames

  LPSOL=LPSOLVER(LP,solver=solver)
 
 (status=LPSOL$status)
 sol=LPSOL$solution
 names(sol)=LP$cnames
 (sol=sol[1:ncX])
 
 fitted = X %*% sol
 resids = y - fitted
 
 list(status=status,coeff=sol,fitted=fitted,resids=resids,
      proctime=LPSOL$proctime,LP=LP)
############################################################################### 
}  # end function my_rq
###############################################################################



###############################################################################
# LP modeling functions written by Joe Atwood
###############################################################################

###############################################################################
ccc=function(nm='X',index=NULL){
 tmp=nm
 if(!is.null(index)){
  #print(paste('length index = ',length(index)))
  (tmp=paste(tmp,'[',index[1],sep=''))
  if(length(index)>=2){
         for(j in 2:length(index)){
   tmp=paste(tmp,index[j],sep=',')
   } # end for
  } # end if    
  tmp=paste(tmp,']',sep='')
 } # end if(!is.null(index))
 tmp
} # end function ccc
###############################################################################
#ccc(index='Y')
#ccc(nm='Y',index=c('Y','Z'))
# ccc('Y',c(1,3,5))
################################################################################


###############################################################################
cvert=function(x) {
  as.character(x)
} # end function cvert
###############################################################################
#cvert(1:5) 
###############################################################################

###############################################################################
createLP=function(sense='min'){
 list(LPsense=sense,nr=0,nc=0,obj = NULL,A=NULL,dirs=NULL,rhs=NULL,
      cnames=NULL,vartypes=NULL,xlower=NULL,xupper=NULL, 
      nvblocks=0,vblocks=list(NULL),SMM='A',ZINDEX=F)
} # end function createLP
###############################################################################
# LP=createLP()
# LP
###############################################################################


###############################################################################
# debug on
#v='X';vindex=expand.grid(i=1);vtype='C';vlower=0.0;vupper=Inf
#
#v='X';vindex=expand.grid(i=1:4);vtype='C';vlower=0.0;vupper=Inf
#
#v='X';vindex=expand.grid(i=1:4,j=1:3);vtype='C';vlower=0.0;vupper=Inf
#vindex=doBy::orderBy(~i+j,data=vindex)
#
#v='X';vindex=expand.grid(i=1:4,j=1:3,k=1:2);vtype='C';vlower=0.0;vupper=Inf
#vindex=doBy::orderBy(~i+j+k,data=vindex)
#
#v='EPXCS';vindex=expand.grid(c=1:C,r=1:R,x=(XS+1):X)
#vtype='C';vlower=-Inf;vupper=Inf
#
# CNAMES=c('CN','SB'); RNAMES=c('IA','NE'); XNAMES=c('LAND','LABOR','FERT') 
# v='EPXCS';vindex=expand.grid(c=CNAMES,r=RNAMES,x=XNAMES)
# vtype='C';vlower=-Inf;vupper=Inf
# #
# vindex=expand.grid(c=CNAMES,r=RNAMES,x=XNAMES[3])
######
# vtype='C'; vupper=Inf
# v='EPND';vindex = expand.grid(c=1:C);vlower = -Inf
# vnames=NULL
# CROPS=c('CN','SB');vnames = expand.grid(c=CROPS) 
#
# vtype='C'; vupper=Inf
###############################################################################
namevarLP=function(v='X',vindex=expand.grid(i=1:5,j=1:3,k=1:2),
         vnames=NULL){
###############################################################################
 if(is.null(vnames)) vnames=vindex 
 vnames=as.matrix(cbind(vnames))
 for(j in 1:ncol(vnames)) vnames[,j]=cvert(vnames[,j])

 (dims=apply(vindex,2,ulong))
#   VAR=array(' ',dim=apply(vindex,2,ulong))
  VAR=array(' ',dim=apply(vindex,2,max))
    
  cnames=''
#  for(icell in 1:nrow(vindex)) cnames[icell]=ccc(v,vindex[icell,])
  for(icell in 1:nrow(vindex)) cnames[icell]=ccc(v,vnames[icell,])
  
   VAR[as.matrix(cbind(vindex))] = cnames                                               
   VAR
###############################################################################  
} # end function namevarLP
###############################################################################
#X=addvarLP(LP)
###############################################################################
# v='X'
# vindex=expand.grid(i=1:5,j=1:3,k=1:2)
# vnames=NULL
# vtype='C'
# vlower=-Inf
# vupper=Inf
#keep=NULL
# 
# (keep=QR>0)
###############################################################################
# v='X'
# vindex=NULL
# vnames=NULL
# vtype='C'
# vlower=-Inf
# vupper=Inf
# keep=NULL
###############################################################################
addvarLP=function(LP,v='X',vindex=NULL,vnames=NULL,vtype='C',
                       vlower=0,vupper=Inf,keep=NULL){
###############################################################################
 if(is.null(vindex)) {VAR=v; cnames=v}
 
 if(!is.null(vindex)){
        if(!is.data.frame(vindex)) vindex=as.data.frame(vindex)
        if(!is.null(vnames)&!is.data.frame(vnames)) vnames=as.data.frame(vnames)
        
        VAR=namevarLP(v=v,vindex=vindex,vnames=vnames)
   
  if(is.null(keep)) {cnames=as.vector(VAR)}
  
  if(!is.null(keep)){(cnames=VAR[keep])}
 } # end  if(!is.null(vindex))  
  
  (nc=length(cnames))
  LP$nc=LP$nc+nc
  LP$cnames=c(LP$cnames,cnames)
  LP$vartypes=c(LP$vartypes,rep(vtype,nc))
  LP$xlower=c(LP$xlower,rep(vlower,nc))
  LP$xupper=c(LP$xupper,rep(vupper,nc))
  LP$nvblocks=LP$nvblocks+1
  LP$vblocks[[LP$nvblocks]]=VAR
  LP$obj=rep(0,LP$nc)
  names(LP$obj)=LP$cnames
  names(LP$xlower)=LP$cnames
  names(LP$xupper)=LP$cnames    
  LP <<- LP
  VAR
###############################################################################  
} # end function addvarLP
###############################################################################


###############################################################################
# End LP modeling functions section written by Joe Atwood
###############################################################################





###############################################################################
# LP solving functions written by Joe Atwood
###############################################################################
###############################################################################
#debug on run off
#LP=LPM;solver=solver
###############################################################################
LPSOLVER=function(LP,solver='clp',printobj=F,time_limit=36000){
  if('cnames' %in% names(LP) == F) LP$cnames = paste('X',(1:LP$nc),sep='')
  if('rnames' %in% names(LP) == F) LP$rnames = paste('R',(1:LP$nr),sep='')
  if('vartypes' %in% names(LP) == F) LP$vartypes = rep('C',LP$nc)
  
  if(solver=='highs'|solver=='HiGHS')           tmpLPM=LP_highs(LP,time_limit=time_limit)  
  if(solver=='clp'|solver=='CLP')               tmpLPM=LP_clp(LP)
  if(solver=='glpk'|solver=='GLPK')             tmpLPM=LP_glpk(LP)
  if(solver=='lpsolve'|solver=='lpSolve')       tmpLPM=LP_lpSolve(LP)
  if(solver=='rglpk'|solver=='Rglpk')           tmpLPM=LP_Rglpk(LP)  
  if(solver=='rsym'|solver=='Rsym')             tmpLPM=LP_Rsymphony(LP,time_limit=time_limit)
  if(solver=='lindo'|solver=='LINDO')           tmpLPM=LP_Lindo(LP)
  if(solver=='gurobi'|solver=='GUROBI')         tmpLPM=LP_gurobi(LP)
  if(solver=='lindo_mi'|solver=='LINDO_MI')     tmpLPM=MILP_Lindo(LP)
  if(solver=='gurobi_mi'|solver=='GUROBI_MI')   tmpLPM=MILP_gurobi(LP)  
  
  if(printobj==T) print(paste(solver,'objval',tmpLPM$objval))
  tmpLPM
} # end function LPSOLVER
###############################################################################
my_seconds=function(){
    time1=proc.time()
    as.numeric(time1[3])
}
###############################################################################
#source('A2SM-scripts.R')
###############################################################################



###############################################################################
# HiGHS LP solver function
LP_highs=function(LP,time_limit=1200){
###############################################################################
#require(highs);require(Matrix)
###############################################################################
if(!"highs" %in% .packages()|!"Matrix" %in% .packages()) {
stop('The packages highs and Matrix must be installed to use the option solver = highs')  
}
###############################################################################
# Note: This code is set up to allow
# the user to measure the time used in
# different procedures. The following
#lines pull initial times from the system
###############################################################################
 time000=my_seconds()
 time00=my_seconds()
#######################################
#Convert LP matrix format (if necessary)
#######################################
 if(LP$SMM!='CRI'|LP$ZINDEX!=F) {
   LP=LP2LP(LP,SMM2='CRI',ZINDEX2=F)
 }
# Compute sparse matrix conversion time used
 SMtime=my_seconds()-time00
#######################################
if('vartypes' %in% names(LP) == F) LP$vartypes = rep('C',length(obj)) 
####################################### 
 
#######################################
# Construct lower and upper bounds on x
# If xbounds are not specified in LP
# we assume 0 <= x <= Inf
#######################################
 if('xlower' %in% names(LP) == F) LP$xlower = rep(0,length(obj))
 if('xupper' %in% names(LP) == F) LP$xupper = rep(Inf,length(obj))
#######################################
# Construct lower and upper bounds on Ax
# If rbounds are not specified in LP
# we construct them from the LP list's
# 'rhs' and 'dirs' vectors
#######################################
 vartypes=LP$vartypes
 xlower=LP$xlower
 xupper=LP$xupper
 xlower[vartypes=='B'|vartypes=='b']=0
 xupper[vartypes=='B'|vartypes=='b']=1
 vartypes[vartypes=='B'|vartypes=='b']='I'
 LP$vartypes=vartypes
 LP$xlower=xlower
 LP$xupper=xupper

 if('rlower' %in% names(LP) == F) {
  rlower     = LP$rhs
  rlower     = ifelse(LP$dirs=='<='|LP$dirs=='<',-Inf,rlower)
  rlower     = ifelse(LP$dirs=='>='|LP$dirs=='>',LP$rhs,rlower)
  LP$rlower  = rlower
 }

 if('rupper' %in% names(LP) == F) {
  rupper     = LP$rhs
  rupper     = ifelse(LP$dirs=='<='|LP$dirs=='<',LP$rhs,rupper)
  rupper     = ifelse(LP$dirs=='>='|LP$dirs=='>',Inf,rupper)
  LP$rupper  = rupper
 }

 (LPsense=LP$LPsense)
 
 (LPsense = ifelse(LPsense=='max'|LPsense=='MAX',T,F))
 # highs blows up if abs of any sparse matrix ra values are <1e-9  
 LP$ra[abs(LP$ra)<2e-9]=2e-9  
 
# ia=LP$ia
# ja=LP$ja
# ra=LP$ra

 ASM <- sparseMatrix(LP$ia,LP$ja,x=LP$ra)

if(!exists('threads')) threads=1 
#print(paste('highs threads =',threads))
time0=my_seconds()

tmp <- highs_solve(L = LP$obj, lower = LP$xlower, upper = LP$xupper,
        A = ASM, lhs = LP$rlower, rhs = LP$rupper,types=LP$vartypes,
        maximum=LPsense,control=list(threads=threads,time_limit=time_limit))

proctime=my_seconds()-time0

 # retrieve the results
 (status=tmp$status)
 (objval=tmp$objective_value)
 (solution=tmp$primal_solution)

 tmp2=tmp$solver_msg
 rcost=tmp2$col_dual
 lhs=tmp2$row_value
 dual=tmp2$row_dual

 LP$solution=solution
 LP$rcost=rcost
 LP$lhs=lhs
 LP$dual=dual
 
 
 TIME=my_seconds()-time000

 names(solution)=LP$cnames
 names(dual)=LP$rnames
 names(rcost)=LP$cnames
 names(lhs)=LP$rnames
 
 list(status=status,objval=objval,solution=solution,dual=dual,rcost=rcost,
      lhs=lhs,TIME=TIME,SMtime=SMtime,
  proctime=proctime,LP=LP)
###############################################################################
} # end function LP_highs
###############################################################################


###############################################################################
# clpAPI LP solver function
LP_clp=function(LP){
###############################################################################
#require(clpAPI)
if(!"clpAPI" %in% .packages()) {
stop('The package clpAPI must be installed and loaded to use option solver = clp')  
}
###############################################################################
# Note: This code is set up to allow
# the user to measure the time used in
# different procedures. The following
#lines pull initial times from the system
###############################################################################
 time000=my_seconds()
 time00=my_seconds()
###############################################################################
#Convert LP matrix format (if necessary)
###############################################################################
if(LP$SMM!='CMO'|LP$ZINDEX!=T){
  LP=LP2LP(LP,SMM2='CMO',ZINDEX2=T)
}  
###############################################################################
typecheck=c(grep('I',LP$vartypes),grep('B',LP$vartypes))
if(length(typecheck)>0) stop('clp cannot solve integer or binary LP problems') 
# Compute sparse matrix conversion time used
 SMtime=my_seconds()-time00
#######################################
# Construct lower and upper bounds on x
# If xbounds are not specified in LP
# we assume 0 <= x <= Inf
#######################################
 if('xlower' %in% names(LP) == F) LP$xlower = rep(0,length(obj))
 if('xupper' %in% names(LP) == F) LP$xupper = rep(Inf,length(obj))
#######################################
# Construct lower and upper bounds on Ax
# If rbounds are not specified in LP
# we construct them from the LP list's
# 'rhs' and 'dirs' vectors
#######################################

 if('rlower' %in% names(LP) == F) {
  rlower     = LP$rhs
  rlower     = ifelse(LP$dirs=='<='|LP$dirs=='<',-Inf,rlower)
  rlower     = ifelse(LP$dirs=='>='|LP$dirs=='>',LP$rhs,rlower)
  LP$rlower  = rlower
 }

 if('rupper' %in% names(LP) == F) {
  rupper     = LP$rhs
  rupper     = ifelse(LP$dirs=='<='|LP$dirs=='<',LP$rhs,rupper)
  rupper     = ifelse(LP$dirs=='>='|LP$dirs=='>',Inf,rupper)
  LP$rupper  = rupper
 }

 nr = LP$nr
 nc = LP$nc
 obj = LP$obj
 (sense = ifelse(LP$LPsense!='max'&LP$LPsense!='MAX',1,-1))

 ia=LP$ia
 ja=LP$ja
 ra=LP$ra

 xlower     = LP$xlower
 xupper     = LP$xupper
 rlower     = LP$rlower
 rupper     = LP$rupper
 dirs       = LP$dirs

lp <- initProbCLP()
 # supress output
 setLogLevelCLP(lp,0)
 # direction of optimization
 setObjDirCLP(lp, sense)  # 1 = min, -1 = max

 # load problem data
 loadProblemCLP(lp, ncols=nc, nrows=nr, ia=ia, ja=ja, ra=ra,
  lb=xlower, ub=xupper,obj_coef=obj,rlb=rlower,rub=rupper)

   scaleModelCLP(lp,3)

  # solve lp problem
 time0=my_seconds()
  solveInitialCLP(lp)
  #primalCLP(lp,ifValP=0)
  #dualCLP(lp,ifValP=0)
  proctime=my_seconds()-time0

 # retrieve the results
 (status=getSolStatusCLP(lp))
 (objval=getObjValCLP(lp))
 (solution=getColPrimCLP(lp))
 (dual=getRowDualCLP(lp))
 (rcost=getColDualCLP(lp))
 (lhs=getRowPrimCLP(lp))

  
  LP$solution=solution
  LP$rcost=rcost
  LP$lhs=lhs
  LP$dual=dual
 # remove problem
 delProbCLP(lp)

 TIME=my_seconds()-time000

 names(solution)=LP$cnames
 names(dual)=LP$rnames
 names(rcost)=LP$cnames
 names(lhs)=LP$rnames
 
 list(status=status,objval=objval,solution=solution,dual=dual,rcost=rcost,
      lhs=lhs,TIME=TIME,SMtime=SMtime,proctime=proctime,LP=LP)
###############################################################################
} # end function LP_CLP
###############################################################################


###############################################################################
LP_glpk=function(LP){
###############################################################################
#require(glpkAPI)
############################################################################### 
if(!"glpkAPI" %in% .packages()) {
stop('The package glpkAPI must be installed and loaded to use option solver = glpk')  
} 
############################################################################### 
time000=my_seconds()
time00=my_seconds()

# Check and convert LP matrix format if necessary
if(LP$SMM!='CRI'|LP$ZINDEX!=T){ 
 LP=LP2LP(LP,SMM2='CRI',ZINDEX2=F)
} 

 SMtime=my_seconds()-time00

 obj=LP$obj
 
 
 if('vartypes'  %in% names(LP) == F) LP$vartypes  = rep('C',length(obj))

 if('xlower' %in% names(LP) == F) LP$xlower = rep(0,length(obj))
 if('xupper' %in% names(LP) == F) LP$xupper = rep(Inf,length(obj))

 if('rlower' %in% names(LP) == F) {
 #if(is.null(LP$rlower)){
  rlower     = LP$rhs
  rlower     = ifelse(LP$dirs=='<='|LP$dirs=='<',-Inf,rlower)
  rlower     = ifelse(LP$dirs=='>='|LP$dirs=='>',LP$rhs,rlower)
  LP$rlower  = rlower
 }

 if('rupper' %in% names(LP) == F) {
 #if(is.null(LP$rupper)){
  rupper     = LP$rhs
  rupper     = ifelse(LP$dirs=='<='|LP$dirs=='<',LP$rhs,rupper)
  rupper     = ifelse(LP$dirs=='>='|LP$dirs=='>',Inf,rupper)
  LP$rupper  = rupper
 }

 (nr = LP$nr)
 (nc = LP$nc)
 obj = LP$obj
 (LPsense = LP$LPsense)

 ia=LP$ia
 ja=LP$ja
 ra=LP$ra
 (nnz=length(ra))

 vartypes   = LP$vartypes
 xlower     = LP$xlower
 xupper     = LP$xupper
 rlower     = LP$rlower
 rupper     = LP$rupper

 lp <- initProbGLPK()
 # direction of optimization
 if(LPsense=='min'|LPsense=='MIN') setObjDirGLPK(lp,GLP_MIN)
 if(LPsense=='max'|LPsense=='MAX') setObjDirGLPK(lp,GLP_MAX)

 # add rows and columns
 addRowsGLPK(lp, nr)
 addColsGLPK(lp, nc)

 setObjCoefsGLPK(lp, j=1:nc, obj_coef=obj)
# setColsBndsGLPK(lp, j=1:nc, lb=xlower, ub=xupper, type=vartypes)
 setColsBndsGLPK(lp, j=1:nc, lb=xlower, ub=xupper)
 setRowsBndsGLPK(lp, c(1:nr), rlower, rupper)

 # load constraint matrix
 loadMatrixGLPK(lp, nnz, ia, ja, ra)
 # suppress terminal output
 termOutGLPK(GLP_OFF)

 time0=my_seconds()
 # solve lp problem
 solveSimplexGLPK(lp)
 proctime=my_seconds()-time0


 # retrieve the results
 (status = getSolStatGLPK(lp))
 (objval = getObjValGLPK(lp))
 solution = getColsPrimGLPK(lp)
 dual=getRowsDualGLPK(lp)
 rcost=getColsDualGLPK(lp)
 lhs=getRowsPrimGLPK(lp) 
 
 
 LP$solution=solution
 LP$rcost=rcost
 LP$lhs=lhs
 LP$dual=dual
 
 # remove problem object
 delProbGLPK(lp)

 TIME=my_seconds()-time000

 names(solution)=LP$cnames
 names(dual)=LP$rnames
 names(rcost)=LP$cnames
 names(lhs)=LP$rnames
 
 list(status=status,objval=objval,solution=solution,
            dual=dual,rcost=rcost,lhs=lhs,
      TIME=TIME,SMtime=SMtime,proctime=proctime,LP=LP)
###############################################################################
} # end function LP_glpk
###############################################################################


###############################################################################
LP_gurobi=function(LP,Threads=1,Cutoff=NULL,time_limit=36000,
  params=list()){
###############################################################################
#require(gurobi)
############################################################################### 
if(!"gurobi" %in% .packages()) {
  stop('The package gurobi must be installed and loaded to use option solver = gurobi')  
}
###############################################################################
 time000=my_seconds()
 model <- list()

if(LP$SMM!='A') {
 LP=LP2LP(LP,SMM2='A',ZINDEX2=F)
}

 model$A          = LP$A
 model$obj        = LP$obj
 (model$modelsense = LP$LPsense)
 model$rhs        = LP$rhs
 model$lb         = LP$xlower
 model$ub         = LP$xupper
 GBsense          = LP$dirs
 GBsense=ifelse(GBsense=='<='|GBsense=='<','<',GBsense)
 GBsense=ifelse(GBsense=='>='|GBsense=='>','>',GBsense)
 GBsense=ifelse(GBsense=='=='|GBsense=='=','=',GBsense)
 #GBsense

 (model$sense      = GBsense)
 (model$vtype      = LP$vartypes)

 if(is.null(Cutoff)){
   Cutoff=ifelse(LP$LPsense=='min'|LP$LPsense=='MIN',1E30,-1E30)
 }


# params <- list(OutputFlag=0)
 params$OutputFlag=0
 params$Threads = Threads
 params$TimeLimit = time_limit

 time0=my_seconds()
 result <- gurobi(model, params)
 proctime=my_seconds()-time0

 (objval=result$objval)
 solution=result$x
 # Clear space
 rm(model, result, params)

 TIME=my_seconds()-time000


 list(objval=objval,solution=solution,
       TIME=TIME,proctime=proctime)
###############################################################################
} # end function LP_Gurobi
###############################################################################


###############################################################################
LP_Lindo=function(LP){
###############################################################################
#require(rLindo)
if(!"rLindo" %in% .packages()) {
  stop('The package rLindo must be installed and loaded to use option solver = lindo')  
}  
###############################################################################  
time000=my_seconds()
  
if(LP$SMM!='CMO'|LP$ZINDEX!='T') {
  LP=LP2LP(LP,SMM2='CMO',ZINDEX2=T)
}  
  
  #Create LINDO enviroment object
  rEnv <- rLScreateEnv()

  #Create LINDO model object
  rModel <- rLScreateModel(rEnv)

  #Disable printing log
  rLSsetPrintLogNull(rModel)

  time00=my_seconds()
#  tmp=Lindo_A2SM(LP$A)
#  tmp=A2SM(LP$A,SMM='CMO',ZINDEX=T)



  (SMtime=my_seconds()-time00 )
  #Define the model data

  (nr=as.integer(LP$nr))
  (nc=as.integer(LP$nc))
  ia=as.integer(LP$ia)
  ja=as.integer(LP$ja)
  ra=LP$ra

  (nNZ=length(ra))


  (LPsense <- ifelse(LP$LPsense=='max'|LP$LPsense=='MAX',LS_MAX,LS_MIN))

  obj=LP$obj

  Ltypes=LP$dirs
  Ltypes=ifelse(Ltypes=='<='|Ltypes=='<','L',Ltypes)
  Ltypes=ifelse(Ltypes=='>='|Ltypes=='>','G',Ltypes)
  Ltypes=ifelse(Ltypes=='=='|Ltypes=='=','E',Ltypes)
  (Ltype=paste(Ltypes,sep='',collapse=''))


  rhs=LP$rhs

  pdLower =  LP$xlower
  pdUpper =  LP$xupper

  if('vartypes' %in% names(LP) == F) LP$vartypes = rep('C',length(obj)) 

  (vartypes=LP$vartypes)
  (Lvartypes=paste(vartypes,sep='',collapse=''))


  rLSloadLPData(model=rModel,nCons=nr,nVars=nc,nObjSense=LPsense,
    dObjConst=0,padC=obj,padB=rhs,pszConTypes=Ltype,nAnnz=nNZ,
    paiAcols=ja,panAcols = NULL,padAcoef=ra,paiArows=ia,
    padL=pdLower,padU=pdUpper)

  #rLSgetErrorMessage(rEnv,2053)

  time0=my_seconds()
  #Solve the model
  rLSoptimize(rModel,LS_METHOD_FREE)
  proctime=my_seconds()-time0


  #Get solution information
  (status=rLSgetDInfo(rModel,LS_DINFO_POBJ)$ErrorCode)
  (objval=rLSgetDInfo(rModel,LS_DINFO_POBJ)$pdResult)
  (solution=rLSgetPrimalSolution(rModel)$padPrimal)
  (dual=rLSgetDualSolution(rModel)$padDual)
  (rcost=rLSgetReducedCosts(rModel)$padReducedCost)
  (slack=rLSgetSlacks(rModel)$padSlack)
  (lhs=ifelse((LP$dirs=='>='|LP$dirs=='>'),slack-rhs,rhs-slack))

  
  LP$solution=solution
  LP$rcost=rcost
  LP$lhs=lhs
  LP$dual=dual
 # #compute lhs directly
 # round(t(SM2A(LP) %*% solution),3)
  
  #Delete the model and environment
  rLSdeleteModel(rModel)
  rLSdeleteEnv(rEnv)

  time1=my_seconds()
  TIME=my_seconds()-time000

  names(solution)=LP$cnames
  names(dual)=LP$rnames
  names(rcost)=LP$cnames
  names(lhs)=LP$rnames
  
 list(status=status,objval=objval,solution=solution,dual=dual,rcost=rcost,
        lhs=lhs,TIME=TIME,SMtime=SMtime,proctime=proctime,LP=LP)
###############################################################################
} # end function LP_LINDO
###############################################################################

###############################################################################
LP_lpSolve=function(LP){
###############################################################################
#require(lpSolveAPI)
###############################################################################
if(!"lpSolveAPI" %in% .packages()) {
  stop('The package lpSolveAPI must be installed and loaded to use option solver = lpsolve')  
}
###############################################################################
time000=my_seconds()

if(LP$SMM!='CMO') LP=LP2LP(LP,SMM2='CMO',ZINDEX2=F)
if(LP$ZINDEX!=F)  LP=LP2LP(LP,SMM2='CMO',ZINDEX2=F)

(LPsense = LP$LPsense)

nr=LP$nr;nc=LP$nc;obj=LP$obj;ia=LP$ia;ja=LP$ja;ra=LP$ra
dirs=LP$dirs; rhs=LP$rhs
###############################################################################
# Construct lower and upper bounds on x
# If xbounds are not specified in LP
# we assume 0 <= x <= Inf
###############################################################################
if('xlower' %in% names(LP) == F) LP$xlower = rep(0,length(obj))
if('xupper' %in% names(LP) == F) LP$xupper = rep(Inf,length(obj))
xlower=LP$xlower;xupper=LP$xupper
###############################################################################
 if('vartypes'  %in% names(LP) == F) LP$vartypes  = rep('C',length(obj))
 vartypes=LP$vartypes
 vartypes[vartypes=='C']='real'
 vartypes[vartypes=='B']='binary'
 vartypes[vartypes=='I']='integer'
###############################################################################
 dirs[dirs=='==']='='
 dirs[dirs=='<']='<='
 dirs[dirs=='>']='>='
###############################################################################
 lp=make.lp(nrow=nr,ncol=nc)
 lp.control(lp,sense=LPsense)

 j=1
 for(j in 1:nc){
  (jdxs=(ja[j]):(ja[j+1]-1))
  set.column(lp,j,ra[jdxs],indices=ia[jdxs])  
 }
 
 set.objfn(lp,obj)
 set.constr.type(lp,dirs)
 set.rhs(lp,rhs)
 set.bounds(lp, lower = xlower, upper = xupper)

 for(j in 1:nc) set.type(lp,j,vartypes[j])

 # solve lp problem
 time0=my_seconds()
 (status=solve(lp))
 (proctime=my_seconds() - time0)

(objval=get.objective(lp))
(solution=get.variables(lp))
(duals=get.dual.solution(lp))
(LHS=get.constraints(lp))

 
 LP$solution=solution
 LP$lhs=LHS
 LP$dual=duals
 
# write.lp(lp,file='LP_MPS_DATA.mps',type='freemps',use.names=c(F,F))
 
 delete.lp(lp)
 
 TIME=my_seconds()-time000

 names(solution)=LP$cnames
 names(duals)=LP$rnames
 names(LHS)=LP$rnames
 
 list(status=status,objval=objval,solution=solution,
   duals=duals,LHS=LHS,TIME=TIME,proctime=proctime,LP=LP)
###############################################################################
} # end function LP_lpSolve
###############################################################################


###############################################################################
LP_Rglpk=function(LP){
###############################################################################
# require(Rglpk); require(Matrix)
###############################################################################
if(!"Rglpk" %in% .packages()|!"Matrix" %in% .packages()) {
  stop('The packages Rglpk and Matrix must be installed and loaded to use the option solver = rglpk')  
}
################################################################################
 time000=my_seconds()
 

 LP=LP2LP(LP,SMM2='CRI',ZINDEX2=F)
 ASM = sparseMatrix(LP$ia, LP$ja, x = LP$ra)
 
 
 (LPsense = ifelse(LP$LPsense=='max'|LP$LPsense=='MAX',T,F)) 
 
 dirs=LP$dirs
 dirs[dirs=='=']='=='
 LP$dirs=dirs
 
 nc=ncol(LP$A)

 if('xlower' %in% names(LP) == F) LP$xlower = rep(0,length(LP$obj))
 if('xupper' %in% names(LP) == F) LP$xupper = rep(Inf,length(LP$obj))
 

 bounds=list(lower=list(ind=1:LP$nc,val=LP$xlower),
             upper=list(ind=1:LP$nc,val=LP$xupper))

 
 time0=my_seconds()

 tmpLP= Rglpk_solve_LP(LP$obj, ASM, LP$dirs, LP$rhs, bounds, 
                 LP$vartypes, max=LPsense) 
 

 proctime=my_seconds()-time0


 (status=as.numeric(tmpLP$status))
 (objval=tmpLP$optimum)
 (solution=tmpLP$solution)
 rcost=tmpLP$solution_dual
 dual=tmpLP$auxiliary$dual
 lhs=tmpLP$auxiliary$primal
 
 LP$solution=solution
 LP$rcost=rcost
 LP$lhs=lhs
 LP$dual=dual
 
 names(solution)=LP$cnames
 names(dual)=LP$rnames
 names(rcost)=LP$cnames
 names(lhs)=LP$rnames
 
 TIME=my_seconds()-time0
 list(status=status,objval=objval,solution=solution,dual=dual,rcost=rcost,
      lhs=lhs,TIME=TIME,proctime=proctime,LP=LP)
###############################################################################   
} # end function LP_Rglpk
###############################################################################

###############################################################################
LP_Rsymphony=function(LP,time_limit=1200){
###############################################################################
# my_require('Rsymphony')
###############################################################################
if(!"Rsymphony" %in% .packages()) {
  stop('The package Rsymphony must be installed and loaded to use the option solver = Rsymphony')  
}
###############################################################################
 time000=my_seconds()
#######################################
nc=LP$nc; nr=LP$nr 
#######################################
(LPsense = ifelse(LP$LPsense=='max'|LP$LPsense=='MAX',T,F))
######################################
#Convert LP matrix format (if necessary)
#######################################
if(LP$SMM!='CRI'|LP$ZINDEX!=F) {
 LP=LP2LP(LP,SMM2='CRI',ZINDEX2=F)
}
#######################################
# Compute sparse matrix conversion time used
SMtime=my_seconds()-time000
#######################################
# convert to sparse matrix triplet
#######################################
ASM = sparseMatrix(LP$ia,LP$ja,x=LP$ra)
#######################################
# put dirs in Rsym format
#######################################
dirs=LP$dirs
dirs[dirs=='=']='=='
dirs[dirs=='<']='<='
dirs[dirs=='>']='>='
LP$dirs=dirs
######################################
# Create Rsym lower_upper bounds object
#####################################
if('xlower' %in% names(LP) == F) LP$xlower = rep(0,length(obj))
if('xupper' %in% names(LP) == F) LP$xupper = rep(Inf,length(obj))
lower=list(ind=1:nc,val=LP$xlower)
upper=list(ind=1:nc,val=LP$xupper)
bounds=list(lower=lower,upper=upper)
#####################################
time0=my_seconds()

#tmpLP=Rsymphony_solve_LP(obj=LP$obj, mat=LP$A, dir=LP$dirs,
#     rhs=LP$rhs,types=LP$vartypes,max = LPsense,bounds = bounds) 
 
 
tmpLP=Rsymphony_solve_LP(obj=LP$obj, mat=ASM, dir=LP$dirs,
  rhs=LP$rhs,types=LP$vartypes,max = LPsense,bounds = bounds,
  time_limit=time_limit) 
 
proctime=my_seconds()-time0

(status=as.numeric(tmpLP$status))
(objval=tmpLP$objval)
(solution=tmpLP$solution)
TIME=my_seconds()-time0
###############################################################################
list(status=status,objval=objval,solution=solution,
      TIME=TIME,proctime=proctime)
###############################################################################
} # end function LP_Rsymphony
###############################################################################


###############################################################################
MILP_Gurobi=function(LPB,Threads=1,Cutoff=NULL,
    TimeLimit=1200,params=list()){
###############################################################################
#require(gurobi)
###############################################################################
if(!"gurobi" %in% .packages()) {
  stop('The package gurobi must be installed and loaded to use the function MILP_Gurobi')  
}
################################################################################
 time000=my_seconds()
 model <- list()

 if(LP$SMM!='A') {
  LP=LP2LP(LP,SMM2='A',ZINDEX2=F)
 }

 model$A          = LPB$A
 model$obj        = LPB$obj
 (model$modelsense = LPB$LPsense)
 model$rhs        = LPB$rhs
 model$lb         = LPB$xlower
 model$ub         = LPB$xupper
 GBsense          = LPB$dirs
 GBsense=ifelse(GBsense=='<='|GBsense=='<','<',GBsense)
 GBsense=ifelse(GBsense=='>='|GBsense=='>','>',GBsense)
 GBsense=ifelse(GBsense=='=='|GBsense=='=','=',GBsense)
 #GBsense

 (model$sense      = GBsense)
 (model$vtype      = LPB$vartypes)

 if(is.null(Cutoff)){
   Cutoff=ifelse(LPB$LPsense=='min'|LPB$LPsense=='MIN',1E30,-1E30)
 }


# params <- list(OutputFlag=0)
 params$OutputFlag=0
 params$Threads=Threads
 params$TimeLimit=TimeLimit

 time0=my_seconds()
 result <- gurobi(model, params)
 proctime=my_seconds()-time0


 (objval=result$objval)
 solution=result$x
 # Clear space
 rm(model, result, params)
 TIME=my_seconds()-time000
 list(objval=objval,solution=solution,
      TIME=TIME,proctime=proctime)
###############################################################################
} # end function MILP_Gurobi
###############################################################################


###############################################################################
MILP_Lindo=function(LPB,TimeLimit=1200){
###############################################################################
#require(rLindo)
###############################################################################
if(!"rLindo" %in% .packages()) {
  stop('The package rLindo must be installed and loaded to use the function MILP_Lindo')  
}
################################################################################
 
 time000=my_seconds()
# tmp=Lindo_A2SM(LPB$A)
# tmp=A2SM(LPB$A,SMM='CMO',ZINDEX = T)

 if(LPB$SMM!='CMO'|LPB$ZINDEX!='T'){
  LPB=LP2LP(LPB,SMM2='CMO',ZINDEX2=T)
 }

 SMtime=my_seconds()-time000

 #Create LINDO environment object
 rEnv <- rLScreateEnv()

 #Create LINDO model object
 rModel <- rLScreateModel(rEnv)

 #Disable printing log
 rLSsetPrintLogNull(rModel)


  #Define the model data

  (nr=as.integer(LPB$nr))
  (nc=as.integer(LPB$nc))
  ia=as.integer(LPB$ia)
  ja=as.integer(LPB$ja)
  ra=LPB$ra

 (nNZ=length(ra))
  (LPsense=LPB$LPsense)

 (LPsense <- ifelse(LPB$LPsense=='max'|LPB$LPsense=='MAX',LS_MAX,LS_MIN))

 obj=LPB$obj

 Ltypes=LPB$dirs
 Ltypes=ifelse(Ltypes=='<='|Ltypes=='<','L',Ltypes)
 Ltypes=ifelse(Ltypes=='>='|Ltypes=='>','G',Ltypes)
 Ltypes=ifelse(Ltypes=='=='|Ltypes=='=','E',Ltypes)
 (Ltype=paste(Ltypes,sep='',collapse=''))

 rhs=LPB$rhs

 pdLower =  LPB$xlower
 pdUpper =  LPB$xupper

 (vartypes=LPB$vartypes)
 (Lvartypes=paste(vartypes,sep='',collapse=''))

 rLSloadLPData(model=rModel,nCons=nr,nVars=nc,nObjSense=LPsense,
     dObjConst=0,padC=obj,padB=rhs,pszConTypes=Ltype,nAnnz=nNZ,
     paiAcols=ja,panAcols = NULL,padAcoef=ra,paiArows=ia,
     padL=pdLower,padU=pdUpper)

 #Load data to the model
 rLSloadVarType(rModel,Lvartypes)

 # load time limits
 rLSsetModelDouParameter(model=rModel,nParameter=LS_DPARAM_MIP_TIMLIM,dValue=TimeLimit)

 #rLSgetErrorMessage(rEnv,2053)


 #Solve the model
 time0=my_seconds()
 rLSsolveMIP(rModel)
 proctime=my_seconds()-time0


 #Get solution information
 (status=rLSgetIInfo(rModel,LS_IINFO_MIP_STATUS)$pnResult) # integer type'd info
 (objval=rLSgetDInfo(rModel,LS_DINFO_MIP_OBJ)$pdResult) # double type'd info

  xvals=rLSgetMIPPrimalSolution(rModel)$padPrimal


 #Delete the model and environment
 rLSdeleteModel(rModel)
 rLSdeleteEnv(rEnv)

TIME=my_seconds()-time000
list(status=status,objval=objval,xvals=xvals,
     TIME=TIME,SMtime=SMtime,proctime=proctime,LPB=LPB)
###############################################################################
} # end function MILP_LINDO
###############################################################################




###############################################################################
# debug on
#LP1=LP1;SMM2='CMO';ZINDEX2='T'
#LP1=LP2;SMM2='A';ZINDEX2='F'
#LP1=LP2;SMM2='CRI';ZINDEX2='F'
###############################################################################
LP2LP=function(LP1,SMM2='CMO',ZINDEX2=F){
###############################################################################
  LP2=LP1
  (SMM1=LP1$SMM)
  (ZINDEX1=LP1$ZINDEX)


  if(SMM1=='A'& SMM2!='A') {
    A2 = A2SM(LP1$A,SMM=SMM2,ZINDEX=ZINDEX2)
  }


  if(SMM1!='A'& SMM2!='A') {
    A1=list(nnz=LP1$nnz,nr = LP1$nr,nc = LP1$nc,ia = LP1$ia,ja = LP1$ja,
        ra = LP1$ra,SMM = SMM1,ZINDEX = ZINDEX1)
    A2 = SM2SM(A1,SMM2=SMM2,ZINDEX2=ZINDEX2)
  }


  if(SMM1!='A' & SMM2=='A'){
   A = SM2A(LP1$A)
  }


 if(SMM2 != 'A'){
  LP2$nnz=A2$nnz
  LP2$nr=A2$nr
  LP2$nc=A2$nc
  LP2$ia=A2$ia
  LP2$ja=A2$ja
  LP2$ra=A2$ra
  LP2$SMM=SMM2
  LP2$ZINDEX=A2$ZINDEX
  LP2=LP2[names(LP2) %in% c('A','ASM')==FALSE]
 }

 if(SMM2 == 'A'){
  LP2$A=A
  LP2$SMM=SMM2
  LP2=LP2[names(LP2)%in%c('nnz','nr','nc','ia','ja','ra','ASM','ZINDEX')==FALSE]
 }

  LP2
} # end function LP2LP
###############################################################################

###############################################################################
#' Write.LP: Writes LP to .csv file format
#'
#' @param LP      A sparse LP object
#' @param fname   .csv file name (include path if desired)
#' @return        = .csv file (can be opened/solved with spreadsheet solver)
###############################################################################
# debugging
#LP=LPB;fname='LPB.csv'
###############################################################################
Write.LP=function(LP,fname='LP.csv'){
###############################################################################
  if(length(grep('LP',search()))>0) detach(LP)
  if(exists('A')) rm(A,nr,nc)

	

  A = SM2A(LP)
  (nc=ncol(A))
  (nr=nrow(A))

	if('cnames' %in% names(LP) == F) LP$cnames=paste('C_',1:nc,sep='')
	if('rnames' %in% names(LP) == F) LP$rnames=paste('R_',1:nc,sep='')	
	
  if('rhs' %in% names(LP)==F&'rupper' %in% names(LP)==F){
   stop('rhs and (rlower,rupper )cannot both be missing in LP object')	
  }
  
 
 (LPsense=LP$LPsense)
 
 obj=LP$obj
# if('cnames' %in% names(LP) == F) LP$cnames = LP$cnames
# if('rnames' %in% names(LP) == F) LP$rnames = LP$rnames   


 if('vartypes' %in% names(LP) == F) LP$vartypes  = rep('C',length(obj))
 if('xlower' %in% names(LP) == F)   LP$xlower = rep(0,length(obj))
 if('xupper' %in% names(LP) == F)   LP$xupper = rep(1E10,length(obj))

 
 if('rlower' %in% names(LP) ==T){if(is.na(LP$rlower[1])) LP$rlower=NULL} 
 if('rupper' %in% names(LP) ==T){if(is.na(LP$rupper[1])) LP$rupper=NULL}
 
 if('rlower' %in% names(LP) == F) {
  rlower     = LP$rhs
  rlower     = ifelse(LP$dirs=='<='|LP$dirs=='<',-1E10,rlower)
  rlower     = ifelse(LP$dirs=='>='|LP$dirs=='>',LP$rhs,rlower)
  LP$rlower  = rlower
 }

 if('rupper' %in% names(LP) == F) {
  rupper     = LP$rhs
  rupper     = ifelse(LP$dirs=='<='|LP$dirs=='<',LP$rhs,rupper)
  rupper     = ifelse(LP$dirs=='>='|LP$dirs=='>',1E10,rupper)
  LP$rupper  = rupper
 }

 solution=0*obj
 if('solution' %in% names(LP) == T) solution=LP$solution 
 
  cnames=LP$cnames
  rnames=LP$rnames
  vartypes=LP$vartypes
  xlower=LP$xlower
  xupper=LP$xupper
  xlower[xlower<=-Inf] = -1E10
  xupper[xupper>=Inf] = 1E10  
 
  dirs=LP$dirs
  rhs=NA
  if('rhs' %in% names(LP) == T) rhs = LP$rhs
  rlower=LP$rlower
  rupper=LP$rupper

  
  M=A
  
  M=cbind(M,'|','lhsA','|')
  (tmp=c(cnames,'|','LHS','|'))

  ncol(M)
  length(tmp)
 
  
  
       
  if(exists('rhs')){
   M=cbind(M,dirs,'|',rhs,'|')	
   tmp=c(tmp,'dirs','|','rhs','|')
   ncol(M)
   length(tmp)
  }
  
  if(exists('rlower')){
    M=cbind(M,rlower,'|',rupper,'|')
   	tmp=c(tmp,'rlower','|','rupper','|')
    ncol(M)
    length(tmp)
  }
            
  M=cbind(rnames,M)
  tmp=c('',tmp)
    
  M=rbind(rep('========',ncol(M)),M)
  M=rbind(tmp,M)
  M=rbind(rep('========',ncol(M)),M)

  M=rbind(M,rep('========',ncol(M)))

  (tmp=c('OBJ',obj,'|','lhsOBJ','=',paste(LPsense)))
  (tmp=c(tmp,rep('',ncol(M)-length(tmp))))

  length(tmp);ncol(M)
  M=rbind(M,tmp)
  
  M=rbind(M,rep('========',ncol(M)))

  (tmp=c('nrows =',nr,'ncols =',nc))
  M=rbind(c(tmp,rep('',ncol(M)-length(tmp))),M)

  
  (tmp=c('xlower',xlower,'|','xlower'))
  M=rbind(M,c(tmp,rep('',ncol(M)-length(tmp))))
  (tmp=c('xupper',xupper,'|','xupper'))
  M=rbind(M,c(tmp,rep('',ncol(M)-length(tmp))))

  M=rbind(M,rep('========',ncol(M)))

#  (tmp=c('',0*obj,'|','xvalues'))
  (tmp=c('xvals',solution,'|','xvalues'))
  
  length(tmp)
  dim(M)
  M=rbind(M,c(tmp,rep('',ncol(M)-length(tmp))))

  M=rbind(M,rep('========',ncol(M)))
  M=rbind(rep('========',ncol(M)),M)

  (tmp=c('VARTYPE',vartypes,'|','vartypes'))
  
  M=rbind(M,c(tmp,rep('',ncol(M)-length(tmp))))
  M=rbind(M,rep('========',ncol(M)))
  #M=rbind(rep('========',ncol(M)),M)

 rm(A,nr,nc)
 write.table(M,file=fname,row.names=F,col.names=F,sep=',')
###############################################################################
} # end function Write.LP
###############################################################################



###############################################################################
# End Section  -- LP solving functions written by Joe Atwood
###############################################################################


