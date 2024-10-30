###############################################################################
# Start with an empty workspace
###############################################################################
rm(list = ls())                    # Remove all objects from workspace
###############################################################################
# The following code dynamically pulls the directory structure from computer.
# The code assumes you are in the main project folder.
###############################################################################
(workdir=paste(getwd(),'/',sep=''))
(datadir=paste(workdir,'data/',sep='')); dir.create(datadir)
(blockdir=paste(datadir,'blockdata/',sep='')); dir.create(blockdir)
(downdir=paste(datadir,'downloads/',sep='')); dir.create(downdir)
(ghcndir=paste(datadir,'ghcnd_all/',sep='')); dir.create(ghcndir)
(hcndir=paste(datadir,'ghcnd_hcn/',sep='')); dir.create(hcndir)
(prcpdir=paste(datadir,'prcp_data/',sep='')); dir.create(prcpdir)
(RDatadir=paste(datadir,'RData/',sep='')); dir.create(RDatadir)
(stationsdir=paste(datadir,'station_info/',sep='')); dir.create(stationsdir)
(tempdir=paste(datadir,'temp_data/',sep='')); dir.create(tempdir)
###############################################################################

source("Misc_functions.R")
source('my_GHCN_functions.R')
###############################################################################
# Source a function that let's us pull function from an external script
# The following functions will also allow you to work from a MAC in that
# the functions will check what operating system you are using and 
# accommodate the machine's operating system.
###############################################################################
#source('LoadFunction.R')
# load functions we will use in this session
#LoadFunction(Miscfile,cleanup,get_os,my_x11,my_require,wait,
    # my_over,PolyfindplotM,distcalc,my_lag,summaryDF,my_convert,
    # my_seconds,length.na,lagMat,jmonth,mo.yr,ulong)
###############################################################################
# Misc functions to directly download and process GHCN data
source('my_GHCN_functions.R')
###############################################################################
# Create vector of objects we want to keep when we "cleanup" the workspace.
###############################################################################
keep=ls(); (keep=c(keep,'keep'))
###############################################################################

###############################################################################
# Save workspace image in .RData file so that we do not need to rebuild each
# time.  Note if you change computers or put another copy of this project's R
# folder in a differnt location you will need to rebuild the workspace
# image WorkSpace.RData file.
###############################################################################
save.image('WorkSpace.RData')
###############################################################################
