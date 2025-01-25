################################################################################
#
# organizeEvents function, reads in binaries denoted to include a bioacoustic
# event in PAMGuard, and organizes them by event.
#
# Emily T. Griffiths
# 2025
#
# emilytgriffiths@ecos.au.dk
#
################################################################################

organizeEvents=function(binaryPath, eventsdf) {
  
  # Required functions
  
  library(PamBinaries)
  library(PAMpal)
  library(plyr)
  library(dplyr)
  library(dbplyr)
  library(stringr)
  library(reshape2)
  library(lubridate)
  library(tictoc)
  library(tuneR)
  library(signal)
  library(seewave)
  library(RSQLite) 
  
  ## Variables ##
  
  # binaryPath    = Path that your binaries are located. NOTE: This script uses the standard terminology for binaries. If you have changed the names in PAMGUARD, you will need to adjust this script.
  # eventsdf      = A three column data.frame which includes the event IDs (e.g. 1, 2 3, etc.), and the event start/end datetimes. It assumes the order: eventID, mindate, maxdate.

  
  ## Check your variables ##
  
  # A path and database are necessary.
  if (missing(binaryPath)) {stop('No binaryPath variable set.')}
  if (missing(eventsdf)) {stop('No events included. Include a data frame that has eventIDs, min and max date.')}
  if (!is.data.frame(eventsdf)) {stop('Event tables needs to be a data.frame.')}
  eventsdf=as.data.frame(eventsdf)  #Format as a basic data frame, because Tibbles don't work. 
  
  binC=list.files(binaryPath, pattern = glob2rx('*Click_Detector_Clicks*.pgdf'), recursive = TRUE, full.names = TRUE) #The regex will need to change if you have altered the name of the detector in PG.

  bindf=as.data.frame(binC)
  bindf$name=basename(binC)
  bindf$datetime=as.POSIXct(gsub("[^\\d]+", "\\1",  bindf$name, perl = TRUE), format ='%Y%m%d%H%M%S', tz='UTC')
  
  for (ev in 1:nrow(eventsdf)) {
    dateRange <- seq(eventsdf[ev,2], eventsdf[ev,3], by="sec")
    
    #Which Bins
    mn=min(abs(difftime(bindf$datetime, min(dateRange))))
    mnDt=which(abs(difftime(bindf$datetime, min(dateRange)))==mn)
    
    mx=min(abs(difftime(bindf$datetime, max(dateRange))))
    mxDt=which(abs(difftime(bindf$datetime, max(dateRange)))==mx)
    
    dtidx=seq(mnDt, mxDt,by=1)
    # Include one extra bin on either side, unless it is the start or the end of the binary list. 
    if (min(dtidx)!=1) {dtidx= append(min(dtidx)-1, dtidx)}
    if (max(dtidx)!=nrow(bindf)) {dtidx=append(dtidx,max(dtidx)+1)}
    
    eventBins = binC[dtidx]  
    
    EvWaveData=list()
    wavData=list()
    EvClckData=NULL
    for (b in 1:length(eventBins)){
      binList=loadPamguardBinaryFile(binC[b])
      
      # Subset to event times
      realDates=lapply(binList$data, function(x) as.POSIXct((x$millis/1000), origin = '1970-01-01', tz = 'UTC'))
      withinEvent=which(realDates >=min(dateRange) & realDates<=max(dateRange))
      
      #Skip extra binary files.
      if (length(withinEvent)==0){
        next
      }
      
      clickData = NULL
      for (cld in withinEvent) {
        if (length(binList$data)==0) {break}
        data = as.data.frame(within(binList$data[[cld]],rm(wave, annotations)))

        # only take clicks classified by PG.
        if (as.numeric(data$type)==0){
          next
        }
        clickData = rbind.fill(clickData, data)
        
        #Extract the waveform
        UID=data$UID
        wave =list(binList$data[[cld]]$wave)
        names(wave)=UID
        wavData = c(wavData,wave)
      }
      
      clickData$datetime=as.POSIXct((clickData$millis/1000), origin = '1970-01-01', tz = 'UTC')
      clickData$eventID = paste0('Event_',eventsdf[ev,1])
      
      EvClckData=rbind(EvClckData, clickData)
      
    }
    EvWaveData[[ev]]=wavData
    
    evtsnms = unique(EvClckData$eventID)
    
    names(EvWaveData)=evtsnms
  }
  
  results = list(EvClckData, EvWaveData)
  
  return(results)
}
  

  

  
  