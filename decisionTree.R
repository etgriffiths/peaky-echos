################################################################################
#
# decisionTree function, which reads in both the click and burst pulse binary 
# files from PAMGuard, and makes a determines if each second may have dolphin 
# bioacoustic activity based on the provided parameters.
#
# Emily T. Griffiths
# 2025
#
# emilytgriffiths@ecos.au.dk
#
################################################################################



decisionTree=function(binaryPath, interval, trainMin, perPosID, whistleMin, clickOnlyConf) {
  
  # Required functions
  
  require(PamBinaries)
  require(PAMpal)
  require(dplyr)
  require(stringr)
  require(reshape2)
  require(lubridate) 
  
  ## Variables ##
  
  # binaryPath    = Path that your binaries are located. NOTE: This script uses the standard terminology for binaries. If you have changed the names in PAMGUARD, you will need to adjust this script.
  # interval      = How many minutes must lapse before a new event starts. Default is 15.
  # trainMin      = Minimum number of click trains per second needed for it to be considered a train. Default is 5.
  # perPosID      = Percent of positive click classifications by PAMGuard per minute over the total number of clicks detected that second. Default is 30%.
  # whistleMin    = Minimum number of whistles/burst pulses detected per second to be considered a burst pulse. Default is 3.
  # clickOnlyConf = Percent of positive click classifications by PAMGuard per minute over the total number of clicks, if no burst pulses are detected. Default is 80%.
  

  ## Check your variables ##
  
  # A path is necessary.
  if (missing(binaryPath)) {stop('No binaryPath variable set.')}
  if (missing(interval)) {
    warning('No event interval set. Set to 15.')
    interval = 15
  }
  if (missing(trainMin)) {
    warning('No event minimum clicks per second (trainMin) set. Set to 5.')
    trainMin = 5
  }
  if (missing(perPosID)) {
    warning('No percent positive ID of clicks per second (perPosID). Set to 30%  (0.3).')
    perPosID = 0.3
  }
  if (missing(whistleMin)) {
    warning('No minimum number of burst pulses set (whistleMin). Set to 3.')
    whistleMin = 3
  }
  if (missing(clickOnlyConf)) {
    warning('No click only percent postive ID of clicks per second (clickOnlyConf). Set to 80% (0.8).')
    clickOnlyConf = 0.8
  }

  
  ## Clicks First ##
  
  binC=list.files(binaryPath, pattern = glob2rx('*Click_Detector_Clicks*.pgdf'), recursive = TRUE, full.names = TRUE) #The regex will need to change if you have altered the name of the detector in PG.
  
  clickEvents=NULL
  for (b in 1:length(binC)){
    
    #Load in click binary file
    binList=loadPamguardBinaryFile(binC[b])
    
    # Isolate click metadata from waveform and annotation data.
    clickData = NULL
    for (cld in 1:length(binList$data)){
      data = as.data.frame(within(binList$data[[cld]],rm(wave, annotations)))
      clickData = rbind(clickData, data)
    }
    
    # Group the clicks by the nearest second, and classification type. 
    clickData$dateSec = round(clickData$date, digits=0)
    groupedSecs=clickData %>% group_by(dateSec, type) %>% dplyr::summarize(count=n())
    
    # Select you classification, as assigned by PAMGuard. This is a setting you can name in PAMGuard, check your click detector for the ID.
    withClass = groupedSecs$dateSec[groupedSecs$type==1]
    if (length(withClass)==0) {next}
    
    ## NOTE ##
    # This script is set up for only one click classifier used in PAMGuard. Therefore 
    # 'type' is either 0 (no classification) or 1 (meets classification criteria). 
    # Modifications will be necessary if more than one classifier is employed in PG.
    ##
    
    # This next section ensures that multiple types were detected within the second,
    # as if only our click type was detected, that is likely the result of an error
    # or recorded noise/static. Remove for animals with a slow ICI or clean recordings.
    DateCounts=table(groupedSecs$dateSec) #All Dates
    withBoth=table(groupedSecs$dateSec)>1 #
    bar=DateCounts[withBoth]
    candidates=as.numeric(names(bar))
    
    candidates = unique(c(candidates, withClass))
    
    # How many clicks were classified, and how many were not.
    myData = groupedSecs[groupedSecs$dateSec %in% candidates,]
    if (nrow(myData)==0) {next}
    dat1=dcast(myData, dateSec ~type)
    names(dat1)=c('dateSec','No', 'Yes')
    dat1$totalClicks = rowSums(dat1[,c("No","Yes")], na.rm=TRUE)
    # Get the percentage of clicks positively classified.
    dat1$PerPos=dat1$Yes/dat1$totalClicks
    
    ## FIRST BRANCH ##
    # Does the data meet the minimum number of clicks necessary to be considered a 
    # train AND is the percentage of positively identified clicks exceed the threshold?
    idx = dat1$Yes >= trainMin & dat1$PerPos >= perPosID
    
    events=dat1[idx,]
    
    events$datetime = as.POSIXct(events$dateSec, origin = '1970-01-01 00:00:00', tz='UTC')
    
    clickEvents = rbind(clickEvents,events)
  }
  
  
  
  ## Whistles/Burst Pulses Second ##
  
  binBP=list.files(binaryPath, pattern = glob2rx('*Whistle_and_Moan*.pgdf'), recursive = TRUE, full.names = TRUE)  #The regex will need to change if you have altered the name of the detector in PG.
  
  #Exctract BP counts
  bpEvents=NULL
  for (b in 1:length(binBP)){
    binList=loadPamguardBinaryFile(binBP[b])
    
    bpData = NULL
    for (wh in 1:length(binList$data)){
      data = as.data.frame(within(binList$data[[wh]],rm(sliceData, contour, contWidth, annotations)))
      bpData = rbind(bpData, data)
    }
    
    bpData$dateSec = round(bpData$date, digits=0)
    woo=bpData %>% group_by(dateSec) %>% dplyr::summarize(burstPulse=n())
    
    #Apply burst pulse min decision tree branch.
    woosub = woo[woo$burstPulse>=whistleMin,]
    woosub$datetime = as.POSIXct(woosub$dateSec, origin = '1970-01-01 00:00:00', tz='UTC')
    
    bpEvents = rbind(bpEvents, woosub)
    
  
  }
  
  #Merge datasets.
  allDetections=merge(clickEvents,bpEvents, by = c('dateSec', 'datetime'), all = TRUE)
  
  # Apply last decision tree branch
  powerEvs = allDetections[allDetections$PerPos >=clickOnlyConf | allDetections$burstPulse >= whistleMin,]
  powerEvs=powerEvs[order(powerEvs$datetime),] 
  #Remove nulls.
  allEvents=powerEvs[complete.cases(powerEvs$dateSec),]
  
  # Create Events
  allEvents$eventID=''
  allEvents$eventID[1]= 1
  counter=1
  
  for (ee in 2:nrow(allEvents)) {
    #Create events based on the minimum time interval.
    timediff=as.numeric(difftime(allEvents$datetime[ee], allEvents$datetime[ee-1], units = 'mins'))
    if(timediff<interval) {
      allEvents$eventID[ee]=counter
    } else {
      counter = counter +1
      allEvents$eventID[ee]=counter
    }
    
  }
  
  return(allEvents)
}