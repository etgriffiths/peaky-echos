################################################################################
#
# ttsplit function, splits the truthed data into training and testing dataset, 
# based on a ratio provided by the user.
#
# Emily T. Griffiths
# 2025
#
# emilytgriffiths@ecos.au.dk
#
################################################################################




ttsplit=function(dataList, ra) {
  
  # Required functions
  library(seewave)
  library(plyr)
  library(dplyr)
  library(lubridate)
  
  ## Variables ##
  
  # dataList  = Output from organizeEvents function, which has both the metadata and wave files organized as a list.
  # ra        = Ratio of how you want to split your data, as in, how much of your data will be included in the training dataset? Default is 35%.
  
  ## Check your variables ##
  
  if (missing(dataList)) {stop('Please provide list with both metadata and waveforms.')}
  if (missing(ra)) {
    warning('No precentage of data to be allocated to training data specified. Set to 35%.')
    ra = .35
  }
  
  testingWaveform=list()
  trainingWaveform=list()
  
  
  classifiedClicks_testing=NULL
  classifiedClicks_training=NULL  
  
  eventNo = names(dataList[[2]])
  for (ev in 1:length(eventNo)) {
    
    eventD = dataList[[1]][dataList[[1]]$eventID==eventNo[ev],]
    eventW = dataList[[2]][[ev]]
    
    
    if (nrow(eventD) < 100) {
      
      classifiedClicks_testing=rbind.fill(classifiedClicks_testing,eventD)
      testingWaveform[[ev]] = eventW
      trainingWaveform[[ev]]=list()
      
    } else {
      
      trainingD=eventD %>%
        group_by(minute(datetime)) %>%
        sample_frac(ra)
      
      cond=eventD$UID %in% trainingD$UID
      
      testingD=subset(eventD, !cond)
      
      trainingW = eventW[cond]
      testingW = eventW[!cond]
      
      if (length(trainingW)!= nrow(trainingD)) {
        stop('There is an issue with your click formating, in that your number of clicks in the metadata does not match the number of clicks waveforms.')
      }
      
      
      
      classifiedClicks_testing=rbind.fill(classifiedClicks_testing,testingD)
      testingWaveform[[ev]] = testingW
      
      classifiedClicks_training=rbind.fill(classifiedClicks_training, trainingD)
      trainingWaveform[[ev]] = trainingW
    }
  }
    
    names(testingWaveform) = eventNo
    names(trainingWaveform) = eventNo
    
    md = list(classifiedClicks_testing, classifiedClicks_training)
    names(md) = c('Testing', 'Training')
    wv = list(testingWaveform, trainingWaveform)
    names(wv)= c('Testing', 'Training')
    
    results = list(md, wv)
    names(results) = c('Metadata', 'Waveforms')
    return(results)
  
}
  