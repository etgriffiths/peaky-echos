################################################################################
#
# Required functions and example script for the entire peaky-echos repository.
#
# Emily T. Griffiths
# 2025
#
# emilytgriffiths@ecos.au.dk
#
################################################################################

library(PamBinaries)
library(PAMpal)
library(stringr)
library(reshape2)
library(lubridate) 
library(plyr)
library(dplyr)
library(dbplyr)
library(tictoc)
library(tuneR)
library(signal)
library(seewave)
library(RSQLite)
library(spectral)

source("./_peaky-echos/decisionTree.R")
source("./_peaky-echos/organizeEvents.R")

path='C:/yourPath'
step1=decisionTree(path)  #Find dolphin events

#Create a list of dolphin events. This example has no supervision.
set.seed(1)
eventsdf = step1 %>%
  dplyr::group_by(eventID) %>%
  dplyr::summarise(minDate=min(step1$datetime), maxDate=max(step1$datetime)) 

step2=organizeEvents(path, eventsdf) #Extract binary information per event from relevant clicks.    

# Applying the broadband noise threshold. This is presented as a script, as there 
# are many variables that will require alteration depending on your project. 
# Decidecade noise levels are calculated using scripts found
# in https://github.com/ices-tools-prod/underwaternoise.

meanBB= 1.28498106978208e-07  # raw value
clipLevel= 175 #dB
f=384000 # Sample rate
loudThres=10*log10(meanBB)+clipLevel+40  #40 dB buffer is adjustable.

eventID = names(step2[[2]])

strongWaveform=list()
strongMD=NULL
for (ev in 1:length(eventID)) {
  
  if (is.null(meanBB)) {
    strongWaveform[[ev]]=list()
    next
  }
  
  event=step2[[2]][[ev]]
  
  if (length(event)==0){
    next
  }
  
  EvMaxPress= NULL
  for (w in 1:length(event)) {
    wav=event[[w]]
    ff=Wave(wav,samp.rate=f,bit=16)
    FFF=seewave::bwfilter(ff@left,f=ff@samp.rate,n=4,from=10000, to = 180000, bandpass=TRUE)
    sp = spec(FFF, f=f, wl=1024, wn='hanning',norm = FALSE, correction = "amplitude", plot=FALSE)
    sp[,2] = 10*log10(sp[,2])+clipLevel
    peakL = as.data.frame(cbind(names(event[w]),max(sp[,2])))
    EvMaxPress = rbind(EvMaxPress, peakL)
  }
  
  EvMaxPress[,2]=as.numeric(EvMaxPress[,2])
  
  clidx = EvMaxPress[EvMaxPress[,2]>=loudThres,]
  cond=step2[[1]]$UID %in% clidx[,1]
  
  evMD= subset(step2[[1]], cond)
  strongMD = rbind(strongMD, evMD)
  strongWaveform[[ev]] = event[cond]
  
}

names(strongWaveform)=eventID

step3 =list(strongMD, strongWaveform)

#Remove excess variables
rm(list=setdiff(ls(), c('path','step1','step2','step3')))
source("./_peaky-echos/ttsplit.R")
source("./_peaky-echos/concatWav.R")

step4 = ttsplit(step3,0.35)  # Split data into training and testing datasets.

#Concatente waveforms in training data.
f=384000 # Sample rate
bit=16  # Bit rate
allWavformTr = list()
eventID=names(step4$Waveforms$Training)
for (ev in 1:length(eventID)) {
  eventW = step4$Waveforms$Training[[ev]]
  UIDs = paste0(eventID[ev],'_',names(eventW))
  names(eventW)=UIDs
  allWavformTr=append(allWavformTr, eventW)
}

fullspec = concatWav(allWavformTr, f, bit)

cFFF=seewave::bwfilter(fullspec@left,f=f,n=4,from=20000, to = 180000, bandpass=TRUE)

# Generate a mean spectra for the training data with a normalized amplitude. FFT and settings are depenedent on your data.
sp = meanspec(cFFF, f=f, wl=512, wn='hanning',norm = TRUE,  plot=TRUE)

# subset you spectra to the frequencies to be used in your template. e.g., where the banding is.
min20 = which.min(abs(sp[,1]-20))
min80 = which.min(abs(sp[,1]-80))
ssp <- data.frame(sp[min20:min80, ])
names(ssp)=c('freq', 'amp')


# Set your starting parameters for the sinusoidal model. 
#I used the equation y(t) = A*sin(Omega*t + Phi) + C where A is the amplitude, Omega the period, Phi the phase shift and C the midline.
A<- (max(ssp$amp)-min(ssp$amp)/2)
C<-((max(ssp$amp)+min(ssp$amp))/2)
omega=pi/4
phi=1

# Create your model
res<- nls(amp ~ A*sin(omega*freq+phi)+C, data=ssp, start=list(A=A,omega=omega,phi=phi,C=C))
# Extract your coefficients.
co <- coef(res)
# Fit your model
fit <- function(x, a, b, c, d) {a*sin(b*x+c)+d}

plot(ssp, type='l')
title('Mean Spectra of Training Data, 20-80 kHz')
# Extract template curve from the plot
sintemp=curve(fit(x, a=co["A"], b=co["omega"], c=co["phi"], d=co["C"]), add=TRUE ,lwd=2, col="steelblue2")

# Interpolate the real data to have as many datapoints as the template.
interp_output<- approx(x=ssp, n=length(sintemp$y))

# Pearson correlation test, as in how well correlated is the template with the real data.
cor.test(interp_output$y, sintemp$y, method = 'pearson')

# Save your template to import as need be.
save(co, sintemp,fit,res, file='./_peaky-echos/templateInfo.rda')
rm(list=setdiff(ls(), c('path', 'step4', 'concatWav', 'f', 'bit')))

load('./_peaky-echos/templateInfo.rda')

eventID = names(step4$Waveforms$Testing)

# Concatenate all events individually within the testing data.
concatEvents = list()
for (ev in 1:length(eventID)) {
  event = step4$Waveforms$Testing[[ev]]
  concat=concatWav(event, f, bit)
  concatEvents=append(concatEvents,concat)
}
names(concatEvents)=eventID

## Fit the modeled template to the testing data.
sumData=list()
for (ev in 1:length(concatEvents)){
  cwav=concatEvents[[ev]]
  
  cFFF=seewave::bwfilter(cwav@left,f=f,n=4,from=20000, to = 80000, bandpass=TRUE)
  sp = meanspec(cFFF, f=f, wl=512, wn='hanning',norm = TRUE, plot=FALSE)
  
  min20 = which.min(abs(sp[,1]-20))
  min80 = which.min(abs(sp[,1]-80))
  
  ssp <- data.frame(sp[min20:min80, ])
  names(ssp)=c('freq', 'amp')

  # Fit modeled values to a curve based on the events data.  This keeps omega or ω,
  # the frequency spacing between the peaks in the spectra, as a fixed measurement.
  # The other variables start at their modeled values.
  resT<- nls(amp ~ A*sin(co['omega']*freq+phi)+C, data=ssp, 
           start=list(A=co['A'],C=co['C'],phi=co['phi']),control = nls.control( warnOnly = TRUE)) 
  
  # Extract coeffients for the event fit.
  coO<- coef(resT)
  
  # Extract your applied curve with the same fixed frequency spacing.
  plot(ssp, type='l')
  sinApp=curve(fit(x, a=coO[1], b=co["omega"], c=coO[3], d=coO[2]), add=TRUE ,lwd=2, col="darkorange")
  
  interp_fit<- approx(x=ssp, n=length(sinApp$y))
  corfit = cor.test(interp_fit$y, sinApp$y, method = 'pearson')

  #Save your data.
  sumData[[ev]]=list()
  sumData[[ev]]$res=summary(resT)
  sumData[[ev]]$corfit=corfit
}
  




