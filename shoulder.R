################################################################################
#
# Shoulder function
#
# Emily T. Griffiths
# 2018
#
# Original version of this function is available in:
# Griffiths, E. T., Archer, F., Rankin, S., Keating, J. L., Keen, E., Barlow, J., & Moore, J. E. (2020). Detection and classification of narrow-band high frequency echolocation clicks from drifting recorders. The Journal of the Acoustical Society of America, 147(5), 3511-3522.
#
# 2025 Update:
# Added lowpass variable option, instead of having it fixed within the function 
# for clicks with needed frequencies below 100 kHz.
#
# emilytgriffiths@ecos.au.dk
#
#
################################################################################


library(tuneR)
library(seewave)

## Variables ##

# spec          = spectrum, produced from the spec() function in Seewave
# f             = sample rate of the recording
# threshold     = Amplitude threshold that you want to find peaks above. Dependent on how the amplitude is assigned in spec()
# amp.slope     = Amplitude slope minimim. Again, dependent on how amplitude was assigned in spec().
# freq.dist     = Minimum distance between peaks in kHz, to avoid getting ripple peaks from jagged spectra on the same peak.
# lowpass       = Lowpass filter to apply to the data, good for NBHF clicks like kogia. In kHz
# plot          = T/F. Do you want a plot?

## Check your variables ##

shoulder=function(spec,f,threshold,amp.slope,freq.dist, lowpass, plot) {
  
  #Find peak Hz
  prim=fpeaks(spec,f=f,nmax=1, plot=F)
  #Find all peaks above your amplitude threshold.
  pks=fpeaks(spec,f=f, threshold = threshold,  amp=c(amp.slope,amp.slope), plot=F)
  #Subset those peaks to remove ones too close the actual peak.
  res=subset(pks,pks[,1] >= prim[1]+freq.dist | pks[,1] <= prim[1]-freq.dist)
  #Applies low pass filter
  if (!missing(lowpass)) {
    res=subset(res,res[,1]>=lowpass)
  }
  #Find Shoulder
  mx=which.max(res[,2])
  shod=res[mx,]
  
  #Plot
  b.pks=rbind(prim,shod)
  
  if (plot) {
    plot(spec, type = "l", xlab = "Frequency (kHz)", ylab = "Amplitde")
    title("Peak and Shoulder Frequencies")
    points(b.pks, col = "red")
    text(b.pks, labels = round(b.pks[, 1], digits = 2),  pos = 1, col = "red")
    abline(h = threshold, col = "red", lty = 2)
  }
  
  results=rbind(shod,prim)
  results
  
}    
  