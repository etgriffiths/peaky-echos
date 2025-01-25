################################################################################
#
# concatenateEvent function, Concatenates the all the clicks from an event into 
# a single waveform.
#
# Emily T. Griffiths
# 2025
#
# emilytgriffiths@ecos.au.dk
#
################################################################################



# Concatenate an event

concatWav = function(wavList, f, bit) {
  
  library(tuneR)
  
  ## Variables ##
  
  # wavList = List of wav spectras to concatenate.
  # f       = sample rate of the recording
  # bit     = bit rate of recording
  
  ## Check your variables ##
  
  if (missing(wavList)) {stop('Please provide list of waveforms/spectra.')}
  if (missing(f)) {stop('Please provide a wav sample rate.')}
  if (missing(bit)) {
    warning('No bit rate specified. Set to 16.')
    bit = 16
  }
  
  
  wav=wavList[[1]]
  concat=Wave(wav,samp.rate=f,bit=bit)
  for (w in 2:length(wavList)){
    wav=wavList[[w]]
    ff=Wave(wav,samp.rate=f,bit=bit)
    
    concat=bind(concat,ff)
  }
  return(concat)
}