
# ---
# title:MMGLib.R
# author: "Alice Rueda"
# date: "February 19, 2019"
# This file contains the functional calls for MMG analysis using accelerometers.
# ---

# Using signal information and pulse indications to detection events
detectMMGEvent <- function(dftmp, percentage = 0.08, ampFactor = 1){
  
  # Adding optional amplification factor to weak signals
  tmpx = dftmp$X.A1.X.*ampFactor
  tmpy = dftmp$X.A1.Y.*ampFactor
  tmpz = dftmp$X.A1.Z.*ampFactor
  
  # Finding threshold, setting at 5% of noise removed maximum magnitude
  tmp <- sqrt(tmpz^2+tmpy^2+tmpx^2)
  
  qplot(tmp)
  
  # Using low pass to remove high frequency noise
  Fs <- 1/0.02
  Fn <- Fs/2
  f <- 0.5      # in Hz
  per <-f/Fn
  bf <- butter(4, per, type="low")
  filtered <- filtfilt(bf, abs(tmp))
  threshold = max(filtered)*percentage
  
  # Event detection from signal
  #mavTremsenData <- movavg(filtered, 30, type="s")
  #eventDetect <- mavTremsenData
  eventDetect <- filtered
  eventDetect[filtered < threshold] <- 0 
  eventDetect[filtered >= threshold] <- 1 
  
  # # Event indication from pulse information
  # diffPulse <- diff(dftmp$X.PULSE, lag=1)
  # pulseStart <- which(diffPulse == 1)
  # eventPulse <- zeros(length(dftmp$X.PULSE))
  # for (i in 1:4){
  #   startIdx <- pulseStart[2*i-1]
  #   stopIdx <- pulseStart[2*i]
  #   eventPulse[startIdx:stopIdx] <- 1
  #   eventDetect[startIdx:stopIdx]<-removeSmallSegment(eventDetect[startIdx:stopIdx])
  # }
  # 
  # # Combing PULSE info with signal event
  # eventDetect <- eventDetect*eventPulse
  
  return(eventDetect)  
}
