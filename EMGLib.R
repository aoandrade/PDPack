# Nonlinear detrending (Step 1)
detrendEMGRaw <- function(rawSignal, t){
  dat <- data.frame(x = t, y = rawSignal)
  print("Start detrending. Please wait.")
  yp = predict(loess(y ~ x, dat, span = 0.1))
  ydetrended <- detrend(rawSignal-yp, tt = 'linear')
  return(as.vector(ydetrended))
}


# Using signal information and pulse indications to detection events
detectEMGBurst <- function(dfSig, pulse, Fs, cutoffFreq=1, percentage = 1, PLOT=FALSE){
  # Finding threshold, setting at 5% of noise removed maximum magnitude
  sigTime <- seq(0,(length(dfSig[,1])-1)/Fs,by=1/Fs)
  
  filtered <- zeros(length(dfSig[,1]))
  for (i in 1:length(dfSig)){
    #sig.detrended<- detrendRaw(sig, sigTime)
    sig.detrended <- filterHP(dfSig[,i], Fs, cutFreq=200)
    s_env <- detectEnv(sig.detrended, Fs, cutoffFreq)
    filtered <- filtered + filterMVC(s_env, Fs, cutoffFreq)
  }
  
  # Finding the threshold using histogram
  h <- hist(filtered)
  threshold <- h$mids[which.max(h$counts)]
  threshold <- threshold*percentage
  
  
  # Using low pass to remove high frequency noise
  #filtered <- filterMVC(tmp, Fs, cutFreq = 2)
  
  # Event detection from signal
  eventDetect <- filtered
  eventDetect[filtered < threshold] <- 0 
  eventDetect[filtered >= threshold] <- 1 
  
  # Event indication from pulse information
  diffPulse <- diff(pulse, lag=1)
  pulseStart <- which(diffPulse == 1)
  eventPulse <- zeros(length(pulse))
  for (i in 1:4){
    startIdx <- pulseStart[2*i-1]
    stopIdx <- pulseStart[2*i]
    eventPulse[startIdx:stopIdx] <- 1
    eventDetect[startIdx:stopIdx]<-removeSmallSegment(eventDetect[startIdx:stopIdx])
  }
  
  # Combing PULSE info with signal event
  eventDetect <- eventDetect*eventPulse
  
  if (PLOT) {
    dfplot<-data.frame(df$time, sig.detrended, pulse*max(sig.detrended), filtered/max(filtered)*max(sig.detrended),eventDetect*max(sig.detrended))
    g<-dygraph(dfplot) %>% dyRangeSelector() %>%
      dyOptions(colors = c('black', 'blue', 'green','red')) %>% 
      dyAxis("x", label="Time (s)") %>% 
      dyAxis("y", label="Magnitude")
    print(g)
  }
  
  return(eventDetect)  
}


# Function to create label vector
createEMGLabels<-function(eventDetect, Fs, minTime=1){
  minTime <- 1 #If the duration of the segment is less than 1 second than it is the same event
  sep <- findSep(eventDetect)
  sepLen <- length(sep)
  Labels <- rep(0, length(eventDetect))
  label <- 0   #The first label is no event
  segCount <-2 #The first segment is no event, so start at 2nd segment
  deltaT <- 1/Fs
  endIdx <- sep[segCount-1]
  while (segCount<sepLen){
    segTime <- sep[segCount]*deltaT #The duration of the segment period 
    startIdx <- endIdx + 1
    endIdx <- startIdx + sep[segCount]-1
    if((segCount %% 2) == 0){  #if even segment
      if (segTime>minTime){
        label <- label + 1
        Labels[startIdx:endIdx]<-label
      }
    } else if (segTime<minTime){  #if odd segment and less than minTine
      #label <- label
      Labels[startIdx:endIdx]<-label   # same label as last
      segCount <- segCount + 1        # also label the next one the same
      startIdx <- endIdx+1
      endIdx <- startIdx + sep[segCount]-1
      Labels[startIdx:endIdx]<-label   # same label as last
    } 
    segCount <- segCount + 1
  }  
  return(Labels)
}