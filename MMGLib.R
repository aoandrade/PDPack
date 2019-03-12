
# ---
# title:MMGLib.R
# author: "Alice Rueda"
# created: "February 19, 2019"
# last updated: "March 12, 2019"
#
# This file contains the functional calls for MMG analysis using the IMU accelerometer signals.
#
# Functions include:
# 1. Detect events based on signal power
# 2. Check if Pulse indication for start and stop of the event period
# 3. Removel small activity segments
# 4. Butterworth filter to remove motion artifacts
# 5. Segment signal dataframe into 4 task dataframes
# ---


# Event Detection ---------------------------------------------------------

library('signal')

# Remove smaller segment
removeSmallSegment<-function(eventDetectRegion){
  #Starting from event detected, ie eventDetectRegion[1]==1
  #eventDetectRegion <- eventDetect[3694:4413]
  tmp <-eventDetectRegion
  len <- length(tmp)
  
  # Finding number of segments
  tmp <- c(tmp, 0)                        # Padding zero for differentiation
  diffEventPulse <- diff(tmp, 1)
  stopSeg <- which(diffEventPulse == -1)  # Dropping edge of the segment
  startSeg <- which(diffEventPulse == 1)  # Raising edge of the segment
  numSegments <- length(stopSeg)
  
  # Finding the longest segment
  # Exit if there is only 1 segment
  if (numSegments == 1){
    return(eventDetectRegion)
    break
  }
  # Handing different starting point and calculate the first segment
  segLength <- zeros(numSegments)
  if (tmp[1]==1) {                        
    segLength[1]<-stopSeg[1]
    nextSeg <- 2
    stopSeg <- stopSeg[-1]    # Removing first falling edge 
  } else {
    nextSeg <- 1
  }
  
  # Finding the rest of the segment length
  for (i in nextSeg:numSegments){
    segLength[i]<-stopSeg[i] - startSeg[i]
  }
  
  # Finding the largest segment
  mainSegment<-which.max(segLength)
  
  # Removing all segments, except the longest one
  eventDetectRegion<-zeros(len)
  startEvent <- startSeg[mainSegment]
  stopEvent <- stopSeg[mainSegment]
  eventDetectRegion[startEvent:stopEvent]<- 1
  return(eventDetectRegion)
}



# Using signal information and pulse indications to detection events
detectMMGEvent <- function(dftmp, percentage = 0.08, PLOT=FALSE, ampFactor=1){
  
  # Adding optional amplification factor , weak signals
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
  
  # Event indication from pulse information
  diffPulse <- diff(dftmp$X.PULSE, lag=1)
  pulseStart <- which(diffPulse == 1)
  totalPulses <- length(pulseStart)
  
  
  if (totalPulses > 6) {         # This part takes care of 7 and 8 pulses
    eventPulse <- zeros(length(dftmp$X.PULSE))
    if (totalPulses == 8) {      # if there are 8 pulses, the pulses come in start-stop pair
      for (i in 1:4) {
        startIdx <- pulseStart[2*i-1]
        stopIdx <- pulseStart[2*i]
        eventPulse[startIdx:stopIdx] <- 1
        eventDetect[startIdx:stopIdx]<-removeSmallSegment(eventDetect[startIdx:stopIdx])
      }
    }
    else {                        # if there are 7 pulses, the first pulse is stop pulse, and then three start-stop pair
      startIdx <- 1
      stopIdx <- pulseStart[1]
      eventPulse[startIdx:stopIdx] <- 1
      eventDetect[startIdx:stopIdx]<-removeSmallSegment(eventDetect[startIdx:stopIdx])
      for (i in 1:3) {
        startIdx <- pulseStart[2*i]
        stopIdx <- pulseStart[2*i+1]
        eventPulse[startIdx:stopIdx] <- 1
        eventDetect[startIdx:stopIdx]<-removeSmallSegment(eventDetect[startIdx:stopIdx])
      }
    }
  }
  else {                          # Any number of pulses other than 7 and 8 will not be ignored
    eventPulse <- ones(length(dftmp$X.PULSE))
  }
   
  # Combing PULSE info with signal event
  eventDetect <- eventDetect*eventPulse
  
  if (PLOT) {
    mag = sqrt(dftmp$X.A1.X.^2+dftmp$X.A1.Y.^2+dftmp$X.A1.Z.^2)
    
    channel = dftmp$X.A1.X.
    dfplot<-data.frame(time = dftmp$X.Time.,channel, dftmp$X.PULSE*max(channel), mag, eventDetect*max(channel))
    dygraph(dfplot, main = "Accelerometer y-axis")%>%dyRangeSelector()%>% 
      dyOptions(colors = c('black', 'blue', 'green','red')) %>% 
      dyAxis("x", label="Time (s)") %>% 
      dyAxis("y", label="Angular velocity")
    
    
  }
  
  return(eventDetect)  
}




# Filters -----------------------------------------------------------------

# Butterworth Lowpass Filter
butterLowPassFilter <- function(channel, cutoff = 2, Fs = 50){
  Fn <- Fs/2
  f <- cutoff
  per <- f/Fn
  bf <- butter(3, per, type = "low")
  
  filtered <- filtfilt(bf, channel)
  
  return(filtered)
}

# Butterworth Highpass Filter
butterHighPassFilter <- function(channel, cutoff = 2, Fs = 50){
  Fn <- Fs/2
  f <- cutoff
  per <- f/Fn
  bf <- butter(3, per, type = "high")
  
  filtered <- filtfilt(bf, channel)
  
  return(filtered)
}



# Partitioning and Segmenting into Tasks Segments -------------------------

# Using the distance of 0 and 1 window length as the separation points
findSep <- function(eventDetect){
  changes <- diff(eventDetect, lag = 1)
  startIdx<-which(changes == 1)#[[1]]
  stopIdx<-which(changes == -1)
  numSegment <- length(startIdx)
  segmentLength <- stopIdx - startIdx
  zeroLength <- startIdx[2:numSegment]-stopIdx[1:numSegment-1]
  zeroLength <- c(zeroLength,(length(eventDetect)-stopIdx[numSegment]))
  sep = startIdx[1]
  total = startIdx[1]
  for (i in 1:numSegment){
    sep = c(sep,segmentLength[i],zeroLength[i])
    total = total +segmentLength[i]+zeroLength[i]
  }
  return(sep)
}

# Partition the signal according to the found separation points
partitionChannel <- function(channel,sep){
  #sep = findSep(eventDetect)
  partChannel <- partition.vector(channel,sep)
  return(partChannel)
}


# Convert Partitioned Data into Dataframe
convertMatrix2DataFrame<-function(mat, colHeader){
  dfmat <- as.data.frame(mat)
  names(dfmat) <- colHeader
  return(dfmat)
}

# Separate the task signals for all channels according to the detect event windowing

partitionDataFrame <- function(dfMMG, eventDetect) {
  
  # Parition according to the eventDetect windowing
  sep <- findSep(eventDetect)
  
  # Number of signals in a dataframe
  ncol <- dim(dfMMG)[2]
  
  # Declare the task matrix to hold segmented signals
  pinch <- matrix(0L, nrow = sep[2], ncol = ncol)
  hand <- matrix(0L, nrow = sep[4], ncol = ncol)
  pron <- matrix(0L, nrow = sep[6], ncol = ncol)
  flex <- matrix(0L, nrow = sep[8], ncol = ncol)
  
  # Partition all channels
  for (j in 1:ncol) {
    partSignal <- partitionChannel(dfMMG[,j],sep)
    
    pinch[,j] <- partSignal$'2'
    hand[,j]  <- partSignal$'4'
    pron[,j]  <- partSignal$'6'
    flex[,j] <- partSignal$'8'
    
  }  
  
  # Convert data matrix into dataframes
  colHeader <- colnames(dfMMG)
  dfpinch <- convertMatrix2DataFrame(pinch,colHeader)
  dfhand <- convertMatrix2DataFrame(hand,colHeader)
  dfpron <- convertMatrix2DataFrame(pron,colHeader)
  dfflex <- convertMatrix2DataFrame(flex,colHeader)
  
  dflist <- list(dfpinch, dfhand, dfpron, dfflex)
  return(dflist)
  
}
