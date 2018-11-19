#clear() # Clear the environment
library(signal)
library(dygraphs)
source("TremsenToolbox.r")


# Function Calls for MVC --------------------------------------------------
# Function to get the MVC from a file
getMVC<-function(filename, pathData, subjectCode, chanList = list(2,3), cutoffFreq=1, percentage = 0.5, PLOT=FALSE){
  filename<-paste(pathData,filename, sep="")
  df <- readWorkbook(filename,
                     sheet = 1,
                     detectDates = TRUE)
  
  Fs <- 1/(df$time[2]-df$time[1]) # sampling frequency in Hz
  
  # define the dataframe for recording
  dfResult <- data.frame(
    subjectCode=as.character(NA),
    filename=as.character(NA), 
    dateYMD=as.character(NA), 
    movementType=as.character(NA),
    channel=as.integer(NA), 
    burstIdx=as.integer(NA),
    medianMVC=as.numeric(NA),
    cutoffFreq=as.integer(NA),
    stringsAsFactors = FALSE
  )
  
  dfNewFile<-dfResult
  
  fprintf("Processing %s \n", filename)
  
  # calculate and record results
  for (i in 1:length(chanList)){
    fprintf("Processing channel number %d \n", chanList[[i]])
    dfNewChannel <- dfNewFile
    chanNum <- chanList[[i]]    # Channel number as provided on chanList with offset 1 for df$time
    dfNewChannel$channel <- chanNum  
    
    # Pre-processing the channel: detrending, getting signal envelope, and filtering
    sig<-df[,chanNum]
    sig<-filterHP(sig, Fs, cutFreq=200)
    sig.detrended <- detrendRaw(sig, df$time)
    s_env<-detectEnv(sig.detrended, Fs, cutFreq = 1)
    sigFilt <- filterMVC(s_env, Fs, cutFreq = 2)
    #sigFilt <- detrend(sigFilt, tt="linear")
    
    # Detect event based on the detrended raw signal
    eventDetect <- NA
    eventDetect = detectMVCEvent(sigFilt, percentage = 0.4, Fs)
    
    #tsig<-(1:length(signal))/Fs
    if (PLOT) {
      dfplot <- data.frame(time=df$time, sig, s_env, sigFilt, eventDetect*max(sig))
      p<-dygraph(dfplot, main = paste("Channel ", chanNum, sep = "")) %>%
        dyOptions(colors = c('black', 'blue', 'green','red')) %>% 
        dyAxis("x", label="Time (s)") %>% 
        dyAxis("y", label="Magnitude")
      print(p)
    }
    
    
    # Partition the signal
    sep <- findSep(eventDetect)
    partSig <- partitionChannel(sigFilt,sep)
    partEvent <- partitionChannel(eventDetect,sep)
    eventIndicator <- zeros(length(partEvent))
    fprintf("Number of segments is %d \n", length(partEvent))
    for (k in 1: length(partEvent)){
      eventIndicator[k] <- round(mean(partEvent[[k]]))
    }
    
    
    # Finding the median value
    medianSig <- 0
    count <- 0
    threshold <- max(sigFilt)*percentage
    for (j in 1:length(eventIndicator)) {
      tmpMedian <- median(abs(partSig[[j]]))
      
      if (eventIndicator[j] == 1) {
        dfNewBurst<-dfNewChannel
        dfNewBurst$medianMVC <- tmpMedian
        medianSig <- medianSig+tmpMedian
        print(medianSig)
        count <- count + 1
        dfNewBurst$burstIdx <- count
        dfResult <- rbind(dfResult, dfNewBurst)
        fprintf("The segment number %d is active \n", j)
      }
    }
    medianSig <- medianSig/count
  }
  #dfSave <- dfResult[2:length(dfResult),]
  
  # Assigning constant values to all entries: filename, subjectCode,movementType, and cutoffFreq
  dfResult$filename<-filename
  str1 <- "CVM_"
  str2 <- "_"
  Type <- str_match(filename, paste(str2, "(.*?)", str2, sep = ''))
  Type <- Type[,2]
  dfResult$movementType <- Type
  Date <- str_match(filename, paste(str1, Type, str2, "(.*?)", str2, sep = ''))
  Date <- Date[,2]
  Date <- paste("20", substr(Date,1,2),"-",substr(Date,3,4), "-",substr(Date,5,6), sep = '')
  dfResult$dateYMD <- Date
  dfResult$cutoffFreq<-cutoffFreq
  dfResult$subjectCode <- subjectCode
  
  # Removing first entry and renumber row index
  dfResult<-dfResult[-1,]
  rownames(dfResult)<-seq(nrow(dfResult))
  
  return(dfResult)
}

# High pass filter
filterHP<-function(sig, Fs, cutFreq=200){
  Fn <- Fs/2
  per <- cutFreq/Fn
  bf <- butter(3, per, type='high')
  s_filt <- filtfilt(bf, sig)   
  return(s_filt)
}

# Detect Evelope
detectEnv<-function(sig, Fs, cutFreq = 1){
  yenv <- abs(seewave::hilbert(sig, f = Fs, fftw = TRUE))
  Fn <- Fs/2
  per <- cutFreq/Fn
  bf <- butter(3, per, type='low')
  s_env <- filtfilt(bf, yenv)   #s_env is the envelope
  return(s_env)
}


# Filter model for MVC
filterMVC<-function(sig, Fs, cutFreq = 2){
  Fn <- Fs/2
  per <- cutFreq/Fn
  bf <- butter(3, per, type='low')
  s_filt <- filtfilt(bf, sig)
  return(s_filt)
}

# Nonlinear detrending 
detrendRaw <- function(rawSignal, t){
  dat <- data.frame(x = t, y = rawSignal)
  print("Start detrending. Please wait.")
  yp = predict(loess(y ~ x, dat, span = 0.1))
  ydetrended <- detrend(rawSignal-yp, tt = 'linear')
  return(ydetrended)
}

# Using signal information and pulse indications to detection MVC events
detectMVCEvent <- function(s_env, percentage = 0.4, Fs){
  # Finding threshold, setting at 5% of noise removed maximum magnitude
  # Forming envelop --------------------------------------------------------
  #yenv <- abs(seewave::hilbert(dataVector, f = Fs, fftw = TRUE))
  #Fn <- Fs/2
  #f <- 1
  #per <- f/Fn
  #bf <- butter(3, per, type='low')
  #s_env <- filtfilt(bf, yenv)   #s_env is the envelope
  
  threshold = (max(s_env)-min(s_env))*percentage+min(s_env)
  minSegmentLength <- Fs*0.5
  
  # Event detection from signal
  eventDetect <- s_env
  eventDetect[s_env < threshold] <- 0 
  eventDetect[s_env >= threshold] <- 1 
  
  eventDetectDiff <- diff(eventDetect,1)
  eventDetectRaiseIdx <- which(eventDetectDiff == 1)
  eventDetectFallIdx <- which(eventDetectDiff == -1)
  for (i in 1:length(eventDetectFallIdx)){
    lengthEvent <- eventDetectFallIdx[i] - eventDetectRaiseIdx[i]
    if (lengthEvent<minSegmentLength){
      eventDetect[eventDetectRaiseIdx[i]:eventDetectFallIdx[i]] <- 0
    }
  }
  
  return(eventDetect)  
}


# Function Calls for Trail Signalss ---------------------------------------

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
detectEvent <- function(dftmp, percentage = 0.08){
  # Finding threshold, setting at 5% of noise removed maximum magnitude
  tmp <- sqrt(dftmp$X.G1.Z.^2+dftmp$X.G1.Y.^2+dftmp$X.G1.X.^2)
  # Using low pass to remove high frequency noise
  Fs <- 1/0.02
  Fn <- Fs/2
  f <- 0.5
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
  eventPulse <- zeros(length(dftmp$X.PULSE))
  for (i in 1:4){
    startIdx <- pulseStart[2*i-1]
    stopIdx <- pulseStart[2*i]
    eventPulse[startIdx:stopIdx] <- 1
    eventDetect[startIdx:stopIdx]<-removeSmallSegment(eventDetect[startIdx:stopIdx])
  }
  
  # Combing PULSE info with signal event
  eventDetect <- eventDetect*eventPulse

  return(eventDetect)  
}

# Segment events, assuming no events at the begining and return the event segments in even list
#MODIFICATION: use PULSE information for demarkation
#MODIFICATION: get rid of the small separations
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

# Partition the signal with return separation indice, provided that sep has been found using findSep()
partitionChannel <- function(channel,sep){
  #sep = findSep(eventDetect)
  partChannel <- partition.vector(channel,sep)
  return(partChannel)
}

# Function to create label vector
createLabels<-function(eventDetect, minTime=2, deltaT=0.02){
  minTime <- 2 #If the duration of the segment is less than 2 seconds than it is the same event
  sep <- findSep(eventDetect)
  sepLen <- length(sep)
  Labels <- rep(0, length(eventDetect))
  label <- 0   #The first label is no event
  segCount <-2 #The first segment is no event, so start at 2nd segment
  deltaT <- df.nonlineardetrended$X.Time.[2]  #To determine the sampling time
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


# This is easier to use than R's built-in acf()
autoCorr <- function(x){
  len = length(x)
  acfx<-rep(0, len)
  for (lags in 0:(len-1)){
    xlag=x[(lags+1):len]
    xtmp=x[1:(len-lags)]
    acfx[lags+1] = sum(t(xtmp)*xlag)
    lastIdx <- len-lags
    xlag[lastIdx] <-0
    xtmp[lastIdx] <-0
  }
  acfx = acfx/acfx[1]  #Normalize the acf
  return(acfx)
}




# Finding peak and valley locations and their magnitudes
findPeakValley<-function(x, deltaT=0.02){
  
  t <-(1:length(x))*deltaT

  xDiff <- diff(x, lag =1)  # ist derivatives
  xDiff <-c(0, xDiff)
  xDiffSign <- xDiff        # if the derivative is +ve, signal is increasing and will be labeled 1. Otherwise, -1
  xDiffSign[xDiff >= 0] <- 1  
  xDiffSign[xDiff < 0] <- -1
  xDDSign <- - c(diff(xDiffSign, lag=1),0)  # "-ve 2nd derivative"
  xMinIdx <- which(xDDSign == -2)           # Valley (min) if "-ve 2nd derivative" is negative, else peak (max)
  xMaxIdx <- which(xDDSign == 2)
  xMinValue <- x[xMinIdx]
  xMaxValue <- x[xMaxIdx]
  PeakValley <- list(xMaxIdx, xMaxValue, xMinIdx, xMinValue)
  return(PeakValley)
}



# Finding interval between peaks, valleys, and peak-valley
findInterval<-function(xMaxIdx,xMinIdx, deltaT=0.02){
  xInterPeakInterval <- diff(xMaxIdx, lag=1)*deltaT#*1000
  xInterValleyInterval <- diff(xMinIdx, lag=1)*deltaT#*1000
  if (length(xMaxIdx)>length(xMinIdx)){
    xCrossInterval <- (xMinIdx - xMaxIdx[1:length(xMaxIdx)-1])*deltaT#*1000
  } else {
    xCrossInterval <- (xMinIdx - xMaxIdx)*deltaT#*1000
  }
  Intervals <- list(xInterPeakInterval,xInterValleyInterval, xCrossInterval)
  return(Intervals)
}




