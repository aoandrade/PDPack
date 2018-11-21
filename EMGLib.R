source("Thresholding.R")

# Functions for MVC --------------------------------------------------
# Normalize the Results using the maximum value of the channel in maxRefTable
getEMGNormResults<-function(Results, maxRefTable){
  normResults<-Results
  dimResults<-size(normResults)
  nChannel <- dimResults[1]
  nBurst <- dimResults[2]
  for (i in 1:nChannel){
    maxValue <- maxRefTable$medianMVC[i]
    for (j in 1:nBurst){
      normResults[[i,j]]$peakValue <- Results[[i,j]]$peakValue/maxValue
      normResults[[i,j]]$valleyValue <- Results[[i,j]]$valleyValue/maxValue
    }
  }
  return(normResults)  
}

# Getting the maximum reference MVC for each channel
getEMGMaxRefTable<-function(filename){
  dfReference <- readWorkbook(filename,
                              sheet = 1,
                              detectDates = TRUE)
  
  dfRef<-dfReference[c(5, 7, 4)]
  refTable<-aggregate(dfRef, by=list(dfRef$movementType, dfRef$channel), FUN=max)
  refTable<-refTable[,3:5]
  
  attach(refTable)
  maxRefTable <- aggregate(medianMVC ~channel, refTable, FUN=max)
  detach(refTable)
  
  return(maxRefTable)
  
}


# Function to get the MVC for each channel and movement type from a CVM*.xlsx file
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
    dfNewChannel$channel <- chanNum-1  
    
    # Pre-processing the channel: detrending, getting signal envelope, and filtering
    sig<-df[,chanNum]
    
    sig<-filterHP(sig, Fs, cutFreq=200)
    sig.detrended <- detrendRaw(sig, df$time)
    s_env<-detectEnv(sig.detrended, Fs, cutoffFreq)
    sigFilt <- filterMVC(s_env, Fs, cutoffFreq)
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


# Nonlinear detrending (Step 1)
detrendRaw <- function(rawSignal, t){
  dat <- data.frame(x = t, y = rawSignal)
  print("Start detrending. Please wait.")
  yp = predict(loess(y ~ x, dat, span = 0.1))
  ydetrended <- detrend(rawSignal-yp, tt = 'linear')
  return(ydetrended)
}


detectMVCEvent <- function(s_env, percentage = 0.5, Fs){
  # Finding threshold, setting at 5% of noise removed maximum magnitude
  
  threshold = (max(s_env)-min(s_env))*percentage+min(s_env)
  minSegmentLength <- Fs*percentage
  
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


# Nonlinear detrending
detrendEMGRaw <- function(rawSignal, t){
  dat <- data.frame(x = t, y = rawSignal)
  print("Start detrending. Please wait.")
  yp = predict(loess(y ~ x, dat, span = 0.1))
  ydetrended <- detrend(rawSignal-yp, tt = 'linear')
  return(as.vector(ydetrended))
}


# Using signal information and pulse indications to detection events
detectEMGBurst <- function(dfSig, cutoffFreq=1, percentage = 1, PLOT=FALSE){
  
  nChannels <- 3
  
  # Finding threshold, setting at 5% of noise removed maximum magnitude
  Fs <- 1/(dfSig$time[2]-dfSig$time[1])
  sigTime <- seq(0,(length(dfSig[,1])-1)/Fs,by=1/Fs)
  
  filtered <- zeros(length(dfSig[,1]))
  
  for (i in 2:nChannels+1){
    #sig.detrended<- detrendRaw(sig, sigTime)
    sig.detrended <- filterHP(dfSig[,i], Fs, cutFreq=200)
    s_env <- detectEnv(sig.detrended, Fs, cutoffFreq)
    filtered <- filtered + filterMVC(s_env, Fs, cutoffFreq)
  }
  
  # Finding the threshold using histogram
  h <- hist(filtered, plot=FALSE)
  threshold <- h$mids[which.max(h$counts)]
  threshold <- threshold*percentage
  
  
  # Using low pass to remove high frequency noise
  #filtered <- filterMVC(tmp, Fs, cutFreq = 2)
  
  # Event detection from signal
  eventDetect <- filtered
  eventDetect[filtered < threshold] <- 0 
  eventDetect[filtered >= threshold] <- 1 
  
  # Check if there is pulse indication
  if (sum(dfSig$pulse) == 0) {
    eventPulse <- ones(length(dfSig$time))
  } else {
    # Event indication from pulse information
    pulse <- dfSig$pulse
    diffPulse <- diff(pulse, lag=1)
    pulseStart <- which(diffPulse == 1)
    eventPulse <- zeros(length(pulse))
    for (i in 1:4){
      startIdx <- pulseStart[2*i-1]
      stopIdx <- pulseStart[2*i]
      eventPulse[startIdx:stopIdx] <- 1
      eventDetect[startIdx:stopIdx]<-removeSmallSegment(eventDetect[startIdx:stopIdx])
    }
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


# Take in EMG signal file, extract max, min, and time difference information for all channels and bursts
getEMGTimeFeatures<-function(df, PLOT=TRUE){
  # Process all channel
  nChannel <- 3
  nBurst <- 4
  deltaT <- df$time[2]-df$time[1]
  Fs <- 1/deltaT
  defaultCutoff <- c(5, 5, 2, 2)
  results <- array(list(), c(nChannel,nBurst))
  
  
  # Return peakTime, valleyTime, peakValue, valleyValue, interPeakInterval, interValleyInterval, crossInterval, Label 
  emgBurst <- detectEMGBurst(df, cutoffFreq=1, percentage = 0.08, PLOT)
  
  # Create labels and partition ----------------------------------------------------
  #Labels <- createEMGLabels(emgBurst,Fs, minTime=1)
  #df$label <- Labels
  
  # Per channel partition: separate the signal into activities
  sep <- findSep(emgBurst)
  for (channel in 2:(nChannel+1)) {
    fprintf("Processing channel %d. \n", channel-1)
    sig <- df[,channel]
    partSig <- partitionChannel(sig,sep)
    
    # Detrend and filter signal based on the type of activity
    detrendedPartSig <- list(partSig[[2]], partSig[[4]], partSig[[6]], partSig[[8]])
    for (i in 1:4){
      fprintf("Burst: %d ", i)
      tmpSig<-detrendedPartSig[[i]]
      sigTime<-seq(0, (length(tmpSig)-1)/Fs, by=1/Fs)
      
      # Detrending the signal
      detrendedPartSig[[i]]<-detrendEMGRaw(tmpSig, sigTime)
      
      # Apply low pass filter
      detrendedPartSig[[i]] <- filterMVC(detrendedPartSig[[i]], Fs, defaultCutoff[i])
      
      if (PLOT){
        Title <- paste("Channel",channel, "burst", i)
        dfplot<-data.frame(time=sigTime, tmpSig, detrendedPartSig[[i]])
        g <- dygraph(dfplot, main = Title) %>% dyRangeSelector() %>%
          dyAxis("x", label = "Time (s)") %>%
          dyAxis("y", label = "Amplitude (uV)")
        print(g)
      }
      
      
      # Locating peaks and valleys ----------------------------------------------
      sigTime<-seq(0, (length(tmpSig)-1)/Fs, by=1/Fs)
      PeakValley <- findPeakValley(detrendedPartSig[[i]])
      
      
      xMaxIdx <- PeakValley[[1]]
      maxTime <- xMaxIdx*deltaT
      xMaxValue <- PeakValley[[2]]
      
      xMinIdx<-PeakValley[[3]]
      minTime <- xMinIdx*deltaT
      xMinValue <- PeakValley[[4]]
      
      maxLocation <- zeros(length(tmpSig))
      maxLocation[xMaxIdx] <- max(tmpSig)
      minLocation <- zeros(length(tmpSig))
      minLocation[xMinIdx] <- min(tmpSig)
      
      if (PLOT){
        dfplot<-data.frame(time=sigTime, detrendedPartSig[[i]], maxLocation, minLocation)
        g<-dygraph(dfplot, main = Title)%>%dyRangeSelector() %>%
          dyAxis("x", label = "Time (s)") %>%
          dyAxis("y", label = "Amplitude (uV)")
        print(g)
          
      }
      
      
      # Finding intervals ----------------------------------------------------------
      
      Intervals <- findInterval(xMaxIdx,xMinIdx, deltaT)
      
      xInterPeakInterval = Intervals[[1]]
      xInterValleyInterval = Intervals[[2]]
      xCrossInterval = Intervals[[3]]
      
      dfPeak <- data.frame(distance = xInterPeakInterval)
      dfValley <- data.frame(distance = xInterValleyInterval)
      dfCross <- data.frame(distance = xCrossInterval)
      #dfPeak$dist <- rep(NA,length(xInterPeakInterval))
      dfPeak$dist<-'Peak'
      dfValley$dist<-'Valley'
      dfCross$dist<-'Cross'
      
      
      if (PLOT){
        interTime<-rbind(dfPeak, dfValley, dfCross)
        gg1<-ggplot(interTime, aes(distance, fill=dist))+geom_density(alpha=0.2)+
          xlab("Time (s)")+ylab("Count")+ggtitle("Histogram")
        gg2<-qplot(x=(xMaxIdx-1)*deltaT, y=xMaxValue)+
          xlab("Time (s)")+ylab("Magnitude")+ggtitle("Peak Values")
        gg3<-qplot(x=(xMinIdx-1)*deltaT, y=xMinValue)+
          xlab("Time (s)")+ylab("Magnitude")+ggtitle("Valley Values")
        gg4<-qplot(x=1:length(xInterPeakInterval), y=xInterPeakInterval)+
          xlab("Count")+ylab("Time (s)")+ggtitle("Inter Peak Interval")
        gg5<-qplot(x=1:length(xInterValleyInterval), y=xInterValleyInterval)+
          xlab("Count")+ylab("Time (s)")+ggtitle("Inter Valley Interval")
        gg6<-qplot(x=1:length(xCrossInterval), y=xCrossInterval)+
          xlab("Count")+ylab("Time (s)")+ggtitle("Cross Interval")
        
        
        grid.arrange(gg1, gg2, gg3, gg4, gg5, gg6, nrow = 2, ncol=3, top = Title)
      }
      
      
      results[[channel-1,i]]<-list(peakTime=maxTime, 
                                   peakValue=xMaxValue, 
                                   valleyTime=minTime,  
                                   valleyValue=xMinValue, 
                                   interPeakInterval=xInterPeakInterval, 
                                   interValleyInterval=xInterValleyInterval, 
                                   crossInterval=xCrossInterval,
                                   burstLabel=i)
    }
    
  }
  
  return(results)
  
}
