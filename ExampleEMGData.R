# Example of extracting information from the EMG signal
# Load all necesssary libraries
library(dplyr)
source("read_Intan_RHD2000_file.R")
source("Thresholding.R")
source("EMGLib.R")



# Declaring parameters
# Pulse information is included in the excel file
pulseExist <- FALSE
# Get the cutoff frequency
cutoffFreq <- dfReference$cutoffFreq[1]


# All about the reference file --------------------------------------------

# Opening the reference file for the MVC values and construct a reference table for normalization
refFilename <- file.choose()  #choose the MVC reference to get the cutoffFreq
dfReference <- readWorkbook(refFilename,
                            sheet = 1,
                            detectDates = TRUE)

dfRef<-dfReference[c(4,5,7)]
refTable<-aggregate(dfRef, by=list(dfRef$channel, dfRef$movementType), FUN=max)

# Get the max value of a channel and the movement 
chan <- 4
movementType <- "flexao"
maxValue <- refTable[which((refTable$movementType==movementType) & (refTable$channel==chan)), 5]





# Working on the actual EMG signal ----------------------------------------


# Opening a selected EMG signal file
filename <- file.choose()
df <- readWorkbook(filename,
                   sheet = 1,
                   detectDates = TRUE)

Fs <- 1/(df$time[2]-df$time[1])



# Check if there is pulse indication
if (sum(df$Pulse) > 0) {
  pulseExist <- TRUE
}

# Event detection and segmentation
dfSig <- df[,2:4]

emgBurst <- NA
emgBurst <- detectEMGBurst(dfSig, df$pulse, Fs, cutoffFreq, percentage = 0.08, PLOT=TRUE)



# Select channel, detrend, filter, detect event, and segment signal
#chan<- list("chan.1", "chan.2", "chan.3")
sig <- df$chan.3
# Preprocessing
#sig<-filterHP(sig, Fs, cutFreq=200)
#sig.detrended <- detrendRaw(sig, df$time)
#s_env<-detectEnv(sig.detrended, Fs, cutoffFreq)
#sigFilt <- filterMVC(s_env, Fs, cutoffFreq)
#sigFilt <- filterMVC(sig.detrended, Fs, 10)

# Create labels and partition ----------------------------------------------------

Labels <- createEMGLabels(emgBurst,Fs, minTime=1)
df$label <- Labels

# Partition the signal
sep <- findSep(emgBurst)
partSig <- partitionChannel(sig,sep)


# Detrend and filter signal based on the type of activity
defaultCutoff <- c(5, 5, 2, 2)
detrendedPartSig <- list(partSig[[2]], partSig[[4]], partSig[[6]], partSig[[8]])
for (i in 1:4){
  tmpSig<-detrendedPartSig[[i]]
  sigTime<-seq(0, (length(tmpSig)-1)/Fs, by=1/Fs)
  
  # Detrending the signal
  detrendedPartSig[[i]]<-detrendEMGRaw(tmpSig, sigTime)
  
  # Apply low pass filter
  detrendedPartSig[[i]] <- filterMVC(detrendedPartSig[[i]], Fs, defaultCutoff[i])
  
  dfplot<-data.frame(time=sigTime, tmpSig, detrendedPartSig[[i]])
  g <- dygraph(dfplot) %>% dyRangeSelector() 
  print(g)
}

# Locating peaks and valleys
x <- detrendedPartSig[[1]]
sigTime<-seq(0, (length(x)-1)/Fs, by=1/Fs)
PeakValley <- findPeakValley(x)
xMaxIdx<- PeakValley[[1]]
xMaxValue <- PeakValley[[2]]
xMinIdx<-PeakValley[[3]]
xMinValue <- PeakValley[[4]]

maxLocation <- zeros(length(x))
maxLocation[xMaxIdx] <- max(x)
minLocation <- zeros(length(x))
minLocation[xMinIdx] <- min(x)

dfplot<-data.frame(time=sigTime, x,maxLocation, minLocation)
dygraph(dfplot)%>%dyRangeSelector()


# Plot Histogram ----------------------------------------------------------
deltaT=1/Fs

Intervals <- findInterval(xMaxIdx,xMinIdx, deltaT)

library(ggplot2)
library(grid)
library(gridExtra)

xInterPeakInterval = Intervals[[1]]
xInterValleyInterval = Intervals[[2]]
xCrossInterval = Intervals[[3]]

dfPeak <-data.frame(distance = xInterPeakInterval)
dfValley<-data.frame(distance = xInterValleyInterval)
dfCross<-data.frame(distance = xCrossInterval)
#dfPeak$dist <- rep(NA,length(xInterPeakInterval))
dfPeak$dist<-'Peak'
dfValley$dist<-'Valley'
dfCross$dist<-'Cross'

interTime<-rbind(dfPeak, dfValley, dfCross)
gg1<-ggplot(interTime, aes(distance, fill=dist))+geom_density(alpha=0.2)+
  xlab("Time (s)")+ylab("Count")+ggtitle("Inter Arriva lHistogram")
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


grid.arrange(gg1, gg2, gg3, gg4, gg5, gg6, nrow = 2, ncol=3, top = "")
























