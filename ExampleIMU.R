# Description -------------------------------------------------------------
# ExampleIMU.R
#
#
# 1. To run this script you should have available a .txt generated by the hardware TREMSEN
# 2. This script illustrates how to estimate relavant information from inertial signals

# Libraries ---------------------------------------------------------------
source("Thresholding.r")
source("TremsenToolbox.r") # Latest version available @ https://github.com/NIATS-UFU/TREMSEN-Toolbox.git


# Load TREMSEN file -------------------------------------------------------

# Choose a TREMSEN (.txt) file 
# strFileName <- file.choose() # If data are available, select the file

strFileName <- "./Data/AFP_E_1.txt" # Sample data

df <- LoadTREMSENFile(strFileName) # Load TREMSEN file

## START working on the signal
#Fs<-1000 #upsampling to 1000Hz
#df<-resampleTremsenData(df,Fs)
dfplot<-data.frame(df$X.Time., df$X.A1.X.)
dygraph(dfplot)%>%dyRangeSelector()



# Detrend and remove DC offset --------------------------------------------
df.nonlineardetrended <- nonLineardetrendTremsenData(df) #Remove both linear and nonlinear trends

# Event detection  ---------------------------------------------------

# Using detrend dataframe for event detection
dftmp<-df.nonlineardetrended
eventDetect = detectEvent(dftmp)

mag = sqrt(dftmp$X.G1.X.^2+dftmp$X.G1.Y.^2+dftmp$X.G1.Z.^2)

channel = df.nonlineardetrended$X.G1.Y.
dfplot<-data.frame(time = df$X.Time.,df$X.PULSE*max(channel), channel, mag, eventDetect*max(channel))
dygraph(dfplot)%>%dyRangeSelector()



# Create labels and add to PULSE.LABEL ----------------------------------------------------

Labels <- createLabels(eventDetect)
df.nonlineardetrended$X.PULSE.LABEL <- Labels



# Filter signal -----------------------------------------------------------


# Define butterworth filter to remove high frequency noise
Fs <- 1/0.02
Fn <- Fs/2
f <- 2
per <-f/Fn
bf <- butter(3, per, type="low")
#freqz(bf)
#zplane(bf)

# Filter out high frequency noise
filtered <- filtfilt(bf, channel)
#dfplot<-data.frame(time = df$X.Time., eventDetect*max(channel), filtered, channel, (channel-filtered))
dfplot<-data.frame(time = df$X.Time., df.nonlineardetrended$X.PULSE*max(channel), channel, filtered, df.nonlineardetrended$X.PULSE.LABEL/4*max(channel))
dygraph(dfplot)%>%dyRangeSelector()



# Partition signal -----------------------------------------------
sep <- findSep(eventDetect)
partChannel <- partitionChannel(filtered,sep)
partOrignal <- partitionChannel(channel,sep)
numSegment <- size(partChannel,2)  # The activities are stored in even number list of partChannel


# Finding peaks and valleys -----------------------------------------------
#Select the even number of partition where events have been detected
x <- partChannel$`8`
xOriginal <-partOrignal$`8`
deltaT = 0.02 #Time between data points

tsig<-c(0:(length(x)-1))*deltaT

PeakValley <- findPeakValley(x)
xMaxIdx<- PeakValley[[1]]
xMaxValue <- PeakValley[[2]]
xMinIdx<-PeakValley[[3]]
xMinValue <- PeakValley[[4]]

maxLocation <- zeros(length(x))
maxLocation[xMaxIdx] <- max(x)
minLocation <- zeros(length(x))
minLocation[xMinIdx] <- min(x)

dfplot<-data.frame(time=tsig, x,maxLocation, minLocation)
dygraph(dfplot)%>%dyRangeSelector()


# Plot Histogram ----------------------------------------------------------


Intervals <- findInterval(xMaxIdx,xMinIdx)

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
  xlab("Time (s)")+ylab("Density")+ggtitle("Kernel density estimation")
gg2<-qplot(x=(1:length(xMaxValue))*deltaT, y=xMaxValue)+
  xlab("Time (s)")+ylab("Magnitude")+ggtitle("Peak Values")
gg3<-qplot(x=(1:length(xMinIdx))*deltaT, y=xMinValue)+
  xlab("Time (s)")+ylab("Magnitude")+ggtitle("Valley Values")
gg4<-qplot(x=1:length(xInterPeakInterval), y=xInterPeakInterval)+
  xlab("Count")+ylab("Time (s)")+ggtitle("Inter Peak Interval")
gg5<-qplot(x=1:length(xInterValleyInterval), y=xInterValleyInterval)+
  xlab("Count")+ylab("Time (s)")+ggtitle("Inter Valley Interval")
gg6<-qplot(x=1:length(xCrossInterval), y=xCrossInterval)+
  xlab("Count")+ylab("Time (s)")+ggtitle("Cross Interval")


grid.arrange(gg1, gg2, gg3, gg4, gg5, gg6, nrow = 2, ncol=3, top = "")

