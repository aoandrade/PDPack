# ---
# title:MMGLib.R
# author: "Alice Rueda"
# created: "November 7, 2018"
# last updated: "September 19, 2019"
#
# This file was originally created to analysis bradykinesia analysis for PDPack to demonstrate how to use
# Thresholding.R, IMULib.R and other library files to analyze and plot the results.
#



setwd('/media/alice/DATA/2019 Adriano/PDPack')
source("Thresholding.R")
source("TREMSENToolbox.R")
source("PlotSignals.R")
library('dygraphs')

# Main function -----------------------------------------------------------

# Loading function and data

strFolder <-  "./Data/"
strFileName <- "AFP_E_1.txt"
#strFileName <- "CAG_D_1.txt"
#strFileName <- "MAS_D_1.txt"
df <- LoadTREMSENFile(paste(strFolder,strFileName,sep = ""))

## START working on the signal
#Fs<-1000 #upsampling to 1000Hz
#df<-resampleTremsenData(df,Fs)
dfplot<-data.frame(df$X.Time., df$X.A1.Z.)
dygraph(dfplot)%>%dyRangeSelector()

gg1<-qplot(x=df$X.Time., df$X.A1.X.)+
  xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A1.X")
gg2<-qplot(x=df$X.Time., df$X.A1.Y.)+
  xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A1.Y")
gg3<-qplot(x=df$X.Time., df$X.A1.Z.)+
  xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A1.Z")
gg4<-qplot(x=df$X.Time., df$X.G1.X.)+
  xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.G1.X")
gg5<-qplot(x=df$X.Time., df$X.G1.Y.)+
  xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.G1.Y")
gg6<-qplot(x=df$X.Time., df$X.G1.Z.)+
  xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.G1.Z")
# gg7<-qplot(x=df$X.Time., df$X.A3.X.)+
#   xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A3.X")
# gg8<-qplot(x=df$X.Time., df$X.A3.Y.)+
#   xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A3.Y")
# gg9<-qplot(x=df$X.Time., df$X.A3.Z.)+
#   xlab("Time (s)")+ylab("Magnitude")+ggtitle("X.A3.Z")


grid.arrange(gg1, gg2, gg3, gg4, gg5, gg6, nrow = 2, ncol=3, top = "") #, gg7, gg8, gg9, nrow = 3, ncol=3, top = "")



# Detrend and remove DC offset --------------------------------------------

df.nonlineardetrended <- nonLineardetrendTremsenData(df) #Remove both linear and nonlinear trends

# Event detection  ---------------------------------------------------

# Using detrend dataframe for event detection
dftmp<-df.nonlineardetrended
eventDetect = detectEvent(dftmp)

mag = sqrt(dftmp$X.A1.X.^2+dftmp$X.A1.Y.^2+dftmp$X.A1.Z.^2)

channel = df.nonlineardetrended$X.A1.Z.
dfplot<-data.frame(time = df$X.Time.,channel, df$X.PULSE*max(channel), mag, eventDetect*max(channel))
dygraph(dfplot, main = "Accelerometer 1 z-axis")%>%dyRangeSelector()%>% 
  dyOptions(colors = c('black', 'blue', 'green','red')) %>% 
  dyAxis("x", label="Time (s)") %>% 
  dyAxis("y", label="Angular velocity")



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


#Plot LPF signal
dfplot<-data.frame(df$X.Time., filtered)
dygraph(dfplot)%>%dyRangeSelector()

# Plot Event Detect with
Pulses <- df$X.PULSE*max(filtered)
Events <- eventDetect*max(filtered)
dfplot<-data.frame(time = df$X.Time.,filtered, Pulses, Events)
dygraph(dfplot, main = "Accelerometer 1 z-axis")%>%dyRangeSelector()%>% 
  dyOptions(colors = c('black', 'blue', 'red')) %>% 
  dyAxis("x", label="Time (s)") %>% 
  dyAxis("y", label="Acceleration (m/s^2)")

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

