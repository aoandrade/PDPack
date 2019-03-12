# ---
# title: "MMGExample"
# author: "Alice Rueda"
# date: "Mar 12, 2019"
# output: html_document
# 
# This is an example of using MMGlib.R functions to:
# 1. Open text data file to read in the dataframe
# 2. Detrend and remove the raw accelerometer signals
# 3. Detect event and construct a window frame to demarcate the start and end of each task
# 4. Remove low frequency motion artifact using a HP Butterworth filter to obtain the MMG signal
# 5. Partition the MMG Signals according to the 4 task windows
# 6. Save the MMG signals into 4 excel files
# ---

<<<<<<< HEAD


source("MMGLib.R")
source("PlotSignals.R")
source("TREMSENToolbox.R")
#install.packages("xlsx")   # install xlsx package if it is not on the package list
library("xlsx")


#subjectCode: get it from the directory
pathData = "/media/alice/DATA/2019 Adriano/2019 Adriano MMG/Data/Coletas 14_11/MMG_Tremsen/COLETAS 14_11/"
signalType = "MMG"
fileType = "txt"




# Read Data and Plot Raw Signal -------------------------------------------

# Get the list of text files in the directory
filename <- sort(list.files(path=pathData, full.names=TRUE, pattern="*.txt"))

# Select which file to work on and plot raw signal
i <- 1
df <- LoadTREMSENFile(filename[i])
ggplotAccelerometer(df)




# Detrend and Remove DC offset --------------------------------------------

df.nonlineardetrended <- nonLineardetrendTremsenData(df) #Remove both linear and nonlinear trends
ggplotAccelerometer(df.nonlineardetrended)





# Event Detection ---------------------------------------------------------

# Using detrend dataframe for event detection
dfdetrended<-df.nonlineardetrended

eventDetect = detectMMGEvent(dfdetrended, percentage = 0.05, PLOT=TRUE)

scaling = 0.5
channel = dfdetrended$X.A1.X.
dfplot1<-data.frame(time = dfdetrended$X.Time., channel, dfdetrended$X.PULSE*max(channel)*scaling, eventDetect*max(channel)*scaling)
dygraph(dfplot1, main = "Accelerometer 1 x-axis") %>%
  dyAxis("x", label="Time (s)") %>% 
  dyAxis("y", label="Angular velocity") %>%
  dyOptions(colors = c('black', 'blue', 'red')) %>% 
  dyLegend(width = 130) %>%
  dyRangeSelector() 





# Remove Motion Artifact --------------------------------------------------

dfMMG <- dfdetrended


for (i in 14:19){
  channel <- dfdetrended[,i]
  dfMMG[,i] <- butterHighPassFilter(channel, cutoff = 2, Fs = 50)  # MMG signals
}

colHeader <- colnames(dfMMG)
dfplot1<-data.frame(time = dfdetrended$X.Time., channel, dfMMG[,i], dfdetrended$X.PULSE*max(channel)*scaling, eventDetect*max(channel)*scaling)
dygraph(dfplot1, main = colHeader[i]) %>%
  dyAxis("x", label="Time (s)") %>% 
  dyAxis("y", label="Angular velocity") %>%
  dyOptions(colors = c('black', 'blue', 'red')) %>% 
  dyLegend(width = 130) %>%
  dyRangeSelector() 





# Partitioning Tasks ------------------------------------------------------

dflist <-partitionDataFrame(dfMMG, eventDetect)
dfpinch <- dflist[[1]]
dfhand <- dflist[[2]]
dfpron <- dflist[[3]]
dfflex <- dflist[[4]]

ggplotAccelerometer(dfpinch)



# Saving MMG Segments -----------------------------------------------------

write.xlsx(dfpinch, paste(pathData,"PinchMMG.xlsx",sep=""))
write.xlsx(dfhand, paste(pathData,"HandMMG.xlsx",sep=""))
write.xlsx(dfpron, paste(pathData,"PronMMG.xlsx",sep=""))
write.xlsx(dfflex, paste(pathData,"FlexMMG.xlsx",sep=""))




