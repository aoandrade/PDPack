# Example to call getEMGTimeFeature, normalizedEMG, and save it to file
source("EMGLib.R")

## Constructing a reference table from the MVC reference file
# Opening the reference file for the MVC values and construct a reference table for normalization
refFilename <- file.choose()  #choose the MVC reference to get the cutoffFreq

maxRefTable<-getEMGMaxRefTable(refFilename)


# Working on the actual EMG signal ----------------------------------------

filename <- file.choose()
df <- readWorkbook(filename,
                   sheet = 1,
                   detectDates = TRUE)

Results<- getEMGTimeFeatures(df, PLOT=FALSE)


# Normalize the Results using the maximum value of the channel in maxRefTable
normResults<-getEMGNormResults(Results, maxRefTable)

