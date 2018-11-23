# Example for calculating MVC for the person and save all the reference values into an excel file

# Load all necesssary libraries
source("read_Intan_RHD2000_file.R")
source("Thresholding.R")
source("EMGLib.R")

#subjectCode: get it from the directory
#pathData ="D:/Adriano/OneDrive/Projetos/2017-Brazil-Canada/DATA/EJRR/09Nov2018/EMG/Excel"
pathData = "/Users/alicerueda/Documents/RStudio/Adriano/DATA/"
# defining subject and filetype
subjectCode = "EJRR"
signalType = "EMG"
fileType = "Excel"

# set working directory to the target files
pathData = paste(pathData,subjectCode,'/',sep = '')
listDir <- list.dirs(pathData, recursive = FALSE)
dateRecorded <- paste(basename(listDir[1]),'/', sep='')
pathData = paste(pathData, dateRecorded, signalType, '/', fileType, '/', sep='')

# list all relevant files
filename <- sort(list.files(path=pathData, full.names=TRUE, pattern="*.xlsx"))
filename <- basename(filename)
filename <- filename[which(grepl("^CVM", filename)==TRUE)]     # Only selecting the files that start with filename "CVM"

# Example for processing one file
#dfSave <- getMVC(filename[[4]], pathData, subjectCode, chanList = list(2,3,4), cutoffFreq=1, percentage = 0.5,PLOT=TRUE)

# Example for processing all files in the directory
dfSave <- lapply(filename, FUN=getMVC, pathData, subjectCode, chanList = list(2,3,4), cutoffFreq=1, percentage = 0.5,PLOT=TRUE)

df<-do.call("rbind",dfSave)

# Save all MVC to an excel file
WriteDF2ExcelFile(out_xlsx=paste(pathData,"MVCReferenceValue.xlsx",sep=""), df, worksheetName="MVC Reference")
