# Example for calculating MVC for the person and save all the reference values into an excel file

# Load all necesssary libraries
libPath  = "D:/Adriano/OneDrive/Projetos/2017-Brazil-Canada/DATA/PROGRAMS/PDPack/"
source("read_Intan_RHD2000_file.R")
source("Thresholding.r")

#subjectCode: get it from the directory
pathData ="D:/Adriano/OneDrive/Projetos/2017-Brazil-Canada/DATA/EJRR/09Nov2018/EMG/Excel"
subjectCode = "EJRR"
#pathData = paste(pathData,subjectCode,sep = '')
setwd(pathData)
filename <- sort(list.files(path=pathData, full.names=TRUE, pattern="*.xlsx"))
filename <- basename(filename)
filename <- filename[which(grepl("^CVM", filename)==TRUE)]     # Only selecting the files that start with filename "CVM"

# Example for processing one file
dfSave <- getMVC(filename[[9]], subjectCode, chanList = list(2,3,4), cutoffFreq=1, percentage = 0.5,PLOT=TRUE)

# Example for processing all files in the directory
#dfSave <- lapply(filename, FUN=getMVC, subjectCode, chanList = list(2,3,4), cutoffFreq=1, percentage = 0.5,PLOT=FALSE)

df<-do.call("rbind", list(dfSave, stringsAsFactors=FALSE))
#df<-t(df)

WriteDF2ExcelFile(out_xlsx="MVCReferenceValue.xlsx", df, worksheetName="MVC Reference")
