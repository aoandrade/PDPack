# open file
# estimate envelope
# detect medium
# save info into table

# DataFrame header: subjectCode, filename, date, movementType, channel, burstIdx, medianMVC, filterInfo

# Load all necesssary libraries
libPath  = "/Users/alicerueda/Documents/RStudio/Adriano/"
setwd(libPath)
source("read_Intan_RHD2000_file.R")
source("Thresholding.r")


#subjectCode: get it from the directory
pathData ="/Users/alicerueda/Documents/RStudio/Adriano/EMGDATA/"
subjectCode = "EJRR"
pathData = paste(pathData,subjectCode,sep = '')
setwd(pathData)
f <- sort(list.files(path=pathData, full.names=TRUE, pattern="*.xlsx"))
filename <- basename(f)
f <- f[which(grepl("^CVM", f)==TRUE)]

DEBUG = FALSE
if (DEBUG) {
  filename <- file.choose()
  filename <- basename(filename)[[4]]
  percentage = 0.4
  PLOT = TRUE
  chanList = list(2,3,4)
  i <- 3
}


dfSave <- getMVC(filename[[7]], chanList = list(2,3,4), percentage = 0.4,PLOT=TRUE)


