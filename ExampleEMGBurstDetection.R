source("EMGLib.R")
source("Thresholding.R")

# Opening a selected EMG signal file
filename <- file.choose()
df <- readWorkbook(filename,
                   sheet = 1,
                   detectDates = TRUE)

# Event detection and segmentation
emgBurst <- detectEMGBurst(df, cutoffFreq=1, percentage = 0.08, PLOT=TRUE)

