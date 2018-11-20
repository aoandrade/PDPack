# Example of extracting information from the EMG signal
# Results contains all the peak and valley locations, and distances between them

source("EMGLib.R")

# Working on the actual EMG signal ----------------------------------------
# Opening a selected EMG signal file
filename <- file.choose()
df <- readWorkbook(filename,
                   sheet = 1,
                   detectDates = TRUE)

Results<- getEMGTimeFeatures(df, PLOT=TRUE)


# Save Results to file












