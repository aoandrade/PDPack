source("read_Intan_RHD2000_file.R")


filename <- file.choose() # File selection
X <- OpenIntanFile(filename) # Open an Intan file

chan <- X$chan.2 # Select the channel to be visualized

# plot it
dfplot <- data.frame(time=X$time, y=chan, yd=X$pulse*max(chan))
dygraph(dfplot)%>%dyRangeSelector()
