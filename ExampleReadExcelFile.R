source("read_Intan_RHD2000_file.R")


# Read Excel file ---------------------------------------------------------

out_xlsx <- file.choose() # Seleção do arquivo Excel
df <- readWorkbook(out_xlsx,
                   sheet = 1,
                   detectDates = TRUE)



y <- df$chan.1 # Select channel to plot

Fs <- 1/(df$time[2]-df$time[1]) # sampling frequency in Hz


# Getting signal blocks/windows
i1 <- which(diff(df$pulse) == 1)

if ( isempty(i1) == TRUE)  {
  
  fprintf("File without annotation \n")
  
  dfplot <- data.frame(time=df$time, y=y)
  dygraph(dfplot)%>%dyRangeSelector()
  
} else {
  
  fprintf("plotting data \n")
  
  wndIndx <- matrix(i1,length(i1)/2,2, byrow= TRUE)
  
  
  
  dfplot <- data.frame(time=df$time, y=y, yd= df$pulse * max(y))
  dygraph(dfplot)%>%dyRangeSelector()%>%
    dyEvent(wndIndx[1:size(wndIndx)[1],1]/Fs, labelLoc = "bottom")%>%
    dyEvent(wndIndx[1:size(wndIndx)[1],2]/Fs, labelLoc = "bottom")
  
}

