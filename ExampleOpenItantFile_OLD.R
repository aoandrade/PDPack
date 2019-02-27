
# Load File ---------------------------------------------------------------

# Avaliable @  https://github.com/aoandrade/NIATS/blob/master/TremsenToolbox/TREMSENToolbox.R
source("TREMSENToolbox.R")


source("read_Intan_RHD2000_file.R")



filename <- file.choose() # Seleção do arquivo
X <- OpenIntanFile(filename) # Open an Intan file

chan <- 1 # coloque aqui o canal a ser visualizado

dfplot <- data.frame(time=X$t_amplifier, y=X$amplifier_data[chan,], yd=X$board_dig_in_data*max(X$amplifier_data[chan,]))
dygraph(dfplot)%>%dyRangeSelector()


# Convert intan files to Excel------------------------------------------------------------------------

# Convert all intan files in the folder where filename is located to an Excel file. 
# File merging is executed so that split data (saved in different intan files) are put together
filename <- file.choose() # Seleção do arquivo .rdh (intan)
ConvertRHD2Excel(filename)



# Read Excel file ---------------------------------------------------------

out_xlsx <- file.choose() # Seleção do arquivo Excel
df <- readWorkbook(out_xlsx,
                       sheet = 1,
                       detectDates = TRUE)


Fs <- 1/(df$time[2]-df$time[1]) # sampling frequency in Hz


# Getting signal blocks/windows
i1 <- which(diff(df$Pulse) == 1)
wndIndx <- matrix(i1,length(i1)/2,2, byrow= TRUE)


# Plots
y <- df$chan.1
dfplot <- data.frame(time=df$time, y=y, yd= df$Pulse * max(y))
dygraph(dfplot)%>%dyRangeSelector()%>%
dyEvent(wndIndx[1:size(wndIndx)[1],1]/Fs, labelLoc = "bottom")%>%
dyEvent(wndIndx[1:size(wndIndx)[1],2]/Fs, labelLoc = "bottom")



# Nonlinear trend detection -----------------------------------------------

loessData <- function(yamp,t) {
  dat <- data.frame(x = t, y = yamp)
  yp = predict(loess(y ~ x, dat, span = 0.1))
  return(yp)
}

seg <- 3 #select the segment (window/signal block)
ysig <- df$chan.1[wndIndx[seg,1]:wndIndx[seg,2]]
tsig <- seq(1,length(ysig), by=1) * Fs

# detrending ...
ydet <- detrend(ysig, tt = 'linear')
yn <- loessData(ydet, tsig)
ydetrended <- detrend(ydet-yn, tt = 'linear')



dfplot <- data.frame(time= tsig,ysig, ydetrended)
dygraph(dfplot) %>% dyOptions(colors = c('red', 'blue'))




# Processing --------------------------------------------------------------
yenv <- abs(seewave::hilbert(ydetrended, f = Fs, fftw = TRUE))


#############################################
Fn <- Fs/2
f <- c(1,10)
per <- f/Fn

bf <- butter(3, per, type='pass')
s_filt <- filtfilt(bf, yenv-mean(yenv))

psx <- pspectrum(s_filt,verbose = TRUE, ntap.init = 10,
                 niter=100, AR=TRUE, x.frqsamp=Fs, plot=FALSE) ##library(psd)
dygraph(data = data.frame(time=psx$freq, psx$spec)) %>%
  dyRangeSelector(dateWindow = c(0,10))

######################################################

Fn <- Fs/2
f <- psx$freq[which.max(psx$spec)]
per <- f/Fn

bf <- butter(3, per, type='low')
s_filt <- filtfilt(bf, yenv-mean(yenv))

dfplot <- data.frame(time=tsig, ydetrended, s_filt)
dygraph(dfplot)%>%dyRangeSelector() %>% 
  dyOptions(colors = c(rgb(0.8,0.8,0.4), 'red'))


# PS analysis -------------------------------------------------------------
s <- s_filt





  #dyAxis("x", label = "Frequency", valueRange = c(40, 60)) 

# Filtering ---------------------------------------------------------------

