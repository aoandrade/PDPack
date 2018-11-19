source("read_Intan_RHD2000_file.R")


# Convert intan files to Excel------------------------------------------------------------------------

# Convert all intan files in the folder where filename is located to an Excel file. 
# File merging is executed so that split data (saved in different intan files) are put together
filename <- file.choose() # Seleção do arquivo .rdh (intan)
ConvertRHD2Excel(filename)
