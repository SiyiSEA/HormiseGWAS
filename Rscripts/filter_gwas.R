###############################description###################################


#################################packages####################################
suppressPackageStartupMessages(library(data.table))
options(bitmapType = "cairo")


#################################main code####################################
arguments <- commandArgs(T)
setpath <- arguments[1]
filename <- arguments[2]
phenotype <- arguments[3]
MAFcol <- arguments[4]
INFOcol <- arguments[5]

setwd(setpath)
print(paste("Working directory set to:", setpath))
GWAresult = fread(filename, header = T, data.table=F, stringsAsFactors=F)

if (!is.na(MAFcol)) {
  MAFcol <- as.numeric(MAFcol)
  print(paste("Filtering the gwas stats based on the MAFcol > 0.01", MAFcol))
  GWAresult  = subset(GWAresult, GWAresult[,MAFcol] > 0.01)
} 

if (!is.na(INFOcol)) {
  INFOcol <- as.numeric(INFOcol)
  print(paste("Filtering the gwas stats based on the INFOcol > 0.8", INFOcol))
  GWAresult  = subset(GWAresult, GWAresult[,INFOcol] > 0.8)
} 

fwrite(GWAresult, file = paste0(phenotype,".LDSCinput"), quote = F, sep = "\t", col.names = T, row.name = F )