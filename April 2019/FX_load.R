# Quandl provides spot FX from BOE (to GBP), FRED (to USD) and ECB (to EUR)
# See https://blog.quandl.com/api-for-currency-data for cryptical symbols



rm(list=ls())



# G10 currencies vs. USD i.e. 9 pairs
#   Name conversion: used in loaded data
asset_vec_FRED_transl<-c("AUDUSD","GBPUSD","USDCAD","EURUSD","USDJPY","NZDUSD","USDNOK","USDSEK","USDCHF")
#   FRED symbols, see https://blog.quandl.com/api-for-currency-data
asset_vec_FRED_quandl<-c("DEXUSAL","DEXUSUK","DEXCAUS","DEXUSEU","DEXJPUS","DEXUSNZ","DEXNOUS","DEXSDUS","DEXSZUS")

series_i_vec<-1:length(asset_vec_FRED_quandl)



# Data
# Daily EURUSD has been included in MDFA-package (18-02-2017)
# But one can refresh the data from Quandl by selecting refresh_data<-T
#   Just plug-in the API key in the corresponding code line

refresh_data<-T


# Libraries: Quandl is required for refreshing the data
if (refresh_data)
  library(Quandl)

library(xts)
#install.packages("devtools")
library(devtools)
# Load MDFA package from github
devtools::install_github("wiaidp/MDFA")
# MDFA package: EURUSD is now part of the data in the package
library(MDFA)


#---------------------------------------------------------------------
# Load data
# Daily EURUSD has been included in MDFA-package (18-02-2017)
# But one can refresh the data from Quandl by selecting refresh_data<-T
#   Just plug-in the API key in the corresponding code line
if (refresh_data)
{
  
  QuandlAPIKey <- "set in here your API key"
  QuandlAPIKey <- "cxvP9tfN13a9tvMmsXag"
  
  Quandl.api_key(QuandlAPIKey)
  
  for (i in 1:length(series_i_vec))#i<-1
  {
  
    FX_pair<-asset_vec_FRED_quandl[series_i_vec[i]]
    
    Data <- Quandl(paste("FRED/",FX_pair,sep=""))
#    Data <- Quandl(paste("ECB/",FX_pair,sep=""))
#    Data <- Quandl("CUR/GBP")
    FX<-Data[,2]
    names(FX)<-Data[,1]
    FX<-as.xts(FX)
    if (i==1)
    {
      FX_mat<-FX
    } else
    {
      FX_mat<-cbind(FX_mat,FX)
    } 
    print(i)
  }
  colnames(FX_mat)<-asset_vec_FRED_transl[series_i_vec]
  tail(FX_mat)
  head(FX_mat)
  save(FX_mat,file=paste(getwd(),"/FX_mat.Rdata",sep=""))
} else
{
# This daily EURUSD from MDFA-package
  tail(FX_mat)
}



