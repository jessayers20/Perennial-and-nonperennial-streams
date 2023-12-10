##############################
#############################
## Combine outputs 
## Make file that goes into model

# Takes the outputs of formated files from scripts (1+2) PRISM datat
# and (3) Streamflow metrics
# merges to one file for input into the RF model

cat("\014") ;  rm(list=ls())

options(max.print=9999999)

library(VSURF); library(varSelRF); library(caret); library(Boruta); library(gamlss)
library(randomForest);library(dplyr);library(data.table);library(foreach)
library(ranger);library(doParallel); registerDoParallel(cores=4); library(ggplot2)
library(magrittr)

# Select low flow metric --------------------------------------------------

## Load data with all FFM observations and associated watershed variables
# Get site ids and comids
met <- read.csv("C:\\Work\\Data\\Functional Flows\\Met files\\Met_Master_Vars_Altered_2022.csv")
met <- rename(met, Year=WaYr)


#key =  newmet[,-c(3,5)]
key=met
key$siteyear <- paste0(key$COMID,"_",key$WaYr)

keyvars <- key %>% 
  dplyr::distinct(siteyear, .keep_all = T)

altered <- read.csv("C:\\Work\\Data\\USGS gages\\California Discharge\\CA_altered.csv",header = T)
altsites = unique(altered$ID)
altered = altered[,c("COMID","ID")]


refs <- read.csv("C:\\Work\\Data\\Functional Flows\\Met files\\met_dryseason_point1_5condays_15percentofyears_reference.csv",header=T)
refs = refs[which(refs$Stat == "Streamclass_annual"),]
refsites=unique(refs$ID)
refs = unique(refs[,c("COMID","ID")])


idcom <- rbind(altered,refs)

# need file that links comids and site ids 
## here


sites <- c(refsites,altsites)
sites = sites[-which(sites %in% c("11297200","11465350"))]

ls.list <- list()

for(s in 1:length(sites))  {
  
  site = sites[s]  #within original temp, extracts into characters separated by / and extracts in the fourth position the characters in the position of 1 to 7
  
  temp1 = read.csv(paste0("C:\\Work\\Data\\Altered Low Flow Metrics Dry Season Consecutive\\",site,"_Annual_low_flow.csv"),header=T) 
  temp1$site = paste0(site)
  
  avg <- mean(temp1$Zero_Flow_Days)
  
  comid = idcom[which(idcom$ID ==site),]$COMID
  
  kvar = met[which(met$COMID == comid),]
  

  # merge annual low flow metrics with climate/watershed variables as input in RF models
  lf <- temp1[,c("site","Year","Zero_Flow_Days")]
  lf$Metric <- paste0("Zero-flow") 
  colnames(lf) <- c("ID","Year","Value","Stat")
  lf2 <- merge(lf,kvar)
  
  lfd <- temp1[,c("site","Year","Zero_Flow_Start_Date_WY")]
  lfd$Metric <- paste0("Zero_Start")
  colnames(lfd) <- c("ID","Year","Value","Stat")
  lfd2 <- merge(lfd,kvar)
  
  
  #} else{
  minmag <- temp1[,c("site","Year","Min_7_Day_Mov_Avg")]
  minmag$Metric <- paste0("Min-7")
  colnames(minmag) <- c("ID","Year","Value","Stat")
  minmag <- merge(minmag,kvar)
  
  minday <- temp1[,c("site","Year","Min_7_Start_Date_WY")]
  minday$Metric <- paste0("Min-7_Start")
  colnames(minday) <- c("ID","Year","Value","Stat")
  minday <- merge(minday,kvar)
  
  #}
  
  lff <- rbind(lf2,lfd2,minmag,minday)
  lff = lff[,c(5,6,4,2,3,1,7:length(lff))]
  
  #set annual scale of zero flow vals
  ipvals <- temp1[,1:2]
  ipvals$Value <- ifelse(ipvals$Zero_Flow_Days >= 5, "Intermittent","Perennial")
  ipvals$Stat <- paste0("Streamclass_annual")
  ipvals$ID <- paste0(site)
  ipvals <- ipvals[c("ID","Year","Value","Stat")]
  classdf <- merge(ipvals,kvar)
  
  lff2 <- rbind(classdf,lff)
  
  ls.list[[s]] <- (lff2)
  rm(lf,lfd,lf2,lfd2,avg,temp1,lff,ipvals,classdf)
  
}

lowdf <- rbindlist(ls.list)

write.csv(lowdf,"C:\\Work\\Data\\Functional Flows\\Met files\\Met_altered_analysis_12_2023.csv",row.names=F)





