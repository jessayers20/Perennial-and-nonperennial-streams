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
met1 <- read.csv("C:\\Work\\Data\\Functional Flows\\Met files\\met_dryseason_point1_5condays_15percentofyears_reference.csv")

# key variables for merging 
key=met
key$siteyear <- paste0(key$COMID,"_",key$Year)

keyvars <- key %>% 
  dplyr::distinct(siteyear, .keep_all = T)


refs = unique(met[,c("COMID","ID")])

ids <- data.frame(unique(rbind(refs,nrefs,altered)))
ids$COMID = as.numeric(ids$COMID)

files <- list.files("D:\\CA_NHDPreds_1950_2023")
comids <- as.numeric(str_remove(files, "\\.csv"))


# need comids for non-reference sites
nonrefs <- read.csv("C:\\Work\\Data\\Functional Flows\\Met files\\Met_altered_analysis_12_2023.csv",header=T)
nonrefs = nonrefs[which(nonrefs$Stat == "Streamclass_annual"),]
nref=unique(nonrefs$ID)
nrefs = unique(nonrefs[,c("COMID","ID")])

ids = rbind(refs,nrefs)

met = read.csv("C:\\Work\\Data\\Functional Flows\\Met files\\Met_Master_Vars_Altered_2023.csv",header=T)

# need file that links comids and site ids 
## here

library(plyr)
#files <- list.files("C:\\Work\\Data\\Altered Low Flow Metrics Dry Season Consecutive\\")
#sites <- unique(ldply(files,function(x) strsplit(x,"_")[[1]][1])$V1)

ls.list <- list()

sites = unique(ids$ID)

for(s in 1:length(sites))  {
  
  site = sites[s]
  #within original temp, extracts into characters separated by / and extracts in the fourth position the characters in the position of 1 to 7

  file = paste0("C:\\Work\\Data\\Altered Low Flow Metrics Dry Season Consecutive\\",site,"_Annual_low_flow.csv") 
  

  if(file.exists(file)){
  
  comid = unique(ids[which(ids$ID == site),]$COMID)[[1]]
    
  temp1 = read.csv(paste0("C:\\Work\\Data\\Altered Low Flow Metrics Dry Season Consecutive\\",site,"_Annual_low_flow.csv"),header=T)
  
  temp1 = temp1[,c(1:6)]
  temp1$site = paste0(site)
  
  avg <- mean(temp1[which(!is.na(temp1$Zero_Flow_Days)),]$Zero_Flow_Days)
  
  kvar = met[which(met$COMID == comid),]
  
  if(nrow(kvar) == 0){
    kvar = met1[which(met1$COMID== comid),]
    kvar = kvar[,-c(2:5)]
    kvar = kvar[which(colnames(kvar) %in% colnames(met))]
  }
  kvar = unique(kvar)
  

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
  lff2$COMID = paste0(comid)
  
  ls.list[[s]] <- (lff2)
  rm(lf,lfd,lf2,lfd2,avg,temp1,lff,ipvals,classdf)
  
  }
}

lowdf <- rbindlist(ls.list)

write.csv(lowdf,"C:\\Work\\Data\\Functional Flows\\Met files\\Met_altered_analysis_1_2024.csv",row.names=F)


# Read in data for reference years
#all.eff <- read.csv("C:/Work/Data/Functional Flows/CA_USGS_Gage_Reference_Screen_December_2023_WRR.csv",header=T); head(all.eff)
# Need to run metfile again bc missing some COMIDS

# Read in data for reference years
all.eff <- read.csv("C:/Work/Data/Functional Flows/CA_Reference_Screen_December_2023_WRR_20years.csv",header=T); head(all.eff)
ffref <- unique(all.eff$ID)


ref.met<-list()

for(s in 1:166){
  
  site = ffref[s]
  #site = 11181390   
  # site = 11480390#"11120550" "11153470" "11220000" "11220500" "11269300" "11413320" "11413323" "11426150" "11433260" "11460100" "11467200" "11472900"
  info <- all.eff[which(all.eff$ID == site),]
  mfile = lowdf[which(lowdf$ID == site),]
  min = info$Min_year
  max = info$Max_year
  
  # check to see if site is in the altered analysis
  analy <- altered[which(altered$ID == site),]
  
  # assign current reference or nonreference
  # assign none if its not in the current analysis
  analy <- ifelse(nrow(analy)==0,analy <- "None",analy <- analy$Ref_status)
  
  if(analy == "Ref"){
    mfile=mfile[which(mfile$Year < 1980),]
  }
  ref.met[[s]]<-mfile
}

refmet <- data.frame(rbindlist(ref.met))
refmet <- refmet[!is.na(refmet$Value),]

setdiff(as.character(ffref),unique(refmet$ID))

