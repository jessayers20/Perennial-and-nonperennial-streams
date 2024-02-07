#######################################
#####     Updated Low Flow        #####
#####         7-day mo avg        #####
#######################################

# This script calculates the 7-day moving average, zero flow and flow duration curve thresholds
# over the California dry season (June 1 to December 31)
# First section provides the output (csv file) of daily discharge and rolling mean

# Second section calculates the number of zero flow days, 7-day minimum, and flow duration curve thresholds
# Calculates the timing of each low flow metric

cat("\014") ;  rm(list=ls())

library(dataRetrieval); library(FlowScreen); library(tidyverse)
library(dplyr); library(lubridate); library(tidyr); library(ggplot2);library(tools)
library(spatial); library(rgeos); library(raster); library(data.table)
library(RColorBrewer); library(rgdal); library(lfstat); library(plyr)

files <- list.files("C:\\Work\\Data\\USGS gages\\California Discharge\\Daily_Qavg_until_12_2023\\")
sites2 <- unique(ldply(files,function(x) strsplit(x,"_")[[1]][1])$V1)

# Read in data for reference years
all.eff <- read.csv("C:/Work/Data/Functional Flows/CA_Reference_Screen_December_2023_WRR_20years.csv",header=T); head(all.eff)
ffref <- unique(all.eff$ID)

met1 <- read.csv("C:\\Work\\Data\\Functional Flows\\Met files\\met_dryseason_point1_5condays_15percentofyears_reference.csv")


altered <- read.csv("C:\\Work\\Data\\USGS gages\\California Discharge\\CA_altered_1980_25years_2.csv",header = T)
altered <- altered[which(altered$Max_1980 >= 2020),]
altered <- altered[which(altered$Npost_1980 > 40),]
alsites = unique(altered$ID)
curnet <- alsites[which(alsites %in% sites2)]

current <- altered[which(altered$Ref_status == "Ref"),]$ID

sites2 = c(ffref,curnet)
sites2 = ffref

# list to store data
seven.list <- list()

for(s in 1:length(sites2)){
  
  #s=20
  site = sites2[s]
  print(site)
  
  
  sf <- read.csv(paste0("C:\\Work\\Data\\USGS gages\\California Discharge\\Daily_Qavg_until_12_2023\\",as.numeric(site),"_Discharge_12_2023.csv"), header=T)
  sf = sf[,1:2]
  
  #calculate stats for the Zero-Flow days for each year and EFF period
  library(lfstat)
  sf$Date <- as.Date(sf$Date, "%Y-%m-%d")
  sf$MD <- substr(sf$Date,6,10)
  sf$year <- as.numeric(substr(sf$Date,1,4))
  sf$Wyear <- water_year(as.Date(sf$Date), origin=10)
  
  colnames(sf) <- c("Date","Discharge","MD","year","Wyear")
  
  #Check for consecutive dates
  min <- min(sf$year);  max <- max(sf$year)
  
  wateryears <- seq(min,max,by=1)
  actyears <- unique(sf$year)
  
  #test to see how many missing observations  
  missing <- length(sf[is.na(sf$Discharge),])
  
  #data frame for rolling mean  
  x <- sf
  
  #Calculates moving average for different window size (n)
  x$roll7Mean <- RcppRoll::roll_mean(x$Discharge, n = 7L, fill = NA)
  x$test <- paste0(missing)
  x$site <- paste0(site)
  
  #remove rows with nas and days that do not have consecutive values
  all.sev = x[!is.na(x$roll7Mean),]
  all.sev <- na.omit(x)
  
  all.sev$year <- as.numeric(substr(all.sev$Date,1,4))
  all.sev$MD <- substr(all.sev$Date,6,10)
  all.sev$Wyear <- water_year(as.Date(all.sev$Date), origin=10)

  
  # Data frame to store site's annual low flow metrics
  intermets <- data.frame(matrix(NA,nrow=length(wateryears),ncol=24))
  
  threshold.list <- list()
  
  for(y in 1:length(wateryears)){
    
    #y=6
    year = wateryears[y]
    
    #extract water year
    daily <- all.sev[which(all.sev$Wyear==year),]
    daily$Day <- seq_len(nrow(daily))
    
    
    if(nrow(daily)>=330){ # no more than 10% of missing data
      
      # to calculate low flow metrics over dry season, need next year's wy 
      nextyear <- all.sev[which(all.sev$Wyear == as.numeric(year)+1),]
      nextyear <- nextyear[order(as.Date(nextyear$Date, format="%Y-%m-%d")),]
      nextyear$Day <- seq_len(nrow(nextyear))+nrow(daily)
      
      #April 1 to Jan 1 for dry season timing 
      #Dry Season -- Calculates number of No Flow days during the dry season
      dry.tim <- as.Date(paste0("",year,"-04-01"),format="%Y-%m-%d")
      dry.end <- as.Date(paste0("",year,"-12-31"),format="%Y-%m-%d")
      
      # get entire ds (into next water year too)
      dry1 <- daily[which(daily$Date >= dry.tim),]
      dry2 <- nextyear[which(nextyear$Date <= dry.end),]
      dry <- rbind(dry1,dry2)
      
      # plots to check data
      library(ggplot2)
      t <- ggplot(daily)+
        geom_line(aes(x=Date,y=Discharge,color="Daily streamflow"))+
        #geom_line(aes(x=Date,y=roll7Mean,color="7-day mov avg"))+
        geom_line(data=dry,aes(x=Date,y=Discharge,color="Dry season"))+
        scale_color_manual(values=c("Daily streamflow"="black","7-day mov avg"="blue",
                                    "Dry season"="red"),name="")+
        theme_light()
      t <- t + theme(legend.position = c(0.8,0.88))
      #print(t)
      
      
      library(dplyr)
      #Calculate sev day min, duration and start date
      SevDmin <- min(dry$roll7Mean)
      Sevdur <- ifelse(nrow(dry)>0, as.numeric(count(dry[which(dry$Discharge == SevDmin),])$n), NA)
      SevDay <- dry[which(dry$roll7Mean == SevDmin),]$Date[1]
      SevStart <- dry[which(dry$roll7Mean == SevDmin),]$Day[1]
      
      # number of zero flow (<0.1) within the dry season (raw)
      fzeros <- ifelse(nrow(dry)>0, as.numeric(nrow(dry[which(dry$Discharge < 0.1),])), NA)
      zdf <- dry[which(dry$Discharge < 0.1),]
      
      # if there are any zero flow days
      if(nrow(zdf)>0){
        
        zeros <- zdf %>%
          mutate(date = as.Date(Date, format = "%Y-%m-%d")) %>%
          arrange(Date) %>%
          group_by(id2 = cumsum(c(T, diff(Date) > 1)), .add = T) %>%
          mutate(num_con_days = ifelse(date == first(Date), last(Date) - Date + 1, 0))
        
        #zeros$Date <- as.Date(zeros$Date, format = "%Y-%m-%d")
        
        first <- zeros[which(zeros$num_con_days >= 5),]
        first <- zdf[1,]$Day
        
        # number of zero flow days after first five consecutive days
        fzeros1 <- nrow(zdf[which(zdf$Day >= first),])
        
        
        
      }  
      else{
        first <- NA
        last <- NA
        duration <- NA
        fzeros1 <- 0
      }
      
      print(fzeros1-fzeros)
      
      intermets[y,1] <- year
      intermets[y,2] <- fzeros1
      intermets[y,3] <- first
      intermets[y,4] <- SevDmin
      intermets[y,5] <- SevStart
      intermets[y,6] <- paste0(SevDay)
      
      # get threshold (discharge) and number of days below threshold (count of days)
    
      
    }
    else{next}
  }
  
  
  
  colnames(intermets) <- c("Year","Zero_Flow_Days","Zero_Flow_Start_Date_WY",
                           "Min_7_Day_Mov_Avg","Zero_date_check", "Min_7_Start_Date_WY")

  write.csv(intermets,paste0("C:\\Work\\Data\\Altered Low Flow Metrics Dry Season Consecutive\\",site,"_Annual_low_flow.csv"),row.names = F)
  rm(sf,intermets,nextyear,dry,dry1,dry2,daily,t,x,threshold.list,all.sev, max,min, missing, s, SevDay, SevDmin, Sevdur,SevStart, first, fzeros, fzeros1,last,site,y,year,wateryears,dry.end,dry.tim,duration,actyears,zdf,zeros)
}




# Calculate SF Class
# need to do separately for reference and non-reference periods
#current netword with comids
altered <- read.csv("C:\\Work\\Data\\USGS gages\\California Discharge\\CA_altered_1980_years.csv",header = T)
curnet <- unique(altered$ID)

curref <- altered[which(altered$Ref_status == "Ref"),]

# file for predictors
tmet <- read.csv("C:\\Work\\Data\\Functional Flows\\Met files\\met_dryseason_point1_5condays_15percentofyears_reference.csv",header=T)
tmet = tmet[which(tmet$Stat == "Streamclass_annual"),]

# remove tmet observations that have current reference gage data (1980-2015)
tmet2 <- tmet[which(tmet$ID %in% curref$ID),]
tmet2 <- tmet2[-which(tmet2$Year >= 1980),]

tmet1 <- tmet[-which(tmet$ID %in% curref$ID),]
tmet <- rbind(tmet1,tmet2)


all.eff <- read.csv("C:\\Work\\Data\\Functional Flows\\CA_USGS_Gage_Reference_Screen_December_2023_WRR.csv",header=T)
tmet = tmet[which(tmet$ID %in% all.eff$ID),]
length(unique(tmet$ID))
ref.list <- list()

#data frame to save values
zavs <- data.frame(matrix(data=NA, nrow=length(sites),ncol=7))

for(s in 1:length(sites)){
  
  # site =11180960
  #s=1
  site = sites[s]
  temp1=read.csv(paste0("C:\\Work\\Data\\Low Flow Metrics Dry Season Consecutive\\",site,"_Annual_low_flow.csv"),header=T)
  temp1$site = paste0(site)
  
  #get ref years
  ref.info = all.eff[which(all.eff$ID == site),]
  
  # if reference gage is in the modern, altered anlaysis
  if(site %in% current){maxyear = 1980}
  temp1 = temp1[which(temp1$Year >= minyear),]
  temp1 = temp1[which(temp1$Year <= maxyear),]
  
  temp1 = temp1[,c(1:6,25)]
  
  # Reference period means 
  amean <- mean(temp1$Zero_Flow_Days)
  amedian <- median(temp1$Zero_Flow_Days)
  
  # number/percent of years with no flow
  occur <- nrow(temp1[which(temp1$Zero_Flow_Days > 0),])
  percent <- occur/nrow(temp1)
  
  # number/percent of years with no flow greater than 15 days
  occur5 <- nrow(temp1[which(temp1$Zero_Flow_Days >= 5),])
  percent5 <- occur5/nrow(temp1)
  
  nyears <- nrow(temp1)
  
  zavs[s,1] <- paste0(site)
  zavs[s,2] <- amean
  zavs[s,3] <- amedian
  zavs[s,4] <- percent
  zavs[s,5] <- percent5
  zavs[s,6] <- nyears
  zavs[s,7] <- occur5
  zavs[s,8] <- ref.info$COMID
  
  
  # list to show when years of reference occur for analysis 
  
  ref.list[[s]] <- temp1
  
  rm(occur,occur5,percent5,percent,amean,amedian,temp1,site,s)
}

zavs=zavs[!is.na(zavs$X1),]
colnames(zavs) <- c("site_no","Mean_Zero","Median_Zero","Percent","Percent5","N_years","N_zero_flow_years","COMID")
zavs$Overall_Class <- ifelse(zavs$Percent5 > 0.15, "Non-perennial","Perennial")

all.ref <- data.frame(rbindlist(ref.list))

counts <- all.ref %>%
  group_by(Year) %>%
  dplyr::summarise(count=n())

analysisplot <- ggplot(data = all.ref, aes(x = Year)) +
  ylab("No. of sites")+
  geom_histogram(bins=length(unique(all.ref$Year)),fill="#31a354",color="darkgreen")+
  theme_light()

ggsave("All_reference_years_altered_removed.png", plot=analysisplot, device="png", path = "C:\\Work\\CEFF Metrics and Models\\Figures\\",
       width = 5, height = 5, dpi=300, limitsize=TRUE)

nonp <- zavs[which(zavs$Overall_Class == "Non-perennial"),]
