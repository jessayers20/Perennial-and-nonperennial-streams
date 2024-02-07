######################################################################
#####                  Model low-flow metrics                    #####
#####                                                           #####
######################################################################

# Reference gages to model low flow metrics using random forests
# minimum 7-day moving average and timing
# zero flow days is not predicted well enough to include in CEFF


cat("\014") ;  rm(list=ls())

options(max.print=9999999)

## Begin modeling
library(randomForest);library(dplyr);library(data.table);library(foreach); library(stats)
library(ranger);library(doParallel); registerDoParallel(cores=4); library(ggplot2); library(caret)
library(ggplot2); library(tidyverse); library(VSURF); library(varSelRF); library(caret); library(Boruta)


## Load data with all FFM observations and associated watershed variables
#drymets <- c("DS_Mag_90","DS_Mag_50","Min-7", "Min-7_Start","DS_Dur_WS","DS_Tim")
timing <- c("Min_7_Date","Zero-Flow","Zero_Date","DS_Tim","DS_Dur_WS","Streamclass_annual")

# file for predictors
met <- read.csv("C:\\Work\\Data\\Functional Flows\\Met files\\Met_altered_analysis_1_2024.csv",header=T)

# To create model, need to get reference sites from original MET file

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


# loop through all functional flows here or just select which metrics using drymets
#for(m in 1:length(drymets)){
#

curmet="Streamclass_annual"
dir.create(paste0("C:\\Work\\CEFF Metrics and Models\\Results\\Model Performance\\Stream Classification\\WRR Update\\" ,curmet,""))
dir.create(paste0("C:\\Work\\CEFF Metrics and Models\\Results\\Variable Selection\\Stream Classification\\WRR Update\\",curmet,""))


tmet<-tmet[which(tmet$Stat==curmet),]

cor.results <- matrix(NA,nrow=length(unique(tmet$ID)),8)



# Begin Modeling ----------------------------------------------------------

curmet.list <- list()

# Begin Modeling ----------------------------------------------------------
# selected a reduced set of predictors for modeling
new <- tmet[,c("ID","Year","Value", "DRAIN_SQKM",
                   
                   "ppt_Oct.wy", "ppt_Nov.wy", "ppt_Dec.wy",
                   "ppt_Oct.nwy", "ppt_Nov.nwy", "ppt_Dec.nwy", "ppt_Jan.wy","ppt_Feb.wy",
                   "ppt_Mar.wy", "ppt_Apr.wy", "ppt_May.wy", "ppt_Jun.wy", "ppt_Jul.wy","ppt_Aug.wy", "ppt_Sep.wy",
                   
                   "tav_Oct.wy" , "tav_Nov.wy" , "tav_Dec.wy",
                   "tav_Oct.nwy" , "tav_Nov.nwy" , "tav_Dec.nwy" , "tav_Jan.wy", "tav_Feb.wy","tav_Mar.wy" ,"tav_Apr.wy",
                   "tav_May.wy","tav_Jun.wy" ,"tav_Jul.wy" , "tav_Aug.wy","tav_Sep.wy",
                   
                   "tav_sum1", "tav_sum2" , "tav_sum3",  "tav_sum4", 
                   "tav_sprg.wy", "tav_fall.wy", "tav_wint.wy", 
                   "tav_ann.wy",
                   
                   "ppt_sprg.wy", "ppt_fall.wy", "ppt_wint.wy",
                   "ppt_sum1", "ppt_sum2", "ppt_sum3", "ppt_sum4",
                   "ppt_ann.wy",
                   
                   
                   "ANN_CONS_WETDAYS", 
                   
                   "ET", "PET","RH_BASIN", "PPTAVG_BASIN" ,"T_AVG_BASIN", "T_MIN_BASIN" , "T_MAX_BASIN","AnnMinPrecip","AnnMaxPrecip",
                   
                   # watershed characteristcs
                   "SLOPE_PCT_30M","ELEV_MEAN_M_BASIN_30M","ELEV_MIN_M_BASIN_30M","ELEV_MAX_M_BASIN_30M","ElvRng","ECO3",
                   "SNOW_PCT_PRECIP", "FST32F_BASIN","LST32F_BASIN",  "CONTACT", "BFI_AVE","TOPWET",             
                   "PERDUN","RFACT","PERHOR",  "DEPTH_WATTAB" ,"SILTAVE", "CLAYAVE" ,"SANDAVE",
                   
                   "KFACT_UP",    "NO10AVE",  "NO200AVE", "NO4AVE" , "OMAVE", "PERMAVE" , "ROCKDEPAVE" , "BDAVE", "AWCAVE","WTDEPAVE" ,          
                   "quarternary", "sedimentary" ,"volcanic", "gneiss", "granitic", "ultramafic",
                   
                   "KO_pct", "CaO_pct", "FeO_pct", "MgO_pct", "P_pct","S_pct","SiO_pct" , "UCS" ,"Lperm", "RECHARGE",  "pmpe",
                   "BDMAX",  "PERMH_WS",  "KRUG_RUNOFF",
                   "HGA" ,"HGB" , "HGC","HGD",
                   
                   "HLR1","HLR2", "HLR3", "HLR4","HLR5", "HLR6","HLR7", "HLR8", "HLR9","HLR10", "HLR11","HLR12", "HLR13", "HLR14","HLR15",
                   "HLR16","HLR17","HLR18","HLR19", "HLR20"
                   
)]




registerDoParallel(cores=4)
sites = unique(tmet$ID)

for(s in 1:length(sites)){
  site = sites[s]
  
  clist<-sites[-s]
  train<-new[new$ID %in% clist,]
  test<-new[new$ID %in% sites[s],]

  # x var needs to be altered for timing metrics (not scaled by DA)
  
 
  
    train$Value <-as.factor(train$Value)
    test$Value <-as.factor(test$Value)
 
  
  
  library(randomForest)
  
  
  rf.modcur <- randomForest(train$Value ~., train[,!(colnames(train) %in% c("ID","Year","Value","DRAIN_SQKM"))], ntrees=2300, mtry=10)
  
  ### Save variable importance
  varimp<-as.data.frame(rf.modcur$importance)
  varimp$Variable = row.names(varimp)
  
  write.csv(varimp, file = paste0("C:\\Work\\CEFF Metrics and Models\\Results\\Variable Selection\\Stream Classification\\WRR Update\\" ,curmet,"\\",unique(test$ID),".csv"),row.names=F)
  
  #if(curmet == "Streamclass_annual"){
    
    predValid <- predict(rf.modcur, test, type = "class")
#    Valid = mean(predValid == test$Value)
    #Valid = match(predValid,test$Value)
    
    ag_data <- data.frame(test$Value, predValid)
    
    table_agree <- confusionMatrix(as.factor(test$Value),as.factor(predValid),mode='everything')
    rf.cor <- table_agree$overall[[1]]

    x = 1:length(test$Value)
    
    R <- rf.cor
    print(R)
    
    FI <- table_agree$table[[2]]/nrow(test) # false negative
    FP <- table_agree$table[[3]]/nrow(test)
    TI <- table_agree$table[[1]]/nrow(test[which(test$Value == "Intermittent"),])
    TP <- table_agree$table[[4]]/nrow(test[which(test$Value == "Perennial"),])
    
    TrueClass <- nrow(test[which(test$Value=="Intermittent"),])/nrow(test)
    PredClass <- nrow(ag_data[which(ag_data$predValid=="Intermittent"),])/nrow(test)
    
    val<-data.frame(ID=test$ID,Year=test$Year,Obs=test$Value,AREA=test$DRAIN_SQKM,#WY=test$Year
                    p50 = predValid) # ,bayescor,bf.mod$yhat.test.mean)
    
    cor.results[s,1] <- paste0(site)
    cor.results[s,2] <-  R
    cor.results[s,3] <- FP
    cor.results[s,4] = FI
    cor.results[s,5] = TI
    cor.results[s,6] = TP
    cor.results[s,7] = TrueClass
    cor.results[s,8] = PredClass
  
  write.csv(val,paste0("C:\\Work\\CEFF Metrics and Models\\Results\\Model Performance\\Stream Classification\\WRR Update\\",curmet,"\\",unique(test$ID),"_valpreds.csv"),row.names=FALSE)
  
  rm(predValid,val)
}

c.test <- data.frame(apply(cor.results, 2, function(x) as.numeric(x)))
colnames(c.test) <- c("ID","Accuracy","False_Perennial","False_Nonperennial","True_Nonperennial","True_Perennial","Overall_Ob_Nonperennial","Overall_Pr_Nonperennial")

c.test$Ob_Class <- ifelse(c.test$Overall_Ob_Nonperennial>= 0.15, "Nonperennial","Perennial")
c.test$Pr_Class <- ifelse(c.test$Overall_Pr_Nonperennial>= 0.15, "Nonperennial","Perennial")
c.test$Match <- c.test$Pr_Class == c.test$Ob_Class

# Overall stats
wrong = c.test[which(c.test$Match == "FALSE"),]
mean(c.test$False_Nonperennial)
mean(c.test$False_Perennial)

# Perennial stats
perennial = c.test[which(c.test$Ob_Class == "Perennial"),]
mean(perennial$Accuracy)
mean(perennial$False_Nonperennial)
mean(perennial$False_Perennial)
pwrong = perennial[which(perennial$Match == "FALSE"),]

# Nonperennial stats
nonperennial = c.test[which(c.test$Ob_Class == "Nonperennial"),]
mean(nonperennial$Accuracy)
mean(nonperennial$False_Perennial)
mean(nonperennial$False_Nonperennial)
nwrong = nonperennial[which(nonperennial$Match == "FALSE"),]

write.csv(c.test,"C:\\Work\\CEFF Metrics and Models\\Results\\Model Performance\\Stream Classification\\WRR Update\\LOOCV_Classification_Accuracy_FP_FN.csv",row.names = F)

# Bind Everything together

library(data.table); library(tidyr)
setwd("C:\\Work\\CEFF Metrics and Models\\Results\\Model Performance\\Stream Classification\\WRR Update\\Streamclass_annual")
file_list <- list.files("C:\\Work\\CEFF Metrics and Models\\Results\\Model Performance\\Stream Classification\\WRR Update\\Streamclass_annual")


myMergedData <-   do.call(rbind, lapply( file_list, read.csv))
head(myMergedData)

myMergedData = myMergedData[which(myMergedData$ID %in% sites),]

#write.csv(myMergedData, "C:\\Work\\CEFF Metrics and Models\\Results\\Model Performance\\Stream Classification\\Overall_Altered_Annual_Obs_Prs_WRR_update.csv")
