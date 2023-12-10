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
met <- read.csv("C:\\Work\\Data\\Functional Flows\\Met files\\Met_altered_analysis_12_2023.csv",header=T)

# To create model, need to get reference sites from original MET file

refmet <- read.csv("C:\\Work\\Data\\Functional Flows\\Met files\\met_dryseason_point1_5condays_15percentofyears_reference.csv",header=T)
refs = refs[which(refs$Stat == "Streamclass_annual"),]
refsites=unique(refs$ID)
refs = unique(refs[,c("COMID","ID")])

# For altered analysis
altered <- read.csv("C:\\Work\\Data\\USGS gages\\California Discharge\\CA_altered.csv",header = T)
altsites = unique(altered$ID)
altered = altered[,c("COMID","ID")]


idcom <- rbind(altered,refs)

# Metrics to predict
drymets <- unique(met$Stat)

# loop through all functional flows here or just select which metrics using drymets
for(m in 1:length(drymets)){
#
m=1
curmet <- drymets[m]

dir.create(paste0("C:\\Work\\CEFF Metrics and Models\\Results\\Model Performance\\Stream Classification\\WRR Update\\" ,curmet,""))
dir.create(paste0("C:\\Work\\CEFF Metrics and Models\\Results\\Variable Selection\\Stream Classification\\WRR Update\\",curmet,""))


tmet<-refmet[which(refmet$Stat==curmet),]

cor.results <- matrix(NA,nrow=length(unique(tmet$ID)),7)



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
                   
                   #"run_ann.wy","run_fall.wy","run_wint.wy",  "run_sprg.wy" , "run_summ.wy",
                   "run_sum1" , "run_sum2",   "run_sum3", "run_sum4" ,     
                   
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

for(s in 1:length(sites)){
  site = sites[s]
  
  clist<-sites[-s]
  train<-new[new$ID %in% clist,]
  test<-met[met$ID %in% sites[s],]
  test<-test[which(test$Stat == curmet),]
  
  # x var needs to be altered for timing metrics (not scaled by DA)
  
  if(curmet %in% timing){
    train$Value <- train$Value
    test$Value <- test$Value
  } else{
    # add one bc often have zero flow vals
    train$Value <- (train$Value+1)/train$DRAIN_SQKM
    test$Value <- (test$Value+1)/test$DRAIN_SQKM
  }
  
  if(curmet=="Streamclass_annual"){
    train$Value <-as.factor(train$Value)
    test$Value <-as.factor(test$Value)
  }
  
  
  library(randomForest)
  
  
  rf.modcur <- randomForest(train$Value ~., train[,!(colnames(train) %in% c("ID","Year","Value","DRAIN_SQKM"))], ntrees=300)
  
  ### Save variable importance
  varimp<-as.data.frame(rf.modcur$importance)
  varimp$Variable = row.names(varimp)
  
  write.csv(varimp, file = paste0("C:\\Work\\CEFF Metrics and Models\\Results\\Variable Selection\\Stream Classification\\WRR Update\\" ,curmet,"\\",unique(test$ID),".csv"),row.names=F)
  
  #if(curmet == "Streamclass_annual"){
    
    predValid <- predict(rf.modcur, test, type = "class")
#    Valid = mean(predValid == test$Value)
    #Valid = match(predValid,test$Value)
    
    ag_data <- data.frame(test$Value, predValid)
    
    table_agree <- confusionMatrix(test$Value,predValid,mode='everything')
    rf.cor <- table_agree$overall[[1]]

    x = 1:length(test$Value)
    
    R <- rf.cor
    print(R)
    
    FP <- table_agree$table[[2]]/nrow(test) # false negative
    FI <- table_agree$table[[3]]/nrow(test)
    TI <- table_agree$table[[1]]/nrow(test[which(test$Value == "Intermittent"),])
    TP <- table_agree$table[[4]]/nrow(test[which(test$Value == "Perennial"),])
    
    val<-data.frame(ID=test$ID,Year=test$Year,Obs=test$Value,AREA=test$DRAIN_SQKM,#WY=test$Year
                    p50 = predValid, cor = R, F_perennial = FP, F_intermittent=FI,T_intermittent=TI, T_perennial=TP) # ,bayescor,bf.mod$yhat.test.mean)
    
  #}
  #else{
  ### Predict to validation sites and save median, 10th, 25th, 75th, & 90th percentiles
  #preds <- predict(rf.modcur, test, type="response",predict.all=T)
  #predp50<-apply(preds$individual,1,median)
  #predp10<-apply(preds$individual,1,function(x) quantile(x,probs=0.1))
  #predp25<-apply(preds$individual,1,function(x) quantile(x,probs=0.25))
  #predp75<-apply(preds$individual,1,function(x) quantile(x,probs=0.75))
  #predp90<-apply(preds$individual,1,function(x) quantile(x,probs=0.9))
  #rf.cor <- cor(predp50,test$Value)^2
  #}
  
  
  #val<-data.frame(ID=test$ID,WY=test$Year,Obs=test$Value,AREA=test$DRAIN_SQKM,
   #               p50=predp50,p10=predp10,p25=predp25,p75=predp75,p90=predp90,
  #                rf.cor ) # ,bayescor,bf.mod$yhat.test.mean)
  
  
  write.csv(val,paste0("C:\\Work\\CEFF Metrics and Models\\Results\\Model Performance\\Stream Classification\\WRR Update\\",curmet,"\\",unique(test$ID),"_valpreds.csv"),row.names=FALSE)
  
  rm(predValid,val)
}


# Bind Everything together

library(data.table); library(tidyr)
setwd("C:\\Work\\CEFF Metrics and Models\\Results\\Model Performance\\Stream Classification\\WRR Update\\Streamclass_annual")
file_list <- list.files("C:\\Work\\CEFF Metrics and Models\\Results\\Model Performance\\Stream Classification\\WRR Update\\Streamclass_annual")


myMergedData <-   do.call(rbind, lapply( file_list, read.csv))
head(myMergedData)

#write.csv(myMergedData, "C:\\Work\\CEFF Metrics and Models\\Results\\Model Performance\\Stream Classification\\Overall_Altered_Annual_Obs_Prs_WRR_update.csv")
