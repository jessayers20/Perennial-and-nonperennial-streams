#####################################################################
#####             Intermittent and Perennial                    #####
#####               Compare reference predictions               #####
#####                      to actual streamclass                #####
#####################################################################

# Model performance of annual scale intermittent and perennial stream classification

cat("\014") ;  rm(list=ls())

options(max.print=9999999)

library(randomForest); library(raster); library(rgdal); library(tidyr)
library(dplyr);library(data.table);library(foreach); library(stats)
library(ranger);library(doParallel); registerDoParallel(cores=4); library(ggplot2)

#State Shapefile
pathfm = "C:/Work/Data/Shapefiles/California State"
usa <- readOGR(pathfm, layer="CA_State_TIGER2016")
usa <- spTransform(usa,CRS("+init=epsg:4326"))
ca = fortify(usa)


# results has altered and reference, need to separate out
ores <- read.csv("C:\\Work\\CEFF Metrics and Models\\Results\\Model Performance\\Stream Classification\\Overall_Altered_Annual_Obs_Prs_WRR_update.csv",header=T)
ores=ores[,-1]
ores <- ores[which(ores$Year >= 1980),]
sites = unique(ores$ID)

# original cutoff network (lost~65 sites)
altered <- read.csv("C:\\Work\\Data\\USGS gages\\California Discharge\\CA_altered.csv",header = T)
altered = altered[which(altered$ID %in% unique(ores$ID)),]

ares <- merge(ores,altered, by="ID")
sites = unique(ares$ID)


comids <- unique(ares$COMID)

# file for predictors
r.results <- matrix(NA,nrow=length(sites),10)
var.list <- list()

sites <- unique(ares$ID)

for(s in 1:length(sites)){
  #dat1 <- read.csv(paste0("C:\\Work\\CEFF Metrics and Models\\Results\\Variable Selection\\model performance vars by site\\" ,curmet,"\\",site,".csv"),header=T)
  site = sites[s]
  
  dat2 <- ares[which(ares$ID == site),]
  dat2 = dat2[which(dat2$Year >= 1980),]
  dat2 = unique(dat2)
  
  dat2$Pre_class <- ifelse(dat2$p50 == "Perennial",1,0)
  dat2$Ob_class <- ifelse(dat2$Obs == "Perennial",1,0)
  
  rfcor = unique(dat2$cor)
  
  var.list[[s]] <- dat2
  
  r.results[s,1] <- site
  r.results[s,2] <- "Streamclass"
  r.results[s,3] <- as.numeric(rfcor)
  
  i_pred <- nrow(dat2[which(dat2$p50 == "Intermittent"),])/nrow(dat2)
  clu <- nrow(dat2[which(dat2$p50 == "Perennial"),])/nrow(dat2)
  
  if(i_pred >= 0.15){
    r.results[s,4] <- paste0("Intermittent")
  } else{r.results[s,4] <- paste0("Perennial")}
  
  act_cat <-   read.csv(paste0("C:\\Work\\Data\\Altered Low Flow Metrics Dry Season Consecutive\\",site,"_Annual_low_flow.csv"))
  act_cat = act_cat[which(act_cat$Year >= 1980),]
  
  
  i_act_cat <- nrow(act_cat[which(act_cat$Zero_Flow_Days >= 5),])/nrow(act_cat)
  
  if(i_act_cat >= 0.15){
    r.results[s,5] <- paste0("Intermittent")
  }else{r.results[s,5] <- paste0("Perennial")}
  
  r.results[s,6]<-ifelse(match(r.results[s,5],r.results[s,4])==1,1,0)
  
  if(nrow(dat2>0)){
    r.results[s,7] <- as.numeric((sum(dat2$Ob_class) - sum(dat2$Pre_class))/nrow(dat2))
    r.results[s,8] <- nrow(dat2)
    
  }
  
  r.results[s,9] <- mean(act_cat$Zero_Flow_Days)
  r.results[s,10] <- median(act_cat$Zero_Flow_Days)
  # var.list[[s]] <- read.csv(paste0("C:\\Work\\CEFF Metrics and Models\\Results\\Variable Selection\\Stream Classification\\",curmet,"\\",site,".csv"),header=T)
  
}

all.results <- data.frame(r.results)
colnames(all.results) <- c("ID","Metric","R","Pre_Cat","Act_Cat","Match","predbias","n_obs","Mean_Zero","Median_Zero")

all.results[is.na(all.results$Match),]$Match <- 0

all.results[,c(3,6:10)] <- apply(all.results[,c(3,6:10)],2,as.numeric)

all.results %>% group_by(Act_Cat) %>%
  summarize(mean_zero = mean(Mean_Zero), med_zero = median(Median_Zero),
            min_z = min(Mean_Zero), max_z=max(Mean_Zero))

all = merge(all.results,altered, by =c("ID"))

write.csv(all, "C:\\Work\\CEFF Metrics and Models\\Results\\Matches_WRR_Updated_Altered_Analysis.csv",row.names = F)

    