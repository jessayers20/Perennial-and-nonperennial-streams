######################################################################
#####             Determine intermittent based                   #####
#####               on new thresholds                            #####
######################################################################

# check # of days on average below threshold
# check X% of the time
# Cluster values based on new thresholds -- 5 consecutive days during the dry sesaon
# uses low/zero flow metrics cacluated from 1_low_flow_updated

# compare to previous kmeans clusters

cat("\014") ;  rm(list=ls())

options(max.print=9999999)

library(VSURF); library(varSelRF); library(caret); library(Boruta); library(gamlss)
library(randomForest);library(dplyr);library(data.table);library(foreach)
library(ranger);library(doParallel); registerDoParallel(cores=4); library(ggplot2)
library(magrittr); library(raster)
library(ggpubr); library(factoextra)


# Select low flow metric --------------------------------------------------
pathfm = "C:/Work/Data/Shapefiles/California State"
usa <- readOGR(pathfm, layer="CA_State_TIGER2016")
usa <- spTransform(usa,CRS("+init=epsg:4326"))
ca = fortify(usa)


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
sites = unique(tmet$ID)

altered$Ref_status <- factor(altered$Ref_status, levels = c("Ref","Nonref"),labels=c("Minimally disturbed","Disturbed"))


ores <- read.csv("C:\\Work\\CEFF Metrics and Models\\Results\\Model Performance\\Stream Classification\\Overall_Altered_Annual_Obs_Prs_WRR_update.csv",header=T)
ores$ID = ores$site
all.alltered = left_join(ores,altered,by="ID")
all.alltered$Class = factor(all.alltered$Obs.Overall.Class, levels=c("Perennial","Intermittent"),labels=c("Perennial","Non-perennial"))

G_values =  
ggplot()+
  geom_polygon(data=ca,aes(x=long,y=lat,group=group),colour="grey10",fill="grey90",size=0.6)+
  geom_point(data=all.alltered, aes(x=LONGITUDE, y=LATITUDE,fill=Ref_status.y,shape=Class),size=3)+
  scale_shape_manual(values=c(21,22),name="Class")+
  scale_fill_manual(values=c("Minimally disturbed"="white","Disturbed"="grey50"), name="Contemporary Gages")+
  scale_color_manual(values=c("black","black"),name="Contemporary Gages")+
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),                                     
        strip.text = element_text(size=rel(0.2),color="white"),
        legend.position = c(0.75,0.8))


ggsave(paste0("Site_Location_Paper_altered_classification.png"), plot=G_values, device="png", path = "C:\\Work\\CEFF Metrics and Models\\Figures\\Altered Analysis",
       width = 6, height = 7, dpi=300)


altered$Ref_status = factor(altered$Ref_status, levels = c("Ref","Nonref"),
                            labels = c("Minimally disturbed","Disturbed"))

G_values =  
  ggplot()+
  geom_polygon(data=ca,aes(x=long,y=lat,group=group),colour="grey10",fill="grey90",size=0.6)+
  # geom_point(data=int, aes(x=dec_long_va, y=dec_lat_va, fill = Percent5),shape=21, color="black",size=4)+
  geom_point(data=altered, aes(x=LONGITUDE, y=LATITUDE,fill=Ref_status),shape=21, color="black",size=3)+
  #   scale_fill_manual(values = c( "#FFFFFF", "#FFEFD1" ,"#FFDFA4" ,"#FFCF77" ,"#FFBF4A", "#FFAF1D" ,"#FF8C00", "#FF4600", "#FF0000","#993404","black"),
  #                   name="Mean zero days (< 0.1cfs)")+
  scale_fill_manual(values=c("Minimally disturbed"="white","Disturbed"="grey50"), name="Modern Gage")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),                                     
        strip.text = element_text(size=rel(0.2),color="white"),
        legend.position = c(0.75,0.8))


ggsave(paste0("Site_Location_altered_paper.png"), plot=G_values, device="png", path = "C:\\Work\\CEFF Metrics and Models\\Figures\\Altered Analysis",
       width = 6, height = 7, dpi=300)


#data frame to save values
zavs <- data.frame(matrix(data=NA, nrow=length(sites),ncol=7))

for(s in 1:length(sites)){
  
  # site=10259000
  site = sites[s]
  temp1=read.csv(paste0("C:\\Work\\Data\\Low Flow Metrics Dry Season Consecutive\\",site,"_Annual_low_flow.csv"),header=T)
  temp1$site = paste0(site)
  
  #get ref years
  #keyvars <- unique(tmet[which(tmet$ID == site),]$Year)
  #temp1 = temp1[which(temp1$Year %in% keyvars),]
  
  # get years of overlapping data for models
  refyears <- unique(tmet[which(tmet$ID == site),]$Year)
  temp1 = temp1[which(temp1$Year %in% refyears),]
  
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
  
  
  rm(occur,occur5,percent5,percent,amean,amedian,temp1,site,s)
}

colnames(zavs) <- c("ID","Mean_Zero","Median_Zero","Percent","Percent5","N_years","N_zero_flow_years")
zavs$ID = as.integer(zavs$ID)


test = zavs[which(zavs$N_zero_flow_years >= 0),]
test = test[which(test$Percent5 < 0.15),]


allper <- left_join(zavs,all.eff,by = "ID")

per = allper[which(allper$Percent5 < 0.15),]
ptest <- per[which(per$Percent5 > 0),]
per$Category <- paste0("Perennial")
int = allper[which(allper$Percent5 >= 0.15),]
itest = int[which(int$Percent5 ==1),]
int$Category <- paste0("Non-perennial")


# Plot of USGS gages to EFF gages
classes_R = seq(0,1,0.1) # Only when I will have negative values I will be able to put -1 as initial value of cuts
n_cuts = length(classes_R)-1
block1 = rev(colorRampPalette(c("white","lightblue","darkblue"),bias=0.5)(n_cuts/2-1))
block3 = colorRampPalette(c("white","orange","red"),bias=0.5)(n_cuts)
block2 = colorRampPalette(c("white","white"))(2)
my.palette_R = c(block3)



allcats <- rbind(per,int)
#write.csv(allcats, "C:\\Work\\Data\\Functional Flows\\Streamflow_classification_5condays_dryseason_altered.csv",row.names=F)

int$category <- cut(int$Percent5, 
                    breaks=c(0, 0.1, 0.2, 0.3,0.4, 0.5,0.6,0.7,0.8,0.9,1), 
                    labels=c("10","20","30","40","50","60","70","80","90","100"))

int$category <- cut(int$Mean, 
                    breaks=c(0, 5, 15, 20, 25, 30, 35, 40, 45, 50, 100, 300), 
                    labels=c("0-5", "5-15", "15-20", "20-25"," 25-30", "30-35"," 35-40", "40-45", "45-50", "50-100"," 100-300"))


R_values =  
  ggplot()+
  geom_polygon(data=ca,aes(x=long,y=lat,group=group),colour="grey10",fill="grey90",size=0.6)+
  geom_point(data=allcats, aes(x=LONGITUDE, y=LATITUDE, fill = Percent5*100,shape=Category), color="black",size=3)+
  scale_shape_manual(values =c(22,21), name="")+
  scale_fill_gradientn(colors =my.palette_R, breaks=round(seq(0,100,10),2), limits=c(0,100),
                       guide = guide_colorbar(barwidth = 1, barheight = 6), oob = scales::squish, name="Percent of years that dry") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),                                     
       # strip.text = element_text(size=rel(0.1),color="white"),
        legend.position = c(0.8,0.77))


ggsave(paste0("Percent_zeroflow_Reference_streams.png"), plot=R_values, device="png", path = "C:\\Work\\CEFF Metrics and Models\\Figures\\Altered Analysis",
      width = 7, height = 8, dpi=300)


test <- allcats[which(allcats$Category == "Perennial"),]
test <- test[which(allcats$N_zero_flow_years > 0),]


#add also annual total flow?
set.seed(123)

res.km <- kmeans(scale(zavs[,c(4)]),2, nstart = 25)
# K-means clusters showing the group of each individuals
res.km$cluster

fviz_cluster(res.km, data = zavs[,c(3,4)],
             palette = c("#2E9FDF", "#00AFBB", "#E7B800","red","green","pink"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)

k.res <- cbind(zavs,res.km$cluster)
write.csv(k.res, "C:\\Work\\Data\\Functional Flows\\Low flow datasets\\SF_Class_Cluster_3_CY_Percent_Zeros_point1",row.names=F)

# For maps of percent zero flow over record -------------------------------

stationsfile$site_no = as.factor(stationsfile$site_no)

clust <- left_join(k.res,stationsfile,by = "site_no")

library(raster)
#Transform lat/long data for EFF data
pts = clust[,c("dec_lat_va","dec_long_va")]
coordinates(pts) = ~dec_long_va+dec_lat_va
projection(pts) <- CRS("+init=epsg:4326")
Long = clust$Long 
Lat = clust$Lat 


C_values =  
  ggplot()+
  geom_polygon(data=ca,aes(x=long,y=lat,group=group),colour="grey10",fill="grey90",size=0.6)+
  geom_point(data=clust,aes(x=dec_long_va, y=dec_lat_va,fill=as.factor(res.km$cluster)),shape=21,color="grey50",size=4)+
  scale_fill_manual(values = c("#fdc086", "#7fc97f", "#beaed4","#fdc086"),
                    name="Cluster (percent years with <0.1 cfs)")+
  #geom_point(data=pos, aes(x=dec_long_va, y=dec_lat_va,fill=category),shape=21, color="black",size=4)+
  #    geom_point(data=none, aes(x=dec_long_va, y=dec_lat_va),fill="white",shape=21, color="grey50",size=2)+
  #   scale_fill_manual(values = c( "#FFFFFF", "#FFEFD1" ,"#FFDFA4" ,"#FFCF77" ,"#FFBF4A", "#FFAF1D" ,"#FF8C00", "#FF4600", "#FF0000","#993404"),
  #                    name="% Zero days (=0 cfs)")+ 
  
  #scale_fill_gradientn(colors =my.palette_R, breaks=round(seq(0,1,0.2),3), limits=c(0,1),
  #                     guide = guide_colorbar(barwidth = 1, barheight = 10), oob = scales::squish, name="Zero Days (% entire record; <0.1cfs)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        strip.text = element_text(size=rel(0.2),color="white"),
        legend.position = c(0.8,0.8))

  ggsave(paste0("Cluster_percent_years_point1.png"), plot=C_values, device="png", path = "C:\\Work\\California Reference Gages\\Figures\\Zero Flow",
         width = 9, height = 9, dpi=300, limitsize=TRUE)
  
  
  library(usmap)
  library(ggplot2)
  
  all_states <- map_data("state")  
  cal <- filter(all_states, region == "california")
  
  
  p <- ggplot()
  p <- p + geom_polygon( data=all_states, aes(x=long, y=lat, group = group),colour="grey30", fill="white",size=1)+
    theme(panel.background = element_rect(color = "black", fill = "white"))
  
  p + geom_polygon(data = cal, aes(x=long, y=lat, group = group),fill="blue")+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          rect = element_blank())  
  
