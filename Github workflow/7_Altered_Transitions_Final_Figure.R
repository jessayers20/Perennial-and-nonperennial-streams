#####################################################################
#####             Intermittent and Perennial                    #####
#####               Compare reference predictions               #####
#####                      to actual streamclass                #####
#####################################################################

# Model performance of annual scale intermittent and perennial stream classification

cat("\014") ;  rm(list=ls())

options(max.print=9999999)

library(randomForest); library(raster); library(rgdal); library(tidyr)
library(dplyr);library(data.table);library(foreach); library(stats); library(sf)
library(ranger);library(doParallel); registerDoParallel(cores=4); library(ggplot2)

#State Shapefile
pathfm = "C:/Work/Data/Shapefiles/California State"
usa <- readOGR(pathfm, layer="CA_State_TIGER2016")
usa <- spTransform(usa,CRS("+init=epsg:4326"))
ca = fortify(usa)


# results has altered and reference, need to separate out
ores <- read.csv("C:\\Work\\CEFF Metrics and Models\\Results\\Matches_WRR_Updated_Altered_Analysis.csv",header = T)

# original cutoff network (lost~65 sites)
altered <- read.csv("C:\\Work\\Data\\USGS gages\\California Discharge\\CA_altered.csv",header = T)
altered = altered[which(altered$ID %in% unique(ores$ID)),]


comids <- unique(ores$COMID)

tmet <- read.csv("C:\\Work\\Data\\Functional Flows\\Met files\\met_dryseason_point1_5condays_15percentofyears_reference.csv",header=T)
tmet = tmet[which(tmet$Stat == "Streamclass_annual"),]


histsites <- ores[which(ores$ID %in% unique(tmet$ID)),]

# Number of sites transitioned to nonreference
histr <- histsites[which(histsites$Ref_status == "Nonref"),]
# Number of sites that are still reference
hist_cur_ref <- histsites[which(histsites$Ref_status == "Ref"),]

# totoal reference sites in the analysis (+4 new reference sites)
cur_ref <- altered[which(altered$Ref_status == "Ref"),]

new <- ores
#new$predbias = new$Overall_sum_match/new$N_Years

new$Transition <- ifelse(new$Match == 1 & new$Act_Cat == "Intermittent","No Change Intermittent",
                         ifelse(new$Match == 1 & new$Act_Cat == "Perennial","No Change Perennial",
                                ifelse(new$Match == 0 & new$Act_Cat == "Intermittent","Perennial to Intermittent",
                                       "Intermittent to Perennial")))

new$Transition = factor(new$Transition, levels = c("No Change Perennial","No Change Intermittent","Intermittent to Perennial","Perennial to Intermittent"),
                        labels = c("No Change","No Change","Intermittent to Perennial","Perennial to Intermittent"))

new = new[which(new$Npost_1980 >= 40),]
library("RColorBrewer")

brewer.pal(10, "RdYlGn")
block1 <- brewer.pal(10, "RdYlGn")[1:5] #oranges
block3 <- brewer.pal(10, "RdYlGn")[6:10] #greens
block2 = colorRampPalette(c("white","white"))(1)
my.palette_R = c(block1,block2, block3)

new$predbias <- new$predbias*100

new$Ref_status.x = factor(new$Ref_status, levels = c("Ref","Nonref"), 
                          labels = c("Reference","Non-reference"))
#new = new[-which(new$site  %in% c("9429100","9429490")),]

new1 = new[which(new$Ref_status.x == "Reference"),]
ntest = new1[which(new1$Transition=="Perennial to Intermittent"),]

new2 = new[which(new$Ref_status.x == "Non-reference"),]
rtest = new2[which(new2$Transition=="Perennial to Intermittent"),]
Ttest = new2[which(new2$Transition=="Intermittent to Perennial"),]



Iplot<-
  ggplot(new, aes(x=LONGITUDE,y=LATITUDE))+
  geom_polygon(data=ca,aes(x=long,y=lat,group=group),colour="grey10",fill="grey90",size=0.6)+
  geom_point(aes(shape=Act_Cat,colour=Act_Cat,fill=predbias), size=3, color="grey50")+#fill=predbias),size=2) +
  facet_grid(Ref_status.x~Transition,switch = "y") +
  scale_shape_manual(values=c(22,21),drop=F,name=NULL) +
  scale_fill_gradientn(colors =my.palette_R, breaks=seq(-100,100,by=25), limits=c(-100,100),
                       guide = guide_colorbar(barwidth = 20, barheight = 1),
                       name="Intermittent (-)                                             Perennial (+)")+
  guides(fill = guide_colourbar(title.position = "bottom",
                                title.hjust = .5,
                                barwidth = 20, barheight = 1))+
  # guides(fill = guide_legend(override.aes = list(shape=21)))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        #axis.moduleticks = element_blank(),
        axis.title = element_blank(),
        strip.background = element_rect(fill="grey30", size=0.3, linetype="solid"), # Pr #308fd8  Ob #4598ba
        panel.border = element_rect(colour = "black", fill=NA),
        strip.text = element_text(size=rel(1),color="white"),
        legend.position = "bottom")

#ggsave("Altered_Gages_Transition_WRR_update.png", plot=Iplot, device="png", path = "C:\\Work\\CEFF Metrics and Models\\Figures\\Altered Analysis",
 #      width = 8, height = 7, dpi=800)




sites=unique(new$ID)
trend.results<-list()

# Calculate trends

for(s in 1:length(sites)){
  #s=6
  station = as.numeric(sites[s])
  
  print(station)
  
  
  ref = read.csv(paste0("C:\\Work\\Data\\Altered Low Flow Metrics Dry Season Consecutive\\",station,"_Annual_low_flow.csv"), header=T)
  
  library(modifiedmk); library(Kendall)
  
  #Subset for predictions or observations
  all = ref
  
  allT <- all[,c(2,4)] 
  
  runmo = colnames(allT)
  
  # if(nrow(allT)>25){
  
  pw<-apply(allT[,1:length(allT)],2,function(x){mkttest(as.vector(x))}) #because zero values in columns
  pw = as.data.frame(pw)#grabs tau from output results
  zval = sapply(pw,"[[",1)
  tau = sapply(pw,"[[",6)
  sl= sapply(pw,"[[",2)
  pval = sapply(pw,"[[",5)
  pos <- ifelse(sl >= 0, "Positive", "Negative")
  five <- ifelse(abs(pval) <= 0.05, "Significant", "Nonsignificant")
  res1 = data.frame(zval,pval,tau,sl,pos,five)
  res1$month = row.names(res1)
  res1$site = paste(station)
  res1$ref_status <- altered[which(altered$ID == station),]$Ref_status
  
  trend.results[[s]] <- res1
  
  write.csv(res1,file=paste0("C:/Work/CEFF Metrics and Models/Results/Altered MK Trends/MK_trends_",station,".csv"), row.names=F)
  #}
  
}

all.trends <- data.frame(rbindlist(trend.results))

all.trends$ID <-all.trends$site

trends <- merge(all.trends,new, by = c("ID"))

trends$Significance = ifelse(trends$five == "Significant" & trends$pos == "Positive", trends$Sig <- "Increasing", 
                             ifelse(trends$five == "Significant" & trends$pos == "Negative", trends$Sig <- "Decreasing", 
                                    "Nonsignificant"))

trends$ref_status = factor(trends$ref_status, levels = c("Ref","Nonref"), labels = c("Reference","Nonreference"))
trends$month = factor(trends$month, levels = c("Min_7_Day_Mov_Avg" , "Zero_Flow_Days"),
                      labels = c("Q7", "Zero flow days"))



Iclass <- trends[which(trends$Act_Cat=="Intermittent"),]
Iclass <- Iclass[which(Iclass$month == "Zero flow days"),]
Pclass <- trends[which(trends$Act_Cat == "Perennial"),]
Pclass <- Pclass[which(Pclass$month == "Q7"),]





bclass <- rbind(Iclass,Pclass)

tplot<-
  ggplot(bclass, aes(x=LONGITUDE,y=LATITUDE))+
  geom_polygon(data=ca,aes(x=long,y=lat,group=group),colour="grey10",fill="grey90",size=0.6)+
  geom_point(aes(shape=Significance,colour=Significance,size=Significance,fill=Significance)) +
  facet_grid(ref_status~month,switch = "y") +
  scale_shape_manual(values=c(24,25,21),drop=F) +
  scale_colour_manual(values=c("black","black","grey60"),drop=F) +
  scale_size_manual(values= c(3,3,0.5),drop=F) +
  scale_fill_manual(values=c("#cb181d","#2171b5","grey60"),drop=F) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        #axis.moduleticks = element_blank(),
        axis.title = element_blank(),
        legend.title=element_blank(),
        strip.background = element_rect(fill="grey30", size=0.3, linetype="solid"), # Pr #308fd8  Ob #4598ba
        panel.border = element_rect(colour = "black", fill=NA),
        strip.text = element_text(size=rel(1),color="white"),
        legend.position = "bottom")

ggsave("Altered_Gages_Trends_Min_7_Zero_Flow_RefNon_classes_WRR_update.png", plot=tplot, device="png", path = "C:\\Work\\CEFF Metrics and Models\\Figures\\Altered Analysis",
       width = 5, height = 6, dpi=800)
