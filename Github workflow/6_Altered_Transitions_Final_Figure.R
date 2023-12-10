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

# original cutoff network (lost~65 sites)
altered <- read.csv("C:\\Work\\Data\\USGS gages\\California Discharge\\CA_altered.csv",header = T)
altered = altered[which(altered$ID %in% unique(ores$ID)),]

ares <- merge(ores,altered, by="ID")

comids <- unique(ares$COMID)

tmet <- read.csv("C:\\Work\\Data\\Functional Flows\\Met files\\met_dryseason_point1_5condays_15percentofyears_reference.csv",header=T)
tmet = tmet[which(tmet$Stat == "Streamclass_annual"),]


histsites <- ares[which(ares$ID %in% unique(tmet$ID)),]
histr <- histsites[which(histsites$Ref_status == "Nonref"),]


ares$Match = ifelse(ares$Obs == ares$p50, 1, 0)
test = new[which(new$site %in% tmet$ID),]

new$predbias = new$Overall_sum_match/new$N_Years

new$Transition <- ifelse(new$Match == 1 & new$Obs.Overall.Class == "Intermittent","No Change Intermittent",
                         ifelse(new$Match == 1 & new$Obs.Overall.Class == "Perennial","No Change Perennial",
                                ifelse(new$Match == 0 & new$Obs.Overall.Class == "Intermittent","Perennial to Intermittent",
                                       "Intermittent to Perennial")))

new$Transition = factor(new$Transition, levels = c("No Change Perennial","No Change Intermittent","Intermittent to Perennial","Perennial to Intermittent"),
                        labels = c("No Change","No Change","Intermittent to Perennial","Perennial to Intermittent"))


library("RColorBrewer")

brewer.pal(10, "RdYlGn")
block1 <- brewer.pal(10, "RdYlGn")[1:5] #oranges
block3 <- brewer.pal(10, "RdYlGn")[6:10] #greens
block2 = colorRampPalette(c("white","white"))(1)
my.palette_R = c(block1,block2, block3)

new$predbias <- new$predbias*100

new$Ref_status.x = factor(new$Ref_status.x, levels = c("Ref","Nonref"), 
                          labels = c("Reference","Non-reference"))
#new = new[-which(new$site  %in% c("9429100","9429490")),]

new1 = new[which(new$Ref_status.x == "Reference"),]
ntest = new1[which(new1$Transition=="Perennial to Intermittent"),]

new2 = new[which(new$Ref_status.x == "Non-reference"),]

master <- read.csv("C:\\Work\\Data\\USGS gages\\usgs_master_list.csv",header=T)
master$site = master$usgs_id
newmast <- left_join(new,master)
newmast = newmast[!is.na(newmast$drain_area),]
newmast$drain_area = newmast$drain_area*2.58999
#new = new[-which(new$site %in% c("9429190",))]

Iplot<-
  ggplot(new, aes(x=LONGITUDE,y=LATITUDE))+
  geom_polygon(data=ca,aes(x=long,y=lat,group=group),colour="grey10",fill="grey90",size=0.6)+
  geom_point(aes(shape=Obs.Overall.Class,colour=Obs.Overall.Class,fill=predbias), size=2)+#fill=predbias),size=2) +
  facet_grid(Ref_status.x~Transition,switch = "y") +
  scale_shape_manual(values=c(24,22),drop=F,name=NULL) +
  scale_colour_manual(values=c("black","black","grey60"),drop=F,name=NULL) +
  # scale_fill_manual(values=c("white","purple","orange"),drop=F) +
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

ggsave("Altered_Gages_Transition.png", plot=Iplot, device="png", path = "C:\\Work\\CEFF Metrics and Models\\Figures\\Altered Analysis",
       width = 8, height = 7, dpi=300)