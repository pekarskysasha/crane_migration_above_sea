library(ARTool) 
library(vegan)
library(grid)
library(plyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(lme4)
library(car)
library(tidyr)
library(effects)
library(ggmap)
library(ggsn)
library(viridis)
library(sf)    
library(spData)
library(spDataLarge)
library(tmap)
library(jtools)
library(performance)
library(flipTime)
library(emmeans) # for post hoc comparisons 
#==========================================================================
# load data: thermal presence in 10 min interval sections ----------------
#==========================================================================
TS<-read.csv("TheramlStats.csv")

# - remove captivity bred
TS <-TS[TS$Individual!=1891513,]
TS <-TS[TS$Individual!=170602,]

TS$Vbgw<-TS$tw+14.4 # Theoretical best-gliding air speed for maximizing gliding range relative to the ground
TS$date <- AsDateTime(TS$date,us.format = FALSE)
TS$dateOnly <- as.Date(format(TS$date, "%Y-%m-%d"))
TS$MonthNum<-as.numeric(format(as.Date(TS$date), "%m"))
TS$Season <- 0
TS$Season[TS$MonthNum>1 & TS$MonthNum<8] <- 1
TS$Season <- factor(TS$Season,labels=c("fall","spring"))
TS$area <- 0
TS$area[TS$lat_start>=42 & TS$lat_start<50] <- 1
TS$area[TS$lat_start>=32 & TS$lat_start<42] <- 2
TS$area[TS$lat_start>=25 & TS$lat_start<32] <- 3
TS$area[TS$lat_start<25] <- 4
TS$area[TS$OverSea==2] <- 6
TS$area[TS$OverSea==1] <- 5
TS$areaGlob <- 0
TS$areaGlob[TS$area>2 & TS$area<5] <- 1
TS$areaGlob[TS$area>4] <- 2
TS$areaGlob2 <- 0
TS$areaGlob2[TS$area>2 & TS$area<5] <- 1
TS$areaGlob2[TS$OverSea==1] <- 2
TS$areaGlob2[TS$OverSea==2] <- 3
TS$areaGlob2 <- factor(TS$areaGlob2,labels=c("North","Desert","Black Sea","Med. sea"))
TS$areaGlob <- factor(TS$areaGlob,labels=c("North","Desert","Sea"))
TS$area <- factor(TS$area,labels=c("Russia","Mid-Latitude","BlackSea-Israel","Sinai","Sahara","Black_sea","Med_sea"))
TS$thermalPresense<-0
TS$thermalPresense[TS$percent_time_in_thermals>0]<-1
TS$day_time <- as.factor(TS$day_time)
#- age
TS$age_categ <- 1
TS$age_categ[TS$age_w=="Subadult"]<-2
TS$age_categ[TS$age_w=="Unknown adult" | TS$age_w=="breeding adult"]<-3
TS$age_categ<- factor(TS$age_categ,labels = c("1-yo birds","2/3-yo birds","adult birds"))

TS$age_categ2 <- 1
TS$age_categ2[TS$age_w=="Unknown adult" | TS$age_w=="breeding adult"]<-2
TS$age_categ2<- factor(TS$age_categ2,labels = c("non-adult","adult birds"))
#- day-night
TS$day_time[TS$time_sice_sunset_h>0] <- 0
TS$day_time[TS$time_sice_sunset_h<=0 & TS$time_sice_sunrize_h>=0] <- 1

#- information
print(paste0("number of 10 min. sections with at least one thermal: ",sum(TS$thermalPresense)))
print(paste0("number of individuals: ",length(unique(TS$Individual))))

#==========================================================================
# [A] Porportion of thermal sections in different regions ----------------
#==========================================================================
#-- general proportion of flapping
ddply(TS, c("areaGlob"),summarise, 
             propThermal = mean(thermalPresense),
             propThermalFlap = 1-propThermal,
             SE=sd(thermalPresense) / sqrt(length(thermalPresense)),
             SD=sd(thermalPresense),
             t=qt(0.95/2 + 0.5, length(thermalPresense)-1),   # tend to 1.96 if sample size is big enough
             CI=t*SE,
             N=length(thermalPresense))


ddply(TS, c("areaGlob","thermalPresense"),summarise, 
      N=length(thermalPresense))


library(ggpattern)

a3 <- ddply(TS, c("areaGlob","Season"),summarise, 
           propThermal = mean(thermalPresense),
           SE=sd(thermalPresense) / sqrt(length(thermalPresense)),
           SD=sd(thermalPresense),
           t=qt(0.95/2 + 0.5, length(thermalPresense)-1),   # tend to 1.96 if sample size is big enough
           CI=t*SE,
           N=length(thermalPresense))

#______FIGURE 2B top_______________________________________________________________________

Desert <- '#e38519'
North <- '#5a7345'
Sea <- '#5262a3'

#  with pattern fill & Color both seas together 
prop<-ggplot(a3, aes(x=areaGlob, y=propThermal, fill=Season)) + 
  geom_col_pattern(
    aes(areaGlob, propThermal, pattern = Season, fill = areaGlob), 
    colour  = 'black',
    alpha  = 0.9,
    pattern_fill    = 'black',
    pattern_colour  = 'black',
    pattern_density= 0.1,
    pattern_spacing = 0.05,
    pattern_key_scale_factor = 0.1,
    position=position_dodge())+ 
  geom_errorbar(aes(ymin=propThermal-SE, ymax=propThermal+SE), width=.2,
                position=position_dodge(.9))+
  ylim(0, 1)+
  scale_fill_manual(values=c(Desert,Sea, North, Sea,Sea,Desert,North))+
  labs(y = "Proportion of sections with thermal soaring")+
  theme_classic()+
  scale_y_continuous(breaks=seq(0,1,0.5), limits=c(0, 1))+
  theme(axis.text.y = element_text(size=17),
        axis.title.x=element_blank(),axis.title.y=element_text(size=15),axis.text.x = element_blank()) 

#==========================================================================
#-- statistical analysis using binomial Generalized Linear Mixed Model (GLMM) 
#==========================================================================
# Full model
m1<-glmer(thermalPresense ~ areaGlob*Season  +
            (1|Individual),
          data=TS,
          family=binomial(link = "logit"))

summary(m1)
library(emmeans)

#====Likelihood Ratio Tests (LRTs)==============================

# Without areaGlob
m3 <- glmer(thermalPresense ~ Season + (1|Individual),
            data=TS, family=binomial(link = "logit"))

# Without Season
m4 <- glmer(thermalPresense ~ areaGlob  + (1|Individual),
            data=TS, family=binomial(link = "logit"))

# Model without the interaction
m2 <- glmer(thermalPresense ~ areaGlob + Season + (1|Individual),
            data=TS, family=binomial(link = "logit"))
# Perform LRTs
anova(m1, m3)  # Test significance of areaGlob
anova(m1, m4)  # Test significance of Season
anova(m1, m2)  # Test significance of interaction (as previously shown)
#====post hoc===================================================


#-(a) Calculate estimated marginal means for areas (North, Desert, Sea): 
#--- What is the average presence of thermal soaring in each area, 
#----------adjusting for the fact that different seasons might have different soaring conditions?
# Estimate marginal means
emm <- emmeans(m1, specs = ~ areaGlob*Season)
# Pairwise comparisons for Season within each areaGlob
pairs(emm, by = "areaGlob",adjust = "Tukey")


#-(b) estimated marginal means for each area, adjusted for season
#--- it gives a warning: Results may be misleading due to involvement in interactions
#---- this is because: the estimated marginal means for 'areaGlob' alone might be 
#-----------misleading because they don't capture the varying effects across different seasons
emm_areas <- emmeans(m1, specs = ~ areaGlob)
pairs(emm_areas)
#===== show the estimated means===============================================
#-(c) pairwise comparisons within each season separately 
#----- compare North, Desert, and Sea within Fall, and then separately compare North, Desert, and Sea within Spring
emm_season <- emmeans(m1, specs = pairwise ~ areaGlob | Season)
summary(emm_season, adjust = "Holm")
emm_season_df <- as.data.frame(emm_season)
emm_season_df$emmean <- exp(emm_season_df$emmean) / (1 + exp(emm_season_df$emmean))  # Back-transform from logit scale
emm_season_df$SE <- exp(emm_season_df$SE) / (1 + exp(emm_season_df$SE)) # Back-transform from logit scale

#-(d) pairwise comparisons within each area separately 
emm_area <- emmeans(m1, specs = pairwise ~ Season | areaGlob)
summary(emm_area, adjust = "Holm")
emm_area_df <- as.data.frame(emm_area)
emm_area_df$emmean <- exp(emm_area_df$emmean) / (1 + exp(emm_area_df$emmean))  # Back-transform from logit scale
emm_area_df$SE <- exp(emm_area_df$SE) / (1 + exp(emm_area_df$SE))  # Back-transform from logit scale

#==========================================================================
# [B] Porportion of thermal sections during day and night ----------------
#==========================================================================
#--(1) general info
#-- proportion of night sections
sum(TS$day_time==0)/nrow(TS)

#-- proportion of night sections Desert
sum(TS$day_time==0 & TS$areaGlob2=="Desert")/sum(TS$areaGlob2=="Desert")
#-- proportion of night sections North
sum(TS$day_time==0 & TS$areaGlob2=="North")/sum(TS$areaGlob2=="North")

#-- proportion of night thermal sections
sum(TS$day_time==0 & TS$thermalPresense==1)/sum(TS$day_time==0)
#-- proportion of day thermal sections
sum(TS$day_time==1 & TS$thermalPresense==1)/sum(TS$day_time==1)

#==========================================================================
# [C] time in thermal in sections (thermalPresense==1 only) in different regions ----------------
#==========================================================================
TS1<-TS[TS$day_time==1,]
TStrip <- aggregate(time_in_thermals_sec~areaGlob +dateOnly+ Individual,data=TS1,FUN=sum)
TStripNum <- aggregate(time_in_thermals_sec~areaGlob +dateOnly+ Individual,data=TS1,FUN=length)
TStrip$NumberOfThermalSection <- TStripNum$time_in_thermals_sec
# calculate proportion time in thermals for trips with at least one thermal
TStrip$PropThermalTime <- TStrip$time_in_thermals_sec/
  (TStripNum$time_in_thermals_sec*600)
TStripWithThermal <- TStrip[TStrip$time_in_thermals_sec>0,]

TStripWithThermal$Individual<-factor(TStripWithThermal$Individual)
#---proportion of time soaring by region in sections with soaring ---------------------------
Sumstat <- ddply(TStripWithThermal, c("areaGlob"),summarise, 
                 thermalPresence = median(PropThermalTime))

mean(Sumstat$thermalPresence[Sumstat$areaGlob!="Sea"])/Sumstat$thermalPresence[Sumstat$areaGlob=="Sea"]    


#______FIGURE 2B bottom_______________________________________________________________________

propsoar<-ggplot(aes(x=areaGlob, y=PropThermalTime, fill=areaGlob),data=TStripWithThermal) +
  geom_boxplot(outlier.shape = NA,alpha = 0.9, colour = "black",fatten = 5)+
  geom_point(shape = 21,alpha = 0.3, size=2, colour = "black",position=position_jitter(0.1))+
  scale_color_manual(values=c(North, Desert, Sea))+
  scale_fill_manual(values=c(North, Desert, Sea))+
  labs(x="Geographic area", y = "Proportion of time soaring")+
  theme_bw()+
  scale_y_continuous(breaks=seq(0,1,0.5), limits=c(0, 1))+
  theme(axis.text.y = element_text(size=17),
        axis.title.x=element_blank(),axis.title.y=element_text(size=15),axis.text.x = element_blank()) 
#_____________________________________________________________________________________
#==========================================================================
#-- statistical analysis using using ART-Anova -  allowed to distinguish significance levels in post-hoc
#==========================================================================
 
m1 = art(PropThermalTime ~ areaGlob  + (1|Individual), data=TStripWithThermal)
anova(m1)
CompM1=contrast(emmeans(artlm(m1, "areaGlob"), ~ areaGlob), method="pairwise", adjust = "tukey",pbkrtest.limit = 3412)

#==========================================================================
# [D] Air speed ------------------------------------------------------------
#==========================================================================
TS %>%
  summarise_each(list(~sum(is.na(.)))) %>%
  gather()

TS0 <-TS%>%
  filter(!is.na(va))
TS0$Individual<-factor(TS0$Individual)
TS0$thermalPresense <-factor(TS0$thermalPresense)

# remove outliers
Q <- quantile(TS0$va, probs=c(.25, .75), na.rm = FALSE)
iqr <- IQR(TS0$va)
up <-  Q[2]+1.5*iqr # Upper Range  
low<- Q[1]-1.5*iqr # Lower Range
low <-0 # in this case 0 is accepted
TS1<- TS0[TS0$va > low & TS0$va < up,]

library(ggstatsplot)


summary_speed<- ddply(TS1, c("thermalPresense"),summarise, 
      N    = length(va),
      Average_va= mean(va),
      sd_va=sd(va),
      Average_vg= mean(vg),
      sd_vg=sd(vg))


ggplot(TS1) +
  geom_histogram(aes(x=va),fill="red", alpha=0.3, position="identity",bins = 50)+
  geom_histogram(aes(x=vg),fill="blue", alpha=0.3, position="identity",bins = 50)

ggplot(TS1) +
  geom_histogram(aes(x=tw),fill="green", alpha=0.7, position="identity",bins = 50)

TS1Agregate<- ddply(TS1, c("Individual","UniqueSectionCounter","thermalPresense"),summarise, 
                    N    = length(va),
                    Average_va= mean(va))

m = art(Average_va ~ thermalPresense  + (1|Individual), data=TS1Agregate)
anova(m)

ggplot(aes(x=thermalPresense, y=Average_va),data=TS1Agregate) +
  geom_boxplot(outlier.shape = NA,alpha = 0.9, colour = "black",fatten = 5)+
  geom_point(aes(fill=Individual),shape = 21,alpha = 0.3, size=2,position=position_jitter(0.1))+
  theme_bw()+
  theme(axis.text.y = element_text(size=17),
        axis.title.x=element_blank(),axis.title.y=element_text(size=15)) 



ggplot(TS1, aes(x=va, y=vg, color=thermalPresense)) +
  geom_point(size=2)
  
#==========================================================================
# [C] Invirometal variables influence on probabilty of thermal cycling in sea----------
#==========================================================================
TS_Sea <- TS[TS$OverSea>0,]
ddply(TS_Sea, c("areaGlob2"),summarise, 
       N = length(Individual),
       Nthermal= sum(thermalPresense))
#==========================================================================
# [C.1] Woodcock graph
#==========================================================================
TS_Sea1 <-TS_Sea
TS_Sea1$thermalPresense <- as.factor(TS_Sea1$thermalPresense)
#-clean nans
TS_Sea1 <- TS_Sea1 %>%
  filter(!is.na(Mean_DeltaT))

ddply(TS_Sea, c("Season"),summarise, 
      N = length(Individual),
      Nthermal= sum(thermalPresense))


TS_Sea1$thermalPresense <- as.factor(TS_Sea1$thermalPresense)
ggplot(TS_Sea1, aes(x=Mean_DeltaT, y=wvel, shape=thermalPresense)) + 
  stat_density_2d(aes(alpha = ..level.., fill=thermalPresense), geom="polygon")+
  #scale_alpha_continuous(range = c(0,1))+
  geom_point(aes(color=thermalPresense),size=2)+
  scale_shape_manual(values=c(4,1))+
  scale_color_manual(values=c("black","red"))+
  scale_fill_manual(values=c("black","red"))+
  labs(x="dalta T", y = "wind speed")+
  theme_classic()+
  theme(text = element_text(size=16),
        axis.text.x = element_text(size=14)) 

# calculate thermal with delta T < 1
length(TS_Sea1$Mean_DeltaT[TS_Sea1$thermalPresense==1 & TS_Sea1$Mean_DeltaT<1])/
  length(TS_Sea1$Mean_DeltaT[TS_Sea1$thermalPresense==1])
hist(TS_Sea1$Mean_DeltaT[TS_Sea1$thermalPresense==1])

#==========================================================================
# [C.2] Binomial GLMM (Table S4)
#==========================================================================
#-(II) create a new database with normalized variables to run the model-----------------------------------
#- (1) Inspect Missing Data
TS_Sea %>%
  summarise_each(list(~sum(is.na(.)))) %>%
  gather()
#-- remove the missing
TS_Sea <- TS_Sea %>%
  filter(!is.na(Mean_DeltaT))

#- (2)-make animal ID facor --------------------------
TS_Sea$Individual<-as.factor(TS_Sea$Individual)

#--(3)- create a clean data frame --------------------------
Thermal_Sea <- as.data.frame(TS_Sea$time_in_thermals_sec)
Thermal_Sea$Individual <- TS_Sea$Individual
Thermal_Sea$age <- TS_Sea$age_categ2
Thermal_Sea$season <- TS_Sea$Season
Thermal_Sea$tail_wind <- TS_Sea$tw
Thermal_Sea$cross_wind <- TS_Sea$sw
Thermal_Sea$wind_speed <- TS_Sea$wvel
Thermal_Sea$BLH <- TS_Sea$Mean_blh
Thermal_Sea$sea_level_pressure <- TS_Sea$Mean_msl
Thermal_Sea$sea_level_pressure[Thermal_Sea$sea_level_pressure>2000]<-Thermal_Sea$sea_level_pressure[Thermal_Sea$sea_level_pressure>2000]/100
Thermal_Sea$delta_T <- TS_Sea$Mean_DeltaT
Thermal_Sea$thermalPresense <- TS_Sea$thermalPresense
Thermal_Sea$sea <- TS_Sea$areaGlob2
names(Thermal_Sea)[names(Thermal_Sea) == "TS_Sea$time_in_thermals_sec"] <- "time_in_thermals_sec"
Thermal_Sea$wind_direction <- TS_Sea$wang
#- (4) Z-Scores Standardizing --------------------------
#---see here: https://www.r-bloggers.com/2018/04/z-is-for-z-scores-and-standardizing/
Thermal_Sea[c(8:9)]<-lapply(Thermal_Sea[c(8:9)], function(x) {
  y<-scale(x, center=TRUE, scale=TRUE)
}
)

#--(III)-- check correlation between variables------------------------------
#== (1) METHOD 1 ====
mydata <-Thermal_Sea[,c(5:9)]
library(corrplot)
mydata.cor = cor(mydata, method = c("pearson"))
corrplot(mydata.cor)

palette = colorRampPalette(c("green", "white", "red")) (20)
heatmap(x = mydata.cor, col = palette, symm = TRUE)


library("Hmisc")
mydata.rcorr = rcorr(as.matrix(mydata))
mydata.rcorr

#== (2) METHOD 2 ====

# this model is run to see the correlation
XX<-glm(time_in_thermals_sec ~   wind_speed + BLH + delta_T + sea_level_pressure,
        data=Thermal_Sea)

vif(XX) # measure of multicolinarity
# varaibles that are above 3 are correlated, below 3 no multicolinirity
#-- season and SeaAirTDiffBLK are almost ---


#== (3) conclusions===
#- BLH and delta_T are correlated

#--(IV)-- run the model for sea crossing together -----------------------------------------------------

#----(1.2) Full model
m1<-glmer(thermalPresense ~ wind_speed + delta_T + sea_level_pressure + age +
            (1|Individual),
          data=Thermal_Sea,
          family=binomial(link = "logit"))

summary(m1)


# proportion of variance explained (Nakagawa and Schielzeth)
#-- marginal: variance explained by the fixed effects
#-- conditional: variance explained by both fixed and random effects 
Variance <- as.data.frame(r2_nakagawa(m1, by_group = FALSE))

#--- plot to see
delta_T <- effect_plot(m1, pred = delta_T, interval = TRUE, plot.points = TRUE)
windSpeed <- effect_plot(m1, pred = wind_speed, interval = TRUE, plot.points = TRUE)
#tail_windBK <- effect_plot(m1, pred = tail_wind, interval = TRUE, plot.points = TRUE)
#cross_windBK<- effect_plot(m1, pred = cross_wind, interval = TRUE, plot.points = TRUE)
sea_level_pressure <- effect_plot(m1, pred = sea_level_pressure, interval = TRUE, plot.points = TRUE)

ggarrange(delta_T,windSpeed,sea_level_pressure,
          heights = c(3, 3,3,3),
          widths = c(3,3,3,3),
          labels = c("(a)", "(b)", "(c)","(d)"),
          ncol = 2, nrow = 2,
          align = "v")

#-- plot to see fixed effect sizes
# see here https://strengejacke.github.io/sjPlot/articles/plot_model_estimates.html
library(sjPlot)
library(sjlabelled)
library(sjmisc)
plot_model(m1, show.values = TRUE, value.offset = .3,transform = NULL, vline.color = "grey",
           axis.labels = c("age","sea level pressure",
                           expression(paste(delta, "T")),"wind speed"))+
  theme_classic()+
  scale_color_manual(values=c("#858581","#080807"))+
  theme(text = element_text(size=20),
        axis.text.x = element_text(size=20)) 
#----(1.2) Null model
m0<-glmer(thermalPresense ~
              (1|Individual),
            data=Thermal_Sea,
            family=binomial(link = "logit"))

#----(1.3)-- compare full and null model---------------------------------------------------------------
anova(m1,m0,test='Chisq')




