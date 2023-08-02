library(ARTool) # ARTANOVA (nonparametric ANOVA - Ingo says, its a nonparametric GLMM)
library(vegan)
library(grid)
library(dplyr)
library(ggplot2)
library(plyr)
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
#===========================================================================================================
#--Read data-------------------------------------------------------------------------------------------
#===========================================================================================================
W<-read.csv("SoarGlide.csv")
# - remove captivity bred
W <-W[W$individual!=1891513,]
W <-W[W$individual!=170602,]

W$Vbgw<-W$tail_wind+14.4 # Theoretical best-gliding air speed for maximizing gliding range relative to the ground
W$soar_Gl_Afic<-(W$distance_gliding_m)/W$time_soaring_sec
W$glide_ratio <-W$distance_gliding_m/((W$sink_speed_m_sec)*-1)
  
W$datetime_start <- AsDateTime(W$datetime_start,us.format = FALSE)
W$Date <- AsDate(W$Date,us.format = FALSE)

W$MonthNum<-as.numeric(format(as.Date(W$datetime_start), "%m"))
W$Season <- 0
W$Season[W$MonthNum>1 & W$MonthNum<8] <- 1
W$Season <- factor(W$Season,labels=c("fall","spring"))

W$area <- 0
W$area[W$lat_start>=42 & W$lat_start<50] <- 1
W$area[W$lat_start>=32 & W$lat_start<42] <- 2
W$area[W$lat_start>=25 & W$lat_start<32] <- 3
W$area[W$lat_start<25] <- 4
W$area[W$sea_land==2] <- 6
W$area[W$sea_land==1] <- 5

W$area2 <- 0
W$area2[W$lat_start>32] <- 1
W$area2[W$lat_start<=32] <- 2
W$area2[W$lat_start>=30.85 & W$lat_start<=36.92] <- 3
W$area2[W$sea_land==2] <- 4
W$area2[W$sea_land==1] <-5
W$area2 <- factor(W$area2,labels=c("North","Desert","Med-SeaLand","MedSea","BKSea"))


W$areaGlob <- 0
W$areaGlob[W$lat_start<=32] <- 1
W$areaGlob[W$sea_land>0] <- 2
W$areaGlob <- factor(W$areaGlob,labels=c("North","Desert","Sea"))

W$SeaOrLand <- 0
W$SeaOrLand[W$area>4] <- 1
W$SeaOrLand <- factor(W$SeaOrLand,labels=c("Land","Sea"))

W$area <- factor(W$area,labels=c("Russia","Mid-Latitude","BlackSea-Israel","Sinai","Sahara","Black_sea","Med_sea"))

#- age
W$age_categ <- 1
W$age_categ[W$age_w=="Subadult"]<-2
W$age_categ[W$age_w=="Unknown adult" | W$age_w=="breeding adult"]<-3
W$age_categ<- factor(W$age_categ,labels = c("1-yo birds","2/3-yo birds","adult birds"))

W$age_categ2 <- 1
W$age_categ2[W$age_w=="Unknown adult" | W$age_w=="breeding adult"]<-2
W$age_categ2<- factor(W$age_categ2,labels = c("non-adult","adult birds"))


#===========================================================================================================
#-- colors ----------------------------------------------------
#===========================================================================================================

Desert <- '#e38519'
North <- '#5a7345'
Sea <- '#5262a3'


grey <- '#404145'
#===========================================================================================================
#-- Tailwind calculated vs. model----------------------------------------------------
#===========================================================================================================
W1 <- W %>%
  filter(!is.na(tail_wind))
W1 <- W1 %>%
  filter(!is.na(tail_wind_calc))

W1 %>% ggplot(aes(x= tail_wind,y=tail_wind_calc)) +
  geom_point(size = 2) + 
  geom_smooth(method = lm, se = FALSE)+
  #xlab("days since year start (fall)") + ylab("crossing distanse (km)")+
  theme_bw()

m=lm(tail_wind ~ tail_wind_calc, data=W1)
summary(m)
#==========================================================================================================
#-- clean and organaize----------------------------------------------------
#===========================================================================================================
#- (1) Inspect Missing Data
W %>%
  summarise_each(list(~sum(is.na(.)))) %>%
  gather()
#-- remove the missing
WRAFI <- W %>%
  filter(!is.na(RAFIcalc))

WRAFI <- WRAFI %>%
  filter(!is.na(falp_rate_climb))

WRAFI <- WRAFI %>%
  filter(!is.na(falp_rate_glide))

WRAFI <- WRAFI %>%
  filter(!is.na(exit_alt_above_terr ))

WRAFI %>%
  summarise_each(list(~sum(is.na(.)))) %>%
  gather()

# If flap ration >1 make it 1
WRAFI$falp_prop_glide[WRAFI$falp_prop_glide>1] <-1
WRAFI$falp_prop_climb[WRAFI$falp_prop_climb>1] <-1

ddply(WRAFI, c("area2"),summarise, N    = length(RAFI))

a <- ddply(WRAFI, c("Date","individual","areaGlob"),summarise, 
           MeanRAFI = mean(RAFIcalc),
           MeanSoarGlideAf = mean(soar_Gl_Afic),
           MeanFlapRateGlide = mean(falp_prop_glide),
           MeanFlapRateSoar = mean(falp_prop_climb),
           MeanClimbRate = mean(climb_rate_m_sec),
           MeanExitAlt = mean(exit_alt_above_terr),
           N = length(exit_alt_above_terr))

a<-a[a$N>1,]

ddply(WRAFI, c("areaGlob"),summarise, 
      M = mean(climb_rate_m_sec),
      SD = sd(climb_rate_m_sec))


ddply(a, c("areaGlob"),summarise, 
      Msoar = mean(MeanFlapRateSoar),
      SDsoar = sd(MeanFlapRateSoar),
      Mglide = mean(MeanFlapRateGlide),
      SDglide = sd(MeanFlapRateGlide))

 
mean(a$MeanFlapRateSoar)
sd(a$MeanFlapRateSoar)
mean(a$MeanFlapRateGlide)
sd(a$MeanFlapRateGlide)

#--Create the figure

RaFi<-ggplot(aes(areaGlob, MeanRAFI, color=areaGlob,fill=areaGlob),data=a) +
  geom_boxplot(outlier.shape = NA,alpha = 0.9, colour = "black",fatten = 5)+
  geom_point(shape = 21,alpha = 0.3, size=2, colour = "black",position=position_jitter(0.1))+
  scale_color_manual(values=c(North, Desert, Sea))+
  scale_fill_manual(values=c(North, Desert, Sea))+
  labs(x="Geographic area", y = "RAFI")+
  theme_bw()+
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=17),
        axis.title.x=element_blank(),axis.title.y=element_text(size=15)) 



SoarGlideAF<-ggplot(aes(areaGlob, MeanSoarGlideAf, color=areaGlob,fill=areaGlob),data=a) +
  geom_boxplot(outlier.shape = NA,alpha = 0.9, colour = "black",fatten = 5)+
  geom_point(shape = 21,alpha = 0.3, size=2, colour = "black",position=position_jitter(0.1))+
  scale_color_manual(values=c(North, Desert, Sea))+
  scale_fill_manual(values=c(North, Desert, Sea))+
  labs(x="Geographic area", y = "Soaring-gliding efficiency")+
  theme_bw()+
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=17),
        axis.title.x=element_blank(),axis.title.y=element_text(size=15)) 



FlapRateGlide<-ggplot(aes(areaGlob, MeanFlapRateGlide,color=areaGlob,fill=areaGlob),data=a) +
  geom_boxplot(outlier.shape = NA,alpha = 0.9, colour = "black",fatten = 5)+
  geom_point(shape = 21,alpha = 0.3, size=2, colour = "black",position=position_jitter(0.1))+
  scale_color_manual(values=c(North, Desert, Sea))+
  scale_fill_manual(values=c(North, Desert, Sea))+
  labs(x="Geographic area", y = "Flap. ratio during gliding")+
  theme_bw()+
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=17),
        axis.title.x=element_blank(),axis.title.y=element_text(size=15)) 


FlapRateSaor<-ggplot(aes(areaGlob, MeanFlapRateSoar, color=areaGlob,fill=areaGlob),data=a) +
  geom_boxplot(outlier.shape = NA,alpha = 0.9, colour = "black",fatten = 5)+
  geom_point(shape = 21,alpha = 0.3, size=2, colour = "black",position=position_jitter(0.1))+
  scale_color_manual(values=c(North, Desert, Sea))+
  scale_fill_manual(values=c(North, Desert, Sea))+
  labs(x="Geographic area", y = "Flap. ratio during climbing")+
  theme_bw()+
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=17),
        axis.title.x=element_blank(),axis.title.y=element_text(size=15)) 

ClimbRate<-ggplot(aes(areaGlob, MeanClimbRate,color=areaGlob,fill=areaGlob),data=a) +
  geom_boxplot(outlier.shape = NA,alpha = 0.9, colour = "black",fatten = 5)+
  geom_point(shape = 21,alpha = 0.3, size=2, colour = "black",position=position_jitter(0.1))+
  scale_color_manual(values=c(North, Desert, Sea))+
  scale_fill_manual(values=c(North, Desert, Sea))+
  labs(x="Geographic area", y = "Climb rate (m/s)")+
  theme_bw()+
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=17),
        axis.title.x=element_blank(),axis.title.y=element_text(size=15)) 

ExitAlt<-ggplot(aes(areaGlob, MeanExitAlt/1000,color=areaGlob,fill=areaGlob),data=a) +
  geom_boxplot(outlier.shape = NA,alpha = 0.9, colour = "black",fatten = 5)+
  geom_point(shape = 21,alpha = 0.3, size=2, colour = "black",position=position_jitter(0.1))+
  scale_color_manual(values=c(North, Desert, Sea))+
  scale_fill_manual(values=c(North, Desert, Sea))+
  labs(x="Geographic area", y = "Exit altit. above ground (km)")+
  theme_bw()+
  theme(axis.text.x = element_text(size=16),axis.text.y = element_text(size=17),
        axis.title.x=element_blank(),axis.title.y=element_text(size=15))  


ggarrange(RaFi,SoarGlideAF,FlapRateSaor,FlapRateGlide, ClimbRate,ExitAlt,
          heights = c(3,3,3,3,3,3),
          widths = c(3,3,3,3,3,3),
          ncol = 2, nrow = 3,
          align = "v",
          legend = 'none')


#############################################################################################
#-- Sampel size: 
  ddply(a, c("areaGlob"),summarise, 
        N = length(MeanRAFI)) 

#-- statistical analysis
library(ARTool) 
library(emmeans)

# (1) RAFI
m = art(MeanRAFI ~ areaGlob + (1|individual), data=a)
anova(m) # use the non-capital anova, because this is the one that works with the package

#               F Df Df.res Pr(>F)  
#1 areaGlob 8.4859  2 252.34 0.0002712 *** 

#-- poshoc without interaction  
contrast(emmeans(artlm(m, "areaGlob"), ~ areaGlob), method="pairwise",pbkrtest.limit = 4194)


# (2) Soaring gliding aficiancy
m = art(MeanSoarGlideAf ~ areaGlob + (1|individual), data=a)
anova(m) # use the non-capital anova, because this is the one that works with the package

#               F Df Df.res Pr(>F)  
#1 areaGlob 1.2714  2 272.99 0.2821 


# (3) Flapping proportion during glide
m = art(MeanFlapRateGlide ~ areaGlob + (1|individual), data=a)
anova(m) # use the non-capital anova, because this is the one that works with the package

#               F Df Df.res Pr(>F)  
#1 areaGlob 20.146  2 272.22 6.9273e-09 ***

#-- poshoc without interaction  
contrast(emmeans(artlm(m, "areaGlob"), ~ areaGlob), method="pairwise",pbkrtest.limit = 4194)


# (4) Flapping proportion during soar
m = art(MeanFlapRateSoar ~ areaGlob + (1|individual), data=a)
anova(m) # use the non-capital anova, because this is the one that works with the package

#               F Df Df.res Pr(>F)  
#1 areaGlob 78.878  2 294.29 < 2.22e-16 ***

#-- poshoc without interaction  
contrast(emmeans(artlm(m, "areaGlob"), ~ areaGlob), method="pairwise",pbkrtest.limit = 4194)


# (5) climb rate

ddply(a, c("areaGlob"),summarise, 
      M = mean(MeanClimbRate),
      med = median(MeanClimbRate),
      SD = sd(MeanRAFI)) 

m = art(MeanClimbRate ~ areaGlob + (1|individual), data=a)
anova(m) # use the non-capital anova, because this is the one that works with the package

#               F Df Df.res Pr(>F)  
#1 20.756  2 266.53 4.1939e-09 ***

#-- poshoc without interaction  
contrast(emmeans(artlm(m, "areaGlob"), ~ areaGlob), method="pairwise",pbkrtest.limit = 4194)


# (5) exit altitude
m = art(MeanExitAlt ~ areaGlob + (1|individual), data=a)
anova(m) # use the non-capital anova, because this is the one that works with the package

#               F Df Df.res Pr(>F)  
#1 areaGlob 17.609  2 234.53 7.4991e-08 ***

#-- poshoc without interaction  
contrast(emmeans(artlm(m, "areaGlob"), ~ areaGlob), method="pairwise",pbkrtest.limit = 4194)




