library(ggplot2)
library(dplyr)
library(plyr)
library(ggpubr)
library(rstatix) # for outlier calculation
library(lubridate) # extract year from date
library(DescTools) # for Dunnett's Test
## mediterranean sea cross ======================================================
ThermalFroBox<-read.csv("AnnotatedTimePointsMedSea.csv")
ThermalFroBox$Date <- as.Date(ThermalFroBox$Date, "%Y-%m-%d")
#--- change to No-thermals events that were highlited manualy as mistake
ToEraze<-read.csv("ToEraze.csv")
TS<-read.csv("TheramlStats.csv")

# - remove captivity bred
TS <-TS[TS$Individual!=1891513,]
TS <-TS[TS$Individual!=170602,]

TS$date <- AsDateTime(TS$date,us.format = FALSE)
TS$dateOnly <- as.Date(format(TS$date, "%Y-%m-%d"))
ToEraze$TimeStart <- AsDateTime(ToEraze$TimeStart,us.format = FALSE)
ToEraze$dateOnly <- as.Date(format(ToEraze$TimeStart, "%Y-%m-%d"))

for (z in 1:nrow(ToEraze)) {
  IND<-TS$Individual==ToEraze$Animal_ID[z] & TS$dateOnly==ToEraze$dateOnl[z]
  temp<-TS[IND,]
    if (sum(temp$percent_time_in_thermals>0 & temp$OverSea>0)==0) { # change to No-thermals
      IND_ThermalFroBox <- ThermalFroBox$indev==ToEraze$Animal_ID[z] & ThermalFroBox$Date==ToEraze$dateOnly[z]
      ThermalFroBox$thermap_presence[IND_ThermalFroBox]=0
      }
}
#----------------------------------------------------------
ThermalFroBox$thermap_presence<-factor(ThermalFroBox$thermap_presence,labels=c("No-thermals","thermals"))
ThermalFroBox$Fall_Spring <- factor(ThermalFroBox$fall,labels=c("fall"))
ThermalFroBox$time_point <- factor(ThermalFroBox$time_point,labels=c("starting point","sea enter","in sea"))
ThermalFroBox <- ThermalFroBox[ThermalFroBox$time_point!="sea enter",]
  
  
  
a <- ddply(ThermalFroBox, c("Date","thermap_presence","time_point"),summarise, N    = length(sstdiff))

temp=ThermalFroBox[ThermalFroBox$daysBack==0 & ThermalFroBox$time_point=="in sea",]
ddply(temp, c("thermap_presence","Date"),summarise, N    = length(sstdiff))
ddply(temp, c("thermap_presence"),summarise, N    = length(sstdiff))

ThermalFroBox <- ThermalFroBox[ThermalFroBox$time_point=="in sea",]
ThermalFroBox %>%
  summarise_each(list(~sum(is.na(.)))) %>%
  gather()

ThermalFroBox <- ThermalFroBox %>%
  filter(!is.na(sstdiff))
#== sample size:

# No-thermals  15; over 9 different migration dates
# thermals    25 ; over 17 different migration dates


# precipitation

Means <- ddply(ThermalFroBox, c("daysBack","thermap_presence"),summarise,
               mean= mean(GPCC),
               sd=sd(GPCC),
               N=length(GPCC),
               se=sd(GPCC)/length(GPCC))

precip<- ggplot(Means, aes(daysBack, mean)) +
  geom_line(aes(group = as.factor(thermap_presence),color = as.factor(thermap_presence)),size = 1)+
  scale_color_manual(values=c("#c70428","#087d8a"))+
  geom_point()+
  geom_errorbar(
    aes(ymin = mean-se, ymax = mean+se, group = as.factor(thermap_presence),color = as.factor(thermap_presence)),
    width = 0.2,size = 1)+
  theme(text = element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey50"))+
  ylab("Daily precipitation (mm/d)")+
  xlab("Number of timescales relative to departure")
  




# air temperature 
Means <- ddply(ThermalFroBox, c("daysBack","thermap_presence"),summarise,
               mean= mean(t2m),
               sd=sd(t2m),
               N=length(t2m),
               se=sd(t2m)/length(t2m))

 t_air<- ggplot(Means, aes(daysBack, mean)) +
   geom_line(aes(group = as.factor(thermap_presence),color = as.factor(thermap_presence)),size = 1)+
   scale_color_manual(values=c("#c70428","#087d8a"))+
   geom_point()+
   geom_errorbar(
     aes(ymin = mean-se, ymax = mean+se, group = as.factor(thermap_presence),color = as.factor(thermap_presence)),
     width = 0.2,size = 1)+
  theme(text = element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey50"))+
  ylab("Air temperature")+
   xlab("Number of timescales relative to departure")
 

# sea-air delta temprature
Means <- ddply(ThermalFroBox, c("daysBack","thermap_presence"),summarise,
               mean= mean(sstdiff),
               sd=sd(sstdiff),
               N=length(sstdiff),
               se=sd(sstdiff)/length(sstdiff))
             
delta_t <- ggplot(Means, aes(daysBack, mean)) +
  geom_line(aes(group = as.factor(thermap_presence),color = as.factor(thermap_presence)),size = 1)+
  scale_color_manual(values=c("#c70428","#087d8a"))+
  geom_point()+
  geom_errorbar(
    aes(ymin = mean-se, ymax = mean+se, group = as.factor(thermap_presence),color = as.factor(thermap_presence)),
    width = 0.2,size = 1)+
  theme(text = element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey50"))+
  ylab("temperture delta (sea-air)")+
  xlab("Number of timescales relative to departure")




# sea level pressure
Means <- ddply(ThermalFroBox, c("daysBack","thermap_presence"),summarise,
               mean= mean(msl),
               sd=sd(msl),
               N=length(msl),
               se=sd(msl)/length(msl))



msl <- ggplot(Means, aes(daysBack, mean)) +
  geom_line(aes(group = as.factor(thermap_presence),color = as.factor(thermap_presence)),size = 1)+
  scale_color_manual(values=c("#c70428","#087d8a"))+
  geom_point()+
  geom_errorbar(
    aes(ymin = mean-se, ymax = mean+se, group = as.factor(thermap_presence),color = as.factor(thermap_presence)),
    width = 0.2,size = 1)+
  theme(text = element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey50"))+
  ylab("Sea level pressure (mb)")+
  xlab("Number of timescales relative to departure")


# Total cloud cover
Means <- ddply(ThermalFroBox, c("daysBack","thermap_presence"),summarise,
               mean= mean(tcc),
               sd=sd(tcc),
               N=length(tcc),
               se=sd(tcc)/length(tcc))

cloadcove <-  ggplot(Means, aes(daysBack, mean)) +
  geom_line(aes(group = as.factor(thermap_presence),color = as.factor(thermap_presence)),size = 1)+
  scale_color_manual(values=c("#c70428","#087d8a"))+
  geom_point()+
  geom_errorbar(
    aes(ymin = mean-se, ymax = mean+se, group = as.factor(thermap_presence),color = as.factor(thermap_presence)),
    width = 0.2,size = 1)+
  theme(text = element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey50"))+
  ylab("Total cloud cover")+
  xlab("Number of timescales relative to departure")


# Tail wind
Means <- ddply(ThermalFroBox, c("daysBack","thermap_presence"),summarise,
               mean= mean(tw),
               sd=sd(tw),
               N=length(tw),
               se=sd(tw)/length(tw))

TailWind <-  ggplot(Means, aes(daysBack, mean),group = as.factor(thermap_presence)) +
  geom_line(aes(group = as.factor(thermap_presence),color = as.factor(thermap_presence)),size = 1)+
  scale_color_manual(values=c("#c70428","#087d8a"))+
  geom_point()+
  geom_errorbar(
    aes(ymin = mean-se, ymax = mean+se, group = as.factor(thermap_presence),color = as.factor(thermap_presence)),
    width = 0.2,size = 1)+
  theme(text = element_text(size=14),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "grey50"))+
  ylab("Tail wind (m/sec)")+
  xlab("Number of timescales relative to departure")



#== plot all together (Figure 4A)



ggarrange(delta_t,msl,cloadcove,precip,TailWind,
          heights = c(3,3,3),
          widths = c(3,3,3),
          ncol = 1, nrow = 5,
          align = "v")


##########################################################################################################
#--------------------------------------------------------------------------------------------------------
# STATISTICS
#--------------------------------------------------------------------------------------------------------
##########################################################################################################

# create the database
InSeaOnly=ThermalFroBox[ThermalFroBox$time_point=="in sea",c(1,2,4,5,26,19,18,31,33)]
InSeaOnly$daysBack1<-factor(InSeaOnly$daysBack,labels=c("t1","t2","t3","t4","t5"))

InSeaOnly$id[InSeaOnly$thermap_presence=="No-thermals"& InSeaOnly$daysBack==0]<-1:15
InSeaOnly$id[InSeaOnly$thermap_presence=="No-thermals"& InSeaOnly$daysBack==1]<-1:15
InSeaOnly$id[InSeaOnly$thermap_presence=="No-thermals"& InSeaOnly$daysBack==-1]<-1:15
InSeaOnly$id[InSeaOnly$thermap_presence=="No-thermals"& InSeaOnly$daysBack==-2]<-1:15
InSeaOnly$id[InSeaOnly$thermap_presence=="No-thermals"& InSeaOnly$daysBack==-3]<-1:15
InSeaOnly$id[InSeaOnly$thermap_presence=="thermals"& InSeaOnly$daysBack==0]<-1:25
InSeaOnly$id[InSeaOnly$thermap_presence=="thermals"& InSeaOnly$daysBack==1]<-1:25
InSeaOnly$id[InSeaOnly$thermap_presence=="thermals"& InSeaOnly$daysBack==-1]<-1:25
InSeaOnly$id[InSeaOnly$thermap_presence=="thermals"& InSeaOnly$daysBack==-2]<-1:25
InSeaOnly$id[InSeaOnly$thermap_presence=="thermals"& InSeaOnly$daysBack==-3]<-1:25
InSeaOnly$daysBack <-factor(InSeaOnly$daysBack)

##################
#--(a)-delta T
##################
#-(1)--Create box plots of the score colored by treatment groups:
ggboxplot(
  InSeaOnly, x = "daysBack", y = "sstdiff",
  color = "thermap_presence", palette = "jco"
)

#-(2)--Chek assumptions
#---(2.1)-Outliers
Outlier <- InSeaOnly %>%
  group_by(thermap_presence, daysBack) %>%
  identify_outliers(sstdiff)
#3 extrime outliers

#---(2.2)-Normality assumption
InSeaOnly %>%
  group_by(thermap_presence, daysBack) %>%
  shapiro_test(sstdiff) 
# thermals in (-3) is not normal
#---(2.3)-Create QQ plot for each cell of design
ggqqplot(InSeaOnly, "sstdiff", ggtheme = theme_bw()) +
  facet_grid(daysBack ~ thermap_presence, labeller = "label_both")

#-(3)--Computation


res.aov <- anova_test(data = InSeaOnly, dv = sstdiff, wid = id,
                      within = c(thermap_presence, daysBack1))

#-(3.1)--Effect of thermals on each time point (which time points are different between groups)
one.way <- InSeaOnly %>%
  group_by(daysBack) %>%
  anova_test(dv = sstdiff, wid = id, within = thermap_presence) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way


#-(3.2)--Effect of time at each  group (on which groups there is time effect)
one.way2 <- InSeaOnly %>%
  group_by(thermap_presence) %>%
  anova_test(dv = sstdiff, wid = id, within = daysBack) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2

#-(3.3)--make Dunnett Test onlt on the significant group
bbThermals<-InSeaOnly[InSeaOnly$thermap_presence=="thermals",]
DunnettTest(sstdiff~daysBack, data=bbThermals,control="0")


##################
#--(b)-sea level pressure
##################
#-(1)--Create box plots of the score colored by treatment groups:
ggboxplot(
  InSeaOnly, x = "daysBack", y = "msl",
  color = "thermap_presence", palette = "jco"
)

#-(2)--Chek assumptions
#---(2.1)-Outliers
Outlier <- InSeaOnly %>%
  group_by(thermap_presence, daysBack) %>%
  identify_outliers(msl)
# NO extrime outliers

#---(2.2)-Normality assumption
InSeaOnly %>%
  group_by(thermap_presence, daysBack) %>%
  shapiro_test(msl) 

#---(2.3)-Create QQ plot for each cell of design
ggqqplot(InSeaOnly, "msl", ggtheme = theme_bw()) +
  facet_grid(daysBack ~ thermap_presence, labeller = "label_both")

#-(3)--Computation


res.aov <- anova_test(data = InSeaOnly, dv = msl, wid = id,
                      within = c(thermap_presence, daysBack1))

#-(3.1)--Effect of thermals on each time point (which time points are different between groups)
one.way <- InSeaOnly %>%
  group_by(daysBack) %>%
  anova_test(dv = msl, wid = id, within = thermap_presence) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way


#-(3.2)--Effect of time at each  group (on which groups there is time effect)
one.way2 <- InSeaOnly %>%
  group_by(thermap_presence) %>%
  anova_test(dv = msl, wid = id, within = daysBack) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2

#-(3.3)--make Dunnett Test only on the significant group
bbThermals<-InSeaOnly[InSeaOnly$thermap_presence=="thermals",]
DunnettTest(msl~daysBack, data=bbThermals,control="0")


##################
#--(b)-cloud cover
##################
#-(1)--Create box plots of the score colored by treatment groups:
ggboxplot(
  InSeaOnly, x = "daysBack", y = "tcc",
  color = "thermap_presence", palette = "jco"
)

#-(2)--Chek assumptions
#---(2.1)-Outliers
Outlier <- InSeaOnly %>%
  group_by(thermap_presence, daysBack) %>%
  identify_outliers(tcc)
# NO extrime outliers

#---(2.2)-Normality assumption
InSeaOnly %>%
  group_by(thermap_presence, daysBack) %>%
  shapiro_test(tcc) 

#---(2.3)-Create QQ plot for each cell of design
ggqqplot(InSeaOnly, "tcc", ggtheme = theme_bw()) +
  facet_grid(daysBack ~ thermap_presence, labeller = "label_both")

#-(3)--Computation


res.aov <- anova_test(data = InSeaOnly, dv = tcc, wid = id,
                      within = c(thermap_presence, daysBack1))

#-(3.1)--Effect of thermals on each time point (which time points are different between groups)
one.way <- InSeaOnly %>%
  group_by(daysBack) %>%
  anova_test(dv = tcc, wid = id, within = thermap_presence) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way


#-(3.2)--Effect of time at each  group (on which groups there is time effect)
one.way2 <- InSeaOnly %>%
  group_by(thermap_presence) %>%
  anova_test(dv = tcc, wid = id, within = daysBack) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2

#-(3.3)--make Dunnett Test only on the significant group
bbThermals<-InSeaOnly[InSeaOnly$thermap_presence=="thermals",]
bbNoThermals<-InSeaOnly[InSeaOnly$thermap_presence=="No-thermals",]
DunnettTest(tcc~daysBack, data=bbThermals,control="0")
DunnettTest(tcc~daysBack, data=bbNoThermals,control="0")

##################
#--(b)-precipitation
##################
#-(1)--Create box plots of the score colored by treatment groups:
ggboxplot(
  InSeaOnly, x = "daysBack", y = "GPCC",
  color = "thermap_presence", palette = "jco"
)

#-(2)--Chek assumptions
#---(2.1)-Outliers
Outlier <- InSeaOnly %>%
  group_by(thermap_presence, daysBack) %>%
  identify_outliers(GPCC)
# NO extrime outliers

#---(2.2)-Normality assumption
InSeaOnly %>%
  group_by(thermap_presence, daysBack) %>%
  shapiro_test(GPCC) 

#---(2.3)-Create QQ plot for each cell of design
ggqqplot(InSeaOnly, "GPCC", ggtheme = theme_bw()) +
  facet_grid(daysBack ~ thermap_presence, labeller = "label_both")

#-(3)--Computation


res.aov <- anova_test(data = InSeaOnly, dv = GPCC, wid = id,
                      within = c(thermap_presence, daysBack1))

#-(3.1)--Effect of thermals on each time point (which time points are different between groups)
one.way <- InSeaOnly %>%
  group_by(daysBack) %>%
  anova_test(dv = GPCC, wid = id, within = thermap_presence) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way


#-(3.2)--Effect of time at each  group (on which groups there is time effect)
one.way2 <- InSeaOnly %>%
  group_by(thermap_presence) %>%
  anova_test(dv = GPCC, wid = id, within = daysBack) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2

#-(3.3)--make Dunnett Test only on the significant group
bbThermals<-InSeaOnly[InSeaOnly$thermap_presence=="thermals",]
DunnettTest(GPCC~daysBack, data=bbThermals,control="0")


##################
#--(b)-tail wind
##################
#-(1)--Create box plots of the score colored by treatment groups:
ggboxplot(
  InSeaOnly, x = "daysBack", y = "tw",
  color = "thermap_presence", palette = "jco"
)

#-(2)--Chek assumptions
#---(2.1)-Outliers
Outlier <- InSeaOnly %>%
  group_by(thermap_presence, daysBack) %>%
  identify_outliers(tw)


#---(2.2)-Normality assumption
InSeaOnly %>%
  group_by(thermap_presence, daysBack) %>%
  shapiro_test(tw) 

#---(2.3)-Create QQ plot for each cell of design
ggqqplot(InSeaOnly, "tw", ggtheme = theme_bw()) +
  facet_grid(daysBack ~ thermap_presence, labeller = "label_both")

#-(3)--Computation


res.aov <- anova_test(data = InSeaOnly, dv = tw, wid = id,
                      within = c(thermap_presence, daysBack1))

#-(3.1)--Effect of thermals on each time point (which time points are different between groups)
one.way <- InSeaOnly %>%
  group_by(daysBack) %>%
  anova_test(dv = tw, wid = id, within = thermap_presence) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way


#-(3.2)--Effect of time at each  group (on which groups there is time effect)
one.way2 <- InSeaOnly %>%
  group_by(thermap_presence) %>%
  anova_test(dv = tw, wid = id, within = daysBack) %>%
  get_anova_table() %>%
  adjust_pvalue(method = "bonferroni")
one.way2

#-(3.3)--make Dunnett Test only on the significant group
bbThermals<-InSeaOnly[InSeaOnly$thermap_presence=="thermals",]
DunnettTest(tw~daysBack, data=bbThermals,control="0")
#--------------------------------------------------------------------------------------------------------
##########################################################################################################

