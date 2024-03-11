#C. sessiliflora/Castilleja hawkmoth day/night exclusion experiment with 2019 and 2012/2013 datasets
#code for analyses in Wenzell, Zhang, Skogen, and Fant, bioRxiv, 2024
#collated for submission 20240310

library(tidyverse)
library(glmmTMB)
library(car)
library(ggpubr)
library(lsmeans)
library(DHARMa)
library(rcompanion)
library(rnaturalearth)


setwd("~/CASE_Rdocs/Cast_pollinators/CastPhD2.2_MS/Submission_1/")


#Read in data for pollinator exclusion experiment
ex4 <- read.csv("Cast_pollExclData_submission20240310.csv")
aggregate(frtSet ~ trt_Open:pop:majorPollinator, ex4, mean)
#remove bagged treatment from 2012/2013 for analysis (almost all zeroes and not included in 2019)
ex5 <- subset(ex4, !trt_Open=="bagged")
aggregate(frtSet ~ trt_Open:pop:majorPollinator, ex5, mean)

#MAP
#plot populations on map
map <- read.csv("Cast_pollExclPop_submission20240310.csv")

world <- ne_countries(scale = "medium", returnclass = "sf")

#Fig 1 map
map <- ggplot(data = world) +
  geom_sf() +
  geom_point(data = map, aes(x = lonW, y = latN, fill= majorPollPeriod, shape= corL), size= 4) +
  scale_fill_manual(values = c("nocturnal"= "#0074A2", "diurnal"= "#FFA500"))+
  scale_shape_manual(values=c('long'=21,'short'=24))+
  geom_text(data = map, aes(x = lonW, y = latN, label= pop), hjust= -.2, vjust= 0.1, size= 3)+
  coord_sf(xlim = c(-115, -83), ylim = c(25, 50), expand = FALSE)+
  theme_classic()
map

#some minor aesthetic edits were made to this map in InkScape

#FLORAL VISITATION DATA

#long- for visualization
vis <- read.csv("Cast_pollExclVisitationDataL_submission20240310.csv")
#wide- for summary table (Supplemental Table S3)
visW <- read.csv("Cast_pollExclVisitationData_submission20240310.csv")

#order of pops
popOrd <- c("LVH", "CQL", "SMP", "SIC", "SBL", "SCL", "SCC", "SILB12", "SILB13", "SNRM")
#order of pollinator groups
pollOrd <- c("hummingbird", "bumblebee", "small_bee", "otherD", "hawkmoth", "otherN")

#count of floral visits
#Fig 2B
vc <- ggplot(vis, aes(x=factor(pop, level=popOrd), y=pgVisit, fill=factor(pollFxGrp2, level= pollOrd)))+
  geom_bar(position= "stack", stat = "identity")+
  scale_fill_manual(values=c('hummingbird'="tomato2", 
                             'small_bee'="orange", 'bumblebee'="goldenrod1", 
                             "otherD" = "coral", 
                             'hawkmoth'="#0074A2", "otherN" = "#80B7C3"))+
  labs(y= "Number of floral visits", x= NULL, title= "B")+
  theme_classic()+
  theme(legend.position= "top", legend.title = element_blank())
vc

#FLORAL DATA FOR COROLLA LENGTH VISUALIZATION

#link pollinator data with floral data
ca <- read.csv("Cast_pollExcl_corL_submission20240310.csv")


#include SCC values and average by flower
sccf <- read.csv("Cast_SCC_floralCorL_perFlr.csv")


sccAv <- group_by(sccf, plant, pop, sp) %>%
  summarise(corL = mean(corL), LatN = mean(LatN), LonW = mean(LonW))

#combine with other populations
caa <- rbind(ca, sccAv)

#add corL category
caa$corLcat <- "long"
caa[caa$sp == "L", "corLcat"] <- "short"
caa[caa$sp == "C", "corLcat"] <- "short"
caa[caa$pop == "SIC", "corLcat"] <- "short"
caa[caa$pop == "SMP", "corLcat"] <- "short"

#Fig 2A
cL <- ggplot(caa, aes(x=factor(pop, level=popOrd), y=corL, pch= corLcat))+
  geom_violin()+
  stat_summary(fun.data=mean_sd, 
               #fun.args = list(mult=1), 
               geom="pointrange", color="black")+
  labs(y= "Corolla length (mm)", x= NULL, title = "A")+
  theme_classic()+
  theme(legend.position = "none")
cL

#####
#VISUALIZE AND ANALYZE EXCLUSION EXPERIMENT DATA

#Run GLMM with all data together, assess if fruit set varies by treatment with interactions of
#corolla length and period of activity of major pollinator

m <- glmmTMB(frtSet ~ trt_Open*corL*majorPollPeriod  + (1|pop)+ (1| plantID) + (1| year), weights=nFlrs, family = "betabinomial", data= ex5)
summary(m)
Anova(m)

#year included as random effect here because one population (SILB) repeated in multiple years
#subsetted datasets (below) havev no repeat years, so year RE not included


#check binomial distribution
m2 <- glmmTMB(frtSet ~ trt_Open*corL*majorPollPeriod  + (1|pop)+ (1| plantID) + (1| year), weights=nFlrs, family = "binomial", data= ex5)

AIC (m, m2) 
#m2 run with binomial, betabinomial performs better/ lower AIC; use m/ betabinomial

#assess model residuals with DHARMa
sim_resid_glmmTMB<- simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#ns- no significant problems/ deviations predicted
#test for zero inflation
testZeroInflation(sim_resid_glmmTMB)
#not significant evidence for zero inflation

#based on significant interactions of corolla length and major pollinator, subset data based on these categories
ex5l <- subset(ex5, corL== "long")
ex5s <- subset(ex5, corL== "short")

#SUBSET BY POLLINATOR CATEGORY
#"Predominant pollinator" based on the pollinator functional group that contributed the greatest 
#number of floral visits to the population during exclusion experiment

#long corolla and nocturnal poll
m <- glmmTMB(frtSet ~ trt_Open  + (1| pop) + (1|plantID), weights=nFlrs, family = "betabinomial", data= subset(ex5l, majorPollPeriod== "nocturnal"))
summary(m)
Anova(m)
DT <- lsmeans(m, pairwise ~ trt_Open, adjust= "tukey")
DT
#day < night < open- all sig

#compact letter display
PT <- DT$contrasts
PT <- as.data.frame(PT)
cldList(p.value ~ contrast, data= PT, threshold = 0.05)

#assess model residuals with DHARMa
sim_resid_glmmTMB<- simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#ns
testZeroInflation(sim_resid_glmmTMB)
#NS

#compare binomial
m1 <- glmmTMB(frtSet ~ trt_Open  + (1|pop) + (1|plantID), weights=nFlrs, family = "binomial", data= subset(ex5l, majorPollPeriod== "nocturnal"))
#compare models based on AIC
AIC(m, m1)
#choose BETABINOMIAL, lower AIC


#short corolla and nocturnal poll 
#no population random effect, because only one population in this dataset

#BETABINOMIAL
m <- glmmTMB(frtSet ~ trt_Open  + (1|plantID), weights=nFlrs, family = "betabinomial", data= subset(ex5s, majorPollPeriod== "nocturnal"))
summary(m)
Anova(m)
#NS- p=0.43

#assess model residuals with DHARMa
sim_resid_glmmTMB<- simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#ns
testZeroInflation(sim_resid_glmmTMB)
#NS

#compare binomial
m1 <- glmmTMB(frtSet ~ trt_Open  + (1|plantID), weights=nFlrs, family = "binomial", data= subset(ex5s, majorPollPeriod== "nocturnal"))
summary(m1)
Anova(m1)

#assess model residuals with DHARMa
sim_resid_glmmTMB<- simulateResiduals(m1, 1000)
plot(sim_resid_glmmTMB)
#ns
testZeroInflation(sim_resid_glmmTMB)
#NS

#compare models by AIC
AIC(m, m1)
#choose binomial, lower AIC (almost identical values)


#DIURNAL POLLINATORS
#long corolla and DIURNAL pollinators

#model won't converge with betabinomial distribution, so use binomial
m <- glmmTMB(frtSet ~ trt_Open  + (1| pop) + (1|plantID), weights=nFlrs, family = "binomial", 
             data= subset(ex5l, majorPollPeriod== "diurnal"))
summary(m)
Anova(m)
#p=0.01
DT <- lsmeans(m, pairwise ~ trt_Open, adjust= "tukey")
DT
#night < open- only significant comparison
#get compact letter display of pairwise differences
PT <- DT$contrasts
PT <- as.data.frame(PT)
cldList(p.value ~ contrast, data= PT, threshold = 0.05)
# ### Groups sharing a letter not significantly different (alpha = 0.05).

#assess model residuals with DHARMa
sim_resid_glmmTMB<- simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#decide to proceed based on inspecting plot and to keep analyses parallel/comparable
testZeroInflation(sim_resid_glmmTMB)

#short corolla and DIURNAL pollinators
m <- glmmTMB(frtSet ~ trt_Open  + (1| pop) + (1|plantID), weights=nFlrs, family = "betabinomial", 
             data= subset(ex5s, majorPollPeriod== "diurnal"))
summary(m)
Anova(m)

DT <- lsmeans(m, pairwise ~ trt_Open, adjust= "tukey")
DT
#night < day = open
#compact letter display
PT <- DT$contrasts
PT <- as.data.frame(PT)
cldList(p.value ~ contrast, data= PT, threshold = 0.05)


#assess model residuals with DHARMa
sim_resid_glmmTMB<- simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
testZeroInflation(sim_resid_glmmTMB)
#NS

#compare binomial
m1 <- glmmTMB(frtSet ~ trt_Open  + (1| pop) + (1|plantID), weights=nFlrs, family = "binomial",
              data= subset(ex5s, majorPollPeriod== "diurnal"))
AIC(m, m1)
#betabinomial model has lower AIC, use betabinomial


#combine pops into corL X pollinator period activity
ex6 <- mutate(ex5,
              popGroup = paste(corL, majorPollPeriod, sep = "_"))
#sample sizes
aggregate(trt_Open ~ popGroup, ex6, length)

ex6$cld <- "a              b              c"
ex6[ex6$popGroup == "short_nocturnal", "cld"] <- "a              a              a"
ex6[ex6$popGroup == "long_diurnal", "cld"] <- "ab            a              b"
ex6[ex6$popGroup == "short_diurnal", "cld"] <- "a              b              a"

#Fig 3
cmb <-   ggplot(data= ex6, aes(x= factor(trt_Open, level= c("day", "night", "open")), 
                               y= frtSet, pch= corL)) +
  geom_boxplot(aes(fill=trt_Open))+
  scale_fill_manual(values = c("day"= "#FFA500", "night"= "#0074A2", "open"= "darkslategray"))+
  facet_wrap(vars(factor(popGroup, level= c("long_nocturnal", "short_nocturnal",
                                            "long_diurnal", "short_diurnal" )), cld),  
             nrow=2)+
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", fill="darkgrey", size= .5)+
  scale_shape_manual(values = c(21, 24))+
  labs(y= "Proportion Fruit Set", x= "Treatment (period open)")+
  theme_classic()+
  theme(strip.background = element_blank(), strip.placement = "outside")
cmb

#Plot average difference in day and night Fruit set, relative to fully open FS by pop

#average by treatment and pop
exm <- group_by(ex5, pop, corL, majorPollPeriod, trt_Open) %>%
  summarise(avgFS = mean(frtSet))

#pivot wider: trt X FS
exmw <- exm %>% 
  pivot_wider(names_from = trt_Open, names_prefix = "FS_", values_from = avgFS) 

#calculate comparisons among treatments on average
exComp <- exmw %>%
  mutate(dayFSpr = FS_day / FS_open,
         nightFSpr = FS_night / FS_open,
         ndNorm = nightFSpr - dayFSpr)

#Plot difference in night FS - day FS (normalized by open FS)
#Fig 2C
ndn <- ggplot(data= exComp, aes(x= factor(pop, level= popOrd), y= ndNorm, shape = corL, col= majorPollPeriod))+
  geom_point(size= 3)+
  scale_color_manual(values=c('nocturnal'="#0074A2", 'diurnal'="#FFA500"))+
  labs(y= "Avg Normalized Night FS - Day FS", x= "Population", title = "C")+
  geom_hline(yintercept=0)+
  theme_classic()+
  theme(legend.position = "bottom", legend.title= element_blank())

ndn

#Plot open-only FS by pop (baseline)
#Fig 2D
op <- ggplot(data= subset(ex6, trt_Open=="open"), aes(x= factor(pop, level= popOrd), y= frtSet, shape= corL))+
  geom_point(size= 1, col= "grey")+
  stat_summary(fun.data=mean_se, 
               fun.args = list(mult=1), 
               geom="pointrange", color="black")+
  labs(y= "Open-treatment fruit set", x= "Population", title = "D")+
  theme_classic()+
  theme(legend.position = "none")

op

#combine FIGURE 2 panels
f2ab <- ggarrange(cL, vc, nrow=2, ncol=1, heights = c(0.7, 1))
f2ab
f2cd <- ggarrange(ndn, op, nrow=2, ncol=1, heights = c(1.2, 1), common.legend = T, legend= "bottom")
f2cd
ggarrange(f2ab, f2cd, nrow=2, ncol=1)


#Analyses using data from Wenzell et al., 2023, Oikos (doi: 10.1111/oik.09708)
#full data and code available here: https://github.com/KWenzell/Castilleja_pollinator_mosaics

#POPULATION AVERAGE FRUIT SET X COROLLA LENGTH AND VISITATION
fs1 <- read.csv("Wenzelletal2023_Oikos_fruitSet.csv")
fs1$year <- as.factor(fs1$year)
str(fs1)
#calculate fruit-flower ratio
fs1 <- fs1 %>% 
  mutate(FrtSet = nFrts / nFlrs)

unique(fs1$popYear)

#REMOVE populations with no observed visitation
fs2 <- subset(fs1, noVisitation== "0")
unique(fs2$popYear)

#long corollas
m <- glmmTMB(FrtSet ~ majorVisPeriod + (1|pop), weights=nFlrs, family = "betabinomial", data= subset(fs2, corLcat== "long"))
summary(m)
Anova(m)

#assess model residuals with DHARMa
sim_resid_glmmTMB<- simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
testZeroInflation(sim_resid_glmmTMB)

#short corollas
m <- glmmTMB(FrtSet ~ majorVisPeriod + (1|pop), weights=nFlrs, family = "betabinomial", data= subset(fs2, corLcat== "short"))
summary(m)
Anova(m)
#NS

#assess model residuals with DHARMa
sim_resid_glmmTMB<- simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#proceed after inspecting plots
testZeroInflation(sim_resid_glmmTMB)

#Fig 4A
fsS <- ggplot(data = fs2, aes(x= corLcat, y= FrtSet, fill= majorVisPeriod, pch= corLcat))+
  geom_violin(aes(fill= majorVisPeriod), trim= T, alpha= 0.9)+
  geom_text(data = NULL, x = 1, y = 1.01, label= "**")+
  geom_text(data = NULL, x = 2, y = 1.01, label= "n.s.")+
  stat_summary(fun.data=mean_se, fun.args = list(mult=1),
               geom="pointrange", color="black", position = position_dodge(width = 0.9))+
  scale_fill_manual(values=c('nocturnal'="#0074A2", 'diurnal'="#FFA500"))+
  labs(y= "Fruit set", x=NULL)+
  theme_classic()

fsS 

#####
####DIVERSITY OF VISITORS X COROLLA LENGTH
myp1 <- read.csv("Wenzelletal2023_Oikos_visitDiv.csv")

#log transform InvSimpD (Inverse Simpson's Diversity Index)
myp2 <- myp1 %>%
  mutate(logInvSimpD = log(ct_InvSimpD))


m <- glmmTMB(logInvSimpD ~ corLcat + (1|pop) + (1|year) + dataset, data= myp2, ziformula = ~1,
             control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(m) 
Anova(m) 

sim_resid_glmmTMB<- simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
testZeroInflation(sim_resid_glmmTMB)
#significant, but accounted for using zero inflation parameter, use of glmmTMB, and further addressed with nonparametric KW test below

#average across observation datasets for Kruskal Wallis test (no fixed effects possible)
mypD <- myp2 %>%
  dplyr::select(popYr, year, pop, dataset, sp, latN, lonW, corLcat, logInvSimpD)
#remove Nas (where 0 visitors recorded)
mypD <- na.omit(mypD)
mypAvgD <- group_by(mypD, pop, year, corLcat, sp, latN, lonW) %>%
  summarise(avgLogInvSimpD = mean(logInvSimpD),
            sdLISD = sd(logInvSimpD))

kruskal.test(avgLogInvSimpD ~ corLcat, data= mypAvgD)

#CHECK LATITUDE
m <- glmmTMB(logInvSimpD ~ latN + (1|year) + dataset, data= myp2, ziformula = ~1,
             control=glmmTMBControl(optimizer=optim, optArgs=list(method="BFGS")))
summary(m) 
Anova(m) 

sim_resid_glmmTMB<- simulateResiduals(m, 1000)
plot(sim_resid_glmmTMB)
#proceed after inspection of plots
testZeroInflation(sim_resid_glmmTMB)
#see above

#check with KW test
kruskal.test(avgLogInvSimpD ~ latN, data= mypAvgD)


#Fig 4B
div <- ggplot(data = myp2, aes(x= corLcat, y= logInvSimpD, pch= corLcat))+
  geom_violin(trim= T)+
  geom_text(data = NULL, x = 1.5, y = 1.3, label= "*")+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5, alpha= 0.5)+
  stat_summary(fun.data=mean_se, fun.args = list(mult=1), 
               geom="pointrange", color="black")+
  labs(x= "Corolla length", y= "Pollinator Diversity")+
  theme_classic()
div

#combine Figure 4 panels
ggarrange(fsS, div, nrow=2, ncol=1, common.legend = T)



#Supplemental Figure S2: boxplots of exclusion experiment results by population
popEx <-   ggplot(data= ex5, 
                 aes(x= factor(trt_Open, level= c("day", "night", "open")), y= frtSet)) +
  geom_boxplot(aes(fill=trt_Open))+
  scale_fill_manual(values = c("#FFA500", "#0074A2", "darkslategray"))+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=.5)+
  facet_wrap(vars(factor(pop, level= c(popOrd)), corL, majorPollPeriod), 
             nrow=3)+
  labs(y= "Proportion Fruit Set", x= "Treatment (period open)")+
  theme_classic()+
  theme(strip.background = element_blank(), strip.placement = "outside")
popEx
