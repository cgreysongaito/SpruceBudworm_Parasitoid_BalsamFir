#############################################################################
#            Coherent whole food web responses to outbreaking budworm
#
#     Christopher J. Greyson-Gaito, Kevin S. McCann, Jochen Fruend, Christopher J. Lucarotti, M. Alex Smith, Eldon S. Eveleigh
#
#
#             R script for statistical analysis and figure plotting
#
###################################################################


# Packages -----

rm(list=ls())

theme_simple <- function () { 
  theme_grey() %+replace% 
    theme(
      axis.line=element_line(colour="black"),
      panel.grid.minor=element_blank(), 
      panel.grid.major=element_blank(),
      panel.background=element_blank(), 
      axis.title=element_text(size=28,face="bold"),
      axis.text.x=element_text(size=24, colour="Black"),
      axis.text.y=element_text(size=24, colour="Black"),
      axis.ticks.length=unit(0.5,"cm"),
      axis.ticks=element_line(size=0.5, colour="Black"),
      panel.border=element_rect(fill=FALSE,size=0.5),
      legend.title=element_text(size=15),
      legend.key=element_blank()
    )
  
}

library(lattice)
library(magrittr)
library(tidyverse); theme_set(theme_simple())
library(grid)
library(gridExtra)
library(nlme)
library(viridis)
library(rms)
library(bipartite)
library(vegan)
library(goeveg)
library(indicspecies)
library(permute)
library(stargazer)

# Data Preparation -----
##Input data as data.table
interactionimport<-read_csv("data/InteractFinal_long_2014Nov10.csv") #can use read_csv and will change duo to NA but still 300 parsing failures in other columns (not important columns for this analysis)
speciesnamesimport<-read_csv("data/TotalSpeciesTable_2014Nov10.csv")

##Data cleaning
interaction<-interactionimport%>%
  mutate(Plot=as.factor(Plot))%>%
  filter(!Plot==4)%>%
  filter(CrownLevel=="mid") %>% #using mid crown because royama 2017 used mid crown and because eveleigh and johns 2014 found that using mid crown samples are the best proxies for estimating budworm and parasitoid densities
  filter(!SpecID.low=="h11")%>% #remove sawfly from dataset
  filter(!SpecID.low=="h00")%>% #remove undetermined herbivores from dataset (so similar to Eveleigh 2007 SI Materials and Methods)
  filter(!SpecID.high %in% c("pA3", "p40","p05", "p07","p04", "p51"))%>% #remove all unresolved taxa #remove diadegma sp. that has freq.low of 2 and uncertain as to the species (either D. pulicalivariae-p94 or Diadegma sp1-p95)
  mutate(Plot=droplevels(Plot))%>%
  mutate(Peak=ifelse(Plot=="3",Year-6,Year))%>% 
  mutate(Peak=Peak-85) #plot 1 and 2 peaked in 85 and plot 3 peaked in 91. Therefore created new variable Peak to standardize the different plots by when they peak

  parashareprep<-interaction%>%
  filter(!grepl("p",SpecID.low))%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  select(Peak,Plot,SBWALT,SpecID.high)%>%
  filter(!is.na(SpecID.high))%>%
  split(.$SBWALT)%>%
  map(function(x) {
    tibble(
      UnSP=unique(x$SpecID.high),
      SBWALT=rep(unique(x$SBWALT),each=length(UnSP))
    )
  })%>%
  bind_rows #produces a list of unique parasitoid taxa that attacked budworm and other caterpillars

parashare<-parashareprep$UnSP[duplicated(parashareprep$UnSP)] #produces a list of the parasitoids that attack both budworm and other caterpillars

create.numSBWALTnparaprep<-function(plot){
  x<-data.frame(expand.grid(c(-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10),c("SBW","ALT")))%>%
    rename(Peak=Var1,SBWALT=Var2)%>%
    mutate(Plot=plot)
  
  if (plot==1){
    x2<-x%>%filter(Peak<5)
  } else if (plot==2){
    x2<-x%>%filter(Peak>0)
  } else {
    x2<-x%>%filter(Peak<4)
  }
}

numSBWALTnparaprep<-c(1,2,3)%>%
  map(create.numSBWALTnparaprep)%>%
  bind_rows()%>%
  mutate(Plot=as.factor(Plot)) #combined with the function create.numSBWALTnparaprep this creates a dataframe with all combinations of Peak (-3 to 10), SBWALT (either SBW or ALT), and Plot (1,2,3) to be ensure all data sets used in analyses have all combinations

speciesnamesshare<-speciesnamesimport%>%
  select(speciesID,origName, Genus, Species)%>%
  filter(grepl("p",speciesID))%>%
  filter(speciesID %in% parashare)%>%
  unite(LatinName, Genus,Species,sep=" ",remove=TRUE)%>%
  rename(SpecID.high=speciesID)%>%
  mutate(LatinName=ifelse(grepl("NA",LatinName),origName,LatinName))%>%
  select(-origName)

speciesnamesshare$LatinName<-gsub("Winthemia","Smidtia", speciesnamesshare$LatinName)
speciesnamesshare$LatinName<-gsub("setifacies","fumipennis", speciesnamesshare$LatinName)
speciesnamesshare$LatinName<-gsub("Ephialtes","Apechthis", speciesnamesshare$LatinName)
speciesnamesshare$LatinName<-gsub("Phaeogenes","Dirophanes", speciesnamesshare$LatinName)
speciesnamesshare$LatinName<-gsub("Macrocentrus","Macrocentrus linearis iridescens", speciesnamesshare$LatinName)
speciesnamesshare$LatinName<-gsub("microgaster","Microgaster", speciesnamesshare$LatinName)
speciesnamesshare$LatinName<-gsub("Glypta sp1","Glypta sp.", speciesnamesshare$LatinName)
speciesnamesshare$LatinName<-gsub("Orgilis","Orgilus", speciesnamesshare$LatinName)
speciesnamesshare$LatinName<-gsub("Synetaeris","Tranosema",speciesnamesshare$LatinName)
speciesnamesshare$LatinName<-gsub("cf. ","",speciesnamesshare$LatinName)
speciesnamesshare$LatinName<-gsub("sp1","sp.",speciesnamesshare$LatinName)

#From speciesnamesshare<- to here, this section creates a list of the latin names of all the parasitoids that attack both budworm and other caterpillars and corrects the latin names if necessary.

propSBWALTbyspeciesdata<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  group_by(SpecID.high)%>%
  summarise(Tot=length(SpecID.high), NSBW=length(SpecID.high[SBWALT=="SBW"]),NALT=length(SpecID.high[SBWALT=="ALT"]))%>%
  mutate(PropSBW=NSBW/Tot)%>%
  mutate(stanerror=sqrt(PropSBW*(1-PropSBW)/Tot), pymax=PropSBW+stanerror, pymin=PropSBW-stanerror)%>% #calculate standard error of the proportion
  merge(speciesnamesshare,by="SpecID.high",all.x=TRUE)%>%
  select(-SpecID.high)%>%
  rename(SpecID.high=LatinName)%>%
  mutate(SpecID.high=as.factor(SpecID.high)) #This section produces a dataset of the Proportion of emergences from budworm for each parasitoid taxon. also gives the proportion standard error and error bars for the standard error.

proplevels<-levels(reorder(propSBWALTbyspeciesdata$SpecID.high,propSBWALTbyspeciesdata$PropSBW)) #produces a list where the parasitoid taxa are sorted by the proportion of emergences from budworm.

speciesnamesprep<-tibble(
  SpecID.high=as.factor(rev(proplevels)),
  numID=seq(1,31,by=1)
) #produces a dataset where each parasitoid taxon has a number. used in producing the bipartite graphs.

#Parasitoid taxa diversity -----

nparaSBWALT<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  mutate(SBWALT=as.factor(ifelse(SpecID.low=="h01","SBW","ALT")))%>%
  group_by(Peak,Plot, SBWALT)%>%
  summarise(npara=length(unique(SpecID.high)))%>%
  ungroup%>%
  right_join(numSBWALTnparaprep,by=c("Peak","Plot","SBWALT"))%>%
  mutate(dummyPeak=Peak+4,npara=as.numeric(npara),Plot=as.factor(Plot),SBWALT=as.factor(SBWALT), npara=ifelse(is.na(npara),0,npara), label=as.factor(ifelse(SBWALT=="SBW", "budworm","other caterpillars")))#Creates a dataframe of the overall number of parasitoid taxa that emerge from budworm and other caterpillars in every year and for each plot.

#Linear mixed effects model of the number of parasitoid taxa from budworm and other caterpillars over time.

#All of the statistical analysis below follows the protocol outlined in Zuur, A., E. N. Ieno, N. Walker, A. A. Saveliev, and G. M. Smith. 2009. Mixed effects models and extensions in ecology with r. First edition. Springer-Verlag New York, New York, New York, United States of America.

#Data Exploration
#1.Outliers in the response and explanatory variables.
op <- par(mfrow = c(1, 2), mar = c(3, 3, 3, 1))
dotchart(nparaSBWALT$npara, main = "Number of Species", group=nparaSBWALT$SBWALT)
dotchart(nparaSBWALT$Peak, main = "Peak")
par(op)
#2. Collinearity of the explanatory variables.
#Not doing as do not have multiple continuous explanatory variables

#3. Relationships between the response variable and the explanatory variables.
xyplot(npara~Peak, data=nparaSBWALT, groups=SBWALT)
#Creating the model
#1.Start with a linear regression model that contains as many explanatory variables and their interactions as possible. 
SBWALTnparafulllinear<-lm(npara~Peak*SBWALT*Plot,data=nparaSBWALT)
nparaSBWALT$resid<-NA
nparaSBWALT$resid<-SBWALTnparafulllinear$residuals

#Test normality and homogeniety of residuals
op <- par(mfrow = c(3, 3), mar = c(5, 4, 1, 2))
plot(SBWALTnparafulllinear, add.smooth = FALSE, which = 1,col=nparaSBWALT$SBWALT)
hist(nparaSBWALT$resid, xlab = "Residuals", main = "")
plot(nparaSBWALT$Peak,nparaSBWALT$resid,xlab="Peak",ylab="Residuals")
plot(nparaSBWALT$Plot,nparaSBWALT$resid,xlab="Plot",ylab="Residuals")
boxplot(nparaSBWALT$resid~nparaSBWALT$SBWALT,xlab="Caterpillar",ylab="Residuals")
plot(nparaSBWALT$Peak[nparaSBWALT$SBWALT=="SBW"],nparaSBWALT$resid[nparaSBWALT$SBWALT=="SBW"],xlab="Peak",ylab="Residuals",main="SBW")
plot(nparaSBWALT$Peak[nparaSBWALT$SBWALT=="ALT"],nparaSBWALT$resid[nparaSBWALT$SBWALT=="ALT"],xlab="Peak",ylab="Residuals",main="ALT")
par(op)

acf(residuals(SBWALTnparafulllinear), na.action = na.pass,
    main = "Auto-correlation plot for residuals")
#Obviously violation of independence of x values (time) but acf graph does not show high correlation between years

#2. Repeat step 1 using the gls function from the nlme package using REML estimation but without any variance structures.
SBWALTnparafullgls<-gls(npara~dummyPeak*SBWALT*Plot,data=nparaSBWALT, method="REML")

#3. Choose an appropriate variance structure depending on graphical analysis above.
SBWALTnparafullvaridentgls<-gls(npara~dummyPeak*SBWALT*Plot,weights=varIdent(form=~1|SBWALT),data=nparaSBWALT, method="REML")
anova(SBWALTnparafullgls,SBWALTnparafullvaridentgls)
#using varIdent does not improve model (AIC actually increases)

SBWALTnparaautocorr<-gls(npara~dummyPeak*SBWALT*Plot,correlation=corAR1(form=~dummyPeak|Plot/SBWALT),data=nparaSBWALT, method="REML")


anova(SBWALTnparafullgls,SBWALTnparaautocorr)#Delta AIC=-5

acf(residuals(SBWALTnparaautocorr, type = "normalized"), na.action = na.pass,
    main = "Auto-correlation plot for residuals")

#reduction in AIC with autocorrelation (greater than 2). small improvement when looking at ACF graph.
#decided to use autocorrelation

#4. Using the gls model with the selected variance structure from 3., try to find the optimal random structure.
#Not doing because Plot is a fixed effect (using reasoning from http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html)

#5. Compare the new gls model with the earlier results using AIC.

#6.Find the optimal ﬁxed component using the likelihood ratio test where ML estimation is used.
SBWALTnparafullml<-gls(npara~dummyPeak*SBWALT*Plot,correlation=corAR1(form=~dummyPeak|Plot/SBWALT),data=nparaSBWALT, method="ML")
SBWALTnparadrop3interactionml<-update(SBWALTnparafullml,.~.-dummyPeak:SBWALT:Plot)
anova(SBWALTnparafullml,SBWALTnparadrop3interactionml)

#drop interaction of dummyPeak,SBWALT, and Plot because not significant

SBWALTnparadropdummyPeakSBWALTml<-update(SBWALTnparadrop3interactionml,.~.-dummyPeak:SBWALT)
anova(SBWALTnparadrop3interactionml,SBWALTnparadropdummyPeakSBWALTml)

SBWALTnparadropdummyPeakPlotml<-update(SBWALTnparadrop3interactionml,.~.-dummyPeak:Plot)
anova(SBWALTnparadrop3interactionml,SBWALTnparadropdummyPeakPlotml)

SBWALTnparadropSBWALTPlotml<-update(SBWALTnparadrop3interactionml,.~.-SBWALT:Plot)
anova(SBWALTnparadrop3interactionml,SBWALTnparadropSBWALTPlotml)

#drop interaction of SBWALT and Plot because least significant

SBWALTnparadropdummyPeakSBWALTml2<-update(SBWALTnparadropSBWALTPlotml,.~.-dummyPeak:SBWALT)
anova(SBWALTnparadropSBWALTPlotml,SBWALTnparadropdummyPeakSBWALTml2)

SBWALTnparadropdummyPeakPlotml2<-update(SBWALTnparadropSBWALTPlotml,.~.-dummyPeak:Plot)
anova(SBWALTnparadropSBWALTPlotml,SBWALTnparadropdummyPeakPlotml2)

#do not drop interaction of dummyPeak and SBWALT. do not drop interaction of dummyPeak and Plot.

#7. Refit the model found in step 6. with REML estimation
SBWALTnparafinalmodelreml<-gls(npara~dummyPeak+SBWALT+Plot+dummyPeak:SBWALT+dummyPeak:Plot,correlation=corAR1(form=~dummyPeak|Plot/SBWALT),data=nparaSBWALT, method="REML")

summary(SBWALTnparafinalmodelreml)

# Calculation of difference in mean parasitoid diversity between budworm and other caterpillars
nparaSBWALT%>%
  group_by(SBWALT)%>%
  summarise(MeanSBWALT=mean(npara))%>%
  spread(SBWALT,MeanSBWALT)

# Create graph of number of parasitoid taxa emerging from budworm and other caterpillars over time,
SBWALTnparaautocorrremlGLS<-Gls(npara~dummyPeak+SBWALT+Plot+dummyPeak:SBWALT+dummyPeak:Plot,correlation=corAR1(form=~dummyPeak|Plot/SBWALT),data=nparaSBWALT, method="REML") #using the function Gls because allows Predict function needed to plot trendline
nparapredict<-Predict(SBWALTnparaautocorrremlGLS, dummyPeak=c(1:14), Plot=seq(from=1,to=3,by=1),SBWALT=c("ALT","SBW"))%>%
  mutate(label=as.factor(ifelse(SBWALT=="SBW", "budworm","other caterpillars")),Peak=dummyPeak-4,Plot=as.factor(Plot))

nparaSBWALTpredict<-left_join(nparaSBWALT,nparapredict,by=c("Plot", "Peak", "SBWALT", "label"))

nparaSBWALTplot<-ggplot(nparaSBWALTpredict)+
  geom_vline(xintercept=0)+
  geom_line(aes(Peak,yhat,linetype=Plot), size=2)+
  geom_point(aes(Peak,npara,shape=Plot), size=5)+
  facet_wrap(~label,scales="fixed",ncol=1)+
  theme(strip.background=element_blank(),strip.text.x=element_text(size=28),
        axis.title.y=element_text(hjust=0.5, vjust=1.5))+
  ylab("Richness of \nparasitoid taxa")+xlab("Years before/after peak")+scale_y_continuous(limits=c(0, 26))+scale_x_continuous(breaks=c(-2,0,2,4,6,8,10))+scale_shape_manual(values=c(16,17,15))+scale_linetype_manual(values=c("dotted","dashed","solid"))

ggsave("figs/nparaSBWALTplot.pdf",plot=nparaSBWALTplot,width=8,height=7) #Figure 1

#Parasitoid community host preference-----

#community level plotting of diet switching using log10 ratios of emergences from budworm and other caterpillars and log10 ratios of total abundances of budworm and other caterpillars

ratioSBWbypeakplot<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  group_by(Peak, Plot)%>%
  summarise(TotSBW=length(SBWALT[SBWALT=="SBW"]),TotALT=length(SBWALT[SBWALT=="ALT"]))%>%
  mutate(ratioSBWpeakplot=TotSBW/TotALT, logratioSBWpeakplot=log10(ratioSBWpeakplot)) #produces a dataset of the ratio (and log 10 ratio) of the total number of budworm and other caterpillars sampled on balsam fir

ratioSBWparacommunity<-interaction %>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  group_by(Peak, Plot)%>%
  summarise(TotSBWEmerge=length(SpecID.high[SBWALT=="SBW"]),TotALTEmerge=length(SpecID.high[SBWALT=="ALT"]))
   #produces a data set of the ratio (and log 10 ratio) of the total number of parasitoid emergences (all taxa combined) from budworm and other caterpillars sampled on balsam fir.

communitydietswitch<-full_join(ratioSBWparacommunity,ratioSBWbypeakplot,by=c("Peak","Plot"))%>%
  filter(TotSBW!=0,TotALT!=0,TotSBWEmerge!=0,TotALTEmerge!=0)%>% #Removed four data points (from a total of 25 data points).
  ungroup()%>%
  mutate(TotSBWEmergepcap=TotSBWEmerge/TotSBW, TotALTEmergepcap=TotALTEmerge/TotALT, ratioSBWEmerge=TotSBWEmerge/TotALTEmerge,logratioSBWEmerge=log10(ratioSBWEmerge),ratioSBWEmergepcap=TotSBWEmergepcap/TotALTEmergepcap,logratioSBWEmergepcap=log10(ratioSBWEmergepcap)) %>%
  mutate(dummyPeak=Peak+4, Plot=as.factor(Plot)) #creates a dummy variable because function gls can not use negative numbers.

###Linear mixed effects model of log10 ratio parasitoid emergence from budworm and other caterpillars to the log10 ratio of total budworm to other caterpillars. Using the diet switching method outlined by Greenwood, J. J. D., and R. A. Elton. 1979. Analysing experiments on frequency-dependent selection by predators. The Journal of Animal Ecology 48:721.
 
#All of the statistical analysis below follows the protocol outlined in Zuur, A., E. N. Ieno, N. Walker, A. A. Saveliev, and G. M. Smith. 2009. Mixed effects models and extensions in ecology with r. First edition. Springer-Verlag New York, New York, New York, United States of America.

#Data Exploration
#1.Outliers in the response and explanatory variables.

dotchart(communitydietswitch$logratioSBWEmerge, main = "Log ratio SBW Emerge Per Capita")
dotchart(communitydietswitch$logratioSBWpeakplot, main = "Log Ratio SBW:ALT")

#2. Collinearity of the explanatory variables.
#Not doing as do not have multiple continuous explanatory variables

#3. Relationships between the response variable and the explanatory variables.
xyplot(logratioSBWEmerge~logratioSBWpeakplot, data=communitydietswitch, groups=Plot) #not percapita
xyplot(TotSBWEmerge~Peak, data=communitydietswitch)
xyplot(TotALTEmerge~Peak, data=communitydietswitch)
#Creating the model
#1.Start with a linear regression model that contains as many explanatory variables and their interactions as possible. 
SBWALTdietswitchfulllinear<-lm(logratioSBWEmerge~logratioSBWpeakplot*Plot,data=communitydietswitch)
communitydietswitch$resid<-NA
communitydietswitch$resid<-SBWALTdietswitchfulllinear$residuals
#Test normality and homogeniety of residuals

plot(SBWALTdietswitchfulllinear, add.smooth = FALSE, which = 1,col=communitydietswitch$Plot)
hist(communitydietswitch$resid, xlab = "Residuals", main = "")
plot(communitydietswitch$Peak,communitydietswitch$resid,xlab="Peak",ylab="Residuals",col=as.factor(communitydietswitch$Plot))
plot(communitydietswitch$Plot,communitydietswitch$resid,xlab="Plot",ylab="Residuals")


acf(residuals(SBWALTdietswitchfulllinear), na.action = na.pass,
    main = "Auto-correlation plot for residuals")
#Obviously violation of independence of x values (time) and repeated measures of the same plot in different years (although acf does not show much correlation between years)

#2. Repeat step 1 using the gls function from the nlme package using REML estimation but without any variance structures.
SBWALTdietswitchgls<-gls(logratioSBWEmerge~logratioSBWpeakplot*Plot,data=communitydietswitch, method="REML")
#3. Choose an appropriate variance structure depending on graphical analysis above.

SBWALTdietswitchautocorr<-gls(logratioSBWEmerge~logratioSBWpeakplot*Plot,data=communitydietswitch,correlation=corAR1(form=~dummyPeak|Plot), method="REML")

anova(SBWALTdietswitchgls,SBWALTdietswitchautocorr)#Delta AIC<2

#check assumptions graphically for model with autocorrelation
communitydietswitch$fmnormresid<-NA
communitydietswitch$fmnormresid<-resid(SBWALTdietswitchautocorr, type="normalized")

plot(SBWALTdietswitchautocorr, add.smooth = FALSE)
hist(communitydietswitch$fmnormresid, xlab = "Normalised Residuals", main = "")
plot(communitydietswitch$Peak,communitydietswitch$fmnormresid,xlab="Peak",ylab="Normalised Residuals")
plot(communitydietswitch$Peak,communitydietswitch$fmnormresid,xlab="Peak",ylab="Residuals",col=as.factor(communitydietswitch$Plot))
plot(communitydietswitch$Plot,communitydietswitch$fmnormresid,xlab="Plot",ylab="Residuals")
acf(communitydietswitch$fmnormresid, na.action = na.pass,
    main = "Auto-correlation plot for normalised residuals")

#AIC lowest for model with autocorr but tiny change in ACF (and ACF does not show large improvement, also residuals are less normally distributed). Changes too small so not including corAR1.

#4. Using the gls model with the selected variance structure from 3., try to find the optimal random structure.
#Not doing because Plot is a fixed effect (using reasoning from https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html)

#6.Find the optimal ﬁxed component using the likelihood ratio test where ML estimation is used.
SBWALTdietswitchmlfull<-gls(logratioSBWEmerge~logratioSBWpeakplot*Plot,data=communitydietswitch, method="ML")

SBWALTdietswitchmlminusinteraction<-gls(logratioSBWEmerge~logratioSBWpeakplot+Plot,data=communitydietswitch, method="ML")
anova(SBWALTdietswitchmlfull,SBWALTdietswitchmlminusinteraction)

# Interaction is significant (p-value = 0.0011, delta AIC = 9.60224). Therefore, keep all fixed effects in the model including the interaction.


#7. Refit the model found in step 6. with REML estimation
SBWALTdietswitchreml<-gls(logratioSBWEmerge~logratioSBWpeakplot*Plot,data=communitydietswitch, method="REML")

summary(SBWALTdietswitchreml)

# average intercept for all plots
avinterceptallplotsDS<-(((3*coef(SBWALTdietswitchreml)[1])+coef(SBWALTdietswitchreml)[3]+coef(SBWALTdietswitchreml)[4])/3)

# average standard error of intercepts for all plots
avSEinterceptallplotsDS <- ((summary(SBWALTdietswitchreml)$tTable[1,2]+summary(SBWALTdietswitchreml)$tTable[3,2]+summary(SBWALTdietswitchreml)$tTable[4,2])/3)

#confidence intervals for average intercept
avSEinterceptallplotsDS*1.96

#Test of whether the average intercept of all plots is different from 0 (test of caterpillar type discrimination)
2*pt(abs((avinterceptallplotsDS-0)/avSEinterceptallplotsDS),df=SBWALTdietswitchreml$dims[[1]]-SBWALTdietswitchreml$dims[[2]],lower.tail = FALSE) #p>0.05 so intercept is  not significantly different from 0


#average slope for all plots
avslopeallplotsDS<-(((3*coef(SBWALTdietswitchreml)[2])+coef(SBWALTdietswitchreml)[5]+coef(SBWALTdietswitchreml)[6])/3)
  
#average standard error of slopes for all plots
avSEslopeallplotsDS<-((summary(SBWALTdietswitchreml)$tTable[2,2]+summary(SBWALTdietswitchreml)$tTable[5,2]+summary(SBWALTdietswitchreml)$tTable[6,2])/3)

#confidence intervals for average slope
avSEslopeallplotsDS*1.96

#Test of whether the average slope of all plots is different from 1 (test of resource tracking)
2*pt(abs((avslopeallplotsDS-1)/avSEslopeallplotsDS),df=SBWALTdietswitchreml$dims[[1]]-SBWALTdietswitchreml$dims[[2]],lower.tail = FALSE) #p<0.05 so b is significantly different from 0

#Create graph of log10 ratio parasitoid emergence from budworm and other caterpillars to the log10 ratio of total budworm to other caterpillars
SBWALTdietswitchremlGLS<-Gls(logratioSBWEmerge~logratioSBWpeakplot*Plot,data=communitydietswitch, method="REML") #using the function Gls because allows Predict function needed to plot trendline

dietswitchpredict<-Predict(SBWALTdietswitchremlGLS, logratioSBWpeakplot=seq(from=-1.5,to=2.1,by=0.1),Plot=seq(from=1,to=3,by=1))%>%
  mutate(Plot=as.factor(Plot))

comswitchlogratioplot<-ggplot(communitydietswitch,aes(logratioSBWpeakplot, logratioSBWEmerge))+
  geom_abline(intercept=0,slope=1,linetype="twodash")+
  geom_abline(intercept=0,slope=0)+
  geom_vline(xintercept=0)+
  geom_line(data=dietswitchpredict,aes(logratioSBWpeakplot,yhat,linetype=Plot),size=2)+
  geom_point(aes(shape=Plot),size=5)+
  theme(axis.title.y=element_text(hjust=0.5, vjust=1.5), 
        legend.text=element_text(size=14),legend.justification=c(1,0), legend.position=c(0.98,0.05),legend.box = "horizontal")+
  coord_fixed(ratio = 1)+
  scale_shape_manual(values=c(16,17,15))+scale_linetype_manual(values=c("dotted","dashed","solid"))+
  ylab("Log10 emergence \nbudworm:other caterpillars")+xlab("Log10 abundance \nbudworm:other caterpillars")

ggsave("figs/communitydietswitchlogratio.pdf", plot=comswitchlogratioplot, width=8,height=8) #Figure 2a

# Per capita parasitoid emergences from budworm and other caterpillars

#Together, the function create.SBWALTparashareprep and the short code that creates the data frame called SBWALTparashareprep create a data frame where there is a row for every combination of peak, caterpillar type and parasitoid taxa but where there aren't rows for when the plots were not sampled (different plots were sampled in different years).
create.SBWALTparashareprep<-function(plot){
  x<-data.frame(expand.grid(c(-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10),c("SBW","ALT")))%>%
  rename(Peak=Var1,SBWALT=Var2)%>%
  mutate(Plot=plot)
  
  if (plot==1){
  x2<-x%>%filter(Peak<5)
  } else if (plot==2){
    x2<-x%>%filter(Peak>0)
  } else {
    x2<-x%>%filter(Peak<4)
  }
}

SBWALTparashareprep<-c(1,2,3)%>%
  map(create.SBWALTparashareprep)%>%
  bind_rows()%>%
  mutate(Plot=as.factor(Plot))

pcapprep <- ratioSBWbypeakplot %>%
  select(Peak,Plot,TotSBW,TotALT) %>%
  rename(SBW=TotSBW,ALT=TotALT)%>%
  gather(SBWALT,Tot,-Peak,-Plot)

SBWALTEM<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  group_by(Peak,Plot)%>%
  summarise(SBW=length(SpecID.low[SBWALT=="SBW"]),ALT=length(SpecID.low[SBWALT=="ALT"]))%>%
  gather(SBWALT,NEM,-Peak,-Plot)%>%
  ungroup%>%
  right_join(SBWALTparashareprep,by=c("Peak","Plot","SBWALT"))%>%
  mutate(Plot=as.factor(Plot),SBWALT=as.factor(SBWALT), NEM=ifelse(is.na(NEM),0,NEM))%>%
  left_join(pcapprep, by=c("Peak","Plot","SBWALT")) %>%
  mutate(NEMpcap=NEM/Tot) %>%
  filter(!is.nan(NEMpcap)) %>%
  #mutate(NEMpcap = ifelse(is.nan(NEMpcap),0,NEMpcap))%>%
 mutate(dummyPeak=Peak+4,label=as.factor(ifelse(SBWALT=="SBW", "budworm","other caterpillars")),SBWALT=as.factor(SBWALT))  #Creates a dataframe of the average number of emergences of parasitoids from budworm or other caterpillars (all parasitoid taxa combined) for each peak and plot. Before taking the average, 0 values are added to any peak, plot, taxa combination that not have a parasitoid emerge from either budworm or other caterpillars.

xyplot(NEMpcap~Peak, data=SBWALTEM, groups=SBWALT)
xyplot(NEM~Peak, data=SBWALTEM, groups=SBWALT)
#Linear mixed effects model of log10 average number of emergences of parasitoids from budworm or other caterpillars over time

#All of the statistical analysis below follows the protocol outlined in Zuur, A., E. N. Ieno, N. Walker, A. A. Saveliev, and G. M. Smith. 2009. Mixed effects models and extensions in ecology with r. First edition. Springer-Verlag New York, New York, New York, United States of America.

#Data Exploration
#1.Outliers in the response and explanatory variables.
op <- par(mfrow = c(2, 2), mar = c(3, 3, 3, 1))
dotchart(SBWALTEM$NEMpcap, main = "Number of Interactions", group=SBWALTEM$SBWALT)
dotchart(SBWALTEM$Peak, main = "Peak")
par(op)
#2. Collinearity of the explanatory variables.
#Not doing as do not have multiple continuous explanatory varaibles

#3. Relationships between the response variable and the explanatory variables.
xyplot(NEMpcap~Peak, data=SBWALTEM, groups=SBWALT)

#Creating the model
#1.Start with a linear regression model that contains as many explanatory variables and their interactions as possible. 
SBWALTEMfulllinear<-lm(NEMpcap~Peak*SBWALT*Plot,data=SBWALTEM)
SBWALTEM$resid<-NA
SBWALTEM$resid<-SBWALTEMfulllinear$residuals
#Test normality and homogeniety of residuals
op <- par(mfrow = c(3, 3), mar = c(5, 4, 1, 2))
plot(SBWALTEMfulllinear, add.smooth = FALSE, which = 1,col=SBWALTEM$SBWALT)
hist(SBWALTEM$resid, xlab = "Residuals", main = "")
plot(SBWALTEM$Peak,SBWALTEM$resid,xlab="Peak",ylab="Residuals")
plot(SBWALTEM$Plot,SBWALTEM$resid,xlab="Plot",ylab="Residuals")
boxplot(SBWALTEM$resid~SBWALTEM$SBWALT,xlab="Caterpillar",ylab="Residuals")
plot(SBWALTEM$Peak[SBWALTEM$SBWALT=="SBW"],SBWALTEM$resid[SBWALTEM$SBWALT=="SBW"],xlab="Peak",ylab="Residuals",main="SBW")
plot(SBWALTEM$Peak[SBWALTEM$SBWALT=="ALT"],SBWALTEM$resid[SBWALTEM$SBWALT=="ALT"],xlab="Peak",ylab="Residuals",main="ALT")
par(op)

#Obviously violation of independence of x values (time), also possibly violation of homogeneity assumption!#

#2. Repeat step 1 using the gls function from the nlme package using REML estimation but without any variance structures.
SBWALTEMfullgls<-gls(NEMpcap~dummyPeak*SBWALT*Plot,data=SBWALTEM, method="REML")

#3. Choose an appropriate variance structure depending on graphical analysis above.
SBWALTEMfixed<-gls(NEMpcap~dummyPeak*SBWALT*Plot,weights=varFixed(~dummyPeak),data=SBWALTEM, method="REML")
SBWALTEMcomb<-gls(NEMpcap~dummyPeak*SBWALT*Plot,weights=varComb(varIdent(form=~1|SBWALT),varFixed(~dummyPeak)),data=SBWALTEM, method="REML")
SBWALTEMident<-gls(NEMpcap~dummyPeak*SBWALT*Plot,weights=varIdent(form=~1|SBWALT),data=SBWALTEM, method="REML")
SBWALTEMpower<-gls(NEMpcap~dummyPeak*SBWALT*Plot,weights=varPower(form=~dummyPeak|SBWALT),data=SBWALTEM, method="REML")

anova(SBWALTEMfullgls,SBWALTEMfixed,SBWALTEMident,SBWALTEMpower,SBWALTEMcomb)

#From AIC values, without any correlation structure is best (-29.64877 for without corr structure versus -28.26646 for with varident)
###Chosen no correlation structure (residuals looked pretty good from above and lowest AIC)
#autocorrelation
SBWALTEMautocorr<-gls(NEMpcap~dummyPeak*SBWALT*Plot,correlation=corAR1(form=~dummyPeak|Plot/SBWALT),data=SBWALTEM, method="REML")
anova(SBWALTEMfullgls,SBWALTEMautocorr)
#ACF
acf(residuals(SBWALTEMfullgls, type = "normalized"), na.action = na.pass,
    main = "Auto-correlation plot for residuals")
acf(residuals(SBWALTEMautocorr, type = "normalized"), na.action = na.pass,
    main = "Auto-correlation plot for residuals")
E2b <- resid(SBWALTEMautocorr, type="normalized")
coplot(E2b ~ dummyPeak | SBWALT,
       ylab = "Normalised residuals", data = SBWALTEM)

#decided to keep corAR1 because does improve model (lower AIC)

#4. Using the gls model with the selected variance structure from 3., try to find the optimal random structure.
#Not doing because Plot is a fixed effect (using reasoning from http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html)

#5. Compare the new gls model with the earlier results using AIC.

#6.Find the optimal ﬁxed component using the likelihood ratio test where ML estimation is used.
SBWALTEMfullml<-gls(NEMpcap~dummyPeak*SBWALT*Plot,correlation=corAR1(form=~dummyPeak|Plot/SBWALT),data=SBWALTEM, method="ML")

SBWALTEMdrop3interactionml<-update(SBWALTEMfullml,.~.-dummyPeak:SBWALT:Plot)
anova(SBWALTEMfullml,SBWALTEMdrop3interactionml)

# do not drop three way interaction because significant

#7. Refit the model found in step 6. with REML estimation
SBWALTEMfinalmodel<-gls(NEMpcap~dummyPeak*SBWALT*Plot,correlation=corAR1(form=~dummyPeak|Plot/SBWALT),data=SBWALTEM, method="REML") 

summary(SBWALTEMfinalmodel)



#Create graph of number of parasitoid emergences from budworm and other caterpillars over time,
SBWALTEMfinalmodelGLS<-Gls(NEMpcap~dummyPeak*SBWALT*Plot,correlation=corAR1(form=~dummyPeak|Plot/SBWALT),data=SBWALTEM, method="REML") #using the function Gls because allows Predict function needed to plot trendline
NEMpredict<-Predict(SBWALTEMfinalmodelGLS,SBWALT=c("ALT","SBW"),dummyPeak=c(1:14),Plot=c(1:3))%>%
  mutate(label=as.factor(ifelse(SBWALT=="SBW", "budworm","other caterpillars")),Plot=as.factor(Plot))

SBWALTEMpredict<-left_join(SBWALTEM,NEMpredict,by=c("SBWALT", "label","Plot","dummyPeak"))

SBWALTEMplot<-SBWALTEMpredict%>%#PLACEHOLDER FIX LINES SO NOT FOR ALL YEARS BUT FOR EACH PLOT
  ggplot+
  geom_vline(xintercept=0)+
  geom_line(aes(Peak,yhat, linetype=Plot), size=2)+
  geom_point(data=SBWALTEM,aes(Peak,NEMpcap,shape=Plot),size=5)+
  facet_wrap(~label,scales="free_y",ncol=1)+
  theme(strip.background=element_blank(),strip.text.x=element_text(size=28), 
        axis.title.y=element_text(hjust=0.5, vjust=1.5))+
  ylab("Per capita \nparasitoid emergences")+xlab("Years before/after peak")+ylim(c(0,NA))+scale_shape_manual(values=c(16,17,15))+scale_linetype_manual(values=c("dotted","dashed","solid"))

ggsave("figs/SBWALTEMplot.pdf",plot=SBWALTEMplot,width=8,height=7) #figure 2b

# Aggregate response-----------------------------------------------------------
numinteractionsspecies<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  group_by(SpecID.high)%>%
  summarise(numinter=length(SpecID.high))%>%
  arrange(desc(numinter)) # Produces a dataframe of the total number of emergences with budworm and other caterpillars for each parasitoid taxon. Then arranged descending from most emergences to lowest emergences.

paralist<-numinteractionsspecies$SpecID.high

drop_top3_sp <- function(x){
communitydietswitchlog_dropsp<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  filter(!SpecID.high %in% x)%>%
  group_by(Peak, Plot)%>%
  summarise(TotSBWEmerge=length(SpecID.high[SBWALT=="SBW"]),TotALTEmerge=length(SpecID.high[SBWALT=="ALT"]))%>%
  full_join(ratioSBWbypeakplot,by=c("Peak","Plot"))%>%
  mutate(ratioSBWEmerge=TotSBWEmerge/TotALTEmerge,logratioSBWEmerge=log10(ratioSBWEmerge)) %>%
  filter(TotSBW!=0,TotALT!=0,TotSBWEmerge!=0,TotALTEmerge!=0)%>%
  mutate(dummyPeak=Peak+4)

SBWALTdietswitchreml_dropsp<-gls(logratioSBWEmerge~logratioSBWpeakplot*Plot,data=communitydietswitchlog_dropsp,method="REML")

avinterceptallplotsDS_dropsp<-(((3*coef(SBWALTdietswitchreml_dropsp)[1])+coef(SBWALTdietswitchreml_dropsp)[3]+coef(SBWALTdietswitchreml_dropsp)[4])/3)

avSEinterceptallplotsDS_dropsp <- ((summary(SBWALTdietswitchreml_dropsp)$tTable[1,2]+summary(SBWALTdietswitchreml_dropsp)$tTable[3,2]+summary(SBWALTdietswitchreml_dropsp)$tTable[4,2])/3)

avslopeallplotsDS_dropsp<-(((3*coef(SBWALTdietswitchreml_dropsp)[2])+coef(SBWALTdietswitchreml_dropsp)[5]+coef(SBWALTdietswitchreml_dropsp)[6])/3)

avSEslopeallplotsDS_dropsp<-((summary(SBWALTdietswitchreml_dropsp)$tTable[2,2]+summary(SBWALTdietswitchreml_dropsp)$tTable[5,2]+summary(SBWALTdietswitchreml_dropsp)$tTable[6,2])/3)


bt<- (avslopeallplotsDS_dropsp-avslopeallplotsDS)/avSEslopeallplotsDS_dropsp

it<-(avinterceptallplotsDS_dropsp-avinterceptallplotsDS)/avSEinterceptallplotsDS_dropsp

ttable<-tibble(
  speciesdrop=paste(x,collapse=" "),
  origbeta=avslopeallplotsDS,
  beta=avslopeallplotsDS_dropsp,
  diffbeta=1-beta/origbeta,
  binterval=avSEslopeallplotsDS_dropsp*1.96,
  btstat=bt,
  bdf=SBWALTdietswitchreml_dropsp$dims[[1]]-SBWALTdietswitchreml_dropsp$dims[[2]],
  bp=pt(abs(bt),df=SBWALTdietswitchreml_dropsp$dims[[1]]-SBWALTdietswitchreml_dropsp$dims[[2]],lower.tail = FALSE),
  bp2=2*bp,
  origint=avinterceptallplotsDS,
  int=avinterceptallplotsDS_dropsp,
  diffint=1-int/origint,
  iinterval=avSEinterceptallplotsDS_dropsp*1.96,
  itstat=it,
  idf=SBWALTdietswitchreml_dropsp$dims[[1]]-SBWALTdietswitchreml_dropsp$dims[[2]],
  ip=pt(abs(it),df=SBWALTdietswitchreml_dropsp$dims[[1]]-SBWALTdietswitchreml_dropsp$dims[[2]],lower.tail = FALSE),
  ip2=2*ip
)
return(ttable)
}

table<-bind_rows(drop_top3_sp("p01"), drop_top3_sp(c("p01","p02")),drop_top3_sp(c("p01","p02","p03"))) %>% #drop the most abundant parasitoid taxon (p01), then drop the two most abundant taxa, then drop the three most abundant taxa. compare coefficient values (average of all three plots), run one sample ttests and combine into one table.
  select(-c(origbeta,origint,diffbeta,diffint,bdf,idf,bp,ip)) %>%
  data.frame()

stargazer(table, type="latex",summary=FALSE,rownames=FALSE)

#graph of per capita emergences from budworm and other caterpillars compared to ratio of budworm and other caterpillars
#Apanteles fumiferanae
communitydietswitchlog_p01<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  filter(SpecID.high=="p01")%>%
  group_by(Peak, Plot)%>%
  summarise(TotSBWEmerge=length(SpecID.high[SBWALT=="SBW"]),TotALTEmerge=length(SpecID.high[SBWALT=="ALT"]))%>%
  full_join(ratioSBWbypeakplot,by=c("Peak","Plot"))%>%
  mutate(ratioSBWEmerge=TotSBWEmerge/TotALTEmerge,logratioSBWEmerge=log10(ratioSBWEmerge)) %>%
  filter(TotSBW!=0,TotALT!=0,TotSBWEmerge!=0,TotALTEmerge!=0) %>%
ggplot(aes(logratioSBWpeakplot, logratioSBWEmerge))+
  geom_abline(intercept=0,slope=1,linetype="twodash")+
  geom_abline(intercept=0,slope=0)+
  geom_vline(xintercept=0)+
  geom_point(aes(shape=Plot),size=5)+
  theme(axis.title.y=element_text(hjust=0.5, vjust=1.5), 
        legend.text=element_text(size=14),legend.justification=c(1,0), legend.position=c(0.98,1),legend.box = "horizontal")+
  coord_fixed(ratio = 1)+
  scale_shape_manual(values=c(16,17,15))+scale_linetype_manual(values=c("dotted","dashed","solid"))+
  ylab("Log10 emergence \nbudworm:other caterpillars")+xlab("Log10 abundance \nbudworm:other caterpillars")

ggsave("figs/communitydietswitchlogratio_p01.pdf", plot=communitydietswitchlog_p01, width=8,height=8) #Supplementary Figure

#Glypta fumiferanae
communitydietswitchlog_p02<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  filter(SpecID.high=="p02")%>%
  group_by(Peak, Plot)%>%
  summarise(TotSBWEmerge=length(SpecID.high[SBWALT=="SBW"]),TotALTEmerge=length(SpecID.high[SBWALT=="ALT"]))%>%
  full_join(ratioSBWbypeakplot,by=c("Peak","Plot"))%>%
  mutate(ratioSBWEmerge=TotSBWEmerge/TotALTEmerge,logratioSBWEmerge=log10(ratioSBWEmerge)) %>%
  filter(TotSBW!=0,TotALT!=0,TotSBWEmerge!=0,TotALTEmerge!=0) %>%
  ggplot(aes(logratioSBWpeakplot, logratioSBWEmerge))+
  geom_abline(intercept=0,slope=1,linetype="twodash")+
  geom_abline(intercept=0,slope=0)+
  geom_vline(xintercept=0)+
  geom_point(aes(shape=Plot),size=5)+
  theme(axis.title.y=element_text(hjust=0.5, vjust=1.5), 
        legend.text=element_text(size=14),legend.justification=c(1,0), legend.position=c(0.2,0.5),legend.box = "horizontal")+
  coord_fixed(ratio = 1)+
  scale_shape_manual(values=c(16,17,15))+scale_linetype_manual(values=c("dotted","dashed","solid"))+
  ylab("Log10 emergence \nbudworm:other caterpillars")+xlab("Log10 abundance \nbudworm:other caterpillars")

ggsave("figs/communitydietswitchlogratio_p02.pdf", plot=communitydietswitchlog_p02, width=8,height=8) #Supplementary Figure

#Meteorus trachynotus
communitydietswitchlog_p03<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  filter(SpecID.high=="p03")%>%
  group_by(Peak, Plot)%>%
  summarise(TotSBWEmerge=length(SpecID.high[SBWALT=="SBW"]),TotALTEmerge=length(SpecID.high[SBWALT=="ALT"]))%>%
  full_join(ratioSBWbypeakplot,by=c("Peak","Plot"))%>%
  mutate(ratioSBWEmerge=TotSBWEmerge/TotALTEmerge,logratioSBWEmerge=log10(ratioSBWEmerge)) %>%
  filter(TotSBW!=0,TotALT!=0,TotSBWEmerge!=0,TotALTEmerge!=0) %>%
  ggplot(aes(logratioSBWpeakplot, logratioSBWEmerge))+
  geom_abline(intercept=0,slope=1,linetype="twodash")+
  geom_abline(intercept=0,slope=0)+
  geom_vline(xintercept=0)+
  geom_point(aes(shape=Plot),size=5)+
  theme(axis.title.y=element_text(hjust=0.5, vjust=1.5), 
        legend.text=element_text(size=14),legend.justification=c(1,0), legend.position=c(0.1,0.6),legend.box = "horizontal")+
  coord_fixed(ratio = 1)+
  scale_shape_manual(values=c(16,17,15))+scale_linetype_manual(values=c("dotted","dashed","solid"))+
  ylab("Log10 emergence \nbudworm:other caterpillars")+xlab("Log10 abundance \nbudworm:other caterpillars")

ggsave("figs/communitydietswitchlogratio_p03.pdf", plot=communitydietswitchlog_p03, width=8,height=8) #Supplementary Figure

#nMDS analysis of parasitoid community over time
turnoverprep<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  group_by(Peak,Plot)%>%
  summarise(Totemerge=length(SpecID.high))#creates a dataframe of the total emergences of all parasitoid taxa combined for each peak and plot (used to standardise the number of emergences per taxon by the total number of emergences because different years had different total number of emergences)

turnover<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  group_by(Peak,Plot, SpecID.high)%>%
  summarise(nemerge=length(SpecID.high))%>%
  ungroup()%>%
  left_join(turnoverprep,by=c("Peak","Plot"))%>%
  mutate(nemergestand=nemerge/Totemerge, bpae=case_when(
    Peak %in% c(-3,-2) ~ "b",
    Peak %in% c(-1,0,1) ~ "p",
    Peak %in% c(2,3) ~ "a",
    Peak %in% c(4,5,6,7,8,9,10) ~ "en"
  ))%>%
  mutate(bpae=as.factor(bpae))
  
turnovercommat<-turnover%>%
  select(-c(nemerge,Totemerge))%>%
  spread(SpecID.high,nemergestand)%>%
  ungroup()#creates a community matrix (but still including the factors) of the dataframe turnover.

factors<-tibble(
  Plot=turnovercommat$Plot,
  Peak=turnovercommat$Peak,
  bpae=turnovercommat$bpae
)#creates a dataframe of the different factors that will be used in the nMDS and the PERMANOVA

turnovercommat%<>%select(-c(Peak,Plot,bpae))#removes the factors from the dataframe to create a true community matrix

turnovercommat[is.na(turnovercommat)]<-0#replaces any NA values with 0. these NA values were created when the long format data was spread into a wide format dataframe and where taxa were not found for a certain Peak or Plot.

dimcheckMDS(turnovercommat, distance = "bray", k = 6, trymax = 30, 
        autotransform=TRUE) #after exmining screeplot from 6 dimensions to 1 dimension, chosen to use 2 dimensions because additional dimensions provide small reductions in stress

turnoverNMDS <- metaMDS(turnovercommat, distance = "bray", k = 2, trymax = 50, 
                        autotransform=TRUE)

stressplot(turnoverNMDS)  ##to test if distances as represented in ordination are correlated with actual distances

stress10<-c(0.07992804 ,0.07981951,0.08191709,0.1194109,0.08585157, 0.08352331,0.0821567,0.1255087,0.07992697,0.1215399)
sd(stress10)

#Create nMDS plot
turnoverdatascores<-as.data.frame(scores(turnoverNMDS))
turnoverdatascores$Plot<-factors$Plot
turnoverdatascores$bpae<-factors$bpae

# function for creating ellipses in nmds plot
#adapted from  http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
veganCovEllipsenew <- function (x, scale = 1, npoints = 100) 
{
  cov <- cov.wt(cbind(x$NMDS1,x$NMDS2),wt=rep(1/length(x$NMDS1),length(x$NMDS1)))$cov 
  center <- c(mean(x$NMDS1),mean(x$NMDS2))
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  bpaeunique<-unique(x$bpae)
  tibble(
    NMDS1=t(center + scale * t(Circle %*% chol(cov)))[,1],
    NMDS2=t(center + scale * t(Circle %*% chol(cov)))[,2],
    bpae=bpaeunique
  )
}


df_ell.bpae<-turnoverdatascores%>%
  split(.$bpae)%>%
  map(veganCovEllipsenew)%>%
  bind_rows()

turnoverNMDSplot<-ggplot(turnoverdatascores)+
  geom_point(aes(NMDS1,NMDS2, shape=Plot, colour=bpae),size=5)+
  geom_path(data=df_ell.bpae, aes(x=NMDS1, y=NMDS2, colour=bpae), size=1.5)+
  theme(axis.title.y=element_text(hjust=0.5, vjust=1.5), 
        legend.text=element_text(size=14),
        legend.title=element_text(size=15),
        legend.justification=c(1,0.99), 
        legend.position=c(1,0.99),
        legend.box = "horizontal")+
  coord_fixed(ratio = 1)+scale_color_viridis(name="Peak", breaks=c("b","p","a", "en"),labels=c("before","during","after", "endemic"), alpha = 1, begin = 0, end = 1, direction = 1, discrete = TRUE, option = "D")+scale_shape_manual(values=c(16,17,15))+scale_x_continuous(breaks=c(-1,0,1))

ggsave("figs/nmdsturnover.pdf",plot=turnoverNMDSplot,width=10,height=10) #Figure 3


# Testing the separation of the parasitoid communities between before, during, and after the budworm peak, using the adonis function for a permanova. - Guide to formulating the permutation structure and permutation test from the following webpage http://thebiobucket.blogspot.ca/2011/04/repeat-measure-adonis-lately-i-had-to.html#more

fit <- adonis(vegdist(turnovercommat)~bpae, factors, permutations=1)

# number of perms
B <- 999

### setting up a dataframe which will be populated by random F Model values:
pop <- rep(NA, B + 1)

# the first entry of the dataframe pop will be the true F. model:
pop[1] <- fit$aov.tab[1, 4]

# set up a the permutation structure where permutations occur withing each plot and where the order of years are maintained for the permutations:
ctrl <- how(within = Within(type = "series",mirror = FALSE), plots = Plots(strata = factors$Plot))

# Number of observations:
nobs <- nrow(turnovercommat)

# loop to populate the dataframe pop with random F model values:
set.seed(123)
for(i in 2:(B+1)){
  idx <- shuffle(nobs, control = ctrl)
  fit.rand <- adonis(vegdist(turnovercommat) ~ bpae[idx], factors,permutations = 1)
  pop[i] <- fit.rand$aov.tab[1, 4]
}

# Return the p-value of the full permanova:
print(pval <- sum(pop >= pop[1]) / (B + 1))


# Topology and per capita interaction strengths---------


speciesnamesprepbi<-speciesnamesshare%>%
  rename(dataID=SpecID.high,SpecID.high=LatinName)%>%
  right_join(speciesnamesprep,by="SpecID.high")%>%
  select(-SpecID.high)%>%
  rename(SpecID.high=dataID,SpecID.highnum=numID) # creates a dataframe that includes the parasitoid code from data collection and the corresponding number for the figure.

bipartitegraphSBWALT<-interaction%>%
  mutate(Peak=as.character(Peak))%>%
  mutate(Peak=gsub("-","minus",Peak))%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  left_join(speciesnamesprepbi,by="SpecID.high")%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  select(SBWALT,SpecID.highnum,Peak)%>%
  frame2webs(varnames=c("SBWALT","SpecID.highnum","Peak")) # creates a dataframe that counts the number of interactions between each parasitoid taxon and budworm or other caterpillars for each year.

## Create bipartite graphs to identify the changing interaction strength distributions between parasitoids, budworm, and other caterpillars

# Calculation of percentage of emergences from budworm and other caterpillas for each bipartite graph (in the main text)
percentSBWALT_fun<-function(year){
  percent<-tibble(
    Total=rowSums(as.matrix(year))[1]+ifelse(is.na(rowSums(as.matrix(year))[2]),0,rowSums(as.matrix(year))[2]),
    SBWpercent=(ifelse(is.na(rowSums(as.matrix(year))[2]),0,rowSums(as.matrix(year))[2])/Total)*100,
    ALTpercent=(rowSums(as.matrix(year))[1]/Total)*100
  )
  return(percent)
}

percentSBWALT_fun(bipartitegraphSBWALT$minus3)

percentSBWALT_fun(bipartitegraphSBWALT$`0`)

percentSBWALT_fun(bipartitegraphSBWALT$`3`)

percentSBWALT_fun(bipartitegraphSBWALT$`10`)

# final creation of bipartite graphs for the main text (Figure 3).
#minus3
pdf("figs/minus3plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'minus3'),ybig=1.6,low.spacing=0.8,high.spacing=0.1,x.lim=c(0,2.6),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.4,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

#0
pdf("figs/0plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'0'),ybig=1.6,low.xoff = 1.2,low.spacing=2,high.spacing=0.2,x.lim=c(0,5.3),labsize=2.5,high.lab.dis=0.05,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

#3
pdf("figs/3plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'3'),ybig=1.6,low.xoff = 1,low.spacing=1.8,high.spacing=0.2,x.lim=c(0,5.61),labsize=2.5,high.lab.dis=0.05,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

#10
pdf("figs/10plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'10'),ybig=1.6,high.spacing=0.1,x.lim=c(0,1.2),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.1,col.low="#C73D73FF", col.high="#160F3BFF",col.interaction="#FECC8FFF", method='normal')
dev.off()


#calculation of median interactions strength to maximum interaction strength
#don't need to calculate per captita emergence because ratio of median per capita emergence to max per capita emergence would cancel out denominator of number of SBW or ALT - so same as using number of emergences
bipartitegraphSBWALTplotprep<-interaction%>%
  mutate(Peak=as.character(Peak))%>%
  mutate(Peak=gsub("-","minus",Peak))%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  left_join(speciesnamesprepbi,by="SpecID.high")%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  select(SBWALT,SpecID.highnum,Peak,Plot)

plot1bipartite<-bipartitegraphSBWALTplotprep %>%
  filter(Plot==1)%>%
  select(SBWALT,SpecID.highnum,Peak)%>%
  frame2webs(varnames=c("SBWALT","SpecID.highnum","Peak")) # creates a dataframe that counts the number of interactions between each parasitoid taxon and budworm or other caterpillars for each year.


plot2bipartite<-bipartitegraphSBWALTplotprep %>%
  filter(Plot==2)%>%
  select(SBWALT,SpecID.highnum,Peak)%>%
  frame2webs(varnames=c("SBWALT","SpecID.highnum","Peak")) # creates a dataframe that counts the number of interactions between each parasitoid taxon and budworm or other caterpillars for each year.

plot3bipartite<-bipartitegraphSBWALTplotprep %>%
  filter(Plot==3)%>%
  select(SBWALT,SpecID.highnum,Peak)%>%
  frame2webs(varnames=c("SBWALT","SpecID.highnum","Peak")) # creates a dataframe that counts the number of interactions between each parasitoid taxon and budworm or other caterpillars for each year.

intdistplot1<-plot1bipartite%>%
  lapply(function(x) {y<-data.frame(t(x))
  if ("ALT" %in% colnames(y)) 
  {mutate(y,ALT=ifelse(ALT==0,NA,ALT),SBW=ifelse(SBW==0,NA,SBW))%>%
      summarise(medALT=median(ALT, na.rm=TRUE),maxALT=max(ALT,na.rm=TRUE),medSBW=median(SBW, na.rm=TRUE),maxSBW=max(SBW, na.rm=TRUE))%>%
      mutate(medmaxratioALT=medALT/maxALT,medmaxratioSBW=medSBW/maxSBW)}
  else
  {mutate(y,SBW=ifelse(SBW==0,NA,SBW))%>%
      summarise(medSBW=median(SBW, na.rm=TRUE),maxSBW=max(SBW, na.rm=TRUE))%>%
      mutate(medmaxratioSBW=medSBW/maxSBW)}})%>%
  bind_rows()%>%
  mutate(Year=names(plot1bipartite), Plot=1)# creates a dataframe of median to maximum interaction strength (number of emergences is a proxy for interaction strength) for parasitoids emerging from budworm (SBW) and other caterpillars (ALT).

intdistplot2<-plot2bipartite%>%
  lapply(function(x) {y<-data.frame(t(x))
  if ("SBW" %in% colnames(y)) 
  {mutate(y,ALT=ifelse(ALT==0,NA,ALT),SBW=ifelse(SBW==0,NA,SBW))%>%
      summarise(medALT=median(ALT, na.rm=TRUE),maxALT=max(ALT,na.rm=TRUE),medSBW=median(SBW, na.rm=TRUE),maxSBW=max(SBW, na.rm=TRUE))%>%
      mutate(medmaxratioALT=medALT/maxALT,medmaxratioSBW=medSBW/maxSBW)}
  else
  {mutate(y,ALT=ifelse(ALT==0,NA,ALT))%>%
      summarise(medALT=median(ALT, na.rm=TRUE),maxALT=max(ALT,na.rm=TRUE))%>%
      mutate(medmaxratioALT=medALT/maxALT)}})%>%
  bind_rows()%>%
  mutate(Year=names(plot2bipartite),Plot=2)# creates a dataframe of median to maximum interaction strength (number of emergences is a proxy for interaction strength) for parasitoids emerging from budworm (SBW) and other caterpillars (ALT).


intdistplot3<-plot3bipartite%>%
  lapply(function(x) {y<-data.frame(t(x))
  if ("SBW" %in% colnames(y)) 
  {mutate(y,ALT=ifelse(ALT==0,NA,ALT),SBW=ifelse(SBW==0,NA,SBW))%>%
      summarise(medALT=median(ALT, na.rm=TRUE),maxALT=max(ALT,na.rm=TRUE),medSBW=median(SBW, na.rm=TRUE),maxSBW=max(SBW, na.rm=TRUE))%>%
      mutate(medmaxratioALT=medALT/maxALT,medmaxratioSBW=medSBW/maxSBW)}
  else
  {mutate(y,ALT=ifelse(ALT==0,NA,ALT))%>%
      summarise(medALT=median(ALT, na.rm=TRUE),maxALT=max(ALT,na.rm=TRUE))%>%
      mutate(medmaxratioALT=medALT/maxALT)}})%>%
  bind_rows()%>%
  mutate(Year=names(plot3bipartite),Plot=3)# creates a dataframe of median to maximum interaction strength (number of emergences is a proxy for interaction strength) for parasitoids emerging from budworm (SBW) and other caterpillars (ALT).

intdistall <- bind_rows(intdistplot1,intdistplot2,intdistplot3) %>%
  select(medmaxratioALT,medmaxratioSBW,Year, Plot) %>%
  rename(ALT=medmaxratioALT,SBW=medmaxratioSBW) %>%
  gather(SBWALT, medmax, -Year, -Plot) %>%
  mutate(Plot=as.factor(Plot))

intdistall$Year<-gsub("minus","-",intdistall$Year)

intdistplot<-ggplot(intdistall)+
  geom_vline(xintercept=0)+
  geom_point(aes(as.numeric(Year),medmax, colour=SBWALT, shape=Plot), size=5)+
  xlab("Years before/after peak")+ylab("Median:maximum \ninteraction strength")+
  theme(legend.text=element_text(size=14),
        legend.title=element_text(size=15),
        legend.justification=c(0,1),
        legend.position=c(0.01,0.99))+
  scale_color_manual(name="caterpillar type", breaks=c("SBW","ALT"),labels=c("budworm","other caterpillars"), values=c("#C73D73FF","#2D1160FF"))+
    scale_x_continuous(breaks=c(-2,0,2,4,6,8,10))+scale_shape_manual(values=c(16,17,15))

ggsave("figs/intdistplot.pdf",intdistplot,width=8,height=7) #Figure 5

# Supplementary Information -----------------------------------------------

###Linear mixed effects model of log10 ratio per capita parasitoid emergence from budworm and other caterpillars to the log10 ratio of total budworm to other caterpillars. Using the diet switching method outlined by Greenwood, J. J. D., and R. A. Elton. 1979. Analysing experiments on frequency-dependent selection by predators. The Journal of Animal Ecology 48:721.

#All of the statistical analysis below follows the protocol outlined in Zuur, A., E. N. Ieno, N. Walker, A. A. Saveliev, and G. M. Smith. 2009. Mixed effects models and extensions in ecology with r. First edition. Springer-Verlag New York, New York, New York, United States of America.

#Data Exploration
#1.Outliers in the response and explanatory variables.
op <- par(mfrow = c(1, 2), mar = c(3, 3, 3, 1))
dotchart(communitydietswitch$logratioSBWEmergepcap, main = "Log ratio SBW Emerge Per Capita")
dotchart(communitydietswitch$logratioSBWpeakplot, main = "Log Ratio SBW:ALT")
par(op)
#2. Collinearity of the explanatory variables.
#Not doing as do not have multiple continuous explanatory variables

#3. Relationships between the response variable and the explanatory variables.
xyplot(logratioSBWEmerge~logratioSBWpeakplot, data=communitydietswitch, groups=Plot) #not percapita
xyplot(logratioSBWEmergepcap~logratioSBWpeakplot, data=communitydietswitch, groups=Plot) # per capita
xyplot(TotSBWEmerge~Peak, data=communitydietswitch)
xyplot(TotALTEmerge~Peak, data=communitydietswitch)
xyplot(TotSBWEmergepcap~Peak, data=communitydietswitch)
xyplot(TotALTEmergepcap~Peak, data=communitydietswitch)
#Creating the model
#1.Start with a linear regression model that contains as many explanatory variables and their interactions as possible. 
SBWALTdietswitchfulllinear<-lm(logratioSBWEmergepcap~logratioSBWpeakplot*Plot,data=communitydietswitch)
communitydietswitch$resid<-NA
communitydietswitch$resid<-SBWALTdietswitchfulllinear$residuals
#Test normality and homogeniety of residuals
op <- par(mfrow = c(2, 2), mar = c(5, 4, 1, 2))
plot(SBWALTdietswitchfulllinear, add.smooth = FALSE, which = 1,col=communitydietswitch$Plot)
hist(communitydietswitch$resid, xlab = "Residuals", main = "")
plot(communitydietswitch$Peak,communitydietswitch$resid,xlab="Peak",ylab="Residuals",col=as.factor(communitydietswitch$Plot))
plot(communitydietswitch$Plot,communitydietswitch$resid,xlab="Plot",ylab="Residuals")
par(op)

acf(residuals(SBWALTdietswitchfulllinear), na.action = na.pass,
    main = "Auto-correlation plot for residuals")
#Obviously violation of independence of x values (time) and repeated measures of the same plot in different years (although acf does not show much correlation between years)

#2. Repeat step 1 using the gls function from the nlme package using REML estimation but without any variance structures.
SBWALTdietswitchgls<-gls(logratioSBWEmergepcap~logratioSBWpeakplot*Plot,data=communitydietswitch, method="REML")
#3. Choose an appropriate variance structure depending on graphical analysis above.

SBWALTdietswitchautocorr<-gls(logratioSBWEmergepcap~logratioSBWpeakplot*Plot,data=communitydietswitch,correlation=corAR1(form=~dummyPeak|Plot), method="REML")

anova(SBWALTdietswitchgls,SBWALTdietswitchautocorr)#Delta AIC<2

#check assumptions graphically for model with autocorrelation
communitydietswitch$fmnormresid<-NA
communitydietswitch$fmnormresid<-resid(SBWALTdietswitchautocorr, type="normalized")

plot(SBWALTdietswitchautocorr, add.smooth = FALSE)
hist(communitydietswitch$fmnormresid, xlab = "Normalised Residuals", main = "")
plot(communitydietswitch$Peak,communitydietswitch$fmnormresid,xlab="Peak",ylab="Normalised Residuals")
plot(communitydietswitch$Peak,communitydietswitch$fmnormresid,xlab="Peak",ylab="Residuals",col=as.factor(communitydietswitch$Plot))
plot(communitydietswitch$Plot,communitydietswitch$fmnormresid,xlab="Plot",ylab="Residuals")
acf(communitydietswitch$fmnormresid, na.action = na.pass,
    main = "Auto-correlation plot for normalised residuals")

#AIC lowest for model with autocorr but tiny change in ACF (and ACF does not show large improvement, also residuals are less normally distributed). Changes too small so not including corAR1.

#4. Using the gls model with the selected variance structure from 3., try to find the optimal random structure.
#Not doing because Plot is a fixed effect (using reasoning from https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html)

#6.Find the optimal ﬁxed component using the likelihood ratio test where ML estimation is used.
SBWALTdietswitchmlfull<-gls(logratioSBWEmergepcap~logratioSBWpeakplot*Plot,data=communitydietswitch, method="ML")

SBWALTdietswitchmlminusinteraction<-gls(logratioSBWEmergepcap~logratioSBWpeakplot+Plot,data=communitydietswitch, method="ML")
anova(SBWALTdietswitchmlfull,SBWALTdietswitchmlminusinteraction)

# Interaction is significant (p-value = 0.0011, delta AIC = 9.60224). Therefore, keep all fixed effects in the model including the interaction.


#7. Refit the model found in step 6. with REML estimation
SBWALTdietswitchreml<-gls(logratioSBWEmergepcap~logratioSBWpeakplot*Plot,data=communitydietswitch, method="REML")

summary(SBWALTdietswitchreml)

# average intercept for all plots
avinterceptallplotsDS<-(((3*coef(SBWALTdietswitchreml)[1])+coef(SBWALTdietswitchreml)[3]+coef(SBWALTdietswitchreml)[4])/3)

# average standard error of intercepts for all plots
avSEinterceptallplotsDS <- ((summary(SBWALTdietswitchreml)$tTable[1,2]+summary(SBWALTdietswitchreml)$tTable[3,2]+summary(SBWALTdietswitchreml)$tTable[4,2])/3)

#average slope for all plots
avslopeallplotsDS<-(((3*coef(SBWALTdietswitchreml)[2])+coef(SBWALTdietswitchreml)[5]+coef(SBWALTdietswitchreml)[6])/3)

#average standard error of slopes for all plots
avSEslopeallplotsDS<-((summary(SBWALTdietswitchreml)$tTable[2,2]+summary(SBWALTdietswitchreml)$tTable[5,2]+summary(SBWALTdietswitchreml)$tTable[6,2])/3)

#confidence intervals for average slope
avSEslopeallplotsDS*1.96

#Test of whether the average slope of all plots is different from 0 (test of resource tracking)
2*pt(abs((avslopeallplotsDS-0)/avSEslopeallplotsDS),df=SBWALTdietswitchreml$dims[[1]]-SBWALTdietswitchreml$dims[[2]],lower.tail = FALSE) #p<0.05 so b is significantly different from 0

#Create graph of log10 ratio parasitoid emergence from budworm and other caterpillars to the log10 ratio of total budworm to other caterpillars
SBWALTdietswitchremlGLS<-Gls(logratioSBWEmergepcap~logratioSBWpeakplot*Plot,data=communitydietswitch, method="REML") #using the function Gls because allows Predict function needed to plot trendline

dietswitchpredict<-Predict(SBWALTdietswitchremlGLS, logratioSBWpeakplot=seq(from=-1.5,to=2.1,by=0.1),Plot=seq(from=1,to=3,by=1))%>%
  mutate(Plot=as.factor(Plot))

comswitchlogratiopcapplot<-ggplot(communitydietswitch,aes(logratioSBWpeakplot, logratioSBWEmergepcap))+
  #geom_abline(intercept=0,slope=1,linetype="twodash")+
  geom_abline(intercept=0,slope=0)+
  geom_vline(xintercept=0)+
  geom_line(data=dietswitchpredict,aes(logratioSBWpeakplot,yhat,linetype=Plot),size=2)+
  geom_point(aes(shape=Plot),size=5)+
  theme(axis.title.y=element_text(hjust=0.5, vjust=1.5), 
        legend.text=element_text(size=14),legend.justification=c(1,0), legend.position=c(0.98,1),legend.box = "horizontal")+
  coord_fixed(ratio = 1)+
  #scale_color_viridis(name="Peak", breaks=c("b","p","a"),labels=c("before","during","after"), alpha = 1, begin = 0, end = 1, direction = 1, discrete = TRUE, option = "D")+
  scale_shape_manual(values=c(16,17,15),guide=FALSE)+scale_linetype_manual(values=c("dotted","dashed","solid"))+
  ylab("Log10 per capita emergence \nbudworm:other caterpillars")+xlab("Log10 abundance \nbudworm:other caterpillars")

ggsave("figs/communitydietswitchlogratiopcap.pdf", plot=comswitchlogratioplotpcap, width=8,height=8) #Figure 1a



#Creation of bipartite graphs for extra years (Figures 6 and 7)
# Calculation of percentage of emergences from budworm and other caterpillas for each bipartite graph (in the supplementary information)
percentSBWALT_fun(bipartitegraphSBWALT$minus2)

percentSBWALT_fun(bipartitegraphSBWALT$minus1)

percentSBWALT_fun(bipartitegraphSBWALT$`1`)

percentSBWALT_fun(bipartitegraphSBWALT$`2`)

percentSBWALT_fun(bipartitegraphSBWALT$`4`)

percentSBWALT_fun(bipartitegraphSBWALT$`5`)

percentSBWALT_fun(bipartitegraphSBWALT$`6`)

percentSBWALT_fun(bipartitegraphSBWALT$`7`)

percentSBWALT_fun(bipartitegraphSBWALT$`8`)

percentSBWALT_fun(bipartitegraphSBWALT$`9`)


# Final creation of biparatite graphs for Supplementary Information
pdf("figs/minus2plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'minus2'),ybig=1.5,low.spacing=0.8,high.spacing=0.1,x.lim=c(0,3.1),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.7,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/minus1plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'minus1'),ybig=1.5,low.spacing=1,high.spacing=0.1,x.lim=c(0,3.4),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.8,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/1plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'1'),ybig=1.5,low.spacing=1,high.spacing=0.1,x.lim=c(0,3.5),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.8,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/2plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'2'),ybig=1.5,low.spacing=0.8,high.spacing=0.1,x.lim=c(0,2.8),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.5,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/4plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'4'),ybig=1.5,low.spacing=0.7,high.spacing=0.1,x.lim=c(0,2.5),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.4,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/5plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'5'),ybig=1.5,low.spacing=0.5,high.spacing=0.1,x.lim=c(0,2),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.2,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/6plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'6'),ybig=1.5,low.spacing=0.4,high.spacing=0.1,x.lim=c(0,1.9),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.1,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/7plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'7'),ybig=1.5,low.spacing=0.5,high.spacing=0.1,x.lim=c(0,2.05),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.1,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/8plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'8'),ybig=1.5,high.spacing=0.1,x.lim=c(0,1.35),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.25,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()


pdf("figs/9plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'9'),ybig=1.6,high.spacing=0.1,x.lim=c(0,1.16),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.05,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

# Useful Extra Analyses ---------------------------------------------------
#Proportion of parasitism of budworm by parasitoid taxa that emerge from both budworm and other caterpillars
allsharedprop<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(!is.na(SpecID.high))%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  group_by(SBWALT)%>%
  summarise(NSha=length(SpecID.low[SpecID.high %in% parashare]),NUnSha=length(SpecID.low[!SpecID.high %in% parashare]),total=length(SpecID.low))%>%
  mutate(PropSh=NSha/total)#90% of all parasitism emergences from budworm were done by taxa that share budworm and other caterpillars - 84% of all parasitism emergences from other caterpillars were done by shared parasitoids

#Chao analysis of interactions - sampling effort
#asking whether sampled enough parasitoids interactions that share sbw and other caterpillars (but dont care about parasitoids that do share)
#probably best to compare chao to percentage of interactions with SBW and ALT with parasitoids that attack both (this dataset)

chaocalc <- interaction %>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  select(SampleID,SpecID.low,SpecID.high)%>%
unite(interaction, SpecID.low,SpecID.high, sep="_" )%>%
  group_by(SampleID,interaction)%>%
  summarise(numinteraction=length(interaction))

SBWALTcode<-interaction %>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  select(SampleID,SpecID.low,SpecID.high)%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT")) %>%
  unite(interaction, SBWALT,SpecID.high, sep="_" )%>%
  group_by(SampleID,interaction)%>%
  summarise(numinteraction=length(interaction))

SBWchaocode<-SBWALTcode%>%
  filter(grepl("SBW",interaction))%>%
  spread(interaction, numinteraction, fill=0)%>%
  ungroup%>%
 select(-SampleID)

ALTchaocode<-SBWALTcode%>%
  filter(grepl("ALT",interaction)) %>%
  spread(interaction, numinteraction, fill=0) %>%
  ungroup%>%
  select(-SampleID)
  

specpool(SBWchaocode)
specpool(ALTchaocode)


#chao on number of other caterpillar species on balsam fir
ALTspecchao<-interaction %>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  select(SampleID,SpecID.low)%>%
  filter(!grepl("h01",SpecID.low))%>%
  group_by(SampleID,SpecID.low)%>%
  summarise(numALTspecies=length(SpecID.low)) %>%
  spread(SpecID.low, numALTspecies, fill=0)%>%
  ungroup%>%
  select(-SampleID)
  
specpool(ALTspecchao)

numinteractionsspecies<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  group_by(SpecID.high)%>%
  summarise(numinter=length(SpecID.high))%>%
  arrange(desc(numinter))

library(scales)
show_col(viridis_pal(option="A")(25))

removedpara <- interactionimport%>%
  mutate(Plot=as.factor(Plot))%>%
  filter(!Plot==4)%>%
  filter(CrownLevel=="mid") %>%
  filter(!SpecID.low=="h11")%>% #remove sawfly from dataset
  filter(!SpecID.low=="h00") #remove undetermined herbivores from dataset (so similar to Eveleigh 2007 SI Materials and Methods)

length(removedpara$SpecID.high) #140452

removedpara %<>% filter(SpecID.high %in% c("pA3", "p40","p05", "p07","p04", "p51"))#remove all unresolved taxa.  #remove diadegma sp. that has freq.low of 2 and uncertain as to the species (either D. pulicalivariae-p94 or Diadegma sp1-p95)

length(removedpara$SpecID.high) #2225

2225/140452 #0.01584171

# find out parasitism rate per year for SBW and ALT
pararate <- interaction %>%
  filter(grepl("h",SpecID.low))%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT")) %>%
  group_by(SBWALT, Peak, Plot) %>%
  summarize(numpara = length(), numtotal = length(SpecID.low))

## check to see how many parasitoids emerge from a single caterpillar (either SBW or ALT)

paraemerge <- interaction %>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  group_by(SBWALT, VialUID) %>%
  summarise(paracount = length(SpecID.high))
  

## Calculation of average number of other caterpillars collected
ALTav <- interaction %>%
    filter(grepl("h",SpecID.low))%>%
    mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
    filter(SBWALT=="ALT") %>%
    group_by(Year) %>%
    summarise(AvALT=length(SBWALT))
mean(ALTav$AvALT)
    
  
  # Notes -------------------------------------------------------------------

#Wintemia fumiferana (Smidtia fumiferana), a major parasitoid of budworm was not included in the shared parasitoid list of parasitoids because found to interact with undetermined herbivores once and with no other other caterpillars (but found to interact with budworm  1241 times) . undetermined herbivores were removed from dataset to follow Eveleigh et al. 2007 SI Methods.


