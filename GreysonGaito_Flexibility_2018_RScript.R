#############################################################################
#            Flexibility of the spruce budworm – parasitoid food web 
#                   on balsam ﬁr
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

# Data Preparation -----
##Input data as data.table
interactionimport<-read_csv("data/InteractFinal_long_2014Nov10.csv") #can use read_csv and will change duo to NA but still 300 parsing failures in other columns (not important columns for this analysis)
speciesnamesimport<-read_csv("data/TotalSpeciesTable_2014Nov10.csv")

##Data cleaning
interaction<-interactionimport%>%
  mutate(Plot=as.factor(Plot))%>%
  filter(!Plot==4)%>%
  filter(!SpecID.low=="h11")%>% #remove sawfly from dataset
  filter(!SpecID.low=="h00")%>% #remove undetermined herbivores from dataset (so similar to Eveleigh 2007 SI Materials and Methods)
  mutate(SpecID.high=ifelse(SpecID.high %in% c("pA3","p40"),"p05",SpecID.high))%>%#lumping unidentified ichnuemonids, unidentified braconids and unidentified ichnuemonoidea together into one group
  mutate(Plot=droplevels(Plot))%>%
  mutate(Peak=ifelse(Plot=="3",Year-6,Year))%>% 
  mutate(Peak=Peak-85)%>% #plot 1 and 2 peaked in 85 and plot 3 peaked in 91. Therefore created new variable Peak to standardize the different plots by when they peak
  filter(CrownLevel=="mid") #using mid crown because royama 2017 used mid crown and because eveleigh and johns 2014 found that using mid crown samples are the best proxies for estimating spruce budworm and parasitoid densities

parashareprep<-interaction%>%
  filter(!grepl("p",SpecID.low))%>%
  filter(SpecID.high!="p51")%>% #remove diadegma sp. that has freq.low of 2 and uncertain as to the species (either D. pulicalivariae-p94 or Diadegma sp1-p95)
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
  bind_rows #produces a list of unique parasitoid taxa that attacked spruce budworm and other caterpillars

parashare<-parashareprep$UnSP[duplicated(parashareprep$UnSP)] #produces a list of the parasitoids that attack both spruce budworm and other caterpillars

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
speciesnamesshare$LatinName<-gsub("undetermined","Und", speciesnamesshare$LatinName)
speciesnamesshare$LatinName<-gsub("Und Ichneumonidae","Und Ichneumonoidea", speciesnamesshare$LatinName)
speciesnamesshare$LatinName[speciesnamesshare$LatinName=="icho"]<-"Und Ichneumonidae sp1"
speciesnamesshare$LatinName<-gsub("microgaster","Microgaster", speciesnamesshare$LatinName)
speciesnamesshare$LatinName<-gsub("Glypta sp1","Glypta sp.", speciesnamesshare$LatinName)
speciesnamesshare$LatinName<-gsub("parasite","parasitoids", speciesnamesshare$LatinName)
speciesnamesshare$LatinName<-gsub("Orgilis","Orgilus", speciesnamesshare$LatinName)
speciesnamesshare$LatinName<-gsub("Synetaeris","Tranosema",speciesnamesshare$LatinName)
speciesnamesshare$LatinName<-gsub("microgaster","Microgaster",speciesnamesshare$LatinName)
speciesnamesshare$LatinName<-gsub("cf. ","",speciesnamesshare$LatinName)
speciesnamesshare$LatinName<-gsub("sp1","sp.",speciesnamesshare$LatinName)

#From speciesnamesshare<- to here, this section creates a list of the latin names of all the parasitoids that attack both spruce budworm and other caterpillars and corrects the latin names if necessary.

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
  mutate(SpecID.high=as.factor(SpecID.high)) #This section produces a dataset of the Proportion of emergences from spruce budworm for each parasitoid taxa. also gives the proportion standard error and error bars for the standard error.

proplevels<-levels(reorder(propSBWALTbyspeciesdata$SpecID.high,propSBWALTbyspeciesdata$PropSBW)) #produces a list where the parasitoid taxa are sorted by the proportion of emergences from spruce budworm.

speciesnamesprep<-tibble(
  SpecID.high=as.factor(rev(proplevels)),
  numID=seq(1,34,by=1)
) #produces a dataset where each parasitoid taxa has a number. used in producing the bipartite graphs.


#The flexible parasitism hypothesis -----

#community level plotting of diet switching using log10 ratios of emergences from spruce budworm and other caterpillars and log10 ratios of total abundances of spruce budworm and other caterpillars

ratioSBWbypeakplot<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  group_by(Peak, Plot)%>%
  summarise(TotSBW=length(SBWALT[SBWALT=="SBW"]),TotALT=length(SBWALT[SBWALT=="ALT"]))%>%
  mutate(ratioSBWpeakplot=TotSBW/TotALT, logratioSBWpeakplot=log10(ratioSBWpeakplot)) #produces a dataset of the ratio (and log 10 ratio) of the total number of spruce budworm and other caterpillars sampled on balsam fir

ratioSBWparacommunity<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  group_by(Peak, Plot)%>%
  summarise(TotSBWEmerge=length(SpecID.high[SBWALT=="SBW"]),TotALTEmerge=length(SpecID.high[SBWALT=="ALT"]))%>%
  mutate(ratioSBWEmerge=TotSBWEmerge/TotALTEmerge,logratioSBWEmerge=log10(ratioSBWEmerge)) #produces a data set of the ratio (and log 10 ratio) of the total number of parasitoid emergences (all species combined) from spruce budworm and other caterpillars sampled on balsam fir.

communitydietswitch<-full_join(ratioSBWparacommunity,ratioSBWbypeakplot,by=c("Peak","Plot"))%>%
  filter(TotSBW!=0,TotALT!=0,TotSBWEmerge!=0,TotALTEmerge!=0)%>% ###CHECK whether should be or for logical question removed years/plots where no SBW collected/attacked and removed years/plots when no ALT collected/attacked. Removed four data points (for a total of 21 data points).
  ungroup()%>%
  mutate(dummyPeak=Peak+4, Plot=as.factor(Plot))%>% #creates a dummy variable because function gls can not use negative numbers.
  mutate(bpa=case_when(
    Peak %in% c(-3,-2)~"b",
    Peak %in% c(-1,0,1)~"p",
    Peak %in% c(2,3,4,5,6,7,8,9,10)~"a"
  ))%>%
  mutate(bpa=as.factor(bpa))

###Linear mixed effects model of log10 ratio parasitoid emergence from spruce budworm and other caterpillars to the log10 ratio of total spruce budworm to other caterpillars. Using the diet switching method outlined by Greenwood, J. J. D., and R. A. Elton. 1979. Analysing experiments on frequency-dependent selection by predators. The Journal of Animal Ecology 48:721.
 
#All of the statistical analysis below follows the protocol outlined in Zuur, A., E. N. Ieno, N. Walker, A. A. Saveliev, and G. M. Smith. 2009. Mixed effects models and extensions in ecology with r. First edition. Springer-Verlag New York, New York, New York, United States of America.

#Data Exploration
#1.Outliers in the response and explanatory variables.
op <- par(mfrow = c(1, 2), mar = c(3, 3, 3, 1))
dotchart(communitydietswitch$logratioSBWEmerge, main = "Log ratio SBW Emerge")
dotchart(communitydietswitch$logratioSBWpeakplot, main = "Log Ratio SBW:ALT")
par(op)
#2. Collinearity of the explanatory variables.
#Not doing as do not have multiple continuous explanatory variables

#3. Relationships between the response variable and the explanatory variables.
xyplot(logratioSBWEmerge~logratioSBWpeakplot, data=communitydietswitch, groups=Plot)

#Creating the model
#1.Start with a linear regression model that contains as many explanatory variables and their interactions as possible. 
SBWALTdietswitchfulllinear<-lm(logratioSBWEmerge~logratioSBWpeakplot*Plot,data=communitydietswitch)
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
SBWALTdietswitchgls<-gls(logratioSBWEmerge~logratioSBWpeakplot*Plot,data=communitydietswitch, method="REML")
#3. Choose an appropriate variance structure depending on graphical analysis above.

SBWALTdietswitchautocorr<-gls(logratioSBWEmerge~logratioSBWpeakplot*Plot,data=communitydietswitch,correlation=corAR1(form=~dummyPeak|Plot), method="REML")

anova(SBWALTdietswitchgls,SBWALTdietswitchautocorr)#Delta AIC<1

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

#AIC lowest for model with autocorr (and ACF shows small improvement). However, changes too small so not including corAR1.

#4. Using the gls model with the selected variance structure from 3., try to find the optimal random structure.
#Not doing because Plot is a fixed effect (using reasoning from https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html)

#6.Find the optimal ﬁxed component using the likelihood ratio test where ML estimation is used.
SBWALTdietswitchmlfull<-gls(logratioSBWEmerge~logratioSBWpeakplot*Plot,data=communitydietswitch, method="ML")

SBWALTdietswitchmlminusinteraction<-gls(logratioSBWEmerge~logratioSBWpeakplot+Plot,data=communitydietswitch, method="ML")
anova(SBWALTdietswitchmlfull,SBWALTdietswitchmlminusinteraction)

# Interaction is significant (p-value = 0.0058, delta AIC = 6.301303). Therefore, keep all fixed effects in the model including the interaction.


#7. Refit the model found in step 6. with REML estimation
SBWALTdietswitchreml<-gls(logratioSBWEmerge~logratioSBWpeakplot*Plot,data=communitydietswitch, method="REML")

summary(SBWALTdietswitchreml)

# average intercept for all plots
avinterceptallplotsDS<-(((3*coef(SBWALTdietswitchreml)[1])+coef(SBWALTdietswitchreml)[3]+coef(SBWALTdietswitchreml)[4])/3)

# average standard error of slopes for all plots
avSEinterceptallplotsDS <- ((summary(SBWALTdietswitchreml)$tTable[1,2]+summary(SBWALTdietswitchreml)$tTable[3,2]+summary(SBWALTdietswitchreml)$tTable[4,2])/3)

#average slope for all plots
avslopeallplotsDS<-(((3*coef(SBWALTdietswitchreml)[2])+coef(SBWALTdietswitchreml)[5]+coef(SBWALTdietswitchreml)[6])/3)
  
#average standard error of slopes for all plots
avSEslopeallplotsDS<-((summary(SBWALTdietswitchreml)$tTable[2,2]+summary(SBWALTdietswitchreml)$tTable[5,2]+summary(SBWALTdietswitchreml)$tTable[6,2])/3)

#Test of whether the average slope of all plots is different from 1 (test of frequency dependent selection)
2*pt(abs((avslopeallplotsDS-1)/avSEslopeallplotsDS),df=SBWALTdietswitchreml$dims[[1]]-SBWALTdietswitchreml$dims[[2]],lower.tail = FALSE) #p>0.05 so b is not significantly different from 1

#Calculate V (measure of frequency independent diet switching) using Greenwood and Elton method to calulate V.
10^(avinterceptallplotsDS/avslopeallplotsDS) #V=1.210321


#Calculating confidence limits of V using the quadratic function given in p731 in GreenWood and Elton (1979). 
result <- function(a,b,c){
  if(delta(a,b,c) > 0){ # first case D>0
    x_1 = (-b+sqrt(delta(a,b,c)))/(2*a)
    x_2 = (-b-sqrt(delta(a,b,c)))/(2*a)
    result = c(x_1,x_2)
  }
  else if(delta(a,b,c) == 0){ # second case D=0
    x = -b/(2*a)
  }
  else {"There are no real roots."} # third case D<0
}

# Constructing delta
delta<-function(a,b,c){
  b^2-4*a*c
} #this function is to find the two roots of the quadratic formula below

#To find a, b, and c of the quadratic formule in Greenwood and Elton (1979), found the confidence intervals of the model parameters and then used the difference between the CI and the parameter in the necessary places in the quadratic formula ###need to deal with confidence intervals correctly

aconfint <- avSEinterceptallplotsDS*1.96
bconfint <- avSEslopeallplotsDS*1.96

#The function result finds the two roots of the quadratic formula and 10 is then raised to the power of each root to convert the roots into the V confidence intervals.
logvconfint<-10^result(((avslopeallplotsDS^2)-(bconfint^2)),2*((avslopeallplotsDS*avinterceptallplotsDS)+((bconfint^2)*mean(communitydietswitch$logratioSBWpeakplot))),(avinterceptallplotsDS^2)-(aconfint^2))

#Create graph of log10 ratio parasitoid emergence from spruce budworm and other caterpillars to the log10 ratio of total spruce budworm to other caterpillars
SBWALTdietswitchremlGLS<-Gls(logratioSBWEmerge~logratioSBWpeakplot*Plot,data=communitydietswitch, method="REML") #using the function Gls because allows Predict function needed to plot trendline

dietswitchpredict<-Predict(SBWALTdietswitchremlGLS, logratioSBWpeakplot=seq(from=-1.5,to=2.1,by=0.1),Plot=seq(from=1,to=3,by=1))%>%
  mutate(Plot=as.factor(Plot))

comswitchlogratioplot<-ggplot(communitydietswitch,aes(logratioSBWpeakplot, logratioSBWEmerge))+
  geom_abline(intercept=0,slope=1,linetype="twodash")+
  geom_abline(intercept=0,slope=0)+
  geom_vline(xintercept=0)+
  geom_line(data=dietswitchpredict,aes(logratioSBWpeakplot,yhat,linetype=Plot),size=2)+
  geom_point(aes(shape=Plot,colour=bpa),size=5)+
  theme(axis.title.y=element_text(hjust=0.5, vjust=1.5), 
        legend.text=element_text(size=14),legend.justification=c(1,0), legend.position=c(0.98,0.02),legend.box = "horizontal")+coord_fixed(ratio = 1)+scale_color_viridis(name="Peak", breaks=c("b","p","a"),labels=c("Before","During","After"), alpha = 1, begin = 0, end = 1, direction = 1, discrete = TRUE, option = "D")+scale_shape_manual(values=c(16,17,15),guide=FALSE)+scale_linetype_manual(values=c("dotted","dashed","solid"))+ylab("Log10 emergence \nspruce budworm:other caterpillars")+xlab("Log10 abundance \nspruce budworm:other caterpillars")


ggsave("figs/communitydietswitchlogratio.pdf", plot=comswitchlogratioplot, width=8,height=8) #Figure 1a

# Absolute abundances of spruce budworm and other caterpillars over time
numSBWALT<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  mutate(SBWALT=as.factor(ifelse(SpecID.low=="h01","SBW","ALT")))%>%
  group_by(Peak,Plot,SBWALT)%>%
  summarise(NSBWALT=length(SpecID.low))%>%
  ungroup%>%
  right_join(numSBWALTnparaprep,by=c("Peak","Plot","SBWALT"))%>%
  mutate(dummyPeak=Peak+4,NSBWALT=as.numeric(NSBWALT),Plot=as.factor(Plot),SBWALT=as.factor(SBWALT), NSBWALT=ifelse(is.na(NSBWALT),0,NSBWALT),label=as.factor(ifelse(SBWALT=="SBW", "spruce budworm","other caterpillars"))) #creates a data frame with the number of spruce budworm and other caterpillars for every year and plot.

##Linear mixed effects model of number of interactions between shared parasitoids and alternative hosts and spruce budworm

#All of the statistical analysis below follows the protocol outlined in Zuur, A., E. N. Ieno, N. Walker, A. A. Saveliev, and G. M. Smith. 2009. Mixed effects models and extensions in ecology with r. First edition. Springer-Verlag New York, New York, New York, United States of America.

#Data Exploration
#1.Outliers in the response and explanatory variables.
op <- par(mfrow = c(1, 2), mar = c(3, 3, 3, 1))
dotchart(numSBWALT$NSBWALT, main = "Number of individuals", group=numSBWALT$SBWALT)
dotchart(numSBWALT$Peak, main = "Peak",group=numSBWALT$SBWALT)
par(op)
#2. Collinearity of the explanatory variables.
#Not doing as do not have multiple continuous explanatory variables

#3. Relationships between the response variable and the explanatory variables.
xyplot(NSBWALT~Peak, data=numSBWALT, groups=SBWALT)

#Creating Model
#1.Start with a linear regression model that contains as many explanatory variables and their interactions as possible. 
SBWALTnumfulllinear<-lm(NSBWALT~Peak*SBWALT*Plot,data=numSBWALT)
numSBWALT$resid<-NA
numSBWALT$resid<-SBWALTnumfulllinear$residuals
#Test normality and homogeniety of residuals
op <- par(mfrow = c(3, 3), mar = c(5, 4, 1, 2))
plot(SBWALTnumfulllinear, add.smooth = FALSE, which = 1,col=numSBWALT$SBWALT)
hist(numSBWALT$resid, xlab = "Residuals", main = "")
plot(numSBWALT$Peak,numSBWALT$resid,xlab="Peak",ylab="Residuals")
plot(numSBWALT$Plot,numSBWALT$resid,xlab="Plot",ylab="Residuals")
boxplot(numSBWALT$resid~numSBWALT$SBWALT,xlab="Caterpillar",ylab="Residuals")
plot(numSBWALT$Peak[numSBWALT$SBWALT=="SBW"],numSBWALT$resid[numSBWALT$SBWALT=="SBW"],xlab="Peak",ylab="Residuals",main="SBW")
plot(numSBWALT$Peak[numSBWALT$SBWALT=="ALT"],numSBWALT$resid[numSBWALT$SBWALT=="ALT"],xlab="Peak",ylab="Residuals",main="ALT")
par(op)

acf(residuals(SBWALTnumfulllinear), na.action = na.pass,
    main = "Auto-correlation plot for residuals")
#Obviously violation of independence of x values (time) and heterogeneity of variation

#2. Repeat step 1 using the gls function from the nlme package using REML estimation but without any variance structures.
SBWALTnumfullgls<-gls(NSBWALT~dummyPeak*SBWALT*Plot,data=numSBWALT, method="REML")
#3. Choose an appropriate variance structure depending on graphical analysis above.

SBWALTnumvarIdent<-gls(NSBWALT~dummyPeak*SBWALT*Plot,weights=varIdent(form=~1|SBWALT),data=numSBWALT, method="REML")

SBWALTnumvarIdentautocorr<-gls(NSBWALT~dummyPeak*SBWALT*Plot,weights=varIdent(form=~1|SBWALT),correlation=corAR1(form=~dummyPeak|Plot/SBWALT),data=numSBWALT, method="REML")

anova(SBWALTnumfullgls,SBWALTnumvarIdent,SBWALTnumvarIdentautocorr)

#check assumptions graphically for model with autocorrelation
numSBWALT$fmnormresid<-NA
numSBWALT$fmnormresid<-resid(SBWALTnumvarIdentautocorr, type="normalized")
plot(SBWALTnumvarIdentautocorr, add.smooth = FALSE,col=numSBWALT$SBWALT)
hist(numSBWALT$fmnormresid, xlab = "Normalised Residuals", main = "")
plot(numSBWALT$Peak,numSBWALT$fmnormresid,xlab="Peak",ylab="Normalised Residuals")
boxplot(numSBWALT$fmnormresid~numSBWALT$SBWALT,xlab="Caterpillar",ylab="Normalised Residuals")
plot(numSBWALT$Peak[numSBWALT$SBWALT=="SBW"],numSBWALT$fmnormresid[numSBWALT$SBWALT=="SBW"],xlab="Peak",ylab="Normalised Residuals",main="SBW")
plot(numSBWALT$Peak[numSBWALT$SBWALT=="ALT"],numSBWALT$fmnormresid[numSBWALT$SBWALT=="ALT"],xlab="Peak",ylab="Normalised Residuals",main="ALT")
coplot(fmnormresid ~ dummyPeak | SBWALT,
       ylab = "Normalised residuals", data = numSBWALT)
acf(numSBWALT$fmnormresid, na.action = na.pass,
    main = "Auto-correlation plot for normalised residuals")

#AIC lowest for model with varIdent and autocorr (and ACF shows marked improvement)
#decided to include autocorrelation due to apriori knowledge about temporal correlation and improvements with autocorrelation

#4. Using the gls model with the selected variance structure from 3., try to find the optimal random structure.
#Not doing because Plot is a fixed effect (using reasoning from http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html)

#5. Compare the new gls model with the earlier results using AIC.


#6.Find the optimal ﬁxed component using the likelihood ratio test where ML estimation is used.
SBWALTnumvarIdentautocorrmlfull<-gls(NSBWALT~dummyPeak*SBWALT*Plot,weights=varIdent(form=~1|SBWALT),correlation=corAR1(form=~dummyPeak|Plot/SBWALT),data=numSBWALT, method="ML")
SBWALTnumvarIdentautocorrdrop3interactionml<-update(SBWALTnumvarIdentautocorrmlfull,.~.-dummyPeak:SBWALT:Plot)

anova(SBWALTnumvarIdentautocorrmlfull,SBWALTnumvarIdentautocorrdrop3interactionml) #three way interaction not significant

SBWALTnumvarIdentautocorrdropdummyPeakSBWALTint<-update(SBWALTnumvarIdentautocorrdrop3interactionml,.~.-dummyPeak:SBWALT)
anova(SBWALTnumvarIdentautocorrdrop3interactionml,SBWALTnumvarIdentautocorrdropdummyPeakSBWALTint)

SBWALTnumvarIdentautocorrdropdummyPeakPlotint<-update(SBWALTnumvarIdentautocorrdrop3interactionml,.~.-dummyPeak:Plot)
anova(SBWALTnumvarIdentautocorrdrop3interactionml,SBWALTnumvarIdentautocorrdropdummyPeakPlotint)

SBWALTnumvarIdentautocorrdropSBWALTPlotint<-update(SBWALTnumvarIdentautocorrdrop3interactionml,.~.-SBWALT:Plot)
anova(SBWALTnumvarIdentautocorrdrop3interactionml,SBWALTnumvarIdentautocorrdropSBWALTPlotint)

#Drop dummy peak and plot interaction
SBWALTnumvarIdentautocorrdropdummyPeakSBWALTint2<-update(SBWALTnumvarIdentautocorrdropdummyPeakPlotint,.~.-dummyPeak:SBWALT)
anova(SBWALTnumvarIdentautocorrdropdummyPeakPlotint,SBWALTnumvarIdentautocorrdropdummyPeakSBWALTint2)

SBWALTnumvarIdentautocorrdropSBWALTPlotint2<-update(SBWALTnumvarIdentautocorrdropdummyPeakPlotint,.~.-SBWALT:Plot)
anova(SBWALTnumvarIdentautocorrdropdummyPeakPlotint,SBWALTnumvarIdentautocorrdropSBWALTPlotint2)
#Drop SBWALT and plot interaction

SBWALTnumvarIdentautocorrdropdummyPeakSBWALTint3<-update(SBWALTnumvarIdentautocorrdropSBWALTPlotint2,.~.-dummyPeak:SBWALT)
anova(SBWALTnumvarIdentautocorrdropSBWALTPlotint2,SBWALTnumvarIdentautocorrdropdummyPeakSBWALTint3)

#keep dummyPeak and SBWALT and their interaction

SBWALTnumvarIdentautocorrdropPlot<-update(SBWALTnumvarIdentautocorrdropSBWALTPlotint2,.~.-Plot)
anova(SBWALTnumvarIdentautocorrdropSBWALTPlotint2,SBWALTnumvarIdentautocorrdropPlot)

#Drop Plot because not significant.

#7. Refit the model found in step 6. with REML estimation
SBWALTnumfinalmodel<-gls(NSBWALT~dummyPeak*SBWALT,weights=varIdent(form=~1|SBWALT),correlation=corAR1(form=~dummyPeak|Plot/SBWALT),data=numSBWALT, method="REML")
summary(SBWALTnumfinalmodel)

#Calculation of mean and coefficient of variation for spruce budworm and other caterpillars over all years.
numSBWALT%>%
  group_by(SBWALT)%>%
  summarise(meannum=mean(NSBWALT),stdevnum=sd(NSBWALT))%>%
  mutate(covnum=stdevnum/meannum)

#Create graph of number of spruce budworm and other caterpillars over time,
SBWALTnumvarIdentautocorrremlGLS<-Gls(NSBWALT~dummyPeak*SBWALT,weights=varIdent(form=~1|SBWALT),correlation=corAR1(form=~dummyPeak|Plot/SBWALT),data=numSBWALT, method="REML") #using the function Gls because allows Predict function needed to plot trendline
SBWALTnumpredict<-Predict(SBWALTnumvarIdentautocorrremlGLS, dummyPeak=c(1:14), SBWALT=c("ALT","SBW"))%>%
  mutate(label=as.factor(ifelse(SBWALT=="SBW", "spruce budworm","other caterpillars")),Peak=dummyPeak-4)

SBWALTnumplot<-SBWALTnumpredict%>%
  ggplot+
  geom_line(aes(Peak,yhat),size=2)+
  geom_point(data=numSBWALT,aes(Peak,NSBWALT,shape=Plot),size=5)+
  facet_wrap(~label,scales="free_y",ncol=1)+
  theme(strip.background=element_blank(),strip.text.x=element_text(size=28), 
        axis.title.y=element_text(hjust=0.5, vjust=1.5),legend.text=element_text(size=14),legend.justification=c(1,0), legend.position=c(0.98,0.2))+
  ylab("Caterpillar abundance")+xlab("Years before/after peak")+ylim(c(0,NA))+scale_x_continuous(breaks=c(-2,0,2,4,6,8,10))+scale_shape_manual(values=c(16,17,15))

ggsave("figs/SBWALTnumplot.pdf",plot=SBWALTnumplot,width=8,height=7) #figure 1b

# Number of parasitoid emergences from spruce budworm and other caterpillars

#Together, the function create.SBWALTparashareprep and the short code that creates the data frame called SBWALTparashareprep create a data frame where there is a row for every combination of peak, caterpillar type and parasitoid taxa but where there aren't rows for when the plots were not sampled (different plots were sampled in different years).
create.SBWALTparashareprep<-function(plot){
  x<-data.frame(expand.grid(c(-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10),c("SBW","ALT"),parashare))%>%
  rename(Peak=Var1,SBWALT=Var2,SpecID.high=Var3)%>%
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

SBWALTEM<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  group_by(Peak,Plot,SpecID.high)%>%
  summarise(SBW=length(SpecID.low[SBWALT=="SBW"]),ALT=length(SpecID.low[SBWALT=="ALT"]))%>%
  gather(SBWALT,NEM,4:5)%>%
  ungroup%>%
  right_join(SBWALTparashareprep,by=c("Peak","Plot","SBWALT","SpecID.high"))%>%
  mutate(Plot=as.factor(Plot),SBWALT=as.factor(SBWALT),NEM=ifelse(is.na(NEM),0,NEM))%>%
  group_by(Peak,Plot,SBWALT)%>%
  summarise(AvNEM=mean(NEM))%>%
  mutate(dummyPeak=Peak+4,label=as.factor(ifelse(SBWALT=="SBW", "spruce budworm","other caterpillars")))  #Creates a dataframe of the average number of emergences of parasitoids from spruce budworm or other caterpillars (all parasitoid taxa combined) for each peak and plot. Before taking the average, 0 values are added to any peak, plot, species combination that not have a parasitoid emerge from either spruce budworm or other caterpillars.
  
#Linear mixed effects model of log10 average number of emergences of parasitoids from spruce budworm or other caterpillars over time

#All of the statistical analysis below follows the protocol outlined in Zuur, A., E. N. Ieno, N. Walker, A. A. Saveliev, and G. M. Smith. 2009. Mixed effects models and extensions in ecology with r. First edition. Springer-Verlag New York, New York, New York, United States of America.

#Data Exploration
#1.Outliers in the response and explanatory variables.
op <- par(mfrow = c(2, 2), mar = c(3, 3, 3, 1))
dotchart(SBWALTEM$AvNEM, main = "Number of Interactions", group=SBWALTEM$SBWALT)
dotchart(SBWALTEM$Peak, main = "Peak")
par(op)
#2. Collinearity of the explanatory variables.
#Not doing as do not have multiple continuous explanatory varaibles

#3. Relationships between the response variable and the explanatory variables.
xyplot(AvNEM~Peak, data=SBWALTEM, groups=SBWALT)

#Creating the model
#1.Start with a linear regression model that contains as many explanatory variables and their interactions as possible. 
SBWALTEMfulllinear<-lm(AvNEM~Peak*SBWALT*Plot,data=SBWALTEM)
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

#Obviously violation of independence of x values (time), also violation of homogeneity assumption!#

#2. Repeat step 1 using the gls function from the nlme package using REML estimation but without any variance structures.
SBWALTEMfullgls<-gls(AvNEM~dummyPeak*SBWALT*Plot,data=SBWALTEM, method="REML")

#3. Choose an appropriate variance structure depending on graphical analysis above.
SBWALTEMfixed<-gls(AvNEM~dummyPeak*SBWALT*Plot,weights=varFixed(~dummyPeak),data=SBWALTEM, method="REML")
SBWALTEMcomb<-gls(AvNEM~dummyPeak*SBWALT*Plot,weights=varComb(varIdent(form=~1|SBWALT),varFixed(~dummyPeak)),data=SBWALTEM, method="REML")
SBWALTEMident<-gls(AvNEM~dummyPeak*SBWALT*Plot,weights=varIdent(form=~1|SBWALT),data=SBWALTEM, method="REML")
SBWALTEMpower<-gls(AvNEM~dummyPeak*SBWALT*Plot,weights=varPower(form=~dummyPeak|SBWALT),data=SBWALTEM, method="REML")

anova(SBWALTEMfullgls,SBWALTEMfixed,SBWALTEMident,SBWALTEMpower,SBWALTEMcomb)

#From AIC values varident is best method to deal with heterogeneity (AIC is 267.2292 versus 357.8575 for without corr structure and 269.3425 for varPower)
#Graphical validation of optimal model
E1 <- resid(SBWALTEMident)
coplot(E1 ~ dummyPeak | SBWALT,
         ylab = "Ordinary residuals", data = SBWALTEM)
E2 <- resid(SBWALTEMident, type="normalized")
coplot(E2 ~ dummyPeak | SBWALT,
       ylab = "Normalised residuals", data = SBWALTEM)

E1b <- resid(SBWALTEMpower)
coplot(E1b ~ dummyPeak | SBWALT,
       ylab = "Ordinary residuals", data = SBWALTEM)
E2b <- resid(SBWALTEMpower, type="normalized")
coplot(E2b ~ dummyPeak | SBWALT,
       ylab = "Normalised residuals", data = SBWALTEM)
##From looking at Normalised residuals of varPower and varIdent there is little difference between varPower and varIdent (but varIdent has lower AIC)
###Chosen varIdent as the variance structure

#autocorrelation
SBWALTEMautocorrvarident<-gls(AvNEM~dummyPeak*SBWALT*Plot,weights=varIdent(form=~1|SBWALT),correlation=corAR1(form=~dummyPeak|Plot/SBWALT),data=SBWALTEM, method="REML")
anova(SBWALTEMident,SBWALTEMautocorrvarident)
#ACF
acf(residuals(SBWALTEMident, type = "normalized"), na.action = na.pass,
    main = "Auto-correlation plot for residuals")
acf(residuals(SBWALTEMautocorrvarident, type = "normalized"), na.action = na.pass,
    main = "Auto-correlation plot for residuals")
E2b <- resid(SBWALTEMautocorrvarident, type="normalized")
coplot(E2b ~ dummyPeak | SBWALT,
       ylab = "Normalised residuals", data = SBWALTEM)

#decided to remove corAR1 because does not improve model

#4. Using the gls model with the selected variance structure from 3., try to find the optimal random structure.
#Not doing because Plot is a fixed effect (using reasoning from http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html)

#5. Compare the new gls model with the earlier results using AIC.

#6.Find the optimal ﬁxed component using the likelihood ratio test where ML estimation is used.
SBWALTEMfullml<-gls(AvNEM~dummyPeak*SBWALT*Plot,weights=varIdent(form=~1|SBWALT),data=SBWALTEM, method="ML")

SBWALTEMdrop3interactionml<-update(SBWALTEMfullml,.~.-dummyPeak:SBWALT:Plot)
anova(SBWALTEMfullml,SBWALTEMdrop3interactionml)

#drop interaction of dummyPeak,SBWALT, and Plot because not significant

SBWALTEMdropdummyPeakSBWALTml<-update(SBWALTEMdrop3interactionml,.~.-dummyPeak:SBWALT)
anova(SBWALTEMdrop3interactionml,SBWALTEMdropdummyPeakSBWALTml)

SBWALTEMdropdummyPeakPlotml<-update(SBWALTEMdrop3interactionml,.~.-dummyPeak:Plot)
anova(SBWALTEMdrop3interactionml,SBWALTEMdropdummyPeakPlotml)

SBWALTEMdropSBWALTPlotml<-update(SBWALTEMdrop3interactionml,.~.-SBWALT:Plot)
anova(SBWALTEMdrop3interactionml,SBWALTEMdropSBWALTPlotml)

#drop interaction of dummyPeak Plot because least significant

SBWALTEMdropdummyPeakSBWALTml2<-update(SBWALTEMdropdummyPeakPlotml,.~.-dummyPeak:SBWALT)
anova(SBWALTEMdropdummyPeakPlotml,SBWALTEMdropdummyPeakSBWALTml2)
 
SBWALTEMdropSBWALTPlotml2<-update(SBWALTEMdropdummyPeakPlotml,.~.-SBWALT:Plot)
anova(SBWALTEMdropdummyPeakPlotml,SBWALTEMdropSBWALTPlotml2)

#drop interaction of SBWALT and Plot because p value is 0.04 and Zuur drops effects at this level of significance.

SBWALTEMdropdummyPeakSBWALTml3<-update(SBWALTEMdropSBWALTPlotml2,.~.-dummyPeak:SBWALT)
anova(SBWALTEMdropSBWALTPlotml,SBWALTEMdropdummyPeakSBWALTml3)

#do not drop dummyPeak and SBWALT and their interaction

SBWALTEMdropPlotml<-update(SBWALTEMdropSBWALTPlotml2,.~.-Plot)
anova(SBWALTEMdropSBWALTPlotml2,SBWALTEMdropPlotml)

# drop Plot because p-value = 0.3747

#7. Refit the model found in step 6. with REML estimation
SBWALTEMfinalmodel<-gls(AvNEM~dummyPeak*SBWALT,weights=varIdent(form=~1|SBWALT),data=SBWALTEM, method="REML") 

summary(SBWALTEMfinalmodel)



#Create graph of number of parasitoid emergences from spruce budworm and other caterpillars over time,
SBWALTEMvaridentremlGLS<-Gls(AvNEM~dummyPeak*SBWALT,weights=varIdent(form=~1|SBWALT),data=SBWALTEM, method="REML") #using the function Gls because allows Predict function needed to plot trendline
AvNEMpredict<-Predict(SBWALTEMvaridentremlGLS, dummyPeak=c(1:14), SBWALT=c("ALT","SBW"))%>%
  mutate(label=as.factor(ifelse(SBWALT=="SBW", "spruce budworm","other caterpillars")),Peak=dummyPeak-4)

SBWALTEMplot<-AvNEMpredict%>%
  ggplot+
  geom_line(aes(Peak,yhat),size=2)+
  geom_point(data=SBWALTEM,aes(Peak,AvNEM,shape=Plot),size=5)+
  facet_wrap(~label,scales="free_y",ncol=1)+
  theme(strip.background=element_blank(),strip.text.x=element_text(size=28), 
        axis.title.y=element_text(hjust=0.5, vjust=1.5))+
  ylab("Number of parasitoid \nemergences")+xlab("Years before/after peak")+ylim(c(0,NA))+scale_shape_manual(values=c(16,17,15),guide=FALSE)

ggsave("figs/SBWALTEMplot.pdf",plot=SBWALTEMplot,width=8,height=7) #figure 1c

# Single species response hypothesis versus aggregate response-----------------------------------------------------------
numinteractionsspecies<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  group_by(SpecID.high)%>%
  summarise(numinter=length(SpecID.high))%>%
  arrange(desc(numinter)) # Produces a dataframe of the total number of emergences with spruce budworm and other caterpillars for each parasitoid taxa. Then arranged descending from most emergences to lowest emergences.

paralist<-numinteractionsspecies$SpecID.high

drop_top3_sp <- function(x){
communitydietswitchlog_dropsp<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  filter(!SpecID.high %in% x)%>%
  group_by(Peak, Plot)%>%
  summarise(TotSBWEmerge=length(SpecID.high[SBWALT=="SBW"]),TotALTEmerge=length(SpecID.high[SBWALT=="ALT"]))%>%
  mutate(ratioSBWEmerge=TotSBWEmerge/TotALTEmerge,logratioSBWEmerge=log10(ratioSBWEmerge))%>%
  full_join(ratioSBWbypeakplot,by=c("Peak","Plot"))%>%
  filter(TotSBW!=0,TotALT!=0,TotSBWEmerge!=0,TotALTEmerge!=0)%>%
  mutate(logratioSBWpeakplot=log10(ratioSBWpeakplot),dummyPeak=Peak+4)

SBWALTdietswitchreml_dropsp<-gls(logratioSBWEmerge~logratioSBWpeakplot*Plot,data=communitydietswitchlog_dropsp,method="REML")

avslopeallplotsDS_dropsp<-(((3*coef(SBWALTdietswitchreml_dropsp)[2])+coef(SBWALTdietswitchreml_dropsp)[5]+coef(SBWALTdietswitchreml_dropsp)[6])/3)

avSEslopeallplotsDS_dropsp<-((summary(SBWALTdietswitchreml_dropsp)$tTable[2,2]+summary(SBWALTdietswitchreml_dropsp)$tTable[5,2]+summary(SBWALTdietswitchreml_dropsp)$tTable[6,2])/3)


t<- (avslopeallplotsDS_dropsp-avslopeallplotsDS)/avSEslopeallplotsDS_dropsp

ttable<-tibble(
  speciesdrop=paste(x,collapse=" "),
  origbeta=avslopeallplotsDS,
  beta=avslopeallplotsDS_dropsp,
  diffbeta=1-beta/origbeta,
  int=avSEslopeallplotsDS_dropsp*1.96,
  tstat=t,
  df=SBWALTdietswitchreml_dropsp$dims[[1]]-SBWALTdietswitchreml_dropsp$dims[[2]],
  p=pt(abs(t),df=SBWALTdietswitchreml_dropsp$dims[[1]]-SBWALTdietswitchreml_dropsp$dims[[2]],lower.tail = FALSE),
  p2=2*p
)
return(ttable)
}

bind_rows(drop_top3_sp("p01"), drop_top3_sp(c("p01","p02")),drop_top3_sp(c("p01","p02","p05")))#drop the most abundant parasitoid taxon (p01), then drop the two most abundant species, then drop the three most abundant species. compare coefficient values (average of all three plots), run one sample ttests and combine into one table.


#nMDS analysis of parasitoid community over time
turnoverprep<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  group_by(Peak,Plot)%>%
  summarise(Totemerge=length(SpecID.high))#creates a dataframe of the total emergences of all parasitoid taxa combined for each peak and plot (used to standardise the number of emergences per species by the total number of emergences because different years had different total number of emergences)

turnover<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  group_by(Peak,Plot, SpecID.high)%>%
  summarise(nemerge=length(SpecID.high))%>%
  ungroup()%>%
  left_join(turnoverprep,by=c("Peak","Plot"))%>%
  mutate(nemergestand=nemerge/Totemerge, bpa=case_when(
    Peak %in% c(-3,-2)~"b",
    Peak %in% c(-1,0,1)~"p",
    Peak %in% c(2,3,4,5,6,7,8,9,10)~"a"
  ))%>%
  mutate(bpa=as.factor(bpa))#creates a data frame of the number of emergences for each parasitoid taxa standardised by the total number of emergences for all species combined for each peak and plot. adds a grouping variable (bpa) which splits the peaks into three groups: before the peak, during the peak, and after the peak.

turnovercommat<-turnover%>%
  select(-c(nemerge,Totemerge))%>%
  spread(SpecID.high,nemergestand)%>%
  ungroup()#creates a community matrix (but still including the factors) of the dataframe turnover.

factors<-tibble(
  Plot=turnovercommat$Plot,
  Peak=turnovercommat$Peak,
  bpa=turnovercommat$bpa
)#creates a dataframe of the different factors that will be used in the nMDS and the PERMANOVA

turnovercommat%<>%select(-c(Peak,Plot,bpa))#removes the factors from the dataframe to create a true community matrix

turnovercommat[is.na(turnovercommat)]<-0#replaces any NA values with 0. these NA values were created when the long format data was spread into a wide format dataframe and where species were not found for a certain Peak or Plot.

dimcheckMDS(turnovercommat, distance = "bray", k = 6, trymax = 30, 
        autotransform=TRUE) #after exmining screeplot from 6 dimensions to 1 dimension, chosen to use 2 dimensions because additional dimensions provide small reductions in stress

turnoverNMDS <- metaMDS(turnovercommat, distance = "bray", k = 2, trymax = 50, 
                        autotransform=TRUE)

stressplot(turnoverNMDS)  ##to test if distances as represented in ordination are correlated with actual distances

stress10<-c(0.09293078,0.0953679,0.101168,0.1065564,0.0929335,0.1400294,0.09614226,0.1065557,0.09536982,0.1066436)
stress5<-c(0.09293078,0.0953679,0.101168,0.1065564,0.0929335)
sd(stress10)
sd(stress5)

#Create nMDS plot
turnoverdatascores<-as.data.frame(scores(turnoverNMDS))
turnoverdatascores$Plot<-factors$Plot
turnoverdatascores$bpa<-factors$bpa

# function for creating ellipses in nmds plot
#adapted from  http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
veganCovEllipsenew <- function (x, scale = 1, npoints = 100) 
{
  cov <- cov.wt(cbind(x$NMDS1,x$NMDS2),wt=rep(1/length(x$NMDS1),length(x$NMDS1)))$cov 
  center <- c(mean(x$NMDS1),mean(x$NMDS2))
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  bpaunique<-unique(x$bpa)
  tibble(
    NMDS1=t(center + scale * t(Circle %*% chol(cov)))[,1],
    NMDS2=t(center + scale * t(Circle %*% chol(cov)))[,2],
    bpa=bpaunique
  )
}


df_ell.bpa<-turnoverdatascores%>%
  split(.$bpa)%>%
  map(veganCovEllipsenew)%>%
  bind_rows()

turnoverNMDSplot<-ggplot(turnoverdatascores)+
  geom_point(aes(NMDS1,NMDS2, shape=Plot, colour=bpa),size=5)+
  geom_path(data=df_ell.bpa, aes(x=NMDS1, y=NMDS2, colour=bpa), size=1.5)+
  theme(axis.title.y=element_text(hjust=0.5, vjust=1.5), 
        legend.text=element_text(size=14),
        legend.title=element_text(size=15),
        legend.justification=c(1,0.99), 
        legend.position=c(1,0.99),
        legend.box = "horizontal")+
  coord_fixed(ratio = 1)+scale_color_viridis(name="Peak", breaks=c("b","p","a"),labels=c("Before","During","After"), alpha = 1, begin = 0, end = 1, direction = 1, discrete = TRUE, option = "D")+scale_shape_manual(values=c(16,17,15))+scale_x_continuous(breaks=c(-1,0,1))

ggsave("figs/nmdsturnover.pdf",plot=turnoverNMDSplot,width=10,height=10) #Figure 2


# Testing the separation of the parasitoid communities between before, during, and after the spruce budworm peak, using the adonis function for a permanova. - Guide to formulating the permutation structure and permutation test from the following webpage http://thebiobucket.blogspot.ca/2011/04/repeat-measure-adonis-lately-i-had-to.html#more

fit <- adonis(vegdist(turnovercommat)~bpa, factors, permutations=1)

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
  fit.rand <- adonis(vegdist(turnovercommat) ~ bpa[idx], factors,permutations = 1)
  pop[i] <- fit.rand$aov.tab[1, 4]
}

# Return the p-value of the full permanova:
print(pval <- sum(pop >= pop[1]) / (B + 1))

# Identify which parasitoid taxa indicate a specific point in time, using a indicator species analysis
turnoverisa<-multipatt(turnovercommat,cluster=factors$bpa,control=ctrl)
summary(turnoverisa)
summary(turnoverisa, indvalcomp = TRUE)#Component A is the probability that the surveyed site belongs to the target site group given that the species has been found and B is probability of the species as indicator of the site group
summary(turnoverisa,alpha=1)#to see all selected species and their associations with each site


# Fluctuating interaction strength distributions ---------

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
  frame2webs(varnames=c("SBWALT","SpecID.highnum","Peak")) # creates a dataframe that counts the number of interactions between each parasitoid taxa and spruce budworm or other caterpillars for each year.

## Create bipartite graphs to identify the changing interaction strength distributions between parasitoids, spruce budworm, and other caterpillars

# Calculation of percentage of emergences from spruce budworm and other caterpillas for each bipartite graph (in the main text)
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
plotweb(as.matrix(bipartitegraphSBWALT$'minus3'),ybig=1.6,low.spacing=0.8,high.spacing=0.1,x.lim=c(0,3.2),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.6,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

#0
pdf("figs/0plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'0'),ybig=1.7,low.xoff = 1.5,low.spacing=2,high.spacing=0.2,x.lim=c(0,6.9),labsize=2.5,high.lab.dis=0.05,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

#3
pdf("figs/3plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'3'),ybig=1.6,low.xoff = 1.5,low.spacing=1.8,high.spacing=0.2,x.lim=c(0,6.7),labsize=2.5,high.lab.dis=0.05,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

#10
pdf("figs/10plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'10'),ybig=1.6,high.spacing=0.1,x.lim=c(0,1.45),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.16,col.low="#C73D73FF", col.high="#160F3BFF",col.interaction="#FECC8FFF", method='normal')
dev.off()


#calculation of median interactions strength to maximum interaction strength

intdist<-bipartitegraphSBWALT%>%
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
  mutate(Year=names(bipartitegraphSBWALT))# creates a dataframe of median to maximum interaction strength (number of emergences is a proxy for interaction strength) for parasitoids emerging from spruce budworm (SBW) and other caterpillars (ALT).

intdist$Year<-gsub("minus","-",intdist$Year)

intdistlong <- intdist %>%
  select(medmaxratioALT,medmaxratioSBW,Year) %>%
  rename(ALT=medmaxratioALT,SBW=medmaxratioSBW) %>%
  gather(SBWALT, medmax, -Year)

intdistplot<-ggplot(intdistlong)+
  geom_point(aes(as.numeric(Year),medmax, colour=SBWALT), size=5)+
  geom_smooth(aes(as.numeric(Year),medmax, colour=SBWALT),method=lm, formula=y~poly(x,2),se=FALSE, size=3)+
  xlab("Years before/after peak")+ylab("Median:maximum \ninteraction strength")+
  theme(legend.text=element_text(size=14),
        legend.title=element_text(size=15),
        legend.justification=c(0,1),
        legend.position=c(0.01,0.99))+
  scale_color_manual(name="caterpillar type", breaks=c("SBW","ALT"),labels=c("spruce budworm","other caterpillars"), values=c("#C73D73FF","#2D1160FF"))+
    scale_x_continuous(breaks=c(-2,0,2,4,6,8,10))

ggsave("figs/intdistplot.pdf",intdistplot,width=8,height=7) #Figure 4

#Parasitoid exodus -----

nparaSBWALT<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  mutate(SBWALT=as.factor(ifelse(SpecID.low=="h01","SBW","ALT")))%>%
  group_by(Peak,Plot, SBWALT)%>%
  summarise(npara=length(unique(SpecID.high)))%>%
  ungroup%>%
  right_join(numSBWALTnparaprep,by=c("Peak","Plot","SBWALT"))%>%
  mutate(dummyPeak=Peak+4,npara=as.numeric(npara),Plot=as.factor(Plot),SBWALT=as.factor(SBWALT), npara=ifelse(is.na(npara),0,npara), label=as.factor(ifelse(SBWALT=="SBW", "spruce budworm","other caterpillars")))#Creates a dataframe of the overall number of parasitoid taxa that emerge from spruce budworm and other caterpillars in every year and for each plot.

#Linear mixed effects model of the number of parasitoid taxa from spruce budworm and other caterpillars over time.

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


anova(SBWALTnparafullgls,SBWALTnparaautocorr)#Delta AIC=-4.0827

acf(residuals(SBWALTnparaautocorr, type = "normalized"), na.action = na.pass,
    main = "Auto-correlation plot for residuals")

#Small reduction in AIC with autocorrelation (still greater than 2). small improvement when looking at ACF graph.
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

# Calculation of difference in mean parasitoid diversity between spruce budworm and other caterpillars
nparaSBWALT%>%
  group_by(SBWALT)%>%
  summarise(MeanSBWALT=mean(npara))%>%
  spread(SBWALT,MeanSBWALT)

# Create graph of number of parasitoid taxa emerging from spruce budworm and other caterpillars over time,
SBWALTnparaautocorrremlGLS<-Gls(npara~dummyPeak+SBWALT+Plot+dummyPeak:SBWALT+dummyPeak:Plot,correlation=corAR1(form=~dummyPeak|Plot/SBWALT),data=nparaSBWALT, method="REML") #using the function Gls because allows Predict function needed to plot trendline
nparapredict<-Predict(SBWALTnparaautocorrremlGLS, dummyPeak=c(1:14), Plot=seq(from=1,to=3,by=1),SBWALT=c("ALT","SBW"))%>%
  mutate(label=as.factor(ifelse(SBWALT=="SBW", "spruce budworm","other caterpillars")),Peak=dummyPeak-4,Plot=as.factor(Plot))
  
nparaSBWALTpredict<-left_join(nparaSBWALT,nparapredict,by=c("Plot", "Peak", "SBWALT", "label"))

nparaSBWALTplot<-ggplot(nparaSBWALTpredict)+
  geom_line(aes(Peak,yhat,linetype=Plot), size=2)+
  geom_point(aes(Peak,npara,shape=Plot), size=5)+
  facet_wrap(~label,scales="fixed",ncol=1)+
  theme(strip.background=element_blank(),strip.text.x=element_text(size=28),
        axis.title.y=element_text(hjust=0.5, vjust=1.5))+
  ylab("Richness of \nparasitoid taxa")+xlab("Years before/after peak")+scale_y_continuous(limits=c(0, 26))+scale_x_continuous(breaks=c(-2,0,2,4,6,8,10))+scale_shape_manual(values=c(16,17,15))+scale_linetype_manual(values=c("dotted","dashed","solid"))

ggsave("figs/nparaSBWALTplot.pdf",plot=nparaSBWALTplot,width=8,height=7) #Figure 5


# Supplementary Information -----------------------------------------------

#Creation of bipartite graphs for extra years (Figures 6 and 7)
# Calculation of percentage of emergences from spruce budworm and other caterpillas for each bipartite graph (in the supplementary information)
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
plotweb(as.matrix(bipartitegraphSBWALT$'minus2'),ybig=1.5,low.spacing=0.8,high.spacing=0.1,x.lim=c(0,3.4),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.7,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/minus1plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'minus1'),ybig=1.5,low.spacing=1,high.spacing=0.1,x.lim=c(0,3.8),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.8,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/1plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'1'),ybig=1.5,low.spacing=1,high.spacing=0.1,x.lim=c(0,3.8),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.8,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/2plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'2'),ybig=1.5,low.spacing=0.8,high.spacing=0.1,x.lim=c(0,3.2),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.5,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/4plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'4'),ybig=1.5,low.spacing=0.8,high.spacing=0.1,x.lim=c(0,2.9),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.5,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/5plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'5'),ybig=1.5,low.spacing=0.6,high.spacing=0.1,x.lim=c(0,2.3),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.4,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/6plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'6'),ybig=1.5,low.spacing=0.6,high.spacing=0.1,x.lim=c(0,2.1),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.2,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/7plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'7'),ybig=1.5,low.spacing=0.6,high.spacing=0.1,x.lim=c(0,2.4),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.2,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/8plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'8'),ybig=1.5,high.spacing=0.1,x.lim=c(0,1.5),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.3,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()


pdf("figs/9plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'9'),ybig=1.5,high.spacing=0.1,x.lim=c(0,1.3),labsize=2.5,high.lab.dis=0.05,low.xoff = 0.2,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

# Useful Extra Analyses ---------------------------------------------------
#Proportion of parasitism of spruce budworm by parasitoid taxa that emerge from both spruce budworm and other caterpillars
allsharedprop<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(!is.na(SpecID.high))%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  group_by(SBWALT)%>%
  summarise(NSha=length(SpecID.low[SpecID.high %in% parashare]),NUnSha=length(SpecID.low[!SpecID.high %in% parashare]),total=length(SpecID.low))%>%
  mutate(PropSh=NSha/total)#90% of all parasitism emergences from spruce budworm were done by species that share spruce budworm and other caterpillars - 84% of all parasitism emergences from other caterpillars were done by shared parasitoids

numinteractionsspecies<-interaction%>%
  filter(grepl("h",SpecID.low))%>%
  filter(SpecID.high %in% parashare)%>%
  group_by(SpecID.high)%>%
  summarise(numinter=length(SpecID.high))%>%
  arrange(desc(numinter))

library(scales)
show_col(viridis_pal(option="A")(25))

# Notes -------------------------------------------------------------------

#Wintemia fumiferana (Smidtia fumiferana), a major parasitoid of spruce budworm was not included in the shared parasitoid list of parasitoids because found to interact with undetermined herbivores once and with no other other caterpillars (but found to interact with spruce budworm  1241 times) . undetermined herbivores were removed from dataset to follow Eveleigh et al. 2007 SI Methods.


