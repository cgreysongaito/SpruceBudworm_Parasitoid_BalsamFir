#############################################################################
#            Parasitoid community responds indiscriminately to fluctuating spruce budworm and other caterpillars on balsam fir
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
library(permute)

# Data Preparation -----
##Input data as data.table
interactionimport<-read_csv("data/SpruceBudworm_OtherCaterpillars_BalsamFir_FoodWeb_EldonEveleigh.csv")

#find parasitoids that could not be identified
interactionun <- interactionimport %>%
  filter(grepl("undetermined", ParasitoidLatinName) | grepl("icho", ParasitoidLatinName)) %>%
  filter(!ParasitoidLatinName %in% c("Dolichogenidea absona","Dolichogenidea renaulti"))

paraund <- unique(interactionun$SpecID.high)

#Remove unidentified parastiods
interaction <-interactionimport %>%
  filter(!SpecID.high %in% paraund)%>% #remove all unresolved taxa
  filter(!SpecID.high %in% c("p51")) %>% #remove diadegma sp. that has freq.low of 2 and uncertain as to the species (either D. pulicalivariae-p94 or Diadegma sp1-p95)
  mutate(Plot = as.factor(Plot))

# Produce a list of the parasitoids that attack budworm
paraSBWprep<-interaction%>%
  filter(!grepl("p",SpecID.low))%>%
  filter(!is.na(SpecID.high))%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  filter(SBWALT!="ALT")

paraSBW <- unique(paraSBWprep$SpecID.high)
#---

speciesnamesSBW <- interaction %>%
  select(SpecID.high, ParasitoidLatinName) %>%
  filter(SpecID.high %in% paraSBW)%>%
  unique() # Lists the latin names of all parasitoid species that attack budworm

#The function, create.numSBWALTnparaprep, and the second section of code creates a dataframe with all combinations of Peak (-3 to 10), SBWALT (either SBW or ALT), and Plot (1,2,3) to be ensure all data sets used in analyses have all combinations

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
  mutate(Plot=as.factor(Plot)) 
#---

# Produces a list where the parasitoid taxa are sorted by the number of emergences from budworm. Each parasitoid taxa is given a number (1 for the taxa with the most emergences from budworm and 48 for the least emergences from budworm).
numSBWALTbyspeciesdata<-interaction%>%
  filter(SpecID.high %in% paraSBW)%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  group_by(SpecID.high)%>%
  summarise(NSBW=length(SpecID.high[SBWALT=="SBW"]),NALT=length(SpecID.high[SBWALT=="ALT"]))%>%
  left_join(speciesnamesSBW,by="SpecID.high") %>%
  select(-SpecID.high)%>%
  rename(SpecID.high=ParasitoidLatinName)%>%
  mutate(SpecID.high=as.factor(SpecID.high))

proplevels<-levels(reorder(numSBWALTbyspeciesdata$SpecID.high,numSBWALTbyspeciesdata$NSBW))

speciesnamesprep<-tibble(
  SpecID.high=as.factor(rev(proplevels)),
  numID=seq(1,48,by=1)
) #produces a dataset where each parasitoid taxon has a number. used in producing the bipartite graphs.
#---

#Parasitoid community host preference-----

#community level plotting of diet switching using log10 ratios of emergences from budworm and other caterpillars and log10 ratios of total abundances of budworm and other caterpillars

ratioSBWbypeakplot<-interaction%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  group_by(Peak, Plot)%>%
  summarise(TotSBW=length(SBWALT[SBWALT=="SBW"]),TotALT=length(SBWALT[SBWALT=="ALT"]))%>%
  mutate(ratioSBWpeakplot=TotSBW/TotALT, logratioSBWpeakplot=log10(ratioSBWpeakplot)) #produces a dataset of the ratio (and log 10 ratio) of the total number of budworm and other caterpillars sampled on balsam fir

ratioSBWparacommunity<-interaction %>%
  filter(SpecID.high %in% paraSBW)%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  group_by(Peak, Plot)%>%
  summarise(TotSBWEmerge=length(SpecID.high[SBWALT=="SBW"]),TotALTEmerge=length(SpecID.high[SBWALT=="ALT"])) #produces a data set of the ratio (and log 10 ratio) of the total number of parasitoid emergences (all taxa combined) from budworm and other caterpillars sampled on balsam fir.

communitydietswitch<-full_join(ratioSBWparacommunity,ratioSBWbypeakplot,by=c("Peak","Plot")) %>%
  filter(TotSBW!=0,TotALT!=0,TotSBWEmerge!=0,TotALTEmerge!=0)%>% #Removed four data points (from a total of 25 data points).
  ungroup()%>%
  mutate(TotSBWEmergepcap=TotSBWEmerge/TotSBW, TotALTEmergepcap=TotALTEmerge/TotALT, ratioSBWEmerge=TotSBWEmerge/TotALTEmerge,logratioSBWEmerge=log10(ratioSBWEmerge),ratioSBWEmergepcap=TotSBWEmergepcap/TotALTEmergepcap,logratioSBWEmergepcap=log10(ratioSBWEmergepcap)) %>%
  mutate(dummyPeak=Peak+4, Plot=as.factor(Plot)) %>% #creates a dummy variable because function gls can not use negative numbers.
  mutate(bpae=case_when(
    Peak %in% c(-3,-2) ~ "b",
    Peak %in% c(-1,0,1) ~ "p",
    Peak %in% c(2,3) ~ "a",
    Peak %in% c(4,5,6,7,8,9,10) ~ "en"
  ))

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

plot(SBWALTdietswitchfulllinear, add.smooth = FALSE, which = 1, col=communitydietswitch$Plot)
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

# Interaction is significant (p-value = 0.0033, delta AIC = 7.42856). Therefore, keep all fixed effects in the model including the interaction.


#7. Refit the model found in step 6. with REML estimation
SBWALTdietswitchreml<-gls(logratioSBWEmerge~logratioSBWpeakplot*Plot,data=communitydietswitch, method="REML")

summary(SBWALTdietswitchreml)

#Create graph of log10 ratio parasitoid emergence from budworm and other caterpillars to the log10 ratio of total budworm to other caterpillars
SBWALTdietswitchremlGLS<-Gls(logratioSBWEmerge~logratioSBWpeakplot*Plot,data=communitydietswitch, method="REML") #using the function Gls because allows Predict function needed to plot trendline

dietswitchpredict<-Predict(SBWALTdietswitchremlGLS, logratioSBWpeakplot=seq(from=-1.5,to=2.1,by=0.1),Plot=seq(from=1,to=3,by=1))%>%
  mutate(Plot=as.factor(Plot))

comswitchlogratioplot<-ggplot()+
  geom_abline(intercept=0,slope=1,linetype="dotdash")+
  geom_abline(intercept=0,slope=0)+
  geom_vline(xintercept=0)+
  geom_ribbon(data=filter(dietswitchpredict, Plot==3), aes(x=logratioSBWpeakplot, ymin=lower, ymax=upper), linetype = 1, size=0.05, alpha=0.3, color = "black")+
  geom_ribbon(data=filter(dietswitchpredict, Plot==2), aes(x=logratioSBWpeakplot, ymin=lower, ymax=upper), linetype = 1, size=0.05, alpha=0.2, color = "black")+
  geom_ribbon(data=filter(dietswitchpredict, Plot==1), aes(x=logratioSBWpeakplot, ymin=lower, ymax=upper), linetype = 1, size=0.05, alpha=0.1, color = "black")+
  geom_line(data=dietswitchpredict, aes(logratioSBWpeakplot,yhat,linetype=Plot),size=1.5)+
  geom_point(data = communitydietswitch, aes(logratioSBWpeakplot, logratioSBWEmerge, shape=Plot, colour=bpae),size=5)+
  theme(axis.title.y=element_text(hjust=0.5, vjust=1.5), 
        legend.text=element_text(size=14),legend.justification=c(1,0), legend.position=c(0.98,0.05),legend.box = "horizontal")+
  coord_fixed(ratio = 1)+
    scale_shape_manual(values=c(16,17,15))+scale_linetype_manual(values=c("solid","dashed","dotted"))+scale_color_manual(name="Peak", breaks=c("b","p","a", "en"),labels=c("before","during","after", "endemic"), values=c("#39568CFF","#B8DE29FF", "#440154FF", "#1F968BFF"))+
  ylab("Relative Budworm Utilization")+xlab("Relative Budworm Frequency")


ggsave("figs/communitydietswitchlogratio.pdf", plot=comswitchlogratioplot, width=8,height=8) #Figure 1

# Test slope and intercept for each plot
plot1 <- lm(logratioSBWEmerge~logratioSBWpeakplot, data = filter(communitydietswitch, Plot == "1"))
plot1_slopet <- (coef(summary(plot1))[2,1]-1)/coef(summary(plot1))[2,2]
2*pt(abs(plot1_slopet),df=5,lower.tail = FALSE)
slope_plot1_allsp <- coef(summary(plot1))[2,1]
int_plot1_allsp <- coef(summary(plot1))[1,1]

plot2 <- lm(logratioSBWEmerge~logratioSBWpeakplot, data = filter(communitydietswitch, Plot == "2"))
plot2_slopet <- (coef(summary(plot2))[2,1]-1)/coef(summary(plot2))[2,2]
2*pt(abs(plot2_slopet),df=5,lower.tail = FALSE)
slope_plot2_allsp <- coef(summary(plot2))[2,1]
int_plot2_allsp <- coef(summary(plot2))[1,1]

plot3 <- lm(logratioSBWEmerge~logratioSBWpeakplot, data = filter(communitydietswitch, Plot == "3"))
plot3_slopet <- (coef(summary(plot3))[2,1]-1)/coef(summary(plot3))[2,2]
2*pt(abs(plot3_slopet),df=5,lower.tail = FALSE)
slope_plot3_allsp <- coef(summary(plot3))[2,1]
int_plot3_allsp <- coef(summary(plot3))[1,1]


# Underlysing causes of parasitoid community host preference

numinteractionsspecies<-interaction%>%
  filter(SpecID.high %in% paraSBW)%>%
  group_by(SpecID.high)%>%
  summarise(numinter=length(SpecID.high))%>%
  arrange(desc(numinter)) # Produces a dataframe of the total number of emergences with BOTH budworm and other caterpillars for each parasitoid taxon. Then arranged descending from most emergences to lowest emergences.

paralist<-numinteractionsspecies$SpecID.high

#The following code removes in turn the parastiod taxa that emerged the most from budworm and other caterpillars, the two taxa that emerged the most, and the three taxa that emerged the most. For each of these subsets of data, the same model as above (parasitoid community host preference) is ran and slops and intercepts calculated for each plot.

slope_tval_plot <- function(plot, origslope){
  (coef(summary(plot))[2,1]-origslope)/coef(summary(plot))[2,2]
}

int_tval_plot <- function(plot, origint){
  (coef(summary(plot))[1,1]-origint)/coef(summary(plot))[1,2]
}

slope_pval_plot <- function(plot, origslope){
  2*pt(abs(slope_tval_plot(plot, origslope)),df=5,lower.tail = FALSE)
}

int_pval_plot <- function(plot, origint){
  2*pt(abs(int_tval_plot(plot, origint)),df=5,lower.tail = FALSE)
}

drop_top3_sp <- function(x){
  communitydietswitchlog_dropsp<-interaction%>%
    filter(SpecID.high %in% paraSBW)%>%
    mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
    filter(!SpecID.high %in% x)%>%
    group_by(Peak, Plot)%>%
    summarise(TotSBWEmerge=length(SpecID.high[SBWALT=="SBW"]),TotALTEmerge=length(SpecID.high[SBWALT=="ALT"]))%>%
    full_join(ratioSBWbypeakplot,by=c("Peak","Plot"))%>%
    mutate(ratioSBWEmerge=TotSBWEmerge/TotALTEmerge,logratioSBWEmerge=log10(ratioSBWEmerge)) %>%
    filter(TotSBW!=0,TotALT!=0,TotSBWEmerge!=0,TotALTEmerge!=0)%>%
    mutate(dummyPeak=Peak+4)

    plot1dropsp <- lm(logratioSBWEmerge~logratioSBWpeakplot, data = filter(communitydietswitchlog_dropsp, Plot == "1"))
    plot2dropsp <- lm(logratioSBWEmerge~logratioSBWpeakplot, data = filter(communitydietswitchlog_dropsp, Plot == "2"))
    plot3dropsp <- lm(logratioSBWEmerge~logratioSBWpeakplot, data = filter(communitydietswitchlog_dropsp, Plot == "3"))

  ttable <-tibble(
    speciesdrop = rep(paste(x,collapse=" "), times=3),
    plot = c("Plot 1", "Plot 2", "Plot 3"),
    slope = c(coef(summary(plot1dropsp))[2,1], coef(summary(plot2dropsp))[2,1], coef(summary(plot3dropsp))[2,1]),
    slopet = c(slope_tval_plot(plot1dropsp, slope_plot1_allsp), slope_tval_plot(plot2dropsp, slope_plot2_allsp), slope_tval_plot(plot3dropsp, slope_plot3_allsp)),
    slopep = c(slope_pval_plot(plot1dropsp, slope_plot1_allsp), slope_pval_plot(plot2dropsp, slope_plot1_allsp), slope_pval_plot(plot3dropsp, slope_plot1_allsp)),
    int = c(coef(summary(plot1dropsp))[1,1], coef(summary(plot2dropsp))[1,1], coef(summary(plot3dropsp))[1,1]),
    intt = c(int_tval_plot(plot1dropsp, int_plot1_allsp), int_tval_plot(plot2dropsp, int_plot2_allsp), int_tval_plot(plot3dropsp, int_plot3_allsp)),
    intp = c(int_pval_plot(plot1dropsp, int_plot1_allsp), int_pval_plot(plot2dropsp, int_plot1_allsp), int_pval_plot(plot3dropsp, int_plot1_allsp)),
    df = c(plot1dropsp$df.residual, plot2dropsp$df.residual, plot3dropsp$df.residual )
    )
  return(ttable)
}

bind_rows(drop_top3_sp("p01"), drop_top3_sp(c("p01","p02")),drop_top3_sp(c("p01","p02","p08"))) %>% #drop the most abundant parasitoid taxon (p01), then drop the two most abundant taxa, then drop the three most abundant taxa. compare coefficient values, run one sample ttests and combine into one table.
  data.frame()
#---

#nMDS analysis of parasitoid community over time
turnoverprep<-interaction%>%
  filter(SpecID.high %in% paraSBW)%>%
  group_by(Peak,Plot)%>%
  summarise(Totemerge=length(SpecID.high))#creates a dataframe of the total emergences of all parasitoid taxa combined for each peak and plot (used to standardise the number of emergences per taxon by the total number of emergences because different years had different total number of emergences)

turnover<-interaction %>%
  filter(SpecID.high %in% paraSBW)%>%
  group_by(Peak,Plot, SpecID.high) %>%
  summarise(nemerge=length(SpecID.high)) %>%
  ungroup() %>%
  left_join(turnoverprep,by=c("Peak","Plot")) %>%
  mutate(nemergestand=nemerge/Totemerge, bpae=case_when(
    Peak %in% c(-3,-2) ~ "b",
    Peak %in% c(-1,0,1) ~ "p",
    Peak %in% c(2,3) ~ "a",
    Peak %in% c(4,5,6,7,8,9,10) ~ "en"
  )) %>%
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

turnovercommatfin <- turnovercommat %>%
  select(-c(Peak,Plot,bpae))#removes the factors from the dataframe to create a true community matrix

turnovercommatfin[is.na(turnovercommatfin)]<-0#replaces any NA values with 0. these NA values were created when the long format data was spread into a wide format dataframe and where taxa were not found for a certain Peak or Plot.

dimcheckMDS(turnovercommatfin, distance = "bray", k = 6, trymax = 30, autotransform=TRUE) #after exmining screeplot from 6 dimensions to 1 dimension, chosen to use 2 dimensions because additional dimensions provide small reductions in stress

set.seed(1223)
turnoverNMDS <- metaMDS(turnovercommatfin, distance = "bray", k = 2, trymax = 100, autotransform=TRUE)

stressplot(turnoverNMDS)  ##to test if distances as represented in ordination are correlated with actual distances

stress10<-c(0.0868891,0.1270212,0.08939534,0.08462844,0.1265077,0.08547444,0.09038554,0.08890437,0.1225755,0.1265186)
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
  coord_fixed(ratio = 1)+scale_color_manual(name="Peak", breaks=c("b","p","a", "en"),labels=c("before","during","after", "endemic"), values=c("#39568CFF","#B8DE29FF", "#440154FF", "#1F968BFF"))+scale_shape_manual(values=c(16,17,15))+scale_x_continuous(breaks=c(-1,0,1))

ggsave("figs/nmdsturnover.pdf",plot=turnoverNMDSplot,width=10,height=10) #Figure 2


# Testing the separation of the parasitoid communities between before, during, and after the budworm peak, using the adonis function for a permanova. - Guide to formulating the permutation structure and permutation test from the following webpage http://thebiobucket.blogspot.ca/2011/04/repeat-measure-adonis-lately-i-had-to.html#more

fit <- adonis(vegdist(turnovercommatfin)~bpae, factors, permutations=1)

# number of perms
B <- 999

### setting up a dataframe which will be populated by random F Model values:
pop <- rep(NA, B + 1)

# the first entry of the dataframe pop will be the true F. model:
pop[1] <- fit$aov.tab[1, 4]

# set up a the permutation structure where permutations occur withing each plot and where the order of years are maintained for the permutations:
ctrl <- how(within = Within(type = "series",mirror = FALSE), plots = Plots(strata = factors$Plot))

# Number of observations:
nobs <- nrow(turnovercommatfin)

# loop to populate the dataframe pop with random F model values:
set.seed(123)
for(i in 2:(B+1)){
  idx <- shuffle(nobs, control = ctrl)
  fit.rand <- adonis(vegdist(turnovercommatfin) ~ bpae[idx], factors,permutations = 1)
  pop[i] <- fit.rand$aov.tab[1, 4]
}

# Return the p-value of the full permanova:
print(pval <- sum(pop >= pop[1]) / (B + 1))


# Food web Topology and per capita interaction strengths---------

speciesnamesprepbi<-speciesnamesSBW %>%
  rename(dataID=SpecID.high,SpecID.high=ParasitoidLatinName)%>%
  mutate(SpecID.high=as.factor(SpecID.high))%>%
  left_join(speciesnamesprep,by="SpecID.high")%>%
  select(-SpecID.high)%>%
  rename(SpecID.high=dataID,SpecID.highnum=numID) # creates a dataframe that includes the parasitoid code from data collection and the corresponding number for the figure.

bipartitegraphSBWALT<-interaction%>%
  mutate(Peak=as.character(Peak))%>%
  mutate(Peak=gsub("-","minus",Peak))%>%
  filter(SpecID.high %in% paraSBW)%>%
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
plotweb(as.matrix(bipartitegraphSBWALT$'minus3'),ybig=1.6,low.spacing=0.7,high.spacing=0.075,x.lim=c(0,2.6),labsize=2.2,high.lab.dis=0.12,low.xoff = 0.6,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

#0
pdf("figs/0plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'0'),ybig=1.6,low.xoff = 1.6,low.spacing=1.6,high.spacing=0.147,x.lim=c(0,5.3),labsize=2.2,high.lab.dis=0.12,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

#3
pdf("figs/3plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'3'),ybig=1.6,low.xoff = 1,low.spacing=1.8,high.spacing=0.18,x.lim=c(0,5.61),labsize=2.2,high.lab.dis=0.12,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

#10
pdf("figs/10plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'10'),ybig=1.6,high.spacing=0.1,x.lim=c(0,1.2),labsize=2.2,high.lab.dis=0.12,low.xoff = 0.1,col.low="#C73D73FF", col.high="#160F3BFF",col.interaction="#FECC8FFF", method='normal')
dev.off()


#calculation of median interactions strength to maximum interaction strength
#don't need to calculate per captita emergence because ratio of median per capita emergence to max per capita emergence would cancel out denominator of number of SBW or ALT - so same as using number of emergences
bipartitegraphSBWALTplotprep<-interaction%>%
  mutate(Peak=as.character(Peak))%>%
  mutate(Peak=gsub("-","minus",Peak))%>%
  filter(SpecID.high %in% paraSBW)%>%
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
  if ("ALT" %in% colnames(y)) {
  mutate(y,ALT=ifelse(ALT==0,NA,ALT),SBW=ifelse(SBW==0,NA,SBW))%>%
      summarise(medALT=median(ALT, na.rm=TRUE),maxALT=max(ALT,na.rm=TRUE),medSBW=median(SBW, na.rm=TRUE),maxSBW=max(SBW, na.rm=TRUE))%>%
      mutate(medmaxratioALT=medALT/maxALT,medmaxratioSBW=medSBW/maxSBW)
  } else {
  mutate(y,SBW=ifelse(SBW==0,NA,SBW))%>%
      summarise(medSBW=median(SBW, na.rm=TRUE),maxSBW=max(SBW, na.rm=TRUE))%>%
      mutate(medmaxratioSBW=medSBW/maxSBW)
    }
  } )%>%
  bind_rows()%>%
  mutate(Year=names(plot1bipartite), Plot=1)# creates a dataframe of median to maximum interaction strength (number of emergences is a proxy for interaction strength) for parasitoids emerging from budworm (SBW) and other caterpillars (ALT).

intdistplot2<-plot2bipartite%>%
  lapply(function(x) {y<-data.frame(t(x))
  if ("SBW" %in% colnames(y)) { 
    mutate(y,ALT=ifelse(ALT==0,NA,ALT),SBW=ifelse(SBW==0,NA,SBW))%>%
      summarise(medALT=median(ALT, na.rm=TRUE),maxALT=max(ALT,na.rm=TRUE),medSBW=median(SBW, na.rm=TRUE),maxSBW=max(SBW, na.rm=TRUE))%>%
      mutate(medmaxratioALT=medALT/maxALT,medmaxratioSBW=medSBW/maxSBW)
  } else {
    mutate(y,ALT=ifelse(ALT==0,NA,ALT))%>%
      summarise(medALT=median(ALT, na.rm=TRUE),maxALT=max(ALT,na.rm=TRUE))%>%
      mutate(medmaxratioALT=medALT/maxALT)
    }
  } )%>%
  bind_rows()%>%
  mutate(Year=names(plot2bipartite),Plot=2)# creates a dataframe of median to maximum interaction strength (number of emergences is a proxy for interaction strength) for parasitoids emerging from budworm (SBW) and other caterpillars (ALT).


intdistplot3<-plot3bipartite%>%
  lapply(function(x) {y<-data.frame(t(x))
  if ("SBW" %in% colnames(y)) {
    mutate(y,ALT=ifelse(ALT==0,NA,ALT),SBW=ifelse(SBW==0,NA,SBW))%>%
      summarise(medALT=median(ALT, na.rm=TRUE),maxALT=max(ALT,na.rm=TRUE),medSBW=median(SBW, na.rm=TRUE),maxSBW=max(SBW, na.rm=TRUE))%>%
      mutate(medmaxratioALT=medALT/maxALT,medmaxratioSBW=medSBW/maxSBW)
  } else {
    mutate(y,ALT=ifelse(ALT==0,NA,ALT))%>%
      summarise(medALT=median(ALT, na.rm=TRUE),maxALT=max(ALT,na.rm=TRUE))%>%
      mutate(medmaxratioALT=medALT/maxALT)
    }
  } )%>%
  bind_rows()%>%
  mutate(Year=names(plot3bipartite),Plot=3)# creates a dataframe of median to maximum interaction strength (number of emergences is a proxy for interaction strength) for parasitoids emerging from budworm (SBW) and other caterpillars (ALT).

intdistall <- bind_rows(intdistplot1,intdistplot2,intdistplot3) %>%
  select(medmaxratioALT,medmaxratioSBW,Year, Plot) %>%
  rename(ALT=medmaxratioALT,SBW=medmaxratioSBW) %>%
  gather(SBWALT, medmax, -Year, -Plot) %>%
  mutate(Plot=as.factor(Plot), SBWALT = as.factor(SBWALT))

intdistall$medmax[intdistall$Year == 3 & intdistall$Plot == 3 & intdistall$SBWALT == "SBW"] <- NA
intdistall$medmax[intdistall$Year == 4 & intdistall$Plot == 2 & intdistall$SBWALT == "SBW"] <- NA
intdistall$medmax[intdistall$Year == 6 & intdistall$Plot == 2 & intdistall$SBWALT == "SBW"] <- NA
intdistall$medmax[intdistall$Year == 7 & intdistall$Plot == 2 & intdistall$SBWALT == "SBW"] <- NA
intdistall$medmax[intdistall$Year == 3 & intdistall$Plot == 3 & intdistall$SBWALT == "ALT"] <- NA
intdistall$medmax[intdistall$Year == 10 & intdistall$Plot == 2 & intdistall$SBWALT == "ALT"] <- NA
#SBW remove 3.3  4.2 6.2 7.2 (8.2 9.2 10.2 already NA) - removes med:max data values where SBW numbers less than 50
#ALT remove 3.3 10.2 - removes med:max data values where ALT numbers less than 50

intdistall$Year<-gsub("minus","-",intdistall$Year)

intdistall %<>% mutate(Year = as.numeric(Year), dummyPeak = Year+4, dummyPeak2 = dummyPeak^2) %>%
  filter(!is.na(medmax))

intdistlmSBW <- lm(medmax ~ dummyPeak+dummyPeak2, data = filter(intdistall, SBWALT == "SBW"))
intdistlmALT <- lm(medmax ~ dummyPeak+dummyPeak2, data = filter(intdistall, SBWALT == "ALT"))
summary(intdistlmSBW)
summary(intdistlmALT)

intdistplot<-ggplot(intdistall)+
  geom_vline(xintercept=0)+
  geom_point(aes(as.numeric(Year),medmax, colour=SBWALT, shape=Plot), size=5)+
  xlab("Years before/after peak")+ylab("Median:maximum \ninteraction strength")+
  theme(legend.text=element_text(size=14),
        legend.title=element_text(size=15),
        legend.justification=c(0,1),
        legend.position=c(0.01,0.99))+
    scale_color_manual(name="caterpillar type", breaks=c("SBW","ALT"),labels=c("budworm","other caterpillars"), values=c("#FA7876FF", "#C73D73FF"))+ 
    scale_x_continuous(breaks=c(-2,0,2,4,6,8,10))+scale_shape_manual(values=c(16,17,15))

ggsave("figs/intdistplot.pdf",intdistplot,width=8,height=7) #Figure 4


# Supporting Information -----------------------------------------------

# Graph of per capita emergences from budworm and other caterpillars compared to ratio of budworm and other caterpillars
#Apanteles fumiferanae
communitydietswitchlog_p01<-interaction%>%
  filter(SpecID.high %in% paraSBW)%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  filter(SpecID.high=="p01")%>%
  group_by(Peak, Plot)%>%
  summarise(TotSBWEmerge=length(SpecID.high[SBWALT=="SBW"]),TotALTEmerge=length(SpecID.high[SBWALT=="ALT"]))%>%
  full_join(ratioSBWbypeakplot,by=c("Peak","Plot"))%>%
  mutate(ratioSBWEmerge=TotSBWEmerge/TotALTEmerge,logratioSBWEmerge=log10(ratioSBWEmerge)) %>%
  filter(TotSBW!=0,TotALT!=0,TotSBWEmerge!=0,TotALTEmerge!=0) %>%
  mutate(bpae=case_when(
    Peak %in% c(-3,-2) ~ "b",
    Peak %in% c(-1,0,1) ~ "p",
    Peak %in% c(2,3) ~ "a",
    Peak %in% c(4,5,6,7,8,9,10) ~ "en"
  )) %>%
  ggplot(aes(logratioSBWpeakplot, logratioSBWEmerge))+
  geom_abline(intercept=0,slope=1,linetype="twodash")+
  geom_abline(intercept=0,slope=0)+
  geom_vline(xintercept=0)+
  geom_point(aes(shape=Plot, colour=bpae),size=5)+
  theme(axis.title.y=element_text(hjust=0.5, vjust=1.5), 
        legend.text=element_text(size=14),legend.justification=c(1,0), legend.position=c(0.32,0.7),legend.box = "horizontal")+
  coord_fixed(ratio = 1)+
  scale_shape_manual(values=c(16,17,15))+scale_linetype_manual(values=c("dotted","dashed","solid"))+scale_color_manual(name="Peak", breaks=c("b","p","a", "en"),labels=c("before","during","after", "endemic"), values=c("#39568CFF","#B8DE29FF", "#440154FF", "#1F968BFF"))+
  ylab("Log10 emergence \nbudworm:other caterpillars")+xlab("Log10 abundance \nbudworm:other caterpillars")

ggsave("figs/communitydietswitchlogratio_p01.pdf", plot=communitydietswitchlog_p01, width=8,height=8) #Figure S1

#Glypta fumiferanae
communitydietswitchlog_p02<-interaction%>%
  filter(SpecID.high %in% paraSBW)%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  filter(SpecID.high=="p02")%>%
  group_by(Peak, Plot)%>%
  summarise(TotSBWEmerge=length(SpecID.high[SBWALT=="SBW"]),TotALTEmerge=length(SpecID.high[SBWALT=="ALT"]))%>%
  full_join(ratioSBWbypeakplot,by=c("Peak","Plot"))%>%
  mutate(ratioSBWEmerge=TotSBWEmerge/TotALTEmerge,logratioSBWEmerge=log10(ratioSBWEmerge)) %>%
  filter(TotSBW!=0,TotALT!=0,TotSBWEmerge!=0,TotALTEmerge!=0) %>%
  mutate(bpae=case_when(
    Peak %in% c(-3,-2) ~ "b",
    Peak %in% c(-1,0,1) ~ "p",
    Peak %in% c(2,3) ~ "a",
    Peak %in% c(4,5,6,7,8,9,10) ~ "en"
  )) %>%
  ggplot(aes(logratioSBWpeakplot, logratioSBWEmerge))+
  geom_abline(intercept=0,slope=1,linetype="twodash")+
  geom_abline(intercept=0,slope=0)+
  geom_vline(xintercept=0)+
  geom_point(aes(shape=Plot, colour=bpae),size=5)+
  theme(axis.title.y=element_text(hjust=0.5, vjust=1.5), 
        legend.text=element_text(size=14),legend.justification=c(1,0), legend.position=c(0.32,0.7),legend.box = "horizontal")+
  coord_fixed(ratio = 1)+
  scale_shape_manual(values=c(16,17,15))+scale_linetype_manual(values=c("dotted","dashed","solid"))+scale_color_manual(name="Peak", breaks=c("b","p","a", "en"),labels=c("before","during","after", "endemic"), values=c("#39568CFF","#B8DE29FF", "#440154FF", "#1F968BFF"))+
  ylab("Log10 emergence \nbudworm:other caterpillars")+xlab("Log10 abundance \nbudworm:other caterpillars")

ggsave("figs/communitydietswitchlogratio_p02.pdf", plot=communitydietswitchlog_p02, width=8,height=8) #Figure S2

#Smidtia fumiferanae
communitydietswitchlog_p08<-interaction%>%
  filter(SpecID.high %in% paraSBW)%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  filter(SpecID.high=="p08")%>%
  group_by(Peak, Plot)%>%
  summarise(TotSBWEmerge=length(SpecID.high[SBWALT=="SBW"]),TotALTEmerge=length(SpecID.high[SBWALT=="ALT"]))%>%
  full_join(ratioSBWbypeakplot,by=c("Peak","Plot"))%>%
  mutate(ratioSBWEmerge=TotSBWEmerge/TotALTEmerge,logratioSBWEmerge=log10(ratioSBWEmerge)) %>%
  filter(TotALT!=0,!is.na(TotSBWEmerge)) %>%
  mutate(bpae=case_when(
    Peak %in% c(-3,-2) ~ "b",
    Peak %in% c(-1,0,1) ~ "p",
    Peak %in% c(2,3) ~ "a",
    Peak %in% c(4,5,6,7,8,9,10) ~ "en"
  )) %>%
  ggplot(aes(logratioSBWpeakplot, log10(TotSBWEmerge)))+
  geom_vline(xintercept=0)+
  geom_point(aes(shape=Plot, colour=bpae),size=5)+
  theme(axis.title.y=element_text(hjust=0.5, vjust=1.5), 
        legend.text=element_text(size=14),legend.justification=c(1,0), legend.position=c(0.32,0.7),legend.box = "horizontal")+
  coord_fixed(ratio = 1)+
  scale_shape_manual(values=c(16,17,15))+scale_linetype_manual(values=c("dotted","dashed","solid"))+scale_color_manual(name="Peak", breaks=c("b","p","a", "en"),labels=c("before","during","after", "endemic"), values=c("#39568CFF","#B8DE29FF", "#440154FF", "#1F968BFF"))+
  ylab("Log10 emergence from budworm")+xlab("Log10 abundance \nbudworm:other caterpillars")

ggsave("figs/communitydietswitchlogratio_p08.pdf", plot=communitydietswitchlog_p08, width=8,height=8) #Figure S3


#Creation of bipartite graphs for extra years (Figures S4 and S5)
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
plotweb(as.matrix(bipartitegraphSBWALT$'minus2'),ybig=1.5,low.spacing=0.8,high.spacing=0.08,x.lim=c(0,3.1),labsize=2.2,high.lab.dis=0.12,low.xoff = 0.7,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/minus1plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'minus1'),ybig=1.5,low.spacing=1,high.spacing=0.075,x.lim=c(0,3.4),labsize=2.0,high.lab.dis=0.12,low.xoff = 0.8,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/1plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'1'),ybig=1.5,low.spacing=1,high.spacing=0.065,x.lim=c(0,3.5),labsize=2.2,high.lab.dis=0.12,low.xoff = 0.8,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/2plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'2'),ybig=1.5,low.spacing=0.7,high.spacing=0.075,x.lim=c(0,2.8),labsize=2.2,high.lab.dis=0.12,low.xoff = 0.6,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/4plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'4'),ybig=1.5,low.spacing=0.7,high.spacing=0.09,x.lim=c(0,2.5),labsize=2.2,high.lab.dis=0.12,low.xoff = 0.4,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/5plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'5'),ybig=1.5,low.spacing=0.5,high.spacing=0.1,x.lim=c(0,2),labsize=2.2,high.lab.dis=0.12,low.xoff = 0.2,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/6plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'6'),ybig=1.5,low.spacing=0.4,high.spacing=0.095,x.lim=c(0,1.9),labsize=2.2,high.lab.dis=0.12,low.xoff = 0.2,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/7plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'7'),ybig=1.5,low.spacing=0.5,high.spacing=0.1,x.lim=c(0,2.05),labsize=2.2,high.lab.dis=0.12,low.xoff = 0.1,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/8plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'8'),ybig=1.5,high.spacing=0.1,x.lim=c(0,1.35),labsize=2.2,high.lab.dis=0.12,low.xoff = 0.2,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

pdf("figs/9plot.pdf",width=7,height=6)
plotweb(as.matrix(bipartitegraphSBWALT$'9'),ybig=1.6,high.spacing=0.1,x.lim=c(0,1.16),labsize=2.2,high.lab.dis=0.12,low.xoff = 0.05,col.low=c("#C73D73FF","#2D1160FF"), col.high="#160F3BFF",col.interaction="#FECC8FFF")
dev.off()

#Calculate number of budworm and other caterpillars in each year and plot (Table S1) (needed in simulations below)
SBWALTabundpara <- full_join(ratioSBWparacommunity,ratioSBWbypeakplot,by=c("Peak","Plot")) %>%
  select(c(Peak, Plot, TotSBW, TotALT, TotSBWEmerge, TotALTEmerge)) %>%
  mutate(paraSBW = TotSBWEmerge/TotSBW, paraALT = TotALTEmerge/TotALT)

summary(SBWALTabundpara)

# Simulation of a community of parasitoids attacking two types of hosts indiscriminately using the sample numbers of budworm and other caterpillars (original code and idea by Brian Van Hezewijk)
paraSp <- c(1:51)
paraAttRates<-round(0.3*0.3^paraSp,4)+0.001 #generate a distribution of plausible attack rates. Rare parasitoids attack 1 in 500 hosts
survRate<-1-sum(paraAttRates[1:50])
paraAttRates[51]<-survRate # set the 51st attack rate to the survival rate
sum(paraAttRates)  # check to see if they all sum to unity

set.seed(12345)

medmax <- function (x) {
  medmaxdata <- tibble(
    Yr = x$Yr[1],
    med=median(x$Freq, na.rm=TRUE),
    max=max(x$Freq, na.rm=TRUE),
    medmax=med/max
  )
  return(medmaxdata)
}

simul_plot <- function(plot){
  totCollA<-round(SBWALTabundpara$TotSBW[SBWALTabundpara$Plot==plot],0) # mean number of budworm sampled each year
  years<- 1:length(totCollA)

  sampYear<-rep(years, times = totCollA)
  dat<-data.frame(sampYear=sampYear,para=0)

  # randomly assign parasitoid attacks to the hosts according to the distribution of species-specific attack rates
  for (i in years){
    pA<-sample(paraSp,totCollA[i],prob=paraAttRates,replace=TRUE)
    dat$para[dat$sampYear==i]<-pA
  }
  dat$para[dat$para==51] <- NA
  
  medmaxdata <- as.data.frame(table(dat$para,dat$sampYear)) %>%
    mutate(Freq = ifelse(Freq == 0, NA, Freq)) %>%
    rename(para = Var1, Yr = Var2) %>%
    mutate(Freq = as.numeric(Freq)) %>%
    select(Yr, Freq) %>%
    split(.$Yr) %>%
    map(medmax) %>%
    bind_rows() %>%
    mutate(Plot = plot)
  return(medmaxdata)
}


rep_coeff <- function() {
  intdistall_simul <- bind_rows(simul_plot(1),simul_plot(2),simul_plot(3)) %>%
    mutate(Yr = as.numeric(Yr), RelYear = ifelse(Plot %in% c(1,3), Yr-4, Yr)) %>%
    select(RelYear, Plot, medmax) %>%
    mutate(Plot=as.factor(Plot), dummyPeak = RelYear + 4, dummyPeak2 = dummyPeak^2)
  intdistall_simul$medmax[intdistall_simul$RelYear == 3 & intdistall_simul$Plot == 3] <- NA
  intdistall_simul$medmax[intdistall_simul$RelYear == 4 & intdistall_simul$Plot == 2] <- NA
  intdistall_simul$medmax[intdistall_simul$RelYear == 6 & intdistall_simul$Plot == 2] <- NA
  intdistall_simul$medmax[intdistall_simul$RelYear == 7 & intdistall_simul$Plot == 2] <- NA
  #SBW remove Year 3 plot 3,  Year 4 Plot 2, Year 6 Plot 2, Year 7 Plot 2 - removed same med:max datavalues in in real data
  intdistlm <- lm(medmax ~ dummyPeak+dummyPeak2, data = intdistall_simul)
  return(coefficients(intdistlm)[3])
}

simulateddata <- replicate(10000,rep_coeff())

simulateddata_plot <- ggplot(data.frame(simulateddata), aes(x=simulateddata)) +
  geom_histogram(color="black", fill="white") + 
  geom_vline(xintercept = 0.007557, size = 1, linetype = "dashed", colour = "black") +
  xlab(quote('Simulated relative'~year^2~' coefficient')) +ylab(quote(Count))

ggsave("figs/simulateddata_plot.png",simulateddata_plot,width=8,height=7) # Figure S6

exceed_count <- length(simulateddata[simulateddata <= 0.007557])

2 * exceed_count / 10000 #p-value (two tailed) - 0.8188

# plot dominance of weak interactions vs sampling year
intdistall_simul <- bind_rows(simul_plot(1),simul_plot(2),simul_plot(3)) %>%
  mutate(Yr = as.numeric(Yr), RelYear = ifelse(Plot %in% c(1,3), Yr-4, Yr)) %>%
  select(RelYear, Plot, medmax) %>%
  mutate(Plot=as.factor(Plot))

intdistall_simul$medmax[intdistall_simul$RelYear == 3 & intdistall_simul$Plot == 3] <- NA
intdistall_simul$medmax[intdistall_simul$RelYear == 4 & intdistall_simul$Plot == 2] <- NA
intdistall_simul$medmax[intdistall_simul$RelYear == 6 & intdistall_simul$Plot == 2] <- NA
intdistall_simul$medmax[intdistall_simul$RelYear == 7 & intdistall_simul$Plot == 2] <- NA

intdistsimulplot <- ggplot(intdistall_simul)+
  geom_vline(xintercept=0)+
  geom_point(aes(RelYear, medmax, shape = Plot), size=5)+
  xlab("Years before/after peak")+ylab("Median:maximum \ninteraction strength")+ylim(c(0,0.65))+
  theme(legend.text=element_text(size=14),
        legend.title=element_text(size=15),
        legend.justification=c(0,1),
        legend.position=c(0.01,0.99))+
  scale_x_continuous(breaks=c(-2,0,2,4), limits = c(-3,5))+scale_shape_manual(values=c(16,17,15))

ggsave("figs/intdistsimulplot.png",intdistsimulplot,width=8,height=6) # Figure S7

# Useful Extra Analyses ---------------------------------------------------
# Determine percentage of emergences for dropped unresolved taxa. Used in last sentence of Laboratory Work section of Methods
interactionundropcomp<-interactionimport%>%
  filter(!is.na(SpecID.high))  #20493

interactiondropcomp<-interaction%>%
  filter(!is.na(SpecID.high)) #18246

( length(interactionundropcomp$SpecID.high) - length(interactiondropcomp$SpecID.high) ) / length(interactionundropcomp$SpecID.high) #0.1096472

#Proportion of parasitism of budworm by parasitoid taxa that emerge from both budworm and other caterpillars
allsharedprop<-interaction%>%
  filter(!is.na(SpecID.high))%>%
  mutate(SBWALT=ifelse(SpecID.low=="h01","SBW","ALT"))%>%
  group_by(SBWALT)%>%
  summarise(NSha=length(SpecID.low[SpecID.high %in% paraSBW]),NUnSha=length(SpecID.low[!SpecID.high %in% paraSBW]),total=length(SpecID.low))%>%
  mutate(PropSh=NSha/total)#100% of all parasitism emergences from budworm were done by taxa that share budworm and other caterpillars - 81% of all parasitism emergences from other caterpillars were done by shared parasitoids. Info used in first paragraph of Statistical Analyses of Methods section

#Chao analysis of interactions - sampling effort
# Asking what percentage of potential (Chao) interactions between the parasitoids and budworm or other caterpillars were sampled. Used in Statistical Analses (Methods section) paragraph 1.

chaocalc <- interaction %>%
  filter(SpecID.high %in% paraSBW)%>%
  select(SampleID,SpecID.low,SpecID.high)%>%
  unite(interaction, SpecID.low,SpecID.high, sep="_" )%>%
  group_by(SampleID,interaction)%>%
  summarise(numinteraction=length(interaction))

SBWALTcode<-interaction %>%
  filter(SpecID.high %in% paraSBW)%>%
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

# Calculate the average abundance in each period (included in Figure 2 caption)
SBWabund<-interaction %>%
  filter(SpecID.low=="h01") %>%
  group_by(Peak, Plot) %>%
  summarise(TotSBW = length(SpecID.low)) %>%
  mutate(bpae=case_when(
    Peak %in% c(-3,-2) ~ "b",
    Peak %in% c(-1,0,1) ~ "p",
    Peak %in% c(2,3) ~ "a",
    Peak %in% c(4,5,6,7,8,9,10) ~ "en"
  )) %>%
  ungroup() %>%
  rbind(tibble(Peak=10,Plot=as.factor(2),TotSBW=0,bpae = "en")) %>%
  group_by(bpae) %>%
  summarise(minSBW=min(TotSBW), maxSBW=max(TotSBW), avSBW=mean(TotSBW))

###Linear mixed effects model of log10 ratio PER CAPITA parasitoid emergence from budworm and other caterpillars to the log10 ratio of total budworm to other caterpillars. Using the diet switching method outlined by Greenwood, J. J. D., and R. A. Elton. 1979. Analysing experiments on frequency-dependent selection by predators. The Journal of Animal Ecology 48:721.

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

ggsave("figs/communitydietswitchlogratiopcap.pdf", plot=comswitchlogratioplotpcap, width=8,height=8)

# Reason why using three parasitoid taxa (used in second paragraph of Parasitoid community host preference section of Methods )
threepara <- numSBWALTbyspeciesdata %>%
  mutate(totemer = NSBW + NALT)

sum(threepara$totemer)
  # Notes -------------------------------------------------------------------

# Undetermined herbivores were removed from dataset to follow Eveleigh et al. 2007 SI Methods.
