setwd("~/Desktop/Research/Dissertation Code and Data/Dissertation-Chapter-3")
library(lme4)
library(car)
library(lmerTest)
library(emmeans)
library(FSA)
library(tidyr)
library(tidyverse)
library(arm)
library(lattice)
library(rv)
library(nnet)
library(reshape2)
library(MASS)
library(ggplot2)
st.err = function(x) {sd(x)/sqrt(length(x))}

## Whole animal metrics ----
data1<-read.csv("WholeAnimal2020Diss.csv")
data1.whole<-subset(data1, id!="2021" & id!="2035") # Removed due to giving birth
mean(data1.whole$CaptureMass)
st.err(data1.whole$CaptureMass)
mean(data1.whole$SVL)
st.err(data1.whole$SVL)
length(data1.whole$id)
aggregate(CaptureMass~Sex, length, data=data1.whole)

# Change in Mass
# Add grams per day measure
gpd.data<-read.csv("gpd.csv") 
data1.whole<-merge(data1.whole, gpd.data, by="id")

aggregate(PercentChangeperDay~treatment, length, data=data1.whole)
a<-aggregate(PercentChangeperDay~treatment, mean, data=data1.whole)
a.se<-aggregate(PercentChangeperDay~treatment, st.err, data=data1.whole)
a$SE<-a.se$PercentChangeperDay
colnames(a)<-c("treatment","Mean","SE")

m1=aov(PercentChangeperDay~treatment, data=data1.whole)
summary(m1)
TukeyHSD(m1)

summary(lm(g.day~ReleaseMass, data=data1.whole))
m2=aov(g.day~treatment, data=data1.whole)
summary(m2)
em.m2<-emmeans(m2, c("treatment"))
summary(em.m2)
contrast(em.m2, 'tukey')

# Tolerance acclimation
data1.whole.TA<-subset(data1.whole, id!="2001" & id!="2004") # 2001 was stressed before CTmax trial and 2004 was injured during tracking
# Change in CTmax
mean(data1.whole.TA$InitialMax)
st.err(data1.whole.TA$InitialMax)
mean(data1.whole.TA$FinalMax)
st.err(data1.whole.TA$FinalMax)
aggregate(ChangeInMax~treatment, length, data=data1.whole.TA)
aggregate(ChangeInMax~treatment, mean, data=data1.whole.TA)
aggregate(ChangeInMax~treatment, st.err, data=data1.whole.TA)
summary(aov(ChangeInMax~treatment, data=data1.whole.TA))
shapiro.test(data1.whole.TA$ChangeInMax)
bartlett.test(ChangeInMax~treatment, data=data1.whole.TA)

# Change in CTmin
mean(data1.whole.TA$InitialMin)
st.err(data1.whole.TA$InitialMin)
mean(data1.whole.TA$FinalMin)
st.err(data1.whole.TA$FinalMin)
aggregate(ChangeInMin~treatment, length, data=data1.whole.TA)
aggregate(ChangeInMin~treatment, mean, data=data1.whole.TA)
aggregate(ChangeInMin~treatment, length, data=data1.whole.TA)
summary(aov(ChangeInMin~treatment, data=data1.whole.TA))
shapiro.test(data1.whole.TA$ChangeInMin)
bartlett.test(ChangeInMin~treatment, data=data1.whole.TA)

# Home Range
aggregate(HR100~treatment, length, data=data1)
HR.m<-aggregate(HR100~treatment, mean, data=data1)
Hr.se<-aggregate(HR100~treatment, st.err, data=data1)
aggregate(HR100~treatment, length, data=data1)
colnames(Hr.se)[2]<-"SE"
HR.m$SE<-Hr.se$SE
colnames(HR.m)[2]<-"HR100"

shapiro.test(log(data1$HR100))
bartlett.test(HR100~treatment, data=data1)

aov.HR<-aov(log(HR100)~treatment, data=data1)
summary(aov.HR)
TukeyHSD(aov.HR)

# Movement
aggregate(SumDailyMovement~treatment, mean, data=data1)
aggregate(SumDailyMovement~treatment, st.err, data=data1)
aggregate(AvgDailyMovement~treatment, length, data=data1)

shapiro.test(log(data1$SumDailyMovement))
bartlett.test(SumDailyMovement~treatment, data=data1)
summary(aov(log(SumDailyMovement)~treatment, data=data1))

# DEE
data1.dee<-subset(data1.whole, treatment!="longitudinal")
data1.dee<-droplevels(data1.dee)

b<-aggregate(DEE~treatment, mean, data=data1.dee)
b.se<-aggregate(DEE~treatment, st.err, data=data1.dee)
b$SE<-b.se$DEE
colnames(b)<-c("treatment","Mean","SE")
aggregate(DEE~treatment, length, data=data1.dee)

t.test(DEE~treatment, data=data1.dee)

## Cort Data ----
CortData<-read.csv("Corticosterone2020Diss.csv")
aggregate(Cort~Time+Treatment, length, data=CortData)
aggregate(Cort~Time+Treatment, range, data=CortData)

shapiro.test(log(CortData$Cort))
bartlett.test(log(Cort)~Treatment, data=CortData)

cort.m<-aggregate(Cort~Time+Treatment, mean, data=CortData)
cort.se<-aggregate(Cort~Time+Treatment, st.err, data=CortData)
cort.m$SE<-cort.se$Cort

m3<-lmer(log(Cort)~Time+Treatment+(1|ID), data=CortData)
em.m3<-emmeans(m3, c("Time"))
summary(em.m3)
contrast(em.m3, 'tukey')
anova(m3, type=3)

## Environmental Temp Data ----
# Air temperatures from ibuttons
All<-read.csv("AirTemps2020Diss.csv")
All$Datetime<-as.POSIXct(strptime(All$Datetime, format("%Y-%m-%d %H:%M:%S")))
All.Day<-subset(All, format(Datetime, "%H:%M:%S")>"06:00:00" & format(Datetime, "%H:%M:%S")<"21:00:00")
All.Night1<-subset(All, format(Datetime, "%H:%M:%S")<"06:00:00")
All.Night2<-subset(All, format(Datetime, "%H:%M:%S")>"21:00:00")
All.Night<-rbind(All.Night1, All.Night2)

# Overall
aggregate(Temperature~Site, mean, data=All)
aggregate(Temperature~Site, sd, data=All)
aggregate(Temperature~Site, st.err, data=All)
All.int<-subset(All, Site=="Intersection")
All.high<-subset(All, Site=="High")
t.test(All.int$Temperature, All.high$Temperature)

# Daytime
aggregate(Temperature~Site, mean, data=All.Day)
aggregate(Temperature~Site, sd, data=All.Day)
aggregate(Temperature~Site, st.err, data=All.Day)
All.Day.int<-subset(All.Day, Site=="Intersection")
All.Day.high<-subset(All.Day, Site=="High")
t.test(All.Day.int$Temperature, All.Day.high$Temperature)

# Nighttime
aggregate(Temperature~Site, mean, data=All.Night)
aggregate(Temperature~Site, sd, data=All.Night)
aggregate(Temperature~Site, st.err, data=All.Night)
All.Night.int<-subset(All.Night, Site=="Intersection")
All.Night.high<-subset(All.Night, Site=="High")
t.test(All.Night.int$Temperature, All.Night.high$Temperature)

## Invert Data ----
# Prey Density (Count per trap)
High<-c(3,1,14,24,13,2)
Low<-c(6,13,1,6,5,6)
t.test(High, Low, var.equal=T)

# Shannon's Diversity
InvertData<-read.csv("Inverts2020Diss.csv")
InvertData$Order<-factor(InvertData$Order, levels = c("Coleoptera", "Hemiptera", "Orthoptera", "Arachnia", "Diptera", "Hymenoptera"))
InvertData$Abundance<-as.numeric(InvertData$Abundance)
InvertData$Site<-as.factor(InvertData$Site)
Inv.High<-subset(InvertData, Site=="Home")
Inv.Low<-subset(InvertData, Site=="Transplant")

# High
H.sum<-sum(Inv.High$Abundance)
abs(sum(log(Inv.High$Abundance/H.sum)*(Inv.High$Abundance/H.sum))) # 0.4764601
# Low
L.sum<-sum(Inv.Low$Abundance)
abs(sum(log(Inv.Low$Abundance/L.sum)*(Inv.Low$Abundance/L.sum))) # 0.9110963

## Light-level data ----
# Modified from Refsnider et al., 2018- https://github.com/songsqian/lizards
a.ll<-read.csv("Light2020Diss.csv")
a.ll.wide<-spread(a.ll, Bins, Lux)
a.ll.wide<-subset(a.ll.wide, ID!="2032" & ID!="2027")
a.ll.wide<-droplevels(a.ll.wide)

PID <- c(3:8)
multinom.better <- multinom(as.matrix(a.ll.wide[,PID])~Treatment,
                            data=a.ll.wide)
summary(multinom.better, corr=F)

New<-as.data.frame(table(a.ll.wide[,"Treatment"]))
colnames(New)<-"Treatment"

pp <- predict(multinom.better, type="probs",
              newdata=New,
              se.fit=T)

pp <- cbind(New, as.data.frame(pp))
pp <- pp[,-2]
pp.molten <- melt(pp, id.var="Treatment")

sim.multinom <- function(M, n.sims=NULL){
  ## M: a multinomial model object of class "multinom"
  ## n.sims: number of Monte Carlo somulations
  library(rv)
  ## a package for random variate simulation and calculation
  if (is.null(n.sims)) n.sims <- getnsims()
  else setnsims(n.sims)
  ## setting simulation numbers to be either user supplied
  ## or rv package default (2500)
  object.class <- class(M)
  if(object.class[1]!="multinom") stop ("Not a multinom object")
  
  summ <- summary(M)
  beta.hat <- as.vector(t(coef(M)))
  V.beta <- vcov(M)
  k <- length(beta.hat)
  beta <- array(NA, c(n.sims, k))
  lbs <- labels(coef(M))
  dmnm <- character()
  for (i in 1:length(lbs[[1]]))
    dmnm <- c(dmnm, paste(lbs[[1]][i], lbs[[2]], sep=":"))
  dimnames(beta) <- list(NULL, dmnm)
  beta <- mvrnorm(n.sims, beta.hat, V.beta)
  return(beta)
}

sim.Better <- rvsims(sim.multinom(multinom.better, 5000))

trc <- as.numeric(pp$Treatment=="Control")
trl <- as.numeric(pp$Treatment=="Longitudinal")
trt <- as.numeric(pp$Treatment=="Transplant")

X <- cbind(trc, trl, trt)

sim.Better <- rvmatrix(sim.Better, nrow=5, ncol=3, byrow=T)

Xb11 <- Reduce('+', X[1,]*sim.Better[1,])
Xb12 <- Reduce('+', X[1,]*sim.Better[2,])
Xb13 <- Reduce('+', X[1,]*sim.Better[3,])
Xb14 <- Reduce('+', X[1,]*sim.Better[4,])
Xb15 <- Reduce('+', X[1,]*sim.Better[5,])

Xb21 <- Reduce('+', X[2,]*sim.Better[1,])
Xb22 <- Reduce('+', X[2,]*sim.Better[2,])
Xb23 <- Reduce('+', X[2,]*sim.Better[3,])
Xb24 <- Reduce('+', X[2,]*sim.Better[4,])
Xb25 <- Reduce('+', X[2,]*sim.Better[5,])

Xb31 <- Reduce('+', X[3,]*sim.Better[1,])
Xb32 <- Reduce('+', X[3,]*sim.Better[2,])
Xb33 <- Reduce('+', X[3,]*sim.Better[3,])
Xb34 <- Reduce('+', X[3,]*sim.Better[4,])
Xb35 <- Reduce('+', X[3,]*sim.Better[5,])

pRV <- rvmatrix(0, nrow=3, ncol=6)

denomsum <- 1+exp(Xb11) + exp(Xb12) + exp(Xb13) +exp(Xb14) +exp(Xb15)
pRV[1,1] <- 1/denomsum
pRV[1,2] <- exp(Xb11)/denomsum
pRV[1,3] <- exp(Xb12)/denomsum
pRV[1,4] <- exp(Xb13)/denomsum
pRV[1,5] <- exp(Xb14)/denomsum
pRV[1,6] <- exp(Xb15)/denomsum

denomsum <- 1+exp(Xb21) + exp(Xb22) + exp(Xb23) +exp(Xb24) +exp(Xb25)
pRV[2,1] <- 1/denomsum
pRV[2,2] <- exp(Xb21)/denomsum
pRV[2,3] <- exp(Xb22)/denomsum
pRV[2,4] <- exp(Xb23)/denomsum
pRV[2,5] <- exp(Xb24)/denomsum
pRV[2,6] <- exp(Xb25)/denomsum

denomsum <- 1+exp(Xb31) + exp(Xb32) + exp(Xb33) + exp(Xb34) + exp(Xb35)
pRV[3,1] <- 1/denomsum
pRV[3,2] <- exp(Xb31)/denomsum
pRV[3,3] <- exp(Xb32)/denomsum
pRV[3,4] <- exp(Xb33)/denomsum
pRV[3,5] <- exp(Xb34)/denomsum
pRV[3,6] <- exp(Xb35)/denomsum

dataPlot <- as.data.frame(cbind(pp.molten, summary(pRV)))

#pdf(paste(plotDIR, "fittedUncas.pdf", sep="/"),
#    height=5, width=5.75)
my.panel <- function(x,y,subscripts, group.number, col, ...){
  myjitter <- c(-0.1,0.1,0)
  panel.dotplot(x,as.numeric(y)+myjitter[group.number],
                cex=0.5, col=col)
  ##    panel.grid()
  panel.segments(dataPlot[subscripts,7],
                 as.numeric(y)+myjitter[group.number],
                 dataPlot[subscripts,11],
                 as.numeric(y)+myjitter[group.number],
                 col=col)
  
}

## Figures ----

# Figure 1 was created in QGIS from data compiled in HomeRange.R

# Figure 2- Invertebrate breakdowns
InvertFig<-ggplot(InvertData, aes(fill=Order, y=Abundance, x=Site))+ 
  geom_bar(position="fill", stat="identity")+
  theme_bw()+
  theme(axis.title.x=element_blank())+
  theme(panel.grid.major=element_line(colour="#FFFFFF"),panel.grid.minor=element_line(colour="#FFFFFF"))+
  theme(axis.text=element_text(size=12,face="bold"), axis.title=element_text(size=14,face="bold"))
#ggsave("InvertBreakdown.png", width=5, height=5, plot=InvertFig)

# Figure 3- Change in mass per day
Fig3<-ggplot(data=data1.whole, aes(x=treatment, y=PercentChangeperDay))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitter(seed=2,width=0.15), color="gray")+
  geom_point(data=a, aes(x=treatment, y=Mean), size=4, shape=18)+
  theme_bw()+
  theme(axis.title.x=element_blank())+
  theme(panel.grid.major=element_line(colour="#FFFFFF"),panel.grid.minor=element_line(colour="#FFFFFF"))+
  theme(axis.text=element_text(size=12,face="bold"), axis.title=element_text(size=14,face="bold"))+
  theme(legend.position="")+
  scale_x_discrete(labels = c('Control','Longitudinal','Transplant'))+
  geom_hline(yintercept=0, linetype="dashed")+
  annotate("text",
           x = c(1, 2, 3),
           y = c(1.8, 1.8, -0.5),
           label = c("A", "AB", "B"), fontface="bold",
           size=6)+
  xlab("Treatment Group")+
  ylab("% Change in Mass Per Day")
#ggsave(path=path, filename="MassFig.jpeg", width=5, height=5, plot=Fig3)

# Figure 4- Daily energy expenditure- Longitudinal was removed due to low sample size
Fig4<-ggplot(data=data1.dee, aes(x=treatment, y=DEE))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitter(seed=2,width=0.15), color="gray")+
  geom_point(data=b, aes(x=treatment, y=Mean), size=4, shape=18)+
  theme_bw()+
  theme(axis.title.x=element_blank())+
  theme(panel.grid.major=element_line(colour="#FFFFFF"),panel.grid.minor=element_line(colour="#FFFFFF"))+
  theme(axis.text=element_text(size=12,face="bold"), axis.title=element_text(size=14,face="bold"))+
  theme(legend.position="")+
  scale_x_discrete(labels = c('Control', 'Transplant'))+
  annotate("text",
           x = c(1, 2),
           y = c(160, 108),
           label = c("A", "B"), fontface="bold",
           size=6)+
  xlab("Treatment Group")+
  ylab("DEE (J-g-d)")
#ggsave(path=path, filename="DEE.jpeg", width=5, height=5, plot=Fig4)

# Figure 5- Stress Response
Fig5<-ggplot(data=CortData, aes(x=Treatment, y=Cort, fill=Time))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(seed=2), color="gray")+
  scale_fill_grey(start=0.4, end=1)+
  geom_point(data=cort.m, aes(x=Treatment, y=Cort, shape=Time), size=4, shape=18, position=position_dodge(.75))+
  theme_bw()+
  theme(axis.title.x=element_blank())+
  theme(panel.grid.major=element_line(colour="#FFFFFF"),panel.grid.minor=element_line(colour="#FFFFFF"))+
  theme(axis.text=element_text(size=12,face="bold"), axis.title=element_text(size=14,face="bold"))+
  theme(legend.position="bottom")+
  xlab("")+
  ylab("Corticosterone (ng/mL)")+
  annotate("text",
           x = c(1, 2, 3),
           y = c(14, 25.5, 11.5),
           label = c("*", "*", "*"),
           size=10)
#ggsave(path=path, filename="Cort.jpeg", width=5, height=5, plot=Fig5)

# Figure 6- Home range size
Fig6<-ggplot(data=data1, aes(x=treatment, y=HR100))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitter(seed=2,width=0.15), color="gray")+
  geom_point(data=HR.m, aes(x=treatment, y=HR100), size=4, shape=18)+
  theme_bw()+
  theme(axis.title.x=element_blank())+
  theme(panel.grid.major=element_line(colour="#FFFFFF"),panel.grid.minor=element_line(colour="#FFFFFF"))+
  theme(axis.text=element_text(size=12,face="bold"), axis.title=element_text(size=14,face="bold"))+
  theme(legend.position="bottom")+
  xlab("")+
  ylab("Home Range (m^2)")+
  annotate("text",
           x = c(1, 2, 3),
           y = c(750,1700, 2200),
           label = c("A", "A", "A"), fontface="bold",
           size=6)
#ggsave(path=path, filename="HR100.jpeg", width=5, height=5, plot=Fig6)

# Figure 7- Proportion of time spent in full sun
# From simulations following Refsnider et al., 2018
dataPlot.1<-subset(dataPlot, variable=="full sun")
Fig7<-ggplot(data=dataPlot.1, aes(x=Treatment, y=mean, ymin=dataPlot.1$'2.5%', ymax=dataPlot.1$'97.5%'))+
  geom_point(size=4, position=position_dodge(width=0.5))+
  theme_bw()+
  theme(axis.title.x=element_blank())+
  theme(panel.grid.major=element_line(colour="#FFFFFF"),panel.grid.minor=element_line(colour="#FFFFFF"))+
  theme(axis.text=element_text(size=12,face="bold"), axis.title=element_text(size=14,face="bold"))+
  theme(legend.position="bottom")+
  geom_linerange(position = position_dodge(.5))+
  ylab("Proportion of Time in Full Sun")+
  annotate("text",
           x = c(1, 2, 3),
           y = c((dataPlot.1[1,11]+(0.002)), (dataPlot.1[2,11]+0.002), (dataPlot.1[3,11]+0.002)),
           label = c("A", "A", "B"),
           size=6)
#ggsave("Fig7.jpeg", width=5, height=5, plot=Fig7)

