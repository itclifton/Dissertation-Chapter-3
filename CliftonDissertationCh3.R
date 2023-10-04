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
library(cowplot)
st.err = function(x) {sd(x)/sqrt(length(x))}
path<-"/Users/ianclifton/Desktop/Research- PhD to Present/PhD/Chapters/Selection/Figures"


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
mean(data1.whole$PercentChangeperDay)
st.err(data1.whole$PercentChangeperDay)

m1=aov(PercentChangeperDay~treatment, data=data1.whole)
summary(m1)
TukeyHSD(m1)

# summary(lm(g.day~ReleaseMass, data=data1.whole))
# m2=aov(g.day~treatment, data=data1.whole)
# summary(m2)
# em.m2<-emmeans(m2, c("treatment"))
# summary(em.m2)
# contrast(em.m2, 'tukey')

# Tolerance acclimation
data1.whole.TA<-subset(data1.whole, id!="2001" & id!="2004") # 2001 was stressed before CTmax trial and 2004 was injured during tracking

# Change in CTmax
mean(data1.whole.TA$InitialMax)
st.err(data1.whole.TA$InitialMax)
summary(aov(InitialMax~treatment, data=data1.whole.TA))
mean(data1.whole.TA$FinalMax)
st.err(data1.whole.TA$FinalMax)
aggregate(ChangeInMax~treatment, length, data=data1.whole.TA)
aggregate(ChangeInMax~treatment, mean, data=data1.whole.TA)
aggregate(ChangeInMax~treatment, st.err, data=data1.whole.TA)
mean(data1.whole.TA$ChangeInMax)
st.err(data1.whole.TA$ChangeInMax)

summary(aov(ChangeInMax~treatment, data=data1.whole.TA))
shapiro.test(data1.whole.TA$ChangeInMax)
bartlett.test(ChangeInMax~treatment, data=data1.whole.TA)

# Change in CTmin
mean(data1.whole.TA$InitialMin)
st.err(data1.whole.TA$InitialMin)
summary(aov(InitialMin~treatment, data=data1.whole.TA))
mean(data1.whole.TA$FinalMin)
st.err(data1.whole.TA$FinalMin)
aggregate(ChangeInMin~treatment, length, data=data1.whole.TA)
aggregate(ChangeInMin~treatment, mean, data=data1.whole.TA)
aggregate(ChangeInMin~treatment, st.err, data=data1.whole.TA)
mean(data1.whole.TA$ChangeInMin)
st.err(data1.whole.TA$ChangeInMin)

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
mean(data1$HR100, na.rm=T)
sd(data1$HR100, na.rm=T)/sqrt(18)

plot(lm(HR100~treatment, data=data1))

aov.HR<-aov(HR100~treatment, data=data1)
summary(aov.HR)
TukeyHSD(aov.HR)

# Movement
aggregate(SumDailyMovement~treatment, mean, data=data1)
aggregate(SumDailyMovement~treatment, st.err, data=data1)
aggregate(SumDailyMovement~treatment, length, data=data1)
mean(data1$SumDailyMovement, na.rm=T)
sd(data1$SumDailyMovement, na.rm=T)/sqrt(18)

plot(lm(SumDailyMovement~treatment, data=data1))
summary(aov(SumDailyMovement~treatment, data=data1))

# DEE
data1.dee<-subset(data1.whole, treatment!="longitudinal")
data1.dee<-droplevels(data1.dee)

b<-aggregate(DEE~treatment, mean, data=data1.dee)
b.se<-aggregate(DEE~treatment, st.err, data=data1.dee)
b$SE<-b.se$DEE
colnames(b)<-c("treatment","Mean","SE")
aggregate(DEE~treatment, length, data=data1.dee)
mean(data1.dee$DEE, na.rm=T)
sd(data1.dee$DEE, na.rm=T)/sqrt(8)

t.test(DEE~treatment, data=data1.dee)

## Cort Data ----
CortData<-read.csv("Corticosterone2020Diss.csv")
aggregate(Cort~Time+Treatment, length, data=CortData)
aggregate(Cort~Time+Treatment, range, data=CortData)

cort.m<-aggregate(Cort~Time+Treatment, mean, data=CortData)
cort.se<-aggregate(Cort~Time+Treatment, st.err, data=CortData)
cort.m$SE<-cort.se$Cort

aggregate(Cort~Time, mean, data=CortData)
aggregate(Cort~Time, st.err, data=CortData)

m3<-lmer(Cort~Time+Treatment+(1|ID), data=CortData)
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
abs(sum(log(Inv.High$Abundance/H.sum)*(Inv.High$Abundance/H.sum))) # 0.4764601 Shannon's
mean(Inv.High$Abundance)
st.err(Inv.High$Abundance)

# Low
L.sum<-sum(Inv.Low$Abundance)
abs(sum(log(Inv.Low$Abundance/L.sum)*(Inv.Low$Abundance/L.sum))) # 0.9110963 Shannon's
mean(Inv.Low$Abundance)
st.err(Inv.Low$Abundance)

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
  geom_point(position=position_jitter(seed=2,width=0.15), color="#636363", size=2.5)+
  geom_point(data=a, aes(x=treatment, y=Mean), size=4, shape=18)+
  theme_classic()+
  theme(axis.title.x=element_blank())+
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
  geom_point(position=position_jitter(seed=2,width=0.15), color="#636363", size=2.5)+
  geom_point(data=b, aes(x=treatment, y=Mean), size=4, shape=18)+
  theme_classic()+
  theme(axis.title.x=element_blank())+
  theme(axis.text=element_text(size=12,face="bold"), axis.title=element_text(size=14,face="bold"))+
  theme(legend.position="")+
  scale_x_discrete(labels = c('Control', 'Transplant'))+
  annotate("text",
           x = c(1, 2),
           y = c(168, 115),
           label = c("A", "B"), fontface="bold",
           size=6)+
  xlab("Treatment Group")+
  ylab(bquote(bold('DEE (J' *~g^-1~day^-1*')')))
#ggsave(path=path, filename="DEE.jpeg", width=5, height=5, plot=Fig4)

# Figure 5- Stress Response
Fig5<-ggplot(data=CortData, aes(x=Treatment, y=Cort, fill=Time))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitterdodge(seed=2), color="#636363", size=2.5)+
  scale_fill_grey(start=0.4, end=1)+
  geom_point(data=cort.m, aes(x=Treatment, y=Cort, shape=Time), size=4, shape=18, position=position_dodge(.75))+
  theme_classic()+
  theme(axis.title.x=element_blank())+
  theme(axis.text=element_text(size=12,face="bold"), axis.title=element_text(size=14,face="bold"))+
  theme(legend.position=c(0.8,0.8))+
  xlab("")+
  ylab("Corticosterone (ng/mL)")+
  annotate("text",
           x = c(1, 2, 3),
           y = c(14, 25.5, 11.5),
           label = c("*", "*", "*"),
           size=10)+
  annotate("text",
           x = c(1, 2, 3),
           y = c(28, 32, 17),
           label = c("A", "A", "A"), fontface="bold",
           size=6)
#ggsave(path=path, filename="Cort.jpeg", width=5, height=5, plot=Fig5)

# Figure 6- Home range size
Fig6<-ggplot(data=data1, aes(x=treatment, y=HR100))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitter(seed=2,width=0.15), color="#636363", size=2.5)+
  geom_point(data=HR.m, aes(x=treatment, y=HR100), size=4, shape=18)+
  theme_classic()+
  theme(axis.title.x=element_blank())+
  theme(axis.text=element_text(size=12,face="bold"), axis.title=element_text(size=14,face="bold"))+
  theme(legend.position="bottom")+
  xlab("")+
  ylab("Home Range (m²)")+
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
  theme_classic()+
  theme(axis.title.x=element_blank())+
  theme(axis.text=element_text(size=12,face="bold"), axis.title=element_text(size=14,face="bold"))+
  theme(legend.position="bottom")+
  geom_linerange(position = position_dodge(.5))+
  ylab("Proportion of Time in Full Sun")+
  annotate("text",
           x = c(1, 2, 3),
           y = c((dataPlot.1[1,11]+(0.002)), (dataPlot.1[2,11]+0.002), (dataPlot.1[3,11]+0.002)),
           label = c("A", "A", "B"), fontface="bold",
           size=6)
#ggsave(path=path, "Fig7.jpeg", width=5, height=5, plot=Fig7)

# Physiological and Behavioral Panel Plot- New Figure 2
Panel=plot_grid(Fig3, Fig4, Fig5, Fig6, Fig7,
                       labels = "AUTO", ncol = 2, nrow=3, align="v")
#ggsave(path=path,"PanelFigure2.jpeg", width=10, height=15, plot=Panel)

## Reviewer Comments ----
aggregate(ReleaseMass~treatment, length, data=data1.whole)
aggregate(ReleaseMass~treatment, mean, data=data1.whole)
aggregate(ReleaseMass~treatment, st.err, data=data1.whole)

aggregate(RecoveredMass~treatment, length, data=data1.whole)
aggregate(RecoveredMass~treatment, mean, data=data1.whole)
aggregate(RecoveredMass~treatment, st.err, data=data1.whole)

aggregate(PercentChangeperDay~treatment, length, data=data1.whole)
aggregate(PercentChangeperDay~treatment, mean, data=data1.whole)
aggregate(PercentChangeperDay~treatment, st.err, data=data1.whole)

cort.1<-c(7.833,
          3.142,
          11.145,
          4.719,
          3.06,
          5.854,
          6.065,
          12.659,
          1.099,
          7.045,
          4.454,
          7.656,
          2.221,
          6.732,
          7.964,
          10.954)
time.1<-c("AM",
          "AM",
          "AM",
          "AM",
          "PM",
          "PM",
          "PM",
          "PM",
          "PM",
          "PM",
          "PM",
          'AM',
          "AM",
          "PM",
          "AM",
          'AM')
summary(aov(as.numeric(cort.1)~time.1, time.2))
aggregate(as.numeric(cort.1)~time.1, mean, data=time.2)

# Requested CORT Supplemental
SuppCORT<-read.csv("CORTSupplemental.csv")

summary(lm(CORT~Time, data=SuppCORT))
s1<-ggplot(aes(x=Time, y=CORT), data=SuppCORT)+
  geom_point(size=2)+
  geom_smooth(method="lm", se=F, color="black")+
  theme_classic()+
  theme(axis.text=element_text(size=20,face="bold"), axis.title=element_text(size=22,face="bold"),
        legend.text=element_text(size=22, face="bold"), legend.title=element_text(size=25, face="bold", hjust=0.5))+
  theme(axis.ticks.length.y=unit(.5, "cm"), axis.ticks.y=element_line(size=1.75), axis.line=element_line(size=1.75))+
  theme(legend.position = c(.75, .2), plot.margin = margin(11, 5.5, 5.5, 5.5, "pt"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 15), breaks=seq(0,15,5))+
  ylab("Corticosterone (ng/mL)")+
  xlab("Time to Draw (s)")

summary(lm(CORT~Tb, data=SuppCORT))
s2<-ggplot(aes(x=Tb, y=CORT), data=SuppCORT)+
  geom_point(size=2)+
  geom_smooth(method="lm", se=F, color="black")+
  theme_classic()+
  theme(axis.text=element_text(size=20,face="bold"), axis.title=element_text(size=22,face="bold"),
        legend.text=element_text(size=22, face="bold"), legend.title=element_text(size=25, face="bold", hjust=0.5))+
  theme(axis.ticks.length.y=unit(.5, "cm"), axis.ticks.y=element_line(size=1.75), axis.line=element_line(size=1.75))+
  theme(legend.position = c(.75, .2), plot.margin = margin(11, 5.5, 5.5, 5.5, "pt"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 15), breaks=seq(0,15,5))+
  ylab("Corticosterone (ng/mL)")+
  xlab("Body Temperature (°C)")

summary(aov(CORT~TimePeriod, data=SuppCORT))
s3<-ggplot(data=SuppCORT, aes(x=TimePeriod, y=CORT))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(position=position_jitter(seed=2,width=0.15), color="#636363", size=2.5)+
  theme_classic()+
  theme(axis.text=element_text(size=20,face="bold"), axis.title=element_text(size=22,face="bold"),
        legend.text=element_text(size=22, face="bold"), legend.title=element_text(size=25, face="bold", hjust=0.5))+
  theme(axis.ticks.length.y=unit(.5, "cm"), axis.ticks.y=element_line(size=1.75), axis.line=element_line(size=1.75))+
  theme(legend.position = c(.75, .2), plot.margin = margin(11, 5.5, 5.5, 5.5, "pt"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 15), breaks=seq(0,15,5))+
  theme(legend.position="")+
  scale_x_discrete(labels = c('AM','PM'))+
  annotate("text",
           x = c(1, 2),
           y = c(12.95, 12.95),
           label = c("A", "A"), fontface="bold",
           size=6)+
  xlab("Time of Day")+
  ylab("Corticosterone (ng/mL)")
PanelSupp=plot_grid(s1, s2, s3,
                labels = "AUTO", ncol = 3, nrow=1, align="v")
ggsave(path=path,"SuppFig1.jpeg", width=15, height=10, plot=PanelSupp)
