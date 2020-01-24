#Biomass meta-analysis
library(metafor)
library(reshape2)
library(gridExtra)
library(lme4)
library(LMERConvenienceFunctions)
library(lmerTest)
library(lsmeans)
library(car)
library(cowplot)
library(vegan)
library(MASS)
library(rgl)
library(lubridate)
library(ggplot2)
library(dplyr)
library(meta)
library(effsize)
library(esc)
library(metaviz)
#library(metagear)


Exotics<-read.csv("exotics.csv", header=T) 

# #justification for use of log response ratio over Hedges'
HG.e <- escalc(n2i = Nc, n1i = Ne, m2i = Mc, m1i = Me,
                  sd2i = Sc, sd1i = Se, data = Exotics, measure = "SMD",
                  append = TRUE)

Exotics$Mc2<-ifelse(Exotics$Mc==0,1/2,Exotics$Mc)
Exotics$Sc2<-ifelse(Exotics$Sc==0,1/2,Exotics$Sc)
Exotics$Me2<-ifelse(Exotics$Me==0,1/2,Exotics$Me)
Exotics$Se2<-ifelse(Exotics$Se==0,1/2,Exotics$Se)

# View(Exotics)
Exotics.kgm <- escalc(n2i = Nc, n1i = Ne, m2i = Mc2, m1i = Me2, 
               sd2i = Sc2, sd1i = Se2, data = Exotics, measure = "ROM", 
               append = TRUE )

# write.csv(Exotics.kgm, file = "ROM.exotics.csv")

# 
# # Check normality of Hedges' vs. log response ratio
 par(mfrow=c(1,2))
# #
Hedges = HG.kgm$yi
# # #
h <- hist(Hedges, breaks = 10, density = 10,
 xlab = "Hedges' (g)")
xfit <- seq(min(Hedges,na.rm=TRUE), max(Hedges,na.rm=TRUE), length = 40)
yfit <- dnorm(xfit, mean = mean(Hedges,na.rm=TRUE), sd = sd(Hedges,na.rm=TRUE))
yfit <- yfit * diff(h$mids[1:2]) * length(Hedges)

lines(xfit, yfit, col = "black", lwd = 2)
# # #
LRR = Exotics.kgm$yi
# # #
r <- hist(LRR, breaks = 10, density = 10, xlab = "Log Response Ratio")
xfit.r <- seq(min(LRR), max(LRR), length = 40)
yfit.r <- dnorm(xfit.r, mean = mean(LRR), sd = sd(LRR))
yfit.r <- yfit.r * diff(r$mids[1:2]) * length(LRR)
# # #
lines(xfit.r, yfit.r, col = "black", lwd = 2)
# #

#random effects model

rand.var=list(~ 1|Experiment/Paper.ID./Case)
 
rma.random <- rma.mv(yi = yi, V = vi, random = rand.var, data = Exotics.kgm, method = "ML")
summary(rma.random)
ran.df<-as.data.frame(coef(summary(rma.random)))

par(mfrow=c(1,1))
funnel(rma.random)
forest(rma.random,slab=paste(Exotics$author, Exotics$year, sep=", "))

