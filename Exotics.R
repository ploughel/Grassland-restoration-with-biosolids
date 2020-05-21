#Exotics meta-analysis
library(metafor)
library(tidyverse)
library(ggplot2)

library(car)
library(vegan)

library(meta)
library(effsize)
library(esc)
library(metaviz)


library(readxl)
library(xlsx)
library(rJava)
library("glmulti")

Exotics <- read_excel("~/Documents/GitHub/Grassland restoration with biosolids/Grassland-restoration-with-biosolids/data_sheetscorr.xlsx", 
                    sheet = "Exotics")%>%
  janitor::clean_names()%>%
  filter(measured_variable=="percent total")%>%
  mutate(mc_0=ifelse(mc ==0, mc+1/2 , mc),
         me_0=ifelse(mc ==0, me+1/2 , me),
         se.coef=round(sum(se, na.rm=TRUE)/sum(me, na.rm=TRUE),2), 
         sc.coef=round(sum(sc, na.rm = TRUE)/sum(mc, na.rm=TRUE),2),
         seimp=round(ifelse(error_in_the_article=="N", me*se.coef, se),2),
         scimp=round(ifelse(error_in_the_article=="N", mc*sc.coef, sc), 2))%>%
  filter(me_0!=0)




ROM.Exotics <- escalc(n2i = nc, n1i = ne, m2i = mc_0, m1i = me_0, 
                    sd2i = scimp, sd1i = seimp, data = Exotics, measure = "ROM", 
                    append = TRUE )

HG.Exotics <- escalc(n2i = nc, n1i = ne, m2i = mc, m1i = me, 
                   sd2i = scimp, sd1i = seimp, data = Exotics, measure = "SMD", 
                   append = TRUE )
dim(ROM.Exotics)


Hedges = HG.Exotics$yi
# #
h <- hist(Hedges, breaks = 10, density = 10,
          xlab = "Hedges' (g)")
xfit <- seq(min(Hedges,na.rm=TRUE), max(Hedges,na.rm=TRUE), length = 40)
yfit <- dnorm(xfit, mean = mean(Hedges,na.rm=TRUE), sd = sd(Hedges,na.rm=TRUE))
yfit <- yfit * diff(h$mids[1:2]) * length(Hedges)
# #
lines(xfit, yfit, col = "black", lwd = 2)
# #
LRR = ROM.Exotics$yi
# #
r <- hist(LRR, breaks = 10, density = 10, xlab = "Log Response Ratio")
xfit.r <- seq(min(LRR), max(LRR), length = 40)
yfit.r <- dnorm(xfit.r, mean = mean(LRR), sd = sd(LRR))
yfit.r <- yfit.r * diff(r$mids[1:2]) * length(LRR)
# #
lines(xfit.r, yfit.r, col = "black", lwd = 2)
#


#random effects model

rand.var=list(~ 1|paper_id/case)

rma.random <- rma.mv(yi = yi, V = vi, random = rand.var, data = ROM.Exotics, method = "ML")
summary(rma.random)
