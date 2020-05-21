#Biomass meta-analysis
library(metafor)
library(tidyverse)
library(ggplot2)

library(car)
library(vegan)

library(meta)
library(effsize)
library(esc)
library(metaviz)
library("glmulti")


library(readxl)

Biomass <- read_excel("~/Documents/GitHub/Grassland restoration with biosolids/Grassland-restoration-with-biosolids/data_sheetscorr.xlsx", 
                      sheet = "Biomass")
#View(Biomass)

ROM.Biomass <- escalc(n2i = Nc, n1i = Ne, m2i = Mc+1/2, m1i = Me, 
                      sd2i = scimp, sd1i = seimp, data = Biomass, measure = "ROM", 
                      append = TRUE )


HG.Biomass <- escalc(n2i = Nc, n1i = Ne, m2i = Mc, m1i = Me, 
                     sd2i = scimp, sd1i = seimp, data = Biomass, measure = "SMD", 
                     append = TRUE )



Hedges = HG.Biomass$yi
# #
h <- hist(Hedges, breaks = 10, density = 10,
          xlab = "Hedges' (g)")
xfit <- seq(min(Hedges,na.rm=TRUE), max(Hedges,na.rm=TRUE), length = 40)
yfit <- dnorm(xfit, mean = mean(Hedges,na.rm=TRUE), sd = sd(Hedges,na.rm=TRUE))
yfit <- yfit * diff(h$mids[1:2]) * length(Hedges)
# #
lines(xfit, yfit, col = "black", lwd = 2)
# #
LRR.0 = ROM.Biomass.0$yi
# #
r <- hist(LRR.0, breaks = 10, density = 10, xlab = "Log Response Ratio")
xfit.r <- seq(min(LRR), max(LRR), length = 40)
yfit.r <- dnorm(xfit.r, mean = mean(LRR), sd = sd(LRR))
yfit.r <- yfit.r * diff(r$mids[1:2]) * length(LRR)
# #
lines(xfit.r, yfit.r, col = "black", lwd = 2)


LRR = ROM.Biomass$yi
# #
r <- hist(LRR, breaks = 10, density = 10, xlab = "Log Response Ratio")
xfit.r <- seq(min(LRR), max(LRR), length = 40)
yfit.r <- dnorm(xfit.r, mean = mean(LRR), sd = sd(LRR))
yfit.r <- yfit.r * diff(r$mids[1:2]) * length(LRR)
# #
lines(xfit.r, yfit.r, col = "black", lwd = 2)
#random effects model

rand.var=list(~ 1|Paper.ID/Case)

rma.random <- rma.mv(yi = yi, V = vi, random = rand.var, data = ROM.Biomass, method = "ML")
summary(rma.random)


#-----------------------# model selection 


rma.glmulti<- function(formula, data, ...)
  rma(formula, vi, data=data, method ="ML",...)

names(ROM.Biomass)
ROM.ai.size <- glmulti(yi ~Biosolid.level..Mg.ha.1.+time+Temp+Precip+Mixture..Y.N.+Burn..Y.N.+
                        Seeded..Y.N.+Multiple.application..Y.N.+Severe.Disturbance+ai, data=ROM.Biomass, method="d",
                      level=1, crit="aicc", fitfunction = rma.glmulti)

ROM.ai.mod <- glmulti(yi ~ Biosolid.level..Mg.ha.1.+time+Temp+Precip+Mixture..Y.N.+Burn..Y.N.+
                       Seeded..Y.N.+Multiple.application..Y.N.+Severe.Disturbance+ai, data=ROM.Biomass, method="h",
                     level=1, crit="aicc", confsetsize=ROM.ai.size, fitfunction = rma.glmulti)

plot(ROM.ai.mod, type="s")




ROM.ai.rma<-rma.mv(yi, vi, mods = ~time+Burn..Y.N.+  Temp+ Severe.Disturbance+Mixture..Y.N.,
                       random = rand.var,data = ROM.Biomass,slab = paste(author, year), method="ML")


#models with interactions

ROM.prod.int1<-rma.mv(yi, vi, mods = ~time*Mixture..Y.N.,
                     random = rand.var,data = ROM.Biomass,slab = paste(author, year), method="ML")

ROM.prod.int2<-rma.mv(yi, vi, mods = ~time*Burn..Y.N.,
                      random = rand.var,data = ROM.Biomass,slab = paste(author, year), method="ML")

ROM.prod.int3<-rma.mv(yi, vi, mods = ~time*Severe.Disturbance,
                      random = rand.var,data = ROM.Biomass,slab = paste(author, year), method="ML")

ROM.prod.int4<-rma.mv(yi, vi, mods = ~Temp*Mixture..Y.N.,
                      random = rand.var,data = ROM.Biomass,slab = paste(author, year), method="ML")

ROM.prod.int5<-rma.mv(yi, vi, mods = ~Temp*Burn..Y.N.,
                       random = rand.var,data = ROM.Biomass,slab = paste(author, year), method="ML")

ROM.prod.int6<-rma.mv(yi, vi, mods = ~Temp*Severe.Disturbance,
                       random = rand.var,data = ROM.Biomass,slab = paste(author, year), method="ML")



mods = list(ROM.ai.rma, ROM.prod.int1,ROM.prod.int2,ROM.prod.int3,ROM.prod.int5,ROM.prod.int6)


# model fit diagnostics: pseudo R2, marginal/conditional R2, AICc, and I2
pseudo.r2 <- c()
marg.r2 <- c()
cond.r2 <- c()
aicc <- c()
I2 <- c()

for (i in 1:length(mods)) {
  # pseudo R2 of proportional reduction of variance explained when fitting reduced model relative to model with random effects only
  pseudo.r2[i] <- (sum(rma.random$sigma2) - sum(mods[[i]]$sigma2)) / sum(rma.random$sigma2)
  
  # marginal and conditional R2 for linear mixed-mods from Nakagawa and Scheilzeth 2013
  # variance from fixed effects
  fe.total <- c()
  for (j in 2:ncol(model.matrix(mods[[i]]))) {
    fe <- mods[[i]]$b[j] * model.matrix(mods[[i]])[, j] 
    fe.total <- c(fe.total, fe)
  }
  v.fix <- var(fe.total)
  # variance from random effects
  v.rand <- sum(mods[[i]]$sigma2)
  # variance from residuals
  v.resid <- var(residuals(mods[[i]]))
  # marginal R2 - total variance explained from fixed effects
  marg.r2[i] <- v.fix / (v.fix + v.rand + v.resid)
  # conditional R2 - total variance explained from fixed and random effects
  cond.r2[i] <- (v.fix + v.rand) / (v.fix + v.rand + v.resid)
  
  # from Wolfgang Viechtbauer (http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate)
  # I2 statistic - amount of heterogeneity relative to the total amount of variance in observed effects (how much heterogeneity contributes to total variance)
  W <- diag(1/ROM.Biomass$vi)
  X <- model.matrix(mods[[i]])
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2[i] <- 100 * sum(mods[[i]]$sigma2) / (sum(mods[[i]]$sigma2) + (mods[[i]]$k - mods[[i]]$p) / sum(diag(P)))
  
  aicc[i] <- fitstats(mods[[i]])[5]
}
# 
cbind(pseudo.r2,marg.r2,cond.r2,aicc,I2)

vif(ROM.ai.rma, table=TRUE)
vif(ROM.prod.int1, table = TRUE) #interaction not significant
vif(ROM.prod.int2, table = TRUE)
vif(ROM.prod.int3, table=TRUE)

ROM.ai.rma

ROM.prod.int2
ROM.prod.int3



#=========================================graphs=========

preds<-predict(ROM.ai.rma, addx = TRUE)

preds<-as.data.frame(preds)



#model with time as variable
ROM.biomass.time<-rma.mv(yi, vi, mods = ~time, 
                         random = rand.var,data = ROM.Biomass,slab = paste(author, year), method="ML")

preds2<-predict(ROM.biomass.time, addx = TRUE)
preds2<-as.data.frame(preds2)


bubble.time2=ggplot(preds2, aes(x=X.time,y=pred))+
  geom_ribbon(aes(ymax=ci.lb, ymin=ci.ub), fill="gray83", alpha=.5) +

  stat_smooth(data=preds2,method="glm",fullrange=F,se=F, size=1, color="black")+ #, linetype=Biome
  geom_hline(yintercept=0, linetype="dashed", size=.5)+
  geom_point(data=ROM.Biomass, aes(x=time, y = yi, size=1/vi, alpha=.5))+
  ylab("Log Response Ratio")+
  xlab("")+
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=10),
        axis.title=element_text(size=10), legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylim(-3,6)

bubble.time2

 
# tiff("prod.time.tiff", width = 16, height= 12, units ='cm', res=600)
# bubble.time
# dev.off()

#plots with interactions
 
#subset by burn
prod.int2a.burn<-rma.mv(yi, vi, mods = ~ time, random = rand.var,data = ROM.Biomass,slab = paste(author, year), subset = Burn..Y.N.=="Y", method = "ML")
preds.b<-predict(prod.int2a.burn, addx = TRUE)
preds.b<-as.data.frame(preds.b)
dim(preds.b)
 
 preds.b$Burn<-rep("Y", times = 25)
 
 prod.int2a.burn.n<-rma.mv(yi, vi, mods = ~ time, random = rand.var,data = ROM.Biomass,slab = paste(author, year), subset = Burn..Y.N.=="N", method = "ML")
 preds.b.n<-predict(prod.int2a.burn.n, addx = TRUE)
 preds.b.n<-as.data.frame(preds.b.n)
 dim(preds.b.n)
 
 preds.b.n$Burn<-rep("N", times = 244)
 
 preds.burn=rbind(preds.b,preds.b.n)
 names(preds.burn)

 bubble.time.Burn=ggplot(preds.burn, aes(x=X.time,y=pred))+
   stat_smooth(data=preds.burn,aes(color= Burn),method="lm",formula=y~x, level =.95,fullrange=F,se=F,size=1)+ #, linetype=Burn 
   geom_ribbon( aes(ymin = preds.burn$ci.lb, ymax = preds.burn$ci.ub, fill = Burn), alpha = .15) +
   #scale_linetype_manual(values=biome.line)+
   geom_hline(yintercept=0, linetype="dashed", size=.5)+
  scale_color_manual(values=c("black","red"))+
   scale_fill_manual(values=c("black","red"))+
   geom_point(data=ROM.Biomass, aes(x=time, y = yi, size=1/vi, alpha=.5,color=Burn..Y.N., group=Burn..Y.N., fill=Burn..Y.N., shape=Burn..Y.N.))+
   ylab("")+
   xlab("Years since restoration")+
   # scale_color_manual(values=Burn.colors)+
   # scale_fill_manual(values=Burn.colors)+
   # 
   theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=10),
         axis.title=element_text(size=15), legend.position = "none")+
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
   ylim(-3,6)
 bubble.time.Burn
 
 # tiff("prod.time.burn.tiff", width = 16, height= 12, units ='cm', res=600)
 # bubble.time.Burn
 # dev.off()
 
 
 # #subset by seeded
 # 
 # prod.int3a.Seeded<-rma.mv(yi, vi, mods = ~ time, random = rand.var,data = ROM.Biomass,slab = paste(author, year), subset = Seeded..Y.N.=="Y", method = "ML")
 # preds.s<-predict(prod.int3a.Seeded, addx = TRUE)
 # preds.s<-as.data.frame(preds.s)
 # dim(preds.s)
 # 
 # preds.s$Seeded<-rep("Y", times = 126)
 # 
 # prod.int3a.Seeded.n<-rma.mv(yi, vi, mods = ~ time, random = rand.var,data = ROM.Biomass,slab = paste(author, year), subset = Seeded..Y.N.=="N", method = "ML")
 # preds.s.n<-predict(prod.int3a.Seeded.n, addx = TRUE)
 # preds.s.n<-as.data.frame(preds.s.n)
 # dim(preds.s.n)
 # 
 # preds.s.n$Seeded<-rep("N", times = 143)
 # 
 # preds.Seeded=rbind(preds.s,preds.s.n)
 # names(preds.Seeded)
 # 
 # bubble.time.Seeded=ggplot(preds.Seeded, aes(x=X.time,y=pred))+
 #   stat_smooth(data=preds.Seeded,aes(color= Seeded, fill=Seeded),method="glm",formula=y~x,fullrange=F,se=F,size=1     )+ #, linetype=Seeded 
 #   geom_ribbon( aes(ymin = preds.Seeded$ci.lb, ymax = preds.Seeded$ci.ub, fill = Seeded), alpha = .15) +
 #   
 #   geom_hline(yintercept=0, linetype="dashed", size=.5)+
 #   scale_color_manual(values=c("darkblue","forestgreen"))+
 #   scale_fill_manual(values=c("darkblue","forestgreen"))+
 #   geom_point(data=ROM.Biomass, aes(x=time, y = yi, size=1/vi, alpha=.5,color=Seeded..Y.N., group=Seeded..Y.N., fill=Seeded..Y.N., shape=Seeded..Y.N.))+
 #   ylab("")+
 #   xlab("")+
 #   # scale_color_manual(values=Seeded.colors)+
 #   # scale_fill_manual(values=Seeded.colors)+
 #   # 
 #   theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=10),
 #         axis.title=element_text(size=10), legend.position = "none")+
 #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
 #         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
 #   ylim(-3,6)
 # bubble.time.Seeded
 # 
 # # tiff("prod.time.seed.tiff", width = 16, height= 12, units ='cm', res=600)
 # # bubble.time.Seeded
 # # dev.off()
 # 
 # library(cowplot)
 # library(ggpubr)
 # 
 # 
 # biomass.plot<-plot_grid(bubble.time2, bubble.time.Burn,bubble.time.Seeded,
 #                       labels = c("a", "b", "c"),
 #                       ncol = 3, nrow = 1)
 # 
 # 
 # # 
 # # tiff("biomass.tiff", width = 25, height= 8, units ='cm', res=600)
 # # biomass.plot
 # # dev.off()
 # # 