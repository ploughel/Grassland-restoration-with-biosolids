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

Biomass <- read_excel("~/Documents/GitHub/Grassland restoration with biosolids/Grassland-restoration-with-biosolids/data_sheets.xlsx", 
                      sheet = "Biomass")


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


# rma.glmulti<- function(formula, data, ...)
#   rma(formula, vi, data=data, method ="ML",...)
# 
# ROM.ai.size <- glmulti(yi ~Biosolid.level..Mg.ha.1.+time+Temp+Precip+Mixture..yes.no.+Burn..Y.N.+
#                         Seeded..Y.N.+Multiple.application..Y.N.+Severe.Disturbance+ai, data=ROM.Biomass, method="d",
#                       level=1, crit="aicc", fitfunction = rma.glmulti)
# 
# ROM.ai.mod <- glmulti(yi ~ Biosolid.level..Mg.ha.1.+time+Temp+Precip+Mixture..yes.no.+Burn..Y.N.+
#                        Seeded..Y.N.+Multiple.application..Y.N.+Severe.Disturbance+ai, data=ROM.Biomass, method="h",
#                      level=1, crit="aicc", confsetsize=ROM.ai.size, fitfunction = rma.glmulti)
# 
# plot(ROM.ai.mod, type="s")




ROM.ai.rma<-rma.mv(yi, vi, mods = ~Burn..Y.N.+ time+ Temp+ Severe.Disturbance+ Mixture..yes.no.+ Seeded..Y.N., 
                       random = rand.var,data = ROM.Biomass,slab = paste(author, year), method="ML")


#models with interactions

ROM.prod.int1<-rma.mv(yi, vi, mods = ~time*Mixture..yes.no.,
                     random = rand.var,data = ROM.Biomass,slab = paste(author, year), method="ML")

ROM.prod.int2<-rma.mv(yi, vi, mods = ~time*Burn..Y.N.,
                      random = rand.var,data = ROM.Biomass,slab = paste(author, year), method="ML")

ROM.prod.int3<-rma.mv(yi, vi, mods = ~time*Seeded..Y.N.,
                      random = rand.var,data = ROM.Biomass,slab = paste(author, year), method="ML")

ROM.prod.int4<-rma.mv(yi, vi, mods = ~time*Severe.Disturbance,
                      random = rand.var,data = ROM.Biomass,slab = paste(author, year), method="ML")

ROM.prod.int5<-rma.mv(yi, vi, mods = ~Temp*Mixture..yes.no.,
                      random = rand.var,data = ROM.Biomass,slab = paste(author, year), method="ML")

ROM.prod.int6<-rma.mv(yi, vi, mods = ~Temp*Burn..Y.N.,
                       random = rand.var,data = ROM.Biomass,slab = paste(author, year), method="ML")

ROM.prod.int7<-rma.mv(yi, vi, mods = ~Temp*Seeded..Y.N.,
                       random = rand.var,data = ROM.Biomass,slab = paste(author, year), method="ML")

ROM.prod.int8<-rma.mv(yi, vi, mods = ~Temp*Severe.Disturbance,
                       random = rand.var,data = ROM.Biomass,slab = paste(author, year), method="ML")



mods = list(ROM.ai.rma, ROM.prod.int1,ROM.prod.int2,ROM.prod.int3,ROM.prod.int5,ROM.prod.int6,ROM.prod.int7,
            ROM.prod.int8)


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


# 
# # write.csv(as.data.frame(coef(summary(ROM.ai.rma))), file="prod.mod.csv")
# df.prod<-read.csv("prod.mod.csv")
# df.prod
# 
# df.ai<-df.prod[2:8,]
# 
# LRR.prod = ggplot(data=df.ai,
#                         aes(x = X, y = estimate, ymin =ci.lb , ymax = ci.ub ))+
#   geom_point(aes(color=X))+
#   geom_hline(yintercept =0, linetype=2)+
#   scale_color_brewer(type = 'div', palette = 4)+
#   xlab('ai')+ ylab("Log Response Ratio")+
#   geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub,col=X),width=0.3,cex=1)+ 
#   #xlim(breaks =rev(mod.lev))+
#   #facet_wrap(~X,scale="free") +
#   theme(plot.title=element_text(size=16,face="bold"),
#         axis.text.x=element_text(face="bold"),
#         axis.title=element_text(size=12,face="bold"),
#         legend.position = c(.8,.7), axis.text.y=element_blank())+
#   theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=20),
#         axis.title=element_text(size=20),legend.title=element_text(size=15), 
#         legend.text=element_text(size=13))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))  +
#   
#   
#   coord_flip()
# LRR.prod

# tiff("prod.tiff", width = 20, height= 12, units ='cm', res=600)
# LRR.prod
# dev.off()



preds<-predict(ROM.ai.rma, addx = TRUE)

names(preds)

preds<-as.data.frame(preds)
names(preds)




#forest plot

forest.prod<-viz_forest(x =preds[c("pred", "se")], variant="classic",summary_label ="Summary Effect", xlab = "Log Response Ratio")
# tiff("prod.forest.tiff", width = 16, height= 25, units ='cm', res=600)
# forest.prod
# dev.off()

bubble.time.auto=ggplot(preds, aes(x=X.time,y=pred))+
  stat_smooth(data=preds,method="auto",formula=y~x,fullrange=T,se=T,size=1, color="black", alpha=0.9)+ #, linetype=Biome 

  geom_hline(yintercept=0, linetype="dashed", size=.5)+
  scale_color_manual(values=c("darkblue","forestgreen"))+
  scale_fill_manual(values=c("darkblue","forestgreen"))+
  geom_hline(yintercept=0, linetype="dashed", size=.5)+
  # stat_smooth(data=preds.biome,aes(x=level.exp,y=ci.lb, color=Biome),method="glm",formula=y~x,fullrange=T,
  #             se=FALSE,size=1,linetype=3)+
  # stat_smooth(data=preds.biome,aes(x=level.exp,y=ci.ub),method="glm",formula=y~x,fullrange=T,
  #             se=FALSE,size=1,linetype=3)+
  # scale_linetype_manual(values=biome.line)+
  geom_point(data=ROM.Biomass, aes(x=time, y = yi, size=1/vi))+
  ylab("Log Response Ratio")+
  xlab("Years Since Restoration")+
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=15),
        axis.title=element_text(size=15), legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylim(-2,6)


bubble.time.auto

# tiff("prod.time.auto.tiff", width = 16, height= 12, units ='cm', res=600)
# bubble.time.auto
# dev.off()
dim(preds)
bubble.time=ggplot(preds, aes(x=X.time,y=pred))+
  stat_smooth(data=preds,method="glm",formula=y~x,fullrange=F,se=T,aes(ymin=ci.lb, ymax=ci.ub), size=1, color="black")+ #, linetype=Biome
  geom_hline(yintercept=0, linetype="dashed", size=.5)+

  # stat_smooth(data=preds,aes(x=X.time,y=ci.lb),method="glm",formula=y~x,fullrange=T,
  #             se=FALSE,size=1,linetype=3)+
  # stat_smooth(data=preds,aes(x=X.time,y=ci.ub),method="glm",formula=y~x,fullrange=T,
  #             se=FALSE,size=1,linetype=3)+
  geom_point(data=ROM.Biomass, aes(x=time, y = yi, size=1/vi))+
  ylab("Log Response Ratio")+
  xlab("")+
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=10),
        axis.title=element_text(size=15), legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylim(-3,6)

bubble.time


# tiff("prod.time.tiff", width = 16, height= 12, units ='cm', res=600)
# bubble.time
# dev.off()

bubble.temp.auto=ggplot(preds, aes(x=X.Temp,y=pred))+
  stat_smooth(data=preds,method="auto",formula=y~x,fullrange=T,se=TRUE,size=1, color="black"     )+ #, linetype=Biome 
  geom_hline(yintercept=0, linetype="dashed", size=.5)+
  # stat_smooth(data=preds.biome,aes(x=level.exp,y=ci.lb, color=Biome),method="glm",formula=y~x,fullrange=T,
  #             se=FALSE,size=1,linetype=3)+
  # stat_smooth(data=preds.biome,aes(x=level.exp,y=ci.ub),method="glm",formula=y~x,fullrange=T,
  #             se=FALSE,size=1,linetype=3)+
  # scale_linetype_manual(values=biome.line)+
  geom_point(data=ROM.Biomass, aes(x=Temp, y = yi, size=1/vi))+
  ylab("Log Response Ratio")+
  xlab(expression('MAT ('*~degree*C*')'))+
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=10),
        axis.title=element_text(size=10), legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylim(-2,6)+
  xlim(0,20)


bubble.temp.auto
# 
# tiff("prod.temp.auto.tiff", width = 16, height= 12, units ='cm', res=600)
# bubble.temp.auto
# dev.off()

bubble.temp=ggplot(preds, aes(x=X.Temp,y=pred))+
  stat_smooth(data=preds,method="lm",formula=y~x,fullrange=F,se=T,size=1, color="black" , level=.95)+ #, linetype=Biome 
  geom_hline(yintercept=0, linetype="dashed", size=.5)+
  geom_point(data=ROM.Biomass, aes(x=Temp, y = yi, size=1/vi))+
  ylab("")+
  xlab(expression('MAT ('*~degree*C*')'))+
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=10),
        axis.title=element_text(size=10), legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylim(-2,6)+
  xlim(0,20)


bubble.temp
 # 
# tiff("prod.temp.auto.tiff", width = 16, height= 12, units ='cm', res=600)
# bubble.temp.auto
# dev.off()
bubble.temp.spline=ggplot(preds, aes(x=X.Temp,y=pred))+
  stat_smooth(data=preds,method="lm",formula=y~splines::bs(x, 2),fullrange=F,se=T,size=1, color="black" , level=.95)+ #, linetype=Biome 
  geom_hline(yintercept=0, linetype="dashed", size=.5)+
  geom_point(data=ROM.Biomass, aes(x=Temp, y = yi, size=1/vi))+
  ylab("")+
  xlab(expression('MAT ('*~degree*C*')'))+
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=10),
        axis.title=element_text(size=10), legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylim(-2,6)+
  xlim(0,20)


bubble.temp.spline
 #plots with interactions
 summary(ROM.prod.int2) # interaction between burn and time

 
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
   
   geom_point(data=ROM.Biomass, aes(x=time, y = yi, size=1/vi, color=Burn..Y.N., group=Burn..Y.N., fill=Burn..Y.N., shape=Burn..Y.N.))+
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
 
 
 #subset by seeded

 prod.int3a.Seeded<-rma.mv(yi, vi, mods = ~ time, random = rand.var,data = ROM.Biomass,slab = paste(author, year), subset = Seeded..Y.N.=="Y", method = "ML")
 preds.s<-predict(prod.int3a.Seeded, addx = TRUE)
 preds.s<-as.data.frame(preds.s)
 dim(preds.s)
 
 preds.s$Seeded<-rep("Y", times = 126)
 
 prod.int3a.Seeded.n<-rma.mv(yi, vi, mods = ~ time, random = rand.var,data = ROM.Biomass,slab = paste(author, year), subset = Seeded..Y.N.=="N", method = "ML")
 preds.s.n<-predict(prod.int3a.Seeded.n, addx = TRUE)
 preds.s.n<-as.data.frame(preds.s.n)
 dim(preds.s.n)
 
 preds.s.n$Seeded<-rep("N", times = 143)
 
 preds.Seeded=rbind(preds.s,preds.s.n)
 names(preds.Seeded)
 
 bubble.time.Seeded=ggplot(preds.Seeded, aes(x=X.time,y=pred))+
   stat_smooth(data=preds.Seeded,aes(color= Seeded, fill=Seeded),method="glm",formula=y~x,fullrange=F,se=F,size=1     )+ #, linetype=Seeded 
   geom_ribbon( aes(ymin = preds.Seeded$ci.lb, ymax = preds.Seeded$ci.ub, fill = Seeded), alpha = .15) +
   
   geom_hline(yintercept=0, linetype="dashed", size=.5)+
   scale_color_manual(values=c("darkblue","forestgreen"))+
   scale_fill_manual(values=c("darkblue","forestgreen"))+
   geom_point(data=ROM.Biomass, aes(x=time, y = yi, size=1/vi, color=Seeded..Y.N., group=Seeded..Y.N., fill=Seeded..Y.N., shape=Seeded..Y.N.))+
   ylab("")+
   xlab("")+
   # scale_color_manual(values=Seeded.colors)+
   # scale_fill_manual(values=Seeded.colors)+
   # 
   theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=10),
         axis.title=element_text(size=10), legend.position = "none")+
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
   ylim(-3,6)
 bubble.time.Seeded
 
 # tiff("prod.time.seed.tiff", width = 16, height= 12, units ='cm', res=600)
 # bubble.time.Seeded
 # dev.off()

 library(cowplot)
 library(ggpubr)
 
 
 biomass.plot<-plot_grid(bubble.time, bubble.time.Burn,bubble.time.Seeded,
                       labels = c("a", "b", "c"),
                       ncol = 3, nrow = 1)
 
 
 
 # tiff("biomass.tiff", width = 25, height= 8, units ='cm', res=600)
 # biomass.plot
 # dev.off()
 