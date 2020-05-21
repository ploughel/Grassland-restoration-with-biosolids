#Cover meta-analysis
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


Cover <- read_excel("~/Documents/GitHub/Grassland restoration with biosolids/Grassland-restoration-with-biosolids/data_sheetscorr.xlsx", 
                    sheet = "Cover")%>%
        janitor::clean_names()

Cover$mc_0 <- ifelse(Cover$mc ==0, Cover$mc+1/2 , Cover$mc)
Cover$me_0 <- ifelse(Cover$mc ==0, Cover$me+1/2 , Cover$me)


ROM.Cover <- escalc(n2i = nc, n1i = ne, m2i = mc_0, m1i = me_0, 
                    sd2i = scimp, sd1i = seimp, data = Cover, measure = "ROM", 
                    append = TRUE )


HG.Cover <- escalc(n2i = nc, n1i = ne, m2i = mc_0, m1i = me_0, 
                   sd2i = scimp, sd1i = seimp, data = Cover, measure = "SMD", 
                   append = TRUE )



Hedges = HG.Cover$yi
# #
h <- hist(Hedges, breaks = 10, density = 10,
          xlab = "Hedges' (g)")
xfit <- seq(min(Hedges,na.rm=TRUE), max(Hedges,na.rm=TRUE), length = 40)
yfit <- dnorm(xfit, mean = mean(Hedges,na.rm=TRUE), sd = sd(Hedges,na.rm=TRUE))
yfit <- yfit * diff(h$mids[1:2]) * length(Hedges)
# #
lines(xfit, yfit, col = "black", lwd = 2)
# #
LRR = ROM.Cover$yi
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

rma.random <- rma.mv(yi = yi, V = vi, random = rand.var, data = ROM.Cover, method = "ML")
summary(rma.random)


#-----------------------# model selection 
# rma.glmulti<- function(formula, data, ...)
#   rma(formula, vi, data=data, method ="ML",...)
# 
# ROM.cover.size <- glmulti(yi ~biosolid_level_mg_ha_1+yeartrans+temp+precip+mixture+burn+
#                          seeded+multiple_application+s_dist+ai, data=ROM.Cover, method="d",
#                       level=1, crit="aicc", fitfunction = rma.glmulti)
# 
# 
# ROM.cover.bf <- glmulti(yi ~ biosolid_level_mg_ha_1+yeartrans+temp+precip+mixture+burn+
#                         seeded+multiple_application+s_dist+ai, data=ROM.Cover, method="h",
#                      level=1, crit="aicc", confsetsize=ROM.cover.size, fitfunction = rma.glmulti)
# 
# 
# plot(ROM.cover.bf, type="s")

#best fit model
ROM.cover.bf<-rma.mv(yi, vi, mods = ~temp+s_dist+burn+multiple_application+seeded+precip, 
                       random = rand.var,data = ROM.Cover,slab = paste(author, year), method="ML")


#model with interactions
ROM.cover.int.1<-rma.mv(yi, vi, mods = ~temp*seeded, 
                      random = rand.var,data = ROM.Cover,slab = paste(author, year), method="ML")

ROM.cover.int.2<-rma.mv(yi, vi, mods = ~temp*burn, 
                       random = rand.var,data = ROM.Cover,slab = paste(author, year), method="ML")

ROM.cover.int.3<-rma.mv(yi, vi, mods = ~temp*multiple_application, 
                       random = rand.var,data = ROM.Cover,slab = paste(author, year), method="ML")

ROM.cover.int.4<-rma.mv(yi, vi, mods = ~temp*s_dist, 
                        random = rand.var,data = ROM.Cover,slab = paste(author, year), method="ML")

ROM.cover.int.5<-rma.mv(yi, vi, mods = ~precip*seeded, 
                       random = rand.var,data = ROM.Cover,slab = paste(author, year), method="ML")

ROM.cover.int.6<-rma.mv(yi, vi, mods = ~precip*burn, 
                       random = rand.var,data = ROM.Cover,slab = paste(author, year), method="ML")

ROM.cover.int.7<-rma.mv(yi, vi, mods = ~precip*multiple_application, 
                       random = rand.var,data = ROM.Cover,slab = paste(author, year), method="ML")
ROM.cover.int.8<-rma.mv(yi, vi, mods = ~precip*s_dist, 
                        random = rand.var,data = ROM.Cover,slab = paste(author, year), method="ML")


mods = list(ROM.cover.bf,ROM.cover.int.1, ROM.cover.int.2,
            ROM.cover.int.3, ROM.cover.int.4,ROM.cover.int.5, ROM.cover.int.6,
            ROM.cover.int.7, ROM.cover.int.8)



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
  W <- diag(1/ROM.Cover$vi)
  X <- model.matrix(mods[[i]])
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2[i] <- 100 * sum(mods[[i]]$sigma2) / (sum(mods[[i]]$sigma2) + (mods[[i]]$k - mods[[i]]$p) / sum(diag(P)))
  
  aicc[i] <- fitstats(mods[[i]])[5]
}
# 

cover.mods=cbind(pseudo.r2,marg.r2,cond.r2,aicc,I2)
cover.mods

vif(ROM.cover.bf, table=TRUE)

vif(ROM.cover.int.1, table=TRUE)
vif(ROM.cover.int.2, table=TRUE)
vif(ROM.cover.int.4, table=TRUE)
vif(ROM.cover.int.5, table=TRUE)
vif(ROM.cover.int.6, table=TRUE) #Vif high
vif(ROM.cover.int.8, table=TRUE)

ROM.cover.bf
ROM.cover.int.1 #temp and seeded significant, but lower aic
ROM.cover.int.2
ROM.cover.int.4
ROM.cover.int.5
ROM.cover.int.8

#burn (main model) and seeded and temperature in interaction model

#=========================================graphs=========
cover.temp.mod<-rma.mv(yi, vi, mods = ~temp, 
                                     random = rand.var,data = ROM.Cover,slab = paste(author, year), method="ML")

preds<-predict(ROM.cover.bf, addx = TRUE)
preds<-as.data.frame(preds)

preds2<-predict(cover.temp.mod, addx = TRUE)
preds2<-as.data.frame(preds2)



bubble.temp=ggplot(preds2, aes(x=X.temp,y=pred))+
  geom_ribbon(aes(ymax=ci.lb, ymin=ci.ub), fill="gray83", alpha=.5) +
  stat_smooth(method="glm",fullrange=F,se=F,size=1, color="black")+ 
  geom_hline(yintercept=0, linetype="dashed", size=.5)+
  
  geom_point(data=ROM.Cover, aes(x=temp, y = yi, size=1/vi, alpha=.5))+
  ylab("")+
  xlab(expression('MAT ('*degree*C*')'))+
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=10),
        axis.title=element_text(size=10), legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylim(-2,5)

bubble.temp

#  tiff("cover.temp.tiff", width = 16, height= 12, units ='cm', res=600)
# bubble.temp
#  dev.off()

#burned
mod1.burn.n<-rma.mv(yi, vi, mods = ~ temp+s_dist+seeded+precip, random = rand.var,data = ROM.Cover,slab = paste(author, year), subset = burn=="N", method = "ML")
mod1.burn<-rma.mv(yi, vi, mods = ~ temp+s_dist+seeded+precip, random = rand.var,data = ROM.Cover,slab = paste(author, year), subset = burn=="Y", method = "ML")



preds.burn<-predict(mod1.burn, addx = TRUE)
preds.b<-as.data.frame(preds.burn)
dim(preds.b)
preds.b$burned<-rep("Y", times = 50)

mean(preds.burn$pred)
preds.nb<-predict(mod1.burn.n, addx = TRUE)

preds.nb<-as.data.frame(preds.nb)
preds.nb$burned<-rep("N", times = 164)
names(preds.nb)
names(preds.b)
preds.burned=rbind(preds.b,preds.nb)



# burn.n=as.data.frame(preds.nb)%>%
#   summarise(mean = mean(pred, na.rm = TRUE),
#             sd = sd(pred, na.rm = TRUE),
#             n = n()) %>%
#   mutate(se = sd / sqrt(n),
#          lower = mean - qt(1 - (0.05 / 2), n - 1) * se,
#          upper = mean + qt(1 - (0.05 / 2), n - 1) * se)
# 
# burn.y=as.data.frame(preds.b)%>%
#   summarise(mean = mean(pred, na.rm = TRUE),
#             sd = sd(pred, na.rm = TRUE),
#             n = n()) %>%
#   mutate(se = sd / sqrt(n),
#          lower = mean - qt(1 - (0.05 / 2), n - 1) * se,
#          upper = mean + qt(1 - (0.05 / 2), n - 1) * se)
# 

bubble.burn=ggplot(preds.burned, aes(x=burned,y=pred, color=burned))+
  geom_boxplot()+
  geom_jitter(data=ROM.Cover, aes(x=burn, y = yi, size=1/vi, color=burn,alpha=.5))+
  geom_hline(yintercept=0, linetype="dashed", size=.5)+
  scale_color_manual(values=c("black", "red"))+
  #geom_point(data=ROM.Cover, aes(x=burn, y = yi, size=1/vi))+
  ylab("Log Response Ratio")+
  xlab("Burn")+
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=10),
        axis.title=element_text(size=10), legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylim(-2,5)

bubble.burn

#   tiff("cover.burn.tiff", width = 16, height= 12, units ='cm', res=600)
# bubble.burn
#dev.off()

mod1.burn.n<-rma.mv(yi, vi,  random = rand.var,data = ROM.Cover,slab = paste(author, year), subset = burn=="N", method = "ML")
mod1.burn<-rma.mv(yi, vi,  random = rand.var,data = ROM.Cover,slab = paste(author, year), subset = burn=="Y", method = "ML")

#seeded

mod1.seeded.n<-rma.mv(yi, vi, mods = ~ temp, random = rand.var,
                      data = ROM.Cover,slab = paste(author, year), subset = seeded=="N", method = "ML")

mod1.seeded<-rma.mv(yi, vi, mods = ~ temp, random = rand.var,
                    data = ROM.Cover,slab = paste(author, year), subset = seeded=="Y", method = "ML")



preds.seeded<-predict(mod1.seeded, addx = TRUE)
preds.b<-as.data.frame(preds.seeded)
dim(preds.b)
preds.b$seeded<-rep("Y", times = 93)

preds.nb<-predict(mod1.seeded.n, addx = TRUE)

preds.nb<-as.data.frame(preds.nb)
preds.nb$seeded<-rep("N", times = 121)
names(preds.nb)
names(preds.b)
preds.seeded=rbind(preds.b,preds.nb)

bubble.seeded=ggplot(preds.seeded, aes(x=seeded,y=pred, color=seeded))+
  geom_boxplot()+
  geom_jitter(data=ROM.Cover, aes(x=seeded, y = yi, size=1/vi, color=seeded,alpha=.5))+

  geom_hline(yintercept=0, linetype="dashed", size=.5)+
  scale_color_manual(values=c("darkblue","forestgreen"))+
  #geom_point(data=ROM.Cover, aes(x=seeded, y = yi, size=1/vi))+
  ylab("")+
  xlab("Seeded")+
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=10),
        axis.title=element_text(size=10), legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylim(-2,5)

bubble.seeded

#   tiff("cover.seeded.tiff", width = 16, height= 12, units ='cm', res=600)
# bubble.seeded
#  dev.off()

mod.seed<-rma.mv(yi, vi,  random = rand.var,data = ROM.Cover,slab = paste(author, year), subset = seeded=="Y", method = "ML")
mod.seed.n<-rma.mv(yi, vi,  random = rand.var,data = ROM.Cover,slab = paste(author, year), subset = seeded=="N", method = "ML")




library(cowplot)
library(ggpubr)


cover.plot<-plot_grid(bubble.burn,bubble.temp, bubble.seeded, 
                      labels = c("a", "b", "c"),
                      ncol = 3, nrow = 1)

cover.plot

#  tiff("cover2.tiff", width = 25, height= 8, units ='cm', res=600)
# cover.plot
#  dev.off()
#  
