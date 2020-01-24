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


BSma.kgm<-read.csv("biomass.csv", header=T) 


#obtain value to multiply means by when sd missing
# BS.se<-BSma.kgm%>%
#   filter(Se!="NA")%>%
#   summarise(value=sum(Se)/sum(Me))
#
#
# BS.sc<-BSma.kgm%>%
#   filter(Se!="NA")%>%
#   summarise(value=sum(Sc)/sum(Mc))



# #justification for use of log response ratio over Hedges'
# HG.kgm <- escalc(n2i = Nc, n1i = Ne, m2i = Mc, m1i = Me,
#                   sd2i = Sc, sd1i = Se, data = BSma.kgm, measure = "SMD",
#                   append = TRUE)

BSma.kgm$Mc2<-ifelse(BSma.kgm$Mc==0,1/2,BSma.kgm$Mc)


head(BSma.kgm)
ROM.kgm <- escalc(n2i = Nc, n1i = Ne, m2i = Mc2, m1i = Me, 
               sd2i = Sc, sd1i = Se, data = BSma.kgm, measure = "ROM", 
               append = TRUE )

# write.csv(ROM.kgm, file="ROM.prod.csv")

# 
# # # Check normality of Hedges' vs. log response ratio
#  par(mfrow=c(1,2))
# # 
# Hedges = HG.kgm$yi
# # #
# h <- hist(Hedges, breaks = 10, density = 10,
#  xlab = "Hedges' (g)")
# xfit <- seq(min(Hedges,na.rm=TRUE), max(Hedges,na.rm=TRUE), length = 40)
# yfit <- dnorm(xfit, mean = mean(Hedges,na.rm=TRUE), sd = sd(Hedges,na.rm=TRUE))
# yfit <- yfit * diff(h$mids[1:2]) * length(Hedges)
# # #
# lines(xfit, yfit, col = "black", lwd = 2)
# # #
# LRR = ROM.kgm$yi
# # #
# r <- hist(LRR, breaks = 10, density = 10, xlab = "Log Response Ratio")
# xfit.r <- seq(min(LRR), max(LRR), length = 40)
# yfit.r <- dnorm(xfit.r, mean = mean(LRR), sd = sd(LRR))
# yfit.r <- yfit.r * diff(r$mids[1:2]) * length(LRR)
# # #
# lines(xfit.r, yfit.r, col = "black", lwd = 2)
# # 

#random effects model

rand.var=list(~ 1|Experiment/Paper.ID./Case)
 
rma.random <- rma.mv(yi = yi, V = vi, random = rand.var, data = ROM.kgm, method = "ML")
summary(rma.random)
ran.df<-as.data.frame(coef(summary(rma.random)))

#transform continuous variables
ROM.kgm$level.log<-log(ROM.kgm$Biosolid.level..Mg.ha.1.)
ROM.kgm$year.log<-log(ROM.kgm$time)
ROM.kgm$temp.sq<-sqrt(ROM.kgm$Temp) 
ROM.kgm$precip.log<-log(ROM.kgm$Precip) 

#log is better than sqrt
names(ROM.kgm)

#models with single covariates
rma.lev<-rma.mv(yi,vi,mods = ~ -1+level.log,random =rand.var, data = ROM.kgm, slab = paste(author, year),method ="ML")
lev.df<-as.data.frame(coef(summary(rma.lev)))
rma.yr<-rma.mv(yi,vi,mods = ~ -1+year.log, random =rand.var,data = ROM.kgm,slab = paste(author, year),method ="ML")
year.df<-as.data.frame(coef(summary(rma.yr)))
rma.temp<-rma.mv(yi,vi,mods = ~ -1+ temp.sq, random =rand.var,data = ROM.kgm,slab = paste(author, year),method ="ML")
temp.df<-as.data.frame(coef(summary(rma.temp)))
rma.precip<-rma.mv(yi,vi,mods = ~ -1+ precip.log, random =rand.var,data = ROM.kgm,slab = paste(author, year),method ="ML")
precip.df<-as.data.frame(coef(summary(rma.precip)))

rma.mix<-rma.mv(yi,vi,mods = ~ -1+Mixture..yes.no., random =rand.var,data = ROM.kgm,slab = paste(author, year),method ="ML")
mix.df<-as.data.frame(coef(summary(rma.mix)))
rma.dis<-rma.mv(yi,vi,mods = ~ -1+Severe.Disturbance, random =rand.var,data = ROM.kgm,slab = paste(author, year),method ="ML")
dis.df<-as.data.frame(coef(summary(rma.dis)))
rma.mult<-rma.mv(yi,vi,mods = ~ -1+Multiple.application..Y.N., random =rand.var,data = ROM.kgm,slab = paste(author, year),method ="ML")
mult.df<-as.data.frame(coef(summary(rma.mult)))
rma.seed<-rma.mv(yi,vi,mods = ~ -1+Seeded..Y.N., random =rand.var,data = ROM.kgm,slab = paste(author, year),method ="ML")
seed.df<-as.data.frame(coef(summary(rma.seed)))
rma.Biome<-rma.mv(yi,vi,mods = ~ -1+Biome, random =rand.var,data = ROM.kgm,slab = paste(author, year),method ="ML")
Biome.df<-as.data.frame(coef(summary(rma.Biome)))

#simple.models<-rbind(ran.df,lev.df,year.df,temp.df,precip.df,mix.df,dis.df,mult.df,seed.df, Biome.df)
#AIC(rma.random,rma.lev,rma.yr, rma.temp,rma.precip,rma.mix,rma.dis,rma.mult, rma.seed, rma.Biome)
#write.csv(simple.models, "simple.models.csv")

# Explore variable correlations ---------------------------------------

# pairs(~ level.log + year.log + Severe.Disturbance +  Multiple.application..Y.N. +Mixture..yes.no.  + Biome+ Seeded..Y.N.+Temp+Precip, data = ROM.kgm)
# library("polycor")
# 
# 
# var<-ROM.kgm%>%
#   dplyr::select(level.log , year.log , Severe.Disturbance ,  Multiple.application..Y.N. ,Mixture..yes.no.  , Biome, Seeded..Y.N.,temp.sq,precip.log)
# 
# correlations <- hetcor(data = var, std.err = FALSE)

#-----------------------# model selection 
library(xlsx)
library(rJava)
library("glmulti")

#model 1 uses biome, not temp/precip
#  model1.size <- glmulti(yi ~ level.log + year.log +   Severe.Disturbance+Multiple.application..Y.N. + Biome+Seeded..Y.N. +Mixture..yes.no., data=ROM.kgm, method="d",
#                     level=1, crit="aicc")
# model1 <- glmulti(yi ~ level.log + year.log +   Severe.Disturbance+Multiple.application..Y.N. +Biome+Seeded..Y.N. +Mixture..yes.no., data=ROM.kgm, method="h",
#                level=1, crit="aicc", confsetsize=model1.size)

# par(mfrow=c(1,1))

#  tiff("prod.import.tiff", width = 16, height= 12, units ='cm', res=300)
 plot(model1, type="s") # Burn, year, disturbance, mixture, seeded most important terms, exclude multiple application, burn, and seeded
#  dev.off()

#best fit from output, excluding multiple application not considered important criteria according to plot
model1.bestfit1<-rma.mv(yi, vi, mods = ~ Severe.Disturbance+Biome+level.log+year.log, random = rand.var,data = ROM.kgm,slab = paste(author, year), method="ML")
# #remove correlations
model1.bestfit1a<-rma.mv(yi, vi, mods = ~Biome+level.log+year.log, random = rand.var,data = ROM.kgm,slab = paste(author, year), method="ML")
model1.bestfit1b<-rma.mv(yi, vi, mods = ~Severe.Disturbance+Biome+year.log, random = rand.var,data = ROM.kgm,slab = paste(author, year), method="ML")

models1 = list(model1.bestfit1, model1.bestfit1a, model1.bestfit1b)
# model fit diagnostics: pseudo R2, marginal/conditional R2, AICc, and I2
pseudo.r2 <- c() 
marg.r2 <- c()
cond.r2 <- c()
aicc <- c()
I2 <- c()

for (i in 1:length(models1)) {
  # pseudo R2 of proportional reduction of variance explained when fitting reduced model relative to model with random effects only
  pseudo.r2[i] <- (sum(rma.random$sigma2) - sum(models1[[i]]$sigma2)) / sum(rma.random$sigma2)
  
  # marginal and conditional R2 for linear mixed-models1 from Nakagawa and Scheilzeth 2013
  # variance from fixed effects
  fe.total <- c()
  for (j in 2:ncol(model.matrix(models1[[i]]))) {
    fe <- models1[[i]]$b[j] * model.matrix(models1[[i]])[, j] 
    fe.total <- c(fe.total, fe)
  }
  v.fix <- var(fe.total)
  # variance from random effects
  v.rand <- sum(models1[[i]]$sigma2)
  # variance from residuals
  v.resid <- var(residuals(models1[[i]]))
  # marginal R2 - total variance explained from fixed effects
  marg.r2[i] <- v.fix / (v.fix + v.rand + v.resid)
  # conditional R2 - total variance explained from fixed and random effects
  cond.r2[i] <- (v.fix + v.rand) / (v.fix + v.rand + v.resid)
  
  # from Wolfgang Viechtbauer (http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate)
  # I2 statistic - amount of heterogeneity relative to the total amount of variance in observed effects (how much heterogeneity contributes to total variance)
  W <- diag(1/ROM.kgm$vi)
  X <- model.matrix(models1[[i]])
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2[i] <- 100 * sum(models1[[i]]$sigma2) / (sum(models1[[i]]$sigma2) + (models1[[i]]$k - models1[[i]]$p) / sum(diag(P)))
  
  aicc[i] <- fitstats(models1[[i]])[5]
}



 models1.df<-cbind(pseudo.r2,marg.r2,cond.r2,aicc,I2)
write.csv(models1.df,file="model.biome.prod.csv")


#
#model 2 uses precip/temp not biom
# model2.size <- glmulti(yi ~ level.log + year.log +   Severe.Disturbance+Multiple.application..Y.N. + temp.sq+precip.log+Seeded..Y.N. +Mixture..yes.no., data=ROM.kgm, method="d",
#                        level=1, crit="aicc")
model2 <- glmulti(yi ~ level.log + year.log +   Severe.Disturbance+Multiple.application..Y.N. + temp.sq+precip.log+Seeded..Y.N. +Mixture..yes.no., data=ROM.kgm, method="h",
                   level=1, crit="aicc", confsetsize=model2.size)

# tiff("prod.import2.tiff", width = 16, height= 12, units ='cm', res=300)
 # plot(model2, type="s") #  exclude multiple application, disturbance, and mixture
# dev.off()

model2.bestfit<-rma.mv(yi, vi, mods = ~ Seeded..Y.N.+level.log+year.log+temp.sq+precip.log, random = rand.var,data = ROM.kgm,slab = paste(author, year), method="ML")

#remove correlations, most resulted in the same models as above
model2.bestfita<-rma.mv(yi, vi, mods = ~ Seeded..Y.N.+level.log+year.log, random = rand.var,data = ROM.kgm,slab = paste(author, year), method="ML")
model2.bestfitb<-rma.mv(yi, vi, mods = ~ level.log+year.log+temp.sq+precip.log, random = rand.var,data = ROM.kgm,slab = paste(author, year),method = "ML")


models2 = list(model2.bestfit, model2.bestfita, model2.bestfitb)
# model fit diagnostics: pseudo R2, marginal/conditional R2, AICc, and I2
pseudo.r2 <- c()
marg.r2 <- c()
cond.r2 <- c()
aicc <- c()
I2 <- c()

for (i in 1:length(models2)) {
  # pseudo R2 of proportional reduction of variance explained when fitting reduced model relative to model with random effects only
  pseudo.r2[i] <- (sum(rma.random$sigma2) - sum(models2[[i]]$sigma2)) / sum(rma.random$sigma2)
  
  # marginal and conditional R2 for linear mixed-models2 from Nakagawa and Scheilzeth 2013
  # variance from fixed effects
  fe.total <- c()
  for (j in 2:ncol(model.matrix(models2[[i]]))) {
    fe <- models2[[i]]$b[j] * model.matrix(models2[[i]])[, j] 
    fe.total <- c(fe.total, fe)
  }
  v.fix <- var(fe.total)
  # variance from random effects
  v.rand <- sum(models2[[i]]$sigma2)
  # variance from residuals
  v.resid <- var(residuals(models2[[i]]))
  # marginal R2 - total variance explained from fixed effects
  marg.r2[i] <- v.fix / (v.fix + v.rand + v.resid)
  # conditional R2 - total variance explained from fixed and random effects
  cond.r2[i] <- (v.fix + v.rand) / (v.fix + v.rand + v.resid)
  
  # from Wolfgang Viechtbauer (http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate)
  # I2 statistic - amount of heterogeneity relative to the total amount of variance in observed effects (how much heterogeneity contributes to total variance)
  W <- diag(1/ROM.kgm$vi)
  X <- model.matrix(models2[[i]])
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2[i] <- 100 * sum(models2[[i]]$sigma2) / (sum(models2[[i]]$sigma2) + (models2[[i]]$k - models2[[i]]$p) / sum(diag(P)))
  
  aicc[i] <- fitstats(models2[[i]])[5]
}

models2.df<-cbind(pseudo.r2,marg.r2,cond.r2,aicc,I2)
     # write.csv(models2.df,file="model.t.p.prod.df.csv")

#include interactions with top model
mod1<-rma.mv(yi, vi, mods = ~Biome*level.log*year.log, random = rand.var,data = ROM.kgm,slab = paste(author, year), method = "ML")


mods = list(mod1, mod2)
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
  W <- diag(1/ROM.kgm$vi)
  X <- model.matrix(mods[[i]])
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2[i] <- 100 * sum(mods[[i]]$sigma2) / (sum(mods[[i]]$sigma2) + (mods[[i]]$k - mods[[i]]$p) / sum(diag(P)))
  
  aicc[i] <- fitstats(mods[[i]])[5]
}
# 
mods.int.df<-cbind(pseudo.r2,marg.r2,cond.r2,aicc,I2)
# write.csv(mods.int.df,file="mods.int.prod.csv")

# par(mfrow=c(1,1))

#Models selected output
model1.bestfit1a
summary(model1.bestfit1a)
summary(mod1)
funnel(mod1)
forest(mod1)


ROM.kgm$Biome

mod1.shrub<-rma.mv(yi, vi, mods = ~ level.log*year.log, random = rand.var,data = ROM.kgm,slab = paste(author, year), subset = Biome=="Shrubland", method = "ML")
mod1.forest<-rma.mv(yi, vi, mods = ~ level.log*year.log, random = rand.var,data = ROM.kgm,slab = paste(author, year), subset = Biome=="Temperate seasonal forest", method = "ML")
mod1.grass<-rma.mv(yi, vi, mods = ~ level.log*year.log, random = rand.var,data = ROM.kgm,slab = paste(author, year), subset = Biome=='Temperate grassland', method = "ML")


# dis.yi.y<-ROM.kgm%>%
#   filter(Severe.Disturbance=="Y")%>%
#   summarise(mean(yi), sd(yi))
# #View(dis.yi.y)
# dis.yi.n<-ROM.kgm%>%
#   filter(Severe.Disturbance=="N")%>%
#   summarise(mean(yi), sd(yi))

#=========================================graphs=========


mod<-c("Biome-Shrubland","Biome-Temperate Grassland", "Biome-Temperate Forest")
mod.lev<-c("Biome-Shrubland","Biome-Temperate Grassland", "Biome-Temperate Forest")
mod.cat<-c("Biome", "Biome","Biome")
df.cat<-as.data.frame(cbind(mod,mod.cat,df))

#write.csv(as.data.frame(coef(summary(mod1))), file="prod.mod.csv")
df.prod<-read.csv("prod.mod.csv")
head(df.prod)

df.biome<-df.prod[1:3,]

biome.colors2<-c("blue4","chartreuse4","orange")
LRR.biome.prod = ggplot(data=df.biome,
                     aes(x = X, y = estimate, ymin =ci.lb , ymax = ci.ub ))+
  geom_point(aes(color=X))+
  geom_hline(yintercept =0, linetype=2)+
  scale_color_manual(values =biome.colors, name = "Biome", labels=c("Temperate Grassland" , "Shrubland","Temperate Seasonal Forest"),
                     limits=c("BiomeTemperate grassland", "intrcpt", "BiomeTemperate seasonal forest"))+
  xlab('Biome')+ ylab("Log Response Ratio")+
  geom_errorbar(aes(ymin=ci.lb, ymax=ci.ub,col=X),width=0.3,cex=1)+ 
  #xlim(breaks =rev(mod.lev))+
  #facet_wrap(~X,scale="free") +
  theme(plot.title=element_text(size=16,face="bold"),
        axis.text.x=element_text(face="bold"),
        axis.title=element_text(size=12,face="bold"),
        legend.position = c(.8,.7), axis.text.y=element_blank())+
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=20),
        axis.title=element_text(size=20),legend.title=element_text(size=15), 
        legend.text=element_text(size=13))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))  +

  
  coord_flip()
LRR.biome.prod

# tiff("prod.biome.tiff", width = 20, height= 12, units ='cm', res=600)
# LRR.biome.prod
# dev.off()


#level
preds<-predict(mod1, addx = TRUE)
#preds<-do.call(cbind.data.frame, preds)

names(preds)

preds<-as.data.frame(preds)
names(preds)
preds$level.exp<-exp(preds[,10]) # convert levels fROM.kgm log
ROM.kgm$level.exp<-exp(ROM.kgm$level.log)

preds$year.exp<-exp(preds[,11]) # convert years fROM.kgm log
ROM.kgm$year.exp<-exp(ROM.kgm$year.log)

dim(ROM.kgm)


#write.csv(preds, file = "lrr.prod.csv")
lrr.prod<-read.csv("lrr.prod.csv", header = TRUE)

biome.colors<-c("blue4","orange"
                ,"chartreuse4")
biome.line<-c("twodash","dashed","solid")

names(lrr.prod)
#forest plot
lrr.prod$author.year<-paste(lrr.prod$author,lrr.prod$year, sep="-")

<<<<<<< HEAD:Datasheets for analysis/biomass.R
forest.prod<-viz_forest(x =coverpreds[c("pred", "se")], variant="classic",study_labels= coverpreds$author.year,summary_label ="Summary Effect", xlab = "Log Response Ratio")
tiff("prod.forest.tiff", width = 16, height= 25, units ='cm', res=600)
forest.prod
dev.off()
=======
forest.prod<-viz_forest(x =lrr.prod[c("pred", "se")], variant="classic",study_labels= lrr.prod$author.year,summary_label ="Summary Effect", xlab = "Log Response Ratio")
# tiff("prod.forest.tiff", width = 16, height= 25, units ='cm', res=600)
# forest.prod
# dev.off()
>>>>>>> 32f4e8f48d171a9b89005bbcb2bea5a866449a68:Datasheets for analysis/wc10/biomass.R

g0 <- ggplot(preds,aes(x=level.exp,y=pred))+ 
  stat_smooth(method="glm",formula=y~x,fullrange=T,se=FALSE,size=1)+
  geom_hline(yintercept=0,  size=1)+
  stat_smooth(data=preds,aes(x=level.exp,y=ci.lb),method="glm",formula=y~x,fullrange=T,
              se=FALSE,size=1,linetype=3,colour="blue")+
  stat_smooth(data=preds,aes(x=level.exp,y=ci.ub),method="glm",formula=y~x,fullrange=T,
              se=FALSE,size=1,linetype=3,colour="blue")+
  geom_point(data=ROM.kgm, aes(x=level.exp, y = yi, size=1/vi))+#, color=Biome, group=Biome, fill=Biome, shape=Biome))+
  ylab("Log Response Ratio")+
  xlab(expression(Level~of~Biosolids~Applied~(Mg~ha^{-1})))+
  #scale_color_manual(values=biome.colors)+
  scale_y_continuous(breaks=seq(-1, 6, 1))+
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=20),
        axis.title=element_text(size=20), legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))  

g0
# tiff("prod.level.tiff", width = 16, height= 12, units ='cm', res=600)
# g0
# dev.off()


#level by biome

preds.s<-predict(mod1.shrub, addx = TRUE)
#preds<-do.call(cbind.data.frame, preds)
preds.s<-as.data.frame(preds.s)
preds.s$level.exp<-exp(preds.s$X.level.log) # convert levels fROM.kgm log
preds.s$Biome<-rep("Shrubland", times = 111)

dim(preds.f)
preds.g<-predict(mod1.grass, addx = TRUE)
#preds<-do.call(cbind.data.frame, preds)
preds.g<-as.data.frame(preds.g)
preds.g$level.exp<-exp(preds.g$X.level.log) # convert levels fROM.kgm log
preds.g$Biome<-rep("Temperate grassland", times = 82)

preds.f<-predict(mod1.forest, addx = TRUE)
#preds<-do.call(cbind.data.frame, preds)
preds.f<-as.data.frame(preds.f)
preds.f$level.exp<-exp(preds.f$X.level.log) # convert levels fROM.kgm log
preds.f$Biome<-rep("Temperate seasonal forest", times = 48)

preds.biome=rbind(preds.g,preds.f, preds.s)

# ROM.s<-ROM.kgm%>%
#   filter(Biome=="Shrubland")
# ROM.g<-ROM.kgm%>%
#   filter(Biome=="Temperate grassland")
# ROM.f<-ROM.kgm%>%
#   filter(Biome=="Temperate seasonal forest")
biome.colors<-c("orange"
                ,"blue4","chartreuse4")
g0.biome <- ggplot(preds,aes(x=level.exp,y=pred, size=1/se))+ 
  stat_smooth(method="glm",formula=y~x,fullrange=T,se=FALSE,size=1)+
  geom_hline(yintercept=0,  size=1)+
  stat_smooth(data=preds,aes(x=level.exp,y=ci.lb),method="glm",formula=y~x,fullrange=T,
              se=FALSE,size=1,linetype=3,colour="black")+
  stat_smooth(data=preds,aes(x=level.exp,y=ci.ub),method="glm",formula=y~x,fullrange=T,
              se=FALSE,size=1,linetype=3,colour="black")+
  geom_point(data=ROM.kgm, aes(x=level.exp, y = yi, size=1/vi, color=Biome, group=Biome, fill=Biome))+
  ylab("Log Response Ratio")+
  xlab(expression(Level~of~Biosolids~Applied~(Mg~ha^{-1})))+
  scale_color_manual(values=biome.colors)+
  scale_y_continuous(breaks=seq(-1, 6, 1))+
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=20),
        axis.title=element_text(size=20))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))  

g0.biome




bubble.level.biome=ggplot(preds.biome, aes(x=level.exp,y=pred))+
  stat_smooth(data=preds.biome,aes(color= Biome, fill=Biome),method="glm",formula=y~x,fullrange=T,se=TRUE,size=1     )+ #, linetype=Biome 
  geom_hline(yintercept=0, linetype="dashed", size=.5)+
  # stat_smooth(data=preds.biome,aes(x=level.exp,y=ci.lb, color=Biome),method="glm",formula=y~x,fullrange=T,
  #             se=FALSE,size=1,linetype=3)+
  # stat_smooth(data=preds.biome,aes(x=level.exp,y=ci.ub),method="glm",formula=y~x,fullrange=T,
  #             se=FALSE,size=1,linetype=3)+
  # scale_linetype_manual(values=biome.line)+
  geom_point(data=ROM.kgm, aes(x=level.exp, y = yi, size=1/vi, color=Biome, group=Biome, fill=Biome, shape=Biome))+
  ylab("Log Response Ratio")+
  xlab(expression(Level~of~Biosolids~Applied~(Mg~ha^{-1})))+
  scale_color_manual(values=biome.colors)+
  scale_fill_manual(values=biome.colors)+
  
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=20),
        axis.title=element_text(size=20), legend.position = "bottom")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
bubble.level.biome

# tiff("prod.level.biome.tiff", width = 16, height= 12, units ='cm', res=600)
# bubble.level.biome
# dev.off()

ROM.kgm$year.exp<-exp(ROM.kgm$year.log)
preds.biome$year.exp<-exp(preds.biome$X.year.log)

year.prod <- ggplot(preds,aes(x=year.exp,y=pred))+ 
  stat_smooth(method="glm",formula=y~x,fullrange=T,se=FALSE,size=1)+
  geom_hline(yintercept=0,  size=1)+
  stat_smooth(data=preds,aes(x=year.exp,y=ci.lb),method="glm",formula=y~x,fullrange=T,
              se=FALSE,size=1,linetype=3,colour="blue")+
  stat_smooth(data=preds,aes(x=year.exp,y=ci.ub),method="glm",formula=y~x,fullrange=T,
              se=FALSE,size=1,linetype=3,colour="blue")+
  geom_point(data=ROM.kgm, aes(x=year.exp, y = yi, size=1/vi))+#, color=Biome, group=Biome, fill=Biome, shape=Biome))+
  ylab("Log Response Ratio")+
  xlab(expression(year~of~Biosolids~Applied~(Mg~ha^{-1})))+
  #scale_color_manual(values=biome.colors)+
  scale_y_continuous(breaks=seq(-1, 6, 1))+
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=20),
        axis.title=element_text(size=20), legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))  

year.prod

# tiff("prod.year.tiff", width = 16, height= 12, units ='cm', res=600)
# year.prod
# 
# dev.off()

# View(ROM.kgm)
bubble.year.biome=ggplot(preds.biome, aes(x=year.exp,y=pred))+
  stat_smooth(data=preds.biome,aes(color= Biome, fill=Biome),method="glm",formula=y~x,fullrange=T,se=TRUE,size=1)+ #, linetype=Biome 
  geom_hline(yintercept=0, linetype="dashed", size=.5)+
  # stat_smooth(data=preds.biome,aes(x=year.exp,y=ci.lb, color=Biome),method="glm",formula=y~x,fullrange=T,
  #             se=FALSE,size=1,linetype=3)+
  # stat_smooth(data=preds.biome,aes(x=year.exp,y=ci.ub),method="glm",formula=y~x,fullrange=T,
  #             se=FALSE,size=1,linetype=3)+
  # scale_linetype_manual(values=biome.line)+
  geom_point(data=ROM.kgm, aes(x=year.exp, y = yi, size=1/vi, color=Biome, group=Biome, fill=Biome, shape=Biome))+
  ylab("Log Response Ratio")+
  xlab("Years since reclamation")+
  scale_color_manual(values=biome.colors)+
  scale_fill_manual(values=biome.colors)+
  
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=20),
        axis.title=element_text(size=20), legend.position = "bottom")+
  ylim(-2,6)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
bubble.year.biome

# tiff("prod.year.biome.tiff", width = 16, height= 12, units ='cm', res=600)
# bubble.year.biome
# dev.off()

dim(mod1)


library(cowplot)
library(ggpubr)

legend.prod<-get_legend(bubble.year.biome)
legend.prod.plot<-as_ggplot(legend.prod)



prod.plot<-plot_grid(bubble.year.biome, bubble.level.biome,
                     labels = c("a", "b"),
                     ncol = 2, nrow = 1)

# tiff("biomes.prod.plot.tiff", width = 28, height= 12, units ='cm', res=600)
# prod.plot
#  dev.off()

# tiff("biomes.prod.legend.tiff", width = 20, height= 12, units ='cm', res=600)
# legend.prod.plot
# dev.off()


