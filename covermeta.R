#Cover meta-analysis
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
library(polycor)
library(glmulti)

cover<- read.csv("~/Documents/GitHub/Reclamation_biosolid/Datasheets for analysis/cover.csv")
cover$Mc2<-ifelse(cover$Mc==0, 1/2, cover$Mc)
cover$Temp2<-ifelse(cover$Temp<=0, 1/2, cover$Temp)
ROM.cover <- escalc(n2i = Nc, n1i = Ne, m2i = Mc2, m1i = Me, 
                   sd2i = scimp, sd1i = seimp, data = cover, measure = "ROM", 
                   append = TRUE )
write.csv(ROM.rich, "ROM.cover.csv")

HG.cover <- escalc(n2i = Nc, n1i = Ne, m2i = Mc2, m1i = Me,
                  sd2i = scimp, sd1i = seimp, data = cover, measure = "SMD",
                  append = TRUE) 
# # Check normality of Hedges' vs. log response ratio
par(mfrow=c(1,2))
#
Hedges = HG.cover$yi
# #
h <- hist(Hedges, breaks = 10, density = 10,
          xlab = "Hedges' (g)")
xfit <- seq(min(Hedges,na.rm=TRUE), max(Hedges,na.rm=TRUE), length = 40)
yfit <- dnorm(xfit, mean = mean(Hedges,na.rm=TRUE), sd = sd(Hedges,na.rm=TRUE))
yfit <- yfit * diff(h$mids[1:2]) * length(Hedges)
# #
lines(xfit, yfit, col = "black", lwd = 2)
text(x= -17, y = 50, labels = "B", xpd = NA)
# #
LRR = ROM.cover$yi
# #
r <- hist(LRR, breaks = 10, density = 10, xlab = "Log Response Ratio")
xfit.r <- seq(min(LRR), max(LRR), length = 40)
yfit.r <- dnorm(xfit.r, mean = mean(LRR), sd = sd(LRR))
yfit.r <- yfit.r * diff(r$mids[1:2]) * length(LRR)
# #
lines(xfit.r, yfit.r, col = "black", lwd = 2)

rand.varc=list(~ 1|Paper.ID./Case)
globaleffcover<- rma.mv(yi = yi, V = vi, random = rand.varc, 
                        data = ROM.cover, method = "ML")
summary(globaleffcover)

##subset: SHRUBLAND

# subset of data with only shrubland to feed to rma
shrubland <- which(ROM.cover$Biome== "Shrubland")
cover.shrubland <- ROM.cover[shrubland, ]
covereffshrub <- rma.mv(yi = yi, V = vi, random = rand.varc, 
                       data = cover.shrubland, method = "ML")
summary(covereffshrub)

##subset: GRASSLAND
# subset of data with only grassland to feed to rma
grassland <- which(ROM.cover$Biome== "Temperate grassland")
cover.grassland <- ROM.cover[grassland, ]
covereffgrass <- rma.mv(yi = yi, V = vi, random = rand.varc, 
                       data = cover.grassland, method = "ML")
summary(covereffgrass)

##subset: FOREST
forest <- which(ROM.cover$Biome== "Temperate seasonal forest")
cover.forest <- ROM.cover[forest, ]
covereffforest<- rma.mv(yi = yi, V = vi, random = rand.varc, 
                       data = cover.forest, method = "ML")
summary(covereffforest)

ROM.cover$level.log<-log(ROM.cover$Biosolid.level)
ROM.cover$year.log<-log(ROM.cover$yeartrans)
ROM.cover$temp.log<-log(ROM.cover$Temp2)
ROM.cover$precip.log<-log(ROM.cover$Precip)

# Moderator analysis \ metaregression with single covariates ---------
#
# generate using rma (random effects model with single moderator 
# and no intercept)
rma.random.c <- rma.mv(yi = yi, V = vi, random = rand.varc, data = ROM.cover, 
                     method = "ML")
summary(rma.random.c)
ran.dfc<-as.data.frame(coef(summary(rma.random.c)))

rma.lev.c<-rma.mv(yi,vi,mods = ~ level.log,random =rand.varc,
                data = ROM.cover, method ="ML")
lev.dfc<-as.data.frame(coef(summary(rma.lev.c)))

rma.yr.c<-rma.mv(yi,vi,mods = ~ year.log, random =rand.varc,
               data = ROM.cover,method ="ML")
yr.dfc<-as.data.frame(coef(summary(rma.yr.c)))

rma.temp.c<-rma.mv(yi,vi,mods = ~ temp.log, random =rand.varc,
                 data = ROM.cover, method ="ML")
temp.dfc<-as.data.frame(coef(summary(rma.temp.c)))

rma.precip.c<-rma.mv(yi,vi,mods = ~  precip.log, random =rand.varc,
                   data = ROM.cover, method ="ML")
precip.dfc<-as.data.frame(coef(summary(rma.precip.c)))

rma.mix.c<-rma.mv(yi,vi,mods = ~ -1 + Mixture, random =rand.varc,
                data = ROM.cover, method ="ML")
mix.dfc<-as.data.frame(coef(summary(rma.mix.c)))

rma.dis.c<-rma.mv(yi,vi,mods = ~ -1 + S.dist, random =rand.varc,
                data = ROM.cover, method ="ML")
dis.dfc<-as.data.frame(coef(summary(rma.dis.c)))

rma.mult.c<-rma.mv(yi,vi,mods = ~ -1 + Multiple.application, random =rand.varc,
                 data = ROM.cover, method ="ML")
mult.dfc<-as.data.frame(coef(summary(rma.mult.c)))

rma.seed.c<-rma.mv(yi,vi,mods = ~ -1 + Seeded, random =rand.varc,
                 data = ROM.cover, method ="ML")
seed.dfc<-as.data.frame(coef(summary(rma.seed.c)))

rma.biome.c<-rma.mv(yi, vi, mods = ~ -1 + Biome, random = rand.varc, 
                    data = ROM.cover, method = "ML")
biome.dfc<-as.data.frame(coef(summary(rma.biome.c)))

simple.cover.models<-rbind(ran.dfc, lev.dfc, yr.dfc, temp.dfc, precip.dfc, mix.dfc, 
                           dis.dfc, mult.dfc, seed.dfc, biome.dfc)
AIC.rma(rma.random.c, rma.lev.c, rma.yr.c, rma.temp.c, rma.precip.c, rma.mix.c, 
    rma.dis.c, rma.mult.c, rma.seed.c, rma.biome.c, correct = TRUE)
write.csv(simple.cover.models, "simple.cover.models.csv")

# Explore variable correlations ---------------------------------------

pairs(~ level.log + year.log + S.dist +  Multiple.application + Mixture + Biome  + Seeded + Temp + precip.log, data = ROM.cover)
library(polycor)


var<-ROM.cover%>%
  dplyr::select(level.log, year.log, S.dist, Multiple.application, Mixture, Biome, Seeded , temp.log, precip.log)

correlations.c <- hetcor(data = var, std.err = FALSE)

# including intercept to get an accurate Q-test for heterogeneity
rma.mv(yi,vi,mods = ~ level.log,random =rand.varc, data = ROM.cover,
       method ="ML")
rma.mv(yi,vi,mods = ~ year.log, random =rand.varc, data = ROM.cover,
       method ="ML")
rma.mv(yi,vi,mods = ~ temp.log, random = rand.varc,data = ROM.cover, 
       method ="ML")
rma.mv(yi,vi,mods = ~  precip.log, random =rand.varc,data = ROM.cover, 
       method ="ML")
rma.mv(yi,vi,mods = ~  Mixture, random =rand.varc, data = ROM.cover,
       method ="ML")
rma.mv(yi,vi,mods = ~ S.dist, random =rand.varc,
       data = ROM.cover, method ="ML")
rma.mv(yi,vi,mods = ~ Multiple.application, random =rand.varc, data = ROM.cover,
       method ="ML")

rma.mv(yi,vi,mods = ~ Seeded, random =rand.varc,
       data = ROM.cover, method ="ML")
rma.mv(yi, vi, mods = ~Biome, random = rand.varc, 
       data = ROM.cover, method = "ML")
AIC(rma.mv(yi, vi, mods = ~ Biome, random =rand.varc,
       data= ROM.cover, method = "ML"))

AIC.rma(rma.random.r, rma.lev.r, rma.yr.r, rma.temp.r, rma.precip.r, rma.mix.r, 
        rma.dis.r, rma.mult.r, rma.seed.r, rma.biome.r, correct = TRUE)

rma.glmulti <- function(formula, data, ...)
  rma(formula, vi, data=data, method="ML", ...)

# run glmulti with meta regression function with biome but no temp and precip
candidatesetc1 <- glmulti(yi ~ S.dist + Multiple.application  + Seeded +
                          Mixture + year.log + level.log + Biome, V = "vi",
                        random = "rand.varc", data = ROM.cover, level = 1, method = "d", 
                        fitfunction = rma.glmulti, crit = "aicc")

# run glmulti with meta regression function
covermodels1 <- glmulti(yi ~ S.dist + Multiple.application + Seeded +
                            Mixture + year.log + level.log + Biome, V = "vi",
                          random = "rand.var", data = ROM.cover, level = 1, method = "h", 
                          fitfunction = rma.glmulti, crit = "aicc", confsetsize = candidatesetc1)

# run glmulti with meta regression function with no biome but temp and precip
candidatesetc2 <- glmulti(yi ~ S.dist + Multiple.application  + Seeded +
                            Mixture + year.log + level.log + temp.log + precip.log, V = "vi",
                          random = "rand.varc", data = ROM.cover, level = 1, method = "d", 
                          fitfunction = rma.glmulti, crit = "aicc")

# run glmulti with meta regression function
covermodels2 <- glmulti(yi ~ S.dist + Multiple.application + Seeded +
                          Mixture + year.log + level.log + temp.log + precip.log, V = "vi",
                        random = "rand.var", data = ROM.cover, level = 1, method = "h", 
                        fitfunction = rma.glmulti, crit = "aicc", confsetsize = candidatesetc2)

print(covermodels1)
bestcovermodel<-summary(covermodels1@objects[[1]])
plot(covermodels1, type="s")

print(covermodels2)
bestcovermodel<-summary(covermodels2@objects[[1]])
plot(covermodels2, type="s")

#Best model from output (precip and temp model), excluding covariates which plot indicate as unimportant
bestcov1<-rma.mv(yi, vi, mods = ~S.dist + Multiple.application + Seeded + Mixture +
          temp.log + precip.log , random = rand.varc,
          data = ROM.cover,slab = paste(author, year), method = "ML")

#subsets of best model removing all correlations greater than 0.3
bestcov2<-rma.mv(yi, vi, mods = ~ S.dist + Multiple.application + Mixture + temp.log, random = rand.varc,
                 data = ROM.cover,slab = paste(author, year), method = "ML")

bestcov3<-rma.mv(yi, vi, mods = ~ Seeded + Multiple.application + temp.log, random = rand.varc,
                 data = ROM.cover,slab = paste(author, year), method = "ML")

bestcov4 <- rma.mv(yi, vi, mods = ~ precip.log + Multiple.application + temp.log, random = rand.varc,
                   data= ROM.cover, slab = paste(author, year), method = "ML")

##Run the best model with correlated variables below 0.3 and interactions 
bestcov3int<- rma.mv(yi, vi, mods = ~ temp.log * Seeded * Multiple.application, random = rand.varc,
                  data = ROM.cover,slab = paste(author, year), method = "ML")


covermodelslist = list(bestcov1, bestcov2, bestcov3, bestcov4, bestcov3int)

# model fit diagnostics: pseudo R2, marginal/conditional R2, AICc, and I2
pseudo.r2c <- c()
marg.r2c <- c()
cond.r2c <- c()
aiccc <- c()
I2c <- c()

for (i in 1:length(covermodelslist)) {
  # pseudo R2 of proportional reduction of variance explained when fitting reduced model relative to model with random effects only
  pseudo.r2c[i] <- (sum(rma.random.c$sigma2) - sum(covermodelslist[[i]]$sigma2)) / sum(rma.random.c$sigma2)
  
  # marginal and conditional R2 for linear mixed-models1 from Nakagawa and Scheilzeth 2013
  # variance from fixed effects
  fe.total <- c()
  for (j in 2:ncol(model.matrix(covermodelslist[[i]]))) {
    fe <- covermodelslist[[i]]$b[j] * model.matrix(covermodelslist[[i]])[, j] 
    fe.total <- c(fe.total, fe)
  }
  v.fix <- var(fe.total)
  # variance from random effects
  v.rand <- sum(covermodelslist[[i]]$sigma2)
  # variance from residuals
  v.resid <- var(residuals(covermodelslist[[i]]))
  # marginal R2 - total variance explained from fixed effects
  marg.r2c[i] <- v.fix / (v.fix + v.rand + v.resid)
  # conditional R2 - total variance explained from fixed and random effects
  cond.r2c[i] <- (v.fix + v.rand) / (v.fix + v.rand + v.resid)
  
  # from Wolfgang Viechtbauer (http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate)
  # I2 statistic - amount of heterogeneity relative to the total amount of variance in observed effects (how much heterogeneity contributes to total variance)
  W <- diag(1/ROM.cover$vi)
  X <- model.matrix(covermodelslist[[i]])
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2c[i] <- 100 * sum(covermodelslist[[i]]$sigma2) / (sum(covermodelslist[[i]]$sigma2) + (covermodelslist[[i]]$k - 
                                                                                            covermodelslist[[i]]$p) / sum(diag(P)))
  
  aiccc[i] <- fitstats(covermodelslist[[i]])[5]
}

models.cov<-cbind(pseudo.r2c,marg.r2c,cond.r2c,aiccc,I2c)
write.csv(models.cov,file="cover models.csv")

funnel(bestcov3int)
forest(bestcovint1)

#nothing significant so I am not sure what to plot
as.data.frame(biome.dfc)
#Additional unused code
#topcover <- weightable(covermodels)
#topcover <- topcover[topcover$aicc <= min(topcover$aicc) + 2,]
#topcover

#eval(metafor:::.glmulti)

#coverweights <- round(coef.glmulti(covermodels, select = nrow(topcover), 
#                                  icmethod = "Burnham"), 7)
#coverweights <- coverweights[, c(1, 4, 5)]
#coverweights
#cummulative metanalysis with forest plot
#forest(cumul(rma.uni(yi, vi, data= ROM.cover)), transf = exp, order=order(ROM.cover$Paper.ID.))

covertable <- data.frame(variate = c(rep(x = "Biosolid level", length(rma.lev.c$b)),
                                        rep(x = "Mixture", length(rma.mix.c$b)), 
                                        rep(x = "Disturbance", length(rma.dis.c$b)), 
                                        rep(x = "Multiple application", length(rma.mult.c$b)),
                                        rep(x = "Seeded", length(rma.seed.c$b)), 
                                        rep(x = "Temperature", length(rma.temp.c$b)), 
                                        rep(x = "Precipitation", length(rma.precip.c$b)),
                                        rep(x = "Biome", length(rma.biome.c$b)),
                                        rep(x = "Year", length(rma.yr.c$b))),
                            LOR = c(rma.lev.c$b, rma.mix.c$b, rma.dis.c$b, rma.mult.c$b,
                                    rma.seed.c$b, rma.temp.c$b, rma.precip.c$b, rma.biome.c$b, rma.yr.c$b),
                            CI.LB = c(rma.lev.c$ci.lb, rma.mix$ci.lb, rma.dis.c$ci.lb, rma.mult.c$ci.lb,
                                      rma.seed.c$ci.lb, rma.temp.c$ci.lb, rma.precip.c$ci.lb,
                                      rma.biome.c$ci.lb, rma.yr.c$ci.lb), 
                            CI.UB = c(rma.lev.c$ci.ub, rma.mix$ci.ub, rma.dis.c$ci.ub, rma.mult.c$ci.ub,
                                      rma.seed.c$ci.ub, rma.temp.c$ci.ub, rma.precip.c$ci.ub,
                                      rma.biome.c$ci.ub, rma.yr.c$ci.ub),
                            lab = c("Biosolid level_int", "Biosolid Level_slope", levels(ROM.cover$Mixture), 
                                    levels(ROM.cover$S.dist), levels(ROM.cover$Multiple.application), 
                                    levels(ROM.cover$Seeded), "Temperature_int", "Temperature_slope",
                                    "Precipitation_int", "Precipitation_slope", levels(ROM.cover$Biome), "Year_int", "Year_slope"),
                            N = c(rep(sum(!is.na(ROM.cover$Biosolid.level)), 2), 
                                  summary(ROM.cover$Mixture[!is.na(ROM.cover$Mixture)]), 
                                  summary(ROM.cover$S.dist[!is.na(ROM.cover$S.dist)]),
                                  summary(ROM.cover$Multiple.application[!is.na(ROM.cover$Multiple.application)]),
                                  summary(ROM.cover$Seeded[!is.na(ROM.cover$Seeded)]),
                                  rep(sum(!is.na(ROM.cover$temp.log)), 2),
                                  rep(sum(!is.na(ROM.cover$precip.log)),2),
                                  summary(ROM.cover$Biome[!is.na(ROM.cover$Biome)]),
                                  rep(sum(!is.na(ROM.cover$year.log)),2)))

coverpreds<-predict(bestcov3int, addx = TRUE)
coverpreds<-as.data.frame(coverpreds)
coverpreds$temp.exp<-exp(coverpreds$X.temp.log)
temp.plot <- ggplot(coverpreds, aes(x = temp.exp, y = pred)) + geom_line() +
  geom_line(aes(x = temp.exp, y = ci.lb), linetype = "dashed") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray65", size = 0.2) +
  geom_line(aes(x = temp.exp, y = ci.ub), linetype = "dashed") + 
  geom_rug(sides = "b") + 
  ylab(label = expression(paste("Effect size (", italic("LRR"),")"))) + 
  xlab(expression(paste("Temperature [",degree,"C]"))) +
  scale_x_continuous(breaks = seq(0, 20, 5), limits = c(0, 20)) + 
  ylim(c(-0.4, 2.7)) + theme_classic() + 
  theme(text = element_text(size = 14))
tiff(file = "Tempc.tiff", width = 16, height= 12, units ='cm', res=600)
temp.plot
dev.off()

preds$level.exp<-exp(coverpreds[,10]) # convert levels fROM.kgm log
ROM.kgm$level.exp<-exp(ROM.cover$level.log)

biome.tab <- subset(covertable, variate == "Biome")
biome.tab$lab <- as.character(biome.tab$lab)
biome.levels <- biome.tab[order(biome.tab$LOR), 5]
biome.plotc <- ggplot(data = subset(covertable, variate == "Biome"), 
                     aes(x = factor(lab, levels = biome.levels), y = LOR)) + 
  geom_point() + geom_errorbar(aes(ymax = CI.UB, ymin = CI.LB, 
                                   width = 0.1)) + geom_text(aes(y = -0.3, 
                                                                 label = paste("n = ", N, sep = "")), size = 3) + 
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray65", size = 0.2) + 
  annotate("text", x = 0.7, y = 2.0, label = " ", 
           fontface = "bold") + scale_x_discrete(name = "", 
                                                 labels = c("Shrubland", "Temperate grassland", "Temperate seasonal forest")) + 
  ylab(label = expression(paste("Effect size (", italic("LRR"),")"))) + 
  coord_cartesian(ylim = c(-0.3, 2.3)) + 
  theme_classic() + theme(text = element_text(size = 14)) 
tiff(file = "biomec.tiff", width = 16, height= 12, units ='cm', res=600)
biome.plotc
dev.off()

seedc.tab <- subset(covertable, variate == "Seeded")
seedc.tab$lab <- as.character(seedc.tab$lab)
seedc.levels <- seedc.tab[order(seedc.tab$LOR), 5]
seed.plotc <- ggplot(data = subset(covertable, variate == "Seeded"), 
                      aes(x = factor(lab, levels = seedc.levels), y = LOR)) + 
  geom_point() + geom_errorbar(aes(ymax = CI.UB, ymin = CI.LB, 
                                   width = 0.1)) + geom_text(aes(y = -0.1, 
                                                                 label = paste("n = ", N, sep = "")), size = 3) + 
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray65", size = 0.2) + 
  annotate("text", x = 0.7, y = 2.0, label = " ", 
           fontface = "bold") + scale_x_discrete(name = "", 
                                                 labels = c("Not seeded", "Seeded")) + 
  ylab(label = expression(paste("Effect size (", italic("LRR"),")"))) + 
  coord_cartesian(ylim = c(-0.1, 2.0)) + 
  theme_classic() + theme(text = element_text(size = 14)) 
tiff(file = "seedc.tiff", width = 16, height= 12, units ='cm', res=600)
seed.plotc
dev.off()


temp.levels <- print.list.rma(predict.rma(rma.temp.c))
temp.levels <- data.frame(temp.levels, 
                         temp = ROM.cover$Temp2[which(!is.na(ROM.cover$Temp2), 
                                                               arr.ind = TRUE)])
temp.plot <- ggplot(temp.levels, aes(x = temp, y = pred)) + geom_line() +
  geom_line(aes(x = temp, y = ci.lb), linetype = "dashed") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray65", size = 0.2) +
  geom_line(aes(x = temp, y = ci.ub), linetype = "dashed") + 
  geom_rug(sides = "b") + 
  ylab(label = expression(paste("Effect size (", italic("LRR"),")"))) + 
  xlab(expression(paste("Temperature [",degree,"C]"))) +
  scale_x_continuous(breaks = seq(0, 20, 5), limits = c(0, 20)) + 
  ylim(c(-0.4, 2.7)) + theme_classic() + 
  theme(text = element_text(size = 17))
tiff(file = "Tempc.tiff", width = 16, height= 12, units ='cm', res=600)
temp.plot
dev.off()


prod.plot<-plot_grid(seedr.plot, bubble.level.biome,
                     labels = c("a", "b"),
                     ncol = 2, nrow = 1)


g0 <- ggplot(coverpreds,aes(x=exp(X.temp.log),y=pred))+ 
  stat_smooth(method="glm",formula=y~x,fullrange=T,se=FALSE,size=1)+
  geom_hline(yintercept=0,  size=1)+
  stat_smooth(data=coverpreds,aes(x=exp(X.temp.log),y=ci.lb),method="glm",formula=y~x,fullrange=T,
              se=FALSE,size=1,linetype=3,colour="blue")+
  stat_smooth(data=coverpreds,aes(x=exp(X.temp.log),y=ci.ub),method="glm",formula=y~x,fullrange=T,
              se=FALSE,size=1,linetype=3,colour="blue")+
  geom_point(data=ROM.cover, aes(x=exp(temp.log), y = yi, size=1/vi, color=Seeded, group=Seeded, fill=Seeded))+
ylab("Log Response Ratio")+
  xlab(expression(Level~of~Biosolids~Applied~(Mg~ha^{-1})))+
  scale_y_continuous(breaks=seq(-1, 6, 1))+
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=20),
        axis.title=element_text(size=20))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))  



#seeded by temp
##subset: Seeded

# subset of data with only grassland to feed to rma
seeded <- which(ROM.cover$Seeded== "Y")
noseeded<- which(ROM.cover$Seeded== "N")
cover.seeded <- ROM.cover[seeded, ]
cover.noseeded <- ROM.cover[noseeded, ]

modcov.seeded <- rma.mv(yi, vi, mods = ~ temp.log, random = rand.varc,data = ROM.cover,slab = paste(author, year), subset = Seeded=="Y", method = "ML")
modcov.unseeded<-rma.mv(yi, vi, mods = ~ temp.log, random = rand.varc,data = ROM.cover,slab = paste(author, year), subset = Seeded=="N", method = "ML")


preds.covs<-predict(modcov.seeded, addx = TRUE)
preds.covn<-predict(modcov.unseeded, addx = TRUE)
#preds<-do.call(cbind.data.frame, preds)
preds.s<-as.data.frame(preds.covs)
preds.n<-as.data.frame(preds.covn)
preds.s$temp.exp<-exp(preds.s$X.temp.log)# convert levels from log
preds.n$temp.exp<-exp(preds.n$X.temp.log)
preds.s$Seeded<-rep("Y", times = 81)
preds.n$Seeded<-rep("N", times = 98)

allpredss=rbind(preds.n, preds.s)

# ROM.s<-ROM.kgm%>%
#   filter(Biome=="Shrubland")
# ROM.g<-ROM.kgm%>%
#   filter(Biome=="Temperate grassland")
# ROM.f<-ROM.kgm%>%
#   filter(Biome=="Temperate seasonal forest")
level.colors<-c("black","red")
g0.level <- ggplot(allpredss,aes(x=temp.exp,y=pred, size=1/se))+ 
  stat_smooth(method="glm",formula=y~x,fullrange=T,se=FALSE,size=1)+
  geom_hline(yintercept=0,  size=1)+
  stat_smooth(data=allpredss,aes(x=temp.exp,y=ci.lb),method="glm",formula=y~x,fullrange=T,
              se=FALSE,size=1,linetype=3,colour="black")+
  stat_smooth(data=allpredss,aes(x=temp.exp,y=ci.ub),method="glm",formula=y~x,fullrange=T,
              se=FALSE,size=1,linetype=3,colour="black")+
  geom_point(data=ROM.cover, aes(x=exp(temp.log), y = yi, size=1/vi, color=Seeded, group=Seeded, fill=Seeded))+
  ylab("Log Response Ratio")+
  xlab(expression(Temperature(Mg~ha^{-1})))+
  scale_color_manual(values=level.colors)+
  scale_y_continuous(breaks=seq(-1, 6, 1))+
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=20),
        axis.title=element_text(size=20))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))  



bubble.temps = ggplot(allpredss, aes(x=temp.exp,y=pred))+
  stat_smooth(data=allpredss,aes(color= Seeded, fill=Seeded),method="glm",formula=y~x,fullrange=T,se=TRUE,size=1     )+ #, linetype=Biome 
  geom_hline(yintercept=0, linetype="dashed", size=.5)+
  # stat_smooth(data=preds.biome,aes(x=level.exp,y=ci.lb, color=Biome),method="glm",formula=y~x,fullrange=T,
  #             se=FALSE,size=1,linetype=3)+
  # stat_smooth(data=preds.biome,aes(x=level.exp,y=ci.ub),method="glm",formula=y~x,fullrange=T,
  #             se=FALSE,size=1,linetype=3)+
  # scale_linetype_manual(values=biome.line)+
  geom_point(data=ROM.cover, aes(x= Temp2, y = yi, size=1/vi, color=Seeded, group=Seeded, fill=Seeded, shape=Seeded))+
  ylab("Log Response Ratio")+
  xlab(expression(paste("Temperature [",degree,"C]")))+
  scale_color_manual(values=level.colors)+
  scale_fill_manual(values=level.colors)+
  
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=20),
        axis.title=element_text(size=20), legend.position = "bottom")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
bubble.temps
tiff(file = "bubble.temps.tiff", width = 16, height= 12, units ='cm', res=600)
bubble.temps
dev.off()

coverf<-ROM.cover
coverf$author.year<-paste(ROM.cover$author, coverf$year, sep = "-")
coverpreds$author.year<-coverf$author.year

##### HEAD:Datasheets for analysis/biomass.R
forest.cover<-viz_forest(x = coverpreds[c("pred", "se")], variant="classic",study_labels= coverpreds$author.year,summary_label ="Summary Effect", xlab = "Log Response Ratio")
tiff("cover.forest.tiff", width = 16, height= 25, units ='cm', res=600)
forest.cover
dev.off()