#Richness meta-analysis
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
library(tidyr)
library(readxl)

richness <- read_excel("~/Documents/GitHub/Grassland restoration with biosolids/Grassland-restoration-with-biosolids/data_sheets.xlsx", 
                    sheet = "Richness")
#richness$Mc2<-ifelse(richness$Control.Mc==0, 1/2, richness$Control.Mc)
#richness$Temp2<-ifelse(richness$Temp<=0, 1/2, richness$Temp)
# #justification for use of log response ratio over Hedges'
 # #HG.rich <- escalc(n2i = Nc, n1i = Ne, m2i = Control.Mc , m1i = Me,
 #                    sd2i = Scimp, sd1i = Seimp, data = richness, measure = "SMD",
 #                   append = TRUE)
 # 
ROM.rich <- escalc(n2i = Nc, n1i = Ne, m2i = Control.Mc + .5, m1i = Me, 
                   sd2i = Scimp, sd1i = Seimp, data = richness, measure = "ROM", 
                   append = TRUE )



# 
# # # Check normality of Hedges' vs. log response ratio
#  par(mfrow=c(1,2))
# Hedges = HG.rich$yi
# h <- hist(Hedges, breaks = 10, density = 10,
#  xlab = "Hedges' (g)")
# xfit <- seq(min(Hedges,na.rm=TRUE), max(Hedges,na.rm=TRUE), length = 40)
# yfit <- dnorm(xfit, mean = mean(Hedges,na.rm=TRUE), sd = sd(Hedges,na.rm=TRUE))
# yfit <- yfit * diff(h$mids[1:2]) * length(Hedges)
# # #
# lines(xfit, yfit, col = "black", lwd = 2)
# text(x= -17, y = 27, labels = "C", xpd = NA)
# # #
# LRR = ROM.rich$yi
# r <- hist(LRR, breaks = 10, density = 10, xlab = "Log Response Ratio")
# xfit.r <- seq(min(LRR), max(LRR), length = 40)
# yfit.r <- dnorm(xfit.r, mean = mean(LRR), sd = sd(LRR))
# yfit.r <- yfit.r * diff(r$mids[1:2]) * length(LRR)
# # #
# lines(xfit.r, yfit.r, col = "black", lwd = 2)
# 
# 
# HG.plot<-plot_grid(plot(h), plot(r),
#                    labels = c("C"),
#                    ncol = 2, nrow = 1)

rand.var=list(~ 1|Paper.ID/Case)
globaleffrich <- rma.mv(yi = yi, V = vi, random = rand.var, 
                        data = ROM.rich, method = "ML")
summary(globaleffrich)

# ROM.rich$level.log<-(ROM.rich$Biosolid.level)
# ROM.rich$year.log<-(ROM.rich$year_restoration)
# ROM.rich$temp.log<-(ROM.rich$Temp2)
# ROM.rich$precip.log<-(ROM.rich$Precip)
# ##subset: SHRUBLAND
# write.csv(ROM.rich, "ROM.rich.csv")
# # subset of data with only grassland to feed to rma
# 
# ROM.rich$level.log<-log(ROM.rich$Biosolid.level)
# ROM.rich$year.log<-log(ROM.rich$yeartrans)
# ROM.rich$temp.log<-log(ROM.rich$Temp2)
# ROM.rich$precip.log<-log(ROM.rich$Precip)

# Moderator analysis \ metaregression with single covariates ---------
#
# generate using rma (random effects model with single moderator 
# and no intercept)

rma.random.r <- rma.mv(yi = yi, V = vi, random = rand.var, data = ROM.rich, 
                       method = "ML")

ran.dfr<-as.data.frame(coef(summary(rma.random.r)))

rma.lev.r<-rma.mv(yi,vi,mods = ~ Biosolid.level,random =rand.var,
                  data = ROM.rich, method ="ML")
lev.dfr<-as.data.frame(coef(summary(rma.lev.r)))

rma.yr.r<-rma.mv(yi,vi,mods = ~ yeartrans, random =rand.var,
                 data = ROM.rich,method ="ML")
yr.dfr<-as.data.frame(coef(summary(rma.yr.r)))

rma.temp.r<-rma.mv(yi,vi,mods = ~ Temp, random =rand.var,
                   data = ROM.rich, method ="ML")
temp.dfr<-as.data.frame(coef(summary(rma.temp.r)))

rma.precip.r<-rma.mv(yi,vi,mods = ~ Precip, random =rand.var,
                     data = ROM.rich, method ="ML")
precip.dfr<-as.data.frame(coef(summary(rma.precip.r)))

rma.mix.r<-rma.mv(yi,vi,mods = ~ -1 + Mixture, random =rand.var,
                  data = ROM.rich, method ="ML")
mix.dfr<-as.data.frame(coef(summary(rma.mix.r)))

rma.dis.r<-rma.mv(yi,vi,mods = ~ -1 + S.disturbance, random =rand.var,
                  data = ROM.rich, method ="ML")
dis.dfr<-as.data.frame(coef(summary(rma.dis.r)))

rma.mult.r<-rma.mv(yi,vi,mods = ~ -1 + multapp, random =rand.var,
                   data = ROM.rich, method ="ML")
mult.dfr<-as.data.frame(coef(summary(rma.mult.r)))

rma.seed.r<-rma.mv(yi,vi,mods = ~ -1 + Seeded, random =rand.var,
                   data = ROM.rich, method ="ML")
seed.dfr<-as.data.frame(coef(summary(rma.seed.r)))

rma.ai.r<-rma.mv(yi,vi,mods = ~ ai, random =rand.var,
                   data = ROM.rich, method ="ML")
ai.dfr<-as.data.frame(coef(summary(rma.ai.r)))

simple.richness.models<- rbind(ran.dfr, lev.dfr, yr.dfr, temp.dfr, precip.dfr, mix.dfr, 
                               dis.dfr, mult.dfr, seed.dfr, ai.dfr)

AIC.rma(rma.random.r, rma.lev.r, rma.yr.r, rma.temp.r, rma.precip.r, rma.mix.r, 
        rma.dis.r, rma.mult.r, rma.seed.r, rma.ai.r, correct = TRUE)
write.csv(simple.richness.models, "simple.richness.models.csv")

# including intercept to get an accurate Q-test for heterogeneity
rma.mv(yi,vi,mods = ~ Biosolid.level,random =rand.var, data = ROM.rich,
       method ="ML")
rma.mv(yi,vi,mods = ~ yeartrans, random =rand.var, data = ROM.rich,
       method ="ML")
rma.mv(yi,vi,mods = ~ Temp, random = rand.var,data = ROM.rich, 
       method ="ML")
rma.mv(yi,vi,mods = ~ Precip, random =rand.var,data = ROM.rich, 
       method ="ML")
rma.mv(yi,vi,mods = ~  Mixture, random =rand.var, data = ROM.rich,
       method ="ML")
rma.mv(yi,vi,mods = ~ S.disturbance, random =rand.var,
                data = ROM.rich, method ="ML")
rma.mv(yi,vi,mods = ~ multapp, random =rand.var, data = ROM.rich,
       method ="ML")
rma.mv(yi,vi,mods = ~ Seeded, random =rand.var,
                 data = ROM.rich, method ="ML")
rma.mv(yi, vi, mods = ~ai, random = rand.var, 
                    data = ROM.rich, method = "ML")

# Explore variable correlations ---------------------------------------

# pairs(~ level.log + year.log + S.disturbance +  multapp + Mixture + biome + Seeded + temp.log + precip.log, data = ROM.rich)
# 
# library(polycor)
# var.rich<-ROM.rich%>%
#   dplyr::select(level.log, year.log, S.disturbance, multapp, Mixture, biome, Seeded, temp.log, precip.log)
# correlations.r <- hetcor(data = var.rich, std.err = FALSE)
# # determine the number of reasonable candidate models for glmulti to run 
# #through using method "d"
rma.glmulti <- function(formula, data, ...)
  rma(formula, vi, data=data, method="ML", ...)

# run glmulti with meta regression function with biome but no temp and precip
candidateset1 <- glmulti(yi ~ S.disturbance + multapp + Seeded + Temp + Burned +
 Mixture + yeartrans + Biosolid.level + ai + Precip, data = ROM.rich, level = 1, method = "d", 
 fitfunction = rma.glmulti, crit = "aicc")

# run glmulti with meta regression function with biome but no temp and precip
richnessmodels1 <- glmulti(yi ~ S.disturbance + multapp + Seeded + Temp + Burned +
                             Mixture + yeartrans + Biosolid.level + ai + Precip, data = ROM.rich, level = 1, method = "h", 
fitfunction = rma.glmulti, crit = "aicc", confsetsize = candidateset1)


print(richnessmodels1)
bestrichnessmodel1<-summary(richnessmodels1@objects[[1]])
plot(richnessmodels1, type="s")



# toprich <- weightable(richnessmodels)
# toprich <- toprich[toprich$aicc <= min(toprich$aicc) + 2,]
# toprich
# 
# eval(metafor:::.glmulti)
# 
# richweights <- round(coef.glmulti(richnessmodels, select = nrow(toprich), 
#                                   icmethod = "Burnham"), 7)
# richweights <- richweights[, c(1, 4, 5)]
# richweights




bestrich1<-rma.mv(yi, vi, mods = ~ S.disturbance + multapp + Seeded + Burned +  yeartrans + Precip, random = rand.var, data = ROM.rich, 
                  slab = paste(author, year))

rich1<-rma.mv(yi, vi, mods = ~  Seeded * yeartrans, random = rand.var, data = ROM.rich, 
       slab = paste(author, year))

rich2<-rma.mv(yi, vi, mods = ~  S.disturbance * yeartrans, random = rand.var, data = ROM.rich, 
       slab = paste(author, year))

rich3 <-rma.mv(yi, vi, mods = ~  multapp * yeartrans, random = rand.var, data = ROM.rich, 
       slab = paste(author, year))

rich4 <- rma.mv(yi, vi, mods = ~  S.disturbance * Precip, random = rand.var, data = ROM.rich, 
       slab = paste(author, year))

rich5  <- rma.mv(yi, vi, mods = ~  Precip * multapp, random = rand.var, data = ROM.rich, 
                 slab = paste(author, year))

rich6  <- rma.mv(yi, vi, mods = ~   Seeded * Precip, random = rand.var, data = ROM.rich, 
                  slab = paste(author, year))

rich7  <- rma.mv(yi, vi, mods = ~   Burned * Precip, random = rand.var, data = ROM.rich, 
                 slab = paste(author, year))

rich8  <- rma.mv(yi, vi, mods = ~   Burned * yeartrans, random = rand.var, data = ROM.rich, 
                 slab = paste(author, year))

# bubble.seed=ggplot(preds.dist, aes(x=S.Disturbance,y=pred,color=S.Disturbance))+
#   geom_jitter(data=ROM.rich, aes(x=S.disturbance, y = yi, size=1/vi, color=S.disturbance, alpha=.5))+
#   geom_boxplot()+
#   geom_hline(yintercept=0, linetype="dashed", size=.5)+
#   scale_color_manual(values=c("black","red"))+
#   ylab("Log Response Ratio")+
#   xlab("Disturbance")+
#   theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=10),
#         axis.title=element_text(size=10), legend.position = c(0.1,1),
#         legend.justification = c(0,1), legend.title = element_blank())+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
#   ylim(-2,2.5) + scale_alpha(guide=FALSE) + scale_size(guide=FALSE) +
#   scale_fill_discrete(labels = c("Not Disturbed", "Disturbed"))

# ggplot(preds.dist, aes(x=S.Disturbance,y=pred,color=S.Disturbance))+
#   geom_jitter(data=ROM.rich, aes(x=S.disturbance, y = yi, size=1/vi, color=S.disturbance, alpha=.5))+
#   geom_boxplot()+
#   geom_hline(yintercept=0, linetype="dashed", size=.5)+
#   scale_color_manual(values=c("black","red"))+
#   ylab("Log Response Ratio")+
#   xlab("Disturbance")+
#   theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=10),
#         axis.title=element_text(size=10), legend.position = "none",
#         legend.justification = c(0,1))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
# #   ylim(-2,2.5) + scale_alpha(guide=FALSE) + scale_size(guide=FALSE)
# 
# ggplot(preds.dist, aes(x=Seeded,y=pred,color=Seeded))+
#   geom_jitter(data=ROM.rich, aes(x=Seeded, y = yi, size=1/vi, color=Seeded, alpha=.5))+
#   geom_boxplot()+
#   geom_hline(yintercept=0, linetype="dashed", size=.5)+
#   scale_color_manual(values=c("black","red"))+
#   ylab("Log Response Ratio")+
#   xlab("Disturbance")+
#   theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=10),
#         axis.title=element_text(size=10), legend.position = "none",
#         legend.justification = c(0,1))+
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))+
#   ylim(-2,2.5) + scale_alpha(guide=FALSE) + scale_size(guide=FALSE)
#plots with interactions

richnessmodelslist = list(bestrich1, rich1, rich2, rich3, rich4, rich5, rich6, rich7,
                          rich8)


# model fit diagnostics: pseudo R2, marginal/conditional R2, AICc, and I2
pseudo.r2r <- c()
marg.r2r <- c()
cond.r2r <- c()
aiccr <- c()
I2r <- c()
for (i in 1:length(richnessmodelslist)) {
  # pseudo R2 of proportional reduction of variance explained when fitting reduced model relative to model with random effects only
  pseudo.r2r[i] <- (sum(rma.random.r$sigma2) - sum(richnessmodelslist[[i]]$sigma2)) / sum(rma.random.r$sigma2)
  
  # marginal and conditional R2 for linear mixed-models1 from Nakagawa and Scheilzeth 2013
  # variance from fixed effects
  fe.total <- c()
  for (j in 2:ncol(model.matrix(richnessmodelslist[[i]]))) {
    fe <- richnessmodelslist[[i]]$b[j] * model.matrix(richnessmodelslist[[i]])[, j] 
    fe.total <- c(fe.total, fe)
  }
  v.fix <- var(fe.total)
  # variance from random effects
  v.rand <- sum(richnessmodelslist[[i]]$sigma2)
  # variance from residuals
  v.resid <- var(residuals(richnessmodelslist[[i]]))
  # marginal R2 - total variance explained from fixed effects
  marg.r2r[i] <- v.fix / (v.fix + v.rand + v.resid)
  # conditional R2 - total variance explained from fixed and random effects
  cond.r2r[i] <- (v.fix + v.rand) / (v.fix + v.rand + v.resid)
  
  # from Wolfgang Viechtbauer (http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate)
  # I2 statistic - amount of heterogeneity relative to the total amount of variance in observed effects (how much heterogeneity contributes to total variance)
  W <- diag(1/ROM.rich$vi)
  X <- model.matrix(richnessmodelslist[[i]])
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2r[i] <- 100 * sum(richnessmodelslist[[i]]$sigma2) / (sum(richnessmodelslist[[i]]$sigma2) + (richnessmodelslist[[i]]$k - 
                                                                                                richnessmodelslist[[i]]$p) / sum(diag(P)))
  
  aiccr[i] <- fitstats(richnessmodelslist[[i]])[5]
}

models.richness<-cbind(pseudo.r2r,marg.r2r, cond.r2r, aiccr,I2r)
models.richness
#Variable diagnostics ---------------------------------------
# hetcor in the polycor allows for the calculation of r between numeric, polyserial correlations between numeric and ordinal, polychoric correlations between ordinal
# correlations <- hetcor(data = dat, std.err = FALSE)


subset.richness <- subset(ROM.rich, select = c(precip.log, S.disturbance,
                                               Seeded, multapp, yi, vi))
#subset.slope <- subset.richness[-which(is.na(subset.richness$vi),)]
subset.richness$ll <- subset.richness$yi - subset.richness$vi
subset.richness$ul <- subset.richness$yi + subset.richness$vi
subset.richness$overlap <- ifelse((subset.richness$yi + subset.richness$vi) *
                                    (subset.richness$yi - subset.richness$vi)
                                    > 0, 1, 0)
subset.richness$pos <- ifelse(subset.richness$overlap > 0 & 
                                subset.richness$vi > 0, 1, 0)
subset.richness$neg <- ifelse(subset.richness$overlap > 0 & 
                                subset.richness$vi < 0, 2, 0)
subset.richness$sign <- as.factor(subset.richness$pos + 
                                    subset.richness$neg)
subset.richness$distances <- ifelse(subset.richness$yi > 0, 
                                    abs(subset.richness$ll), 
                                    abs(subset.richness$ul))
subset.richness$distances[which(subset.richness$sign == 0)] <- 0
richness.se <- arrange(subset.richness, distances)
richness.se$index <- c(1:nrow(richness.se))

richness_forest <- ggplot(data = richness.se, aes(x = yi, y = index)) + 
  geom_point() + geom_errorbarh(aes(xmax = ul, xmin = ll)) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray65", size = 0.2) + 
  ylab("") + xlab("Log Response Ratio") + coord_cartesian(ylim = c(0, 120)) + 
  theme_classic() + theme(text = element_text(size = 12))
print(richness_forest)

richnesstable <- data.frame(variate = c(rep(x = "Biosolid level", length(rma.lev.r$b)),
                                rep(x = "Mixture", length(rma.mix.r$b)), 
                                rep(x = "Disturbance", length(rma.dis.r$b)), 
                                rep(x = "Multiple application", length(rma.mult.r$b)),
                                rep(x = "Seeded", length(rma.seed.r$b)), 
                                rep(x = "Temperature", length(rma.temp.r$b)), 
                                rep(x = "Precipitation", length(rma.precip.r$b)),
                                rep(x = "Biome", length(rma.biome.r$b)),
                                rep(x = "Year", length(rma.yr.r$b))),
                    LOR = c(rma.lev.r$b, rma.mix.r$b, rma.dis.r$b, rma.mult.r$b, 
                            rma.seed.r$b, rma.temp.r$b, rma.precip.r$b, rma.biome.r$b, rma.yr.r$b),
                    CI.LB = c(rma.lev.r$ci.lb, rma.mix.r$ci.lb, rma.dis.r$ci.lb, rma.mult.r$ci.lb,
                              rma.seed.r$ci.lb, rma.temp.r$ci.lb, rma.precip.r$ci.lb,
                               rma.biome.r$ci.lb, rma.yr.r$ci.lb), 
                    CI.UB = c(rma.lev.r$ci.ub, rma.mix.r$ci.ub, rma.dis.r$ci.ub, rma.mult.r$ci.ub,
                             rma.seed.r$ci.ub, rma.temp.r$ci.ub, rma.precip.r$ci.ub,
                              rma.biome.r$ci.ub, rma.yr.r$ci.ub),
                    lab = c("Biosolid Level_int", "Biosolid Level_slope", levels(ROM.rich$Mixture), 
                            levels(ROM.rich$S.disturbance), levels(ROM.rich$multapp), 
                            levels(ROM.rich$Seeded), "Temperature_int", "Temperature_slope",
                            "Precipitation_int", "Precipitation_slope", levels(ROM.rich$biome), "Year_int", "Year_slope"),
                    N = c(rep(sum(!is.na(ROM.rich$Biosolid.level)),2), 
                          summary(ROM.rich$Mixture[!is.na(ROM.rich$Mixture)]), 
                          summary(ROM.rich$S.disturbance[!is.na(ROM.rich$S.disturbance)]),
                          summary(ROM.rich$multapp[!is.na(ROM.rich$multapp)]), 
                          summary(ROM.rich$Seeded[!is.na(ROM.rich$Seeded)]),
                          rep(sum(!is.na(ROM.rich$temp.log)),2),
                          rep(sum(!is.na(ROM.rich$precip.log)),2),
                          summary(ROM.rich$biome[!is.na(ROM.rich$biome)]),
                          rep(sum(!is.na(ROM.rich$year.log)),2)))



biome.tab <- subset(table, variate == "Biome")
biome.tab$lab <- as.character(biome.tab$lab)
biome.levels <- biome.tab[order(biome.tab$LOR), 5]
mixture.tab <- subset(table, variate == "Mixture")
mixture.tab$lab <- as.character(mixture.tab$lab)
mixture.levels <- mixture.tab[order(mixture.tab$LOR), 5]

biome.plot <- ggplot(data = subset(richnesstable, variate == "Biome"), 
                     aes(x = factor(lab, levels = biome.levels), y = LOR)) + 
  geom_point() + geom_errorbar(aes(ymax = CI.UB, ymin = CI.LB, 
                                   width = 0.1)) + geom_text(aes(y = -1.2, 
                                                                 label = paste("n = ", N, sep = "")), size = 3) + 
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray65", size = 0.2) + 
  annotate("text", x = 0.7, y = 0.5, label = "b", 
           fontface = "bold") + scale_x_discrete(name = "", 
                                                 labels = c("Temperate seasonal forest", "Shrubland", "Temperate grassland")) + 
  ylab(label = expression(paste("Effect size (", italic("LRR"),")"))) + 
  coord_cartesian(ylim = c(-1.2, 0.5)) + 
  theme_classic() + theme(text = element_text(size = 14)) 
png(file = "biome.jpg")
biome.plot
dev.off()

mult.tab <- subset(richnesstable, variate == "Multiple application")
mult.tab$lab <- as.character(mult.tab$lab)
mult.levels <- mult.tab[order(mult.tab$LOR), 5]
mult.plot <- ggplot(data = subset(table, variate == "Multiple application"), 
                    aes(x = factor(lab, levels = mult.levels), y = LOR)) + 
  geom_point() + geom_errorbar(aes(ymax = CI.UB, ymin = CI.LB, width = 0.1)) + 
  geom_text(aes(y = -0.8, label = paste("n = ", N, sep = "")), size = 3) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray65", size = 0.2) + 
  annotate("text", x = 0.7, y = 0.3, label = "b", fontface = "bold") + 
  scale_x_discrete(name = "", labels = c("Yes", "No")) + 
  ylab(label = expression(paste("Effect size (", italic("LRR"),")"))) + 
  coord_cartesian(ylim = c(-0.8, 0.3)) + theme_classic() + 
  theme(text = element_text(size = 14)) 
png(file = "mult.jpg")
mult.plot
dev.off()

dist.tab <- subset(richnesstable, variate == "Disturbance")
dist.tab$lab <- as.character(dist.tab$lab)
dist.levels <- dist.tab[order(dist.tab$LOR), 5]
distr.plot <- ggplot(data = subset(richnesstable, variate == "Disturbance"), 
                    aes(x = factor(lab, levels = dist.levels), y = LOR)) + 
  geom_point() + geom_errorbar(aes(ymax = CI.UB, 
                                   ymin = CI.LB, width = 0.1)) + geom_text(aes(y = -0.6, 
                                                                               label = paste("n = ", N, sep = "")), size = 3) + 
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "gray65", size = 0.2) + annotate("text",
                                                      x = 0.7, y = 0.5, label = "b", fontface = "bold") + 
  scale_x_discrete(name = "", labels = c("Severe disturbance absent", "Severe disturbance present")) + 
  ylab(label = expression(paste("Effect size (", italic("LRR"),")"))) + 
  coord_cartesian(ylim = c(-0.6, 0.2)) + 
  theme_classic() + theme(text = element_text(size = 14)) 
tiff(file = "distr.tiff", width = 16, height= 12, units ='cm', res=600)
distr.plot
dev.off()

seed.tab <- subset(richnesstable, variate == "Seeded")
seed.tab$lab <- as.character(seed.tab$lab)
seed.levels <- seed.tab[order(seed.tab$LOR), 5]
seedr.plot <- ggplot(data = subset(richnesstable, variate == "Seeded"), 
                    aes(x = factor(lab, levels = seed.levels), y = LOR)) + 
  geom_point() + geom_errorbar(aes(ymax = CI.UB, ymin = CI.LB, 
                                   width = 0.1)) + geom_text(aes(y = -1.7, label = paste("n = ", N, sep = "")), size = 5) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray65", size = 0.2) + 
  annotate("text", x = 0.7, y = 0.5, label = "", fontface = "bold") + 
  scale_x_discrete(name = "", labels = c("Seeded", "Not Seeded")) + 
  ylab(label = expression(paste("Log Response Ratio"))) + 
  coord_cartesian(ylim = c(-1.7, 1.7)) + 
  theme_classic() + theme(text = element_text(size = 20)) 
tiff(file = "seedr.tiff", width = 16, height= 12, units ='cm', res=600)
seedr.plot
dev.off()

mixture.tab <- subset(richnesstable, variate == "Mixture")
mixture.tab$lab <- as.character(mixture.tab$lab)
mixture.levels <- mixture.tab[order(mixture.tab$LOR), 5]
mixr.plot <- ggplot(data = subset(richnesstable, variate == "Mixture"), 
                     aes(x = factor(lab, levels = mixture.levels), y = LOR)) + 
  geom_point() + geom_errorbar(aes(ymax = CI.UB, ymin = CI.LB, 
                                   width = 0.1)) + geom_text(aes(y = -0.9, label = paste("n = ", N, sep = "")), size = 3) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray65", size = 0.2) + 
  annotate("text", x = 0.7, y = 0.5, label = "b", fontface = "bold") + 
  scale_x_discrete(name = "", labels = c("Mixed", "Not Mixed")) + 
  ylab(label = expression(paste("Effect size (", italic("LRR"),")"))) + 
  coord_cartesian(ylim = c(-0.9, 0.3)) + 
  theme_classic() + theme(text = element_text(size = 14)) 
tiff(file = "mixr.tiff", width = 16, height= 12, units ='cm', res=600)
mixr.plot
dev.off()

precip.levels <- print.list.rma(predict.rma(rma.precip.r))
precip.levels <- data.frame(precip.levels, 
                          precip = ROM.rich$Precip[which(!is.na(ROM.rich$Precip), 
                                                       arr.ind = TRUE)])
precipr.plot <- ggplot(precip.levels, aes(x = precip, y = pred)) + geom_line() +
  geom_line(aes(x = precip, y = ci.lb), linetype = "dashed") + 
  geom_point(data = ROM.rich, aes(x = Precip, y = yi, size=1/vi)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray65", size = 0.2) +
  geom_line(aes(x = precip, y = ci.ub), linetype = "dashed") + 
  geom_rug(sides = "b") + 
  ylab(label = expression(paste("Log Response Ratio"))) + 
  xlab(expression(paste("Precipitation [mm]"))) +
  scale_x_continuous(breaks = seq(35, 105, 30), limits = c(35, 105)) + 
  ylim(c(-1.7, 1.7)) + theme_classic() + 
  theme(text = element_text(size = 20))
tiff(file = "precipr.tiff", width = 16, height= 12, units ='cm', res=600)
precipr.plot
dev.off()


seedprecipr<-plot_grid(seedr.plot, precipr.plot,
                     labels = c("a", "b"),
                     ncol = 2, nrow = 1)

tiff("seedprecipr.tiff", width = 28, height= 12, units ='cm', res=600)
seedprecipr
dev.off()

tab.levels <- print.list.rma(predict.rma(rma.lev))
tab.levels <- data.frame(tab.levels, 
                         bio.levels = ROM.rich$level.log[which(!is.na(ROM.rich$level.log), 
                                                               arr.ind = TRUE)])
levels.plot <- ggplot(tab.levels, aes(x = bio.levels, y = pred)) + geom_line() +
  geom_line(aes(x = bio.levels, y = ci.lb), linetype = "dashed") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray65", size = 0.2) +
  geom_line(aes(x = bio.levels, y = ci.ub), linetype = "dashed") + 
  geom_rug(sides = "b") + 
  ylab(label = expression(paste("Effect size (", italic("LRR"),")"))) + 
  xlab("Biosolid Level (log transformed") + 
  scale_x_continuous(breaks = seq(0, 6, 0.5), limits = c(0, 6)) + 
  ylim(c(-0.6, 0.6)) + theme_classic() + 
  theme(text = element_text(size = 14))

year.levels <- print.list.rma(predict.rma(rma.yr))
year.levels <- data.frame(year.levels, 
                          year.levels = ROM.rich$year.log[which(!is.na(ROM.rich$year.log), 
                                                                arr.ind = TRUE)])
year.plot <- ggplot(year.levels, aes(x = year.levels, y = pred)) + geom_line() +
  geom_line(aes(x = year.levels, y = ci.lb), linetype = "dashed") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray65", size = 0.2) +
  geom_line(aes(x = year.levels, y = ci.ub), linetype = "dashed") + 
  geom_rug(sides = "b") + 
  ylab(label = expression(paste("Effect size (", italic("LRR"),")"))) + 
  xlab("Years (log transformed") + 
  scale_x_continuous(breaks = seq(0, 3, 0.5), limits = c(0, 3)) + 
  ylim(c(-0.5, 0.25)) + theme_classic() + 
  theme(text = element_text(size = 14))

temp.levels <- print.list.rma(predict.rma(rma.temp.r))
temp.levels <- data.frame(temp.levels, 
                          temp.levels = ROM.rich$Temp[which(!is.na(ROM.rich$Temp), 
                                                            arr.ind = TRUE)])
temp.plot <- ggplot(temp.levels, aes(x = temp.levels, y = pred)) + geom_line() +
  geom_line(aes(x = temp.levels, y = ci.lb), linetype = "dashed") + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray65",
             size = 0.2) + geom_line(aes(x = temp.levels, y = ci.ub), 
                                     linetype = "dashed") + geom_rug(sides = "b") + 
  ylab(label = expression(paste("Effect size (", italic("LRR"),")"))) + 
  xlab("Temperature") + scale_x_continuous(breaks = seq(-1, 18, 1), 
                                           limits = c(-1, 18)) + ylim(c(-0.6, 0.7)) + theme_classic() +
  theme(text = element_text(size = 14))

richf<-ROM.rich
richf$author.year<-paste(ROM.rich$author, richf$year, sep = "-")


##### HEAD:Datasheets for analysis/biomass.R
forest.rich<-viz_forest(x = richf[c("pred", "se")], variant="classic",study_labels= richf$author.year,summary_label ="Summary Effect", xlab = "Log Response Ratio")
tiff("rich.forest.tiff", width = 16, height= 25, units ='cm', res=600)
forest.rich
dev.off()

seedprecipr<-plot_grid(seedr.plot, precipr.plot,
                     labels = c("a", "b"),
                     ncol = 2, nrow = 1)
tiff("seedprecipr.tiff", width = 16, height= 25, units ='cm', res=600)
seedprecipr
dev.off()
