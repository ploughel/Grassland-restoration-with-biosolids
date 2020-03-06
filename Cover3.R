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


Cover <- read_excel("~/Documents/GitHub/Grassland restoration with biosolids/Grassland-restoration-with-biosolids/data_sheets.xlsx", 
                    sheet = "Cover")


ROM.Cover <- escalc(n2i = Nc, n1i = Ne, m2i = Mc+1/2, m1i = Me, 
                    sd2i = scimp, sd1i = seimp, data = Cover, measure = "ROM", 
                    append = TRUE )


HG.Cover <- escalc(n2i = Nc, n1i = Ne, m2i = Mc, m1i = Me, 
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

rand.var=list(~ 1|Paper.ID/Case)

rma.random <- rma.mv(yi = yi, V = vi, random = rand.var, data = ROM.Cover, method = "ML")
summary(rma.random)
plot(residuals(rma.random))

rma.random.hg <- rma.mv(yi = yi, V = vi, random = rand.var, data = HG.Cover, method = "ML")
plot(residuals(rma.random.hg))

#-----------------------# model selection 
# rma.glmulti<- function(formula, data, ...)
#   rma(formula, vi, data=data, method ="ML",...)

# ROM.cover.size <- glmulti(yi ~Biosolid.level..Mg.ha.1.+yeartrans+Temp+Precip+Mixture+Burn+
#                          Seeded+Multiple.application+S.dist+ai, data=ROM.Cover, method="d",
#                       level=1, crit="aicc", fitfunction = rma.glmulti)
# 
# 
# ROM.cover.bf <- glmulti(yi ~ Biosolid.level..Mg.ha.1.+yeartrans+Temp+Precip+Mixture+Burn+
#                         Seeded+Multiple.application+S.dist+ai, data=ROM.Cover, method="h",
#                      level=1, crit="aicc", confsetsize=ROM.cover.size, fitfunction = rma.glmulti)
# 

plot(ROM.cover.bf, type="s")

#best fit model
ROM.cover.bf<-rma.mv(yi, vi, mods = ~Temp+Seeded+Burn+Multiple.application+Precip+S.dist, 
                       random = rand.var,data = ROM.Cover,slab = paste(author, year), method="ML")


#model with interactions
ROM.cover.int.1<-rma.mv(yi, vi, mods = ~Temp*Seeded, 
                      random = rand.var,data = ROM.Cover,slab = paste(author, year), method="ML")


ROM.cover.int.2<-rma.mv(yi, vi, mods = ~Temp*Burn, 
                       random = rand.var,data = ROM.Cover,slab = paste(author, year), method="ML")

ROM.cover.int.3<-rma.mv(yi, vi, mods = ~Temp*Multiple.application, 
                       random = rand.var,data = ROM.Cover,slab = paste(author, year), method="ML")


ROM.cover.int.4<-rma.mv(yi, vi, mods = ~Precip*Seeded, 
                       random = rand.var,data = ROM.Cover,slab = paste(author, year), method="ML")

ROM.cover.int.5<-rma.mv(yi, vi, mods = ~Precip*Burn, 
                       random = rand.var,data = ROM.Cover,slab = paste(author, year), method="ML")

ROM.cover.int.6<-rma.mv(yi, vi, mods = ~Precip*Multiple.application, 
                       random = rand.var,data = ROM.Cover,slab = paste(author, year), method="ML")


mods = list(ROM.cover.bf,ROM.cover.int.1, ROM.cover.int.2,
            ROM.cover.int.3, ROM.cover.int.4,ROM.cover.int.5, ROM.cover.int.6)



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

cbind(pseudo.r2,marg.r2,cond.r2,aicc,I2)

vif(ROM.cover.int.1, table=TRUE)
ROM.cover.bf
ROM.cover.int.1




#=========================================graphs=========



# write.csv(as.data.frame(coef(summary(ROM.cover.bf))), file="cover.mod.csv")
df.cover<-read.csv("cover.mod.csv")
df.cover

df.ai<-df.cover[2:6,]

LRR.ai.cover = ggplot(data=df.ai,
                        aes(x = X, y = estimate, ymin =ci.lb , ymax = ci.ub ))+
  geom_point(aes(color=X))+
  geom_hline(yintercept =0, linetype=2)+
  #scale_color_brewer(type = 'div', palette = 1)+
  xlab('ai')+ ylab("Log Response Ratio")+
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
LRR.ai.cover

# tiff("cover.ai.tiff", width = 20, height= 12, units ='cm', res=600)
# LRR.ai.cover
# dev.off()


cover.temp.mod<-rma.mv(yi, vi, mods = ~Temp, 
                                     random = rand.var,data = ROM.Cover,slab = paste(author, year), method="ML")

preds<-predict(ROM.cover.bf, addx = TRUE)
preds<-as.data.frame(preds)

preds2<-predict(cover.temp.mod, addx = TRUE)
preds2<-as.data.frame(preds2)

head(preds)
head(preds2)
#forest plot

forest.cover<-viz_forest(x =preds[c("pred", "se")], variant="classic",summary_label ="Summary Effect", xlab = "Log Response Ratio")


# tiff("cover.forest.tiff", width = 16, height= 25, units ='cm', res=600)
# forest.cover
# dev.off()

bubble.temp=ggplot(preds, aes(x=X.Temp,y=pred))+
  stat_smooth(data=preds,method="glm",formula=y~x,fullrange=T,se=T,size=1, color="black"     )+ 
  
  geom_hline(yintercept=0, linetype="dashed", size=.5)+
  scale_color_manual(values=c("darkblue","forestgreen"))+
  scale_fill_manual(values=c("darkblue","forestgreen"))+
  geom_point(data=ROM.Cover, aes(x=Temp, y = yi, size=1/vi, alpha=.5))+
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


#seeded

mod1.seed.n<-rma.mv(yi, vi, mods = ~ Temp, random = rand.var,data = ROM.Cover,slab = paste(author, year), subset = Seeded=="N", method = "ML")
mod1.seed<-rma.mv(yi, vi, mods = ~ Temp, random = rand.var,data = ROM.Cover,slab = paste(author, year), subset = Seeded=="Y", method = "ML")


preds.seed<-predict(mod1.seed, addx = TRUE)
preds.s<-as.data.frame(preds.seed)

preds.s$Seeded<-rep("Y", times = 87)

preds.ns<-predict(mod1.seed.n, addx = TRUE)
preds.ns<-as.data.frame(preds.ns)
preds.ns$Seeded<-rep("N", times = 127)


preds.seeded=rbind(preds.s,preds.ns)

names(preds.seeded)


bubble.seed=ggplot(preds.seeded, aes(x=Seeded,y=pred,color=Seeded))+
geom_boxplot()+
geom_hline(yintercept=0, linetype="dashed", size=.5)+
  scale_color_manual(values=c("darkblue","forestgreen"))+
  
  #geom_point(data=ROM.Cover, aes(x=Seeded, y = yi, size=1/vi))+
  ylab("")+
  xlab("Seeded")+
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=10),
        axis.title=element_text(size=10), legend.position = "bottom")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylim(-2,4.1)

bubble.seed




#burned
mod1.burn.n<-rma.mv(yi, vi, mods = ~ Temp, random = rand.var,data = ROM.Cover,slab = paste(author, year), subset = Burn=="N", method = "ML")
mod1.burn<-rma.mv(yi, vi, mods = ~ Temp, random = rand.var,data = ROM.Cover,slab = paste(author, year), subset = Burn=="Y", method = "ML")


preds.burn<-predict(mod1.burn, addx = TRUE)
preds.b<-as.data.frame(preds.burn)
dim(preds.b)
preds.b$burned<-rep("Y", times = 50)

preds.nb<-predict(mod1.burn.n, addx = TRUE)

preds.nb<-as.data.frame(preds.nb)
preds.nb$burned<-rep("N", times = 164)
names(preds.nb)
names(preds.b)
preds.burned=rbind(preds.b,preds.nb)

names(ROM.Cover)


bubble.burn=ggplot(preds.burned, aes(x=burned,y=pred, color=burned))+
  geom_boxplot()+
  geom_jitter(data=ROM.Cover, aes(x=Burn, y = yi, size=1/vi, color=Burn,alpha=.5))+
  geom_hline(yintercept=0, linetype="dashed", size=.5)+
  scale_color_manual(values=c("black", "red"))+
  #geom_point(data=ROM.Cover, aes(x=Burn, y = yi, size=1/vi))+
  ylab("Log Response Ratio")+
  xlab("Burned")+
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=10),
        axis.title=element_text(size=10), legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylim(-2,5)

bubble.burn

#   tiff("cover.burn.tiff", width = 16, height= 12, units ='cm', res=600)
# bubble.burn
#  dev.off()


#plots with interactions

head(ROM.cover)
prod.int2.Seeded<-rma.mv(yi, vi, mods = ~ Temp, random = rand.var,data = ROM.Cover,slab = paste(author, year), subset = Seeded=="Y", method = "ML")
preds.s<-predict(prod.int2.Seeded, addx = TRUE)
preds.s<-as.data.frame(preds.s)
dim(preds.s)

preds.s$Seeded<-rep("Y", Temps = 87)

prod.int2.Seeded.n<-rma.mv(yi, vi, mods = ~ Temp, random = rand.var,data = ROM.Cover,slab = paste(author, year), subset = Seeded=="N", method = "ML")
preds.s.n<-predict(prod.int2.Seeded.n, addx = TRUE)
preds.s.n<-as.data.frame(preds.s.n)
dim(preds.s.n)

preds.s.n$Seeded<-rep("N", Temps = 124)

preds.Seeded=rbind(preds.s,preds.s.n)
names(preds.Seeded)

bubble.Temp.Seeded=ggplot(preds.Seeded, aes(x=X.Temp,y=pred))+
  stat_smooth(data=preds.Seeded,aes(color= Seeded, fill=Seeded),method="glm",formula=y~x,fullrange=F,se=F,size=1     )+ #, linetype=Seeded 
  geom_ribbon( aes(ymin = preds.Seeded$ci.lb, ymax = preds.Seeded$ci.ub, fill = Seeded), alpha = .15) +
  geom_hline(yintercept=0, linetype="dashed", size=.5)+
  geom_point(data=ROM.Cover, aes(x=Temp, y = yi, size=1/vi, alpha=.5,color=Seeded, group=Seeded, fill=Seeded, shape=Seeded))+
  ylab("")+
  xlab(expression('MAT ('*degree*C*')'))+
  scale_color_manual(values=c("darkblue","forestgreen"))+
scale_fill_manual(values=c("darkblue","forestgreen"))+
  # 
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=10),
        axis.title=element_text(size=10), legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylim(-2,5)

bubble.Temp.Seeded

bubble.Temp.Seeded.raw=ggplot(ROM.Cover, aes(x=Temp,y=yi, size=1/vi, shape=Seeded, fill=Seeded, color=Seeded))+
  stat_smooth(aes(color= Seeded, fill=Seeded),method="glm",fullrange=F,se=T,size=1)+ 
  #geom_ribbon( aes(ymin = preds.Seeded$ci.lb, ymax = preds.Seeded$ci.ub, fill = Seeded), alpha = .15) +
  geom_hline(yintercept=0, linetype="dashed", size=.5)+
  geom_point()+
  ylab("")+
  xlab(expression('MAT ('*degree*C*')'))+
  scale_color_manual(values=c("darkblue","forestgreen"))+
  scale_fill_manual(values=c("darkblue","forestgreen"))+
  # 
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=10),
        axis.title=element_text(size=10), legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylim(-2,5)
 
bubble.Temp.Seeded.raw

library(cowplot)
library(ggpubr)


cover.plot<-plot_grid(bubble.burn,bubble.temp, bubble.Temp.Seeded, 
                      labels = c("a", "b", "c"),
                      ncol = 3, nrow = 1)

cover.plot

#  tiff("cover.tiff", width = 25, height= 8, units ='cm', res=600)
# cover.plot
#  dev.off()
 