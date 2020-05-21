library(readxl)
library(xlsx)

overall <- read_excel("~/Documents/GitHub/Grassland restoration with biosolids/Grassland-restoration-with-biosolids/overall_effect_size.xlsx")%>%
  janitor::clean_names()

head(overall)

overall.plot=ggplot(overall, aes(x=variable,y=lrr))+
  geom_point( fill="black", size=3) +
  geom_errorbar(aes(ymax=lower, ymin=upper),width=.25 )+
  geom_hline(yintercept=0, linetype="dashed", size=.5)+
  ylab("Log Response Ratio")+
  scale_x_discrete(name ="", labels=c("Exotic species abundance (%)","Shannon diversity","Species richness",
                                    "Total vegetative cover" , "ANPP"),
        breaks=c("exotics","shannon diversity","richness","cover","productivity"),
        limits=c("exotics","shannon diversity","richness","cover","productivity"))  +
  coord_flip()+
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=10),
        axis.title=element_text(size=10), legend.position = "none")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  ylim(-1,2)
overall.plot
