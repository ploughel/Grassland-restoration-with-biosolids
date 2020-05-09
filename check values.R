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
library(janitor)


library(readxl)

Biomass <- read_excel("~/Documents/GitHub/Grassland restoration with biosolids/Grassland-restoration-with-biosolids/data_sheetscorr.xlsx", 
                      sheet = "Biomass")%>%
  janitor::clean_names()
  

Cover <- read_excel("~/Documents/GitHub/Grassland restoration with biosolids/Grassland-restoration-with-biosolids/data_sheetscorr.xlsx", 
                      sheet = "Cover")%>%
  janitor::clean_names()

Richness <- read_excel("~/Documents/GitHub/Grassland restoration with biosolids/Grassland-restoration-with-biosolids/data_sheetscorr.xlsx", 
                      sheet = "Richness")%>%
  janitor::clean_names()

Diversity <- read_excel("~/Documents/GitHub/Grassland restoration with biosolids/Grassland-restoration-with-biosolids/data_sheetscorr.xlsx", 
                      sheet = "Diversity(H)")%>%
  janitor::clean_names()

names(Biomass)
names(Richness)
names(Cover)
names(Diversity)

colnames(Cover)[16] <- "time"
colnames(Cover)[26] <- "latitude"
colnames(Cover)[27] <- "longitude"




colnames(Biomass)[19] <- "mixture"
colnames(Biomass)[20] <- "s_dist"
colnames(Biomass)[21] <- "burn"
colnames(Biomass)[22] <- "seeded"
colnames(Biomass)[23] <- "multiple_application"
colnames(Biomass)[26] <- "latitude"
colnames(Biomass)[27] <- "longitude"


colnames(Richness)[16] <- "time"
colnames(Richness)[22] <- "s_dist"
colnames(Richness)[23] <- "burn"
colnames(Richness)[24] <- "multiple_application"
colnames(Richness)[5] <- "biosolid_level_mg_ha_1"


colnames(Diversity)[5] <- "biosolid_level_mg_ha_1"
colnames(Diversity)[16] <- "time"
colnames(Diversity)[22] <- "s_dist"
colnames(Diversity)[23] <- "burn"
colnames(Diversity)[24] <- "multiple_application"





Biomass_meta<-Biomass%>%
  select(paper_id, year,case, author, biosolid_level_mg_ha_1, time, temp, precip, 
                mixture, s_dist, burn, seeded, 
         multiple_application, country, latitude,longitude, ai)



Cover_meta<-Cover%>%
  select(paper_id, year,case, author, biosolid_level_mg_ha_1, time, temp, precip, 
                mixture, s_dist, burn, seeded, 
                multiple_application, country, latitude,longitude, ai)

Richness_meta<-Richness%>%
  dplyr::select(paper_id, year,case, author, biosolid_level_mg_ha_1, time, temp, precip, 
                mixture, s_dist, burn, seeded, 
                multiple_application, country, latitude,longitude, ai)

Diversity_meta<-Diversity%>%
  select(paper_id, case, author,year, biosolid_level_mg_ha_1, time, temp, precip, 
         mixture, s_dist, burn, seeded, 
         multiple_application, country, latitude,longitude, ai)




b.c<-Biomass_meta %>%
  filter(author %in% Cover_meta$author)

b.r<-Biomass_meta %>%
  filter(author %in% Richness_meta$author)

b.d<-Biomass_meta %>%
  filter(author %in% Diversity_meta$author)

c.b<-Cover_meta %>%
  filter(author %in% Biomass_meta$author)

c.r<-Cover_meta %>%
  filter(author %in% Richness_meta$author)

c.d<-Cover_meta %>%
  filter(author %in% Diversity_meta$author)

r.b<-Richness_meta %>%
  filter(author %in% Biomass_meta$author)

r.c<-Richness_meta %>%
  filter(author %in% Cover_meta$author)

r.d<-Richness_meta %>%
  filter(author %in% Diversity_meta$author)

d.b<-Diversity_meta %>%
  filter(author %in% Biomass_meta$author)

d.r<-Diversity_meta %>%
  filter(author %in% Richness_meta$author)

d.c<-Diversity_meta %>%
  filter(author %in% Cover_meta$author)



library("compareDF")
library("htmlTable")
biomass.cover<-compare_df(b.c, c.b, c("paper_id", "year", "biosolid_level_mg_ha_1","time","latitude","longitude"))
create_output_table(biomass.cover, limit=500)

biomass.rich<-compare_df(b.r, r.b, c("paper_id", "year", "biosolid_level_mg_ha_1","time","latitude","longitude"))
create_output_table(biomass.rich, limit=500)

biomass.div<-compare_df(d.b, b.d, c("paper_id", "year","case", "biosolid_level_mg_ha_1","time","latitude","longitude"))
create_output_table(biomass.rich, limit=500)

cover.div<-compare_df(c.d, d.c, c("paper_id", "year","case", "biosolid_level_mg_ha_1","time","latitude","longitude"))
create_output_table(cover.div, limit=500)

cover.rich<-compare_df(c.r, r.c, c("paper_id", "year","case", "biosolid_level_mg_ha_1","time","latitude","longitude"))
create_output_table(cover.rich, limit=500)

div.rich<-compare_df(r.d, d.r, c("paper_id", "year","case", "biosolid_level_mg_ha_1","time","latitude","longitude"))
