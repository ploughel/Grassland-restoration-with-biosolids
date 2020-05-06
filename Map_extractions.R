#Map and climate zones
#install.packages("cowplot")
library(dplyr)
library(ggplot2)
library(rgdal)
library(raster)
library(ggplot2)
library(sp)
library("ggspatial")
library(readxl)
library(maps)
library("kgc")


#upload data


Richness <- read_excel("~/Documents/GitHub/Grassland restoration with biosolids/Grassland-restoration-with-biosolids/data_sheets.xlsx", 
                    sheet = "Richness")
Exotics <- read_excel("~/Documents/GitHub/Grassland restoration with biosolids/Grassland-restoration-with-biosolids/data_sheets.xlsx", 
                       sheet = "Exotics")
Cover <- read_excel("~/Documents/GitHub/Grassland restoration with biosolids/Grassland-restoration-with-biosolids/data_sheets.xlsx", 
                       sheet = "Cover")

Biomass <- read_excel("~/Documents/GitHub/Grassland restoration with biosolids/Grassland-restoration-with-biosolids/data_sheets.xlsx", 
                       sheet = "Biomass")



#

world<-map_data('world')


names
Biomass$lat<-Biomass$`Coordinates (Latitude)`
Biomass$long<-Biomass$`Coordinate (longitude)`


Cover$lat<-Cover$`Coordinates (Latitude)`
Cover$long<-Cover$`Coordinate (longitude)`

Richness$lat<-Richness$Latitude
Richness$long<-Richness$Longitude

Exotics$lat<-Exotics$`Coordinates (Latitude)`
Exotics$long<-Exotics$`Coordinate (longitude)`



list<-rep("Productivity",length(Biomass$author))
list.p.lat.long<-Biomass[,29:30]
papers<-Biomass[,1:4]

BH.map <- cbind(list, list.p.lat.long,papers)

names(Cover)
list<-rep("Cover",length(Cover$author))
list.c.lat.long<-Cover[,30:31]
papers.c<-Cover[,1:4]

cover.map <- cbind(list,
                   list.c.lat.long,papers.c  )

summary(cover.map)
list<-rep("Richness",length(Richness$author))
list.r.lat.long<-Richness[,29:30]
papers<-Richness[,1:4]

rich.map <- cbind(list, list.r.lat.long,papers)

dim(Exotics)
list<-rep("Exotics",length(Exotics$author))
list.e.lat.long<-Exotics[,28:29]
papers<-Exotics[,1:4]

e.map <- cbind(list, list.e.lat.long,papers)

names(cover.map)

Map<-rbind(BH.map,rich.map,e.map, cover.map)

names(Map)
p <- ggplot(legend=FALSE) +
  geom_polygon( data=world, aes(x=long, y=lat,group=group), fill="grey") +
  #theme(panel.background = theme_blank()) +
  #theme(panel.grid.major = theme_blank()) +
  #theme(panel.grid.minor = theme_blank()) +
  #theme(axis.text.x = theme_blank(),axis.text.y = theme_blank()) +
  #theme(axis.ticks = theme_blank()) +
  xlab("Longitude") + ylab("Latitude")+
  geom_point(data=Map,aes(x=long,y=lat),colour="black")+
  ggtitle("")+
  ylab("Latitude")+
  xlab("Longitude")+
  theme(axis.ticks.length=unit(.5, "cm"), axis.text=element_text(size=20),
        axis.title=element_text(size=20))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
p

r <- getData("worldclim",var="bio",res=10)
names(r)
r <- r[[c(1,12)]]
names(r) <- c("Temp","Precip")

lats <- Map$lat
lons <- Map$long

coords <- data.frame(x=lons,y=lats)

points <- SpatialPoints(coords, proj4string = r@crs)
values <- extract(r,points)

df <- cbind.data.frame(coordinates(points),values)
df$Temp<-df$Temp/10
df$Precip<-df$Precip/10



df$author.year<-as.factor(paste(Map$author,Map$year, sep="-"))
df$Variable<-as.factor(Map$list)

# write.csv(df, file ="tem.pre.all.csv")


#Extract ai
r=raster("~/Desktop/ai_et0.tif")
#plot(r)

xy<- cbind(lons,lats)

points <- SpatialPoints(xy, proj4string = r@crs)
values <- extract(r,points) 
plot(points)


df <- cbind.data.frame(coordinates(points),values)

df.xy<-cbind(extract(r, xy, df=T), xy)
 # write.csv(df.xy,file="ai_all.csv")

pal <- colorRampPalette(c("firebrick", "yellow", "forestgreen", "cyan","royalblue4", "purple4"))

plot(r, col=pal(50))
points(xy, pch=19, cex=.65)


