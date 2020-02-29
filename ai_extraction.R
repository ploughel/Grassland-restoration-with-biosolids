library(rgdal)
library(raster)
library(ggplot2)
library(sp)
library("ggspatial")

r=raster("~/Desktop/ai_et0.tif")
#plot(r)

xy<-read.csv("biomass.coord.csv", header=T)

points <- SpatialPoints(xy, proj4string = r@crs)
values <- extract(r,points) 
plot(points)


df <- cbind.data.frame(coordinates(points),values)

df.xy<-cbind(extract(r, xy, df=T), xy)
  # write.csv(df.xy,file="ai_biomass2.csv")
table(is.na(df.xy$ai_et0))

#cover

xy.cover<-read.csv("cover.coord.csv", header = T)


points.c <- SpatialPoints(xy.cover, proj4string = r@crs)

values.c <- extract(r,points.c) #buffer is in m at 1000km at 100km still have 150-T and 91-F

plot(points.c)


df <- cbind.data.frame(coordinates(points.c),values.c)

df.xy.c<-cbind(extract(r, xy.cover, df=T), xy.cover)


# write.csv(df.xy.c,file="ai_cover.csv")

#richness

library(readxl)



Richness <- read_excel("~/Documents/GitHub/Grassland restoration with biosolids/Grassland-restoration-with-biosolids/data_sheets.xlsx", 
                    sheet = "Richness")

names(Richness)

lat<-Richness$Latitude
long<-Richness$Longitude

xy.richness<-cbind(long,lat)

points.r <- SpatialPoints(xy.richness, proj4string = r@crs)

values.r <- extract(r,points.r) 

plot(points.r)


df <- cbind.data.frame(coordinates(points.r),values.r)

df.xy.r<-cbind(extract(r, xy.richness, df=T), xy.richness)


 #write.csv(df.xy.r,file="ai_richness.csv")

#exotics

xy.exotics<-read.csv("exotics.coord.csv", header = T)


points.e <- SpatialPoints(xy.exotics, proj4string = r@crs)

values.e <- extract(r,points.e) 

plot(points.e)


df <- cbind.data.frame(coordinates(points.e),values.e)

df.xy.e<-cbind(extract(r, xy.exotics, df=T), xy.exotics)


# write.csv(df.xy.e,file="ai_exotics.csv")

df<-rbind(df.xy,df.xy.c,df.xy.e,df.xy.r)
pal <- colorRampPalette(c("firebrick", "yellow", "forestgreen", "cyan","royalblue4", "purple4"))

plot(r, col=pal(50))
points(xy.exotics, pch=19)
