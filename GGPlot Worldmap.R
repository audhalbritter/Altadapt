###############################################
### DRAW WORLD MAP WITH WORLDCLIM ELEV DATA
###############################################

#### Libraries ####
library(raster); library(ggplot2);
#library(plyr)

# colours
col.pop <- "#E69F00" # orange
col.garden <- "#009E73" # green
col.molec <- "#CC79A7" # pink
col.chamber <- "#0072B2" # blue

#### Coord pop and site ####
# sites
data.1.ID <- DAT.S$StudyID
keep.these.ID <- dat2$StudyID
dat.site <- DAT.S[data.1.ID %in% keep.these.ID,] # match all sites with the sites included in dat2
dat.site <- dat.site[ !duplicated(dat.site$StudyID,fromLast=FALSE), ]
dat.site <- dat.site %>% select(studysite, latitude_site, longitude_site)

# populations
data.1.ID <- DAT.T$studyunit
keep.these.ID <- dat2$studyunit
dat.pop <- DAT.T[data.1.ID %in% keep.these.ID,] # match all pops with the pops included in dat2
dat.pop <- dat.pop %>% select(lat.pop, long.pop)

# import mol genetic coords
mol <- read.csv("CoordMolGen.csv", header=TRUE, sep = ";")
molec <- mol %>% select(-StudyID)


# Create new factor combining Climate chambers and Greenhouses
new.studysite <- mutate(dat.site, studysite2 =  ifelse(studysite=="climate_chamber","Climate chamber", ifelse(studysite == "greenhouse", "Climate chamber", "Garden")))
new.studysite <- rbind(new.studysite, molec)
Studytype <- new.studysite$studysite2

#### Get Wolrdclim elevation data ####
elev <- getData('worldclim', var='alt', res=2.5)

#### WORLDMAP ####
elev.spdf <- as(elev, "SpatialPixelsDataFrame")
elev.df <- as.data.frame(elev.spdf)
# replace all values > 6000 with 6000 (everything above is not interesting)
elev.df$alt[elev.df$alt > 6000] <- 6000
range(elev.df)

# merge pop and site data set together
worldData <- rbind(new.studysite, data.frame(studysite = "Population", latitude_site = dat.pop$lat.pop, longitude_site = dat.pop$long.pop, studysite2 = "Population"))
worldData <- worldData[with(worldData, order(studysite2, decreasing = TRUE)), ]
worldData$studysite2 <- factor(worldData$studysite2)

# plot world map
w.map <- ggplot() + 
  geom_raster(data = elev.df, aes(x=x, y=y, fill = alt)) + 
  coord_equal() +
  labs(x="", y = "", fill = "Elevation\n", colour = "Study type", shape = "Study type", alpha = "Study type") + 
  scale_y_continuous(limits = c(-60, 85)) +
  scale_color_manual(values = c(col.chamber, col.garden, col.molec, col.pop)) +
  scale_shape_manual(values = c(17, 17, 15, 16)) +
  scale_alpha_manual(values = c(0.9, 0.9, 0.9, 0.3)) +
  scale_fill_gradient(low = "grey0", high = "grey100", limits=c(-416,4000)) + 
  geom_point(aes(x=longitude_site, y=latitude_site, colour  = studysite2, shape = studysite2, alpha = studysite2 ), data = worldData, size=2) + 
  annotate("text", x = -190, y = 80, label = "A", size= 4) +
  theme(legend.position = "none")
w.map



#### EUROPE MAP ####
# Subset EU
dat.site.eu <- subset(worldData, (longitude_site >-10 & longitude_site < 25) & (latitude_site > 36 & latitude_site < 66))

e <- extent(-10,25,36,66)
elev.eu <- crop(elev, e)

# To convert your RasterLayer to a data.frame, you need to convert it to
# a SpatialPixelsDataFrame first
elev.eu.spdf <- as(elev.eu, "SpatialPixelsDataFrame")
elev.eu.df <- as.data.frame(elev.eu.spdf)

# plot Europe map
eu.map <- ggplot() +
  geom_raster(data = elev.eu.df, aes(x=x, y=y, fill = alt)) +
  coord_equal() +
  labs(x="", y = "", fill = "Elevation", colour = "Study type", shape = "Study type", alpha = "Study type") +
  theme(legend.position="none") +
  #scale_y_continuous(limits = c(-60, 85)) +
  scale_color_manual(values = c(col.chamber, col.garden, col.molec, col.pop)) +
  scale_shape_manual(values = c(17, 17, 15, 16)) +
  scale_alpha_manual(values = c(0.9, 0.9, 0.9, 0.3)) +
  scale_fill_gradient(low = "grey0", high = "grey100", limits=c(-416,4000)) + 
  geom_point(aes(x=longitude_site, y=latitude_site, colour  = studysite2, shape = studysite2, alpha = studysite2), data = dat.site.eu, size=2) +
  annotate("text", x = -9, y = 65, label = "C", size= 4)
eu.map



#### USA MAP ####
# Subset US
dat.site.us <- subset(worldData, (longitude_site >-130 & longitude_site < -70) & (latitude_site > 18 & latitude_site < 70))

us <- extent(-130,-70,18,70)
elev.us <- crop(elev, us)

# To convert your RasterLayer to a data.frame, you need to convert it to
# a SpatialPixelsDataFrame first
elev.us.spdf <- as(elev.us, "SpatialPixelsDataFrame")
elev.us.df <- as.data.frame(elev.us.spdf)

# plot US map
us.map <- ggplot() +
  geom_raster(data = elev.us.df, aes(x=x, y=y, fill = alt)) +
  coord_equal() +
  labs(x="", y = "", fill = "Elevation", colour = "Study type", shape = "Study type", alpha = "Study type") +
  theme(legend.position="none") +
  scale_color_manual(values = c(col.chamber, col.garden, col.molec, col.pop)) +
  scale_shape_manual(values = c(17, 17, 15, 16)) +
  scale_alpha_manual(values = c(0.9, 0.9, 0.9, 0.3)) +
  scale_fill_gradient(low = "grey0", high = "grey100", limits=c(-416,4000)) + 
  geom_point(aes(x=longitude_site, y=latitude_site, colour  = studysite2, shape = studysite2, alpha = studysite2), data = dat.site.us, size=2) +
  annotate("text", x = -127, y = 67, label = "B", colour = "white", size= 4)
us.map

library(grid)
library(gridExtra)

LegendPlot <- ggplot() +
  geom_raster(data = elev.us.df, aes(x=x, y=y, fill = alt)) +
  coord_equal() +
  labs(x="", y = "", fill = "Elevation", colour = "Study type", shape = "Study type", alpha = "Study type") +
  scale_fill_gradient(low = "grey0", high = "grey100", limits=c(-416,4000))
  
guides(colour = guide_legend(override.aes = list(size = 2)))
legend_t <- cowplot::get_legend(LegendPlot)



pdf ("Fig.1_Worldmap_GGPlot5.pdf", width=10, height=8, pointsize=6, onefile=TRUE, paper="special")
grid.arrange(w.map,legend_t,us.map,eu.map, layout_matrix = rbind(c(1,1,1,1,1,2),c(3,3,3,4,4,4)))
dev.off()



#### GREY WORLD MAP ####
"ggmap", "maptools", "maps"
library(ggmap)
library(maptools)
library(maps)

mp <- NULL
mapWorld <- borders("world", colour="gray50", fill="gray50") # create a layer of borders
mp <- ggplot() + mapWorld

#Now Layer the cities on top
mp <- mp+ geom_point(aes(x=pop.x, y=pop.y) ,color=col.green, pch=16, size=3)  + geom_point(aes(x=site.x, y=site.y) ,color=ifelse(dat.site$studysite=="Garden",col.red, col.yellow), pch=ifelse(dat.site$studysite=="Garden",15,17), size=3) + theme(axis.title.x = element_blank()) + theme(axis.title.y = element_blank())
mp