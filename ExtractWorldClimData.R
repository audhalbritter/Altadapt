### EXTRACT BIOCLIM VARIABLES FOR POP AND SITE COORDINATES

# READ IN SITE DATA
setwd("/")
site <- read.csv("160210_Altadapt_Study.csv", header=TRUE, fill=TRUE, sep=";")
head(site)
str(site)

# READ IN TRAIT DATA FOR POPS
setwd("/")
pop <- read.csv("160210_Altadapt_Traits.csv", header=TRUE, fill=TRUE, sep=";")
head(pop)
str(pop)


long.site <- site$longitude_site
lat.site <- site$latitude_site
coord <- as.data.frame(cbind(long.site, lat.site))
#coord <- na.omit(coord)
# possible to draw points on the map
points(coord$long.site, coord$lat.site, col="black", pch=1)



setwd("/Volumes/SILVER/Worldclim Data/bio1-9_30s_bil") # all data


# EXTRACT WORLDCLIM VARS
# fun accepts na.rm=T, but then you need to use function(x)...
# Buffer defines a radius around the point in which a value is included or not.
# If df=TRUE, the results are returned as a data frame
# Method: simple returns a sigle value for the cell point. bilinear returns interpolated values from the values of the four nearest raster cells. Only useful on a flat area.
bio1 <- extract(climate, coord, method="simple", buffer=NULL, fun=mean, df=TRUE) 
head(bio1)



dat2 <-bio[,-1]
temp <- dat2[,c(1:3,12:19)]
prec <- dat2[,c(4:11)]

temp <- temp/10
head(temp)

# Temperature
for (i in 1:length(temp)){
  for (j in 1:dim(temp)[1]){
    if (temp[j,i] >= 60000) temp[j,i] <- temp[j,i] - 65535	#negative values are expressed as 65535-x
  }
  temp[,i] <- temp[,i]/10							#temperature in raw data is multiplied by 10
}

# Precipitation
for (i in 1:length(prec)){
  for (j in 1:dim(prec)[1]){
    if (prec[j,i] >= 60000) prec[j,i] <- prec[j,i] - 65535	#negative values are expressed as 65535-x
  }}

dat2 <- cbind(dat,temp,prec)

write.table(dat2, "trait_bioclim2.csv", sep=",", row.names = F)
write.table(dat2, "study_bioclim2.csv", sep=",", row.names = F)

write.table(dat2, "trait_all_bioclim.csv", sep=",", row.names = F)




# GET CLIMATE DATA - DOES ONLY WORK FOR HIGER RESOLUTION OT TILES
bioclim <- getData('worldclim', var='bio', res=2.5)
plot(bioclim$bio1, main="Annual Mean Temperature")
bioclim <- getData('worldclim', var='bio', res=0.5, lon=c(-175, 171), lat=c(-45, 69))






# extract coords for Frei and some other sites

site <- read.csv("Pop_Frei.csv", header=TRUE, sep=";")
head(site)
coord <- as.data.frame(cbind(site$Long, site$Lat))
points(coord$V1, coord$V2, col="black", pch=1)


setwd("G:/Worldclim Data/bio1-9_30s_bil")
# create a list of .bil files that exist in the wd
files <- list.files(pattern='\\.bil$')
# combine all list elements into a stack
bioclim <- stack(files)
# quick plot to make sure nothing went drastically wrong
plot(bioclim$bio_1)

# extract data
bio1 <- extract(bioclim, coord, method="simple", buffer=NULL, fun=mean, df=TRUE) 

data.new <- cbind(site, bio1)
head(data.new)
write.table(dat.new, "Frei_bioclim.csv", sep=",", row.names = F)

