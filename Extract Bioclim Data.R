### EXTRACT BIOCLIM VARIABLES FOR POP AND SITE COORDINATES
getwd()

# READ IN SITE DATA
setwd("/Users/audhalbritter/Dropbox/121005_analysis/ANALYSIS")
site <- read.csv("160210_Altadapt_Study.csv", header=TRUE, fill=TRUE, sep=";")
site <- read.csv("NewCoords.csv", header=TRUE, fill=TRUE, sep=";")
head(site)
str(site)

long.site <- site$longitude_site
lat.site <- site$latitude_site
coord.s <- as.data.frame(cbind(long.site, lat.site))
#coord.s <- na.omit(coord.s)
# possible to draw points on the map to check the locations
points(coord.s$long.site, coord.s$lat.site, col="blue", pch=16)


# EXTRACT WORLDCLIM VARS
# set wd where the data is stored
setwd("/Volumes/SILVER/Worldclim Data/bio1-9_30s_bil")

library("raster")
# create a list of .bil files that exist in the wd
files <- list.files(pattern='\\.bil$')
# combine all list elements into a stack
bioclim <- stack(files)
# quick plot to make sure nothing went drastically wrong
plot(bioclim$bio_1)
dim(bioclim)

# EXTRACT
# fun accepts na.rm=T, but then you need to use function(x)...
# Buffer defines a radius around the point in which a value is included or not.
# If df=TRUE, the results are returned as a data frame
# Method: simple returns a single value for the cell point. bilinear returns interpolated values from the values of the four nearest raster cells. Only useful on a flat area.

# extract bioclim data for sites
bio.s <- extract(bioclim, coord.s, method="simple", buffer=NULL, fun=mean, df=TRUE) 
head(bio.s)
dat2 <-bio.s[,-1]
dat2 <- cbind(dat2[,c(1, 12:19, 2:3)]/10, dat2[,c(4:11)])

site.new <- cbind(site, dat2)
setwd("/Users/audhalbritter/Dropbox/121005_analysis/ANALYSIS")
write.table(site.new, "160210_Study_Bioclim.csv", sep=",", row.names = FALSE)

dd <-bio.s[,-1]
dd <- cbind(dd[,c(1, 12:19, 2:3)]/10, dd[,c(4:11)])
dd <- cbind(site, dd)
write.table(dd, "bioclim.new.txt", row.names = FALSE)

# READ IN TRAIT DATA FOR POPS
setwd("/Users/audhalbritter/Dropbox/121005_analysis/ANALYSIS")
pop <- read.csv("160210_Altadapt_Traits.csv", header=TRUE, fill=TRUE, sep=";")
head(pop)
str(pop)
pop$long.pop

long.pop <- pop$long.pop
lat.pop <- pop$lat.pop
coord.p <- as.data.frame(cbind(long.pop, lat.pop))
# coord.p <- na.omit(coord.p)
# possible to draw points on the map to check the locations
points(coord.p$long.pop, coord.p$lat.pop, col="red", pch=16)


# extract bioclim data for pops
bio.p <- extract(bioclim, coord.p, method="simple", buffer=NULL, fun=mean, df=TRUE) 
head(bio.p)
dat2 <-bio.p[,-1] # remove first column
dat2 <- cbind(dat2[,c(1, 12:19, 2:3)]/10, dat2[,c(4:11)]) # devide temp variables by 10 and add prec variables

pop.new <- cbind(pop, dat2)
setwd("/Users/audhalbritter/Dropbox/121005_analysis/ANALYSIS")
write.table(pop.new, "160210_Trait_Bioclim.csv", sep=",", row.names = FALSE)





### NOT USED
dat2 <-bio.s[,-1]
ddd <- ifelse(dat2$bio_1 > 6000, dat2$bio_1-65535, dat2$bio_1)
ddd <- ifelse(dat2 > 6000, dat2-65535, dat2)

head(dat2)
temp <- temp[,c(1:3)]/10

for (i in 1:length(dat2)){
  for (j in 1:dim(dat2)[1]){
    if (dat2[j,i] >= 60000) dat2[j,i] <- dat2[j,i] - 65535	#negative values are expressed as 65535-x
    eles(dat2[j,i])
  }
}


  



# GET CLIMATE DATA - DOES ONLY WORK FOR HIGER RESOLUTION OT TILES
bioclim <- getData('worldclim', var='bio', res=2.5)
plot(bioclim$bio1, main="Annual Mean Temperature")
bioclim <- getData('worldclim', var='bio', res=0.5, lon=c(-175, 171), lat=c(-45, 69))


