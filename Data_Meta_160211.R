###########################################################################
######################### META ANALYSIS ALTADAPT ########################## 16-02-11
###########################################################################

#Libraries
library(geosphere); library(lme4); library(arm); library(multcomp)

########################
### IMPORT DATA
########################
# IMPORT TRAIT DATA

DAT.T <- read.csv("160210_Trait_Bioclim.csv",header=T, fill=TRUE, sep=";")
DAT.T <- subset(DAT.T, DAT.T$EXCLUDE==0) # exclude the studies we don't want
# Create new variable for each studyunit (Trait within species within Site) and ID (species within StudyID)
DAT.T$studyunit <- paste(DAT.T$SiteID, DAT.T$species, DAT.T$Trait, sep="_")
DAT.T$ID <- paste(DAT.T$SiteID, DAT.T$species, sep="_")
DAT.T <- DAT.T[order(DAT.T$studyunit),]

# Calculate mean and range elev, lat and long and mean temperature (bio_1) and spatial distance for each studyunit
new.means2 <- get.means(DAT.T)

# IMPORT STUDY DATA
DAT.S <- read.table("160210_Study_Bioclim.csv",header=T, fill=TRUE, sep=";")
DAT.S <- DAT.S[order(paste(DAT.S$SiteID, DAT.S$species, sep="_")),]
DAT.S <- subset(DAT.S, DAT.S$excluded!="1") # exclude "bad" studies
DAT.S <- droplevels(DAT.S)
DAT.S$ID <- paste(DAT.S$SiteID, DAT.S$species, sep="_") # create ID for species within StudyID
# add site to bioclim vars to distinguish between pop and site vars
colnames(DAT.S)[33:51] <- c(paste(colnames(DAT.S)[33:51], "site", sep="_"))
length(unique(DAT.T$ID)) # check length of data set
setdiff(DAT.S$ID,DAT.T$ID) # check if T and S have the same StudyIDs

# CREATE META.STUDYUNIT
meta.studyunit <- unique(DAT.T[,c("studyunit", "ID", "StudyID", "SiteID", "species", "Trait", "study_type")])
head(meta.studyunit)
meta.studyunit$studyunit==rownames(new.means2) # check if the two data sets have the same order!
#meta.studyunit <- cbind(meta.studyunit, new.means2) # combine meta data with calculated means
# if line above does not work, the next two lines do.
meta.studyunit[,(ncol(meta.studyunit)+1):(ncol(meta.studyunit)+ncol(new.means2))] <- new.means2[,1:ncol(new.means2)]
colnames(meta.studyunit)[8:16] <- c("mean.elev", "r.elev", "mean.lat", "mean.long", "mean.temp", "r.temp", "dist.dd", "dist.km", "mean.dist.km")

# add TRAIT_CAT1-4 and elev.pop, traitmean, bio_1
meta.studyunit[,17:23] <- DAT.T[match(meta.studyunit$studyunit,DAT.T$studyunit),c("TRAIT_CAT1", "TRAIT_CAT2", "TRAIT_CAT3", "TRAIT_CAT4", "elev", "traitmean", "bio_1")]
# add studysite, family, breedingsystem, longevity, growthform, elevation_site, generation, region and intro
meta.studyunit[,24:38] <- DAT.S[match(meta.studyunit$ID,DAT.S$ID),c("studysite", "Coord", "coor.var", "family", "generation", "duration_exp", "breedingsystem", "longevity", "growthform", "intro", "country_of_origin", "country_of_exp", "region", "elevation_site", "breed")]

# CALCULATE SLOPES
## standardise trait values within studyunits to mean 0 sd 1 and extract regression coefs of trait on mean temperature (bio_1)
dat1 <- get.slopes(DAT.T, meta.studyunit)
dat2 <- dat1

# EXCLUDING DATA
# Excluding studies with mean.dist.km < 1km and dist.km > 1000km
dat2 <- subset(dat2, dat2$mean.dist.km > 1) # Gonz2009JEco, Haut2009JPlE, Sund1995Scan
dat2 <- subset(dat2, dat2$dist.km < 1000) # Ginw2004aSil, Ginw2004bSil, Gomo2011AnPS, Larw2010Silv, Nels1967NePh
dat2 <- subset(dat2, dat2$coor.var=="1")# only used when taking temp data (we only want studies where coord differes for pops)
# Exclude data studies from Rehfeld, where coords are only extracted for some regions
dat2 <- subset(dat2, dat2$Coord!="extracted_rough_2")

# delete Studies without a slope
dat2 <- subset(dat2, dat2$slopes!="NA")
length(unique(dat2$StudyID)) # 80

# Excluding climate chamber and greenhouse experiments
#dat2 <- subset(dat2, dat2$studysite=="Garden")

#write.csv(dat2, file = "dat2_160212.csv",  row.names=FALSE)



# CORRELATION BETWEEN ELEV AND BIO_1
# exclude studies with only 2 pops
dat.new <- dat2[dat2$sample.size > 2,]
range(dat.new$cor.elev.bio_1, na.rm=TRUE)
sort(dat.new$cor.elev.bio_1)
ddd <- dat.new$cor.elev.bio_1 < -0.7
length(ddd[ddd==TRUE])
dat07 <- subset(dat.new, cor.elev.bio_1 < -0.7)
mean(dat07$cor.elev.bio_1)
sd(dat07$cor.elev.bio_1)/sqrt(length(dat07$cor.elev.bio_1))


head(dat2)
dat3 <- subset(dat2, dat2$TRAIT_CAT1 == "FITNESS")
dat3$TRAIT_CAT2 <- factor(dat3$TRAIT_CAT2)
tab <- table(dat3$studyunit, dat3$TRAIT_CAT2)
head(tab)
length(which(tab[,1]>0))
length(which(tab[,2]>0))
length(which(tab[,3]>0))
length(which(tab[,4]>0))
length(which(tab[,5]>0))
length(which(tab[,6]>0))
length(which(tab[,7]>0))
length(which(tab[,8]>0))
length(which(tab[,9]>0))
length(which(tab[,10]>0))
length(which(tab[,11]>0))
length(which(tab[,13]>0))



# SUMMARIES
dd <- read.csv("dat2.csv", header=TRUE, sep=";")
dd <- dat2
unique(dd$StudyID)

by(dd, dd$StudyID, function(id){
  study <- unique(id$study_type)
  site <- unique(id$studysite)
  sp <- unique(id$species)
  gf <- unique(id$growthform)
  ss <- mean(id$sample.size)
  return(gf)
})

length(unique(dat2$StudyID))
length(unique(dat2$SiteID))
length(unique(dat2$studyunit))
length(unique(dat2$species))

length(unique(dat2$family))
sort(table(dat2$family))
grow.f <- t(table(dat2$growthform, dat2$family))

# growthform
ddd <- unique(dat2[,c("species", "growthform")])
dim(ddd)
dd <- count(ddd$growthform)
dd$freq*100/80

# Intro
ddd <- unique(dat2[,c("species", "intro")])
ddd[with(ddd, order(intro)), ]
dim(ddd)
library(plyr)
count (ddd$intro)

# Regions
ddd <- unique(dat2[,c("StudyID", "region")])
dd <- count(ddd$region)
dd$freq*100/80

# study_type
ddd <- unique(dat2[,c("StudyID", "study_type")])
count(ddd$study_type)
ddd[with(ddd, order(StudyID)), ]


ddd <- unique(dat2[,c("StudyID", "studysite")])
count(ddd$studysite)

# Temperature range
(sort(dat2$r.temp))
ddd <- na.omit(dat2$r.temp)
mean(ddd)
sd(ddd)

(sort(dat2$r.elev))
(sort(dat2$dist.km))


names(LA2)
LA2[,59:61] <- meta.studyunit[match(LA2$ID,meta.studyunit$ID),c("r.elev", "r.temp", "dist.km")]
(sort(LA2$r.temp))
mean(LA2$r.temp)
sd(LA2$r.temp)
(sort(LA2$r.elev))
(sort(LA2$dist.km))


### Make a summary table
meta.dat2 <- unique(dat2[,c("study_type", "growthform", "species", "StudyID", "studysite")])
meta.dat2 <- meta.dat2[with(meta.dat2, order(study_type, growthform, species)),]
head(meta.dat2)
write.csv(meta.dat2, "meta.dat.csv", row.names = FALSE)




#### COMPARE LINEAR AND QUADRATIC MODEL (R^2) ####
neudat <- DAT.T[DAT.T$StudyID %in% dat2$StudyID,]
compare.lin.quad <- by(neudat, neudat$studyunit, function(x){
  standardized.traitmean <- (x$traitmean - mean(x$traitmean, na.rm=TRUE))/sd(x$traitmean, na.rm=TRUE)
  # get r2 for comparison of linear and quadratic regression
  fit1 <- tryCatch(lm(standardized.traitmean ~ bio_1, x), error=function(e) NA)
  fit2 <- tryCatch(lm(standardized.traitmean ~ bio_1+I(bio_1^2), x), error=function(e) NA)
  p.value <- tryCatch(anova(fit1, fit2)$"Pr(>F)"[2], error=function(e) NA)
  sample.size <- tryCatch((length(resid(lm(standardized.traitmean~bio_1,x)))), error=function(e) NA)
  aic.value <- tryCatch(AIC(fit1)-AIC(fit2), error=function(e) NA)
  cbind(p.value, sample.size, aic.value)
})
comp <- t(sapply(compare.lin.quad, I))
colnames(comp) <- c("p.value", "sample.size", "aic.value")
comp2 <- data.frame(comp)
names(comp2)
plot(comp2$p.value, col=ifelse(comp2$p.value < 0.05,"red", "black"))
quad <- subset(comp2, comp2$p.value < 0.05)
dim(quad)[1]*100/dim(comp2)[1]
# in 54 out of 1048 cases, the quadratic model is better!!! 5.2%