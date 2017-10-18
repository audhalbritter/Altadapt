###########################################################################
######################### META ANALYSIS ALTADAPT ########################## 16-02-11
###########################################################################

#Libraries
library("geosphere"); library("arm"); library("multcomp"); 
library("tidyverse"); library("lubridate"); library("broom")

source("Functions Metaanalysis.R")

########################
### IMPORT DATA
########################
# IMPORT TRAIT DATA
DAT.T <- read.csv("160210_Trait_Bioclim.csv",header=T, fill=TRUE, sep=";")
DAT.T <- DAT.T %>%
  filter(EXCLUDE == 0) %>% # exclude non relevant studies
  # Create new variable for each studyunit (Trait within species within Site) and ID (species within StudyID)
  mutate(studyunit = paste(SiteID, species, Trait, sep="_")) %>% 
  mutate(ID = paste(SiteID, species, sep="_")) %>% # ID by SiteID and species to join with study data
  arrange(studyunit)


# Calculate mean and range elev, lat and long and mean temperature (bio_1) and spatial distance for each studyunit
new.means2 <- get.means(DAT.T)
new.means2 <- data.frame(new.means2)
new.means2 <- rownames_to_column(new.means2, "studyunit")

# IMPORT STUDY DATA
DAT.S <- read.table("160210_Study_Bioclim.csv",header=T, fill=TRUE, sep=";")
DAT.S <- DAT.S %>%
  mutate(ID = paste(SiteID, species, sep="_")) %>% # create ID for species within StudyID
  arrange(ID) %>% 
  filter(excluded != "1") %>% # exclude "bad" studies
  # add site to bioclim vars to distinguish between pop and site vars
  rename(bio_1_site = bio_1, bio_2_site = bio_2, bio_3_site = bio_3, bio_4_site = bio_4, bio_5_site = bio_5, bio_6_site = bio_6, bio_7_site = bio_7, bio_8_site = bio_8, bio_9_site = bio_9, bio_10_site = bio_10, bio_11_site = bio_11, bio_12_site = bio_12, bio_13_site = bio_13, bio_14_site = bio_14, bio_15_site = bio_15, bio_16_site = bio_16, bio_17_site = bio_17, bio_18_site = bio_18, bio_19_site = bio_19) %>% 
  droplevels()
#colnames(DAT.S)[33:51] <- paste(grep("bio", colnames(DAT.S), value = TRUE), "site", sep = "_")

length(unique(DAT.T$ID)) # check length of data set
setdiff(DAT.S$ID,DAT.T$ID) # check if T and S have the same StudyIDs


# Select TRAIT_CAT1-4 and elev.pop, traitmean
DAT.Trait <- DAT.T %>% 
  select(studyunit, Year, TRAIT_CAT1, TRAIT_CAT2, TRAIT_CAT3, TRAIT_CAT4) %>% 
  distinct()

DAT.Site <- DAT.S %>% 
  select(ID, studysite, Coord, coor.var, family, generation, duration_exp, breedingsystem, longevity, growthform, intro, country_of_origin, country_of_exp, region, elevation_site, breed)

# CREATE META.STUDYUNIT
meta.studyunit <- DAT.T %>% 
  distinct(studyunit, ID, StudyID, SiteID, species, Trait, study_type) %>% 
  left_join(new.means2, by = c("studyunit")) %>%  # left join new.means2
  # merge trait and site data
  left_join(DAT.Site, by = "ID") %>% 
  left_join(DAT.Trait, by = "studyunit")


# CALCULATE SLOPES
## standardise trait values within studyunits to mean 0 sd 1 and extract regression coefs of trait on mean temperature (bio_1)
dat.slope <- get.slopes2(DAT.T)
dat1 <- meta.studyunit %>% left_join(dat.slope, by = "studyunit")
dat2 <- dat1


# EXCLUDING DATA
dat2 <- dat2 %>% 
  # Excluding studies with mean.dist.km < 1km; Gonz2009JEco, Haut2009JPlE, Sund1995Scan
  filter(mean.dist.km > 1) %>% 
  # Excluding studies with dist.km > 1000km; Ginw2004aSil, Ginw2004bSil, Gomo2011AnPS, Larw2010Silv, Nels1967NePh
  filter(dist.km < 1000) %>% 
  # only used when taking temp data (we only want studies where coord differes for pops)
  filter(coor.var == "1") %>% 
  # Exclude data studies from Rehfeld, where coords are only extracted for some regions
  filter(Coord != "extracted_rough_2") %>% 
  # delete Studies without a slope
  filter(!is.na(estimate))
head(dat2)
length(unique(dat2$StudyID)) # 70


save(dat2, file = "dat2.Rdata")


# subsetting data for height, biomass and phenology, create a variable var
# create variables for absolute slope and latitude
dat4 <- dat2 %>% 
  filter(TRAIT_CAT1 == "SIZE" & TRAIT_CAT2 == "HEIGHT" | TRAIT_CAT1 == "SIZE" & TRAIT_CAT2 == "BIOMASS" & TRAIT_CAT3 == "ABOVEGR" | TRAIT_CAT1 == "PHENOLOGY" & TRAIT_CAT2 %in% c("FLOWERING", "LEAF_BUD") & TRAIT_CAT4 == c("first", "peak")) %>%
  mutate(var = ifelse(TRAIT_CAT1 == "SIZE" & TRAIT_CAT2 == "HEIGHT", "height",
                      ifelse(TRAIT_CAT1 == "SIZE" & TRAIT_CAT2 == "BIOMASS" & TRAIT_CAT3 == "ABOVEGR", "biomass", "phenology"))) %>% 
  rename(slopes = estimate, sample.size = n) %>% 
  mutate(abs.slopes = abs(slopes)) %>% 
  mutate(abs.lat = abs(mean.lat))
head(dat4)

save(dat4, file = "dat4.Rdata")

# Excluding climate chamber and greenhouse experiments
noGarden <- dat4 %>% filter(studysite == "Garden")

fit1 <- lm(abs(slopes) ~ sample.size, phenology)
fit2 <- lm(abs(slopes) ~ sample.size + I(sample.size^2), phenology)
AIC(fit1, fit2)
summary(fit2)


# CORRELATION BETWEEN ELEV AND BIO_1
# exclude studies with only 2 pops

CorBioElev <- DAT.T %>% 
  group_by(StudyID, studyunit) %>% 
  filter(!is.na(bio_1)) %>%
  mutate(n = n()) %>% 
  filter(n > 2) %>% 
  summarise(corBioElev = cor(bio_1, elev))

dat2 %>% 
  inner_join(CorBioElev, by = "studyunit") %>% 
  summarise(n = n(), mean = mean(corBioElev), se = sd(corBioElev)/sqrt(n))

DAT.T %>% filter(StudyID %in% c("Frei2014GlCB"))
DAT.T %>% filter(StudyID %in% c("Nels1967NePh", "Park2003ConG", "Rehf1990BoGa")) %>%
  left_join(ss, by = "StudyID") %>% 
  filter(sample.size > 2) %>% 
  ggplot(aes(x = elev, y = bio_1, color = StudyID)) +
  geom_point() +
  facet_wrap(~region)


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



# RESULTS
dat2 %>% 
  select(r.temp, dist.km, r.elev) %>% 
  gather(key = var, value = value) %>% 
  group_by(var) %>% 
  summarise(min = min (value, na.rm = TRUE), max = max (value, na.rm = TRUE), mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE))



dat4 %>% select(species, growthform) %>% distinct(species, growthform) %>% group_by(growthform) %>% summarise(n = n())
dat2 %>% distinct(StudyID, intro) %>% group_by(intro) %>% summarise(n = n(), percent = n * 100 / 71)
dat2 %>% distinct(StudyID, region) %>% group_by(region) %>% summarise(n = n(), percent = n * 100 / 70)

# Sample size
dat2 %>% distinct(StudyID, n) %>% group_by(n) %>% summarise(nr = n()) %>% filter(n < 10) %>% summarise(sum = sum(nr))

dat4 %>% 
  group_by(var) %>% 
  summarise(median = median(sample.size))

# proportion sign. p value
dat4 %>%
  group_by(var) %>% 
  summarise(n = n(), ss = sum(p.value < 0.05), percent = ss * 100 / n)


# mean and CI
dat4 %>% group_by(var) %>% summarise(n = n(), mean = mean(slopes), se = sd(slopes)/sqrt(n), CI.low = mean - 1.96*se, CI.high = mean + 1.96*se)


dd <- read.csv("dat2.csv", header=TRUE, sep=";")
dd <- dat2
unique(dd$StudyID)

dd <- dat2
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
ddd <- unique(dat4[,c("StudyID", "study_type")])
plyr::count(ddd$study_type)
ddd[with(ddd, order(StudyID)), ]

unique(dat4[,c("StudyID", "SiteID", "study_type", "studysite")])
ddd <- unique(dat2[,c("SiteID", "studysite")])
plyr::count(ddd$studysite)
table(dat2$SiteID, dat2$studysite)

library(dplyr)
library(tidyr)
dat2 %>% 
  select(StudyID, SiteID, studysite) %>% 
  group_by(StudyID, studysite) %>% 
  summarise(n = n()) %>% 
  spread(key = studysite, value = n) %>% 
  print(n = 80)

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
dd <- dat2 %>% select(study_type, growthform, StudyID, species, studysite) %>% 
  distinct() %>% arrange(study_type, growthform, StudyID, species) %>% 
  mutate(studysite = plyr::mapvalues(studysite, c("Garden", "climate_chamber", "greenhouse"), c("G", "CC", "GH"))) %>% 
  select(-study_type)

write.csv(dd, "dd.csv", row.names = FALSE)

# check significance
dat4 %>%
  mutate(significant = ifelse(p.slope < 0.05, "sign", "nonsign")) %>% 
  group_by(var, significant) %>% 
  summarise(n = n()) %>% 
  mutate(freq = n / sum(n) * 100)


#### COMPARE LINEAR AND QUADRATIC MODEL (R^2) ####
DAT.T %>% 
  group_by(StudyID, studyunit) %>% 
  filter(!is.na(traitmean)) %>%
  mutate(n = n()) %>% 
  filter(n > 2) %>% 
  mutate(std.traitmean = scale(traitmean)) %>% 
  group_by(n, add = TRUE) %>% # add n for sample size
  do({
    # fit model
    fit1 <- lm(std.traitmean ~ bio_1, data = .)
    fit2 <- lm(std.traitmean ~ bio_1 + I(bio_1^2), data = .)
    anova(fit1, fit2)
    #delta.aic <- AIC(fit1) - AIC(fit2)
    #res <- cbind(p.value, delta.aic)
    #return(res)

  })


neudat <- DAT.T[DAT.T$StudyID %in% dat4$StudyID,]
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
# in 52 out of 1042 cases, the quadratic model is better!!! 4.99%