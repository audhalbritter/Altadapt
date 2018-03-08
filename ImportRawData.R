###########################################################################
######################### META ANALYSIS ALTADAPT ########################## 
###########################################################################

#### READ IN RAW DATA ####
#### CALCULATE SLOPES ####
#### MERGE WITH METADATA ####

#Libraries
library("geosphere"); library("arm"); library("multcomp"); 
library("tidyverse"); library("lubridate"); library("broom")

source("Functions Metaanalysis.R")

########################
### IMPORT DATA
########################
# IMPORT TRAIT DATA
DAT.T <- read.csv("data/160210_Trait_Bioclim.csv",header=T, fill=TRUE, sep=";")
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
DAT.S <- read.table("data/160210_Study_Bioclim.csv",header=T, fill=TRUE, sep=";")
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

#write.csv(dat2, file = "FinalSharedData/TraitDifferentiation.csv")
#save(dat2, file = "dat2.Rdata")