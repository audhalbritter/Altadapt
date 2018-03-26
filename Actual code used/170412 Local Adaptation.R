#### LOCAL ADAPTATION ####

library("lme4")
library("tidyverse")
library("MuMIn")
library("cowplot")

# Only inflorescence fecundity
#Plue2011PlED site 4-7
# ifelse(s %in% c(4, 5), "low", "high")
# ifelse(p == 1, "low", "high")
# Site 4 and 5 = 600m
# Site 6 and 7 = 1200m
# Pop 1 = 600m
# Pop 2 = 1200m

Plue2011 <- DAT.T %>%
  filter(SiteID %in% c("Plue2011PlED.4", "Plue2011PlED.5", "Plue2011PlED.6", "Plue2011PlED.7")) %>% 
  filter(Trait %in% c("Seed_fecundity_2009", "Seed_fecundity_2010")) %>% 
  mutate(study_type = "reciprocal.transplant") %>% 
  mutate(data_type = "transplant") %>% 
  # add speceis and site info
  left_join(DAT.S, by = c("ID", "StudyID", "SiteID", "species", "data_type")) %>% 
  
  # extract last character from site and pop
  mutate(LastSiteID = substr(as.character(SiteID), nchar(as.character(SiteID)), nchar(as.character(SiteID)))) %>% 
  mutate(LastPopID = substr(as.character(PopID), nchar(as.character(PopID)), nchar(as.character(PopID)))) %>% 
  mutate(SiteID = ifelse(LastSiteID %in% c("4", "5"), "l", "h")) %>% 
  mutate(PopID = ifelse(LastPopID == "1", "l", "h")) %>% select(-LastSiteID, -LastPopID)



#Frei2014GlCB all  
Frei2014 <- DAT.T %>%
  # exclude site at 600m
  filter(SiteID %in% c("Frei2014GlCB.2", "Frei2014GlCB.3", "Frei2014GlCB.5", "Frei2014GlCB.6", "Frei2014GlCB.8", "Frei2014GlCB.9")) %>% 
  mutate(study_type = "reciprocal.transplant") %>%
  # add speceis and site info
  left_join(DAT.S, by = c("ID", "StudyID", "SiteID", "species", "data_type")) %>% 
  
  # get everything before the first m for SiteID
  mutate(SiteID = sapply(strsplit(as.character(PopID), split = "m"), "[", 1)) %>% 
  mutate(SiteID = ifelse(SiteID == "1200", "l", "h")) %>% 
  
  # get last four strings for PopID
  mutate(PopID = substr(as.character(PopID), nchar(as.character(PopID))-4+1, nchar(as.character(PopID)))) %>% 
  mutate(PopID = ifelse(PopID == "1200", "l", "h")) 



# SUBSET DATA
LA <- DAT.T %>%
  #mutate(ID = factor(ID), StudyID = factor(StudyID), SiteID = factor(SiteID), species = factor(species), data_type = factor(data_type)) %>% 
  # only use reciprocal transplant experiments
  filter(study_type == "reciprocal.transplant") %>% 
  # add speceis and site info
  left_join(DAT.S, by = c("ID", "StudyID", "SiteID", "species", "data_type")) %>% 
  
  # Rbind Plue and Frei
  bind_rows(Plue2011) %>% 
  bind_rows(Frei2014) %>% 
  
  # add r.elev, mean.dist, etc.
  left_join(new.means2, by = "studyunit") %>% 
  
  # Excluding studies with mean.dist.km < 1km; Gonz2009JEco
  filter(mean.dist.km > 1) %>% 
  # Exclude studies where pops do not differ in coordinates
  filter(StudyID != "Gaut1998NewP") %>% 
  # delete rows with no traitmean
  filter(!is.na(traitmean)) %>% 
  # only variables: biomass, fitness and survival
  filter(TRAIT_CAT1 == "SIZE" & TRAIT_CAT2 == "BIOMASS" & TRAIT_CAT3 == "ABOVEGR" | TRAIT_CAT1 == "FITNESS" | TRAIT_CAT1 == "VITAL_RATES"& TRAIT_CAT2 == "SURVIVAL") %>%
  mutate(var = ifelse(TRAIT_CAT1 == "SIZE" & TRAIT_CAT2 == "BIOMASS" & TRAIT_CAT3 == "ABOVEGR", "biomass",
                      ifelse(TRAIT_CAT1 == "FITNESS", "fitness", "survival"))) 

# SYMPATRY- ALLOPATRY
# create new factor for expressing whether a population is growing in its own habitat ("sympatric" = 1) or not ("allopatric" = 0)
LA <- LA %>% 
  # calculate elev between site and pop (only useful for Gardens) and make it a factor (High/Low)
  mutate(elev.diff = elevation_site - elev) %>% 
  
  mutate(symp = substr(SiteID, nchar(as.character(SiteID)), nchar(as.character(SiteID))) ==
           substr(PopID, nchar(as.character(PopID)), nchar(as.character(PopID)))) %>% 
  mutate(symp = factor(symp + 0)) %>% # symp = 1, allo = 0
  # direciton of transplant: elev.diff > 0 pop has been transplanted up
  mutate(direction = ifelse(symp == 1, NA,
                            ifelse(elev.diff > 0, "up", "down"))) %>% 
  mutate(direction = factor(direction)) %>% 
  mutate(symp = plyr::mapvalues(symp, c("0", "1"), c("foreign site", "home site"))) %>% 
  mutate(symp = factor(symp, levels = (c("home site", "foreign site"))))




# Create site information: which site is low, mid and high. factor.site2 is for the second analysis to check interaction symp * site (needs to exclude mid elevation sites, but for Halbritter studies high elevation sites are different for 3 species)
site.info <- LA %>% 
  select(StudyID, SiteID, elevation_site, species) %>% 
  unique() %>% 
  group_by(StudyID) %>% 
  mutate(low.site = min(elevation_site), high.site = max(elevation_site)) %>% 
  mutate(factor.site = ifelse(elevation_site == low.site, "low", ifelse(elevation_site == high.site, "high", "mid"))) %>%
  mutate(factor.site2 = factor.site) %>% 
  mutate(factor.site2 = ifelse(StudyID == "Halb0000subm" & species == "Plantago_lanceolata" & factor.site == "mid", "high", factor.site2)) %>% 
  mutate(factor.site2 = ifelse(StudyID == "Halb0000subm" & species == "Plantago_lanceolata" & factor.site == "high", "mid", factor.site2)) %>% 
  mutate(factor.site2 = ifelse(StudyID == "Halb0001subm" & species %in% c("Plantago_lanceolata", "Senecio_viscosus", "Medicago_lupulina") & factor.site == "mid", "high", factor.site2)) %>% 
  mutate(factor.site2 = ifelse(StudyID == "Halb0001subm" & species %in% c("Plantago_lanceolata", "Senecio_viscosus", "Medicago_lupulina") & factor.site == "high", "mid", factor.site2)) %>% 
  mutate(factor.site = ifelse(StudyID == "Plue2011PlED" & elevation_site %in% c(626, 600), "low", ifelse(StudyID == "Plue2011PlED" & elevation_site %in% c(1251, 1235), "high", factor.site))) %>% 
  mutate(factor.site2 = ifelse(StudyID == "Plue2011PlED" & elevation_site %in% c(626, 600), "low", ifelse(StudyID == "Plue2011PlED" & elevation_site %in% c(1251, 1235), "high", factor.site2))) %>% select(-low.site, -high.site)


# add population as factor (High/Mid/Low)
elev_factor <- read.table("data/elev_factor.csv", sep=";", header=TRUE)
elev_factor <- elev_factor %>% distinct(StudyID, elev, elevation_site, factor.pop)

# -----------------------------------------------
# Version 1
#st.traitmean <- by(LA, LA$studyunit, function(x){
  #standardized.traitmean <- scale(x$traitmean)
  #standardized.traitmean <- (x$traitmean - mean(x$traitmean, na.rm=TRUE))/sd(x$traitmean, na.rm=TRUE)
#})

# Check if studyunit are the same
#st <- plyr::ldply(st.traitmean, data.frame)
#st <- setNames(st, c("studyunit", "st.traitmean"))
#st$studyunit == LA$studyunit
#setdiff(st$studyunit, LA$studyunit)

# join standardized traitmean
#st.traitmean <- as.vector(unlist(st.traitmean))
#LA <- cbind(LA, st.traitmean)
# -----------------------------------------------

# join site info
LA <- LA %>% 
  #left_join(standardized.traitmean, by = "studyunit") %>% 
  left_join(site.info, by = c("StudyID", "SiteID", "species", "elevation_site")) %>% 
  left_join(elev_factor, by = c("StudyID", "elev", "elevation_site")) %>% 
  mutate(factor.pop = factor(factor.pop, levels = c("low", "mid", "high"))) %>% 
  mutate(factor.site = factor(factor.site, levels = c("low", "mid", "high"))) %>% 
  mutate(factor.site2 = factor(factor.site2, levels = c("low", "mid", "high"))) %>% 
  # Version 2: Standardize traitmean
  group_by(studyunit) %>% 
  mutate(st.traitmean = (traitmean - mean(traitmean, na.rm=TRUE))/sd(traitmean, na.rm=TRUE)) %>% 
  mutate(sample.size = n())
  

LA <- LA %>% 
  select(-SEQ, -data_type, -TRAIT_CAT1, -TRAIT_CAT2, -TRAIT_CAT3, -TRAIT_CAT4, -YearOfExp, -Data.entered.by, -EXCLUDE, -bio_1, -bio_2, -bio_3, -bio_4, -bio_5, -bio_6, -bio_7, -bio_8, -bio_9, -bio_10, -bio_11, -SEQ1, -extracted_by, -lead_author, -year, -journal, -studytype, -studysite, -no.sites, -no.traits, -Coord, -coor.var, -LATLONG_extr, -breedingsystem, -excluded, -why_exclude, -bio_1_site, -bio_2_site, -bio_3_site, -bio_4_site, -bio_5_site, -bio_6_site, -bio_7_site, -bio_8_site, -bio_9_site, -bio_10_site, -bio_11_site, -bio_12_site, -bio_13_site, -bio_14_site, -bio_15_site, -bio_16_site, -bio_17_site, -bio_18_site, -bio_19_site, -mean.elev, -r.elev, -mean.lat, -mean.long, -mean.temp, -mean.seasonality, -dist.dd, -dist.km, -mean.dist.km)

write.csv(LA, file = "FinalSharedData/LocalAdaptation.csv")

s### RESULTS
# count nr of studies
LA %>% group_by(var) %>% distinct(StudyID) %>% summarise(n())
summary(LA)
LA %>% group_by(var) %>%  summarise(median(sample.size))


####################################
#### TEST FOR LOCAL ADAPTATION ####
####################################

# BIOMASS
biomass <- LA %>% filter(var == "biomass") %>% filter(!is.na(st.traitmean))
fit1 <- lmer((st.traitmean) ~ symp * r.temp + (1|PopID)+(1|species)+(1|StudyID), biomass, na.action = "na.fail", REML = FALSE) 
#fit2 <- lmer((st.traitmean) ~ 1 + (1|PopID)+(1|species)+(1|StudyID), biomass, na.action = "na.fail", REML = FALSE) # better model
#modsel(list(fit1, fit2),1000)
model.set <- dredge(fit1, rank = "AICc", extra = "R^2")
model.set
averaged.model <- model.avg(model.set)
res <- data.frame(summary(averaged.model)$coefmat.full)
BiomassLA <- res %>% 
  rownames_to_column(var = "Variable") %>% 
  setNames(., c("Variable", "Estimate", "StError", "AdjSE", "Zvalue", "Pvalue")) %>% 
  mutate(CI.low = Estimate - 1.96 * StError) %>% 
  mutate(CI.high = Estimate + 1.96 * StError) %>% 
  mutate(Estimate = round(Estimate, 2), CI = paste(round(CI.low, 2), round(CI.high, 2), sep = " - "), Zvalue = round(Zvalue, 2), Pvalue = round(Pvalue, 3)) %>% 
  select(Variable, Estimate, CI, Zvalue, Pvalue) %>%
  #mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "sympforeign site"), c("Intercept", "Own-Foreign")))
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "sympforeign site", "r.temp", "r.temp:sympforeign site"), c("Intercept", "Own-Foreign", "T Range", "Own-Foreign * T Range")))

summary(fit1)
fix.check(fit2)
anova(fit1, fit2)
plot(biomass$symp, biomass$st.traitmean)

# FITNESS
fitness <- LA %>% filter(var == "fitness") %>% filter(!is.na(st.traitmean))
hist((fitness$st.traitmean))
fit1 <- lmer((st.traitmean) ~ symp * r.temp + (1|PopID)+(1|species)+(1|StudyID), fitness, na.action = "na.fail", REML = FALSE) # better model
#fit2 <- lmer((st.traitmean) ~ 1 + (1|PopID)+(1|species)+(1|StudyID), fitness, na.action = "na.fail", REML = FALSE)
#modsel(list(fit1, fit2),1000)
model.set <- dredge(fit1, rank = "AICc", extra = "R^2")
model.set
averaged.model <- model.avg(model.set)
res <- data.frame(summary(averaged.model)$coefmat.full)
FitnessLA <- res %>% 
  rownames_to_column(var = "Variable") %>% 
  setNames(., c("Variable", "Estimate", "StError", "AdjSE", "Zvalue", "Pvalue")) %>% 
  mutate(CI.low = Estimate - 1.96 * StError) %>% 
  mutate(CI.high = Estimate + 1.96 * StError) %>% 
  mutate(Estimate = round(Estimate, 2), CI = paste(round(CI.low, 2), round(CI.high, 2), sep = " - "), Zvalue = round(Zvalue, 2), Pvalue = round(Pvalue, 3)) %>% 
  select(Variable, Estimate, CI, Zvalue, Pvalue) %>%
  #mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "sympforeign site"), c("Intercept", "Own-Foreign")))
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "sympforeign site", "r.temp", "r.temp:sympforeign site"), c("Intercept", "Own-Foreign", "T Range", "Own-Foreign * T Range")))

summary(fit1)
fix.check(fit1)

plot(fitness$symp, fitness$st.traitmean)


# SURVIVAL
survival <- LA %>% filter(var == "survival") %>% filter(!is.na(st.traitmean))
hist((survival$st.traitmean))
fit1 <- lmer((st.traitmean) ~ symp * r.temp + (1|PopID)+(1|species)+(1|StudyID), survival, na.action = "na.fail", REML = FALSE) # better model
#fit2 <- lmer((st.traitmean) ~ 1 + (1|PopID)+(1|species)+(1|StudyID), survival, na.action = "na.fail", REML = FALSE)
#modsel(list(fit1, fit2),1000)
model.set <- dredge(fit1, rank = "AICc", extra = "R^2")
model.set
averaged.model <- model.avg(model.set)
res <- data.frame(summary(averaged.model)$coefmat.full)
SurvivalLA <- res %>% 
  rownames_to_column(var = "Variable") %>% 
  setNames(., c("Variable", "Estimate", "StError", "AdjSE", "Zvalue", "Pvalue")) %>% 
  mutate(CI.low = Estimate - 1.96 * StError) %>% 
  mutate(CI.high = Estimate + 1.96 * StError) %>% 
  mutate(Estimate = round(Estimate, 2), CI = paste(round(CI.low, 2), round(CI.high, 2), sep = " - "), Zvalue = round(Zvalue, 2), Pvalue = round(Pvalue, 3)) %>% 
  select(Variable, Estimate, CI, Zvalue, Pvalue) %>%
  #mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "sympforeign site"), c("Intercept", "Own-Foreign")))
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "sympforeign site", "r.temp", "r.temp:sympforeign site"), c("Intercept", "Own-Foreign", "T Range", "Own-Foreign * T Range")))

summary(fit1)
fix.check(fit1)
plot(survival$symp, survival$st.traitmean)

ResultsLA <- BiomassLA %>% 
  bind_rows(FitnessLA, SurvivalLA) %>% 
  #mutate(Trait = c(rep("Biomass", 2), rep("Rep. output", 2), rep("Survival", 2))) %>% 
  mutate(Trait = c(rep("Biomass", 4), rep("Rep. output", 4), rep("Survival", 4))) %>% 
  select(Trait, Variable, Estimate, CI, Zvalue, Pvalue)


aggregate(st.traitmean ~ symp, biomass, mean)
aggregate(st.traitmean ~ symp, biomass, sd)
aggregate(st.traitmean ~ symp, fitness, mean)
aggregate(st.traitmean ~ symp, fitness, sd)
aggregate(st.traitmean ~ symp, survival, mean)
aggregate(st.traitmean ~ symp, survival, sd)
coef(fit1)

aggregate(st.traitmean ~ symp, LA.biomass, mean)
aggregate(st.traitmean ~ symp, LA.fitness, mean)
aggregate(st.traitmean ~ symp, LA.survival, mean)

################################################
#### DIFFERENCE BETWEEN HIGH AND LOW POPS ####
################################################
LA %>% group_by(var, breed) %>% summarise(n())

# BIOMASS
biomass2  <-  biomass %>% filter(factor.site2 != "mid")
fit1 <- lmer((st.traitmean) ~ symp*factor.site2 + (1|PopID)+(1|species)+(1|StudyID), biomass2, na.action = "na.fail", REML = FALSE) ### BEST MODEL
model.set <- dredge(fit1, rank = "AICc")
model.set
averaged.model <- model.avg(model.set)
res <- data.frame(summary(averaged.model)$coefmat.full)
BiomassHL <- res %>% 
  rownames_to_column(var = "Variable") %>% 
  setNames(., c("Variable", "Estimate", "StError", "AdjSE", "Zvalue", "Pvalue")) %>% 
  mutate(CI.low = Estimate - 1.96 * StError) %>% 
  mutate(CI.high = Estimate + 1.96 * StError) %>% 
  mutate(Estimate = round(Estimate, 2), CI = paste(round(CI.low, 2), round(CI.high, 2), sep = " - "), Zvalue = round(Zvalue, 2), Pvalue = round(Pvalue, 3)) %>% 
  select(Variable, Estimate, CI, Zvalue, Pvalue) %>% 
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "factor.site2high", "sympforeign site", "factor.site2high:sympforeign site"), c("Intercept", "Transplant site", "Own-foreign", "Transplant site * Own-foreign")))

# Testing symp for low and high sites
biomass3  <-  biomass %>% filter(factor.site2 == "high")
fit1 <- lmer((st.traitmean) ~ symp + (1|PopID)+(1|species)+(1|StudyID), biomass3, na.action = "na.fail", REML = FALSE) ### BEST MODEL
model.set <- dredge(fit1, rank = "AICc")
model.set
averaged.model <- model.avg(model.set)
data.frame(summary(averaged.model)$coefmat.full)


# 
BiomassSurvival  <-  LA %>% 
  filter(StudyID %in% c("Bast2015JPlE", "Rice1991Oeco", "Will1995Ecol")) %>% 
  filter(var != "fitness") %>% 
  filter(factor.site2 != "mid") %>% 
  select(StudyID, var, symp, factor.site2, st.traitmean)

ggplot(BiomassSurvival, aes(x = symp, y = st.traitmean, color = var, shape = factor.site2)) +
  geom_jitter() +
  facet_grid(~StudyID)


summary(fit1)
fix.check(fit1)
ggplot(biomass2, aes(x = symp, y = st.traitmean)) + geom_boxplot() + facet_wrap(~ factor.site2)



# FITNESS
fitness2  <-  fitness %>% filter(factor.site2 != "mid")
fit1 <- lmer((st.traitmean) ~ symp*factor.site2 + (1|PopID)+(1|species)+(1|StudyID), fitness2, na.action = "na.fail", REML = FALSE)
model.set <- dredge(fit1, rank = "AICc", extra = "R^2")
model.set
averaged.model <- model.avg(model.set)
res <- data.frame(summary(averaged.model)$coefmat.full)
FitnessHL <- res %>% 
  rownames_to_column(var = "Variable") %>% 
  setNames(., c("Variable", "Estimate", "StError", "AdjSE", "Zvalue", "Pvalue")) %>% 
  mutate(CI.low = Estimate - 1.96 * StError) %>% 
  mutate(CI.high = Estimate + 1.96 * StError) %>% 
  mutate(Estimate = round(Estimate, 2), CI = paste(round(CI.low, 2), round(CI.high, 2), sep = " - "), Zvalue = round(Zvalue, 2), Pvalue = round(Pvalue, 3)) %>% 
  select(Variable, Estimate, CI, Zvalue, Pvalue) %>% 
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "factor.site2high", "sympforeign site", "factor.site2high:sympforeign site"), c("Intercept", "Transplant site", "Own-foreign", "Transplant site * Own-foreign")))

ggplot(fitness2, aes(x = symp, y = st.traitmean)) + geom_boxplot() + facet_wrap(~ factor.site2)

summary(fit3)
fix.check(fit3)

fit1 <- lmer((st.traitmean) ~ symp*factor.site2 + (1|PopID)+(1|species)+(1|StudyID), fitness2, na.action = "na.fail", REML = FALSE)
fit2 <- lmer((st.traitmean) ~ factor.site2 + (1|PopID)+(1|species)+(1|StudyID), fitness2, na.action = "na.fail", REML = FALSE)
fit3 <- lmer((st.traitmean) ~ symp + (1|PopID)+(1|species)+(1|StudyID), fitness2, na.action = "na.fail", REML = FALSE)
fit4 <- lmer((st.traitmean) ~ 1 + (1|PopID)+(1|species)+(1|StudyID), fitness2, na.action = "na.fail", REML = FALSE)
modsel(list(fit1, fit2, fit3, fit4), 1000)

fit1 <- lmer((st.traitmean) ~ symp*factor.site2 + (1|PopID)+(1|species)+(1|StudyID), fitness2, na.action = "na.fail", REML = FALSE)
xx <- dredge(fit1, rank = "AICc")
mod.sel(fit1, fit2)
r.squaredLR(fit1)
r.squaredGLMM(fit1)

# SURVIVAL
survival2  <-  survival %>% filter(factor.site2 != "mid")
fit1 <- lmer((st.traitmean) ~ symp*factor.site2 + (1|PopID)+(1|species)+(1|StudyID), survival2, na.action = "na.fail", REML = FALSE)
model.set <- dredge(fit1, rank = "AICc")
model.set
averaged.model <- model.avg(model.set)
res <- data.frame(summary(averaged.model)$coefmat.full)
SurvivalHL <- res %>% 
  rownames_to_column(var = "Variable") %>% 
  setNames(., c("Variable", "Estimate", "StError", "AdjSE", "Zvalue", "Pvalue")) %>% 
  mutate(CI.low = Estimate - 1.96 * StError) %>% 
  mutate(CI.high = Estimate + 1.96 * StError) %>% 
  mutate(Estimate = round(Estimate, 2), CI = paste(round(CI.low, 2), round(CI.high, 2), sep = " - "), Zvalue = round(Zvalue, 2), Pvalue = round(Pvalue, 3)) %>% 
  select(Variable, Estimate, CI, Zvalue, Pvalue) %>% 
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "factor.site2high", "sympforeign site", "factor.site2high:sympforeign site"), c("Intercept", "Transplant site", "Own-foreign", "Transplant site * Own-foreign")))
ggplot(survival2, aes(x = symp, y = st.traitmean)) + geom_boxplot() + facet_wrap(~ factor.site2)


summary(fit3)
fix.check(fit3)


ResultsHL <- BiomassHL %>% 
  bind_rows(FitnessHL, SurvivalHL) %>% 
  mutate(Trait = c(rep("Biomass", 4), rep("Rep. output", 4), rep("Survival", 4))) %>% 
  select(Trait, Variable, Estimate, StError, Pvalue)





#### TEST THE EFFCET OF OTHER VARIABLES ####
# COUNT NUMBERS
LA %>% filter(var == "biomass", !is.na(breed), !is.na(generation)) %>% 
  select(studyunit, breed, longevity, growthform, intro) %>%
  gather(key = variable, value = value, -studyunit) %>% 
  arrange(variable) %>% 
  group_by(variable, value) %>% 
  summarise(n = n()) %>% print(n = 21)

LA %>% filter(var == "fitness", !is.na(breed), !is.na(generation)) %>% 
  select(studyunit, breed, longevity, growthform, intro) %>%
  gather(key = variable, value = value, -studyunit) %>% 
  arrange(variable) %>% 
  group_by(variable, value) %>% 
  summarise(n = n()) %>% print(n = 21)

LA %>% filter(var == "survival", !is.na(breed), !is.na(generation)) %>% 
  select(studyunit, breed, longevity, growthform, intro) %>%
  gather(key = variable, value = value, -studyunit) %>% 
  arrange(variable) %>% 
  group_by(variable, value) %>% 
  summarise(n = n()) %>% print(n = 21)


unique(LA$species)

# RESULTS
LA %>% 
  select(r.temp, dist.km, r.elev) %>% 
  gather(key = variable, value = value) %>% 
  group_by(variable) %>% 
  summarise(min = min (value, na.rm = TRUE), max = max (value, na.rm = TRUE), mean = mean(value, na.rm = TRUE), sd = sd(value, na.rm = TRUE))



LA %>% group_by(var) %>% summarise(n())



# BREEDING SYSTEM
Bio <- biomass %>% filter(!is.na(breed))
modBioBreed<- lmer(st.traitmean ~ symp*breed + (1|PopID)+(1|species)+(1|StudyID), Bio, na.action = "na.fail", REML = FALSE)
model.set1 <- dredge(modBioBreed, rank = "AICc")
averaged.model1 <- model.avg(model.set1)
res1 <- data.frame(summary(averaged.model1)$coefmat.full)
resBioBreed <- res1 %>% 
  rownames_to_column(var = "Variable") %>% 
  setNames(., c("Variable", "Estimate", "StError", "AdjSE", "Zvalue", "Pvalue")) %>% 
  mutate(CI.low = Estimate - 1.96 * StError) %>% 
  mutate(CI.high = Estimate + 1.96 * StError) %>% 
  mutate(Estimate = round(Estimate, 2), CI = paste(round(CI.low, 2), round(CI.high, 2), sep = " - "), Zvalue = round(Zvalue, 2), Pvalue = round(Pvalue, 3)) %>% 
  select(Variable, Estimate, CI, Zvalue, Pvalue) %>% 
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "sympforeign site", "breedmixed_mating", "breedoutcrossing", "breedmixed_mating:sympforeign site", "breedoutcrossing:sympforeign site"), c("Intercept", "OF", "BS - mixed mating", "BS - outcrossing", "OF*BS - mixed mating", "OF*BS - ourcrossing")))


Fit <- fitness %>% filter(!is.na(breed))
modFitBreed<- lmer(st.traitmean ~ symp*breed + (1|PopID)+(1|species)+(1|StudyID), Fit, na.action = "na.fail", REML = FALSE)
model.set2 <- dredge(modFitBreed, rank = "AICc")
averaged.model2 <- model.avg(model.set2)
res2 <- data.frame(summary(averaged.model2)$coefmat.full)
resFitBreed <- res2 %>% 
  rownames_to_column(var = "Variable") %>% 
  setNames(., c("Variable", "Estimate", "StError", "AdjSE", "Zvalue", "Pvalue")) %>% 
  mutate(CI.low = Estimate - 1.96 * StError) %>% 
  mutate(CI.high = Estimate + 1.96 * StError) %>% 
  mutate(Estimate = round(Estimate, 2), CI = paste(round(CI.low, 2), round(CI.high, 2), sep = " - "), Zvalue = round(Zvalue, 2), Pvalue = round(Pvalue, 3)) %>% 
  select(Variable, Estimate, CI, Zvalue, Pvalue) %>% 
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "sympforeign site", "breedmixed_mating", "breedoutcrossing", "breedmixed_mating:sympforeign site", "breedoutcrossing:sympforeign site"), c("Intercept", "OF", "BS - mixed mating", "BS - outcrossing", "OF*BS - mixed mating", "OF*BS - ourcrossing")))


Sur <- survival %>% filter(!is.na(breed))
modSurBreed<- lmer(st.traitmean ~ symp*breed + (1|PopID)+(1|species)+(1|StudyID), Sur, na.action = "na.fail", REML = FALSE)
model.set3 <- dredge(modSurBreed, rank = "AICc")
averaged.model3 <- model.avg(model.set3)
res3 <- data.frame(summary(averaged.model3)$coefmat.full)
resSurBreed <- res3 %>% 
  rownames_to_column(var = "Variable") %>% 
  setNames(., c("Variable", "Estimate", "StError", "AdjSE", "Zvalue", "Pvalue")) %>% 
  mutate(CI.low = Estimate - 1.96 * StError) %>% 
  mutate(CI.high = Estimate + 1.96 * StError) %>% 
  mutate(Estimate = round(Estimate, 2), CI = paste(round(CI.low, 2), round(CI.high, 2), sep = " - "), Zvalue = round(Zvalue, 2), Pvalue = round(Pvalue, 3)) %>% 
  select(Variable, Estimate, CI, Zvalue, Pvalue) %>% 
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "sympforeign site", "breedmixed_mating", "breedmixed_mating:sympforeign site"), c("Intercept", "OF", "BS - mixed mating", "OF*BS - mixed mating")))


ggplot(data=(dat3), aes(x=symp, y=st.traitmean, fill = breed)) + geom_boxplot()
dat3 %>% group_by(symp, breed) %>% summarise(mean = mean(st.traitmean)) %>% spread(key = symp, value = mean)


# GROWTHFORM
Bio <- biomass %>% filter(!is.na(growthform))
modBioGrowth<- lmer(st.traitmean ~ symp*growthform + (1|PopID)+(1|species)+(1|StudyID), Bio, na.action = "na.fail", REML = FALSE)
model.set4 <- dredge(modBioGrowth, rank = "AICc")
averaged.model4 <- model.avg(model.set4)
res4 <- data.frame(summary(averaged.model4)$coefmat.full)
resBioGrowth <- res4 %>% 
  rownames_to_column(var = "Variable") %>% 
  setNames(., c("Variable", "Estimate", "StError", "AdjSE", "Zvalue", "Pvalue")) %>% 
  mutate(CI.low = Estimate - 1.96 * StError) %>% 
  mutate(CI.high = Estimate + 1.96 * StError) %>% 
  mutate(Estimate = round(Estimate, 2), CI = paste(round(CI.low, 2), round(CI.high, 2), sep = " - "), Zvalue = round(Zvalue, 2), Pvalue = round(Pvalue, 3)) %>% 
  select(Variable, Estimate, CI, Zvalue, Pvalue) %>% 
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "sympforeign site", "growthformgrass", "growthformmoss", "growthformtree", "growthformgrass:sympforeign site", "growthformmoss:sympforeign site", "growthformtree:sympforeign site"), c("Intercept", "OF", "GF - grass", "GF - moss", "GF - tree", "OF*GF - grass", "OF*GF - moss", "OF*GF - tree")))


Fit <- fitness %>% filter(!is.na(growthform))
modFitGrowth<- lmer(st.traitmean ~ symp*growthform + (1|PopID)+(1|species)+(1|StudyID), Fit, na.action = "na.fail", REML = FALSE)
model.set5 <- dredge(modFitGrowth, rank = "AICc")
averaged.model5 <- model.avg(model.set5)
res5 <- data.frame(summary(averaged.model5)$coefmat.full)
resFitGrowth <- res5 %>% 
  rownames_to_column(var = "Variable") %>% 
  setNames(., c("Variable", "Estimate", "StError", "AdjSE", "Zvalue", "Pvalue")) %>% 
  mutate(CI.low = Estimate - 1.96 * StError) %>% 
  mutate(CI.high = Estimate + 1.96 * StError) %>% 
  mutate(Estimate = round(Estimate, 2), CI = paste(round(CI.low, 2), round(CI.high, 2), sep = " - "), Zvalue = round(Zvalue, 2), Pvalue = round(Pvalue, 3)) %>% 
  select(Variable, Estimate, CI, Zvalue, Pvalue) %>% 
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "sympforeign site", "growthformgrass", "growthformgrass:sympforeign site"), c("Intercept", "OF", "GF - grass", "OF*GF - grass")))


Sur <- survival %>% filter(!is.na(growthform))
modSurGrowth<- lmer(st.traitmean ~ symp*growthform*r.temp + (1|PopID)+(1|species)+(1|StudyID), Sur, na.action = "na.fail", REML = FALSE)
model.set5 <- dredge(modSurGrowth, rank = "AICc")
averaged.model5 <- model.avg(model.set5)
res5 <- data.frame(summary(averaged.model5)$coefmat.full)
resSurGrowth <- res5 %>% 
  rownames_to_column(var = "Variable") %>% 
  setNames(., c("Variable", "Estimate", "StError", "AdjSE", "Zvalue", "Pvalue")) %>% 
  mutate(CI.low = Estimate - 1.96 * StError) %>% 
  mutate(CI.high = Estimate + 1.96 * StError) %>% 
  mutate(Estimate = round(Estimate, 2), CI = paste(round(CI.low, 2), round(CI.high, 2), sep = " - "), Zvalue = round(Zvalue, 2), Pvalue = round(Pvalue, 3)) %>% 
  select(Variable, Estimate, CI, Zvalue, Pvalue) %>% 
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "sympforeign site", "growthformforb", "growthformgrass", "growthformforb:sympforeign site", "growthformgrass:sympforeign site"), c("Intercept", "OF", "GF - forb", "GF - grass", "OF*GF - forb", "OF*GF - grass")))


# LONGEVITY
Bio <- biomass %>% filter(!is.na(longevity))
modBioLong<- lmer(st.traitmean ~ symp*longevity + (1|PopID)+(1|species)+(1|StudyID), Bio, na.action = "na.fail", REML = FALSE)
model.set6 <- dredge(modBioLong, rank = "AICc")
averaged.model6 <- model.avg(model.set6)
res6 <- data.frame(summary(averaged.model6)$coefmat.full)
resBioLong <- res6 %>% 
  rownames_to_column(var = "Variable") %>% 
  setNames(., c("Variable", "Estimate", "StError", "AdjSE", "Zvalue", "Pvalue")) %>% 
  mutate(CI.low = Estimate - 1.96 * StError) %>% 
  mutate(CI.high = Estimate + 1.96 * StError) %>% 
  mutate(Estimate = round(Estimate, 2), CI = paste(round(CI.low, 2), round(CI.high, 2), sep = " - "), Zvalue = round(Zvalue, 2), Pvalue = round(Pvalue, 3)) %>% 
  select(Variable, Estimate, CI, Zvalue, Pvalue) %>% 
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "sympforeign site", "longevityperennial", "longevityperennial:sympforeign site"), c("Intercept", "OF", "LO - perennial", "OF*LO - perennial")))


Fit <- fitness %>% filter(!is.na(longevity))
modFitLong<- lmer(st.traitmean ~ symp*longevity + (1|PopID)+(1|species)+(1|StudyID), Fit, na.action = "na.fail", REML = FALSE)
model.set7 <- dredge(modFitLong, rank = "AICc")
averaged.model7 <- model.avg(model.set7)
res7 <- data.frame(summary(averaged.model7)$coefmat.full)
resFitLong <- res7 %>% 
  rownames_to_column(var = "Variable") %>% 
  setNames(., c("Variable", "Estimate", "StError", "AdjSE", "Zvalue", "Pvalue")) %>% 
  mutate(CI.low = Estimate - 1.96 * StError) %>% 
  mutate(CI.high = Estimate + 1.96 * StError) %>% 
  mutate(Estimate = round(Estimate, 2), CI = paste(round(CI.low, 2), round(CI.high, 2), sep = " - "), Zvalue = round(Zvalue, 2), Pvalue = round(Pvalue, 3)) %>% 
  select(Variable, Estimate, CI, Zvalue, Pvalue) %>% 
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "sympforeign site", "longevityperennial", "longevityperennial:sympforeign site"), c("Intercept", "OF", "LO - perennial", "OF*LO - perennial")))


Sur <- survival %>% filter(!is.na(longevity))
modSurLong<- lmer(st.traitmean ~ symp*longevity + (1|PopID)+(1|species)+(1|StudyID), Sur, na.action = "na.fail", REML = FALSE)
model.set8 <- dredge(modSurLong, rank = "AICc")
averaged.model8 <- model.avg(model.set8)
res8 <- data.frame(summary(averaged.model8)$coefmat.full)
resSurLong <- res8 %>% 
  rownames_to_column(var = "Variable") %>% 
  setNames(., c("Variable", "Estimate", "StError", "AdjSE", "Zvalue", "Pvalue")) %>% 
  mutate(CI.low = Estimate - 1.96 * StError) %>% 
  mutate(CI.high = Estimate + 1.96 * StError) %>% 
  mutate(Estimate = round(Estimate, 2), CI = paste(round(CI.low, 2), round(CI.high, 2), sep = " - "), Zvalue = round(Zvalue, 2), Pvalue = round(Pvalue, 3)) %>% 
  select(Variable, Estimate, CI, Zvalue, Pvalue) %>% 
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "sympforeign site", "longevityperennial", "longevityperennial:sympforeign site"), c("Intercept", "OF", "LO - perennial", "OF*LO - perennial")))



resBioBreed %>% 
  bind_rows(resFitBreed, resSurBreed, resBioGrowth, resFitGrowth, resSurGrowth, resBioLong, resFitLong, resSurLong) %>% 
  mutate(Trait = c(rep("Biomass", 6), rep("Fitness", 6), rep("Survival", 4), rep("Biomass", 8), rep("Fitness", 4), rep("Survival", 6), rep("Biomass", 4), rep("Fitness", 4), rep("Survival", 4))) %>% 
  mutate(LifeHistory = c(rep("Breeding system", 16), rep("Growthform", 18), rep("Longevity", 12))) %>% 
  select(Trait, LifeHistory, Variable, Estimate, CI, Zvalue, Pvalue) %>% 
  mutate(Star = ifelse(Pvalue <= 0.05, "***", Pvalue))



####################################################
#### Fig. 4 ADAPTATION ####
####################################################

library("cowplot")

# CLIMATE DISTANCE
# define range
yylim <- c(min(d6$std.traitmean), max(d5$std.traitmean))

p1 <- ggplot(d4, aes(x = clim.diff, y = std.traitmean)) +
  geom_vline(xintercept = 0, colour = "grey", linetype = "longdash") +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, colour = "steelblue", linetype = "solid", size = 1) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), colour = "steelblue", linetype = "longdash", size = 1) +
  labs(x = "", y = "", title = "Biomass") +
  ylim(yylim[1], yylim[2]) +
  annotate("text", x = -8, y = 4, label = "C", size = 5) +
  panel_border(colour = "black", remove = FALSE) +
  theme(text = element_text(size = 12), plot.title = element_text(size=12))

p2 <- ggplot(d5, aes(x = clim.diff, y = std.traitmean)) +
  geom_vline(xintercept = 0, colour = "grey", linetype = "longdash") +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, colour = "steelblue", linetype = "solid", size = 1) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), colour = "steelblue", linetype = "longdash", size = 1) +
  labs(x = "", y = "", title = "Reproductive output") +
  ylim(yylim[1], yylim[2]) +
  annotate("text", x = -8, y = 4, label = "B", size = 5) +
  panel_border(colour = "black", remove = FALSE) +
  theme(text = element_text(size = 12), plot.title = element_text(size=12))

p3 <- ggplot(d6, aes(x = clim.diff, y = std.traitmean)) +
  geom_vline(xintercept = 0, colour = "grey", linetype = "longdash") +
  geom_point() +
  stat_smooth(method = "lm", formula = y ~ x, colour = "steelblue", linetype = "solid", size = 1) +
  stat_smooth(method = "lm", formula = y ~ x + I(x^2), colour = "steelblue", linetype = "longdash", size = 1) +
  labs(x = "", y = "", title = "Survival") +
  ylim(yylim[1], yylim[2]) +
  annotate("text", x = -5, y = 3, label = "upward", size= 4.5) +
  annotate("text", x = 10, y = 3, label = "downward transplant", size= 4.5) +
  annotate("text", x = -8, y = 4, label = "A", size = 5) +
  panel_border(colour = "black", remove = FALSE) +
  theme(text = element_text(size = 12), plot.title = element_text(size=12))


# LOCAL ADAPTATION
# define range
yylim <- c(min(biomass$st.traitmean, na.rm = TRUE), max(biomass$st.traitmean, na.rm = TRUE))

p4 <- LA %>%
  filter(var == "biomass") %>% 
  group_by(var, symp) %>% 
  summarise(n = n(), mean = mean(st.traitmean), se = sd(st.traitmean)/sqrt(n), CI.low = mean - 1.96*se, CI.high = mean + 1.96*se) %>% 
  mutate(symp = plyr::mapvalues(symp, c("home site", "foreign site"), c("own", "foreign"))) %>%
  ggplot(aes(x = symp, y = mean, ymin = CI.low, ymax = CI.high)) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "longdash") +
  geom_point(size = 3) +
  geom_errorbar(width = 0) +
  labs(x = "", y = "", title = "Biomass") +
  annotate("text", x = 0.6, y = 0.7, label = "C", size = 5) +
  panel_border(colour = "black", remove = FALSE) +
  ylim(-0.7, 0.8)

p5 <- LA %>%
  filter(var == "fitness") %>% 
  group_by(var, symp) %>% 
  summarise(n = n(), mean = mean(st.traitmean), se = sd(st.traitmean)/sqrt(n), CI.low = mean - 1.96*se, CI.high = mean + 1.96*se) %>% 
  mutate(symp = plyr::mapvalues(symp, c("home site", "foreign site"), c("own", "foreign"))) %>%
  ggplot(aes(x = symp, y = mean, ymin = CI.low, ymax = CI.high)) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "longdash") +
  geom_point(size = 3) +
  geom_errorbar(width = 0) +
  labs(x = "", y = "", title = "Reproductive output") +
  annotate("text", x = 0.6, y = 0.7, label = "B", size = 5) +
  panel_border(colour = "black", remove = FALSE) +
  ylim(-0.7, 0.8)

p6 <- LA %>%
  filter(var == "survival") %>% 
  group_by(var, symp) %>% 
  summarise(n = n(), mean = mean(st.traitmean), se = sd(st.traitmean)/sqrt(n), CI.low = mean - 1.96*se, CI.high = mean + 1.96*se) %>% 
  mutate(symp = plyr::mapvalues(symp, c("home site", "foreign site"), c("own", "foreign"))) %>%
  ggplot(aes(x = symp, y = mean, ymin = CI.low, ymax = CI.high)) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "longdash") +
  geom_point(size = 3) +
  geom_errorbar(width = 0) +
  labs(x = "", y = "", title = "Survival") +
  annotate("text", x = 0.6, y = 0.7, label = "A", size = 5) +
  ylim(-0.7, 0.8) +
  panel_border(colour = "black", remove = FALSE) +
  annotate("text", x = 2.1, y = 0.7, label = "OF*", size= 5)


p4.old <- ggplot(biomass, aes(x = symp, y = st.traitmean)) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "longdash") +
  geom_point() +
  labs(x = "", y = "St. traitmean") +
  scale_x_discrete(labels=c("home site" = "own", "foreign site" = "foreign")) +
  ylim(yylim[1], yylim[2]) 

p5.old <- ggplot(fitness, aes(x = symp, y = st.traitmean)) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "longdash") +
  geom_boxplot() +
  labs(x = "", y = "") +
  ylim(yylim[1], yylim[2]) +
  scale_x_discrete(labels=c("home site" = "own", "foreign site" = "foreign"))

p6.old <- ggplot(survival, aes(x = symp, y = st.traitmean)) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "longdash") +
  geom_boxplot() +
  labs(x = "", y = "") +
  ylim(yylim[1], yylim[2]) +
  scale_x_discrete(labels=c("home site" = "own", "foreign site" = "foreign")) +
  annotate("text", x = 1.9, y = yylim[2]*0.9, label = "*", size= 6)


# SYMP * SITE

p7 <- biomass %>%
  filter(factor.site2 !="mid") %>% 
  mutate(symp.site = paste(symp, factor.site2, sep = "_")) %>% 
  mutate(symp = plyr::mapvalues(symp, c("home site", "foreign site"), c("own", "foreign"))) %>%
  mutate(factor.site2 = factor(factor.site2, levels = c("low", "high"))) %>% 
  group_by(var, symp.site) %>% 
  summarise(n = n(), mean = mean(st.traitmean), se = sd(st.traitmean)/sqrt(n), CI.low = mean - 1.96*se, CI.high = mean + 1.96*se) %>% 
  #mutate(symp.site.nr = c(2.1, 1.9, 1.1, 0.9)) %>% 
  mutate(symp.site.nr = c(1.1, 1.9, 2.1, 0.9)) %>% 
  mutate(elevation = c("high elevation", "low elevation", "high elevation", "low elevation")) %>%
  ggplot(aes(x = symp.site.nr, y = mean, ymin = CI.low, ymax = CI.high, shape = elevation, linetype = elevation)) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "longdash") +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(width = 0, linetype = "solid") +
  labs(x = "", y = "") +
  #scale_x_continuous(breaks = c(1,2), labels = c("own", "foreign")) +
  scale_x_continuous(breaks = c(1,2), labels = c("low", "high")) +
  scale_shape_manual(name = "Elevation of origin", values = c(17, 16)) +
  scale_linetype_manual(name = "Elevation of origin", values = c("dashed", "solid")) +
  annotate("text", x = 2, y = 0.9, label = "OF***", size= 5) +
  annotate("text", x = 2, y = 0.7, label = "TS***", size= 5) +
  annotate("text", x = 1.85, y = 0.5, label = "OF x TS***", size= 5) +
  annotate("text", x = 0.95, y = 0.9, label = "F", size = 5) +
  ylim(-1, 1) +
  panel_border(colour = "black", remove = FALSE) +
  theme(legend.position = c(0.35, 0.2))


p8 <- fitness %>%
  filter(factor.site2 !="mid") %>% 
  mutate(symp.site = paste(symp, factor.site2, sep = "_")) %>% 
  mutate(symp = plyr::mapvalues(symp, c("home site", "foreign site"), c("own", "foreign"))) %>%
  mutate(factor.site2 = factor(factor.site2, levels = c("low", "high"))) %>% 
  group_by(var, symp.site) %>% 
  summarise(n = n(), mean = mean(st.traitmean), se = sd(st.traitmean)/sqrt(n), CI.low = mean - 1.96*se, CI.high = mean + 1.96*se) %>% 
  mutate(symp.site.nr = c(1.1, 1.9, 2.1, 0.9)) %>% 
  mutate(elevation = c("high elevation", "low elevation", "high elevation", "low elevation")) %>%
  ggplot(aes(x = symp.site.nr, y = mean, ymin = CI.low, ymax = CI.high, shape = elevation, linetype = elevation)) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "longdash") +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(width = 0, linetype = "solid") +
  labs(x = "", y = "") +
  scale_x_continuous(breaks = c(1,2), labels = c("low", "high")) +
  scale_shape_manual(values = c(17, 16)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  ylim(-1, 1) +
  annotate("text", x = 0.95, y = 0.9, label = "E", size = 5) +
  panel_border(colour = "black", remove = FALSE) +
  theme(legend.position = "none")


p9 <- survival %>%
  filter(factor.site2 !="mid") %>% 
  mutate(symp.site = paste(symp, factor.site2, sep = "_")) %>% 
  mutate(symp = plyr::mapvalues(symp, c("home site", "foreign site"), c("own", "foreign"))) %>%
  mutate(factor.site2 = factor(factor.site2, levels = c("low", "high"))) %>% 
  group_by(var, symp.site) %>% 
  summarise(n = n(), mean = mean(st.traitmean), se = sd(st.traitmean)/sqrt(n), CI.low = mean - 1.96*se, CI.high = mean + 1.96*se) %>% 
  mutate(symp.site.nr = c(1.1, 1.9, 2.1, 0.9)) %>% 
  mutate(elevation = c("high elevation", "low elevation", "high elevation", "low elevation")) %>%
  ggplot(aes(x = symp.site.nr, y = mean, ymin = CI.low, ymax = CI.high, shape = elevation, linetype = elevation)) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "longdash") +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(width = 0, linetype = "solid") +
  labs(x = "", y = "") +
  scale_x_continuous(breaks = c(1,2), labels = c("low", "high")) +
  scale_shape_manual(values = c(17, 16)) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  annotate("text", x = 1.9, y = 0.9, label = "OF*", size= 5) +
  annotate("text", x = 0.95, y = 0.9, label = "D", size = 5) +
  ylim(-1, 1) +
  panel_border(colour = "black", remove = FALSE) +
  theme(legend.position = "none")


p7.old <- biomass %>% filter(factor.site2 !="mid") %>% 
  mutate(factor.site2 = factor(factor.site2, levels = c("low", "high"))) %>% 
  ggplot(aes(x = symp, y = st.traitmean)) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "longdash") +
  geom_boxplot() +
  labs(x = "", y = "") +
  ylim(yylim[1], yylim[2]) +
  facet_wrap(~ factor.site2) +
  scale_x_discrete(labels=c("home site" = "own", "foreign site" = "foreign", "home site" = "own", "foreign site" = "foreign")) +
  annotate("text", x = 1.8, y = yylim[2]*0.9, label = "*", size= 6)

p8.old <- fitness %>% filter(factor.site2 !="mid") %>% 
  mutate(factor.site2 = factor(factor.site2, levels = c("low", "high"))) %>% 
  ggplot(aes(x = symp, y = st.traitmean)) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "longdash") +
  geom_boxplot() +
  labs(x = "", y = "") +
  ylim(yylim[1], yylim[2]) +
  facet_wrap(~ factor.site2) +
  scale_x_discrete(labels=c("home site" = "own", "foreign site" = "foreign", "home site" = "own", "foreign site" = "foreign"))

p9.old <- survival %>% filter(factor.site2 !="mid") %>% 
  mutate(factor.site2 = factor(factor.site2, levels = c("low", "high"))) %>% 
  ggplot(aes(x = symp, y = st.traitmean)) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "longdash") +
  geom_boxplot() +
  labs(x = "", y = "") +
  ylim(yylim[1], yylim[2]) +
  facet_wrap(~ factor.site2) +
  scale_x_discrete(labels=c("home site" = "own", "foreign site" = "foreign", "home site" = "own", "foreign site" = "foreign")) +
  annotate("text", x = 1.8, y = yylim[2]*0.9, label = "*", size= 7)

library("cowplot")
library("gridExtra")
library(grid)

AdaptationPlot <- grid.arrange(p6, p5, p4, p9, p8, p7, layout_matrix = rbind(c(1,2,3),c(4,5,6)), bottom = textGrob("Elevation of transplant site", vjust = 0.1, gp = gpar(fontsize = 15, font = 8)), left = textGrob("Standardized trait mean", rot = 90, vjust = 1, gp = gpar(fontsize = 15, font = 8)))
ggsave("Fig4_Adaptation4.pdf", AdaptationPlot, height = 6, width = 10)


AdaptationAppendixPlot <- grid.arrange(p3, p2, p1, layout_matrix = rbind(c(1,2,3)), bottom = textGrob("Climate difference", vjust = 0.1, gp = gpar(fontsize = 15, font = 8)), left = textGrob("Standardized trait mean", rot = 90, vjust = 1, gp = gpar(fontsize = 15, font = 8)))
ggsave("Fig_AdaptationAppendix.pdf", AdaptationAppendixPlot, height = 6)



#### Correlation between Biomass and Fitness ####

biomass <- subset(LA2, LA2$TRAIT_CAT1 == "SIZE"& LA2$TRAIT_CAT2 == "BIOMASS"& LA2$TRAIT_CAT3 == "ABOVEGR")
fitness <- subset(LA2, LA2$TRAIT_CAT1 == "FITNESS")
#fitness2 <- subset(LA2, LA2$TRAIT_CAT1 == "FITNESS" & LA2$TRAIT_CAT2=="CUMFIT")
#fitness <- rbind(fitness, fitness2)
survival <- subset(LA2, LA2$TRAIT_CAT1 == "VITAL_RATES"& LA2$TRAIT_CAT2 == "SURVIVAL")
biomass$StudyID
fitness$StudyID

biofit <- biomass[biomass$ID %in% fitness$ID,]
biofit$st.fitness <- biofit[match(biofit$ID, fitness$ID), "st.traitmean"]
plot(biofit$st.traitmean, biofit$st.fitness)
fit <- lm(st.fitness~st.traitmean, biofit)
summary(fit)
abline(fit)

biosur <- biomass[biomass$ID %in% survival$ID,]
biosur$st.survival <- biosur[match(biosur$ID, survival$ID), "st.traitmean"]
plot(biosur$st.traitmean, biosur$st.fitness)
fit <- lm(st.survival~st.traitmean, biosur)
summary(fit)
abline(fit)


### Count how many cases home pops have higher value
biomass %>% group_by(StudyID, Trait, species, symp) %>% summarise(mean = mean(traitmean)) %>% spread(key = symp, value = mean) %>% mutate(diff = `home site` - `foreign site`, adapt = ifelse(diff > 0, 1, 0)) %>% ungroup() %>% summarise(sum = sum(adapt))
# 18 * 100 / 28 = 64%

fitness %>% group_by(StudyID, Trait, species, symp) %>% summarise(mean = mean(traitmean)) %>% spread(key = symp, value = mean) %>% mutate(diff = `home site` - `foreign site`, adapt = ifelse(diff > 0, 1, 0)) %>% ungroup() %>% summarise(sum = sum(adapt))
# 22 * 100 / 33 = 67%

survival %>% group_by(StudyID, Trait, species, symp) %>% summarise(mean = mean(traitmean)) %>% spread(key = symp, value = mean) %>% mutate(diff = `home site` - `foreign site`, adapt = ifelse(diff > 0, 1, 0)) %>% ungroup() %>% summarise(sum = sum(adapt))
# 9 * 100 / 9 = 100%

