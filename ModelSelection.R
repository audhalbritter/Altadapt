##############################################################################################
###### MODEL SELECTION WITH MODEL AVERAGING ######
##############################################################################################

library("lme4")
library("MuMIn")
library("tibble")
library("cowplot")

#  change the default "na.omit" to prevent models from being fitted to different datasets in case of missing values.
options(na.action = "na.fail") # can also be put in the model
options(na.action = "na.omit") # change back
# Alternatively put it in the modle


# subsetting data for height, biomass and phenology, create a variable var
# create variables for absolute slope and latitude
dat2 <- read.csv(file = "FinalSharedData/TraitDifferentiation.csv")

dat4 <- dat2 %>% 
  filter(TRAIT_CAT1 == "SIZE" & TRAIT_CAT2 == "HEIGHT" | TRAIT_CAT1 == "SIZE" & TRAIT_CAT2 == "BIOMASS" & TRAIT_CAT3 == "ABOVEGR" | TRAIT_CAT1 == "PHENOLOGY" & TRAIT_CAT2 %in% c("FLOWERING", "LEAF_BUD") & TRAIT_CAT4 == c("first", "peak")) %>%
  mutate(var = ifelse(TRAIT_CAT1 == "SIZE" & TRAIT_CAT2 == "HEIGHT", "height",
                      ifelse(TRAIT_CAT1 == "SIZE" & TRAIT_CAT2 == "BIOMASS" & TRAIT_CAT3 == "ABOVEGR", "biomass", "phenology"))) %>% 
  rename(slopes = estimate, sample.size = n) %>% 
  mutate(abs.slopes = abs(slopes)) %>% 
  mutate(abs.lat = abs(mean.lat))


#######################################################
#### HEIGHT ####
#######################################################

hh <- dat4 %>% filter(var == "height") %>% 
  filter(!is.na(breed)) %>%  # remove NA's
  mutate(Year = factor(Year))


# UNTRANSFORMED SLOPES
# Define Model
modHeight_UT <- lmer(slopes ~ r.temp + abs.lat + dist.km + breed + growthform + intro + (1|species)+(1|StudyID)+(1|Year), data = hh, REML = FALSE, weights = sample.size, na.action = "na.fail")

# check model assumptions
fix.check(modHeight_UT)


# Model selection and averageing mixed effects models, produce table with estimate, Pvalue, Importance
ResHeight_UT <- ModelAverage(modHeight_UT, printFullTable = TRUE, print95Table = TRUE, percent.thresh = 0.95, "Height_UT")


# Model selection and averageing mixed effects models, produce table with 95% weighted models
ResH_UT <- ModelSelectionTable(modHeight_UT, 0.95, "Height_UT")


# Produce plot with estimates
PlotHeight_UT <- PlotEstimates3(modHeight_UT, percent.thresh = 0.95)


# ABSOLUTE SLOPES
# Define Model
modHeight_ABS <- lmer(abs.slopes ~ r.temp + abs.lat + dist.km + breed + growthform + intro + (1|species)+(1|StudyID)+(1|Year), data = hh, REML = FALSE, weights = sample.size, na.action = "na.fail")

# check model assumptions
fix.check(modHeight_ABS)


# Model selection and averageing mixed effects models, produce table with estimate, Pvalue, Importance
ResHeight_ABS <- ModelAverage(modHeight_ABS, printFullTable = TRUE, print95Table = TRUE, percent.thresh = 0.95, "Height_ABS")


# Model selection and averageing mixed effects models, produce table with 95% weighted models
ResH_ABS <- ModelSelectionTable(modHeight_ABS, 0.95, "Height_ABS")


# Produce plot with estimates
PlotHeight_ABS <- PlotEstimates3(modHeight_ABS, percent.thresh = 0.95)


#######################################################
#### BIOMASS ####
#######################################################

bb <- dat4 %>% filter(var == "biomass") %>% 
  filter(!is.na(breed)) %>%  # remove NA's
  mutate(Year = factor(Year))


# UNTRANSFORMED SLOPES
# Define model
modBiomass_UT <- lmer(slopes ~ r.temp + abs.lat + dist.km + breed + growthform + intro + (1|species)+(1|StudyID)+(1|Year), data = bb, REML = FALSE, weights = sample.size, na.action = "na.fail")

# check model assumptions
fix.check(modBiomass_UT)

# Model selection and averageing mixed effects models, produce table with estimate, Pvalue, Importance
ResBiomass_UT <- ModelAverage(modBiomass_UT, printFullTable = TRUE, print95Table = TRUE, percent.thresh = 0.95, "Biomass_UT")


# Model selection and averageing mixed effects models, produce table with 95% weighted models
ResB_UT <- ModelSelectionTable(modBiomass_UT, 0.95, "Biomass_UT")


# Produce plot with estimates
PlotBiomass_UT <- PlotEstimates3(modBiomass_UT, percent.thresh = 0.95)


# ABSOLUTE SLOPES
# Define model
modBiomass_ABS <- lmer(abs.slopes ~ r.temp + abs.lat + dist.km + breed + growthform + intro + (1|species)+(1|StudyID)+(1|Year), data = bb, REML = FALSE, weights = sample.size, na.action = "na.fail")

# check model assumptions
fix.check(modBiomass_ABS)


# Model selection and averageing mixed effects models, produce table with estimate, Pvalue, Importance
ResBiomass_ABS <- ModelAverage(modBiomass_ABS, printFullTable = TRUE, print95Table = TRUE, percent.thresh = 0.95, "Biomass_ABS")


# Model selection and averageing mixed effects models, produce table with 95% weighted models
ResB_ABS <- ModelSelectionTable(modBiomass_ABS, 0.95, "Biomass_ABS")


# Produce plot with estimates
PlotBiomass_ABS <- PlotEstimates3(modBiomass_ABS, percent.thresh = 0.95)




#######################################################
#### PHENOLOGY ####
#######################################################

pe <- dat4 %>% filter(var == "phenology") %>% 
  filter(!is.na(breed)) %>%  # remove NA's
  mutate(Year = factor(Year))


# UNTRANSFORMED SLOPES
# Define model
modPhenologyEmergence_UT <- lmer(slopes ~ r.temp + abs.lat + dist.km + breed + growthform + intro + (1|species)+(1|StudyID)+(1|Year), data = pe, REML = FALSE, weights = sample.size, na.action = "na.fail")

# check model assumptions
fix.check(modPhenologyEmergence_UT)


# Model selection and averageing mixed effects models, produce table with estimate, Pvalue, Importance
ResPhenologyEmergence_UT <- ModelAverage(modPhenologyEmergence_UT, printFullTable = TRUE, print95Table = TRUE, percent.thresh = 0.95, "Phenology_UT")


# Model selection and averageing mixed effects models, produce table with 95% weighted models
ResP_UT <- ModelSelectionTable(modPhenologyEmergence_UT, 0.95, "Phenology_UT")


# Produce plot with estimates
PlotPhenology_UT <- PlotEstimates3(modPhenologyEmergence_UT, percent.thresh = 0.95)



# ABSOLUTE SLOPES
# Define model
modPhenologyEmergence_ABS <- lmer(abs.slopes ~ r.temp + abs.lat + dist.km + breed + growthform + intro + (1|species)+(1|StudyID)+(1|Year), data = pe, REML = FALSE, weights = sample.size, na.action = "na.fail")

# check model assumptions
fix.check(modPhenologyEmergence_ABS)


# Model selection and averageing mixed effects models, produce table with estimate, Pvalue, Importance
ResPhenologyEmergence_ABS <- ModelAverage(modPhenologyEmergence_ABS, printFullTable = TRUE, print95Table = TRUE, percent.thresh = 0.95, "Phenology_ABS")


# Model selection and averageing mixed effects models, produce table with 95% weighted models
ResP_ABS <- ModelSelectionTable(modPhenologyEmergence_ABS, 0.95, "Phenology_ABS")


# Produce plot with estimates
PlotPhenology_ABS <- PlotEstimates3(modPhenologyEmergence_ABS, percent.thresh = 0.95)






#######################################################
#### RESULT TABLE ####
#######################################################

# Main
ResultEstimate <- ResHeight_UT %>% 
  bind_rows(ResHeight_ABS, ResBiomass_UT, ResBiomass_ABS, ResPhenologyEmergence_UT, ResPhenologyEmergence_ABS) %>% 
  select(Trait, Variable, Estimate, CI, Zvalue, Pvalue, Importance) %>% 
  mutate(Variable = plyr:: mapvalues(Variable, c("(Intercept)", "T Range", "Abs. Latitude", "dist.km", "Native", "Grass", "Tree", "Mixed Mating", "Outcrossing"), c("Intercept", "T Range", "Abs. Latitude", "Distance", "Introduction status:Native", "Growthform:Grass", "Growthform:Tree", "Breeding System:Mixed Mating", "Breeding Sytem:Outcrossing"))) %>% 
  mutate(Variable = factor(Variable, levels = c("Intercept", "T Range", "Abs. Latitude", "Distance", "Introduction status:Native", "Growthform:Grass", "Growthform:Tree", "Breeding System:Mixed Mating", "Breeding Sytem:Outcrossing"))) %>% 
  mutate(Trait = factor(Trait, levels = c("Height_UT", "Height_ABS", "Biomass_UT", "Biomass_ABS","Phenology_UT", "Phenology_ABS"))) %>% 
  arrange(Trait, Variable) %>% 
  mutate(CI = paste("(", CI, ")", sep = ""))




#######################################################
#### MODEL TABLE S4 ####
#######################################################

ResModelSelection <- ResH_UT %>% 
  bind_rows(ResH_ABS, ResB_UT, ResB_ABS, ResP_UT, ResP_ABS) %>% 
  select(Trait, Model, df, delta, R.square, weight)

# count models
ResModelSelection %>%   
  group_by(Trait) %>% 
  summarise(n())



#######################################################
#### PLOT ESTIMATE ####
#######################################################
# Plot Estimates

ParameterEstimate <- ResHeight_UT %>% 
  bind_rows(ResHeight_ABS, ResBiomass_UT, ResBiomass_ABS, ResPhenologyEmergence_UT, ResPhenologyEmergence_ABS)

ParameterEstimatePlot <- ParameterEstimate %>% 
  filter(Variable != c("(Intercept)")) %>% 
  filter(Variable != c("T Range")) %>% 
  separate(col = CI, into = c("CI.low", "CI.high"), sep = " - ", convert = TRUE) %>% 
  mutate(Variable = plyr::mapvalues(Variable, c("Outcrossing", "Mixed Mating", "Tree", "Grass", "Native", "Abs. Latitude", "Distance"), c("Outcrossing", "Mixed Mating", "Tree", "Grass", "Native", "Abs. Latitude", "Max Distance"))) %>% 
  mutate(Variable = factor(Variable, levels = c("Outcrossing", "Mixed Mating", "Tree", "Grass", "Native", "Max Distance", "Abs. Latitude"))) %>% 
  separate(col = Trait, into = c("Trait", "Slope"), sep = "_") %>% 
  mutate(Slope = plyr::mapvalues(Slope, c("ABS", "UT"), c("Absolute slope", "Slope"))) %>% 
  mutate(Slope = factor(Slope, levels = c("Slope", "Absolute slope"))) %>% 
  mutate(Trait = factor(Trait, levels = c("Height", "Biomass", "Phenology"))) %>% 
  ggplot(aes(y = Variable, x = Estimate, xmin = CI.low, xmax = CI.high)) + 
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
  geom_point() +
  labs(x = "Parameter estimate", y = "") +
  geom_errorbarh(height = 0) +
  panel_border(colour = "black", remove = FALSE) +
  facet_grid(Slope ~ Trait, scales = "free_x") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))

ggsave("FinalFigures/Fig3ParameterEstimatePlot.pdf", ParameterEstimatePlot, height = 7, width = 10, dpi = 300)


summaryMean <- dat4 %>% 
  group_by(var) %>% 
  summarise(n = n(), mean = mean(abs.slopes), se = sd(abs.slopes)/sqrt(n), CI.low = mean - 1.96 * se, CI.high = mean + 1.96 * se)






# Mean and 95% CI
Fig2EffectSizePlot <- dat4 %>% 
  select(var, slopes, p.value) %>% 
  filter(var != "phenology duration") %>% 
  setNames(., c("Trait", "Slope", "Pvalue")) %>% 
  mutate(Trait = plyr::mapvalues(Trait, c("height", "biomass", "phenology"), c("Height", "Biomass", "Phenology"))) %>%
  mutate(Trait = factor(Trait, level = c("Height", "Biomass", "Phenology"))) %>%
  mutate(Pvalue = ifelse(Pvalue < 0.05, "sign.", "non sign.")) %>% 
  #gather(key = Slope, value = Value, -Trait, -Pvalue) %>% 
  #mutate(Slope = factor(Slope, levels = c("Untransformed slope", "Absolute slope"))) %>% 
  ggplot(aes(x = Trait, y = Slope)) + 
  geom_violin(draw_quantiles = c(0.5), position = position_dodge(width = 0.60), fill = "grey80") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_jitter(width = 0.2, aes(shape = Pvalue)) +
  labs(x = "", y = "Temperature effect (regression slope)") +
  scale_shape_manual(values = c(1,17)) +
  stat_summary(fun.data = "mean_cl_normal", colour = "black", alpha = 0.8, size = 1.2, shape = 16, position = position_dodge(width = 0.60)) +
  theme(legend.position = "none")

ggsave("FinalFigures/Fig2EffectSizePlot.pdf", Fig2EffectSizePlot, width = 8, height = 6, dpi = 300)


# Calculate mean, se and CI for slopes
dat4 %>% 
  select(var, slopes, abs.slopes) %>% 
  gather(key = slope, value = value, -var) %>% 
  group_by(var, slope) %>% 
  summarise(n = n(), mean = mean(value), se = sd(value)/sqrt(n), min(value), max(value)) %>% 
  mutate(ci.low = mean - 1.96 * se, ci.high = mean + 1.96 * se)


dd <- dat4 %>% 
  select(var, slopes, abs.slopes, StudyID, species, Year) %>% 
  gather(key = slope, value = value, slopes, abs.slopes) %>% 
  group_by(var, slope) %>% 
  do(fit = lmer(value ~ 1 + (1|species)+(1|StudyID)+(1|Year), data = .))

tidy(dd, fit) %>% 
  filter(!is.na(std.error)) %>% 
  mutate(ci.low = estimate - 1.96 * std.error, ci.high = estimate + 1.96 * std.error)

