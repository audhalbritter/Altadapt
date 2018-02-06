##########################################################################################################
###### MODEL SELECTION WITH MODEL AVERAGING ######
##########################################################################################################

library("lme4")
library("MuMIn")
library("tibble")
library("cowplot")

#  change the default "na.omit" to prevent models from being fitted to different datasets in case of missing values.
options(na.action = "na.fail") # can also be put in the model
options(na.action = "na.omit") # change back
# Alternatively put it in the modle


#######################################################
#### HEIGHT ####
#######################################################

hh <- dat4 %>% filter(var == "height") %>% 
  filter(!is.na(breed)) %>%  # remove NA's
  mutate(Year = factor(Year))

hist(log(hh$abs.slopes)) 

# By Hand
#res <- dredge(modAbsHeight, fixed = "r.temp", rank = "AICc", extra = "R^2")
#res <- dredge(modAbsHeight, fixed = "r.temp", rank = "AICc")


# RAW SLOPES
modHeight <- lmer(slopes ~ r.temp + mean.temp + mean.seasonality + dist.km + breed + growthform + intro + (1|species)+(1|StudyID)+(1|Year), data = hh, REML = FALSE, weights = sample.size, na.action = "na.fail")

# check model assumptions
fix.check(modHeight)

ResHeight <- ModelAverage(modHeight, printFullTable = TRUE, print95Table = TRUE, percent.thresh = 0.95)
ResHeight
#PlotEstimates(ResHeight)

model.set <- dredge(modHeight, fixed = "r.temp", rank = "AICc", extra = "R^2")
mm <- data.frame(model.set)
mm$cumsum <- cumsum(mm$weight)
ResH <- mm %>% 
  filter(cumsum < 0.95) %>% 
  mutate(breed = ifelse(!is.na(breed), "BS", "")) %>% 
  mutate(growthform = ifelse(!is.na(growthform), "GF", "")) %>% 
  mutate(intro = ifelse(!is.na(intro), "IS", "")) %>% 
  mutate(mean.seasonality = ifelse(!is.na(mean.seasonality), "MS", "")) %>% 
  mutate(mean.temp = ifelse(!is.na(mean.temp), "MT", "")) %>% 
  mutate(dist.km = ifelse(!is.na(dist.km), "DI", "")) %>% 
  mutate(Model = paste(breed, growthform, intro, dist.km, mean.seasonality, mean.temp, sep = "+")) %>% 
  mutate(delta = round(delta, 2), weight = round(weight, 3), R.square = round(R.2, 3)) %>% 
  select(Model, df, delta, weight, R.square) %>% 
  mutate(Model = gsub("^\\+*|\\+*$", "", Model)) %>% # removes + at start and end
  mutate(Model = gsub("\\+{2,}", "+", Model)) %>%  # removes 2-x + and replace by +; does not touch +
  mutate(Trait = "Height")

PlotHeight <- PlotEstimates3(modHeight, percent.thresh = 0.95)


# Raw slopes without sample size
modHeight.ww <- lmer(slopes ~ r.temp + mean.temp + mean.seasonality + breed + growthform + intro + (1|species)+(1|StudyID)+(1|Year), data = hh, REML = FALSE, na.action = "na.fail")
ResHeight.ww <- ModelAverage(modHeight.ww, printFullTable = FALSE, print95Table = TRUE, percent.thresh = 0.95)

model.set <- dredge(modHeight.ww, fixed = "r.temp", rank = "AICc")
mm <- data.frame(model.set)
mm$cumsum <- cumsum(mm$weight)
ResHww <- mm %>% 
  filter(cumsum < 0.95) %>% 
  mutate(breed = ifelse(!is.na(breed), "BS", "")) %>% 
  mutate(growthform = ifelse(!is.na(growthform), "GF", "")) %>% 
  mutate(intro = ifelse(!is.na(intro), "IS", "")) %>% 
  mutate(mean.seasonality = ifelse(!is.na(mean.seasonality), "MS", "")) %>% 
  mutate(mean.temp = ifelse(!is.na(mean.temp), "MT", "")) %>% 
  mutate(Model = paste(breed, growthform, intro, mean.seasonality, mean.temp, sep = "+")) %>% 
  mutate(delta = round(delta, 2), weight = round(weight, 3)) %>% 
  select(Model, df, delta, weight) %>% 
  mutate(Model = gsub("^\\+*|\\+*$", "", Model)) %>% # removes + at start and end
  mutate(Model = gsub("\\+{2,}", "+", Model)) %>%  # removes 2-x + and replace by +; does not touch +
  mutate(Trait = "Height withough weights")

# Sample size > 3
hh3 <- dat4 %>% filter(var == "height") %>% 
  filter(!is.na(breed)) %>%  # remove NA's
  mutate(Year = factor(Year)) %>% 
  filter(sample.size > 3)

modHeight3 <- lmer(slopes ~ r.temp + mean.temp + mean.seasonality + breed + growthform + intro + (1|species)+(1|StudyID)+(1|Year), data = hh3, REML = FALSE, weights = sample.size, na.action = "na.fail")

# check model assumptions
fix.check(modHeight3)

ResHeight3 <- ModelAverage(modHeight3, printFullTable = TRUE, print95Table = TRUE, percent.thresh = 0.95)



# Sample size > 5
hh5 <- dat4 %>% filter(var == "height") %>% 
  filter(!is.na(breed)) %>%  # remove NA's
  mutate(Year = factor(Year)) %>% 
  filter(sample.size > 5)

modHeight5 <- lmer(slopes ~ r.temp + mean.temp + mean.seasonality + breed + growthform + intro + (1|species)+(1|StudyID)+(1|Year), data = hh5, REML = FALSE, weights = sample.size, na.action = "na.fail")

# check model assumptions
fix.check(modHeight5)

ResHeight5 <- ModelAverage(modHeight5, printFullTable = TRUE, print95Table = TRUE, percent.thresh = 0.95)




#######################################################
#### BIOMASS ####
#######################################################

bb <- dat4 %>% filter(var == "biomass") %>% 
  filter(!is.na(breed)) %>%  # remove NA's
  mutate(Year = factor(Year))

hist(bb$slopes) 

# RAW SLOPES
modBiomass <- lmer(slopes ~ r.temp + mean.temp + mean.seasonality + dist.km + breed + growthform + intro + (1|species)+(1|StudyID)+(1|Year), data = bb, REML = FALSE, weights = sample.size, na.action = "na.fail")
fix.check(modBiomass)

ResBiomass <- ModelAverage(modBiomass, printFullTable = TRUE, print95Table = TRUE, percent.thresh = 0.95)

model.set <- dredge(modBiomass, fixed = "r.temp", rank = "AICc", extra = "R^2")
mm <- data.frame(model.set)
mm$cumsum <- cumsum(mm$weight)
ResB <- mm %>% 
  filter(cumsum < 0.95) %>% 
  mutate(breed = ifelse(!is.na(breed), "BS", "")) %>% 
  mutate(growthform = ifelse(!is.na(growthform), "GF", "")) %>% 
  mutate(intro = ifelse(!is.na(intro), "IS", "")) %>% 
  mutate(mean.seasonality = ifelse(!is.na(mean.seasonality), "MS", "")) %>% 
  mutate(mean.temp = ifelse(!is.na(mean.temp), "MT", "")) %>% 
  mutate(dist.km = ifelse(!is.na(dist.km), "DI", "")) %>% 
  mutate(Model = paste(breed, growthform, intro, dist.km, mean.seasonality, mean.temp, sep = "+")) %>% 
  mutate(delta = round(delta, 2), weight = round(weight, 3), R.square = round(R.2, 3)) %>% 
  select(Model, df, delta, weight, R.square) %>% 
  mutate(Model = gsub("^\\+*|\\+*$", "", Model)) %>% # removes + at start and end
  mutate(Model = gsub("\\+{2,}", "+", Model)) %>%  # removes 2-x + and replace by +; does not touch +
  mutate(Trait = "Biomass")


# Raw slope without sample size
modBiomass.ww <- lmer(slopes ~ r.temp + mean.temp + mean.seasonality + breed + growthform + intro + (1|species)+(1|StudyID)+(1|Year), data = bb, REML = FALSE, na.action = "na.fail")
fix.check(modBiomass.ww)

ResBiomass.ww <- ModelAverage(modBiomass.ww, printFullTable = TRUE, print95Table = TRUE, percent.thresh = 0.95)

model.set <- dredge(ResBiomass.ww, fixed = "r.temp", rank = "AICc")
mm <- data.frame(model.set)
mm$cumsum <- cumsum(mm$weight)
ResBww <- mm %>% 
  filter(cumsum < 0.95) %>% 
  mutate(breed = ifelse(!is.na(breed), "BS", "")) %>% 
  mutate(growthform = ifelse(!is.na(growthform), "GF", "")) %>% 
  mutate(intro = ifelse(!is.na(intro), "IS", "")) %>% 
  mutate(mean.seasonality = ifelse(!is.na(mean.seasonality), "MS", "")) %>% 
  mutate(mean.temp = ifelse(!is.na(mean.temp), "MT", "")) %>% 
  mutate(Model = paste(breed, growthform, intro, mean.seasonality, mean.temp, sep = "+")) %>% 
  mutate(delta = round(delta, 2), weight = round(weight, 3)) %>% 
  select(Model, df, delta, weight) %>% 
  mutate(Model = gsub("^\\+*|\\+*$", "", Model)) %>% # removes + at start and end
  mutate(Model = gsub("\\+{2,}", "+", Model)) %>%  # removes 2-x + and replace by +; does not touch +
  mutate(Trait = "Biomass withouth weights")


# Sample size > 3
bb3 <- dat4 %>% filter(var == "biomass") %>% 
  filter(!is.na(breed)) %>%  # remove NA's
  mutate(Year = factor(Year)) %>% 
  filter(sample.size > 3)

modBiomass3 <- lmer(slopes ~ r.temp + mean.temp + mean.seasonality + breed + growthform + intro + (1|species)+(1|StudyID)+(1|Year), data = bb3, REML = FALSE, weights = sample.size, na.action = "na.fail")

# check model assumptions
fix.check(modBiomass3)

ResBiomass3 <- ModelAverage(modBiomass3, printFullTable = TRUE, print95Table = TRUE, percent.thresh = 0.95)


# Sample size > 5
bb5 <- dat4 %>% filter(var == "biomass") %>% 
  filter(!is.na(breed)) %>%  # remove NA's
  mutate(Year = factor(Year)) %>% 
  filter(sample.size > 5)

modBiomass5 <- lmer(slopes ~ r.temp + mean.temp + mean.seasonality + breed + growthform + intro + (1|species)+(1|StudyID)+(1|Year), data = bb5, REML = FALSE, weights = sample.size, na.action = "na.fail")

# check model assumptions
fix.check(modBiomass5)

ResBiomass5 <- ModelAverage(modBiomass5, printFullTable = TRUE, print95Table = TRUE, percent.thresh = 0.95)



#######################################################
#### PHENOLOGY ####
#######################################################

pe <- dat4 %>% filter(var == "phenology") %>% 
  filter(!is.na(breed)) %>%  # remove NA's
  mutate(Year = factor(Year))

hist(pe$slopes)

# RAW SLOPES
modPhenologyEmergence <- lmer(slopes ~ r.temp + mean.temp + mean.seasonality + dist.km + breed + growthform + intro + (1|species)+(1|StudyID)+(1|Year), data = pe, REML = FALSE, weights = sample.size, na.action = "na.fail")
fix.check(modPhenologyEmergence)

ResPhenologyEmergence <- ModelAverage(modPhenologyEmergence, printFullTable = TRUE, print95Table = TRUE, percent.thresh = 0.95)

model.set <- dredge(modPhenologyEmergence, fixed = "r.temp", rank = "AICc", extra = "R^2")
mm <- data.frame(model.set)
mm$cumsum <- cumsum(mm$weight)
ResPE <- mm %>% 
  filter(cumsum < 0.95) %>% 
  mutate(breed = ifelse(!is.na(breed), "BS", "")) %>% 
  mutate(growthform = ifelse(!is.na(growthform), "GF", "")) %>% 
  mutate(intro = ifelse(!is.na(intro), "IS", "")) %>% 
  mutate(mean.seasonality = ifelse(!is.na(mean.seasonality), "MS", "")) %>% 
  mutate(mean.temp = ifelse(!is.na(mean.temp), "MT", "")) %>% 
  mutate(dist.km = ifelse(!is.na(dist.km), "DI", "")) %>% 
  mutate(Model = paste(breed, growthform, intro, dist.km, mean.seasonality, mean.temp, sep = "+")) %>% 
  mutate(delta = round(delta, 2), weight = round(weight, 3), R.square = round(R.2, 3)) %>% 
  select(Model, df, delta, weight, R.square) %>% 
  mutate(Model = gsub("^\\+*|\\+*$", "", Model)) %>% # removes + at start and end
  mutate(Model = gsub("\\+{2,}", "+", Model)) %>%  # removes 2-x + and replace by +; does not touch +
  mutate(Trait = "Phenology")


# RAW slopes withouth weight
modPhenologyEmergence.ww <- lmer(slopes ~ r.temp + mean.temp + mean.seasonality + breed + growthform + intro + (1|species)+(1|StudyID)+(1|Year), data = pe, REML = FALSE, na.action = "na.fail")
fix.check(modPhenologyEmergence.ww)

ResPhenologyEmergence.ww <- ModelAverage(modPhenologyEmergence.ww, printFullTable = TRUE, print95Table = TRUE, percent.thresh = 0.95)

model.set <- dredge(modPhenologyEmergence.ww, fixed = "r.temp", rank = "AICc")
mm <- data.frame(model.set)
mm$cumsum <- cumsum(mm$weight)
ResPEww <- mm %>% 
  filter(cumsum < 0.95) %>% 
  mutate(breed = ifelse(!is.na(breed), "BS", "")) %>% 
  mutate(growthform = ifelse(!is.na(growthform), "GF", "")) %>% 
  mutate(intro = ifelse(!is.na(intro), "IS", "")) %>% 
  mutate(mean.seasonality = ifelse(!is.na(mean.seasonality), "MS", "")) %>% 
  mutate(mean.temp = ifelse(!is.na(mean.temp), "MT", "")) %>% 
  mutate(Model = paste(breed, growthform, intro, mean.seasonality, mean.temp, sep = "+")) %>% 
  mutate(delta = round(delta, 2), weight = round(weight, 3)) %>% 
  select(Model, df, delta, weight) %>% 
  mutate(Model = gsub("^\\+*|\\+*$", "", Model)) %>% # removes + at start and end
  mutate(Model = gsub("\\+{2,}", "+", Model)) %>%  # removes 2-x + and replace by +; does not touch +
  mutate(Trait = "Phenology without weights")


# Sample size > 3
pe3 <- dat4 %>% filter(var == "phenology") %>% 
  filter(!is.na(breed)) %>%  # remove NA's
  mutate(Year = factor(Year)) %>% 
  filter(sample.size > 3)

modPhenology3 <- lmer(slopes ~ r.temp + mean.temp + mean.seasonality + breed + growthform + intro + (1|species)+(1|StudyID)+(1|Year), data = pe3, REML = FALSE, weights = sample.size, na.action = "na.fail")

# check model assumptions
fix.check(modPhenology3)

ResPhenology3 <- ModelAverage(modPhenology3, printFullTable = TRUE, print95Table = TRUE, percent.thresh = 0.95)



# Sample size > 5
pe5 <- dat4 %>% filter(var == "phenology") %>% 
  filter(!is.na(breed)) %>%  # remove NA's
  mutate(Year = factor(Year)) %>% 
  filter(sample.size > 5)

modPhenology5 <- lmer(slopes ~ r.temp + mean.temp + mean.seasonality + breed + growthform + intro + (1|species)+(1|StudyID)+(1|Year), data = pe5, REML = FALSE, weights = sample.size, na.action = "na.fail")

# check model assumptions
fix.check(modPhenology5)

ResPhenology5 <- ModelAverage(modPhenology5, printFullTable = TRUE, print95Table = TRUE, percent.thresh = 0.95)



#######################################################
#### RESULT TABLE ####
#######################################################

# Main
ResHeight %>% 
  bind_rows(ResBiomass, ResPhenologyEmergence) %>% 
  mutate(Trait = c(rep("Height", 10), rep("Biomass", 10), rep("Phenology", 10))) %>% 
  select(Trait, Variable, Estimate, CI, Zvalue, Pvalue, Importance) %>% 
  mutate(Variable = plyr:: mapvalues(Variable, c("(Intercept)", "T Range", "Mean T", "dist.km", "Native", "T Seasonality","Grass", "Tree", "Mixed Mating", "Outcrossing"), c("Intercept", "T Range", "Mean T", "Distance", "Introduction status:Native", "T Seasonality","Growthform:Grass", "Growthform:Tree", "Breeding System:Mixed Mating", "Breeding Sytem:Outcrossing"))) %>% 
  mutate(Variable = factor(Variable, levels = c("Intercept", "T Range", "Mean T", "T Seasonality", "Distance", "Introduction status:Native", "Growthform:Grass", "Growthform:Tree", "Breeding System:Mixed Mating", "Breeding Sytem:Outcrossing"))) %>% 
  mutate(Trait = factor(Trait, levels = c("Height", "Biomass", "Phenology"))) %>% 
  arrange(Trait, Variable) %>% 
  mutate(CI = paste("(", CI, ")", sep = ""))

# Additional tests
ResHeight.ww %>% 
  bind_rows(ResBiomass.ww, ResPhenologyEmergence.ww, ResHeight3, ResHeight5, ResBiomass3, ResBiomass5, ResPhenology3, ResPhenology5) %>% 
  mutate(Trait = c(rep("Height", 9), rep("Biomass", 9), rep("Phenology", 9), rep("Height", 9), rep("Height", 9), rep("Biomass", 9), rep("Biomass", 9), rep("Phenology", 9), rep("Phenology", 7))) %>% 
  mutate(Test = c(rep("withouth weight", 9), rep("withouth weight", 9), rep("withouth weight", 9), rep("sample size > 3", 9), rep("sample size > 5", 9), rep("sample size > 3", 9), rep("sample size > 5", 9), rep("sample size > 3", 9), rep("sample size > 5", 7))) %>% 
  select(Test, Trait, Variable, Estimate, CI, Zvalue, Pvalue, Importance) %>% 
  mutate(Variable = plyr:: mapvalues(Variable, c("(Intercept)", "T Range", "Mean T", "Native", "T Seasonality","Grass", "Tree", "Mixed Mating", "Outcrossing"), c("Intercept", "T Range", "Mean T", "Introduction status:Native", "T Seasonality","Growthform:Grass", "Growthform:Tree", "Breeding System:Mixed Mating", "Breeding Sytem:Outcrossing"))) %>% 
  mutate(Variable = factor(Variable, levels = c("Intercept", "T Range", "Mean T", "T Seasonality", "Introduction status:Native", "Growthform:Grass", "Growthform:Tree", "Breeding System:Mixed Mating", "Breeding Sytem:Outcrossing"))) %>% 
  mutate(Trait = factor(Trait, levels = c("Height", "Biomass", "Phenology"))) %>% 
  arrange(Test, Trait, Variable)%>% 
  mutate(CI = paste("(", CI, ")", sep = ""))



#######################################################
#### MODEL TABLE 2 ####
#######################################################

ResH %>% 
  bind_rows(ResB, ResPE) %>% 
  select(Trait, Model, df, delta, R.square, weight) %>% 
  group_by(Trait) %>% 
  summarise(n())




#######################################################
#### PLOT IMPORTANCE ####
#######################################################

# ROW_BIND PARAMETER TABLE
ImportancePlot <- ResHeight %>% 
  bind_rows(ResBiomass, ResPhenologyEmergence) %>% 
  filter(Variable != "(Intercept)") %>% 
  mutate(Trait = c(rep("Height", 9), rep("Biomass", 9), rep("Phenology", 9))) %>% 
  mutate(Trait = factor(Trait, level = c("Height", "Biomass", "Phenology"))) %>% 
  group_by(Trait, Category) %>% 
  summarise(Importance = mean(Importance)) %>% 
  mutate(Category = factor(Category, level = c("T Range", "Mean T", "T Seasonality", "dist.km", "Breeding System", "Growthform", "Introduction Status"))) %>%
  arrange(Trait, -Importance) %>% 
  ggplot(aes(x = Category, y = Importance, fill = Trait)) +
  geom_bar(stat = "identity", position = "dodge")
ImportancePlot
ggsave("Fig3ImportancePlot.pdf", ImportancePlot)


ResHeight.ww %>% 
  bind_rows(ResBiomass.ww, ResPhenologyEmergence.ww, ResHeight3, ResHeight5, ResBiomass3, ResBiomass5, ResPhenology3, ResPhenology5) %>% 
  mutate(Trait = c(rep("Height", 9), rep("Biomass", 9), rep("Phenology", 9), rep("Height", 9), rep("Height", 9), rep("Biomass", 9), rep("Biomass", 9), rep("Phenology", 9), rep("Phenology", 7))) %>% 
  mutate(Test = c(rep("withouth weight", 9), rep("withouth weight", 9), rep("withouth weight", 9), rep("sample size > 3", 9), rep("sample size > 5", 9), rep("sample size > 3", 9), rep("sample size > 5", 9), rep("sample size > 3", 9), rep("sample size > 5", 7))) %>% 
  filter(Variable != "(Intercept)") %>% 
  mutate(Trait = factor(Trait, level = c("Height", "Biomass", "Phenology"))) %>% 
  group_by(Test, Trait, Category) %>% 
  summarise(Importance = mean(Importance)) %>% 
  mutate(Category = factor(Category, level = c("T Range", "Mean T", "T Seasonality", "Breeding System", "Growthform", "Introduction Status"))) %>%
  arrange(Test, Trait, -Importance) %>% 
  ggplot(aes(x = Category, y = Importance, fill = Trait)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  facet_wrap(~ Test)



#######################################################
#### PLOT ESTIMATE ####
#######################################################
# Plot Estimates
ParameterEstimatePlot <- ResHeight %>% 
  bind_rows(ResBiomass, ResPhenologyEmergence) %>% 
  mutate(Trait = c(rep("Height", 10), rep("Biomass", 10), rep("Phenology", 10))) %>% 
  filter(Variable != c("(Intercept)")) %>% 
  filter(Variable != c("T Range")) %>% 
  separate(col = CI, into = c("CI.low", "CI.high"), sep = " - ", convert = TRUE) %>% 
  mutate(Variable = plyr::mapvalues(Variable, c("Outcrossing", "Mixed Mating", "Tree", "Grass", "Native", "T Seasonality", "Mean T", "Distance"), c("Outcrossing", "Mixed Mating", "Tree", "Grass", "Native", "T Seasonality", "Mean T", "Max Distance"))) %>% 
  mutate(Variable = factor(Variable, levels = c("Outcrossing", "Mixed Mating", "Tree", "Grass", "Native", "Max Distance", "T Seasonality", "Mean T"))) %>% 
  mutate(Trait = factor(Trait, levels = c("Height", "Biomass", "Phenology"))) %>% 
  ggplot(aes(y = Variable, x = Estimate, xmin = CI.low, xmax = CI.high)) + 
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
  geom_point() +
  labs(x = "Parameter estimate", y = "") +
  geom_errorbarh(height = 0) +
  panel_border(colour = "black", remove = FALSE) +
  facet_wrap(~ Trait) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black"))
  
ggsave("Fig3BParameterEstimatePlot.pdf", ParameterEstimatePlot, height = 5)


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
  #gather(key = Slope, value = Value, -Trait) %>% 
  ggplot(aes(x = Trait, y = Slope)) + 
  geom_violin(draw_quantiles = c(0.5), position = position_dodge(width = 0.60), fill = "grey80") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_jitter(width = 0.2, aes(shape = Pvalue)) +
  labs(x = "", y = "Temperature effect (regression slope)") +
  scale_shape_manual(values = c(1,17)) +
  stat_summary(fun.data = "mean_cl_normal", colour = "black", alpha = 0.8, size = 1.2, shape = 16, position = position_dodge(width = 0.60)) +
  theme(legend.position = "none")

ggsave("Fig2EffectSizePlot.pdf", Fig2EffectSizePlot)


# Check different sample sizes
smaller5 <- dat4 %>% filter(r.temp < 5)
s5 <- smaller5 %>% 
  select(var, slopes, p.value) %>% 
  setNames(., c("Trait", "Slope", "Pvalue")) %>% 
  mutate(Trait = plyr::mapvalues(Trait, c("height", "biomass", "phenology"), c("Height", "Biomass", "Phenology"))) %>%
  mutate(Trait = factor(Trait, level = c("Height", "Biomass", "Phenology"))) %>%
  mutate(Pvalue = ifelse(Pvalue < 0.05, "sign.", "non sign.")) %>% 
  ggplot(aes(x = Trait, y = Slope)) + 
  geom_violin(draw_quantiles = c(0.5), position = position_dodge(width = 0.60), fill = "grey80") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_jitter(width = 0.2, aes(shape = Pvalue)) +
  labs(x = "", y = "Temperature effect (regression slope)") +
  scale_shape_manual(values = c(1,17)) +
  stat_summary(fun.data = "mean_cl_normal", colour = "black", alpha = 0.8, size = 1.2, shape = 16, position = position_dodge(width = 0.60)) +
  annotate("text", x = 2, y = 1.4, label = "Temperature range < 5°C", size= 6) +
  theme(legend.position = "none")


larger5 <- dat4 %>% filter(r.temp > 5)
l5 <- larger5 %>% 
  select(var, slopes, p.value) %>% 
  setNames(., c("Trait", "Slope", "Pvalue")) %>% 
  mutate(Trait = plyr::mapvalues(Trait, c("height", "biomass", "phenology"), c("Height", "Biomass", "Phenology"))) %>%
  mutate(Trait = factor(Trait, level = c("Height", "Biomass", "Phenology"))) %>%
  mutate(Pvalue = ifelse(Pvalue < 0.05, "sign.", "non sign.")) %>% 
  ggplot(aes(x = Trait, y = Slope)) + 
  geom_violin(draw_quantiles = c(0.5), position = position_dodge(width = 0.60), fill = "grey80") +
  geom_hline(yintercept = 0, color = "grey") +
  geom_jitter(width = 0.2, aes(shape = Pvalue)) +
  labs(x = "", y = "Temperature effect (regression slope)") +
  scale_shape_manual(values = c(1,17)) +
  stat_summary(fun.data = "mean_cl_normal", colour = "black", alpha = 0.8, size = 1.2, shape = 16, position = position_dodge(width = 0.60)) +
  annotate("text", x = 2, y = 0.5, label = "Temperature range > 5°C", size= 6) +
  theme(legend.position = "none")

EffectSizePlot5 <- plot_grid(s5, l5, label_size = 12, nrow = 2, ncol = 1)



dd <- dat4 %>% 
  mutate(Five = ifelse(r.temp < 5, "high", "low")) %>% 
  filter(var == "phenology", Five == "low") %>% distinct(StudyID)

summary(lm(slopes ~ r.temp, dd))

ggplot(dat4, aes(x = r.temp, y = slopes)) + geom_point() + facet_wrap(~ var) + geom_vline(xintercept = 5)



modAbsNullHeight <- lmer(log(abs.slopes) ~ 1 + (1|species)+(1|StudyID)+(1|Year), data = hh, REML = FALSE, weights = sample.size, na.action = "na.fail")
summary(modAbsNullHeight)
exp(-1.6600) -1.96*exp(0.2428)
exp(-1.6600) +1.96*exp(0.2428)

