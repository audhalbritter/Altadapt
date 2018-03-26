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


#######################################################
#### HEIGHT ####
#######################################################

hh <- dat4 %>% filter(var == "height") %>% 
  filter(!is.na(breed)) %>%  # remove NA's
  mutate(Year = factor(Year))

hist(log(hh$slopes)) 

# By Hand
#res <- dredge(modAbsHeight, fixed = "r.temp", rank = "AICc", extra = "R^2")
#res <- dredge(modAbsHeight, fixed = "r.temp", rank = "AICc")

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

hist(bb$slopes) 

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


model.set <- dredge(modBiomass_ABS, fixed = "r.temp", rank = "AICc", extra = "R^2")
mm <- data.frame(model.set)
mm$cumsum <- cumsum(mm$weight)
mm95 <- mm %>% filter(cumsum < percent.thresh)
averaged.model <- model.avg(model.set, cumsum(weight) <= 0.95)
res <- data.frame(summary(averaged.model)$coefmat.full)







#######################################################
#### PHENOLOGY ####
#######################################################

pe <- dat4 %>% filter(var == "phenology") %>% 
  filter(!is.na(breed)) %>%  # remove NA's
  mutate(Year = factor(Year))

hist(pe$slopes)

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

write_csv(ResultEstimate, path = "FinalResults/ResultEstimate.csv", col_names = TRUE)

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
#### MODEL TABLE S4 ####
#######################################################

ResModelSelection <- ResH_UT %>% 
  bind_rows(ResH_ABS, ResB_UT, ResB_ABS, ResP_UT, ResP_ABS) %>% 
  select(Trait, Model, df, delta, R.square, weight)
  
write_csv(ResModelSelection, path = "FinalResults/ResModelSelection.csv", col_names = TRUE)

# count models
ResModelSelection %>%   
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





#######################################################
#### PLOT ESTIMATE ####
#######################################################
# Plot Estimates

ParameterEstimate <- ResHeight_UT %>% 
  bind_rows(ResHeight_ABS, ResBiomass_UT, ResBiomass_ABS, ResPhenologyEmergence_UT, ResPhenologyEmergence_ABS)
save(ParameterEstimate, file = "ParameterEstimate.Rdata")

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

