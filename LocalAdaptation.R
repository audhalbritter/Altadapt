##############################################################################################
###### LOCAL ADAPTATION ANALYSIS ######
##############################################################################################

# Load libraries
library("cowplot")
library("gridExtra")
library("grid")

# Read in data
LA <- read.csv(file = "FinalSharedData/LocalAdaptation.csv")

# BIOMASS
biomass <- LA %>% filter(var == "biomass") %>% filter(!is.na(st.traitmean))
fit1 <- lmer((st.traitmean) ~ symp * r.temp + (1|PopID)+(1|species)+(1|StudyID), biomass, na.action = "na.fail", REML = FALSE) 
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
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "symphome site", "r.temp", "r.temp:symphome site"), c("Intercept", "Own-Foreign", "T Range", "Own-Foreign * T Range")))


# FITNESS
fitness <- LA %>% filter(var == "fitness") %>% filter(!is.na(st.traitmean))
hist((fitness$st.traitmean))
fit1 <- lmer((st.traitmean) ~ symp * r.temp + (1|PopID)+(1|species)+(1|StudyID), fitness, na.action = "na.fail", REML = FALSE) # better model
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
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "symphome site", "r.temp", "r.temp:symphome site"), c("Intercept", "Own-Foreign", "T Range", "Own-Foreign * T Range")))


# SURVIVAL
survival <- LA %>% filter(var == "survival") %>% filter(!is.na(st.traitmean))
hist((survival$st.traitmean))
fit1 <- lmer((st.traitmean) ~ symp * r.temp + (1|PopID)+(1|species)+(1|StudyID), survival, na.action = "na.fail", REML = FALSE) # better model
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
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "symphome site", "r.temp", "r.temp:symphome site"), c("Intercept", "Own-Foreign", "T Range", "Own-Foreign * T Range")))

ResultsLA <- BiomassLA %>% 
  bind_rows(FitnessLA, SurvivalLA) %>% 
  mutate(Trait = c(rep("Biomass", 4), rep("Rep. output", 4), rep("Survival", 4))) %>% 
  select(Trait, Variable, Estimate, CI, Zvalue, Pvalue)



################################################
#### DIFFERENCE BETWEEN HIGH AND LOW POPS ####
################################################

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
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "factor.site2low", "symphome site", "factor.site2low:symphome site"), c("Intercept", "Transplant site", "Own-foreign", "Transplant site * Own-foreign")))



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
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "factor.site2low", "symphome site", "factor.site2low:symphome site"), c("Intercept", "Transplant site", "Own-foreign", "Transplant site * Own-foreign")))



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
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "factor.site2low", "symphome site", "factor.site2low:symphome site"), c("Intercept", "Transplant site", "Own-foreign", "Transplant site * Own-foreign")))



ResultsHL <- BiomassHL %>% 
  bind_rows(FitnessHL, SurvivalHL) %>% 
  mutate(Trait = c(rep("Biomass", 4), rep("Rep. output", 4), rep("Survival", 4))) %>% 
  select(Trait, Variable, Estimate, Pvalue)





#### TEST THE EFFCET OF OTHER VARIABLES ####

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
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "symphome site", "breedmixed_mating", "breedoutcrossing", "breedmixed_mating:symphome site", "breedoutcrossing:symphome site"), c("Intercept", "OF", "BS - mixed mating", "BS - outcrossing", "OF*BS - mixed mating", "OF*BS - ourcrossing")))


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
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "symphome site", "breedmixed_mating", "breedoutcrossing", "breedmixed_mating:symphome site", "breedoutcrossing:symphome site"), c("Intercept", "OF", "BS - mixed mating", "BS - outcrossing", "OF*BS - mixed mating", "OF*BS - ourcrossing")))


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
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "symphome site", "breedmixed_mating", "breedmixed_mating:symphome site"), c("Intercept", "OF", "BS - mixed mating", "OF*BS - mixed mating")))


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
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "symphome site", "growthformgrass", "growthformmoss", "growthformtree", "growthformgrass:symphome site", "growthformmoss:symphome site", "growthformtree:symphome site"), c("Intercept", "OF", "GF - grass", "GF - moss", "GF - tree", "OF*GF - grass", "OF*GF - moss", "OF*GF - tree")))


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
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "symphome site", "growthformgrass", "growthformgrass:symphome site"), c("Intercept", "OF", "GF - grass", "OF*GF - grass")))


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
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "symphome site", "growthformforb", "growthformgrass", "growthformforb:symphome site", "growthformgrass:symphome site"), c("Intercept", "OF", "GF - forb", "GF - grass", "OF*GF - forb", "OF*GF - grass")))


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
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "symphome site", "longevityperennial", "longevityperennial:symphome site"), c("Intercept", "OF", "LO - perennial", "OF*LO - perennial")))


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
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "symphome site", "longevityperennial", "longevityperennial:symphome site"), c("Intercept", "OF", "LO - perennial", "OF*LO - perennial")))


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
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "symphome site", "longevityperennial", "longevityperennial:symphome site"), c("Intercept", "OF", "LO - perennial", "OF*LO - perennial")))



resBioBreed %>% 
  bind_rows(resFitBreed, resSurBreed, resBioGrowth, resFitGrowth, resSurGrowth, resBioLong, resFitLong, resSurLong) %>% 
  mutate(Trait = c(rep("Biomass", 6), rep("Fitness", 6), rep("Survival", 4), rep("Biomass", 8), rep("Fitness", 4), rep("Survival", 6), rep("Biomass", 4), rep("Fitness", 4), rep("Survival", 4))) %>% 
  mutate(LifeHistory = c(rep("Breeding system", 16), rep("Growthform", 18), rep("Longevity", 12))) %>% 
  select(Trait, LifeHistory, Variable, Estimate, CI, Zvalue, Pvalue) %>% 
  mutate(Star = ifelse(Pvalue <= 0.05, "***", Pvalue))



####################################################
#### Fig. 4 ADAPTATION ####
####################################################

# define range
yylim <- c(min(biomass$st.traitmean, na.rm = TRUE), max(biomass$st.traitmean, na.rm = TRUE))

p4 <- LA %>%
  filter(var == "biomass") %>% 
  group_by(var, symp) %>% 
  summarise(n = n(), mean = mean(st.traitmean), se = sd(st.traitmean)/sqrt(n), CI.low = mean - 1.96*se, CI.high = mean + 1.96*se) %>% 
  mutate(symp = plyr::mapvalues(symp, c("home site", "foreign site"), c("own", "foreign"))) %>%
  mutate(symp = factor(symp, level = c("own", "foreign"))) %>% 
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
  mutate(symp = factor(symp, level = c("own", "foreign"))) %>% 
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
  mutate(symp = factor(symp, level = c("own", "foreign"))) %>% 
  ggplot(aes(x = symp, y = mean, ymin = CI.low, ymax = CI.high)) +
  geom_hline(yintercept = 0, colour = "grey", linetype = "longdash") +
  geom_point(size = 3) +
  geom_errorbar(width = 0) +
  labs(x = "", y = "", title = "Survival") +
  annotate("text", x = 0.6, y = 0.7, label = "A", size = 5) +
  ylim(-0.7, 0.8) +
  panel_border(colour = "black", remove = FALSE) +
  annotate("text", x = 2.1, y = 0.7, label = "OF*", size= 5)



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




AdaptationPlot <- grid.arrange(p6, p5, p4, p9, p8, p7, layout_matrix = rbind(c(1,2,3),c(4,5,6)), bottom = textGrob("Elevation of transplant site", vjust = 0.1, gp = gpar(fontsize = 15, font = 8)), left = textGrob("Standardized trait mean", rot = 90, vjust = 1, gp = gpar(fontsize = 15, font = 8)))
ggsave("FinalFiugres/Fig4_Adaptation4.pdf", AdaptationPlot, height = 6, width = 10)
