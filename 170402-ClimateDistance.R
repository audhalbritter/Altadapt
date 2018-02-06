#####################
# CLIMATE DISTANCE #
#####################

# select important variables from site
site.CD <- DAT.S %>%
  select(ID, studysite, bio_1_site, elevation_site)

# use only studies from dat2
CD <- DAT.T %>%
  filter(StudyID %in% dat2$StudyID) %>% 
  # add studysite, elev and bio_1 from the sites
  left_join(site.CD, by = "ID") %>% 
  # calculate elev and clim.diff between site and pop (only useful for Gardens) and make it a factor (High/Low)
  mutate(elev.diff = elevation_site - elev) %>% 
  mutate(clim.diff = bio_1_site - bio_1) %>% 
  mutate(clim.diff2 = clim.diff^2) %>% 
  mutate(rel.elev = factor(ifelse(elev.diff > 0, -1, 1))) %>% 
  filter(studysite == "Garden") %>% 
  filter(!is.na(elev.diff)) %>% 
  filter(!is.na(clim.diff)) %>% 
  # standardize traitmean by studyunit
  group_by(studyunit) %>% 
  mutate(std.traitmean = scale(traitmean)[,1]) %>% 
  mutate(abs.std.traitmean = abs(std.traitmean))



clim <- CD %>% 
  filter(TRAIT_CAT1 == "SIZE" & TRAIT_CAT2 == "BIOMASS" & TRAIT_CAT3 == "ABOVEGR" | TRAIT_CAT1 == "FITNESS" | TRAIT_CAT1 == "VITAL_RATES" & TRAIT_CAT2 == "SURVIVAL") %>%
  mutate(var = ifelse(TRAIT_CAT1 == "SIZE" & TRAIT_CAT2 == "BIOMASS" & TRAIT_CAT3 == "ABOVEGR", "Biomass",
                      ifelse(TRAIT_CAT1 == "FITNESS", "Fitness", "Survival")))

clim %>% 
  group_by(var) %>%
  distinct(SiteID) %>% 
  summarise(n())

clim %>% 
  group_by(var) %>%
  filter(StudyID == "Alex2010JBio") %>% 
  select(StudyID, studyunit)
  summarise(n())


# ABOVEGROUND BIOMASS
d4 <- clim %>% filter(var == "Biomass", !is.na(std.traitmean))
BiomassFull <- lmer(std.traitmean ~ clim.diff + clim.diff2 + (1|species), REML = FALSE, d4, na.action = "na.fail")

model.set <- dredge(BiomassFull, rank = "AICc")
model.set
averaged.model <- model.avg(model.set)
res <- data.frame(summary(averaged.model)$coefmat.full)
BiomassCD <- res %>% 
  rownames_to_column(var = "Variable") %>% 
  setNames(., c("Variable", "Estimate", "StError", "AdjSE", "Zvalue", "Pvalue")) %>% 
  mutate(CI.low = Estimate - 1.96 * StError) %>% 
  mutate(CI.high = Estimate + 1.96 * StError) %>% 
  mutate(Estimate = round(Estimate, 2), CI = paste(round(CI.low, 2), round(CI.high, 2), sep = " - "), Zvalue = round(Zvalue, 2), Pvalue = round(Pvalue, 3)) %>% 
  select(Variable, Estimate, CI, Zvalue, Pvalue) %>%
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "clim.diff", "clim.diff2"), c("Intercept", "Climate difference", "Climate difference^2")))



# Count the number of species per studies!
library("tidyr")
library("dplyr")
d4 %>%
  select(StudyID, species) %>%
  group_by(StudyID) %>%
  summarise(n = n_distinct(species))
# mainly one study, only a few have more


# REPRODUCTIVE OUTPUT
d5 <- clim %>% filter(var == "Fitness", !is.na(std.traitmean))
FintessFull <- lmer(std.traitmean ~ clim.diff + clim.diff2 + (1|species), REML = FALSE, d5, na.action = "na.fail")

model.set <- dredge(FintessFull, rank = "AICc")
model.set
averaged.model <- model.avg(model.set)
res <- data.frame(summary(averaged.model)$coefmat.full)
FitnessCD <- res %>% 
  rownames_to_column(var = "Variable") %>% 
  setNames(., c("Variable", "Estimate", "StError", "AdjSE", "Zvalue", "Pvalue")) %>% 
  mutate(CI.low = Estimate - 1.96 * StError) %>% 
  mutate(CI.high = Estimate + 1.96 * StError) %>% 
  mutate(Estimate = round(Estimate, 2), CI = paste(round(CI.low, 2), round(CI.high, 2), sep = " - "), Zvalue = round(Zvalue, 2), Pvalue = round(Pvalue, 3)) %>% 
  select(Variable, Estimate, CI, Zvalue, Pvalue) %>%
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "clim.diff", "clim.diff2"), c("Intercept", "Climate difference", "Climate difference^2")))



### SURVIVAL
d6 <- clim %>% filter(var == "Survival", !is.na(std.traitmean))
SurvivalFull <- lmer(std.traitmean ~ clim.diff + clim.diff2 + (1|species), REML = FALSE, d6, na.action = "na.fail")

model.set <- dredge(SurvivalFull, rank = "AICc")
model.set
averaged.model <- model.avg(model.set)
res <- data.frame(summary(averaged.model)$coefmat.full)
SurvivalCD <- res %>% 
  rownames_to_column(var = "Variable") %>% 
  setNames(., c("Variable", "Estimate", "StError", "AdjSE", "Zvalue", "Pvalue")) %>% 
  mutate(CI.low = Estimate - 1.96 * StError) %>% 
  mutate(CI.high = Estimate + 1.96 * StError) %>% 
  mutate(Estimate = round(Estimate, 2), CI = paste(round(CI.low, 2), round(CI.high, 2), sep = " - "), Zvalue = round(Zvalue, 2), Pvalue = round(Pvalue, 3)) %>% 
  select(Variable, Estimate, CI, Zvalue, Pvalue) %>%
  mutate(Variable = plyr::mapvalues(Variable, c("(Intercept)", "clim.diff", "clim.diff2"), c("Intercept", "Climate difference", "Climate difference^2")))



hist(d6$std.traitmean)
fit0 <- lmer(std.traitmean ~ 1 + (1|species), d6)
fit1 <- lmer(std.traitmean ~ clim.diff + (1|species), d6)
fit2 <- lmer(std.traitmean ~ clim.diff + I(clim.diff^2) + (1|species), d6)
modsel(list(fit0, fit1, fit2),1000)
summary(fit2)
fix.check(fit0)
unique(d6$Trait)
hist(sqrt(d6$st.traitmean))
plot(resid(mod01))


ResultsCD <- BiomassCD %>% 
  bind_rows(FitnessCD, SurvivalCD) %>% 
  mutate(Trait = c(rep("Biomass", 3), rep("Rep. output", 3), rep("Survival", 3))) %>% 
  select(Trait, Variable, Estimate, StError, Pvalue)



### FIG.4 CLIMATE DISTANCE
pdf ("Fig4_Climate Distance.pdf", width=10, height=4, pointsize=5, onefile=TRUE, paper="special")
par(mfrow=c(1,3), mar=c(1,2,2,0.5), mgp=c(3,1,0), oma=c(3,3,1,1), cex=1.7)

fit1 <- lm(std.traitmean ~ clim.diff, d4)
fit2 <- lm(std.traitmean ~ clim.diff + I(clim.diff^2), d4)
plot(std.traitmean ~ clim.diff, d4, ylab="", xlab="", axes=FALSE, main = "a) Aboveground biomass")
box(); axis(1, labels=TRUE); axis(2,labels=TRUE)
abline(v=0, col="darkgrey", lty=2)
abline(fit1, col="black")
lines(sort(d4$clim.diff), fitted(fit2)[order(d4$clim.diff)], col='black', lty=2, type='l')
#f <- summary(fit1)$fstatistic
legend("topright", c("linear", "quadratic"), col=c("black"), lty=c(1,2), bty="n")

fit1 <- lm(std.traitmean ~ clim.diff, d5)
fit2 <- lm(std.traitmean ~ clim.diff + I(clim.diff^2), d5)
f <- summary(fit1)$fstatistic
plot(std.traitmean ~ clim.diff, d5, ylab="", xlab="", axes=FALSE, main = "b) Reproductive output")
box(); axis(1, labels=TRUE); axis(2,labels=TRUE)
abline(v=0, col="darkgrey", lty=2)
abline(fit1, col="black")
lines(sort(d5$clim.diff), fitted(fit2)[order(d5$clim.diff)], col='black', lty=2, type='l')



fit1 <- lm(std.traitmean ~ clim.diff, d6)
fit2 <- lm(std.traitmean ~ clim.diff + I(clim.diff^2), d6)
anova(fit1, fit2)
AIC(fit1, fit2)
summary(fit1)
f <- summary(fit1)$fstatistic
plot(std.traitmean ~ clim.diff, d6, ylab="", xlab="", axes=FALSE, main = "c) Survival")
box(); axis(1, labels=TRUE); axis(2,labels=TRUE)
abline(v=0, col="darkgrey", lty=2)
abline(fit1, col="black")
lines(sort(d6$clim.diff), fitted(fit2)[order(d6$clim.diff)], col='black', lty=2, type='l')
text(-4,2.7, "upward")
text(10,2.7, "downward Transplant")

mtext("Standardized mean", 2, outer=TRUE, line=1, cex=2)
mtext("Temperature difference in °C", 1, outer=TRUE, line=1.5, cex=2)
dev.off()




pdf ("Fig4_ABS_Climate Distance.pdf", width=10, height=4, pointsize=5, onefile=TRUE, paper="special")
par(mfrow=c(1,2), mar=c(1,4,2,0.5), mgp=c(3,1,0), oma=c(3,1,1,1), cex=1.7)

fit1 <- lm(abs.st.traitmean ~ clim.diff, d5)
fit2 <- lm(abs.st.traitmean ~ clim.diff + I(clim.diff^2), d5)
anova(fit1, fit2)
AIC(fit1, fit2)
summary(fit1)
f <- summary(fit1)$fstatistic
plot(abs.st.traitmean ~ clim.diff, d5, ylab="Abs. stand. reproductive output", xlab="", axes=FALSE)
box(); axis(1, labels=TRUE); axis(2,labels=TRUE)
abline(v=0, col="darkgrey", lty=2)
abline(fit1, col="black")
lines(sort(d5$clim.diff), fitted(fit2)[order(d5$clim.diff)], col='black', lty=2, type='l')
legend("topleft", c("linear", "quadratic"), col=c("black"), lty=c(1,2), bty="n")

# Null model is the best, and species has almost no variance, so lm is probably ok
fit0 <- lmer(abs.st.traitmean ~ 1 + (1|species), d5)
fit1 <- lmer(abs.st.traitmean ~ clim.diff + (1|species), d5)
fit2 <- lmer(abs.st.traitmean ~ clim.diff + I(clim.diff^2) + (1|species), d5)
modsel(list(fit0, fit1, fit2),1000)
summary(fit1)


fit1 <- lm(abs.st.traitmean ~ clim.diff, d6)
fit2 <- lm(abs.st.traitmean ~ clim.diff + I(clim.diff^2), d6)
anova(fit1, fit2)
AIC(fit1, fit2)
summary(fit1)
f <- summary(fit1)$fstatistic
plot(abs.st.traitmean ~ clim.diff, d6, ylab="Abs. stand. survival", xlab="", axes=FALSE)
box(); axis(1, labels=TRUE); axis(2,labels=TRUE)
abline(v=0, col="darkgrey", lty=2)
abline(fit1, col="black")
lines(sort(d6$clim.diff), fitted(fit2)[order(d6$clim.diff)], col='black', lty=2, type='l')
text(-4,2.7, "upwards")
text(10,2.7, "downwards Transplant")

# no variance for species, lm ok
fit0 <- lmer(abs.st.traitmean ~ 1 + (1|species), d6)
fit1 <- lmer(abs.st.traitmean ~ clim.diff + (1|species), d6)
fit2 <- lmer(abs.st.traitmean ~ clim.diff + I(clim.diff^2) + (1|species), d6)
modsel(list(fit0, fit1, fit2),1000)
summary(fit0)


mtext("Temperature difference", 1, outer=TRUE, line=1, cex=2)
#mtext("Standardized traitmean", 2, outer=TRUE, line=1, cex=2)
dev.off()
