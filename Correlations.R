### MODEL WITH ONE FULL MODEL
library("ggplot2")
library("GGally")
library("cowplot")
library("lme4")
library("broom")
library("knitr")


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



### CORRELATIONS AMONG CONTINOUS VARIABLES
dat4 %>% select(abs.lat, mean.temp, mean.seasonality, dist.km) %>% ggpairs()
### absolute latitude is highly correlated with mean.temp and mean.seasonality
### maybe exclude dist.km, because it does not seem so important
ddd <- dat4 %>% 
  select(abs.lat, mean.temp, mean.seasonality, dist.km)

kable(round(cor(ddd, method = "pearson"), 2))


#### COUNT NUMBER OF OBSERVATIONS PER FACTOR LEVEL ####
d1 <- dat4 %>% filter(var == "biomass", !is.na(breed), !is.na(generation)) %>% 
  select(studyunit, region, breed, generation, longevity, growthform, intro, studysite, family) %>%
  gather(key = variable, value = value, -studyunit) %>% 
  arrange(variable) %>% 
  group_by(variable, value) %>% 
  summarise(n = n())

d2 <- dat4 %>% filter(var == "phenology emergence", !is.na(breed), !is.na(generation)) %>% 
  select(studyunit, region, breed, generation, longevity, growthform, intro, studysite, family) %>%
  gather(key = variable, value = value, -studyunit) %>% 
  arrange(variable) %>% 
  group_by(variable, value) %>% 
  summarise(n = n())

d3 <- dat4 %>% filter(var == "phenology duration", !is.na(breed), !is.na(generation)) %>% 
  select(studyunit, region, breed, generation, longevity, growthform, intro, studysite, family) %>%
  gather(key = variable, value = value, -studyunit) %>% 
  arrange(variable) %>% 
  group_by(variable, value) %>% 
  summarise(n = n())

d4 <- dat4 %>% filter(var == "height", !is.na(breed), !is.na(generation)) %>% 
  select(studyunit, region, breed, generation, longevity, growthform, intro, studysite, family) %>%
  gather(key = variable, value = value, -studyunit) %>% 
  arrange(variable) %>% 
  group_by(variable, value) %>% 
  summarise(n = n()) %>% 
  full_join(d1, by = c("variable", "value")) %>% 
  full_join(d2, by = c("variable", "value")) %>% 
  full_join(d3, by = c("variable", "value")) %>% 
  arrange(variable, value) %>% 
  setNames(., c("Variable", "Level", "Height", "Biomass", "Phenology emergence", "Phenology duration"))


d4 %>% print(n = 42)
### not enough data for family, longevity, generation and region
### studysite ok, if climate chamber and greenhouse are merged


# Chi-squared Test of Independence (or Fisher for small sample size)
# breed, growthform and intro
dd <- dat4 %>% filter(var == "height") %>% 
  filter(!is.na(breed))
dd <- droplevels(dd)
#pairs(dd)
tbl  <-  table(dd$growthform, dd$intro) 
tbl                 # the contingency table
chisq.test(tbl, correct = FALSE) 
chisq.test(tbl, correct = TRUE) # correct = TRUE if 2x2 table
fisher.test(tbl)

dd <- dat4 %>% filter(var == "biomass") %>% 
  filter(!is.na(breed))
dd <- droplevels(dd)
tbl  <-  table(dd$breed, dd$growthform) 
tbl                 # the contingency table
chisq.test(tbl, correct = FALSE) 


### ANOVAS: relashion between coninous and categorical variables
# y needs to be changed

# abs.lat
lin_fit.lat <- function(dat) {
  fit <- lm(abs.lat ~ value, dat)
  tidy(fit)
}
res <- dat4 %>% select(abs.lat, region, breed, generation, longevity, growthform, intro, studysite) %>% 
  gather(key = variables, value = value, -abs.lat) %>% 
  mutate(value = factor(value), variables = factor(variables)) %>% 
  group_by(variables) %>% 
  do(lin_fit.lat(.))
res %>% 
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(t.test = round(statistic, 2)) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "< 0.001", round(p.value, 3))) %>% 
  select(variables, term, estimate, t.test, p.value) %>% 
  print(n = 24)

# r.temp
lin_fit.rtemp <- function(dat) {
  fit <- lm(r.temp ~ value, dat)
  tidy(fit)
}
res <- dat4 %>% select(r.temp, breed, growthform, intro) %>% 
  gather(key = variables, value = value, -r.temp) %>% 
  mutate(value = factor(value), variables = factor(variables)) %>% 
  group_by(variables) %>% 
  do(lin_fit.rtemp(.))
res %>% 
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(t.test = round(statistic, 2)) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "< 0.001", round(p.value, 3))) %>% 
  select(variables, term, estimate, t.test, p.value) %>% 
  print(n = 24)


# dist.km
lin_fit.dist <- function(dat) {
  fit <- lm(dist.km ~ value, dat)
  tidy(fit)
}
res <- dat4 %>% select(dist.km, region, breed, generation, longevity, growthform, intro, studysite) %>% 
  gather(key = variables, value = value, -dist.km) %>% 
  mutate(value = factor(value), variables = factor(variables)) %>% 
  group_by(variables) %>% 
  do(lin_fit.dist(.))
res %>% 
  mutate(estimate = round(estimate, 2)) %>% 
  mutate(t.test = round(statistic, 2)) %>% 
  mutate(p.value = ifelse(p.value < 0.001, "< 0.001", round(p.value, 3))) %>% 
  select(variables, term, estimate, t.test, p.value) %>% 
  print(n = 24)
