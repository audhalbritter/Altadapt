#### USA VA EUROPE

hh.us <- dat4 %>% filter(var == "phenology", region == "NAmerica") %>% 
  filter(!is.na(breed)) %>%  # remove NA's
  mutate(Year = factor(Year))

hh.eu <- dat4 %>% filter(var == "phenology", region == "Europe") %>% 
  filter(!is.na(breed)) %>%  # remove NA's
  mutate(Year = factor(Year))

modHeight <- lmer(slopes ~ r.temp + mean.temp + mean.seasonality + dist.km + breed + growthform + (1|species)+(1|StudyID)+(1|Year), data = hh.eu, REML = FALSE, weights = sample.size, na.action = "na.fail")

model.set <- dredge(modHeight, rank = "AICc", extra = "R^2")
model.set
averaged.model <- model.avg(model.set)
res <- data.frame(summary(averaged.model)$coefmat.full)
res

ggplot(dat4, aes(y = slopes, x = abs.lat)) +
  geom_point() +
  facet_wrap(~ region)

