#####################
# CLIMATE DISTANCE #
#####################

# use only studies from dat2
CD <- DAT.T[DAT.T$StudyID %in% dat2$StudyID,]
# add studysite, elev and bio_1 from the sites
site.CD <- DAT.S %>% select(ID, studysite, bio_1_site, elevation_site)
CD <- CD %>% left_join(site.CD, by = "ID")
head(CD)

# calculate elev and clim.diff between site and pop (only useful for Gardens) and make it a factor (High/Low)
CD$elev.diff <- CD$elevation_site - CD$elev
CD$clim.diff <- CD$bio_1_site - CD$bio_1
CD$rel.elev <- as.factor(ifelse(CD$elev.diff > 0, -1, 1))

CD <- subset(CD, CD$traitmean!="NA")
st.traitmean <- by(CD, CD$studyunit, function(x){
  standardized.traitmean <- (x$traitmean - mean(x$traitmean, na.rm=TRUE))/sd(x$traitmean, na.rm=TRUE)
})
st.traitmean <- as.vector(unlist(st.traitmean))
clim <- cbind(CD, st.traitmean)
head(clim)

clim <- subset(clim, clim$studysite=="Garden")
clim <- subset(clim, clim$elev.diff!="NA")
clim <- subset(clim, clim$clim.diff!="NA")
clim <- subset(clim, clim$st.traitmean!="NA")
clim$abs.st.traitmean <- abs(clim$st.traitmean) # traitmean

d4 <- subset(clim, clim$TRAIT_CAT1 == "SIZE"& clim$TRAIT_CAT2 == "BIOMASS"& TRAIT_CAT3 == "ABOVEGR")
d5 <- subset(clim, clim$TRAIT_CAT1 == "FITNESS")
d6 <- subset(clim, clim$TRAIT_CAT1 == "VITAL_RATES" & clim$TRAIT_CAT2 == "SURVIVAL")
v1 <- c(rep(1, nrow(d4)), rep(2, nrow(d5)), rep(3, nrow(d6)))
d7 <- cbind(rbind(d4,d5,d6),v1)

length(unique(d4$StudyID))
length(unique(d5$StudyID))
length(unique(d6$StudyID))

# Aboveground Biomass
fit0 <- lm(st.traitmean ~ 1, d4)
fit1 <- lm(st.traitmean ~ clim.diff, d4)
fit2 <- lm(st.traitmean ~ clim.diff + I(clim.diff^2), d4)
anova(fit1, fit2)
AIC(fit1, fit2)
anova(fit0,fit1)
AIC(fit0,fit1)
summary(fit1)

plot(st.traitmean ~ clim.diff, d4, ylab="Standardized mean", xlab="", axes=FALSE, main = "Aboveground Biomass")
box(); axis(1, labels=TRUE); axis(2,labels=TRUE)
abline(v=0, col="darkgrey", lty=2)
abline(fit1, col="black")
lines(sort(d4$clim.diff), fitted(fit2)[order(d4$clim.diff)], col='black', lty=2, type='l')
f <- summary(fit1)$fstatistic
text(10,0.3, bquote("F"["1, 252"] ~ "=" ~ .(round(f[1], 2))))

fit0 <- lmer(st.traitmean ~ 1 + (1|species), d4)
fit1 <- lmer(st.traitmean ~ clim.diff + (1|species), d4)
fit2 <- lmer(st.traitmean ~ clim.diff + I(clim.diff^2) + (1|species), d4)
modsel(list(fit0, fit1, fit2),1000)
summary(fit1)

# Count the number of species per studies!
library("tidyr")
library("plyr")
library("dplyr")
d4 %>%
  select(StudyID, species) %>%
  group_by(StudyID) %>%
  summarise(n = n_distinct(species))
# mainly one study, only a few have more


# Reproductive Output
fit0 <- lm(st.traitmean ~ 1, d5)
fit1 <- lm(st.traitmean ~ clim.diff, d5)
fit2 <- lm(st.traitmean ~ clim.diff + I(clim.diff^2), d5)
anova(fit1, fit2)
AIC(fit1, fit2)
anova(fit0,fit1)
AIC(fit0,fit1)
summary(fit1)
plot(st.traitmean ~ clim.diff, d5, ylab="Standardized mean", xlab="", axes=FALSE, main = "Reproductive output")
box(); axis(1, labels=TRUE); axis(2,labels=TRUE)
abline(v=0, col="darkgrey", lty=2)
abline(fit1, col="black")
lines(sort(d5$clim.diff), fitted(fit2)[order(d5$clim.diff)], col='black', lty=2, type='l')
f <- summary(fit1)$fstatistic
text(10,0.3, bquote("F"["1, 252"] ~ "=" ~ .(round(f[1], 2))))
text(-4,3, "upwards")
text(5,3, "downwards Transplant")

fit0 <- lmer(st.traitmean ~ 1 + (1|species), d5)
fit1 <- lmer(st.traitmean ~ clim.diff + (1|species), d5)
fit2 <- lmer(st.traitmean ~ clim.diff + I(clim.diff^2) + (1|species), d5)
modsel(list(fit0, fit1, fit2),1000)
summary(fit2)


### SURVIVAL
fit0 <- lm(st.traitmean ~ 1, d6)
fit1 <- lm(st.traitmean ~ clim.diff, d6)
fit2 <- lm(st.traitmean ~ clim.diff + I(clim.diff^2), d6)
anova(fit1, fit2)
AIC(fit1, fit2)
summary(fit1)
anova(fit0, fit1)
AIC(fit0, fit1)

plot(st.traitmean ~ clim.diff, d6, ylab="Standardized mean", xlab="", axes=FALSE, main = "Survival")
box(); axis(1, labels=TRUE); axis(2,labels=TRUE)
abline(v=0, col="darkgrey", lty=2)
abline(fit1, col="black")
lines(sort(d6$clim.diff), fitted(fit2)[order(d6$clim.diff)], col='black', lty=2, type='l')


hist(sqrt(d6$st.traitmean+10))
fit0 <- lmer(sqrt(st.traitmean+10) ~ 1 + (1|species), d6)
fit1 <- lmer(sqrt(st.traitmean+10) ~ clim.diff + (1|species), d6)
fit2 <- lmer(sqrt(st.traitmean+10) ~ clim.diff + I(clim.diff^2) + (1|species), d6)
modsel(list(fit0, fit1, fit2),1000)
summary(fit2)
fix.check(fit0)
unique(d6$Trait)
hist(sqrt(d6$st.traitmean))
plot(resid(mod01))

### FIG.4 CLIMATE DISTANCE
pdf ("Fig4_Climate Distance.pdf", width=10, height=4, pointsize=5, onefile=TRUE, paper="special")
par(mfrow=c(1,3), mar=c(1,2,2,0.5), mgp=c(3,1,0), oma=c(3,3,1,1), cex=1.7)

fit1 <- lm(st.traitmean ~ clim.diff, d4)
fit2 <- lm(st.traitmean ~ clim.diff + I(clim.diff^2), d4)
plot(st.traitmean ~ clim.diff, d4, ylab="", xlab="", axes=FALSE, main = "a) Aboveground biomass")
box(); axis(1, labels=TRUE); axis(2,labels=TRUE)
abline(v=0, col="darkgrey", lty=2)
abline(fit1, col="black")
lines(sort(d4$clim.diff), fitted(fit2)[order(d4$clim.diff)], col='black', lty=2, type='l')
#f <- summary(fit1)$fstatistic
legend("topright", c("linear", "quadratic"), col=c("black"), lty=c(1,2), bty="n")

fit1 <- lm(st.traitmean ~ clim.diff, d5)
fit2 <- lm(st.traitmean ~ clim.diff + I(clim.diff^2), d5)
f <- summary(fit1)$fstatistic
plot(st.traitmean ~ clim.diff, d5, ylab="", xlab="", axes=FALSE, main = "b) Reproductive output")
box(); axis(1, labels=TRUE); axis(2,labels=TRUE)
abline(v=0, col="darkgrey", lty=2)
abline(fit1, col="black")
lines(sort(d5$clim.diff), fitted(fit2)[order(d5$clim.diff)], col='black', lty=2, type='l')



fit1 <- lm(abs.st.traitmean ~ clim.diff, d6)
fit2 <- lm(abs.st.traitmean ~ clim.diff + I(clim.diff^2), d6)
anova(fit1, fit2)
AIC(fit1, fit2)
summary(fit1)
f <- summary(fit1)$fstatistic
plot(abs.st.traitmean ~ clim.diff, d6, ylab="", xlab="", axes=FALSE, main = "c) Survival")
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
