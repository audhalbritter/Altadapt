# RELATIONSHIP AND FIGURE: MAT VS. ELEVATION
neudat <- DAT.T[DAT.T$StudyID %in% dat2$StudyID,]
neudat$region <- DAT.S[match(neudat$ID,DAT.S$ID),c(26)]

pdf ("Fig_MAT_ELEV.pdf", width=10, height=7, pointsize=5, onefile=TRUE, paper="special")
par(mfrow=c(2,4), mar=c(2,2,0.5,0.5), mgp=c(3,1,0), oma=c(3,3,0.5, 0.5), cex=2)
plot(bio_1~elev, neudat, subset=neudat$region=="Africa", ylim=c(-5,30))
box(); axis(1); axis(2)
abline(lm(bio_1~elev, neudat, subset=neudat$region=="Africa"))
text(95,29,"a) Africa")

plot(bio_1 ~ elev, neudat, subset=neudat$region=="Asia", ylim=c(-5,30), axes=FALSE)
box(); axis(1); axis(2, labels=FALSE)
abline(lm(bio_1~elev, neudat, subset=neudat$region=="Asia"))
text(2500,29,"b) Asia")

plot(bio_1 ~ elev, neudat, subset=neudat$region=="Europe", ylim=c(-5,30), axes=FALSE)
box(); axis(1); axis(2, labels=FALSE)
abline(lm(bio_1~elev, neudat, subset=neudat$region=="Europe"))
text(700,29,"c) Europe")

plot(bio_1 ~ elev, neudat, subset=neudat$region=="MAmerica", ylim=c(-5,30), axes=FALSE)
box(); axis(1); axis(2, labels=FALSE)
abline(lm(bio_1~elev, neudat, subset=neudat$region=="MAmerica"))
text(1600,29,"d) Middle America")

plot(bio_1 ~ elev, neudat, subset=neudat$region=="NAmerica", ylim=c(-10,30), axes=FALSE)
box(); axis(1); axis(2)
abline(lm(bio_1~elev, neudat, subset=neudat$region=="NAmerica"))
text(1400,29,"e) North America")

plot(bio_1 ~ elev, neudat, subset=neudat$region=="Oceania", ylim=c(-5,30), axes=FALSE)
box(); axis(1); axis(2, labels=FALSE)
abline(lm(bio_1~elev, neudat, subset=neudat$region=="Oceania"))
text(500,29,"f) Oceania")

plot(bio_1 ~ elev, neudat, subset=neudat$region=="SAmerica", ylim=c(-5,30), axes=FALSE)
box(); axis(1); axis(2, labels=FALSE)
abline(lm(bio_1~elev, neudat, subset=neudat$region=="SAmerica"))
text(1000,29,"g) South America")

mtext("Elevation in m a.s.l", side=1, line=1.5, outer=TRUE, cex=2.5)
mtext("Mean annual temperature in °C", side=2, line=1.5, outer=TRUE, cex=2.5)
dev.off()

# plot all together
by(neudat, neudat$region, function(z){
  fit <- lm(bio_1~elev, z)
  plot(bio_1~elev, z, ylim=c(0,30))
  text(500,29, unique(z$region))
  abline(fit)
  summary(fit)
  cor(z$bio_1, z$elev)
})

# correlation for Asia and Oceania (remove NAs)
asia <- subset(neudat, neudat$region=="Asia")
asia <- subset(asia, asia$bio_1!="NA")
cor(asia$bio_1, asia$elev)


# Compare linear and quadratic models (R^2)
compare.lin.quad <- by(DAT.T, DAT.T$studyunit, function(x){
  standardized.traitmean <- (x$traitmean - mean(x$traitmean, na.rm=TRUE))/sd(x$traitmean, na.rm=TRUE)
  # get r2 for comparison of linear and quadratic regression
  fit1 <- tryCatch(lm(standardized.traitmean ~ bio_1, x), error=function(e) NA)
  fit2 <- tryCatch(lm(standardized.traitmean ~ bio_1+I(bio_1^2), x), error=function(e) NA)
  p.value <- tryCatch(anova(fit1, fit2)$"Pr(>F)"[2], error=function(e) NA)
  sample.size <- tryCatch((length(resid(lm(standardized.traitmean~bio_1,x)))), error=function(e) NA)
  aic.value <- tryCatch(AIC(fit1)-AIC(fit2), error=function(e) NA)
  cbind(p.value, sample.size, aic.value)
})
comp <- t(sapply(compare.lin.quad, I))
colnames(comp) <- c("p.value", "sample.size", "aic.value")
#comp2<- cbind(meta.studyunit, comp)
meta.studyunit[,41:43] <- comp[,1:3]
colnames(meta.studyunit)[41:43] <- c("p.value", "sample.size", "aic.value")

comp2 <- meta.studyunit
names(comp2)
plot(comp2$p.value, col=ifelse(comp2$p.value < 0.05,"red", "black"))
quad <- subset(comp2, comp2$p.value < 0.05)
# in 66 out of 1354 cases, the quadratic model is better!!! 4.9%

# If I use only dat2 data set then it is 5.8%
compare.lin.quad <- by(neudat, neudat$studyunit, function(x){
  standardized.traitmean <- (x$traitmean - mean(x$traitmean, na.rm=TRUE))/sd(x$traitmean, na.rm=TRUE)
  # get r2 for comparison of linear and quadratic regression
  fit1 <- tryCatch(lm(standardized.traitmean ~ bio_1, x), error=function(e) NA)
  fit2 <- tryCatch(lm(standardized.traitmean ~ bio_1+I(bio_1^2), x), error=function(e) NA)
  p.value <- tryCatch(anova(fit1, fit2)$"Pr(>F)"[2], error=function(e) NA)
  sample.size <- tryCatch((length(resid(lm(standardized.traitmean~bio_1,x)))), error=function(e) NA)
  aic.value <- tryCatch(AIC(fit1)-AIC(fit2), error=function(e) NA)
  cbind(p.value, sample.size, aic.value)
})
comp <- t(sapply(compare.lin.quad, I))
colnames(comp) <- c("p.value", "sample.size", "aic.value")
comp2 <- as.data.frame(comp)
quad <- subset(comp2, comp2$p.value < 0.05)
dim(comp2)
dim(quad)

100*53/908




# CLIMATE DISTANCE
neudat <- DAT.T[DAT.T$StudyID %in% dat2$StudyID,]
neudat$bio_1_site <- DAT.S[match(neudat$ID,DAT.S$ID),c(34)]
neudat$studysite <- DAT.S[match(neudat$ID,DAT.S$ID),c(9)]
clim <- neudat

st.traitmean <- by(neudat, neudat$studyunit, function(x){
  standardized.traitmean <- (x$traitmean - mean(x$traitmean, na.rm=TRUE))/sd(x$traitmean, na.rm=TRUE)
})
st.traitmean <- as.vector(unlist(st.traitmean))
clim <- cbind(neudat, st.traitmean)


clim <- subset(clim, clim$bio_1!="NA")
clim <- subset(clim, clim$studysite=="Garden")
clim$clim.diff <- clim$bio_1_site - clim$bio_1 # calculate difference in temperature between site and pop

clim$abs.st.traitmean <- (clim$st.traitmean) # traitmean
clim$abs.clim.diff <- (clim$clim.diff) # climate difference

d5 <- subset(clim, clim$TRAIT_CAT1 == "FITNESS")
plot(abs.st.traitmean ~ abs.clim.diff, d5, ylab="Abs. standardized traitmean", xlab="Difference in MAT of site and population")
fit0 <- lm(abs.st.traitmean ~ 1, d5)
fit1 <- lm(abs.st.traitmean ~ abs.clim.diff, d5)
d5$abs.clim.diff.2 <- d5$abs.clim.diff^2
fit2 <- lm(abs.st.traitmean ~ abs.clim.diff + abs.clim.diff.2, d5)
anova(fit1, fit2)
summary(fit1)
abline(fit1, col="grey")
f <- summary(fit1)$fstatistic
text(9,2, bquote("F"["1, 303"] ~ "=" ~ .(f[1])))


### Data summary
dat2 <- droplevels(dat2)
length(unique(dat2$StudyID)) # nr of studies
length(unique(dat2$SiteID)) # sites
sort(unique(dat2$SiteID)) # sites
unique(dat2$species) # sp
unique(dat2$growthform)
unique(dat2$family) # family
table(dat2$species, dat2$growthform)
table(dat2$species, dat2$intro)
table(dat2$StudyID, dat2$region)
table(dat2$StudyID, dat2$studysite)
table(dat2$StudyID, dat2$study_type)
mean(dat2$r.temp, na.rm=TRUE)
sd(dat2$r.temp, na.rm=TRUE)
mean(dat2$dist.km, na.rm=TRUE)
sd(dat2$dist.km, na.rm=TRUE)

# Population data
range(dat2$r.temp, na.rm=TRUE) # some observations in Rice1991Oeco are < 1, but only because some observations are missing.
sort(dat2$r.temp)
sort(dat2$r.elev) # some observations in Vita2013Oeco are < 150m, but only because some observations are missing.
sort(dat2$dist.km)




trans <- subset(dat2, dat2$study_type=="reciprocal.transplant")
head(trans)
range(trans$dist.km)
unique(trans$StudyID)



