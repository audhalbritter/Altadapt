#############################
 ### FUNCTIONS ALTADAPT ####
#############################

### PREPARE THE DATA

# CALCULATE MEANS
# range elev, lat and long and mean temperature (bio_1) and max and mean spatial distance for each studyunit
# input DAT.T
get.means <- function (dd){
  new.means <- by(dd, dd$studyunit, function(x){
    mean.elev.pop <- mean(x$elev)
    r.elev.pop <- max(x$elev)-min(x$elev)
    mean.lat.pop <- mean(x$lat.pop)
    mean.long.pop <- mean(x$long.pop)
    mean.temp.pop <- mean(x$bio_1)
    r.temp.pop <- max(x$bio_1)-min(x$bio_1)
    distance.dd <- max(dist(cbind(x$lat.pop, x$long.pop)))
    distance.km <- max(distm(cbind(x$long.pop, x$lat.pop), fun=distHaversine), na.rm=TRUE)/1000
    mean.distance.km <- mean(distm(cbind(x$long.pop, x$lat.pop), fun=distHaversine), na.rm=TRUE)/1000
    cbind(mean.elev.pop, r.elev.pop, mean.lat.pop, mean.long.pop, mean.temp.pop, r.temp.pop, distance.dd, distance.km, mean.distance.km)
  })
  new.means2<-t(sapply(new.means,I))
  colnames(new.means2) <- c("mean.elev", "r.elev", "mean.lat", "mean.long", "mean.temp", "r.temp", "dist.dd", "dist.km", "mean.dist.km")
  return(new.means2)
}
#ddd <- get.means(DAT.T)


# CALCULATE SLOPES
## standardise trait values within studyunits to mean 0 sd 1 and extract regression coefs of trait on mean temperature (bio_1)
get.slopes <- function(dd, meta){
  result <- by(dd, dd$studyunit, function(x){
    standardized.traitmean <- (x$traitmean - mean(x$traitmean, na.rm=TRUE))/sd(x$traitmean, na.rm=TRUE)
    
    # get slopes for bio_1, bio_1^2 and elev
    slopes <- tryCatch(as.numeric(coef(lm(standardized.traitmean ~ bio_1, x))[2]), error=function(e)NA)
    slopes2 <- tryCatch(as.numeric(coef(lm(standardized.traitmean ~ bio_1+I(bio_1^2), x))[2]), error=function(e)NA)
    slopes.elev <- tryCatch(as.numeric(coef(lm(standardized.traitmean ~ elev, x))[2]), error=function(e)NA)
    slopes.elev2 <- tryCatch(as.numeric(coef(lm(standardized.traitmean ~ elev+I(elev^2), x))[2]), error=function(e)NA)
    
    # get sample size
    sample.size <- tryCatch((length(resid(lm(standardized.traitmean~bio_1,x)))), error=function(e) NA)
    
    # get p-values
    p.slope <- tryCatch(as.numeric(anova(lm(standardized.traitmean ~ bio_1, x))[1,5]), error=function(e) NA)
    p.slope2 <- tryCatch(as.numeric(anova(lm(standardized.traitmean ~ bio_1+I(bio_1^2), x))[1,5]), error=function(e) NA)
    p.slope.elev <- tryCatch(as.numeric(anova(lm(standardized.traitmean ~ elev, x))[1,5]), error=function(e) NA)
    p.slope.elev2 <- tryCatch(as.numeric(anova(lm(standardized.traitmean ~ elev+I(elev^2), x))[1,5]), error=function(e) NA)
    
    # calculate correlation between elev and MAT
    cor.elev.bio_1 <- tryCatch(cor(x$elev, x$bio_1), error=function(e) NA)
    
    cbind(slopes, slopes2, slopes.elev, slopes.elev2, sample.size, p.slope, p.slope2, p.slope.elev, p.slope.elev2, cor.elev.bio_1)
  })
  result2 <- t(sapply(result, I))
  colnames(result2) <- c("slopes", "slopes2", "slopes.elev", "slopes.elev2", "sample.size", "p.slope", "p.slope2", "p.slope.elev", "p.slope.elev2", "cor.elev.bio_1")
  meta[,39:48] <- result2[,1:10]
  colnames(meta)[39:48] <- c("slopes", "slopes2", "slopes.elev", "slopes.elev2", "sample.size", "p.slope", "p.slope2", "p.slope.elev", "p.slope.elev2", "cor.elev.bio_1")
  return(meta)
}
#meta2 <- get.slopes(DAT.T, meta.studyunit)




bars <- function(x,y,z,c){for (k in 1:length(y)) for (i in c(-1, 1)) arrows(x[k], y[k], x[k], y[k]+i*z[k], angle=90, length=0, col=c)}
# x and y are the location of the points, z is the standard error, c the colour
bars2 <- function(x,y,z1,z2,c){for (k in 1:length(y)){ arrows(x[k], y[k], x[k], z1[k], angle=90, length=0, col=c)
                                                       arrows(x[k], y[k], x[k], z2[k], angle=90, length=0, col=c)}}   #for CI


overdisp <- function(mod){
  rp <- resid(mod)
  rdf = as.numeric(attributes(summary(mod))$dims[2] - attr(logLik(mod), "df"))
  print(paste("P=",round(pchisq(sum(rp^2), df=rdf,lower.tail=FALSE),5)))	#notsig, indicating no overdispersion
  print(paste("Dispersion =",round(sum(rp^2)/rdf,3)))}


fix.check <- function(mod){		#function to produce model-checking plots for the fixed effects of an lmer model
  par(mfrow = c(2,2))
  plot(fitted(mod),resid(mod))	#should have no pattern
  abline(h=0)
  print(anova(lm(fitted(mod)~resid(mod))))	#should be non-significant
  qqnorm(resid(mod), ylab="Residuals")		#should be approximately straight line
  qqline(resid(mod))
  plot(density(resid(mod)))					#should be roughly normally distributed
  rug(resid(mod))}


#function for QAICc. NB, phi is the scaling parameter from the quasi-family model. If using e.g. a poisson family, phi=1 and QAICc returns AICc, or AIC if QAICc=FALSE.
QAICc <- function(mod, scale, QAICc = TRUE) {
  ll <- as.numeric(logLik(mod))
  df <- attr(logLik(mod), "df")
  n <- length(resid(mod))
  if (QAICc)
    qaic = as.numeric(-2 * ll/scale + 2 * df + 2 * df * (df + 1)/(n - df - 1))
  else qaic = as.numeric(-2 * ll/scale + 2 * df)
  qaic
}


## code for model selection. First fit mod01, then run this code.
modsel <- function(mods,x){	
  phi=1
  dd <- data.frame(Model=1:length(mods), K=1, QAIC=1)
  for(j in 1:length(mods)){
    dd$K[j] = attr(logLik(mods[[j]]),"df")
    dd$QAIC[j] = QAICc(mods[[j]],phi)
  }
  dd$delta.i <- dd$QAIC - min(dd$QAIC)
  dd <- subset(dd,dd$delta.i<x)
  dd$re.lik <- round(exp(-0.5*dd$delta.i),3)
  sum.aic <- sum(exp(-0.5*dd$delta.i))
  wi <- numeric(0)
  for (i in 1:length(dd$Model)){wi[i] <- round(exp(-0.5*dd$delta.i[i])/sum.aic,3)}; dd$wi<-wi
  print(dds <- dd[order(dd$QAIC), ])
  assign("mstable",dd,envir=.GlobalEnv)
}


# Vertical histogram of the same data
vertical.hist <- function(dat, col.bar){
  hs <- hist(dat, breaks=20, plot=FALSE)
  plot (NA, type='n', axes=FALSE, yaxt='n',
        xlab='', ylab=NA, main=NA, cex.lab=0.7,
        xlim=c(0,max(hs$counts)),
        ylim=c(1,length(hs$counts)))
  #axis (1, cex.axis=0.7)
  arrows(rep(0,length(hs$counts)),1:length(hs$counts),
         hs$counts,1:length(hs$counts),
         length=0,angle=0, lwd=3, col="lightgray")
}
# par(old.par)
# invisible ()


# CREATE PLOT WITH MULTIPLE FIGURES
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
