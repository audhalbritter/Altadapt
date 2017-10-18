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
    mean.seasonality <- mean(x$bio_4)
    distance.dd <- max(dist(cbind(x$lat.pop, x$long.pop)))
    distance.km <- max(distm(cbind(x$long.pop, x$lat.pop), fun=distHaversine), na.rm=TRUE)/1000
    mean.distance.km <- mean(distm(cbind(x$long.pop, x$lat.pop), fun=distHaversine), na.rm=TRUE)/1000
    cbind(mean.elev.pop, r.elev.pop, mean.lat.pop, mean.long.pop, mean.temp.pop, r.temp.pop, mean.seasonality, distance.dd, distance.km, mean.distance.km)
  })
  new.means2<-t(sapply(new.means,I))
  colnames(new.means2) <- c("mean.elev", "r.elev", "mean.lat", "mean.long", "mean.temp", "r.temp", "mean.seasonality", "dist.dd", "dist.km", "mean.dist.km")
  return(new.means2)
}
#ddd <- get.means(DAT.T)

get.means2 <- function(dd){
  res <- dd %>% 
    group_by(studyunit) %>% 
    mutate(mean.elev = mean(elev)) %>% 
    mutate(r.elev = max(elev) - min(elev)) %>%
    mutate(mean.lat = mean(lat.pop)) %>%
    mutate(mean.long = mean(long.pop)) %>%
    mutate(mean.temp = mean(bio_1)) %>%
    mutate(r.temp = max(bio_1) - min(bio_1)) %>% 
    mutate(mean.seasonality = mean(bio_4)) %>% 
    # fix these
    distance.dd <- max(dist(cbind(x$lat.pop, x$long.pop)))
    distance.km <- max(distm(cbind(x$long.pop, x$lat.pop), fun=distHaversine), na.rm=TRUE)/1000
    mean.distance.km <- mean(distm(cbind(x$long.pop, x$lat.pop), fun=distHaversine), na.rm=TRUE)/1000

  return(res)
}



# CALCULATE SLOPES
## standardise trait values within studyunits to mean 0 sd 1 and extract regression coefs of trait on mean temperature (bio_1)
get.slopes2 <- function(dd){
  dd <- dd %>% 
    group_by(studyunit) %>% 
    filter(!is.na(traitmean), !is.na(bio_1)) %>%
    mutate(n = n()) %>% 
    filter(n > 2) %>% 
    mutate(std.traitmean = scale(traitmean)) %>% 
    group_by(n, add = TRUE) %>% # add n for sample size
    do({
      # fit model
      fit <- lm(std.traitmean ~ bio_1, data = .)
      
      # extract coefficients
      res <- tidy(fit) %>% 
        filter(term == "bio_1") %>% 
        select(-term, -std.error, -statistic)
  })
  return(dd)
}
 

 
get.slopes <- function(dd){
  result <- by(dd, dd$studyunit, function(x){
    standardized.traitmean <- (x$traitmean - mean(x$traitmean, na.rm=TRUE))/sd(x$traitmean, na.rm=TRUE) #scale(traitmean)
    
    # get slopes for bio_1, bio_1^2 and elev
    slopes <- tryCatch(as.numeric(coef(lm(standardized.traitmean ~ bio_1, x))[2]), error=function(e)NA)
    slopes2 <- tryCatch(as.numeric(coef(lm(standardized.traitmean ~ bio_1+I(bio_1^2), x))[2]), error=function(e)NA)
    slopes.elev <- tryCatch(as.numeric(coef(lm(standardized.traitmean ~ elev, x))[2]), error=function(e)NA)
    slopes.elev2 <- tryCatch(as.numeric(coef(lm(standardized.traitmean ~ elev+I(elev^2), x))[2]), error=function(e)NA)
    
    # Pearson Correlation Coeff
    correlation <- cor(standardized.traitmean, bio_1, method = c("pearson"))
    
    # get sample size
    sample.size <- tryCatch((length(resid(lm(standardized.traitmean~bio_1,x)))), error=function(e) NA)
    
    # get p-values
    p.slope <- tryCatch(as.numeric(anova(lm(standardized.traitmean ~ bio_1, x))[1,5]), error=function(e) NA)
    p.slope2 <- tryCatch(as.numeric(anova(lm(standardized.traitmean ~ bio_1+I(bio_1^2), x))[1,5]), error=function(e) NA)
    p.slope.elev <- tryCatch(as.numeric(anova(lm(standardized.traitmean ~ elev, x))[1,5]), error=function(e) NA)
    p.slope.elev2 <- tryCatch(as.numeric(anova(lm(standardized.traitmean ~ elev+I(elev^2), x))[1,5]), error=function(e) NA)
    
    # calculate correlation between elev and MAT
    cor.elev.bio_1 <- tryCatch(cor(x$elev, x$bio_1), error=function(e) NA)
    
    cbind(correlation, slopes, slopes2, slopes.elev, slopes.elev2, sample.size, p.slope, p.slope2, p.slope.elev, p.slope.elev2, cor.elev.bio_1)
  })
  result2 <- t(sapply(result, I))
  colnames(result2) <- c("correlation", "slopes", "slopes2", "slopes.elev", "slopes.elev2", "sample.size", "p.slope", "p.slope2", "p.slope.elev", "p.slope.elev2", "cor.elev.bio_1")
  return(result2)
}





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


### Calculate collinearity (Variance inflation factor)
vif.mer <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  nam <- names(fixef(fit))
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}




# Function for model selection and averageing mixed effects models
# Input: full model
# printTable = TRUE: prints dredge table
# threshold in percent the cumulative sum that should be used for the model averaging
ModelAverage <- function(mod, printFullTable, print95Table, percent.thresh){
  model.set <- dredge(mod, fixed = "r.temp", rank = "AICc")
  mm <- data.frame(model.set)
  mm$cumsum <- cumsum(mm$weight)
  mm95 <- mm %>% filter(cumsum < percent.thresh)
  if(printFullTable == TRUE){
    print(model.set)
  }
  if(print95Table == TRUE){
    print(mm95)
  }
  averaged.model <- model.avg(model.set, cumsum(weight) <= percent.thresh)
  res <- data.frame(summary(averaged.model)$coefmat.full)
  
  # Importance
  imp <- data.frame(importance(averaged.model))
  imp <- imp %>%
    rownames_to_column(var = "Category") %>% 
    setNames(., c("Category", "Importance")) %>%
    mutate(Category = plyr::mapvalues(Category, c("r.temp", "mean.temp", "mean.seasonality", "dist.km", "growthform", "intro", "breed"), c("T Range", "Mean T", "T Seasonality", "Distance", "Growthform", "Introduction Status", "Breeding System"))) %>% 
    mutate(Importance = round(Importance, 2))
  
  res2 <- res %>% 
    rownames_to_column(var = "Variable") %>% 
    setNames(., c("Variable", "Estimate", "StError", "AdjSE", "Zvalue", "Pvalue")) %>% 
    select(-AdjSE) %>% 
    mutate(Category = Variable) %>% 
    mutate(Category = plyr::mapvalues(Category, c("r.temp", "mean.temp", "mean.seasonality", "dist.km", "growthformgrass", "growthformtree", "intronative", "breedmixed_mating", "breedoutcrossing"), c("T Range", "Mean T", "T Seasonality", "Distance", "Growthform", "Growthform", "Introduction Status", "Breeding System", "Breeding System"))) %>% 
    mutate(Variable = plyr::mapvalues(Variable, c("r.temp", "mean.temp", "mean.seasonality", "dist.km", "growthformgrass", "growthformtree", "intronative", "breedmixed_mating", "breedoutcrossing"), c("T Range", "Mean T", "T Seasonality", "Distance", "Grass", "Tree", "Native", "Mixed Mating", "Outcrossing"))) %>% 
    mutate(CI.low = Estimate - 1.96 * StError) %>% 
    mutate(CI.high = Estimate + 1.96 * StError) %>% 
    mutate(Estimate = round(Estimate, 2), CI.low = round(CI.low, 2), CI.high = round(CI.high, 2), Zvalue = round(Zvalue, 2), Pvalue = round(Pvalue, 3)) %>% 
    mutate(CI = paste(CI.low, CI.high, sep = " - ")) %>% 
    select(Category, Variable, Estimate, CI, Zvalue, Pvalue) %>% 
    #filter(!Variable == "(Intercept)") %>% # remove intercept
    # add importance for ranking
    left_join(imp, by = "Category") %>% 
    mutate(Variable = factor(Variable, levels = Variable))
  return(res2)
}


# Function to print 95 cumulative selection table
# Input: full model
# threshold in percent the cumulative sum that should be used for the model averaging
ModelSelectionTable <- function(mod, percent.thresh){
  model.set <- dredge(mod, fixed = "r.temp", rank = "AICc")
  mm <- data.frame(model.set)
  mm$cumsum <- cumsum(mm$weight)
  res2 <- mm %>% 
    filter(cumsum < percent.thresh) %>% 
    mutate(breed = ifelse(!is.na(breed), "BS", "")) %>% 
    mutate(growthform = ifelse(!is.na(growthform), "GF", "")) %>% 
    mutate(intro = ifelse(!is.na(intro), "IS", "")) %>% 
    mutate(mean.seasonality = ifelse(!is.na(mean.seasonality), "MS", "")) %>% 
    mutate(mean.temp = ifelse(!is.na(mean.temp), "MT", "")) %>% 
    mutate(Model = paste(breed, growthform, intro, mean.seasonality, mean.temp, sep = "+")) %>% 
    mutate(delta = round(delta, 2), weight = round(weight, 3)) %>% 
    select(Model, df, delta, weight) %>% 
    mutate(Model = gsub("^\\+*|\\+*$", "", Model)) %>% # removes + at start and end
    mutate(Model = gsub("\\+{2,}", "+", Model)) # removes 2-x + and replace by +; does not touch +
  return(res2)
}



# Function for plotting output from ModelAverage
# Input: output form ModelAverage
PlotEstimates <- function(dat){
  dat %>% 
    filter(!Variable == "(Intercept)") %>%# remove intercept
    ggplot(aes(x = Variable, y = Estimate, ymin = CI.low, ymax = CI.high)) + 
    geom_point() +
    labs(x = "") +
    geom_hline(yintercept = 0, color = "grey") +
    geom_errorbar(width=0.25) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}



PlotEstimates3 <- function(mod, percent.thresh){
  model.set <- dredge(mod, fixed = "r.temp", rank = "AICc")
  mm <- data.frame(model.set)
  mm$cumsum <- cumsum(mm$weight)
  mm95 <- mm %>% filter(cumsum < percent.thresh)

  averaged.model <- model.avg(model.set, cumsum(weight) <= percent.thresh)
  res <- data.frame(summary(averaged.model)$coefmat.full)
  
  res2 <- res %>% 
    rownames_to_column(var = "Variable") %>% 
    setNames(., c("Variable", "Estimate", "StError", "AdjSE", "Zvalue", "Pvalue")) %>% 
    select(-AdjSE) %>% 
    mutate(Category = Variable) %>% 
    mutate(Category = plyr::mapvalues(Category, c("r.temp", "mean.temp", "mean.seasonality", "dist.km", "growthformgrass", "growthformtree", "intronative", "breedmixed_mating", "breedoutcrossing"), c("T Range", "Mean T", "T Seasonality", "Distance", "Growthform", "Growthform", "Introduction Status", "Breeding System", "Breeding System"))) %>% 
    mutate(Variable = plyr::mapvalues(Variable, c("r.temp", "mean.temp", "mean.seasonality", "dist.km", "growthformgrass", "growthformtree", "intronative", "breedmixed_mating", "breedoutcrossing"), c("T Range", "Mean T", "T Seasonality", "Distance", "Grass", "Tree", "Native", "Mixed Mating", "Outcrossing"))) %>% 
    mutate(CI.low = Estimate - 1.96 * StError) %>% 
    mutate(CI.high = Estimate + 1.96 * StError) %>% 
    mutate(Estimate = round(Estimate, 2), CI.low = round(CI.low, 2), CI.high = round(CI.high, 2), Zvalue = round(Zvalue, 2), Pvalue = round(Pvalue, 3)) %>% 
    select(Category, Variable, Estimate, CI.low, CI.high) %>% 
    filter(!Variable == "(Intercept)") %>% # remove intercept
    mutate(Variable = factor(Variable, levels = c("Mean T", "T Seasonalyit", "Distance", "Native", "Grass", "Tree", "Mixed Mating", "Outcrossing")))
  
  res2 %>% 
    filter(!Variable == "T Range") %>% 
    ggplot(aes(y = Variable, x = Estimate, xmin = CI.low, xmax = CI.high)) + 
    geom_point() +
    labs(x = "Parameter estimate", y = "") +
    geom_vline(xintercept = 0, color = "grey", linetype = "dashed") +
    geom_errorbarh(height = 0) 
}
