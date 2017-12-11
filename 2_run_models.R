##----README----
# Script torun GLMMs on full data set

##----AUTHOR----
## Kristen Dybala, kdybala@gmail.com


dat = read.csv('data/final_transformed_data.csv') 
summary(dat)
str(dat)
dat$geoid10 = as.factor(dat$geoid10)
dat$locationID = as.factor(dat$locationID)
dat$ecoregionID = as.factor(dat$ecoregionID)
dat$geoid10Tract = as.factor(dat$geoid10Tract)

##----baseline model----
## all candidate variables included on all scales
## quadratic effect on largest scale, nested smaller scales
baseline = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control=lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))
save(baseline, file='output/baseline_model.RData')

##----examine spatial autocorrelation----
## use spline correlogram to examine correlation (and 95% bootstrap confidence
## intervals) in residuals of baseline model (based on Zuur chapter 21)
## Note: spline correlogram used here estimates spatial dependence as a 
## continuous function of distance rather than binning into pre-determined 
## distance classes.

library(rgdal)

p = SpatialPointsDataFrame(
  coords=data.frame(long=dat$longitude, lat=dat$latitude), 
  proj4string=CRS('+proj=longlat'), data=dat)
#convert to UTM so distances are in km
p = spTransform(p, CRS('+proj=utm +zone=10 +datum=WGS84 +units=km')) 
coords = as.data.frame(coordinates(p))
rm(p)

baseline.res = residuals(baseline, type='pearson')
rm(baseline, dat)

memory.limit(size=16000) #requires 64-bit R
sp.corr <- ncf::spline.correlog(
  x = coords$long, y = coords$lat, z = baseline.res, xmax=5)
summary(sp.corr)
save(sp.corr, file='output/autocorr_results.RData')

custom.plot = function (x, xmax=0, text = TRUE, conversion = 1, ...) {
  obj <- x
  xmax <- ifelse(xmax == 0, obj$max.distance, xmax)
  x <- round(obj$real$x, 1)
  y <- round(obj$real$y, 2)
  xul <- round(quantile(obj$boot$boot.summary$x.intercept, 
                        probs = c(0.025, 0.975), na.rm = TRUE), 1)
  yul <- round(quantile(obj$boot$boot.summary$y.intercept, 
                        probs = c(0.025, 0.975), na.rm = TRUE), 2)
  plot(obj$real$predicted$x/conversion, obj$real$predicted$y, 
       xlim = c(0, xmax/conversion), type = "n", ylab = "Correlation", ...)
  if (!is.null(obj$boot$boot.summary$predicted$x)) {
    polygon(c(obj$boot$boot.summary$predicted$x/conversion, 
              rev(obj$boot$boot.summary$predicted$x/conversion)),
            c(obj$boot$boot.summary$predicted$y["0.025", ], 
              rev(obj$boot$boot.summary$predicted$y["0.975", ])),
            border=NA, col='gray90')
    lines(obj$boot$boot.summary$predicted$x / conversion, 
          obj$boot$boot.summary$predicted$y["0.5", ], col='black')
  }
  lines(c(0, max(obj$real$predicted$x)), c(0, 0), lty='dashed')
}
custom.plot(sp.corr, xmax=5, conversion=1, xlab='Distance (km)', las=1, 
            ylim=c(-.15, .15))
custom.plot(sp.corr, xmax=1, conversion=1, xlab='Distance (km)', las=1, 
            ylim=c(-.15, .15))
custom.plot(sp.corr, xmax=.5, conversion=1, xlab='Distance (km)', las=1, 
            ylim=c(-.15, .15))


##----ROUND 1----
load('output/baseline_model.RData')

namedList <- function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)),deparse)[-1]
  if (is.null(nm <- names(L))) nm <- snm
  if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
  setNames(L,nm)
}

# start with baseline, dropping 2500 scale for each variable, one at a time

##----mod1a: natdiv----
mod1a1 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1a2 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1a3 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv1000 + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1a4 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv500 + I(natdiv500^2) + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1a5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv500 + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1a6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

(mod1atab = AICcmodavg::aictab(
  namedList(baseline, mod1a1, mod1a2, mod1a3, mod1a4, mod1a5, mod1a6)))

# keep form a2 only; no support for a6 dropping natdiv
list = ls(pattern = 'mod1')
list = c('baseline',list)
save(list=list, file='output/models_round1.RData')

##----mod1b: popd----
mod1b1 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + popd1000new + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1b2 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd1000 + I(sqrt.popd1000^2) + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
    family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1b3 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd1000 + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew +
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
    family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1b4 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) +
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1b5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd500 + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1b6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optCtrl=list(optimizer='bobyqa', maxfun=100000)))

(mod1btab = AICcmodavg::aictab(
  namedList(baseline,mod1b1,mod1b2,mod1b3,mod1b4,mod1b5,mod1b6)))

## 1b: keep forms b5, b4; no support for b6 dropping popd
list = ls(pattern = 'mod1')
list = c('baseline',list)
save(list=list, file='output/models_round1.RData')

##----mod1c: POC----
mod1c1 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc2500 + poc1000new + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1c2 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc1000 + I(poc1000^2) + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1c3 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc1000 + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1c4 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc500 + I(poc500^2) +
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1c5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) +
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc500 + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1c6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

(mod1ctab = AICcmodavg::aictab(
  namedList(baseline, mod1c1, mod1c2, mod1c3, mod1c4, mod1c5, mod1c6)))

## 1c: keep forms c5, c4, c3, c6
list = ls(pattern = 'mod1')
list = c('baseline',list)
save(list=list, file='output/models_round1.RData')

##----mod1d: income----
mod1d1 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income2500diff + income1000diffnew + income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1d2 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income1000diff + I(income1000diff^2) + income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1d3 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income1000diff + income500diffnew + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1d4 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income500diff + I(income500diff^2) +
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1d5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income500diff + 
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1d6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    
    yrbuilt2500 + I(yrbuilt2500^2) + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

(mod1dtab = AICcmodavg::aictab(
  namedList(baseline, mod1d1, mod1d2, mod1d3, mod1d4, mod1d5, mod1d6)))

## 1d: keep forms d4, d2, d6
list = ls(pattern = 'mod1')
list = c('baseline',list)
save(list=list, file='output/models_round1.RData')

##----mod1e: yrbuilt----
mod1e1 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt2500 + yrbuilt1000new + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1e2 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt1000 + I(yrbuilt1000^2) + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1e3 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt1000 + yrbuilt500new + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1e4 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt500 + I(yrbuilt500^2) +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id),
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev +
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    yrbuilt500 + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod1e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv2500 + I(natdiv2500^2) + natdiv1000new + natdiv500new + 
    sqrt.popd2500 + I(sqrt.popd2500^2) + popd1000new + popd500new + 
    poc2500 + I(poc2500^2) + poc1000new + poc500new + 
    income2500diff + I(income2500diff^2) + income1000diffnew + 
    income500diffnew + 
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

(mod1etab = AICcmodavg::aictab(
  namedList(baseline, mod1e1, mod1e2, mod1e3, mod1e4, mod1e5, mod1e6)))

## 1e: keep forms e6, e5
list = ls(pattern = 'mod1')
list = c('baseline',list)
save(list=list, file='output/models_round1.RData')

## SUMMARY:
## 1a: keep form a2 only; no support for a6 dropping natdiv
## 1b: keep forms b5, b4; no support for b6 dropping popd
## 1c: keep forms c5, c4, c3, c6
## 1d: keep forms d4, d2, d6
## 1e: keep forms e6, e5

rm(list = ls(pattern = 'mod1'))

##----ROUND 2----
## all combinations of retained model forms within each set of round 1 models

##----a2b5----
##----a2b5c5----
##----a2b5c5d4----
mod2a2b5c5d4e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    poc500 + 
    income500diff + I(income500diff^2) +
    
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b5c5d4e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    poc500 + 
    income500diff + I(income500diff^2) + 
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

##----a2b5c5d2----
mod2a2b5c5d2e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    poc500 + 
    income1000diff + I(income1000diff^2) + income500diffnew + 
    
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b5c5d2e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    poc500 + 
    income1000diff + I(income1000diff^2) + income500diffnew + 
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

##----a2b5c5d6----
mod2a2b5c5d6e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    poc500 + 
    
    
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b5c5d6e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    poc500 + 
    
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

save(list=ls(pattern='mod2'), file='output/models_round2.RData')

##----a2b5c4----
##----a2b5c4d4----
mod2a2b5c4d4e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    poc500 + I(poc500^2) + 
    income500diff + I(income500diff^2) +
    
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b5c4d4e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    poc500 + I(poc500^2) + 
    income500diff + I(income500diff^2) + 
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

##----a2b5c4d2----
mod2a2b5c4d2e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    poc500 + I(poc500^2) + 
    income1000diff + I(income1000diff^2) + income500diffnew + 
    
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b5c4d2e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    poc500 + I(poc500^2) + 
    income1000diff + I(income1000diff^2) + income500diffnew + 
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

##----a2b5c4d6----
mod2a2b5c4d6e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    poc500 + I(poc500^2) + 

    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b5c4d6e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    poc500 + I(poc500^2) + 
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

save(list=ls(pattern='mod2'), file='output/models_round2.RData')

##----a2b5c3----
##----a2b5c3d4----
mod2a2b5c3d4e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    poc1000 + poc500new + 
    income500diff + I(income500diff^2) +
    
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b5c3d4e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    poc1000 + poc500new + 
    income500diff + I(income500diff^2) + 
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

##----a2b5c3d2----
mod2a2b5c3d2e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    poc1000 + poc500new + 
    income1000diff + I(income1000diff^2) + income500diffnew + 
    
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b5c3d2e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    poc1000 + poc500new +
    income1000diff + I(income1000diff^2) + income500diffnew + 
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

##----a2b5c3d6----
mod2a2b5c3d6e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    poc1000 + poc500new + 
    
    
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b5c3d6e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    poc1000 + poc500new + 
    
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

save(list=ls(pattern='mod2'), file='output/models_round2.RData')

##----a2b5c6----
##----a2b5c6d4----
mod2a2b5c6d4e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    
    income500diff + I(income500diff^2) +

    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b5c6d4e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    
    income500diff + I(income500diff^2) + 
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

##----a2b5c6d2----
mod2a2b5c6d2e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    
    income1000diff + I(income1000diff^2) + income500diffnew + 

    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b5c6d2e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    
    income1000diff + I(income1000diff^2) + income500diffnew + 
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

##----a2b5c6d6----
mod2a2b5c6d6e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    
    
    
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b5c6d6e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + 
    
    
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

save(list=ls(pattern='mod2'), file='output/models_round2.RData')

##----a2b4----
##----a2b4c5----
##----a2b4c5d4----
mod2a2b4c5d4e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    poc500 + 
    income500diff + I(income500diff^2) +

    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b4c5d4e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    poc500 + 
    income500diff + I(income500diff^2) +
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

##----a2b4c5d2----
mod2a2b4c5d2e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    poc500 + 
    income1000diff + I(income1000diff^2) + income500diffnew + 

    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b4c5d2e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    poc500 + 
    income1000diff + I(income1000diff^2) + income500diffnew + 
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

##----a2b4c5d6----
mod2a2b4c5d6e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    poc500 + 
    
    
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b4c5d6e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    poc500 + 
    
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

save(list=ls(pattern='mod2'), file='output/models_round2.RData')

##----a2b4c4----
##----a2b4c4d4----
mod2a2b4c4d4e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    poc500 + I(poc500^2) + 
    income500diff + I(income500diff^2) +
    
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b4c4d4e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    poc500 + I(poc500^2) + 
    income500diff + I(income500diff^2) + 
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

##----a2b4c4d2----
mod2a2b4c4d2e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    poc500 + I(poc500^2) + 
    income1000diff + I(income1000diff^2) + income500diffnew + 

    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b4c4d2e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    poc500 + I(poc500^2) + 
    income1000diff + I(income1000diff^2) + income500diffnew + 
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

##----a2b4c4d6----
mod2a2b4c4d6e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    poc500 + I(poc500^2) + 
    
    
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b4c4d6e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    poc500 + I(poc500^2) + 
    
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

save(list=ls(pattern='mod2'), file='output/models_round2.RData')

##----a2b4c3----
##----a2b4c3d4----
mod2a2b4c3d4e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    poc1000 + poc500new + 
    income500diff + I(income500diff^2) +

    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b4c3d4e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    poc1000 + poc500new + 
    income500diff + I(income500diff^2) + 
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

##----a2b4c3d2----
mod2a2b4c3d2e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    poc1000 + poc500new + 
    income1000diff + I(income1000diff^2) + income500diffnew + 
    
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b4c3d2e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    poc1000 + poc500new +
    income1000diff + I(income1000diff^2) + income500diffnew + 
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

##----a2b4c3d6----
mod2a2b4c3d6e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    poc1000 + poc500new + 
    
    
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b4c3d6e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    poc1000 + poc500new + 
    
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optCtrl=list(maxfun=100000)))

save(list=ls(pattern='mod2'), file='output/models_round2.RData')

##----a2b4c6----
##----a2b4c6d4----
mod2a2b4c6d4e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    
    income500diff + I(income500diff^2) +

    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b4c6d4e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    
    income500diff + I(income500diff^2) + 
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

##----a2b4c6d2----
mod2a2b4c6d2e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    
    income1000diff + I(income1000diff^2) + income500diffnew + 
    
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b4c6d2e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    
    income1000diff + I(income1000diff^2) + income500diffnew + 
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

##----a2b4c6d6----
mod2a2b4c6d6e6 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    
    
    
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))

mod2a2b4c6d6e5 = lme4::glmer(
  richness ~ hours + I(hours^2) + dist + I(dist^2) + yday + I(yday^2) + 
    time + I(time^2) + temp + I(temp^2) + elev + 
    natdiv1000 + I(natdiv1000^2) + natdiv500new + 
    sqrt.popd500 + I(sqrt.popd500^2) + 
    
    
    yrbuilt500 +
    (1|observer_id) + (1|locationID) + (1|geoid10Tract) + (1|geoid10) + 
    (1|ecoregionID) + (1|sample_id), 
  family = poisson, data = dat, 
  control = lme4::glmerControl(optimizer='bobyqa', optCtrl=list(maxfun=100000)))


##----all round 2 models----
save(list=ls(pattern='mod2'), file='output/models_round2.RData')
