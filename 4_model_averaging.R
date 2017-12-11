
master = read.csv('data/final_compiled_data.csv')
dat = read.csv('data/final_transformed_data.csv')
load('data/scaling.RData')

## ----model list----
load('output/models_round2.RData')
modlist = lapply(ls(pattern='mod2'), get)
names(modlist) = ls(pattern='mod2')
set = AICcmodavg::confset(modlist) ## models in the "confidence set"
nrow(set$table) #32
conflist = mget(x=as.character(set$table$Modnames))
topmod = AICcmodavg::aictab(modlist)$Modnames[1]


## ----FIND BEST CASE EFFORT & BIOGEOGRAPHIC VALUES----
## for model averaging, hold survey/biogeography values at best case, 
## all others at mean (=0 for centered/scaled vars)

## find best effort vars

## hours
newdat = data.frame(
  hours=c(min(dat$hours), seq(-0.60, 2.05, .05), max(dat$hours)), 
  dist=0, yday=0, time=0, temp=0, elev=0,
  natdiv2500 = 0, natdiv1000 = 0, natdiv500=0, 
    natdiv1000new = 0, natdiv500new = 0,
  sqrt.popd2500 = 0, sqrt.popd1000 = 0, sqrt.popd500 = 0, 
    popd1000new = 0, popd500new = 0,
  poc2500 = 0, poc1000 = 0, poc500 = 0, poc1000new=0, poc500new = 0,
  income2500diff = 0, income1000diff = 0, income500diff = 0, 
    income1000diffnew = 0, income500diffnew = 0,
  yrbuilt2500 = 0, yrbuilt1000=0, yrbuilt500=0, 
    yrbuilt1000new=0, yrbuilt500new=0)

x = predict(get(as.character(topmod)), newdat, re.form=NA)
std.hours = newdat$hours[which(x==max(x))]
(std.hours * (2 * vars[names(vars)=='hours']) + 
    means[names(means)=='hours'])
#2.2 hours

## dist
newdat = data.frame(
  hours=0, dist=c(min(dat$dist), seq(-1, 0.8, .05), max(dat$dist)), 
  yday=0, time=0, temp=0, elev=0,
  natdiv2500 = 0, natdiv1000 = 0, natdiv500=0, 
    natdiv1000new = 0, natdiv500new = 0,
  sqrt.popd2500 = 0, sqrt.popd1000 = 0, sqrt.popd500 = 0, 
    popd1000new = 0, popd500new = 0,
  poc2500 = 0, poc1000 = 0, poc500 = 0, poc1000new=0, poc500new = 0,
  income2500diff = 0, income1000diff = 0, income500diff = 0, 
    income1000diffnew = 0, income500diffnew = 0,
  yrbuilt2500 = 0, yrbuilt1000=0, yrbuilt500=0, 
    yrbuilt1000new=0, yrbuilt500new=0)

x = predict(get(as.character(topmod)), newdat, re.form=NA)
std.dist = newdat$dist[which(x==max(x))]
(std.dist * (2 * vars[names(vars)=='dist']) + means[names(means)=='dist']) 
#0.807 km

## yday
newdat = data.frame(
  hours=0, dist=0, yday=c(min(dat$yday), seq(-.65, .95, .05), max(dat$yday)), 
  time=0, temp=0, elev=0,
  natdiv2500 = 0, natdiv1000 = 0, natdiv500=0, 
    natdiv1000new = 0, natdiv500new = 0,
  sqrt.popd2500 = 0, sqrt.popd1000 = 0, sqrt.popd500 = 0, 
    popd1000new = 0, popd500new = 0,
  poc2500 = 0, poc1000 = 0, poc500 = 0, poc1000new=0, poc500new = 0,
  income2500diff = 0, income1000diff = 0, income500diff = 0, 
    income1000diffnew = 0, income500diffnew = 0,
  yrbuilt2500 = 0, yrbuilt1000=0, yrbuilt500=0, 
    yrbuilt1000new=0, yrbuilt500new=0)

x = predict(get(as.character(topmod)), newdat, re.form=NA)
std.yday = newdat$yday[which(x==max(x))]
(std.yday * (2 * vars[names(vars)=='yday']) + means[names(means)=='yday']) 
#day 120

## time
newdat = data.frame(
  hours=0, dist=0, yday=0, 
  time=c(min(dat$time), seq(-0.75, 1, .05), max(dat$time)), 
  temp=0, elev=0,
  natdiv2500 = 0, natdiv1000 = 0, natdiv500=0, 
    natdiv1000new = 0, natdiv500new = 0,
  sqrt.popd2500 = 0, sqrt.popd1000 = 0, sqrt.popd500 = 0, 
    popd1000new = 0, popd500new = 0,
  poc2500 = 0, poc1000 = 0, poc500 = 0, poc1000new=0, poc500new = 0,
  income2500diff = 0, income1000diff = 0, income500diff = 0, 
    income1000diffnew = 0, income500diffnew = 0,
  yrbuilt2500 = 0, yrbuilt1000=0, yrbuilt500=0, 
    yrbuilt1000new=0, yrbuilt500new=0)

x = predict(get(as.character(topmod)), newdat, re.form=NA)
std.time = newdat$time[which(x==max(x))]
(std.time * (2 * vars[names(vars)=='time']) + means[names(means)=='time'])
# hour 5 (a.k.a. 5am)

## temp
newdat = data.frame(
  hours=0, dist=0, yday=0, time=0, 
  temp=c(min(dat$temp), seq(-1.1, 1.1, .05), max(dat$temp)), 
  elev=0,
  natdiv2500 = 0, natdiv1000 = 0, natdiv500=0, 
    natdiv1000new = 0, natdiv500new = 0,
  sqrt.popd2500 = 0, sqrt.popd1000 = 0, sqrt.popd500 = 0, 
    popd1000new = 0, popd500new = 0,
  poc2500 = 0, poc1000 = 0, poc500 = 0, poc1000new=0, poc500new = 0,
  income2500diff = 0, income1000diff = 0, income500diff = 0, 
    income1000diffnew = 0, income500diffnew = 0,
  yrbuilt2500 = 0, yrbuilt1000=0, yrbuilt500=0, 
    yrbuilt1000new=0, yrbuilt500new=0)

x = predict(get(as.character(topmod)), newdat, re.form=NA)
std.temp = newdat$temp[which(x==max(x))]
(std.temp * (2 * vars[names(vars)=='temp']) + means[names(means)=='temp'])
# 11.935

## elev
newdat = data.frame(
  hours = 0, dist = 0, yday = 0, time = 0, temp = 0, 
  elev = c(min(dat$elev), seq(-0.25, 2.65, .05), max(dat$elev)),
  natdiv2500 = 0, natdiv1000 = 0, natdiv500=0, 
    natdiv1000new = 0, natdiv500new = 0,
  sqrt.popd2500 = 0, sqrt.popd1000 = 0, sqrt.popd500 = 0, 
    popd1000new = 0, popd500new = 0,
  poc2500 = 0, poc1000 = 0, poc500 = 0, poc1000new=0, poc500new = 0,
  income2500diff = 0, income1000diff = 0, income500diff = 0, 
    income1000diffnew = 0, income500diffnew = 0,
  yrbuilt2500 = 0, yrbuilt1000=0, yrbuilt500=0, 
    yrbuilt1000new=0, yrbuilt500new=0)

x = predict(get(as.character(topmod)), newdat, re.form=NA)
std.elev = newdat$elev[which(x==max(x))]
(std.elev * (2 * vars[names(vars)=='elev']) + means[names(means)=='elev'])
# -1.87m (minimum)


## ----MODEL AVERAGING: EFFORT & BIOGEOGRAPHIC VARS----
custom.modavgPred <- function(cand.set, newdata, ...) {
  aic <- AICcmodavg::aictab(cand.set = cand.set, sort = FALSE)
  ##extract fitted values and se
  fit <- as.data.frame(lapply(X = cand.set, FUN = function(i) {
    AICcmodavg::predictSE(i, se.fit = TRUE, newdata = newdata, 
                          type = 'link', level = 0)
  }))
  newdata$link = rowSums(sweep(
    fit[, seq(1, ncol(fit), 2)], 2, aic$AICcWt, '*'))
  newdata$link.se = sqrt(rowSums(sweep(
    fit[, seq(2,ncol(fit), 2)]^2 + 
      (fit[, seq(1, ncol(fit), 2)] - newdata$link)^2, 2, aic$AICcWt, '*')
  ))
  if (fam.link.mer(cand.set[[1]])$link=='log') {
    newdata$pred = exp(newdata$link)
    newdata$lcl = exp(newdata$link-2*newdata$link.se)
    newdata$ucl = exp(newdata$link+2*newdata$link.se)
  }
  return(newdata)
}

##----HOURS----
newdat = data.frame(
  hours=c(min(dat$hours), seq(-0.55, 2.05, .05), max(dat$hours)), 
  dist=std.dist, yday=std.yday, time=std.time, temp=std.temp, elev=std.elev,
  natdiv2500 = 0, natdiv1000 = 0, natdiv500=0, 
    natdiv1000new = 0, natdiv500new = 0,
  sqrt.popd2500 = 0, sqrt.popd1000 = 0, sqrt.popd500 = 0, 
    popd1000new = 0, popd500new = 0,
  poc2500 = 0, poc1000 = 0, poc500 = 0, poc1000new=0, poc500new = 0,
  income2500diff = 0, income1000diff = 0, income500diff = 0, 
    income1000diffnew = 0, income500diffnew = 0,
  yrbuilt2500 = 0, yrbuilt1000=0, yrbuilt500=0, 
    yrbuilt1000new=0, yrbuilt500new=0)

system.time(pred.hours <- custom.modavgPred(cand.set=conflist, newdata=newdat))
pred.hours$var = 'hours'
#back-transform and convert to minutes
pred.hours$value = (pred.hours$hours * (2 * vars[names(vars)=='hours']) + 
                      means[names(means)=='hours'])

##----DIST----
newdat = data.frame(
  hours=std.hours, 
  dist=c(min(dat$dist), seq(-1,0.8,.05), max(dat$dist)), 
  yday=std.yday, time=std.time, temp=std.temp, elev=std.elev,
  natdiv2500 = 0, natdiv1000 = 0, natdiv500=0, 
    natdiv1000new = 0, natdiv500new = 0,
  sqrt.popd2500 = 0, sqrt.popd1000 = 0, sqrt.popd500 = 0, 
    popd1000new = 0, popd500new = 0,
  poc2500 = 0, poc1000 = 0, poc500 = 0, poc1000new=0, poc500new = 0,
  income2500diff = 0, income1000diff = 0, income500diff = 0, 
    income1000diffnew = 0, income500diffnew = 0,
  yrbuilt2500 = 0, yrbuilt1000=0, yrbuilt500=0, 
    yrbuilt1000new=0, yrbuilt500new=0)

system.time(pred.dist <- custom.modavgPred(cand.set=conflist, newdata=newdat))
pred.dist$var = 'dist'
pred.dist$value = pred.dist$dist * (2 * vars[names(vars)=='dist']) + 
  means[names(means)=='dist']

##----YDAY----
newdat = data.frame(
  hours=std.hours, dist=std.dist, 
  yday=c(min(dat$yday), seq(- .65, .95, .05), max(dat$yday)), 
  time=std.time, temp=std.temp, elev=std.elev,
  natdiv2500 = 0, natdiv1000 = 0, natdiv500=0, 
    natdiv1000new = 0, natdiv500new = 0,
  sqrt.popd2500 = 0, sqrt.popd1000 = 0, sqrt.popd500 = 0, 
    popd1000new = 0, popd500new = 0,
  poc2500 = 0, poc1000 = 0, poc500 = 0, poc1000new=0, poc500new = 0,
  income2500diff = 0, income1000diff = 0, income500diff = 0, 
    income1000diffnew = 0, income500diffnew = 0,
  yrbuilt2500 = 0, yrbuilt1000=0, yrbuilt500=0, 
    yrbuilt1000new=0, yrbuilt500new=0)

system.time(pred.yday <- custom.modavgPred(cand.set=conflist, newdata=newdat))
pred.yday$var = 'yday'
pred.yday$value = pred.yday$yday * (2 * vars[names(vars)=='yday']) + 
  means[names(means)=='yday']


##----TIME----
newdat = data.frame(
  hours=std.hours, dist=std.dist, yday=std.yday, 
  time=c(min(dat$time), seq(-0.75, 1, .05), max(dat$time)), 
  temp=std.temp, elev=std.elev,
  natdiv2500 = 0, natdiv1000 = 0, natdiv500=0, 
    natdiv1000new = 0, natdiv500new = 0,
  sqrt.popd2500 = 0, sqrt.popd1000 = 0, sqrt.popd500 = 0, 
    popd1000new = 0, popd500new = 0,
  poc2500 = 0, poc1000 = 0, poc500 = 0, poc1000new=0, poc500new = 0,
  income2500diff = 0, income1000diff = 0, income500diff = 0, 
    income1000diffnew = 0, income500diffnew = 0,
  yrbuilt2500 = 0, yrbuilt1000=0, yrbuilt500=0, 
    yrbuilt1000new=0, yrbuilt500new=0)

system.time(pred.time <- custom.modavgPred(cand.set=conflist, newdata=newdat))
pred.time$var = 'time'
pred.time$value = pred.time$time * (2 * vars[names(vars)=='time']) + 
  means[names(means)=='time']

##----TEMP----
newdat = data.frame(
  hours=std.hours, dist=std.dist, yday=std.yday, time=std.time, 
  temp=c(min(dat$temp), seq(-1.1, 1.1, .05), max(dat$temp)), elev=std.elev,
  natdiv2500 = 0, natdiv1000 = 0, natdiv500=0, 
    natdiv1000new = 0, natdiv500new = 0,
  sqrt.popd2500 = 0, sqrt.popd1000 = 0, sqrt.popd500 = 0, 
    popd1000new = 0, popd500new = 0,
  poc2500 = 0, poc1000 = 0, poc500 = 0, poc1000new=0, poc500new = 0,
  income2500diff = 0, income1000diff = 0, income500diff = 0, 
    income1000diffnew = 0, income500diffnew = 0,
  yrbuilt2500 = 0, yrbuilt1000=0, yrbuilt500=0, 
    yrbuilt1000new=0, yrbuilt500new=0)

system.time(pred.temp <- custom.modavgPred(cand.set=conflist, newdata=newdat))
pred.temp$var = 'temp'
pred.temp$value = pred.temp$temp * (2 * vars[names(vars)=='temp']) + 
  means[names(means)=='temp']

##----ELEV----
newdat = data.frame(
  hours=std.hours, dist=std.dist, yday=std.yday, time=std.time, temp=std.temp, 
  elev=c(min(dat$elev), seq(-0.25, 2.65, .05), max(dat$elev)),
  natdiv2500 = 0, natdiv1000 = 0, natdiv500=0, 
    natdiv1000new = 0, natdiv500new = 0,
  sqrt.popd2500 = 0, sqrt.popd1000 = 0, sqrt.popd500 = 0, 
    popd1000new = 0, popd500new = 0,
  poc2500 = 0, poc1000 = 0, poc500 = 0, poc1000new=0, poc500new = 0,
  income2500diff = 0, income1000diff = 0, income500diff = 0, 
    income1000diffnew = 0, income500diffnew = 0,
  yrbuilt2500 = 0, yrbuilt1000=0, yrbuilt500=0, 
    yrbuilt1000new=0, yrbuilt500new=0)

system.time(pred.elev <- custom.modavgPred(cand.set=conflist, newdata=newdat))
pred.elev$var = 'elev'
pred.elev$value = pred.elev$elev * (2 * vars[names(vars)=='elev']) + 
  means[names(means)=='elev']

## ----COMPILE EFFORT & BIOGEOGRAPHIC VARS----
predeffort = do.call(rbind.data.frame, lapply(ls(pattern='pred.'), get))
write.csv(predeffort, 'output/predicted_values_effortvars.csv', row.names=F)

## ----CANDIDATE VARIABLES----
## use same range of values for all 3 scales

##----NATDIV----
summary(dat[,c('natdiv500','natdiv1000','natdiv2500')])
newdat = data.frame(
  hours=std.hours, dist=std.dist, yday=std.yday, time=std.time, temp=std.temp, 
  elev=std.elev,
  natdiv2500 = c(min(dat$natdiv1000), seq(-0.4,2.5,.05), max(dat$natdiv1000)),
  natdiv1000 = c(min(dat$natdiv1000), seq(-0.4,2.5,.05), max(dat$natdiv1000)), 
  natdiv500 = c(min(dat$natdiv1000), seq(-0.4,2.5,.05), max(dat$natdiv1000)), 
  natdiv1000new = 0, natdiv500new = 0,
  sqrt.popd2500 = 0, sqrt.popd1000 = 0, sqrt.popd500 = 0, 
    popd1000new = 0, popd500new = 0,
  poc2500 = 0, poc1000 = 0, poc500 = 0, poc1000new=0, poc500new = 0,
  income2500diff = 0, income1000diff = 0, income500diff = 0, 
    income1000diffnew = 0, income500diffnew = 0,
  yrbuilt2500 = 0, yrbuilt1000=0, yrbuilt500=0, 
    yrbuilt1000new=0, yrbuilt500new=0)

system.time(pred.natdiv <- custom.modavgPred(cand.set=conflist, newdata=newdat))
pred.natdiv$var = 'natdiv'
pred.natdiv$value = pred.natdiv$natdiv1000 * 
  (2 * vars[names(vars)=='natdiv500']) + means[names(means)=='natdiv500']


##----POPD----
## popd500 = popd overall (only one model in candidate set)
summary(dat[,c('sqrt.popd2500','sqrt.popd1000','sqrt.popd500')])
newdat = data.frame(
  hours=std.hours, dist=std.dist, yday=std.yday, time=std.time, temp=std.temp, 
  elev=std.elev,
  natdiv2500 = 0, natdiv1000 = 0, natdiv500=0, 
    natdiv1000new = 0, natdiv500new = 0,
  sqrt.popd2500 = c(min(dat$sqrt.popd500), seq(-0.85,3.8,.05), 
                    max(dat$sqrt.popd500)), 
  sqrt.popd1000 = c(min(dat$sqrt.popd500), seq(-0.85,3.8,.05), 
                    max(dat$sqrt.popd500)), 
  sqrt.popd500 = c(min(dat$sqrt.popd500), seq(-0.85,3.8,.05), 
                   max(dat$sqrt.popd500)), 
    popd1000new = 0, popd500new = 0,
  poc2500 = 0, poc1000 = 0, poc500 = 0, poc1000new=0, poc500new = 0,
  income2500diff = 0, income1000diff = 0, income500diff = 0, 
    income1000diffnew = 0, income500diffnew = 0,
  yrbuilt2500 = 0, yrbuilt1000=0, yrbuilt500=0, 
    yrbuilt1000new=0, yrbuilt500new=0)

system.time(pred.popd <- custom.modavgPred(cand.set=conflist, newdata=newdat))
pred.popd$var = 'popd'
pred.popd$value = (pred.popd$sqrt.popd500 * 
                     (2 * vars[names(vars)=='sqrt.popd500']) + 
                     means[names(means)=='sqrt.popd500'])^2

##----poc----
summary(dat[,c('poc2500','poc1000','poc500')])
newdat = data.frame(
  hours=std.hours, dist=std.dist, yday=std.yday, time=std.time, temp=std.temp, 
  elev=std.elev,
  natdiv2500 = 0, natdiv1000 = 0, natdiv500=0, natdiv1000new = 0, 
    natdiv500new = 0,
  sqrt.popd2500 = 0, sqrt.popd1000 = 0, sqrt.popd500 = 0, 
    popd1000new = 0, popd500new = 0,
  poc2500 = c(min(dat$poc500), seq(-0.75,1.2,0.05), max(dat$poc500)), 
  poc1000 = c(min(dat$poc500), seq(-0.75,1.2,0.05), max(dat$poc500)), 
  poc500 = c(min(dat$poc500), seq(-0.75,1.2,0.05), max(dat$poc500)), 
    poc1000new=0, poc500new = 0,
  income2500diff = 0, income1000diff = 0, income500diff = 0, 
    income1000diffnew = 0, income500diffnew = 0,
  yrbuilt2500 = 0, yrbuilt1000=0, yrbuilt500=0, 
    yrbuilt1000new=0, yrbuilt500new=0)

system.time(pred.poc <- custom.modavgPred(cand.set=conflist, newdata=newdat))
pred.poc$var = 'poc'
pred.poc$value = pred.poc$poc500 * (2 * vars[names(vars)=='poc500']) + 
  means[names(means)=='poc500']


##----INCOME----
summary(dat[,c('income2500diff','income1000diff','income500diff')])
newdat = data.frame(
  hours=std.hours, dist=std.dist, yday=std.yday, time=std.time, temp=std.temp, 
  elev=std.elev,
  natdiv2500 = 0, natdiv1000 = 0, natdiv500=0, 
    natdiv1000new = 0, natdiv500new = 0,
  sqrt.popd2500 = 0, sqrt.popd1000 = 0, sqrt.popd500 = 0, 
    popd1000new = 0, popd500new = 0,
  poc2500 = 0, poc1000 = 0, poc500 = 0, poc1000new=0, poc500new = 0,
  income2500diff = c(min(dat$income500diff), seq(-1.3,2.6,.05), 
                     max(dat$income500diff)),
  income1000diff = c(min(dat$income500diff), seq(-1.3,2.6,.05), 
                     max(dat$income500diff)), 
  income500diff = c(min(dat$income500diff), seq(-1.3,2.6,.05), 
                    max(dat$income500diff)), 
    income1000diffnew = 0, income500diffnew = 0,
  yrbuilt2500 = 0, yrbuilt1000=0, yrbuilt500=0, 
    yrbuilt1000new=0, yrbuilt500new=0)

system.time(pred.income <- custom.modavgPred(cand.set=conflist, newdata=newdat))
pred.income$var = 'income'
pred.income$value = pred.income$income500diff * 
  (2 * vars[names(vars)=='income500diff']) + 
  means[names(means)=='income500diff']


##----YRBUILT----
summary(dat[,c('yrbuilt2500','yrbuilt1000','yrbuilt500')])
newdat = data.frame(
  hours=std.hours, dist=std.dist, yday=std.yday, time=std.time, temp=std.temp, 
  elev=std.elev,
  natdiv2500 = 0, natdiv1000 = 0, natdiv500=0, 
    natdiv1000new = 0, natdiv500new = 0,
  sqrt.popd2500 = 0, sqrt.popd1000 = 0, sqrt.popd500 = 0, 
    popd1000new = 0, popd500new = 0,
  poc2500 = 0, poc1000 = 0, poc500 = 0, poc1000new=0, poc500new = 0,
  income2500diff = 0, income1000diff = 0, income500diff = 0, 
    income1000diffnew = 0, income500diffnew = 0,
  yrbuilt2500 = c(min(dat$yrbuilt1000), seq(-0.75,1,.05), max(dat$yrbuilt1000)), 
  yrbuilt1000 = c(min(dat$yrbuilt1000), seq(-0.75,1,.05), max(dat$yrbuilt1000)), 
  yrbuilt500 = c(min(dat$yrbuilt1000), seq(-0.75,1,.05), max(dat$yrbuilt1000)), 
    yrbuilt1000new=0, yrbuilt500new=0)

system.time(pred.yrbuilt <- custom.modavgPred(cand.set=conflist, newdata=newdat))
pred.yrbuilt$var = 'yrbuilt'
pred.yrbuilt$value = pred.yrbuilt$yrbuilt500 * 
  (2 * vars[names(vars)=='yrbuilt500']) + means[names(means)=='yrbuilt500']


## ----COMPILE CANDIDATE VARIABLES----
predvars = do.call(rbind.data.frame, lapply(ls(pattern='pred.'), get))
write.csv(predvars, 'output/predicted_values_candidatevars.csv', row.names=F)


##----POC & INCOME----
## predict richness for a grid of income & poc data to make contour plot
newdat = expand.grid(
  hours=std.hours, dist=std.dist, yday=std.yday, time=std.time, temp=std.temp, 
  elev=std.elev,
  natdiv2500 = 0, natdiv1000 = 0, natdiv500=0, 
    natdiv1000new = 0, natdiv500new = 0,
  sqrt.popd2500 = 0, sqrt.popd1000 = 0, sqrt.popd500 = 0, 
    popd1000new = 0, popd500new = 0,
  poc500 = seq(0, 1, 0.025), poc1000new=0, poc500new = 0,
  income500diff = seq(-60, 160, 5), income1000diffnew = 0, income500diffnew = 0,
  yrbuilt2500 = 0, yrbuilt1000=0, yrbuilt500=0, 
    yrbuilt1000new=0, yrbuilt500new=0)

newdat$poc500 = 
  (newdat$poc500 - means[names(means)=='poc500']) / 
  (2 * vars[names(vars)=='poc500'])
newdat$poc2500 = newdat$poc500
newdat$poc1000 = newdat$poc500

newdat$income500diff = 
  (newdat$income500diff - means[names(means)=='income500diff']) / 
  (2 * vars[names(vars)=='income500diff'])
newdat$income1000diff = newdat$income500diff
newdat$income2500diff = newdat$income500diff

system.time(pred.minin2 <- custom.modavgPred(cand.set=conflist, newdata=newdat))
pred.minin2$var = 'poc-income'
pred.minin2$value.poc = pred.minin2$poc500 * 
  (2 * vars[names(vars)=='poc500']) + means[names(means)=='poc500']
pred.minin2$value.income = pred.minin2$income500diff * 
  (2 * vars[names(vars)=='income500diff']) + means[names(means)=='income500diff']
write.csv(pred.minin2, 
          'output/predicted_values_poc&income_grid.csv', row.names=F)

