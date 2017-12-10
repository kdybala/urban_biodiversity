##----READ ME----
# Script for examining correlations among variables & spatial scales,
# scaling variables, and creating nested variables

##----AUTHOR----
## Kristen Dybala, kdybala@gmail.com


dat = read.csv('data/final_compiled_data.csv')
# 16,949 observations of 51 variables

##----CORRELATIONS AMONG VARIABLES----
# remove repeat observations at same unique coordinates
tmp = dat[-which(duplicated(dat$locationID)), ] #4094

# custom correlation function
customcor = function(data, mapping, method='spearman', sizeRange=c(3,4), ...) {
  ct = suppressWarnings(cor.test(x=eval(mapping$x, data), 
                                 y=eval(mapping$y, data), method=method))
  rt <- round(unname(ct$estimate), digits=2)
  cex = max(sizeRange)
  # helper function to calculate a useable size
  percent_of_range <- function(percent, range) {
    percent * diff(range) + min(range, na.rm = TRUE)
  }
  GGally::ggally_text(
    label = as.character(rt), mapping = ggplot2::aes(), 
    size = I(percent_of_range(cex * abs(rt), sizeRange)), ...) 
}

customsmooth = function(data, mapping, ...) {
  ggplot2::ggplot(data, mapping) + 
    ggplot2::geom_smooth(method='gam', formula = y ~ s(x))
}

GGally::ggpairs(tmp[, c(
  'popd2500', 'housed2500', 'imperv2500', 'seminatural2500', 'natdiv2500', 
  'yrbuilt2500', 'poc2500', 'hs2500', 'income2500diff')], 
  upper=list(continuous=customcor), lower=list(continuous=customsmooth))
ggplot2::ggsave('figs/correlations2500.png')

GGally::ggpairs(tmp[, c(
  'popd1000', 'housed1000', 'imperv1000', 'seminatural1000', 'natdiv1000', 
  'yrbuilt1000', 'poc1000', 'hs1000', 'income1000diff')], 
  upper=list(continuous=customcor), lower=list(continuous=customsmooth))
ggplot2::ggsave('figs/correlations1000.png')

GGally::ggpairs(tmp[, c(
  'popd500', 'housed500', 'imperv500', 'seminatural500', 'natdiv500', 
  'yrbuilt500', 'poc500', 'hs500', 'income500diff')], 
  upper=list(continuous=customcor), lower=list(continuous=customsmooth))
ggplot2::ggsave('figs/correlations500.png')

## DECISION: drop housed, imperv, seminatural, hs

dat = dat[, -grep('housed|imperv|seminatural|hs', names(dat))]

## sqrt-transformed population density
dat$sqrt.popd500 = sqrt(dat$popd500)
dat$sqrt.popd1000 = sqrt(dat$popd1000)
dat$sqrt.popd2500 = sqrt(dat$popd2500)


##----CORRELATIONS AMONG SCALES----
# remove repeat observations at same unique coordinates
tmp = dat[-which(duplicated(dat$locationID)), ] #4094

## all highly correlated
GGally::ggpairs(tmp[, c('popd500', 'popd1000', 'popd2500')],
                upper=list(continuous=customcor), 
                lower=list(continuous=customsmooth))
GGally::ggpairs(tmp[, c('natdiv500', 'natdiv1000', 'natdiv2500')],
                upper=list(continuous=customcor), 
                lower=list(continuous=customsmooth))
GGally::ggpairs(tmp[, c('poc500', 'poc1000', 'poc2500')],
                upper=list(continuous=customcor), 
                lower=list(continuous=customsmooth))
GGally::ggpairs(tmp[, c('income500diff', 'income1000diff', 'income2500diff')],
                upper=list(continuous=customcor), 
                lower=list(continuous=customsmooth))
GGally::ggpairs(tmp[, c('yrbuilt500', 'yrbuilt1000', 'yrbuilt2500')],
                upper=list(continuous=customcor), 
                lower=list(continuous=customsmooth))

##----SCALE VARIABLES----
## find scaling parameters (mean/sd) on 500m scale
i = c(9:12, 17:18, grep('\\D500', names(dat)))
means = colMeans(dat[, i], na.rm=T)
vars = apply(dat[, i], 2, sd, na.rm=T)
save(means, vars, file='analysis/data/scaling.RData')
# save original means/sd for later back-transformation

cols1000 = grep('1000', names(dat))
cols2500 = grep('2500', names(dat))

## scale according to means & 2*sd on 500m scale
dat.scaled = dat
dat.scaled[, i] = scale(dat.scaled[, i], center=means, scale=2*vars)
dat.scaled[, cols1000] = scale(
  dat.scaled[, cols1000], center = means[7:length(means)], 
  scale = 2 * vars[7:length(means)])
dat.scaled[, cols2500] = scale(
  dat.scaled[, cols2500], center = means[7:length(means)], 
  scale = 2 * vars[7:length(means)])
summary(dat.scaled)


##----CREATE NESTED VARIABLES----
dat.scaled$popd500new = dat.scaled$popd500 - dat.scaled$popd1000
dat.scaled$popd1000new = dat.scaled$popd1000 - dat.scaled$popd2500

dat.scaled$natdiv500new = dat.scaled$natdiv500 - dat.scaled$natdiv1000
dat.scaled$natdiv1000new = dat.scaled$natdiv1000 - dat.scaled$natdiv2500

dat.scaled$poc500new = dat.scaled$poc500 - dat.scaled$poc1000
dat.scaled$poc1000new = dat.scaled$poc1000 - dat.scaled$poc2500

dat.scaled$income500diffnew = dat.scaled$income500diff - dat.scaled$income1000diff
dat.scaled$income1000diffnew = dat.scaled$income1000diff - dat.scaled$income2500diff

dat.scaled$yrbuilt500new = dat.scaled$yrbuilt500 - dat.scaled$yrbuilt1000
dat.scaled$yrbuilt1000new = dat.scaled$yrbuilt1000 - dat.scaled$yrbuilt2500

## check correlations among nested vars
tmp = dat.scaled[-which(duplicated(dat.scaled$locationID)),]
GGally::ggpairs(tmp[,c('popd500new', 'popd1000new', 'popd2500')],
                upper=list(continuous=customcor), 
                lower=list(continuous=customsmooth))
GGally::ggpairs(tmp[, c('natdiv500new', 'natdiv1000new', 'natdiv2500')],
                upper=list(continuous=customcor), 
                lower=list(continuous=customsmooth))
GGally::ggpairs(tmp[, c('poc500new', 'poc1000new', 'poc2500')],
                upper=list(continuous=customcor), 
                lower=list(continuous=customsmooth))
GGally::ggpairs(tmp[, c('income500diffnew', 'income1000diffnew', 'income2500diff')],
                upper=list(continuous=customcor), 
                lower=list(continuous=customsmooth))
GGally::ggpairs(tmp[, c('yrbuilt500new', 'yrbuilt1000new', 'yrbuilt2500')],
                upper=list(continuous=customcor), 
                lower=list(continuous=customsmooth))

write.csv(dat.scaled, 'analysis/data/final_transformed_data.csv', row.names=F)
