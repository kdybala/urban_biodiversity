
master = read.csv('data/final_compiled_data.csv')

## ----FIG 1. Map of survey locations----
coords = master[,c('latitude','longitude','geoid10','cityname')]
coords = coords[-which(duplicated(coords$geoid10)),] 
#just keep one lat/long per city

citystats = plyr::ddply(master, plyr::.(geoid10, cityname), plyr::summarize, 
                        loc = length(unique(locationID)), 
                        observer = length(unique(observer_id)), 
                        surveys = length(unique(sample_id)))

citystats = merge(citystats,coords, by=c('geoid10','cityname'))

library(ggplot2)
theme_custom = theme_classic() + theme(
    legend.title = element_blank(), legend.position = c(0,0), 
    legend.justification = c(0,0), legend.text = element_text(size=12), 
    legend.key.height = unit(0.5,'cm'), legend.key.width = unit(0.5,'cm'),
    plot.title = element_text(size=20, face='bold', vjust=1),
    axis.text = element_text(size=16, color='black'), 
    axis.title = element_text(size=16, face='bold', vjust=0),
    strip.text = element_text(size=14, face='bold'), 
    strip.background = element_blank(),
    axis.line.x=element_line(color='black'),
    axis.line.y=element_line(color='black'))

theme_remove_all <- theme(
  axis.text=element_blank(), axis.line=element_blank(), 
  axis.ticks=element_blank(), axis.title=element_blank())

map = borders('state', colour='gray50', fill='gray95')

m = ggplot(data=citystats, aes(x=longitude, y=latitude, size=surveys)) + 
  map + coord_fixed(ratio=1.2) + geom_point(alpha=0.5) +
  xlab(NULL) + ylab(NULL) + scale_x_continuous(breaks=NULL) + 
  scale_y_continuous(breaks=NULL) +
  theme_custom + theme(legend.position='none')

plot_top = ggplot(master, aes(longitude)) +
  geom_density(color='gray50', fill='gray95', adjust=4/5) + 
  theme_custom + theme_remove_all

plot_right = ggplot(master, aes(latitude)) +
  geom_density(color='gray50', fill='gray90', adjust=.75) + coord_flip() + 
  theme_custom + theme_remove_all

g = ggplotGrob(m)
g = gtable::gtable_add_cols(g, unit(0.5,"in"))
g = gtable::gtable_add_grob(g, ggplotGrob(plot_right), t=nrow(g), l=ncol(g), b=1)
g = gtable::gtable_add_rows(g, unit(0.5,"in"))
g = gtable::gtable_add_grob(g, ggplotGrob(plot_top), t=nrow(g), l=7, r=1)

png('figs/Figure1_map.png', width=10, height=5, units='in', res=300)
grid::grid.newpage()
grid::grid.draw(g)
dev.off()


## ----FIG 2. PLOT CANDIDATE VARIABLES----
predvars = read.csv('output/predicted_values_candidatevars.csv')
predvars$var = factor(predvars$var, levels=c(
  'popd', 'natdiv', 'poc', 'income', 'yrbuilt'))

dat = master[,c('locationID', 'popd500', 'natdiv1000', 'poc500', 
                'income500diff', 'yrbuilt500')]
dat = dat[-which(duplicated(dat)),]
dat = reshape2::melt(dat, id.vars='locationID', variable.name='var', 
                     value.name='value')
dat$var = gsub('500|1000|diff|\\.med', '', dat$var)
dat$var = factor(dat$var, levels=c(
  'popd', 'natdiv', 'poc', 'income', 'yrbuilt'))

sdat = plyr::ddply(dat, plyr::.(var), plyr::summarize, q25=quantile(value,.25), 
                   median=median(value), q75=quantile(value,.75))
sdat = reshape2::melt(sdat, id.vars=c('var'))

a = ggplot(predvars[predvars$var=='popd',], aes(x=value)) + 
  geom_ribbon(fill='gray40', alpha=0.5, aes(ymin=lcl, ymax=ucl)) + 
  geom_line(aes(y=pred)) + 
  geom_density(data=dat[dat$var=='popd',], aes(y=(..scaled..)*6), adjust=2, 
               fill='gray90', color='gray50') + 
  geom_text(data=sdat[sdat$variable=='median' & sdat$var=='popd',], 
            aes(x=value, y=2), label='|') + ylim(0,40) +
  xlab('Population density (individuals/ha)') + ylab("Species richness") + 
  ggtitle('A') + theme_custom + 
  scale_x_sqrt(limits=c(0,350), breaks=seq(0,300,100), expand=c(0,0))

b = ggplot(predvars[predvars$var=='natdiv',], aes(x=value)) + 
  geom_ribbon(fill='gray40', alpha=0.5, aes(ymin=lcl, ymax=ucl)) + 
  geom_line(aes(y=pred)) + 
  geom_density(data=dat[dat$var=='natdiv',], aes(y=(..scaled..)*6), adjust=2, 
               fill='gray90', color='gray50') + 
  geom_text(data=sdat[sdat$variable=='median' & sdat$var=='natdiv',], 
            aes(x=value, y=2), label='|') + 
  xlab('Diversity of natural land cover types') + ylab(NULL) + theme_custom + 
  ggtitle('B') + 
  theme(axis.line.y=element_blank(), axis.ticks.y=element_blank()) +
  scale_y_continuous(limits=c(0,40), labels=NULL) + 
  scale_x_continuous(limits=c(1,6.5), expand=c(0,0))

c = ggplot(predvars[predvars$var=='poc',], aes(x=value)) + 
  geom_ribbon(fill='gray40', alpha=0.5, aes(ymin=lcl, ymax=ucl)) + 
  geom_line(aes(y=pred)) + 
  geom_density(data=dat[dat$var=='poc',], aes(y=(..scaled..)*6), adjust=2, 
               fill='gray90', color='gray50') + 
  geom_text(data=sdat[sdat$variable=='median' & sdat$var=='poc',], 
            aes(x=value, y=2), label='|') + ylim(0,40) +
  xlab('Proportion people of color') + ylab("Species richness") + 
  theme_custom + ggtitle('C') +
  scale_x_continuous(expand=c(0,0), limits=c(0,1), breaks=seq(0,1,.2))

d = ggplot(predvars[predvars$var=='income',], aes(x=value)) + 
  geom_ribbon(fill='gray40', alpha=0.5, aes(ymin=lcl, ymax=ucl)) + 
  geom_line(aes(y=pred)) + 
  geom_density(data=dat[dat$var=='income',], aes(y=(..scaled..)*6), adjust=2, 
               fill='gray90', color='gray50') + 
  geom_text(data=sdat[sdat$variable=='median' & sdat$var=='income',], 
            aes(x=value, y=2), label='|') + 
  xlab('Relative med household income (thousands)') + ylab(NULL) + 
  theme_custom + ggtitle('D') + 
  scale_x_continuous(limits=c(-65,170), breaks=seq(-50,150,50), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,40), labels=NULL) + 
  theme(axis.line.y=element_blank(), axis.ticks.y=element_blank())

e = ggplot(predvars[predvars$var=='yrbuilt',], aes(x=value)) + 
  geom_ribbon(fill='gray40', alpha=0.5, aes(ymin=lcl, ymax=ucl)) + 
  geom_line(aes(y=pred)) + 
  geom_density(data=dat[dat$var=='yrbuilt',], aes(y=(..scaled..)*6), adjust=2, 
               fill='gray90', color='gray50') + 
  geom_text(data=sdat[sdat$variable=='median' & sdat$var=='yrbuilt',], 
            aes(x=value, y=2), label='|') + ylim(0,40) +
  xlab('Median year built') + ylab("Species richness") + theme_custom + 
  ggtitle('E') + 
  scale_x_continuous(breaks=c(1939,1960,1980,2005), 
                     labels=c('<1939','1960','1980','>2005'))

cowplot::plot_grid(a,b,c,d,e, ncol=2, align='hv')
ggsave('figs/Figure2.png', height=7.5, width=6.5)

## ----FIG 3.----
load('output/models_round1.RData')

round1tab = as.data.frame(
  rbind(mod1atab, mod1btab, mod1ctab, mod1dtab, mod1etab))
round1tab$scale = 'No\neffect'
round1tab$scale[grep('baseline|1$', round1tab$Modnames)] = '2.5 km'
round1tab$scale[grep('2$|3$', round1tab$Modnames)] = '1 km'
round1tab$scale[grep('4$|5$', round1tab$Modnames)] = '0.5 km'
round1tab$scale = factor(round1tab$scale, levels=c(
  'No\neffect', '0.5 km', '1 km', '2.5 km'))

round1tab$form = 'No effect'
round1tab$form[grep('baseline|2$|4$', round1tab$Modnames)] = 'Quadratic'
round1tab$form[grep('1$|3$|5$', round1tab$Modnames)] = 'Linear'
round1tab$form = factor(round1tab$form, levels=c(
  'Quadratic', 'Linear', 'No effect'))

round1tab$variable = NA
round1tab$variable[grep('mod1a', round1tab$Modnames)] = 
  'B. Diversity of natural land cover types'
round1tab$variable[grep('mod1b', round1tab$Modnames)] = 
  'A. Population density'
round1tab$variable[grep('mod1c', round1tab$Modnames)] = 
  'C. Proportion people of color'
round1tab$variable[grep('mod1d', round1tab$Modnames)] = 
  'D. Median household income'
round1tab$variable[grep('mod1e', round1tab$Modnames)] = 
  'E. Median year built'
round1tab$variable[round1tab$Modnames=='baseline'] = c(
  'B. Diversity of natural land cover types','A. Population density',
  'C. Proportion people of color','D. Median household income','E. Median year built')

round1tab = round1tab[order(round1tab$variable, round1tab$scale, round1tab$form),]

theme_custom = theme_classic() +
  theme(legend.title=element_blank(), legend.position=c(0.98,1), 
        legend.justification=c(1,1), 
        legend.text=element_text(size=10),
        axis.text=element_text(size=11, color='black'), 
        axis.title=element_text(size=11, vjust=1, face='plain'), 
        axis.line.y = element_line(color="black"), 
        axis.line.x = element_line(color="black"),
        strip.background = element_blank(), 
        strip.text=element_text(hjust=0, face='bold', size=11),
        panel.spacing = unit(1,'lines'))

ggplot(round1tab, aes(x=scale, y=AICcWt, fill=form)) + geom_col() + 
  ylim(0,1) + xlab(NULL) + ylab('Akaike weight') + 
  scale_fill_manual(values=c('grey80','grey50','grey30'), na.value='grey30') + 
  facet_wrap(~variable, scales = 'free_x') + theme_custom
ggsave('figs/Figure3.png', height=6, width=10)


## ----FIG S1. PLOT EFFORT & BIOGEOGRAPHIC VARIABLES----
predeffort = read.csv('output/predicted_values_effortvars.csv')
predeffort$var = factor(predeffort$var, levels=c(
  'hours', 'dist', 'yday', 'time', 'temp', 'elev'))

dat = master[,c('locationID','yday','hours','dist','time','elev','temp')]
dat$hours = dat$hours/60
dat$elev = dat$elev-1.873
dat = dat[-which(duplicated(dat)),]
dat = reshape2::melt(
  dat, id.vars='locationID', variable.name='var', value.name='value')
dat$var = factor(dat$var, levels=c('hours','dist','yday','time','temp','elev'))
dat$value[dat$var=='hours'] = dat$value[dat$var=='hours']*60

sdat = plyr::ddply(dat, plyr::.(var), plyr::summarize, q25=quantile(value, .25), 
                   median=median(value), q75=quantile(value, .75))
sdat = reshape2::melt(sdat, id.vars=c('var'))


theme_custom = theme_classic() + theme(
  legend.title = element_blank(), legend.position = c(0.98, 1), 
  legend.justification = c(1, 1), legend.text = element_text(size = 9), 
  axis.text = element_text(size = 10, color='black'), 
  axis.title = element_text(size = 10, vjust = 1, face = 'plain'), 
  axis.line.y = element_line(color = "black"), 
  axis.line.x = element_line(color = "black"), 
  plot.title = element_text(hjust = 0, size = 10), 
  plot.margin = unit(c(0, 0.5, 0.5, 0), 'lines'), 
  strip.background = element_blank(), 
  strip.text = element_text(hjust = 0, face = 'bold', size = 10))


a = ggplot(predeffort[predeffort$var=='hours', ], aes(x=value)) + 
  geom_ribbon(fill='gray40', alpha=0.5, aes(ymin=lcl, ymax=ucl)) + 
  geom_line(aes(y=pred)) + 
  geom_density(data=dat[dat$var=='hours',], aes(y=(..scaled..)*6), adjust=2, 
               fill='gray90', color='gray50') + 
  geom_text(data=sdat[sdat$variable=='median' & sdat$var=='hours',], 
            aes(x=value, y=2), label='|') + 
  xlab('Survey length (hours)') + ylab('Species richness') + theme_custom + 
  ggtitle('A') + scale_y_continuous(limits=c(0,40)) + 
  scale_x_continuous(expand=c(0,0), limits=c(0,3))

b = ggplot(predeffort[predeffort$var=='dist',], aes(x=value)) + 
  geom_ribbon(fill='gray40', alpha=0.5, aes(ymin=lcl, ymax=ucl)) + 
  geom_line(aes(y=pred)) + 
  geom_density(data=dat[dat$var=='dist',], 
               aes(y=(..scaled..)*6), adjust=2, fill='gray90', color='gray50') + 
  geom_text(data=sdat[sdat$variable=='median' & sdat$var=='dist',], 
            aes(x=value, y=2), label='|') + 
  xlab('Survey distance (km)') + ylab(NULL) + theme_custom + 
  theme(axis.line.y=element_blank(), axis.ticks.y=element_blank()) +
  scale_x_continuous(expand=c(0,0), limits=c(0,1), breaks=seq(0,1,0.2)) + 
  scale_y_continuous(limits=c(0,40), labels=NULL) + ggtitle('B') 

c = ggplot(predeffort[predeffort$var=='yday',], aes(x=value)) + 
  geom_ribbon(fill='gray40', alpha=0.5, aes(ymin=lcl, ymax=ucl)) + 
  geom_line(aes(y=pred)) + 
  geom_density(data=dat[dat$var=='yday',], aes(y=(..scaled..)*6), adjust=2, 
               fill='gray90', color='gray50') + 
  geom_text(data=sdat[sdat$variable=='median' & sdat$var=='yday',], 
            aes(x=value, y=2), label='|') + ylim(0,40) +
  xlab('Day of year') + ylab("Species richness") + theme_custom + ggtitle('C') +
  scale_x_continuous(expand=c(0,0))

d = ggplot(predeffort[predeffort$var=='time',], aes(x=value)) + 
  geom_ribbon(fill='gray40', alpha=0.5, aes(ymin=lcl, ymax=ucl)) + 
  geom_line(aes(y=pred)) + 
  geom_density(data=dat[dat$var=='time',], aes(y=(..scaled..)*6), adjust=2, 
               fill='gray90', color='gray50') + 
  geom_text(data=sdat[sdat$variable=='median' & sdat$var=='time',], 
            aes(x=value, y=2), label='|') + 
  xlab('Time of day (hour)') + ylab(NULL) + theme_custom + ggtitle('D') + 
  scale_x_continuous(limits=c(5,20), breaks=seq(5,20,5), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,40), labels=NULL) + 
  theme(axis.line.y=element_blank(), axis.ticks.y=element_blank())

e = ggplot(predeffort[predeffort$var=='temp',], aes(x=value)) + 
  geom_ribbon(fill='gray40', alpha=0.5, aes(ymin=lcl, ymax=ucl)) + 
  geom_line(aes(y=pred)) + 
  geom_density(data=dat[dat$var=='temp',], aes(y=(..scaled..)*6), adjust=2, 
               fill='gray90', color='gray50') + 
  geom_text(data=sdat[sdat$variable=='median' & sdat$var=='temp',], 
            aes(x=value, y=2), label='|') + ylim(0,40) +
  xlab('Mean temp (C)') + ylab('Species richness') + theme_custom + 
  ggtitle('E') + scale_x_continuous(expand=c(0,0))

f = ggplot(predeffort[predeffort$var=='elev',], aes(x=value)) + 
  geom_ribbon(fill='gray40', alpha=0.5, aes(ymin=lcl, ymax=ucl)) + 
  geom_line(aes(y=pred)) + 
  geom_density(data=dat[dat$var=='elev',], aes(y=(..scaled..)*6), adjust=2, 
               fill='gray90', color='gray50') + 
  geom_text(data=sdat[sdat$variable=='median' & sdat$var=='elev',], 
            aes(x=value, y=2), label='|') + 
  xlab('Elevation (m)') + ylab(NULL) + theme_custom + ggtitle('F') + 
  scale_x_continuous(limits=c(-5,2100), breaks=seq(0,2000,500), expand=c(0,0)) +
  scale_y_continuous(limits=c(0,40), labels=NULL) + 
  theme(axis.line.y=element_blank(), axis.ticks.y=element_blank())

cowplot::plot_grid(a,b,c,d,e,f, ncol=2, align='hv')
ggsave('figs/FigureS1.png', height=7.5, width=6.5)


##----FIG S2. POC & INCOME----
pred.pocincome = read.csv('output/predicted_values_poc&income_grid.csv')

tmp = master[-which(duplicated(master$locationID)),]

ggplot(pred.pocincome, aes(value.income, value.poc)) + 
  geom_raster(aes(fill=pred), interpolate=T) + 
  geom_point(data=tmp, aes(income500diff, poc500), size=1) + 
  geom_contour(aes(z=pred), color='white') +
  xlab('Relative median household income (thousands)') + 
  ylab('Proportion people of color') + 
  theme_custom + theme(legend.position='none') + 
  scale_fill_gradientn(colors=terrain.colors(10)) + xlim(-60,160)
ggsave('figs/FigureS2.png', height=4, width=6)

