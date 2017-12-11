
make.table = function(modlist) {
  aic = as.data.frame(AICcmodavg::aictab(modlist))
  
  vars = plyr::rbind.fill(lapply(modlist, function(x) {
    res = as.data.frame(t(lme4::fixef(x)))
  }))
  vars$modname = names(modlist)
  
  tmp = reshape2::melt(vars, id.vars='modname', 
                       variable.name='parm', value.name='est')
  tmp = tmp[-which(is.na(tmp$est)), ] #4168
  tmp$form = 1
  tmp$est = NULL
  tmp$parm = gsub('I\\(', '', tmp$parm)
  tmp$parm = gsub('\\^2\\)', '', tmp$parm)
  tmp$parm = gsub('sqrt.|diff|new', '', tmp$parm)
  
  tmp$parm = factor(tmp$parm, levels=c(
    '(Intercept)','dist','hours','time','yday','temp','elev',
    'popd500', 'popd1000', 'popd2500', 'natdiv500', 'natdiv1000', 'natdiv2500',
    'poc500', 'poc1000', 'poc2500', 'income500', 'income1000', 'income2500',
    'yrbuilt500', 'yrbuilt1000', 'yrbuilt2500'))
  table = reshape2::dcast(tmp, modname~parm, value.var='form', sum, fill=0)
  table[, 2:ncol(table)] = apply(table[, 2:ncol(table)], 2, function(x) {
    x = gsub('1', 'L', x)
    x = gsub('2', 'Q', x)
  })
  table = merge(table, aic[, c(1:4, 6)], by.x='modname', by.y='Modnames', 
                all.x=T, sort=F)
  table = table[order(-table$AICcWt),]
  return(table)
}

##-----SUMMARIZE ROUND 1 MODEL TABLES----

load('output/models_round1.RData')
rm(mod1atab, mod1btab, mod1ctab, mod1dtab, mod1etab)

modlista = lapply(ls(pattern='mod1a|baseline'), get)
names(modlista) = ls(pattern='mod1a|baseline')
tablea = make.table(modlista)
tablea$variable = 'natdiv'

modlistb = lapply(ls(pattern='mod1b|baseline'), get)
names(modlistb) = ls(pattern='mod1b|baseline')
tableb = make.table(modlistb)
tableb$variable = 'popd'

modlistc = lapply(ls(pattern='mod1c|baseline'), get)
names(modlistc) = ls(pattern='mod1c|baseline')
tablec = make.table(modlistc)
tablec$variable = 'poc'

modlistd = lapply(ls(pattern='mod1d|baseline'), get)
names(modlistd) = ls(pattern='mod1d|baseline')
tabled = make.table(modlistd)
tabled$variable = 'income'

modliste = lapply(ls(pattern='mod1e|baseline'), get)
names(modliste) = ls(pattern='mod1e|baseline')
tablee = make.table(modliste)
tablee$variable = 'yrbuilt'

tablea = tablea[, c(1, 12:14, 24:ncol(tablea))]
tableb = tableb[, c(1, 9:11, 24:ncol(tableb))]
tablec = tablec[, c(1, 15:17, 24:ncol(tablec))]
tabled = tabled[, c(1, 18:20, 24:ncol(tabled))]
tablee = tablee[, c(1, 21:23, 24:ncol(tablee))]

names(tablea) = gsub('natdiv', 'scale', names(tablea))
names(tableb) = gsub('popd', 'scale', names(tableb))
names(tablec) = gsub('poc', 'scale', names(tablec))
names(tabled) = gsub('income', 'scale', names(tabled))
names(tablee) = gsub('yrbuilt', 'scale', names(tablee))

table = rbind(tableb, tablea, tablec, tabled, tablee)
table$Delta_AICc = round(table$Delta_AICc, digits=2)
table$AICcWt = round(table$AICcWt, digits=2)
write.csv(table, 'output/results_tableS2.csv', row.names=F)


##----SUMMARIZE ROUND 2 MODELS----

load('output/models_round2.RData')
modlist = lapply(ls(pattern='mod2'), get)
names(modlist) = ls(pattern='mod2')

table = make.table(modlist)
table = table[, c(1, 9:ncol(table))]

table2 = table
table2$Delta_AICc = round(table2$Delta_AICc, digits=2)
table2$AICcWt = round(table2$AICcWt, digits=2)
write.csv(table2, 'output/results_tableS3.csv', row.names=F)

table3 = reshape2::melt(table, id.vars=c(
  'modname', 'K', 'AICc', 'Delta_AICc', 'AICcWt' ))
table3$value[table3$value %in% c('L','Q')] = 1
table3$value = as.numeric(table3$value)
table3$variable = gsub('500|1000|2500', '', table3$variable)
table3$variable = factor(table3$variable, levels=c(
  'popd', 'natdiv', 'poc', 'income', 'yrbuilt'))
table3 = table3[-which(duplicated(table3)),]
table3 = reshape2::dcast(table3, ...~variable, value.var='value', sum)
(wt = colSums(sweep(table3[,6:ncol(table3)], 1, table3$AICcWt, '*')))
## popd & natdiv: 1
## poc: 0.95
## income: 0.91
## yrbuilt: 0.29
write.csv(wt, 'output/results_tableS3_wtsum.csv', row.names=F)


