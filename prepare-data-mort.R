## Restrict to data on numbers of deaths among adults; aggregate the 45-54 and
## 55+ age groups into a single group

data.raw = readRDS('Data/dunn-table2.rds')
data.raw = subset(data.raw, outcome=='death' & !(agegr %in% c('0-4', '5-14')) & cd4cat != 'all')
data.raw$agegr = as.character(data.raw$agegr)
data.raw$agegr[which(data.raw$agegr %in% c('45-54', '55+'))] = '45+'
data.tmp = aggregate(data.raw[,c('x', 'pys')], by=list(data.raw$agegr, data.raw$cd4cat), FUN=sum)

data.mort = list(
  deaths = t(array(data.tmp$x,   c(4,6), dimnames=list(unique(data.raw$agegr), unique(data.raw$cd4cat)))),
  pyears = t(array(data.tmp$pys, c(4,6), dimnames=list(unique(data.raw$agegr), unique(data.raw$cd4cat))))
)

data.mort$deaths = data.mort$deaths[6:1,]
data.mort$pyears = data.mort$pyears[6:1,]

rm(data.raw, data.tmp)
