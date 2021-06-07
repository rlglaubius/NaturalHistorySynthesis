library(data.table)
library(dplyr)

data.long = data.table(read.csv('Data/cascade_cd4init_agegr.csv', stringsAsFactors=FALSE))
data.long$sex = factor(data.long$sex,   levels=c('Male', 'Female'))
data.long$age = factor(data.long$agegr,
                       levels=c('15-19', '20-24', '25-29', '30-34', '35-39', '40-44', '45-49', '50-54', '55-59'),
                       labels=c('15-24', '15-24', '25-34', '25-34', '35-44', '35-44', '45+',    '45+',  '45+'))
data.long$cd4 = factor(data.long$cd4cat, levels=c('500+', '350-499', '200-349', '100-199', '50-99', '<50'))
data.aggr = by(data.long, list(data.long$cd4, data.long$age, data.long$sex), function(obs) {sum(obs$count)})
data.dist = array(data.aggr, dim(data.aggr), dimnames(data.aggr))

rm(data.long, data.aggr)
