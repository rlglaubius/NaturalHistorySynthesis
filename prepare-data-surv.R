library(haven)
library(lubridate)

## Sites:
##  3 = Kisesa, Tanzania
##  4 = Masaka, Uganda
##  6 = Rakai, Uganda
##  7 = South African miners
## 21 = Rwanda MCH

## Sex:
## 1 = male
## 2 = female

data.surv = as.data.frame(read_dta("Data/todd-survival-postsc.dta"))

data.surv$sero = decimal_date(as.Date(data.surv$seroconv_date, origin="1960-01-01"))
data.surv$exit = decimal_date(as.Date(data.surv$date,          origin="1960-01-01"))
data.surv$fpos = decimal_date(as.Date(data.surv$frstpos_date,  origin="1960-01-01"))
data.surv$agegr = cut(data.surv$infecage, breaks=c(15,25,35,45,100))
data.surv = subset(data.surv, site!=7) # Exclude South African miners

data.surv$int_age = floor(data.surv$infecage)
data.surv$wait = data.surv$exit - data.surv$sero
