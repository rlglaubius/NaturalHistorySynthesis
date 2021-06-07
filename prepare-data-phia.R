library(readxl)

## We load the PHIA data in long and wide formats. The wide format is used for
## inference, but drops confidence intervals on the data. Those CIs are retained
## in the long format so we can visualize them later.

data.phia.long = as.data.frame(read_excel('Data/phia-cd4-data.xlsx', sheet='Combined', col_names=TRUE))
data.phia.long$Percent = 0.01 * data.phia.long$Percent # convert from percentages to proportions
data.phia.long$LowerCL = 0.01 * data.phia.long$LowerCL
data.phia.long$UpperCL = 0.01 * data.phia.long$UpperCL

data.phia = reshape(data.phia.long,
                    timevar = "CD4",
                    idvar = c("Country", "Alpha", "Age", "Year"),
                    drop = c("N", "LowerCL", "UpperCL", "StdErr"),
                    direction = "wide")

## Add up sample sizes across CD4 categories
data.phia.size = plyr::ddply(data.phia.long, .(Country, Alpha, Age, Year), function(df) {
  data.frame(N = sum(df$N, na.rm=TRUE))
})
data.phia = dplyr::left_join(data.phia, data.phia.size, by=c("Country", "Alpha", "Age", "Year"))
data.phia[is.na(data.phia)] = 0.0 # fill in CD4 categories with no data (i.e., no respondents found)
colnames(data.phia) = gsub("Percent.", "", colnames(data.phia))

rm(data.phia.size)
