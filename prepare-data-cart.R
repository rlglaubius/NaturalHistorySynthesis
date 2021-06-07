library(readxl)

data.cart = as.data.frame(read_excel('Data/iedea_cd4_categories_v2.xlsx', na="NA"))
data.cart = data.cart[,c('country', 'recart_yr', 'gender', 'n_lower_50', 'n_50_to_100', 'n_100_to_200', 'n_200_to_350', 'n_350_to_500', 'n_greater_500')]
colnames(data.cart) = c('country', 'year', 'gender', 'CD4 < 50', 'CD4 50-99', 'CD4 100-199', 'CD4 200-349', 'CD4 350-500', 'CD4 >= 500')
