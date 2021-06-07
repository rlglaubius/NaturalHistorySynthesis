library(readxl)


load.spec.defaults = function() {
  spec.param = list(
    mort = data.matrix(read_excel('Data/nathist-params-defaults.xlsx', range="B4:I10",  sheet="Estimates", col_names=sprintf("X%d", 1:8))),
    dist = data.matrix(read_excel('Data/nathist-params-defaults.xlsx', range="B25:I31", sheet="Estimates", col_names=sprintf("X%d", 1:8))),
    prog = data.matrix(read_excel('Data/nathist-params-defaults.xlsx', range="B15:I20", sheet="Estimates", col_names=sprintf("X%d", 1:8))))
  
  spec.param = lapply(spec.param, function(par) {
    stages = cd4.names[1:nrow(par)]
    arr = array(dim=c(nrow(par), 4, 2), dimnames=list(stages, age.names, sex.names))
    arr[,,1] = par[,1:4]
    arr[,,2] = par[,5:8]
    return(arr)
  })
  
  spec.param$dist = 0.01 * spec.param$dist # Convert from probabilities to proportions
  spec.param$art = data.matrix(read_excel('Data/nathist-params-defaults.xlsx', range="B34:B34", sheet="Estimates", col_names="X1"))
  
  return(spec.param)
}