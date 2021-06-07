pretty.upper = function(values, base=10, na.rm=FALSE) {
  m = max(values, na.rm=na.rm)
  a = floor(log(m, base))
  return(base^a * ceiling(m / base^a))
}

pretty.scale = function(values, base=1000, na.rm=FALSE) {
  upper = pretty.upper(values, na.rm=na.rm)
  m = floor(log(upper, base))
  return(base^m)
}
