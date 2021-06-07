dfisk = function(x, scale, shape, log=FALSE) {
  lnfx = log(shape) + (shape - 1.0) * log(x) - shape * log(scale) - 2.0 * log1p((x / scale)^shape)
  return(if(log == TRUE) lnfx else exp(lnfx))
}

pfisk = function(q, scale, shape, log.p=FALSE) {
  lnFx = -log(1.0 + (q / scale)^(-shape))
  return(if(log.p == TRUE) lnFx else exp(lnFx))
}

qfisk = function(p, scale, shape, log.p=FALSE) {
  lnQx = log(scale) + (-1.0 / shape) * log(1.0 / p - 1.0)
  return(if(log.p == TRUE) lnQx else exp(lnQx))
}

rfisk = function(n, scale, shape) {
  return(qfisk(runif(n), scale, shape))
}
