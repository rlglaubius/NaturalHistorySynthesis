
plot.trend.dist = function(imis.out, lwd, col) {
  draws = as.data.frame(imis.out$resample)
  map.ind = which.max(imis.out$prior + imis.out$lhood)

  eval.x = seq(0, 1200, 10)
  eval.y = lapply(0:3, function(a) {
    med = draws$dist.2 * (1.0 - a * draws$dist.3)
    shp = draws$dist.1 + 1
    sapply(eval.x, function(x) {dfisk(x, scale=med, shape=shp)})
  })
  cred.y = lapply(eval.y, function(vals) {apply(vals, 2, function(y) {quantile(y, c(0.025, 0.975))})})

  xlim = range(eval.x)
  ylim = c(0, pretty.upper(sapply(cred.y, max)))

  par(las=1, mar=c(2.5,3.5,1.0,1.0), cex=1)
  plot(xlim, ylim, cex=0, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
  for (a in 1:4) {polygon(c(eval.x, rev(eval.x)), c(cred.y[[a]][1,], rev(cred.y[[a]][2,])), col=sprintf("%s30", col[a]), border=NA)}
  for (a in 1:4) {lines(eval.x, eval.y[[a]][map.ind,], col=col[a], lwd=lwd)}
  axis(1, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(1), labels=axTicks(1))
  axis(2, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=ylim, labels=rep("", 2))
  legend('topright', inset=c(0.025,0.00), box.lty=0, bg=NA, legend=age.names, lty=1, lwd=lwd, col=col, xpd=1)
  title(xlab=expression(paste("CD4 cell count (cells/", mm^3, ")")), line=1.50)
  title(ylab="Density", line=1.75)
}


plot.trend.mort = function(imis.out, lwd, col) {
  draws = as.data.frame(imis.out$resample)
  map.ind = which.max(imis.out$prior + imis.out$lhood)

  eval.x = seq(0, 1200, 10)
  eval.y = lapply(0:3, function(a) {
    base = draws$mort.2 * (1.0 + a * draws$mort.3)
    sapply(eval.x, function(x) {base * draws$mort.1^x})
  })
  cred.y = lapply(eval.y, function(vals) {apply(vals, 2, function(y) {quantile(y, c(0.025, 0.975))})})

  xlim = c(0, 600)
  ylim = c(0, 1)

  par(las=1, mar=c(2.5,3.5,1.0,1.0), cex=1)
  plot(xlim, ylim, cex=0, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
  for (a in 1:4) {polygon(c(eval.x, rev(eval.x)), c(cred.y[[a]][1,], rev(cred.y[[a]][2,])), col=sprintf("%s30", col[a]), border=NA)}
  for (a in 1:4) {lines(eval.x, eval.y[[a]][map.ind,], col=col[a], lwd=lwd)}
  axis(1, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(1), labels=axTicks(1))
  axis(2, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(2), labels=100*axTicks(2))
  legend('topright', inset=c(0.025,0.00), box.lty=0, bg=NA, legend=age.names, lty=1, lwd=lwd, col=col, xpd=1)
  title(xlab=expression(paste("CD4 cell count (cells/", mm^3, ")")), line=1.50)
  title(ylab="Mortality rate (deaths/100 PY)", line=1.75)
}


plot.trend.prog = function(imis.out, lwd, col) {
  conditional.median = function(x, scale, shape) {
    return(qfisk(1.0 - 0.5 * (1.0 - pfisk(x, scale, shape)), scale, shape))
  }

  draws = as.data.frame(imis.out$resample)
  map.ind = which.max(imis.out$prior + imis.out$lhood)

  eval.x = seq(0,20,0.1)
  eval.y = lapply(0:3, function(a) {
    dist.shape = draws$dist.1 + 1
    dist.scale = draws$dist.2 * (1.0 - a * draws$dist.3)
    c0 = conditional.median(500, dist.scale, dist.shape)
    prog.scale = with(draws, prog.2 * (1.0 - a * prog.3))
    prog.shape = draws$prog.1
    sapply(eval.x, function(x) {ifelse(x <= prog.scale, c0 * (1.0 - x / prog.scale)^(1.0 + prog.shape), 0.0)})
  })
  cred.y = lapply(eval.y, function(vals) {apply(vals, 2, function(y) {quantile(y, c(0.025, 0.975), na.rm=TRUE)})})

  xlim = range(eval.x)
  ylim = c(0, pretty.upper(sapply(cred.y, max)))

  par(las=1, mar=c(2.5,3.5,1.0,1.0), cex=1)
  plot(xlim, ylim, cex=0, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
  for (a in 1:4) {polygon(c(eval.x, rev(eval.x)), c(cred.y[[a]][1,], rev(cred.y[[a]][2,])), col=sprintf("%s30", col[a]), border=NA)}
  for (a in 1:4) {lines(eval.x, eval.y[[a]][map.ind,], col=col[a], lwd=lwd)}
  axis(1, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(1), labels=axTicks(1))
  axis(2, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(2), labels=axTicks(2))
  legend('topright', inset=c(0.025,0.00), box.lty=0, bg=NA, legend=age.names, lty=1, lwd=lwd, col=col, xpd=1)
  title(xlab="Years since seroconversion", line=1.50)
  title(ylab=expression(paste("CD4 cell count (cells/", mm^3, ")")), line=1.75)
}
