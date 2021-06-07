library(openxlsx)

write.params.xlsx = function(model, xlsx.file) {
  rnames = c("CD4>500", "350-500", "250-349", "200-249", "100-199", "50-99", "CD4<50")
  cnames = c('Male 15-24', 'Male 25-34', 'Male 35-44', 'Male 45+', 'Female 15-24', 'Female 25-34', 'Female 35-44', 'Female 45+')
  titles = list(mort = 'Annual rate of HIV-related mortality when not on ART',
                prog = 'Annual rate of progression to next lower CD4 category',
                dist = 'Distribution of new infections by CD4 count (percent)',
                init = 'Allocation method for new ART patients')
  
  mort = cbind(model$mort[,,1], model$mort[,,2])
  colnames(mort) = cnames
  rownames(mort) = rnames
  
  prog = cbind(model$prog[,,1], model$prog[,,2])
  colnames(prog) = cnames
  rownames(prog) = rnames[1:6]
  
  ## We lose some precision when going from Excel to Spectrum, so the resulting
  ## distribution may not quite add up to 100%. To avoid this, we round
  ## percentages to 2 decimal places then put the residual in CD4>500
  dist = round(100 * cbind(model$dist[,,1], model$dist[,,2]), 2)
  dist[1,] = dist[1,] + 100 - colSums(dist)
  colnames(dist) = cnames
  rownames(dist) = rnames
  
  init = rbind(model$art, 1.0 - model$art)
  rownames(init) = c('Expected mortality', 'Proportion eligible')
  
  wb = createWorkbook()
  addWorksheet(wb, "Estimates")
  
  ## Mortality block
  addStyle(wb, 1, createStyle(numFmt="0.000"), cols=2:9, rows=4:10, gridExpand=TRUE)
  writeData(wb, 1, titles$mort, startCol=1, startRow=1, colNames=FALSE, rowNames=FALSE)
  writeData(wb, 1, mort, startCol=1, startRow=3, colNames=TRUE, rowNames=TRUE)
  
  ## Progression block
  addStyle(wb, 1, createStyle(numFmt="0.000"), cols=2:9, rows=15:20, gridExpand=TRUE)
  writeData(wb, 1, titles$prog, startCol=1, startRow=12, colNames=FALSE, rowNames=FALSE)
  writeData(wb, 1, prog, startCol=1, startRow=14, colNames=TRUE, rowNames=TRUE)
  
  ## Initial CD4 count block
  addStyle(wb, 1, createStyle(numFmt="0.00"), cols=2:9, rows=24:31, gridExpand=TRUE)
  writeData(wb, 1, titles$dist, startCol=1, startRow=22, colNames=FALSE, rowNames=FALSE)
  writeData(wb, 1, dist, startCol=1, startRow=24, colNames=TRUE, rowNames=TRUE)

  ## ART initiation weights
  addStyle(wb, 1, createStyle(numFmt="0.000"), cols=2, rows=35:36, gridExpand=TRUE)
  writeData(wb, 1, titles$init, startCol=1, startRow=33, colNames=FALSE, rowNames=FALSE)
  writeData(wb, 1, init, startCol=1, startRow=34, colNames=TRUE, rowNames=TRUE)
  
  saveWorkbook(wb, xlsx.file, overwrite=TRUE)
}

gen.param.table = function(param.mode, param.post, fmt="%0.2f", col.names, row.names) {
  lower = apply(param.post, 1:2, function(x) {quantile(x, c(0.025))})
  upper = apply(param.post, 1:2, function(x) {quantile(x, c(0.975))})
  if (any(lower > param.mode) | any(upper < param.mode)) {warning("Bounds do not include all point estimates")}
  fmtstr = sprintf("%s (%s-%s)", fmt, fmt, fmt)
  rval = matrix(sprintf("%0.1f (%0.1f-%0.1f)", param.mode, lower, upper), ncol=ncol(param.mode))
  colnames(rval) = col.names
  rownames(rval) = row.names
  return(rval)
}

sex.age.names = apply(expand.grid(age.names, sex.names), 1, function(val) {sprintf("%s %s", val[2], val[1])})

write.params.xlsx(imis.out$param[[post.mode.ind]], "param-spectrum-format.xlsx")

param.table = list(
  dist = array(sapply(imis.out$param, function(param) {array(param$dist, c(7,8))}), c(7,8,3000), dimnames=list(cd4.names, sex.age.names, 1:3000)),
  mort = array(sapply(imis.out$param, function(param) {array(param$mort, c(7,8))}), c(7,8,3000), dimnames=list(cd4.names, sex.age.names, 1:3000)),
  prog = array(sapply(imis.out$param, function(param) {array(param$prog, c(6,8))}), c(6,8,3000), dimnames=list(cd4.names[1:6], sex.age.names, 1:3000)))
write.csv(gen.param.table(100 * array(imis.out$param[[post.mode.ind]]$dist, c(7,8)), 100*param.table$dist, "%0.1f", sex.age.names, cd4.names), "param-post-dist.csv")
write.csv(gen.param.table(100 * array(imis.out$param[[post.mode.ind]]$prog, c(6,8)), 100*param.table$prog, "%0.1f", sex.age.names, cd4.names[1:6]), "param-post-prog.csv")
write.csv(gen.param.table(100 * array(imis.out$param[[post.mode.ind]]$mort, c(7,8)), 100*param.table$mort, "%0.1f", sex.age.names, cd4.names), "param-post-mort.csv")
write.csv(rbind(
  imis.out$resample[post.mode.ind,],
  apply(imis.out$resample, 2, function(hyper) {quantile(hyper, c(0.025, 0.975))})),
  "param-post-hyper.csv")

