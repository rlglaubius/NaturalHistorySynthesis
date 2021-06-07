source("Utils/load-spec-defaults.R")

## Plot a comparison of all-cause mortality rate data (with 95% confidence
## intervals) to posterior distributions of mortality rates (with 95% credible
## intervals)
plot.post.mort = function(imis.rval, mort.data, tiff.name) {
  spec.param = load.spec.defaults()
  map.ind = which.max(imis.rval$prior + imis.rval$lhood)
  
  plot.data.mort = function(mort.raw, mort.adj, age) {
    ## Confidence interval calculation follows the code for poisson.test
    data.mu = mort.data$deaths / mort.data$pyears
    data.ci = rbind(qgamma(0.025, mort.data$deaths[,age]    ) / mort.data$pyears[,age],
                    qgamma(0.975, mort.data$deaths[,age] + 1) / mort.data$pyears[,age])
    
    raw.flat = sapply(mort.raw, function(param) {param[,age]})
    raw.cred = apply(raw.flat, 1, function(rate) {quantile(rate, c(0.025, 0.975))})
    raw.mode = mort.raw[[map.ind]]
    
    adj.flat = sapply(mort.adj, function(param) {param[,age]})
    adj.cred = apply(adj.flat, 1, function(rate) {quantile(rate, c(0.025, 0.975))})
    adj.mode = mort.adj[[map.ind]]
    
    spec.flat = spec.param$mort[,age,1]
    spec.flat = wpp2019mort[,age,1] + c(spec.flat[1:2], (2*spec.flat[3]+spec.flat[4])/3, spec.flat[5:7])

    xdata = 1:6 - 0.225
    xraw  = 1:6 - 0.075
    xadj  = 1:6 + 0.075
    xspec = 1:6 + 0.225
    
    col = c('#000000', brewer.pal(8,'Dark2'))
    par(las=1, mar=c(2.5,3.0,1.0,1.0), cex=1)
    plot(c(0.5,6.5), c(0,2), cex=0, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
    axis(1, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=0:6 + 0.5, labels=rep("",7))
    axis(1, mgp=c(0.1, 0.3, 0.0), tck=-0.00, at=1:6, labels=c(">500", "350-500", "200-349", "100-199", "50-99", "<50"), cex.axis=0.95)
    axis(2, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(2), labels=100*axTicks(2))
    points(xdata, data.mu[,age], pch=15, col=col[1], xpd=1)
    arrows(xdata, data.ci[1,], xdata, data.ci[2,], col=col[1], code=3, angle=90, length=0.01)
    points(xraw, raw.mode[,age], pch=16, col=col[3], xpd=1)
    arrows(xraw, raw.cred[1,], xraw, raw.cred[2,], col=col[3], code=3, angle=90, length=0.01)
    points(xadj, adj.mode[,age], pch=16, col=col[2], xpd=1)
    arrows(xadj, adj.cred[1,], xadj, adj.cred[2,], col=col[2], code=3, angle=90, length=0.01)
    points(xspec, spec.flat, col=col[5], pch=4, xpd=1)
    if (age == 1) {
      legend('topleft', inset=c(0.025,0.00), box.lty=0, bg=NA,
             legend=c('Observed rates', 'Model fit (unadjusted)', 'Model fit (adjusted)', 'Spectrum'),
             pch=c(15,16,16,4), lty=c(1,1,1,NA), col=col[c(1,3,2,5)], xpd=1)
    }
    title(ylab='Mortality rate, deaths/100 PY', line=1.75)
    title(xlab=expression(paste("CD4 cell count (cells/", mm^3, ")")), line=1.5)
    title(main=sprintf('Age %s', age.names[age]))
  }
  
  prop.f =  499/3497 ## 14% of CASCADE patients were females (dunn2008jid)
  prop.m = 2998/3497 ## 86% male

  ## Aggregate rates for 200-250 and 250-350 categories
  ## raw: Unadjusted for mortality underascertainment in dunn2008jid  
  mort.raw = mapply(function(param, r) {
    param$mort = param$mort / r # factor out the mortality rate ratio (mort.4)
    param.agg = array(NA, c(6,4,2))
    param.agg[,,1] = wpp2019mort[,,1] + rbind(param$mort[1:2,,1], (2*param$mort[3,,1] + param$mort[4,,1]) / 3, param$mort[5:7,,1])
    param.agg[,,2] = wpp2019mort[,,2] + rbind(param$mort[1:2,,2], (2*param$mort[3,,2] + param$mort[4,,2]) / 3, param$mort[5:7,,2])
    return(prop.m * param.agg[,,1] + prop.f * param.agg[,,2])
  }, imis.rval$param, data.frame(imis.rval$resample)$mort.4, SIMPLIFY=FALSE)

  ## adj: Adjusted for mortality underascertainment
  mort.adj = lapply(imis.rval$param, function(param) {
      param.agg = array(NA, c(6,4,2))
      param.agg[,,1] = wpp2019mort[,,1] + rbind(param$mort[1:2,,1], (2*param$mort[3,,1] + param$mort[4,,1]) / 3, param$mort[5:7,,1])
      param.agg[,,2] = wpp2019mort[,,2] + rbind(param$mort[1:2,,2], (2*param$mort[3,,2] + param$mort[4,,2]) / 3, param$mort[5:7,,2])
      return(prop.m * param.agg[,,1] + prop.f * param.agg[,,2])})
  
  tiff(tiff.name, units="mm", w=180, h=2*60, pointsize=8, compression="lzw", res=600)
  layout(matrix(1:4, nrow=2, byrow=TRUE))
  plot.data.mort(mort.raw, mort.adj, 1) # 15-24
  plot.data.mort(mort.raw, mort.adj, 2) # 25-34
  plot.data.mort(mort.raw, mort.adj, 3) # 35-44
  plot.data.mort(mort.raw, mort.adj, 4) # 45+
  dev.off()
}

## Plot a comparison of all-cause mortality rate data (with 95% confidence
## intervals) to posterior distributions of mortality rates (with posterior
## predictive distributions)
plot.pred.mort = function(imis.rval, mort.data, tiff.name, show.raw=FALSE, show.spec=FALSE) {
  src.names = c("500+", "350-499", "200-349", "100-199", "50-99", "<50") # CD4 category labels in data set
  cat.names = c(">500", "350-500", "200-349", "100-199", "50-99", "<50") # Standardized CD4 category names

  tiny.names = c("data", "raw", "adj", "spec")
  nice.names = c("Observed rates", "Model fit (unadjusted)", "Model fit (adjusted)", "Spectrum")
  
  map.ind = which.max(imis.rval$prior + imis.rval$lhood)
  
  ## Helper function. Takes Spectrum parameters in array format, aggregates
  ## mortality rates 200-250 and 250-350 CD4 categories, adds background
  ## mortality, then averages mortality rates by sex
  aggregate.mort = function(param) {
    prop.f =  499/3497 ## 14% of CASCADE patients were females (dunn2008jid)
    prop.m = 2998/3497 ## 86% male
    param.agg = array(NA, dim=c(6,4,2), dimnames=list(CD4=cat.names, Age=age.names, Sex=sex.names))
    param.agg[,,1] = wpp2019mort[,,1] + rbind(param$mort[1:2,,1], (2*param$mort[3,,1] + param$mort[4,,1]) / 3, param$mort[5:7,,1])
    param.agg[,,2] = wpp2019mort[,,2] + rbind(param$mort[1:2,,2], (2*param$mort[3,,2] + param$mort[4,,2]) / 3, param$mort[5:7,,2])
    return(prop.m * param.agg[,,1] + prop.f * param.agg[,,2])
  }

  ## Reformat the data for plotting
  mort.data.list = lapply(mort.data, as.data.frame.table)
  colnames(mort.data.list$deaths) = c("CD4", "Age", "Deaths")
  colnames(mort.data.list$pyears) = c("CD4", "Age", "PYears")
  mort.data.long = dplyr::left_join(mort.data.list$deaths, mort.data.list$pyears, by=c("CD4", "Age"))
  mort.data.long$CD4 = factor(mort.data.long$CD4, levels=src.names, labels=cat.names)
  mort.data.rate = mutate(mort.data.long,
                          Value = Deaths / PYears,
                          Lower = qgamma(0.025, Deaths    ) / PYears, # see poisson.test code
                          Upper = qgamma(0.975, Deaths + 1) / PYears)
  
  ## Reformat the unadjusted posterior parameter sets for plotting and prediction
  mort.post.raw.list = mapply(function(param, r) {
    param$mort = param$mort / r # factor out the mortality rate ratio (mort.4)
    param.avg = aggregate.mort(param)
    return(as.data.frame.table(param.avg, responseName="Value"))
  }, imis.rval$param, data.frame(imis.rval$resample)$mort.4, SIMPLIFY=FALSE)
  mort.post.raw.long = dplyr::bind_rows(mort.post.raw.list, .id="Sample")
  mort.post.raw.long$Sample = as.numeric(mort.post.raw.long$Sample)
  
  mort.post.raw.rate = plyr::ddply(mort.post.raw.long, .(CD4, Age), function(df) {
    ptEst = df$Value[which(df$Sample == map.ind)]
    bound = quantile(df$Value, c(0.025, 0.975))
    data.frame(Value = ptEst, Lower = bound[1], Upper = bound[2])
  })
  
  ## Reformat the adjusted posterior parameter sets for plotting and prediction
  mort.post.adj.list = lapply(imis.rval$param, function(param) {
    param.avg = aggregate.mort(param)
    return(as.data.frame.table(param.avg, responseName="Value"))
  })
  mort.post.adj.long = dplyr::bind_rows(mort.post.adj.list, .id="Sample")
  mort.post.adj.long$Sample = as.numeric(mort.post.adj.long$Sample)
  
  mort.post.adj.rate = plyr::ddply(mort.post.adj.long, .(CD4, Age), function(df) {
    ptEst = df$Value[which(df$Sample == map.ind)]
    bound = quantile(df$Value, c(0.025, 0.975))
    data.frame(Value = ptEst, Lower = bound[1], Upper = bound[2])
  })
  
  ## Reformat the Spectrum defaults for plotting
  spec.param = load.spec.defaults()
  spec.avg = aggregate.mort(spec.param)
  mort.spec.rate = as.data.frame.table(spec.avg, responseName="Value")
  mort.spec.rate$Lower = NA
  mort.spec.rate$Upper = NA

  mort.each.rate = dplyr::bind_rows(list(
    data = mort.data.rate,
    raw  = mort.post.raw.rate,
    adj  = mort.post.adj.rate,
    spec = mort.spec.rate
  ), .id="Source")
  mort.each.rate$Source = factor(mort.each.rate$Source, level=tiny.names, labels=nice.names)
  
  ## Sample the posterior predictive distribution of all-cause mortality rates
  mort.pred.raw.aug = dplyr::left_join(mort.post.raw.long, mort.data.long, by=c("CD4", "Age"))       # augment modeled rates with observed person-years
  mort.pred.raw.aug$Deaths = with(mort.pred.raw.aug, rpois(nrow(mort.pred.raw.aug), PYears * Value)) # predict numbers of deaths
  mort.pred.raw.aug$Value  = with(mort.pred.raw.aug, Deaths / PYears)                                # replace true rates with rates estimated from synthetic data
  
  mort.pred.adj.aug = dplyr::left_join(mort.post.adj.long, mort.data.long, by=c("CD4", "Age"))
  mort.pred.adj.aug$Deaths = with(mort.pred.adj.aug, rpois(nrow(mort.pred.adj.aug), PYears * Value))
  mort.pred.adj.aug$Value  = with(mort.pred.adj.aug, Deaths / PYears)
  
  ## Align violin plots by preparing dummy replicates of the data and Spectrum
  ## defaults. We just copy the original information once per posterior sample
  ## so that the resulting violin plots are empty.
  data.dummy = dplyr::bind_rows(lapply(1:nrow(imis.rval$resample), function(j) {mort.data.rate}), .id="Sample")
  data.dummy$Lower = NULL
  data.dummy$Upper = NULL
  data.dummy$Sample = as.numeric(data.dummy$Sample)
  
  spec.dummy = dplyr::bind_rows(lapply(1:nrow(imis.rval$resample), function(j) {mort.spec.rate}), .id="Sample")
  spec.dummy$Lower = NULL
  spec.dummy$Upper = NULL
  spec.dummy$Sample = as.numeric(spec.dummy$Sample)
  
  violin.data = dplyr::bind_rows(list(data = data.dummy,
                                      raw  = mort.pred.raw.aug,
                                      adj  = mort.pred.adj.aug,
                                      spec = spec.dummy), .id="Source")
  violin.data$Source = factor(violin.data$Source, level=tiny.names, labels=nice.names)
  
  show.mask = c(TRUE, show.raw, TRUE, show.spec)
  mort.each.rate = subset(mort.each.rate, Source %in% nice.names[show.mask])
  violin.data = subset(violin.data, Source %in% nice.names[show.mask])
  col = c("#000000", brewer.pal(8,"Dark2")[c(2,1,4)])[show.mask]
  pch = c(15, 16, 16, 4)[show.mask]

  ggplot(mort.each.rate, aes(x=CD4, y=100*Value, color=Source)) +
    geom_point(position=position_dodge(width=0.6), size=1, stroke=1, aes(color=Source, shape=Source)) +
    geom_pointrange(position=position_dodge(width=0.6), fatten=2, aes(ymin=100*Lower, ymax=100*Upper, color=Source, shape=Source), show.legend=FALSE) +
    geom_violin(data=violin.data, scale="width", position=position_dodge(0.6), color=NA, alpha=0.5, aes(x=CD4, y=100*Value, fill=Source), show.legend=FALSE) +
    facet_wrap(~Age, nrow=2) +
    scale_color_manual(values=col) +
    scale_fill_manual(values=col) +
    scale_shape_manual(values=pch) +
    xlab(expression(paste("CD4 cell count (cells/", mm^3, ")"))) +
    ylab("Mortality rate, deaths/100 PY") +
    ylim(0,200) +
    theme_bw() +
    theme(legend.position = "top",
          legend.margin = margin(t=0, b=0, unit="cm"),
          plot.margin = margin(t=0, b=0, l=0.05, r=0.05, unit="cm"),
          panel.border = element_rect(fill=NA, color="#000000"),
          strip.background = element_rect(fill="grey85", color="#000000"),
          text = element_text(size=10),
          axis.title = element_text(size=rel(1.0)),
          axis.text.x = element_text(color="#000000", size=rel(1.0)),
          axis.text.y = element_text(color="#000000", size=rel(1.0)))
  ggsave(tiff.name, compression="lzw", dpi=600, units="mm", width=180, height=2*60)
}
