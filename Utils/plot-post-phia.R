## Use EPP-ASM to calculate CD4 distributions in untreated PLHIV at each PHIA
## for each posterior sample
calculate.phia.sims = function(imis.rval, pjnz.info, phia.data) {
  phia.post = lapply(names(pjnz.info), function(alpha.code) {
    phia.year = phia.data$Year[which(data.phia$Alpha == alpha.code)][1]
    lapply(imis.rval$param, function(param) {
      cd4.dist = eval.cd4.dist(pjnz.info[[alpha.code]], param, phia.year)
      dimnames(cd4.dist) = list(CD4=rownames(cd4.dist), Age=colnames(cd4.dist))
      return(cd4.dist)
    })
  })
  names(phia.post) = names(pjnz.info)
  return(phia.post)
}

plot.post.phia = function(imis.rval, post.sims, phia.data, tiff.name) {
  src.names = c("CD4 >=500", "CD4 350-499", "CD4 200-349", "CD4 100-199", "CD4 50-99", "CD4 <50")
  cat.names = c(">500", "350-500", "200-349", "100-199", "50-99", "<50")
  
  tiny.names = c("data", "post")
  nice.names = c("Data", "Model fit")
  
  ## Reformat posterior CD4 counts for plotting
  map.ind = which.max(imis.rval$prior + imis.rval$lhood)
  post.sims.list = lapply(post.sims, function(country) {
    post.long = dplyr::bind_rows(lapply(country, as.data.frame.table), .id="Sample")
    post.long$Sample = as.numeric(post.long$Sample)
    return(post.long)
  })
  post.sims.long = dplyr::bind_rows(post.sims.list, .id="Alpha")
  post.sims.summ = plyr::ddply(post.sims.long, .(Alpha, CD4, Age), function(df) {
    ptEst = df$Freq[which(df$Sample == map.ind)]
    bound = quantile(df$Freq, c(0.025, 0.975))
    data.frame(Percent = ptEst, LowerCL = bound[1], UpperCL = bound[2])
  })
  post.sims.summ$CD4 = factor(post.sims.summ$CD4, levels=cat.names)
  post.sset = subset(post.sims.summ, Age=="15-49")

  phia.cols = c("Alpha", "Age", "CD4", "Percent", "LowerCL", "UpperCL")
  phia.data$CD4 = factor(phia.data$CD4, levels=src.names, labels=cat.names)
  phia.sset = subset(phia.data, Age=="15-49" & Alpha %in% metadata$alpha)[,phia.cols]
  phia.sset$Percent[is.na(phia.sset$Percent)] = 0.0
  
  aggr.long = dplyr::bind_rows(list(data = phia.sset, post = post.sset), .id="Source")
  aggr.long$Country = plyr::mapvalues(aggr.long$Alpha, from=metadata$alpha, to=metadata$country)
  aggr.long$Country[which(aggr.long$Alpha=="SWZ")] = "Eswatini"
  aggr.long$Country[which(aggr.long$Alpha=="CIV")] = "Côte d'Ivoire"
  aggr.long$Source = factor(aggr.long$Source, levels=tiny.names, labels=nice.names)
  
  ggplot(aggr.long, aes(x=CD4, y=100*Percent, ymin=100*LowerCL, ymax=100*UpperCL, group=Source)) +
    geom_pointrange(aes(shape=Source, color=Source), size=0.25) +
    facet_wrap(~Country, nrow=2) +
    xlab("CD4 cell count") + ylab("Share of PLHIV off ART, %") +
    scale_color_manual(values=c("#000000", brewer.pal(8,"Dark2")[1])) +
    scale_shape_manual(values=c(15,16)) +
    theme_bw() +
    theme(legend.position = "top",
          legend.margin = margin(t=0, b=0, unit="cm"),
          plot.margin = margin(t=0, b=0, l=0.05, r=0.05, unit="cm"),
          panel.border = element_rect(fill=NA, color="#000000"),
          strip.background = element_rect(fill="grey85", color="#000000"),
          text = element_text(size=10),
          axis.title = element_text(size=rel(1.0)),
          axis.text.x = element_text(color="#000000", size=rel(0.8), angle=45, hjust=1),
          axis.text.y = element_text(color="#000000", size=rel(0.8)))
  ggsave(tiff.name, compression="lzw", dpi=600, units="mm", width=180, height=1.5*60)
}

plot.pred.phia = function(imis.rval, post.sims, phia.data, tiff.name) {
  src.names = c("CD4 >=500", "CD4 350-499", "CD4 200-349", "CD4 100-199", "CD4 50-99", "CD4 <50")
  cat.names = c(">500", "350-500", "200-349", "100-199", "50-99", "<50")
  
  tiny.names = c("data", "post")
  nice.names = c("Data", "Model fit")
  
  ## Reformat posterior CD4 counts for plotting
  map.ind = which.max(imis.rval$prior + imis.rval$lhood)
  post.sims.list = lapply(post.sims, function(country) {
    post.long = dplyr::bind_rows(lapply(country, as.data.frame.table), .id="Sample")
    post.long$Sample = as.numeric(post.long$Sample)
    return(post.long)
  })
  post.sims.long = dplyr::bind_rows(post.sims.list, .id="Alpha")
  post.sims.summ = plyr::ddply(post.sims.long, .(Alpha, CD4, Age), function(df) {
    ptEst = df$Freq[which(df$Sample == map.ind)]
    bound = quantile(df$Freq, c(0.025, 0.975))
    data.frame(Percent = ptEst, LowerCL = bound[1], UpperCL = bound[2])
  })
  post.sims.summ$CD4 = factor(post.sims.summ$CD4, levels=cat.names)
  post.sset = subset(post.sims.summ, Age=="15-49")
  
  phia.cols = c("Alpha", "Age", "CD4", "Percent", "LowerCL", "UpperCL")
  phia.data$CD4 = factor(phia.data$CD4, levels=src.names, labels=cat.names)
  phia.sset = subset(phia.data, Age=="15-49" & Alpha %in% metadata$alpha)[,phia.cols]
  phia.sset$Percent[is.na(phia.sset$Percent)] = 0.0
  
  aggr.long = dplyr::bind_rows(list(data = phia.sset, post = post.sset), .id="Source")
  aggr.long$Country = plyr::mapvalues(aggr.long$Alpha, from=metadata$alpha, to=metadata$country)
  aggr.long$Country[which(aggr.long$Alpha=="SWZ")] = "Eswatini"
  aggr.long$Country[which(aggr.long$Alpha=="CIV")] = "Côte d'Ivoire"
  aggr.long$Source = factor(aggr.long$Source, levels=tiny.names, labels=nice.names)
  
  ## get PHIA sample sizes
  phia.size = plyr::ddply(phia.data, .(Alpha, Age), function(df) {data.frame(N = sum(df$N, na.rm=TRUE))})
  
  ## Prepare synthetic data sets. We drop 15-49 because we will reconstruct this
  ## age group by summing over constituent five year age groups
  post.sims.wide = reshape2::dcast(subset(post.sims.long, Age != "15-49"), Alpha + Sample + Age ~ CD4, value.var="Freq")
  post.sims.augm = dplyr::left_join(post.sims.wide, phia.size, by=c("Alpha", "Age")) # augment with numbers of PHIA respondents by country and age group
  
  pred.raw = rmultinomial(nrow(post.sims.augm), size=post.sims.augm$N, prob=post.sims.augm[,cat.names])
  colnames(pred.raw) = cat.names
  pred.wide = cbind(post.sims.augm[,c("Alpha", "Sample", "Age", "N")], pred.raw)
  pred.flat = reshape2::melt(pred.wide, id.vars=c("Alpha", "Sample", "Age"), measure.vars=cat.names, variable.name="CD4", value.name="Count")
  
  ## Aggregate up to 15-49
  pred.aggr = plyr::ddply(pred.flat, .(Alpha, CD4, Sample), function(df) {
    data.frame(Age = "15-49", Count = sum(df$Count))
  })
  
  ## Calculate proportions
  pred.prop = plyr::ddply(pred.aggr, .(Alpha, Sample), function(df) {
    data.frame(CD4 = df$CD4, Percent = df$Count / sum(df$Count), Source="post")
  })
  
  ## To align violin plots of the posterior predictive distribution with point
  ## estimates in the figure, we need to have a dummy dataset for "data" that has
  ## the same size as the pred.prop
  data.dummy = dplyr::bind_rows(lapply(1:nrow(imis.rval$resample), function(j) {phia.sset}), .id="Sample")
  data.dummy$Lower = NULL
  data.dummy$Upper = NULL
  data.dummy$Source = "data"
  data.dummy$Sample = as.numeric(data.dummy$Sample)
  
  violin.data = dplyr::bind_rows(list(data=data.dummy, post=pred.prop))
  violin.data$Source = factor(violin.data$Source, level=tiny.names, labels=nice.names)
  violin.data$Country = plyr::mapvalues(violin.data$Alpha, from=metadata$alpha, to=metadata$country)
  violin.data$Country[which(violin.data$Alpha=="SWZ")] = "Eswatini"
  violin.data$Country[which(violin.data$Alpha=="CIV")] = "Côte d'Ivoire"
  
  ggplot(aggr.long, aes(x=CD4, y=100*Percent, color=Source)) +
    geom_pointrange(position=position_dodge(width=0.6), fatten=2, aes(ymin=100*LowerCL, ymax=100*UpperCL, shape=Source, color=Source)) +
    geom_violin(data=violin.data, scale="width", position=position_dodge(0.6), color=NA, alpha=0.5, aes(x=CD4, y=100*Percent, fill=Source), show.legend=FALSE) +
    facet_wrap(~Country, nrow=2) +
    xlab("CD4 cell count") + ylab("Share of PLHIV off ART, %") +
    scale_color_manual(values=c("#000000", brewer.pal(8,"Dark2")[1])) +
    scale_fill_manual(values=c("#000000", brewer.pal(8,"Dark2")[1])) +
    scale_shape_manual(values=c(15,16)) +
    theme_bw() +
    theme(legend.position = "top",
          legend.margin = margin(t=0, b=0, unit="cm"),
          plot.margin = margin(t=0, b=0, l=0.05, r=0.05, unit="cm"),
          panel.border = element_rect(fill=NA, color="#000000"),
          strip.background = element_rect(fill="grey85", color="#000000"),
          text = element_text(size=10),
          axis.title = element_text(size=rel(1.0)),
          axis.text.x = element_text(color="#000000", size=rel(0.8), angle=45, hjust=1),
          axis.text.y = element_text(color="#000000", size=rel(0.8)))
  ggsave(tiff.name, compression="lzw", dpi=600, units="mm", width=180, height=1.5*60)
  
  ## Check how many posterior predictive intervals contain the data
  pred.summ = plyr::ddply(pred.prop, .(Alpha, CD4), function(df) {
    q = quantile(df$Percent, c(0.0, 0.025, 0.5, 0.975, 1.0))
    names(q) = c("Min", "Q02.5", "Median", "Q97.5", "Max")
    return(q)
  })
  pred.summ.augm = dplyr::left_join(pred.summ, phia.sset, by=c("Alpha", "CD4"))
  
  qr100.numer = with(pred.summ.augm, sum(Percent >= Min   & Percent <= Max))
  qr095.numer = with(pred.summ.augm, sum(Percent >= Q02.5 & Percent <= Q97.5))
  
  cat(sprintf("Posterior predictive interval coverage\n"))
  cat(sprintf("\t100%% interval: %0.0f%% (%d/%d)\n", 100 * qr100.numer / nrow(pred.summ.augm), qr100.numer, nrow(pred.summ.augm)))
  cat(sprintf("\t 95%% interval: %0.0f%% (%d/%d)\n", 100 * qr095.numer / nrow(pred.summ.augm), qr095.numer, nrow(pred.summ.augm)))
}


