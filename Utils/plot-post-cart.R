## Use EPP-ASM to calculate CD4 distributions at ART initiation
## for each posterior sample
calculate.cart.sims = function(imis.rval, pjnz.info, cart.data) {
  countries = intersect(names(pjnz.info), cart.data$country)
  cart.post = lapply(countries, function(alpha.code) {
    lapply(imis.rval$param, function(param) {
      cd4.dist = eval.art.dist(pjnz.info[[alpha.code]], param, 2000:2020)
      cd4.dist.arr = array(NA, c(6, 21, 2), dimnames=list(CD4=rownames(cd4.dist$m), Year=colnames(cd4.dist$m), Sex=sex.names))
      cd4.dist.arr[,,1] = cd4.dist$m
      cd4.dist.arr[,,2] = cd4.dist$f
      return(cd4.dist.arr)
    })
  })
  names(cart.post) = countries
  return(cart.post)
}

plot.post.cart = function(imis.rval, post.sims, cart.data, tiff.name) {
  countries = intersect(metadata$alpha, cart.data$country)

  src.names = c("CD4 >= 500", "CD4 350-500", "CD4 200-349", "CD4 100-199", "CD4 50-99", "CD4 < 50")
  cat.names = c(">500", "350-500", "200-349", "100-199", "50-99", "<50")
  cat.group = c("CD4>350", "CD4 200-350", "CD4 50-200", "CD4<50")

  tiny.names = c("data", "post")
  nice.names = c("Data", "Model fit")

  map.ind = which.max(imis.rval$prior + imis.rval$lhood)

  ## Reformat data for plotting
  cart.wide = subset(cart.data, country %in% countries)
  cart.long = reshape2::melt(cart.wide,
                             id.vars = c("country", "year", "gender"),
                             measure.vars = src.names,
                             variable.name = "CD4",
                             value.name = "Freq")
  cart.long = dplyr::rename(cart.long, Alpha = country, Sex = gender, Year = year)
  cart.long$Stage = factor(plyr::mapvalues(cart.long$CD4, from=src.names, to=rep(cat.group, c(2, 1, 2, 1))), levels=cat.group)
  cart.aggr = plyr::ddply(cart.long, .(Year, Stage), function(df) {
    data.frame(Freq = sum(df$Freq)) # aggregate frequencies by compact CD4 category, sex, and country
  })
  cart.prop = plyr::ddply(cart.aggr, .(Year), function(df) {
    data.frame(Stage = df$Stage, Share = df$Freq / sum(df$Freq)) # convert from frequencies to proportions
  })

  ## Calculate weights by country and sex needed for aggregating simulation
  ## outputs. These weights are the number of people observed by year, sex, and
  ## country, summed across CD4 categories at ART initiation
  cart.weight = plyr::ddply(cart.long, .(Alpha, Sex, Year), function(df) {
    data.frame(Weight = sum(df$Freq))
  })

  ## Reformat posterior simulations for plotting
  post.list = lapply(post.sims, function(country) {
    sims.list = lapply(country, function(sim) {
      dn = dimnames(sim)
      sim.aggr = array(NA, dim=c(4, 21, 2), dimnames=list(Stage=cat.group, Year=dn$Year, Sex=dn$Sex))
      sim.aggr[,,1] = rowsum(sim[,,1], c(1,1,2,3,3,4))
      sim.aggr[,,2] = rowsum(sim[,,2], c(1,1,2,3,3,4))
      return(as.data.frame.table(sim.aggr, responseName="Share"))
    })
    sims.flat = dplyr::bind_rows(sims.list, .id="Sample")
    sims.flat$Sample = as.numeric(sims.flat$Sample)
    return(sims.flat)
  })
  post.long = dplyr::bind_rows(post.list, .id="Alpha")
  post.long$Year = as.numeric(as.character(post.long$Year)) # convert Year from factor to number
  post.augm = dplyr::left_join(post.long, cart.weight, by=c("Alpha", "Sex", "Year")) # Augment with data weights
  post.augm$Weight[is.na(post.augm$Weight)] = 0.0 # give zero weight when there are no data
  post.aggr = plyr::ddply(post.augm, .(Year, Stage, Sample), function(df) {
    data.frame(Freq = sum(df$Weight * df$Share)) # aggregate across sex and country, weighted by numbers in data set
  })
  post.prop = plyr::ddply(post.aggr, .(Year, Sample), function(df) {
    data.frame(Stage = df$Stage, Share = df$Freq / sum(df$Freq)) # calculate distributions across CD4 categories
  })
  post.summ = plyr::ddply(post.prop, .(Year, Stage), function(df) {
    ptEst = df$Share[which(df$Sample == map.ind)]
    bound = quantile(df$Share, c(0.025, 0.975), na.rm=TRUE)
    data.frame(Share=ptEst, Lower=bound[1], Upper=bound[2], row.names=NULL)
  })

  both.prop = dplyr::bind_rows(list(data = cart.prop, post = post.summ), .id="Source")
  both.prop$Source = factor(both.prop$Source, levels=tiny.names, labels=nice.names)
  both.prop$Stage = factor(both.prop$Stage, levels=rev(cat.group))

  ggplot(subset(both.prop, Year>=2005 & Year<=2014), aes(x=Year, y=100*Share, ymin=100*Lower, ymax=100*Upper, color=Source, fill=Source)) +
    geom_line() +
    geom_ribbon(alpha=0.2, color=NA) +
    facet_wrap(~Stage, nrow=1) +
    scale_color_manual(values=c("#000000", brewer.pal(8,"Dark2")[1:2])) +
    scale_fill_manual(values=c("#000000", brewer.pal(8,"Dark2")[1:2])) +
    xlab("Year") + ylab("Share of ART initiators, %") +
    xlim(2004,2015) +
    theme_bw() +
    theme(legend.position = "right",
          legend.margin = margin(t=0, b=0, unit="cm"),
          plot.margin = margin(t=0, b=0, l=0.05, r=0.05, unit="cm"),
          panel.border = element_rect(fill=NA, color="#000000"),
          strip.background = element_rect(fill="grey85", color="#000000"),
          text = element_text(size=10),
          axis.title = element_text(size=rel(1.0)),
          axis.text.x = element_text(color="#000000", size=rel(0.8)),
          axis.text.y = element_text(color="#000000", size=rel(0.8)))
  ggsave(tiff.name, compression="lzw", dpi=600, units="mm", width=180, height=60)
  
  cat(sprintf("\nObserved proportion initiating at CD4<50\n"))
  print(subset(both.prop, Year %in% c(2005, 2014) & Source=="Data"  & Stage=="CD4<50"))
  cat(sprintf("\nModeled proportion initiating at CD4<50\n"))
  print(subset(both.prop, Year %in% c(2005, 2014) & Source=="Model fit" & Stage=="CD4<50"))
  cat(sprintf("\nObserved proportion initiating at CD4>350 during 2005-2012\n"))
  print(range(subset(both.prop, Year >= 2005 & Year <= 2012 & Source=="Data"  & Stage=="CD4>350")$Share))
  cat(sprintf("Modeled proportion initiating at CD4>350 during 2005-2012\n"))
  print(range(subset(both.prop, Year >= 2005 & Year <= 2012 & Source=="Model fit" & Stage=="CD4>350")$Share))
}

plot.pred.cart = function(imis.rval, post.sims, cart.data, tiff.name) {
  countries = intersect(metadata$alpha, cart.data$country)

  src.names = c("CD4 >= 500", "CD4 350-500", "CD4 200-349", "CD4 100-199", "CD4 50-99", "CD4 < 50")
  cat.names = c(">500", "350-500", "200-349", "100-199", "50-99", "<50")
  cat.group = c("CD4>350", "CD4 200-350", "CD4 50-200", "CD4<50")

  tiny.names = c("data", "post")
  nice.names = c("Data", "Model fit")

  map.ind = which.max(imis.rval$prior + imis.rval$lhood)

  ## Reformat data for plotting
  cart.wide = subset(cart.data, country %in% countries)
  cart.long = reshape2::melt(cart.wide,
                             id.vars = c("country", "year", "gender"),
                             measure.vars = src.names,
                             variable.name = "CD4",
                             value.name = "Freq")
  cart.long = dplyr::rename(cart.long, Alpha = country, Sex = gender, Year = year)
  cart.long$Stage = factor(plyr::mapvalues(cart.long$CD4, from=src.names, to=rep(cat.group, c(2, 1, 2, 1))), levels=cat.group)
  cart.aggr = plyr::ddply(cart.long, .(Year, Stage), function(df) {
    data.frame(Freq = sum(df$Freq)) # aggregate frequencies by compact CD4 category, sex, and country
  })
  cart.prop = plyr::ddply(cart.aggr, .(Year), function(df) {
    data.frame(Stage = df$Stage, Share = df$Freq / sum(df$Freq)) # convert from frequencies to proportions
  })

  ## Calculate weights by country and sex needed for aggregating simulation
  ## outputs. These weights are the number of people observed by year, sex, and
  ## country, summed across CD4 categories at ART initiation
  cart.weight = plyr::ddply(cart.long, .(Alpha, Sex, Year), function(df) {
    data.frame(Weight = sum(df$Freq))
  })

  ## Reformat posterior simulations for plotting
  post.list = lapply(post.sims, function(country) {
    sims.list = lapply(country, function(sim) {
      dn = dimnames(sim)
      sim.aggr = array(NA, dim=c(4, 21, 2), dimnames=list(Stage=cat.group, Year=dn$Year, Sex=dn$Sex))
      sim.aggr[,,1] = rowsum(sim[,,1], c(1,1,2,3,3,4))
      sim.aggr[,,2] = rowsum(sim[,,2], c(1,1,2,3,3,4))
      return(as.data.frame.table(sim.aggr, responseName="Share"))
    })
    sims.flat = dplyr::bind_rows(sims.list, .id="Sample")
    sims.flat$Sample = as.numeric(sims.flat$Sample)
    return(sims.flat)
  })
  post.long = dplyr::bind_rows(post.list, .id="Alpha")
  post.long$Year = as.numeric(as.character(post.long$Year)) # convert Year from factor to number
  post.augm = dplyr::left_join(post.long, cart.weight, by=c("Alpha", "Sex", "Year")) # Augment with data weights
  post.augm$Weight[is.na(post.augm$Weight)] = 0.0 # give zero weight when there are no data
  post.aggr = plyr::ddply(post.augm, .(Year, Stage, Sample), function(df) {
    data.frame(Freq = sum(df$Weight * df$Share)) # aggregate across sex and country, weighted by numbers in data set
  })
  post.prop = plyr::ddply(post.aggr, .(Year, Sample), function(df) {
    data.frame(Stage = df$Stage, Share = df$Freq / sum(df$Freq)) # calculate distributions across CD4 categories
  })
  post.summ = plyr::ddply(post.prop, .(Year, Stage), function(df) {
    ptEst = df$Share[which(df$Sample == map.ind)]
    bound = quantile(df$Share, c(0.025, 0.975), na.rm=TRUE)
    data.frame(Share=ptEst, Lower=bound[1], Upper=bound[2], row.names=NULL)
  })

  ## Generate synthetic data sets
  temp.wide = reshape2::dcast(subset(post.long, Year >= 2005 & Year <= 2014), Alpha + Sample + Sex + Year ~ Stage, value.var="Share")
  temp.augm = dplyr::left_join(temp.wide, cart.weight, by=c("Alpha", "Sex", "Year")) # Augment with sample sizes by sex, year, and country
  temp.augm$Weight[is.na(temp.augm$Weight)] = 0 # Missing values correspond to no participants

  pred.raw = rmultinomial(nrow(temp.augm), size=temp.augm$Weight, prob=temp.augm[,cat.group])
  colnames(pred.raw) = cat.group
  pred.wide = cbind(temp.augm[,c("Alpha", "Sample", "Sex", "Year")], pred.raw)
  pred.flat = reshape2::melt(pred.wide, id.vars=c("Alpha", "Sample", "Sex", "Year"), measure.vars=cat.group, variable.name="Stage", value.name="Count")

  ## Aggregate across country and sex
  pred.aggr = plyr::ddply(pred.flat, .(Sample, Year, Stage), function(df) {
    data.frame(Count = sum(df$Count))
  })

  ## Calculate proportions
  pred.prop = plyr::ddply(pred.aggr, .(Sample, Year), function(df) {
    data.frame(Stage = df$Stage, Share = df$Count / sum(df$Count))
  })

  ## Calculate 95% interquantile range of posterior predictions
  pred.summ = plyr::ddply(pred.prop, .(Year, Stage), function(df) {
    bound = quantile(df$Share, c(0.025, 0.975))
    data.frame(Lower=bound[1], Upper=bound[2])
  })

  ## Join the posterior point estimate into the predictive interval
  pred.summ = dplyr::left_join(pred.summ, post.summ[,c("Year", "Stage", "Share")], by=c("Year", "Stage"))

  both.prop = dplyr::bind_rows(list(data = cart.prop, post = pred.summ), .id="Source")
  both.prop$Source = factor(both.prop$Source, levels=tiny.names, labels=nice.names)
  both.prop$Stage = factor(both.prop$Stage, levels=rev(cat.group))

  ggplot(subset(both.prop, Year>=2005 & Year<=2014), aes(x=Year, y=100*Share, ymin=100*Lower, ymax=100*Upper, color=Source, fill=Source)) +
    geom_line() +
    geom_ribbon(alpha=0.2, color=NA) +
    facet_wrap(~Stage, nrow=1) +
    scale_color_manual(values=c("#000000", brewer.pal(8,"Dark2")[1:2])) +
    scale_fill_manual(values=c("#000000", brewer.pal(8,"Dark2")[1:2])) +
    xlab("Year") + ylab("Share of ART initiators, %") +
    xlim(2004,2015) +
    theme_bw() +
    theme(legend.position = "right",
          legend.margin = margin(t=0, b=0, unit="cm"),
          plot.margin = margin(t=0, b=0, l=0.05, r=0.05, unit="cm"),
          panel.border = element_rect(fill=NA, color="#000000"),
          strip.background = element_rect(fill="grey85", color="#000000"),
          text = element_text(size=10),
          axis.title = element_text(size=rel(1.0)),
          axis.text.x = element_text(color="#000000", size=rel(0.8)),
          axis.text.y = element_text(color="#000000", size=rel(0.8)))
  ggsave(tiff.name, compression="lzw", dpi=600, units="mm", width=180, height=60)

  cat(sprintf("\nObserved proportion initiating at CD4<50\n"))
  print(subset(both.prop, Year %in% c(2005, 2014) & Source=="Data"  & Stage=="CD4<50"))
  cat(sprintf("\nModeled proportion initiating at CD4<50\n"))
  print(subset(both.prop, Year %in% c(2005, 2014) & Source=="Model fit" & Stage=="CD4<50"))
  cat(sprintf("\nObserved proportion initiating at CD4>350 during 2005-2012\n"))
  print(range(subset(both.prop, Year >= 2005 & Year <= 2012 & Source=="Data"  & Stage=="CD4>350")$Share))
  cat(sprintf("Modeled proportion initiating at CD4>350 during 2005-2012\n"))
  print(range(subset(both.prop, Year >= 2005 & Year <= 2012 & Source=="Model fit" & Stage=="CD4>350")$Share))
}
