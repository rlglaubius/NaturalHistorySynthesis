source("Utils/load-spec-defaults.R")

## Plot a comparison of initial CD4 category (with 95% confidence intervals) to
## posterior distributions of categories (with 95% credible intervals)
plot.post.dist = function(imis.rval, data.dist, tiff.name) {
  
  spec.param = load.spec.defaults()
  
  plot.data.dist = function(dist.post, dist.mode, sex, age) {
    data.conf = t(multinomialCI(data.dist[,age,sex], 0.05))
    
    post.flat = sapply(dist.post, function(param) {param[,age,sex]})
    post.cred = apply(post.flat, 1, function(prop) {quantile(prop, c(0.025, 0.975))})
    
    spec.dist = spec.param$dist[,age,sex]
    spec.flat = c(spec.dist[1:2], sum(spec.dist[3:4]), spec.dist[5:7]) # combine 200-249 and 250-349
    
    xdata = 1:6 - 0.15
    xpost = 1:6
    xspec = 1:6 + 0.15

    col = c('#000000', brewer.pal(8,'Dark2'))
    par(las=1, mar=c(2.5,3.0,1.0,1.0), cex=1)
    plot(c(0.5,6.5), c(0,1), cex=0, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
    axis(1, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=0:6 + 0.5, labels=rep("",7))
    axis(1, mgp=c(0.1, 0.3, 0.0), tck=-0.00, at=1:6, labels=c(">500", "350-500", "200-349", "100-199", "50-99", "<50"), cex.axis=0.95)
    axis(2, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(2), labels=100*axTicks(2))
    points(xdata, data.dist[,age,sex] / sum(data.dist[,age,sex]), pch=15, col=col[1], xpd=1)
    arrows(xdata, data.conf[1,], xdata, data.conf[2,], col=col[1], code=3, angle=90, length=0.01)
    points(xpost, dist.mode[,age,sex], pch=16, col=col[2], xpd=1)
    arrows(xpost, post.cred[1,], xpost, post.cred[2,], col=col[2], code=3, angle=90, length=0.01)
    points(xspec, spec.flat, pch=4, col=col[5], xpd=1)
    if (age == 1 & sex == 2) {
      legend('topright', inset=c(0.025,0.00), box.lty=0, bg=NA,
             legend=c('Observed distribution', 'Model fit', 'Spectrum'),
             pch=c(15,16,4), lty=c(1,1,NA), col=col[c(1,2,5)], xpd=1)
    }
    title(xlab=expression(paste("CD4 cell count (cells/", mm^3, ")")), line=1.5)
    title(ylab="Population, %", line=1.75)
    title(main=sprintf('Age %s %s', age.names[age], sex.names[sex]))
  }
  
  dist.post = lapply(imis.rval$param, function(param) {
    param.agg = array(NA, c(6,4,2))
    param.agg[,,1] = rbind(param$dist[1:2,,1], colSums(param$dist[3:4,,1]), param$dist[5:7,,1]) # aggregate 200-250 and 250-350
    param.agg[,,2] = rbind(param$dist[1:2,,2], colSums(param$dist[3:4,,2]), param$dist[5:7,,2])
    return(param.agg)
  })
  
  map.ind = which.max(imis.rval$prior + imis.rval$lhood)
  
  ## Compare CD4 counts at seroconversion to the data. Figure shows estimates
  ## for males, since inputs are the same for both sexes
  tiff(tiff.name, units="mm", w=180, h=3*60, pointsize=8, compression="lzw", res=600)
  layout(matrix(1:8, nrow=4, byrow=FALSE))
  plot.data.dist(dist.post, dist.post[[map.ind]], 1, 1) # 15-24 males
  plot.data.dist(dist.post, dist.post[[map.ind]], 1, 2) # 25-34
  plot.data.dist(dist.post, dist.post[[map.ind]], 1, 3) # 35-44
  plot.data.dist(dist.post, dist.post[[map.ind]], 1, 4) # 45+
  plot.data.dist(dist.post, dist.post[[map.ind]], 2, 1) # 15-24 females
  plot.data.dist(dist.post, dist.post[[map.ind]], 2, 2) # 25-34
  plot.data.dist(dist.post, dist.post[[map.ind]], 2, 3) # 35-44
  plot.data.dist(dist.post, dist.post[[map.ind]], 2, 4) # 45+
  dev.off()
}

## Plot a comparison of initial CD4 category (with 95% confidence intervals) to
## posterior distributions of categories (with posterior predictive
## distributions)
plot.pred.dist = function(imis.rval, dist.data, tiff.name) {
  src.names = c("500+", "350-499", "200-349", "100-199", "50-99", "<50") # CD4 category labels in data set
  cat.names = c(">500", "350-500", "200-349", "100-199", "50-99", "<50") # Standardized CD4 category names
  
  tiny.names = c("data", "post", "spec")
  nice.names = c("Data", "Model fit", "Spectrum")
  
  ## Reformat the data set and calculate sample sizes by age and sex
  dist.data.long = as.data.frame.table(data.dist)
  colnames(dist.data.long) = c("CD4", "Age", "Sex", "Count")
  dist.data.long$CD4 = factor(dist.data.long$CD4, levels=src.names, labels=cat.names)
  
  dist.data.size = plyr::ddply(dist.data.long, .(Sex, Age), function(df) {
    data.frame(Count = sum(df$Count))
  })

  dist.data.prop = plyr::ddply(dist.data.long, .(Age, Sex), function(df) {
    share = df$Count / sum(df$Count)
    bound = multinomialCI(df$Count, alpha=0.05)
    data.frame(CD4=df$CD4, Share=share, Lower=bound[,1], Upper=bound[,2])
  })
  
  ## Reformat the posterior parameter sets
  map.ind = which.max(imis.rval$prior + imis.rval$lhood)
  dist.post = lapply(imis.rval$param, function(param) {
    param.agg = array(NA, c(6,4,2))
    param.agg[,,1] = rbind(param$dist[1:2,,1], colSums(param$dist[3:4,,1]), param$dist[5:7,,1]) # aggregate 200-250 and 250-350 categories
    param.agg[,,2] = rbind(param$dist[1:2,,2], colSums(param$dist[3:4,,2]), param$dist[5:7,,2])
    dimnames(param.agg) = list(CD4=cat.names, Age=age.names, Sex=sex.names)
    return(param.agg)
  })
  
  dist.post.list = lapply(dist.post, function(dist) {
    dist.long = as.data.frame.table(dist)
  })
  dist.post.long = dplyr::bind_rows(dist.post.list, .id="Sample")
  dist.post.long$Sample = as.numeric(dist.post.long$Sample)
  
  # Calculate parameter point estimates and credible intervals
  dist.post.prop = plyr::ddply(dist.post.long, .(Age, Sex, CD4), function(df) {
    share = df$Freq[which(df$Sample == map.ind)]
    bound = quantile(df$Freq, c(0.025, 0.975))
    data.frame(Share = share, Lower = bound[1], Upper = bound[2])
  })
  
  ## Load and format Spectrum defaults
  spec.param = load.spec.defaults()
  dist.spec = array(NA, c(6,4,2))
  dist.spec[,,1] = rbind(spec.param$dist[1:2,,1], colSums(spec.param$dist[3:4,,1]), spec.param$dist[5:7,,1])
  dist.spec[,,2] = rbind(spec.param$dist[1:2,,2], colSums(spec.param$dist[3:4,,2]), spec.param$dist[5:7,,2])
  dimnames(dist.spec) = list(CD4=cat.names, Age=age.names, Sex=sex.names)
  
  dist.spec.prop = as.data.frame.table(dist.spec, responseName="Share")
  dist.spec.prop$Lower = NA
  dist.spec.prop$Upper = NA
  
  dist.each.prop = dplyr::bind_rows(list(
    data = dist.data.prop,
    post = dist.post.prop,
    spec = dist.spec.prop), .id="Source")
  dist.each.prop$Source = factor(dist.each.prop$Source, level=tiny.names, labels=nice.names)
  
  dist.pred.list = lapply(dist.post.list, function(dist.long) {
    ## Augment posterior sample with observed sample sizes
    dist.wide = reshape(dist.long, timevar="CD4", idvar=c("Sex", "Age"), direction="wide")
    dist.augm = dplyr::left_join(dist.wide, dist.data.size, by=c("Sex", "Age"))

    ## Generate a synthetic data set by drawing counts in each CD4 category
    ## conditional on sex and age group
    pred.list = apply(dist.augm, 1, function(group) {
      cols = sprintf("Freq.%s", cat.names)
      freq = rmultinom(1, size=group['Count'], prob=group[cols])
      data.frame(Sex=group['Sex'],
                 Age=group['Age'],
                 CD4=cat.names,
                 Count=freq, row.names=NULL)
    })
    pred.flat = dplyr::bind_rows(pred.list)
  })

  ## Summarize synthetic data sets
  dist.pred = dplyr::bind_rows(dist.pred.list, .id="Sample")
  dist.pred$CD4 = factor(dist.pred$CD4, levels=cat.names)
  dist.pred$Sample = as.integer(dist.pred$Sample)

  dist.pred.prop = plyr::ddply(dist.pred, .(Sample, Age, Sex), function(df) {
    data.frame(CD4=df$CD4, Share = df$Count / sum(df$Count), Source="post")
  })

  ## To align violin plots of the posterior predictive distribution with point
  ## estimates in the figure, we need to have dummy datasets for "data" and
  ## "spec" that have the same size as the dist.pred.prop
  data.dummy = dplyr::bind_rows(lapply(1:nrow(imis.rval$resample), function(j) {dist.data.prop}), .id="Sample")
  data.dummy$Lower = NULL
  data.dummy$Upper = NULL
  data.dummy$Source = "data"
  data.dummy$Sample = as.numeric(data.dummy$Sample)

  spec.dummy = dplyr::bind_rows(lapply(1:nrow(imis.rval$resample), function(j) {dist.spec.prop}), .id="Sample")
  spec.dummy$Lower = NULL
  spec.dummy$Upper = NULL
  spec.dummy$Source = "spec"
  spec.dummy$Sample = as.numeric(spec.dummy$Sample)

  violin.data = dplyr::bind_rows(list(data=data.dummy, post=dist.pred.prop, spec=spec.dummy))
  violin.data$Source = factor(violin.data$Source, level=tiny.names, labels=nice.names)
  
  ggplot(dist.each.prop, aes(x=CD4, y=100*Share, color=Source)) +
    geom_point(position=position_dodge(width=0.6), size=1, stroke=1, aes(color=Source, shape=Source)) +
    geom_pointrange(position=position_dodge(width=0.6), fatten=2, aes(ymin=100*Lower, ymax=100*Upper, color=Source, shape=Source), show.legend=FALSE) +
    geom_violin(data=violin.data, scale="width", position=position_dodge(0.6), color=NA, alpha=0.5, aes(x=CD4, y=100*Share, fill=Source), show.legend=FALSE) +
    facet_grid(Sex~Age) +
    scale_color_manual(values=c("#000000", brewer.pal(8,"Dark2")[c(1,4)])) +
    scale_fill_manual(values=c("#000000", brewer.pal(8,"Dark2")[c(1,4)])) +
    scale_shape_manual(values=c(15, 16, 4)) +
    xlab(expression(paste("CD4 cell count (cells/", mm^3, ")"))) +
    ylab("Population, %") +
    theme_bw() +
    theme(legend.position = "top",
          legend.margin = margin(t=0, b=0, unit="cm"),
          plot.margin = margin(t=0, b=0, l=0.05, r=0.05, unit="cm"),
          panel.border = element_rect(fill=NA, color="#000000"),
          strip.background = element_rect(fill="grey85", color="#000000"),
          text = element_text(size=10),
          axis.title = element_text(size=rel(1.0)),
          axis.text.x = element_text(color="#000000", size=rel(1.0), angle=45, hjust=1),
          axis.text.y = element_text(color="#000000", size=rel(1.0)))
  ggsave(tiff.name, compression="lzw", dpi=600, units="mm", width=180, height=1.5*60)
}
