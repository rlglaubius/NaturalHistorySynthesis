library(ggfortify)

source("Utils/load-spec-defaults.R")

## Use cohort model to calculate survival after seroconversion for each
## posterior sample. This returns two survival curves for each sample,
## corresponding to adjusted and unadjusted HIV-related mortality rates.
calculate.surv.sims = function(imis.rval, pjnz.info) {
  draws.raw = data.frame(imis.rval$resample)
  draws.raw$mort.4 = 1 # turn off mortality rate ratio
  
  param.raw = apply(draws.raw, 1, calc.parameters)
  param.adj = imis.rval$param
  
  surv.raw = lapply(param.raw, function(p) {eval.surv(pjnz.info, p)})
  surv.adj = lapply(param.adj, function(p) {eval.surv(pjnz.info, p)})
  
  return(list(raw = surv.raw, adj = surv.adj))
}

plot.post.surv = function(imis.rval, post.sims, surv.data, tiff.name) {
  map.ind = which.max(imis.rval$prior + imis.out$lhood)
  
  surv.age.names = c("15-24", "25-34", "35-44", "45-54")
  
  spec.param = load.spec.defaults()
  surv.spec = 0.05 * eval.surv(spec.info$UGA, spec.param)$prop_surv
  
  plot.data.surv = function(surv.raw, surv.adj, age) {
    surv.flat.raw = sapply(surv.raw, function(s) {0.05 * s$prop_surv[age,]})
    surv.flat.adj = sapply(surv.adj, function(s) {0.05 * s$prop_surv[age,]})
    
    surv.summ.adj = apply(surv.flat.adj, 1, function(flat) {quantile(flat, c(0.025, 0.975), na.rm=TRUE)})
    surv.summ.raw = apply(surv.flat.raw, 1, function(flat) {quantile(flat, c(0.025, 0.975), na.rm=TRUE)})
    
    col = c('#000000', brewer.pal(8,'Dark2')[c(1,2,4)])
    par(las=1, mar=c(2.5,3.0,1.0,1.0), cex=1)
    plot(c(0,20), c(0,1), cex=0, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
    axis(1, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(1), labels=axTicks(1), cex.axis=0.9)
    axis(2, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(2), labels=100*axTicks(2))
    lines(survfit(Surv(wait, death) ~ 1, data=subset(surv.data, agegr==levels(agegr)[age] & infecage < 55)), col=col[1])
    polygon(c(0:20, 20:0), c(surv.summ.raw[1,], rev(surv.summ.raw[2,])), col=sprintf("%s60", col[3]), border=NA)
    polygon(c(0:20, 20:0), c(surv.summ.adj[1,], rev(surv.summ.adj[2,])), col=sprintf("%s60", col[2]), border=NA)
    
    lines(0:20, surv.flat.adj[,map.ind], col=col[2], lwd=1.5)
    lines(0:20, surv.flat.raw[,map.ind], col=col[3], lwd=1.5)
    lines(0:20, surv.spec[age,], col=col[4], lwd=1.5)
    if (age == 1) {
      legend('bottomleft', inset=c(0.025,0.00), box.lty=0, bg=NA,
             legend=c("Observed survival", "Model (adjusted)", "Model (unadjusted)", "Spectrum"),
             lty=1, lwd=1.5, col=col, xpd=1)
    }
    title(ylab="Surviving population, %", line=1.75)
    title(xlab="Years since seroconversion", line=1.50)
    title(main=sprintf("Age %s", surv.age.names[age]))
    
    cat(sprintf("Median survival, age %s\n", surv.age.names[age]))
    cat(sprintf("\tAdjusted\t%0.1f (%0.1f-%0.1f) years\n",
                approx(surv.flat.adj[,map.ind], 0:20, 0.5, method="linear")$y,
                approx(surv.summ.adj[1,], 0:20, 0.5, method="linear")$y,
                approx(surv.summ.adj[2,], 0:20, 0.5, method="linear")$y))
    cat(sprintf("\tUnadjusted\t%0.1f (%0.1f-%0.1f) years\n",
                approx(surv.flat.raw[,map.ind], 0:20, 0.5, method="linear")$y,
                approx(surv.summ.raw[1,], 0:20, 0.5, method="linear")$y,
                approx(surv.summ.raw[2,], 0:20, 0.5, method="linear")$y))
  }
  
  tiff(tiff.name, units="mm", w=180, h=2*60, pointsize=8, compression="lzw", res=600)
  layout(matrix(1:4, nrow=2, byrow=TRUE))
  plot.data.surv(post.sims$raw, post.sims$adj, 1) # 15-24
  plot.data.surv(post.sims$raw, post.sims$adj, 2) # 25-34
  plot.data.surv(post.sims$raw, post.sims$adj, 3) # 35-44
  plot.data.surv(post.sims$raw, post.sims$adj, 4) # 45+
  dev.off()
}


plot.pred.surv = function(imis.rval, post.sims, surv.data, tiff.name, show.raw=FALSE, show.spec=FALSE) {
  map.ind = which.max(imis.rval$prior + imis.out$lhood)
  
  surv.age.names = c("15-24", "25-34", "35-44", "45-54")
  
  tiny.names = c("data", "adj", "raw", "spec")
  nice.names = c("Observed survival", "Model fit (adjusted)", "Model fit (unadjusted)", "Spectrum")
  
  ## Analyze data and reformat Kaplan-Meier curves for plotting
  surv.data.sset = subset(surv.data, infecage < 55)
  surv.data.sset = haven::zap_labels(surv.data.sset) # haven labels break some functions, so we remove them
  surv.data.sset$Age = factor(surv.data.sset$agegr, levels=levels(surv.data.sset$agegr), labels=surv.age.names)
  surv.data.km = survfit(Surv(wait, death) ~ Age, data=surv.data.sset)
  
  ## Workaround so that we can plot Kaplan-Meier curves alongside smooth survival
  ## curves, and ensure that plotted K-M curves start at t=0
  surv.data.dummy = data.frame(time=0.0, surv=1.0, lower=1.0, upper=1.0, strata=surv.age.names)
  surv.data.frame = fortify(surv.data.km)
  surv.data.augm = rbind(surv.data.dummy, surv.data.frame[,c("time", "surv", "lower", "upper", "strata")])
  colnames(surv.data.augm) = c("Year", "Alive", "Lower", "Upper", "Age")
  colnames(surv.data.dummy) = colnames(surv.data.augm)
  
  ## Reformat posterior simulations for plotting
  surv.post = lapply(post.sims, function(mort) { # mort is "raw" or "adj"
    mort.list = lapply(mort, function(s) {
      arr = 0.05 * s$prop_surv
      dimnames(arr) = list(Age=rownames(arr), Year=colnames(arr))
      return(as.data.frame.table(arr, responseName="Alive"))
    })
    mort.flat = dplyr::bind_rows(mort.list, .id="Sample")
    mort.flat$Sample = as.numeric(mort.flat$Sample)
    return(mort.flat)
  })
  
  surv.summ.list = lapply(surv.post, function(mort) {
    summ = plyr::ddply(mort, .(Age, Year), function(df) {
      ptEst = df$Alive[which(df$Sample==map.ind)]
      bound = quantile(df$Alive, c(0.025, 0.975))
      data.frame(Alive=ptEst, Lower=bound[1], Upper=bound[2])
    })
    summ$Year = as.numeric(as.character(summ$Year)) - 1997
    return(summ)
  })
  
  ## Calculate Spectrum survival
  spec.param = load.spec.defaults()
  surv.spec = 0.05 * eval.surv(spec.info$UGA, spec.param)$prop_surv
  dimnames(surv.spec) = list(Age=rownames(surv.spec), Year=colnames(surv.spec))
  surv.spec.long = as.data.frame.table(surv.spec, responseName="Alive")
  surv.spec.long$Year = as.numeric(as.character(surv.spec.long$Year)) - 1997
  
  ## Calculate posterior predictive intervals
  
  simulate.person = function(sex, age, surv.prob) {
    ## Our cohort model returns the proportion surviving annually after
    ## seroconversion. We sample the number of whole years survived from a
    ## discrete distribution, then calculate the time of death within that year as
    ## a uniform random variate. This is equivalent to assuming the survival
    ## function is piecewise linear. Note that we truncate survival curves and
    ## plots at 20 years after seroconversion.
    cdf = 1.0 - c(surv.prob[age - 15 + 1, sex, ], 0.0)
    yrs = 0:20
    return(sample(yrs, 1, prob=diff(cdf)) + runif(1, 0, 1))
  }
  
  surv.pred.adj.sims = sapply(post.sims$adj, function(surv.sim) {
    apply(surv.data.sset, 1, function(person) {
      sex = as.numeric(person['sex'])
      age = as.numeric(person['int_age'])
      simulate.person(sex, age, surv.sim$size_cohort)
    })
  })
  colnames(surv.pred.adj.sims) = sprintf("pred.%04d", 1:3000)
  
  surv.pred.raw.sims = sapply(post.sims$raw, function(surv.sim) {
    apply(surv.data.sset, 1, function(person) {
      sex = as.numeric(person['sex'])
      age = as.numeric(person['int_age'])
      simulate.person(sex, age, surv.sim$size_cohort)
    })
  })
  colnames(surv.pred.raw.sims) = sprintf("pred.%04d", 1:3000)
  
  surv.pred.sims = list(raw = surv.pred.raw.sims,
                        adj = surv.pred.adj.sims)
  
  ## For each mortality model (raw or adj) and age group, calculate a Kaplan-Meier
  ## curve for each synthetic data set. Interpolate to get values at annual
  ## intervals
  surv.pred.func = lapply(surv.pred.sims, function(mort) {
    age.list = lapply(surv.age.names, function(age) {
      age.ind = which(surv.data.sset$Age == age)
      apply(mort[age.ind,], 2, function(surv.times) {
        sf = survfit(Surv(surv.times) ~ 1)
        x = c(0, sf$time, 100)
        y = c(1, sf$surv,   0)
        return(approx(x, y, 0:20, method="linear")$y)
      })
    })
    names(age.list) = surv.age.names
    return(age.list)
  })
  
  ## Calculate summary statistics of posterior predicted Kaplan-Meier curves for
  ## each mortality model and age group
  surv.pred.summ = lapply(surv.pred.func, function(mort) {
    age.list = lapply(mort, function(age) {
      surv.stat = apply(age, 1, function(prop.alive) {quantile(prop.alive, c(0.025, 0.975))})
      surv.stat.mtrx = data.matrix(surv.stat)
      dimnames(surv.stat.mtrx) = list(Stat=c("Lower", "Upper"), Year=0:20)
      surv.stat.long = as.data.frame.table(surv.stat.mtrx, responseName="Alive")
      surv.stat.wide = reshape2::dcast(surv.stat.long, Year ~ Stat, value.var="Alive")
      return(surv.stat.wide)
    })
    return(dplyr::bind_rows(age.list, .id="Age"))
  })
  
  surv.pred = dplyr::bind_rows(surv.pred.summ, .id="Source")
  surv.pred$Year = as.numeric(as.character(surv.pred$Year))
  surv.pred$Alive = NA
  surv.pred$Source = factor(surv.pred$Source, levels=tiny.names, labels=nice.names)
  
  surv.each = dplyr::bind_rows(list(data = surv.data.dummy,
                                    raw  = surv.summ.list$raw,
                                    adj  = surv.summ.list$adj,
                                    spec = surv.spec.long), .id="Source")
  surv.each$Source = factor(surv.each$Source, levels=tiny.names, labels=nice.names)
  
  show.mask = c(TRUE, TRUE, show.raw, show.spec)
  surv.each = subset(surv.each, Source %in% nice.names[show.mask])
  surv.pred = subset(surv.pred, Source %in% nice.names[show.mask])
  col = c("#000000", brewer.pal(8,"Dark2")[c(1,2,4)])[show.mask]

  ggplot(surv.data.augm, aes(x=Year, y=100*Alive)) +
    geom_step() +
    geom_step(data=surv.data.augm, aes(x=Year, y=100*Lower), linetype="dashed", size=0.3, show.legend=FALSE) +
    geom_step(data=surv.data.augm, aes(x=Year, y=100*Upper), linetype="dashed", size=0.3, show.legend=FALSE) +
    geom_line(data=surv.each, aes(color=Source)) +
    geom_ribbon(data=surv.each, aes(ymin=100*Lower, ymax=100*Upper, fill=Source), alpha=0.6, show.legend=FALSE) +
    geom_ribbon(data=surv.pred, aes(ymin=100*Lower, ymax=100*Upper, fill=Source), alpha=0.4, show.legend=FALSE) +
    scale_color_manual(values=col) +
    scale_fill_manual(values=col) +
    xlim(0,20) +
    xlab("Years since seroconversion") +
    ylab("Surviving population, %") +
    facet_wrap(~Age, nrow=2) +
    theme_bw() +
    theme(legend.position = "top",
          legend.margin = margin(t=0, b=0, unit="cm"),
          plot.margin = margin(t=0, b=0, l=0.05, r=0.05, unit="cm"),
          panel.border = element_rect(fill=NA, color="#000000"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill="grey85", color="#000000"),
          text = element_text(size=10),
          axis.title = element_text(size=rel(1.0)),
          axis.text.x = element_text(color="#000000", size=rel(1.0)),
          axis.text.y = element_text(color="#000000", size=rel(1.0)))
  ggsave(tiff.name, compression="lzw", dpi=600, units="mm", width=180, height=2*60)
}

plot.median.surv = function(imis.rval, post.sims, tiff.name) {
  ## Load median survival estimates
  med.surv.data = read.csv("Data/median-survival.csv")
  med.surv.data$Upper[is.na(med.surv.data$Upper)] = Inf # replace unbounded upper limits

  map.ind = which.max(imis.rval$prior + imis.rval$lhood)
  
  ## Calculate modeled median survival times
  surv.post = lapply(post.sims, function(mort) { # mort is "raw" or "adj"
    mort.list = lapply(mort, function(s) {
      arr = 0.05 * s$prop_surv
      dimnames(arr) = list(Age=rownames(arr), Year=colnames(arr))
      return(as.data.frame.table(arr, responseName="Alive"))
    })
    mort.flat = dplyr::bind_rows(mort.list, .id="Sample")
    mort.flat$Sample = as.numeric(mort.flat$Sample)
    return(mort.flat)
  })
  
  surv.summ.list = lapply(surv.post, function(mort) {
    summ = plyr::ddply(mort, .(Age, Year), function(df) {
      ptEst = df$Alive[which(df$Sample==map.ind)]
      bound = quantile(df$Alive, c(0.025, 0.975))
      data.frame(Alive=ptEst, Lower=bound[1], Upper=bound[2])
    })
    summ$Year = as.numeric(as.character(summ$Year)) - 1997
    return(summ)
  })
  
  med.surv.post = lapply(surv.summ.list, function(mort) {
    med = plyr::ddply(mort, .(Age), function(df) {
      data.frame(Value = approx(df$Alive, 0:20, 0.5, method="linear")$y,
                 Lower = approx(df$Lower, 0:20, 0.5, method="linear")$y,
                 Upper = approx(df$Upper, 0:20, 0.5, method="linear")$y)
    })
  })

  # tiny.names = c("todd2007aids", "babiker2000lancet", "adj", "raw")
  # nice.names = c("Todd et al.", "Babiker et al.", "Model (adjusted)", "Model (unadjusted)")
  tiny.names = c("raw", "adj", "babiker2000lancet", "todd2007aids")
  nice.names = c("Model (unadjusted)", "Model (adjusted)", "Babiker et al.", "Todd et al.")

  med.surv.each = rbind(med.surv.data, dplyr::bind_rows(med.surv.post, .id="Source"))
  med.surv.each$Source = factor(med.surv.each$Source, levels=tiny.names, labels=nice.names)
  med.surv.each$Age = gsub("45\\+", "45-54", med.surv.each$Age) # Show 45+ survival from Todd 2007 AIDS with 45-54 estimates from other sources
  
  ggplot(med.surv.each, aes(x=Source, y=Value, ymin=Lower, ymax=Upper, group=Source)) +
    geom_pointrange(aes(color=Source, shape=Source), size=0.25) +
    facet_wrap(~Age, ncol=1, strip.position="right") +
    ylab("Median survival, years") +
    scale_color_manual(values=c(brewer.pal(8,"Dark2")[2:1], "#000000", "#000000")) +
    scale_shape_manual(values=c(16,16,15,15)) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "none",
          plot.margin = margin(t=0, b=0, l=0.05, r=0.05, unit="cm"),
          panel.border = element_rect(fill=NA, color="#000000"),
          strip.background = element_rect(fill="grey85", color="#000000"),
          text = element_text(size=10),
          panel.spacing = unit(0, "lines"),
          axis.title.x = element_text(size=rel(0.8)),
          axis.title.y = element_blank(),
          axis.text.x = element_text(color="#000000", size=rel(0.8)),
          axis.text.y = element_text(color="#000000", size=rel(0.8)))
  ggsave(tiff.name, compression="lzw", dpi=600, units="mm", width=80, height=60)
  
  print(med.surv.each)
}


  
  
  