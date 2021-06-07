source("Utils/load-spec-defaults.R")

plot.post.prog = function(imis.rval, tiff.name) {
  col = c('#000000', brewer.pal(8,'Dark2'))[c(1, 2, 3, 5)]
  pch = c(15, 16, 16, 4)
  spec.param = load.spec.defaults()
  map.ind = which.max(imis.rval$prior + imis.rval$lhood)
  
  plot.prog.rate = function(prog.post, sex, age, legend.flag) {
    post.flat = sapply(prog.post, function(param) {param[,age,sex]})
    post.cred = apply(post.flat, 1, function(rate) {quantile(rate, c(0.025, 0.975))})
    spec.flat = spec.param$prog[,age,sex]
    
    xpost = 1:6 - 0.075
    xspec = 1:6 + 0.075
    par(las=1, mar=c(2.5,3.0,1.0,1.0), cex=1)
    plot(c(0.5,6.5), c(0.0,1.2), cex=0, ann=FALSE, axes=FALSE, xaxs='i', yaxs='i')
    points(xspec, spec.flat, pch=pch[4], col=col[4])
    points(xpost, post.flat[,map.ind], pch=pch[2], col=col[2])
    arrows(xpost, post.cred[1,], xpost, post.cred[2,], col=col[2], code=3, angle=90, length=0.01)
    if (legend.flag) {legend("topright", inset=c(0.025,0.00), box.lty=0, bg=NA,
                             legend=c("Model fit", "Spectrum"),
                             pch=pch[c(2,4)], lty=c(1,NA), col=col[c(2,4)], xpd=1)}
    axis(1, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=0:6 + 0.5, labels=rep("",7))
    axis(1, mgp=c(0.1, 0.3, 0.0), tck=-0.00, at=1:6, labels=c(">500", "350-500", "250-349", "200-249", "100-199", "50-99"), cex.axis=0.95)
    axis(2, mgp=c(0.1, 0.3, 0.0), tck=-0.01, at=axTicks(2), labels=100*axTicks(2))
    title(xlab=expression(paste("CD4 cell count (cells/", mm^3, ")")), line=1.5)
    title(ylab='Progression rate, events/100 PY', line=1.75)
    title(main=sprintf('Age %s', age.names[age]))
  }
  
  prog.post = lapply(imis.rval$param, function(param) {param$prog})
  
  ## Compare progression rates by CD4 cell count
  tiff(tiff.name, units="mm", w=180, h=2*60, pointsize=8, compression="lzw", res=600)
  layout(matrix(1:4, nrow=2, byrow=TRUE))
  plot.prog.rate(prog.post, 1, 1, legend=TRUE ) # 15-24 males
  plot.prog.rate(prog.post, 1, 2, legend=FALSE) # 25-34
  plot.prog.rate(prog.post, 1, 3, legend=FALSE) # 35-44
  plot.prog.rate(prog.post, 1, 4, legend=FALSE) # 45+
  dev.off()
}

plot.pred.prog = function(imis.rval, tiff.name) {
  map.ind = which.max(imis.rval$prior + imis.rval$lhood)
  spec.param = load.spec.defaults()
  cat.names = c(">500", "350-500", "250-349", "200-249", "100-199", "50-99")
  
  tiny.names = c("post", "spec")
  nice.names = c("Model fit", "Spectrum")
  
  post.prog.list = lapply(imis.rval$param, function(param) {
    prog = param$prog
    dimnames(prog) = list(CD4=cat.names, Age=age.names, Sex=sex.names)
    return(as.data.frame.table(prog, responseName="Rate"))
  })
  post.prog.flat = dplyr::bind_rows(post.prog.list, .id="Sample")
  post.prog.flat$Sample = as.numeric(post.prog.flat$Sample)
  post.prog.flat$CD4 = factor(post.prog.flat$CD4, levels=cat.names)
  
  post.prog.summ = plyr::ddply(post.prog.flat, .(CD4, Age, Sex), function(df) {
    ptEst = df$Rate[which(df$Sample == map.ind)]
    bound = quantile(df$Rate, c(0.025, 0.975))
    data.frame(Rate=ptEst, Lower=bound[1], Upper=bound[2])
  })
  
  spec.prog = spec.param$prog
  dimnames(spec.prog) = list(CD4=cat.names, Age=age.names, Sex=sex.names)
  spec.prog = as.data.frame.table(spec.prog, responseName="Rate")
  spec.prog$CD4 = factor(spec.prog$CD4, levels=cat.names)
  spec.prog$Lower = NA
  spec.prog$Upper = NA
  
  prog.each = dplyr::bind_rows(list(post = post.prog.summ, spec=spec.prog), .id="Source")
  prog.each$Source = factor(prog.each$Source, levels=tiny.names, labels=nice.names)
  
  ggplot(subset(prog.each, Sex=="Male"), aes(x=CD4, y=100*Rate, color=Source)) +
    geom_point(position=position_dodge(width=0.3), size=1, stroke=1, aes(color=Source, shape=Source)) +
    geom_pointrange(position=position_dodge(width=0.3), fatten=2, aes(ymin=100*Lower, ymax=100*Upper, color=Source, shape=Source), show.legend=FALSE) +
    facet_wrap(~Age, nrow=2) +
    scale_color_manual(values=brewer.pal(8,"Dark2")[c(1,4)]) +
    scale_shape_manual(values=c(16, 4)) +
    xlab(expression(paste("CD4 cell count (cells/", mm^3, ")"))) +
    ylab("Progression rate, events/100 PY") +
    ylim(0,120) +
    theme_bw() +
    theme(legend.position = "top",
          legend.margin = margin(t=0, b=0, unit="cm"),
          plot.margin = margin(t=0, b=0, l=0.05, r=0.05, unit="cm"),
          panel.border = element_rect(fill=NA, color="#000000"),
          strip.background = element_rect(fill="grey85", color="#000000"),
          text = element_text(size=10),
          axis.title = element_text(size=rel(1.0)),
          axis.text.x = element_text(color="#000000", size=rel(1.0)), # angle=45, hjust=1),
          axis.text.y = element_text(color="#000000", size=rel(1.0)))
  ggsave(tiff.name, compression="lzw", dpi=600, units="mm", width=180, height=2*60)
}









