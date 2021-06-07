library(ggplot2)

plot.prior.post = function(post.sample, tiff.name) {
  num.breaks = 30
  hyper.name = c("astart", "dist.1", "dist.2", "dist.3", "mort.1", "mort.2", "mort.3", "mort.4", "prog.1", "prog.2", "prog.3")
  hyper.expr = c("omega",  "psi[1]", "psi[2]", "psi[3]", "varphi[1]", "varphi[2]", "varphi[3]", "varphi[4]", "theta[1]", "theta[2]", "theta[3]")
  hyper.brks = list(astart = seq(0, 1, 1/(num.breaks-1)),
                    dist.1 = seq(0, 20, 20/(num.breaks-1)),
                    dist.2 = seq(200, 1100, 900/(num.breaks-1)),
                    dist.3 = seq(0, 1/3, (1/3)/(num.breaks-1)),
                    mort.1 = seq(0.94, 1.00, 0.06/(num.breaks-1)),
                    mort.2 = seq(0, 4, 4/(num.breaks-1)),
                    mort.3 = seq(0, 8, 8/(num.breaks-1)),
                    mort.4 = seq(0, 8, 8/(num.breaks-1)),
                    prog.1 = seq(0, 8, 8/(num.breaks-1)),
                    prog.2 = seq(0, 150, 150/(num.breaks-1)),
                    prog.3 = seq(0, 1/3, (1/3)/(num.breaks-1)))

  hyper.wide = dplyr::bind_rows(list(
    Prior     = as.data.frame(sample.prior(nrow(post.sample))),
    Posterior = as.data.frame(post.sample)),
    .id="Distribution")
  hyper.long = reshape(hyper.wide,
                       varying=colnames(post.sample),
                       v.names="Value",
                       timevar="hyper",
                       times=colnames(post.sample),
                       direction="long")
  hyper.long$Distribution = factor(hyper.long$Distribution, levels=c("Prior", "Posterior"))
  hyper.long$hyper = factor(hyper.long$hyper, levels=hyper.name, labels=hyper.expr)
  
  ## Configure breaks across multiple histograms, per
  ## https://stackoverflow.com/questions/17271968/different-breaks-per-facet-in-ggplot2-histogram
  hist.list = mapply(function(x, b) {geom_histogram(data=x, breaks=b, aes(x=Value, fill=Distribution), alpha=0.8, position="identity")},
                     plyr::dlply(hyper.long, .(hyper)), hyper.brks)
  
  ggplot(hyper.long, aes(x=Value, fill=Distribution)) +
    hist.list +
    facet_wrap(~hyper, scales="free_x", label="label_parsed") +
    xlab("Parameter value") +
    ylab("Frequency") +
    theme_bw() +
    theme(legend.position = "top",
          legend.margin = margin(t=0, b=0, unit="cm"),
          plot.margin = margin(t=0, b=0, l=0.05, r=0.05, unit="cm"),
          panel.border = element_rect(fill=NA, color="#000000"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          strip.background = element_rect(fill="grey85", color="#000000"),
          text = element_text(size=10),
          axis.title = element_text(size=rel(1.0)),
          axis.text.x = element_text(color="#000000", size=rel(0.8)),
          axis.text.y = element_text(color="#000000", size=rel(0.8)))
  ggsave(tiff.name, dpi=600, compression="lzw", units="mm", width=180, height=135)
}
