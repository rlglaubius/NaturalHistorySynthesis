library(ggplot2)
library(mc2d)
library(MultinomialCI)
library(plyr)
library(survival)
library(RColorBrewer)
library(Rcpp)

load("fit-ws.RData")
sourceCpp("cohort-sim.cpp")

source("Utils/plotutils.R")

source("Utils/plot-prior-post.R")  # Prior and posterior hyperparameter distribution comparison
source("Utils/plot-post-trends.R") # Continuous natural history model
source("Utils/plot-post-dist.R")   # Initial CD4 cell count distribution
source("Utils/plot-post-mort.R")   # All-cause mortality rates
source("Utils/plot-post-prog.R")   # HIV disease progression rates
source("Utils/plot-post-surv.R")   # Survival after seroconversion
source("Utils/plot-post-phia.R")   # CD4 counts in untreated PLHIV PHIA respondents
source("Utils/plot-post-cart.R")   # CD4 counts at ART initiation

## Precalculate survival curves and CD4 distributions in the PHIAs and at ART
## initiation for each posterior sample from IMIS
surv.sims = calculate.surv.sims(imis.out, spec.info$UGA)
phia.sims = calculate.phia.sims(imis.out, spec.info, data.phia)
cart.sims = calculate.cart.sims(imis.out, spec.info, data.cart)

plot.prior.post(imis.out$resample, "fitted-params.tiff")

tiff("fitted-trends.tiff", units="mm", w=80, h=3*60, pointsize=8, compression="lzw", res=600)
layout(matrix(1:3, nrow=3, byrow=TRUE))
plot.trend.dist(imis.out, lwd=1.5, col=brewer.pal(8,'Dark2')[1:4]); legend('topleft', inset=c(-0.18,-0.05), title=expression(bold('A')), legend=c(''), xpd=1, box.lty=0, bg=NA, cex=1.25)
plot.trend.mort(imis.out, lwd=1.5, col=brewer.pal(8,'Dark2')[1:4]); legend('topleft', inset=c(-0.18,-0.05), title=expression(bold('B')), legend=c(''), xpd=1, box.lty=0, bg=NA, cex=1.25)
plot.trend.prog(imis.out, lwd=1.5, col=brewer.pal(8,'Dark2')[1:4]); legend('topleft', inset=c(-0.18,-0.05), title=expression(bold('C')), legend=c(''), xpd=1, box.lty=0, bg=NA, cex=1.25)
dev.off()

## Goodness-of-fit plots showing credible intervals
plot.post.dist(imis.out, data.dist, "plot-post-dist.tiff")
plot.post.mort(imis.out, data.mort, "plot-post-mort.tiff")
plot.post.surv(imis.out, surv.sims, data.surv, "plot-post-surv.tiff")
plot.post.phia(imis.out, phia.sims, data.phia.long, "plot-post-phia.tiff")
plot.post.cart(imis.out, cart.sims, data.cart, "plot-post-cart.tiff")
plot.post.prog(imis.out, "plot-post-prog.tiff")
plot.median.surv(imis.out, surv.sims, "plot-post-median-surv.tiff")

## Goodness-of-fit plots showing posterior predictive intervals
plot.pred.dist(imis.out, data.dist, "plot-pred-dist.tiff")
plot.pred.mort(imis.out, data.mort, "plot-pred-mort.tiff", show.raw=TRUE, show.spec=TRUE)
plot.pred.surv(imis.out, surv.sims, data.surv, "plot-pred-surv.tiff", show.raw=TRUE, show.spec=TRUE)
plot.pred.phia(imis.out, phia.sims, data.phia.long, "plot-pred-phia.tiff")
plot.pred.cart(imis.out, cart.sims, data.cart, "plot-pred-cart.tiff")
plot.pred.prog(imis.out, "plot-pred-prog.tiff")



