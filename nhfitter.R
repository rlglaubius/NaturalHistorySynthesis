library(first90)
library(plyr)
library(Rcpp)

source('Utils/fisk.R')
source('Utils/logimis.R')

source('nhfitter-utils.R')
sourceCpp('cohort-sim.cpp')   # load code that simulates untreated cohorts of PLHIV

## Data for likelihood evaluation
source('prepare-data-surv.R') # loads data.surv: ALPHA network data
source('prepare-data-mort.R') # loads data.mort: CASCADE mortality data
source('prepare-data-dist.R') # loads data.dist: CASCADE CD4 cell counts at seroconversion data
source('prepare-data-phia.R') # loads data.phia: PHIA data on CD4 among untreated respondents living with HIV
source('prepare-data-cart.R') # loads data.cart: IeDEA data on CD4 at ART initiation

source('prepare-spec-info.R')

## Background mortality assumed when fitting mortality rates to CASCADE data
## Europe 1985-1990 central mortality rates by five-year age group, from WPP 2019
wpp2019mort = array(NA, dim=c(6, 4, 2), dimnames=list(c('>500', '350-500', '200-350', '100-200', '50-100', '<50'),
                                                      c('15-24', '25-34', '35-44', '45-54'),
                                                      c('male', 'female')))
wpp2019mort[,,1] = matrix(c(0.5 * (0.00093 + 0.00147),  # average 15-19 and 20-24, male
                            0.5 * (0.00178 + 0.00222),  # average 25-29 and 30-34
                            0.5 * (0.00294 + 0.00411),  # average 35-39 and 40-44
                            0.5 * (0.00663 + 0.01039)), # average 45-49 and 50-54
                          nrow=6, ncol=4, byrow=TRUE)
wpp2019mort[,,2] = matrix(c(0.5 * (0.00037 + 0.00048),  # average 15-19 and 20-24, female
                            0.5 * (0.00057 + 0.00079),  # average 25-29 and 30-34
                            0.5 * (0.00117 + 0.00176),  # average 35-39 and 40-44
                            0.5 * (0.00282 + 0.00439)), # average 45-49 and 50-54
                          nrow=6, ncol=4, byrow=TRUE)

prior = function(hyper) {
  theta = hyper.frame(hyper)
  
  ## Priors on mortality rates are designed to reflect estimates used in
  ## Thembisa
  ssize = 256
  m = mean(c(0.9882, 0.9876, 0.9899, 0.9887, 0.9856))
  shape1 = m * (ssize - 2) + 1
  shape2 = ssize - shape1
  
  ## Depletion time (number of years until CD4=0) prior for a 15-24 year-old
  ## seroconverter. We assume a mode of 15 years and a median of 20 years, as
  ## this yields a wide prior with <5% probability of reaching CD4=0 within 8
  ## years after seroconversion, which seems reasonable since it usually takes
  ## at least 8 years to reach CD4<200, let alone CD4=0.
  mode.depl = 15 # mode
  half.depl = 20 # median
  logmu = log(half.depl)
  logsd = sqrt(log(half.depl) - log(mode.depl))

  ## Lodi 2011 CID intercept for age < 25 is 24.167^2 = 584 
  alpha = sqrt(585)
  
  ## Mortality rates  
  pmort = cbind(mort.1 = dbeta( theta$mort.1, shape1=shape1, shape2=shape2,  log=TRUE),
                mort.2 = dgamma(theta$mort.2, shape=0.7 * 4 + 1, scale=0.25, log=TRUE),
                mort.3 = dexp(  theta$mort.3, 1.0,                           log=TRUE),
                mort.4 = dgamma(theta$mort.4, shape=2.0, scale=1.0,          log=TRUE))
  term.mort = rowSums(pmort)
  
  ## Progression rates
  ## We model CD4 trends with a polynomial of order prog.1+1. We assume this is
  ## a low-order polynomial with order at least 1
  pprog = cbind(prog.1 = dexp(  theta$prog.1, 1.0,          log=TRUE),
                prog.2 = dlnorm(theta$prog.2, logmu, logsd, log=TRUE),
                prog.3 = dunif( theta$prog.3,  0, 1/3,      log=TRUE))
  term.prog = rowSums(pprog)
  
  ## CD4 at seroconversion
  ## We model this with a Fisk distribution. A shape parameter dist.1 makes it
  ## unlikely that the distribution has mode 0.
  pdist = cbind(
    dist.1 = dexp(  theta$dist.1, 1/3,                      log=TRUE),
    dist.2 = dgamma(theta$dist.2, shape=alpha, scale=alpha, log=TRUE),
    dist.3 = dunif( theta$dist.3, 0, 1/3,                   log=TRUE))
  term.dist = rowSums(pdist)
  
  ## Prior justification: Spectrum generally allocates ART based on the number
  ## of people and the expected number of deaths in each ART-eligible CD4
  ## category. A user-input value between 0 and 1, usually 0.5, controls how
  ## much weight is given to the number of people versus the expected number of
  ## deaths. There is no real-world data on this weighting, so we assume a
  ## non-informative prior.
  term.astart = dunif(theta$astart, 0, 1, log=TRUE)
  
  return(term.dist + term.mort + term.prog + term.astart)
}

sample.prior = function(n) {
  ## See comments in 'prior' for justification of these choices
  ssize = 256
  m = mean(c(0.9882, 0.9876, 0.9899, 0.9887, 0.9856))
  shape1 = m * (ssize - 2) + 1
  shape2 = ssize - shape1
  
  mode.depl = 15 # mode
  half.depl = 20 # median
  logmu = log(half.depl)
  logsd = sqrt(log(half.depl) - log(mode.depl))
  
  alpha = sqrt(585)
  
  cbind(
    mort.1 = rbeta(n, shape1=shape1, shape2=shape2),   # phi[1]; mortality shape parameter
    mort.2 = rgamma(n, shape=0.7 * 4 + 1, scale=0.25), # phi[2]; mortality at age 15-24 when CD4=0
    mort.3 = rexp(n, 1.0),                             # phi[3]; mortality age effect
    mort.4 = rgamma(n, shape=2, scale=1),              # phi[4]; mortality scale factor
    
    prog.1 = rexp(n, 1.0),                             # theta[1]; progression shape parameter
    prog.2 = rlnorm(n, logmu, logsd),                  # theta[2]; CD4 depletion time for 15-24 year-olds
    prog.3 = runif(n, 0, 1/3),                         # theta[3]; progression age effect
    
    dist.1 = rexp(n, 1/3),                             # psi[1]; initial CD4 shape parameter
    dist.2 = rgamma(n, shape=alpha, scale=alpha),      # psi[2]; initial CD4 median among seroconverters aged 15-24
    dist.3 = runif(n, 0, 1/3),                         # psi[3]; initial CD4 age effect

    astart = runif(n, 0, 1))                           # omega; ART initiation weight (0, eligibility; 1, mortality)
}

likelihood.dist = function(param) {
  lhood.helper = function(param) {
    param.agg = array(NA, c(6,4,2))
    param.agg[,,1] = rbind(param$dist[1:2,,1], colSums(param$dist[3:4,,1]), param$dist[5:7,,1]) # aggregate 200-250 and 250-350
    param.agg[,,2] = rbind(param$dist[1:2,,2], colSums(param$dist[3:4,,2]), param$dist[5:7,,2])
    sum(c(dmultinom(data.dist[,1,1], prob=param.agg[,1,1], log=TRUE),
          dmultinom(data.dist[,2,1], prob=param.agg[,2,1], log=TRUE),
          dmultinom(data.dist[,3,1], prob=param.agg[,3,1], log=TRUE),
          dmultinom(data.dist[,4,1], prob=param.agg[,4,1], log=TRUE),
          dmultinom(data.dist[,1,2], prob=param.agg[,1,2], log=TRUE),
          dmultinom(data.dist[,2,2], prob=param.agg[,2,2], log=TRUE),
          dmultinom(data.dist[,3,2], prob=param.agg[,3,2], log=TRUE),
          dmultinom(data.dist[,4,2], prob=param.agg[,4,2], log=TRUE)))
  }
  return(sapply(param, lhood.helper))
}

likelihood.mort = function(param, scale) {
  lhood.helper = function(theta) {
    prop.f =  499/3497 ## 14% of CASCADE patients were females (dunn2008jid)
    prop.m = 2998/3497 ## 86% male
    ## Assume excess mortality rate at CD4 200-350 is the weighted average of 200-250 and 250-350
    param.agg = array(NA, c(6,4,2))
    param.agg[,,1] = wpp2019mort[,,1] + rbind(theta$mort[1:2,,1], (2*theta$mort[3,,1] + theta$mort[4,,1]) / 3, theta$mort[5:7,,1])
    param.agg[,,2] = wpp2019mort[,,2] + rbind(theta$mort[1:2,,2], (2*theta$mort[3,,2] + theta$mort[4,,2]) / 3, theta$mort[5:7,,2])
    return(sum(dpois(data.mort$deaths, (prop.m * param.agg[,,1] + prop.f * param.agg[,,2]) * data.mort$pyears, log=TRUE)))
  }
  ## Since we assume mortality may be underascertained in CASCADE all-cause
  ## mortality data, we need to remove the scale factor when calculating the
  ## CASCADE likelihood
  pmort = mapply(function(p,r) {p$mort = p$mort / r; return(p)}, param, scale, SIMPLIFY=FALSE)
  return(sapply(pmort, lhood.helper))
}

likelihood.surv = function(param) {
  alpha.dead = as.data.table(subset(data.surv, death==1 & infecage < 60))
  alpha.lost = as.data.table(subset(data.surv, death==0 & infecage < 60))

  surv.helper = function(surv) {
    helper = function(data, model) {
      wait = data$wait
      wmin = floor(wait)
      wmax = wmin + 1

      imin = cbind(data$int_age - 14, data$sex, wmin + 1)
      imax = cbind(data$int_age - 14, data$sex, wmax + 1)

      ymin = model[imin]
      ymax = model[imax]

      return(sum(log(ymin * (wmax - wait) + ymax * (wait - wmin)), na.rm=TRUE))
    }

    ## term_dead and term_lost respectively accumulate the mortality density at
    ## time of death or the proportion surviving at time of censoring for each
    ## individual
    term.dead = helper(alpha.dead, surv$mort_density)
    term.lost = helper(alpha.lost, surv$size_cohort)
    return(term.dead + term.lost)
  }

  lhood = rep(-Inf, length(param))
  surv.model = lapply(param, function(p) {eval.surv(spec.info$UGA, p)})
  valid = which(sapply(surv.model, function(surv) {all(is.finite(surv$size_cohort)) & !any(surv$size_cohort < 0) & !any(surv$size_cohort > 1)}))
  if (length(valid) > 0) {
    lhood[valid] = sapply(surv.model[valid], surv.helper)
  }

  ## We weight of the survival data proportional to the number of countries
  ## evaluated, since we previously evaluated the survival likelihood once per
  ## country, but this was very slow and differences in background mortality
  ## weren't enough to yield meaningfully change the marginal likelihood across
  ## countries.
  return(lhood * length(spec.info))
  
  # lhood = matrix(-Inf, nrow=length(param), ncol=length(spec.info))
  # for (k in 1:length(spec.info)) {
  #   surv.model = lapply(param, function(p) {eval.surv(spec.info[[k]], p)})
  #   valid = which(sapply(surv.model, function(surv) {all(is.finite(surv$size_cohort)) & !any(surv$size_cohort < 0) & !any(surv$size_cohort > 1)}))
  #   if (length(valid) > 0) {
  #     lhood[valid,k] = sapply(surv.model[valid], surv.helper)
  #   }
  # }
  # return(rowSums(lhood))
}

likelihood.phia = function(param) {
  lhood = matrix(-Inf, nrow=length(param), ncol=length(spec.info))
  for (k in 1:length(spec.info)) {
    data.sset = subset(data.phia, Alpha==metadata$alpha[k])
    data.dist = data.matrix(data.sset[1:7,10:5]) # get data by 5-year age (15-49) and CD4 category (in decreasing order)
    data.size = data.sset$N[1:7]
    data.year = data.sset$Year[1]

    ## Evaluate each parameter set
    calc.phia = lapply(param, function(p) {eval.cd4.dist(spec.info[[k]], p, data.year)})
    
    ## Filter out parameter sets that produce invalid epidemics
    valid = which(sapply(calc.phia, function(dist) {all(is.finite(dist)) & all(dist >= 0)}))
    if (length(valid) > 0) {
      ## Evaluate the likelihood of the valid parameter sets
      lhood[valid,k] = sapply(calc.phia[valid], function(dist) {
        sum(dmultinom(data.dist[1,] * data.size[1], prob=dist[,1], log=TRUE), # 15-19
            dmultinom(data.dist[2,] * data.size[2], prob=dist[,2], log=TRUE), # 20-24
            dmultinom(data.dist[3,] * data.size[3], prob=dist[,3], log=TRUE), # 25-29
            dmultinom(data.dist[4,] * data.size[4], prob=dist[,4], log=TRUE), # 30-34
            dmultinom(data.dist[5,] * data.size[5], prob=dist[,5], log=TRUE), # 35-39
            dmultinom(data.dist[6,] * data.size[6], prob=dist[,6], log=TRUE), # 40-44
            dmultinom(data.dist[7,] * data.size[7], prob=dist[,7], log=TRUE)) # 45-49
      })
    }
  }
  
  return(rowSums(lhood))
}

likelihood = function(hyper) {
  theta = hyper.frame(hyper)

  ## rates >= 10 events/yr are considered invalid because they can cause
  ## degenerate behavior in our model simulators
  param = apply(theta, 1, calc.parameters)
  valid = which(sapply(param, function(p) {all(p$dist >= 0) & all(p$prog >= 0) & all(p$prog < 10) & all(p$mort >= 0) & all(p$mort < 10)}))
  lhood = matrix(-Inf, nrow=nrow(theta), ncol=4)
  if (length(valid) > 0) {
    lhood[valid,1] = likelihood.dist(param[valid])
    lhood[valid,2] = likelihood.mort(param[valid], theta$mort.4[valid])
    lhood[valid,3] = likelihood.surv(param[valid])
    lhood[valid,4] = likelihood.phia(param[valid])
  }
  return(rowSums(lhood))
}

sex.names = c('Male', 'Female')
age.names = c('15-24', '25-34', '35-44', '45+')
cd4.names = c(">500","350-500", "250-349", "200-249", "100-199", "50-99", "<50")

if (exists("fit.model") && fit.model == TRUE) { # must set "fit.model=TRUE" before running script
  imis.out = IMIS.log(prior, likelihood, sample.prior, B = 1000, B.re = 3000, number_k = 400, D = 0)
  imis.out$prior = prior(imis.out$resample)
  imis.out$lhood = likelihood(imis.out$resample)
  imis.out$param = apply(imis.out$resample, 1, calc.parameters)
  post.mode.ind = which.max(imis.out$prior + imis.out$lhood)
  
  cat(sprintf("# params    %8.0f\n", ncol(imis.out$resample)))
  cat(sprintf("Lhood(Dist) %8.1f\n", likelihood.dist(imis.out$param[post.mode.ind])))
  cat(sprintf("Lhood(Mort) %8.1f\n", likelihood.mort(imis.out$param[post.mode.ind], imis.out$resample[post.mode.ind,4])))
  cat(sprintf("Lhood(Surv) %8.1f\n", likelihood.surv(imis.out$param[post.mode.ind])))
  cat(sprintf("Lhood(PHIA) %8.1f\n", likelihood.phia(imis.out$param[post.mode.ind])))
  cat(sprintf("Lhood       %8.1f\n", imis.out$lhood[post.mode.ind]))
  cat(sprintf("AIC         %8.1f\n", 2*ncol(imis.out$resample) - 2*imis.out$lhood[post.mode.ind]))
  
  save.image("fit-ws.RData") # save the workspace for later analysis and visualization
  source('write-params.R')
  source('visualize-fit.R')
} else {
  warning("Model fitting not performed")
}
