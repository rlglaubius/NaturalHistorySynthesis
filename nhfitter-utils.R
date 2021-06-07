hyper.names = c("mort.1", "mort.2", "mort.3", "mort.4", "prog.1", "prog.2", "prog.3", "dist.1", "dist.2", "dist.3", "astart")

## Take model hyperparameters as an unlabeled matrix or vector and return an
## annotated data frame
hyper.frame = function(hyper) {
  if (is.null(dim(hyper))) {
    theta = data.frame(t(hyper))
  } else {
    theta = data.frame(hyper)
  }
  colnames(theta) = hyper.names
  return(theta)
}

## Convert a vector of model hyperparameters to Spectrum input parameters
calc.parameters = function(hyper) {
  ## Median CD4 count among seroconverters with CD4 > x cells/mm3
  conditional.median = function(x, scale, shape) {
    return(qfisk(1.0 - 0.5 * (1.0 - pfisk(x, scale, shape)), scale, shape))
  }
  
  ## recast the hyperparameter vector as a dataframe so we know which parameter
  ## is which
  p = hyper.frame(hyper)

  CD4 = 1:7
  age = c('15-24', '25-34', '35-44', '45+')
  sex = c('male', 'female')
  
  ## CD4 stage extents
  upper = c(Inf, 500, 350, 250, 200, 100, 50)
  lower = c(500, 350, 250, 200, 100,  50,  0)
  width = upper - lower
  width[1] = 185
  
  ## CD4 cell counts at seroconversion
  dist.shape = rep(p$dist.1 + 1.0, 8)                    # shape parameter by sex and age (indices 1-4, males; 5-8, females)
  dist.scale = rep(p$dist.2 * (1.0 - 0:3 * p$dist.3), 2) # scale = median CD4 count at seroconversion

  dist = array(NA, c(7,4,2), dimnames=list(CD4, age, sex))
  dist[,1,1] = pfisk(upper, shape=dist.shape[1], scale=dist.scale[1]) - pfisk(lower, shape=dist.shape[1], scale=dist.scale[1])
  dist[,2,1] = pfisk(upper, shape=dist.shape[2], scale=dist.scale[2]) - pfisk(lower, shape=dist.shape[2], scale=dist.scale[2])
  dist[,3,1] = pfisk(upper, shape=dist.shape[3], scale=dist.scale[3]) - pfisk(lower, shape=dist.shape[3], scale=dist.scale[3])
  dist[,4,1] = pfisk(upper, shape=dist.shape[4], scale=dist.scale[4]) - pfisk(lower, shape=dist.shape[4], scale=dist.scale[4])
  dist[,1,2] = pfisk(upper, shape=dist.shape[5], scale=dist.scale[5]) - pfisk(lower, shape=dist.shape[5], scale=dist.scale[5])
  dist[,2,2] = pfisk(upper, shape=dist.shape[6], scale=dist.scale[6]) - pfisk(lower, shape=dist.shape[6], scale=dist.scale[6])
  dist[,3,2] = pfisk(upper, shape=dist.shape[7], scale=dist.scale[7]) - pfisk(lower, shape=dist.shape[7], scale=dist.scale[7])
  dist[,4,2] = pfisk(upper, shape=dist.shape[8], scale=dist.scale[8]) - pfisk(lower, shape=dist.shape[8], scale=dist.scale[8])
  
  ## HIV-related mortality rates by CD4 cell count
  mort.shape = rep(p$mort.1, 8)
  mort.rate  = rep(p$mort.4 * p$mort.2 * (1.0 + 0:3 * p$mort.3), 2)
  
  upper[1] = 685 # Need a finite upper bound on stage CD4>500 when calculating mortality 
  mort = array(NA, c(7,4,2), dimnames=list(CD4, age, sex))
  mort[,1,1] = mort.rate[1] * ((mort.shape[1]^upper - mort.shape[1]^lower) / (upper - lower)) / log(mort.shape[1])
  mort[,2,1] = mort.rate[2] * ((mort.shape[2]^upper - mort.shape[2]^lower) / (upper - lower)) / log(mort.shape[2])
  mort[,3,1] = mort.rate[3] * ((mort.shape[3]^upper - mort.shape[3]^lower) / (upper - lower)) / log(mort.shape[3])
  mort[,4,1] = mort.rate[4] * ((mort.shape[4]^upper - mort.shape[4]^lower) / (upper - lower)) / log(mort.shape[4])
  mort[,1,2] = mort.rate[5] * ((mort.shape[5]^upper - mort.shape[5]^lower) / (upper - lower)) / log(mort.shape[5])
  mort[,2,2] = mort.rate[6] * ((mort.shape[6]^upper - mort.shape[6]^lower) / (upper - lower)) / log(mort.shape[6])
  mort[,3,2] = mort.rate[7] * ((mort.shape[7]^upper - mort.shape[7]^lower) / (upper - lower)) / log(mort.shape[7])
  mort[,4,2] = mort.rate[8] * ((mort.shape[8]^upper - mort.shape[8]^lower) / (upper - lower)) / log(mort.shape[8])

  ## CD4 cell count declines
  prog.shape = rep(1.0 / (p$prog.1 + 1.0), 8)
  prog.scale = rep(p$prog.2 * (1.0 - 0:3 * p$prog.3), 2) # time when CD4 count would hit zero

  ## We calculate stage duration for, say, CD4 200-250 as time between reaching
  ## CD4=250 and reaching CD4=200. We don't have an upper bound for CD4>500, so
  ## we use the median CD4 count among seroconverters who start out at CD4>500.
  ## We could use a fixed upper limit instead, but that can attenuate some of
  ## the effects of age on disease progression. By contrast, we do use a fixed
  ## upper bound for the mortality rate calculation above because mortality
  ## rates at CD4>500 are so low that the effect of doing something more
  ## elaborate is negligible.
  prog = array(NA, c(6,4,2), dimnames=list(1:6, age, sex))
  upper[1] = conditional.median(500, dist.scale[1], dist.shape[1]); prog[,1,1] = 1.0 / (prog.scale[1] * ((upper[1:6] / upper[1])^prog.shape[1] - (lower[1:6] / upper[1])^prog.shape[1]))
  upper[1] = conditional.median(500, dist.scale[2], dist.shape[2]); prog[,2,1] = 1.0 / (prog.scale[2] * ((upper[1:6] / upper[1])^prog.shape[2] - (lower[1:6] / upper[1])^prog.shape[2]))
  upper[1] = conditional.median(500, dist.scale[3], dist.shape[3]); prog[,3,1] = 1.0 / (prog.scale[3] * ((upper[1:6] / upper[1])^prog.shape[3] - (lower[1:6] / upper[1])^prog.shape[3]))
  upper[1] = conditional.median(500, dist.scale[4], dist.shape[4]); prog[,4,1] = 1.0 / (prog.scale[4] * ((upper[1:6] / upper[1])^prog.shape[4] - (lower[1:6] / upper[1])^prog.shape[4]))
  upper[1] = conditional.median(500, dist.scale[5], dist.shape[5]); prog[,1,2] = 1.0 / (prog.scale[1] * ((upper[1:6] / upper[1])^prog.shape[5] - (lower[1:6] / upper[1])^prog.shape[5]))
  upper[1] = conditional.median(500, dist.scale[6], dist.shape[6]); prog[,2,2] = 1.0 / (prog.scale[2] * ((upper[1:6] / upper[1])^prog.shape[6] - (lower[1:6] / upper[1])^prog.shape[6]))
  upper[1] = conditional.median(500, dist.scale[7], dist.shape[7]); prog[,3,2] = 1.0 / (prog.scale[3] * ((upper[1:6] / upper[1])^prog.shape[7] - (lower[1:6] / upper[1])^prog.shape[7]))
  upper[1] = conditional.median(500, dist.scale[8], dist.shape[8]); prog[,4,2] = 1.0 / (prog.scale[4] * ((upper[1:6] / upper[1])^prog.shape[8] - (lower[1:6] / upper[1])^prog.shape[8]))
  
  return(list(mort=mort, prog=prog, dist=dist, art=p$astart))
}

## Simulate untreated survival in seroconverter cohorts
## fp     Spectrum projection inputs
## param  Natural history parameters
eval.surv = function(fp, param) {
  y_init = fp$ss$proj_start
  y_span = fp$ss$PROJ_YEARS
  y_year = y_init:(y_init + y_span - 1)
  
  y_eval_init = 1997 # median time of entry in ALPHA Network cohorts
  y_eval_span = 21
  y_eval = y_eval_init:(y_eval_init + y_eval_span - 1)
  
  ## Insert a new infection cohort at 1997 (1 person per sex and single age
  ## 15,16,...,60). For ART-naive cohort simulation, assume new infections are
  ## uniformly distributed within four age groups (15-24, 25-34, 35-44, 45+)
  fp$new_inf = array(0, dim=c(66, 2, y_span), dimnames=list(15:80, c('male', 'female'), y_init:(y_init + y_span - 1)))
  fp$new_inf[15:60-14,,y_eval_init - y_init + 1] = 1

  ## Reaggregate natural history parameters for EPP-ASM age groups
  fp$cd4_initdist = param$dist[,c(1,1,1,2,2,3,3,4,4),]
  fp$cd4_mort     = param$mort[,c(1,1,1,2,2,3,3,4,4),]
  fp$cd4_prog     = param$prog[,c(1,1,1,2,2,3,3,4,4),]
  
  ## Calculate the number surviving in each age cohort by year
  hivpop = array(cohort_sim(fp), c(7, 66, 2, y_span), dimnames=list(1:7, 15:80, c('male', 'female'), y_year))
  
  ## Survival in ten-year age bands, used for visualizing survival curves
  prop_surv = array(0, c(4, y_eval_span), dimnames=list(c('15-24', '25-34', '35-44', '45-54'), y_eval))
  for (k in 1:y_eval_span) {
    y_ind = k + y_eval_init - y_init
    prop_surv[1,k] = sum(hivpop[,15:24+(k-1)-14,,y_ind])
    prop_surv[2,k] = sum(hivpop[,25:34+(k-1)-14,,y_ind])
    prop_surv[3,k] = sum(hivpop[,35:44+(k-1)-14,,y_ind])
    prop_surv[4,k] = sum(hivpop[,45:54+(k-1)-14,,y_ind])
  }
  
  ## Cache mortality rates
  mort_hiv = param$mort[,rep(1:4, c(10, 10, 10, 36)),]
  mort_nat = -log(fp$Sx)
  
  ## Mortality density and surviving proportions are kept by sex and single age
  ## (15-59) at infection for likelihood evaluation
  mort_density = array(0.0, c(7, 45, 2, y_eval_span), dimnames=list(1:7, 15:59, c('male', 'female'), y_eval))
  size_cohort  = array(0.0, c(7, 45, 2, y_eval_span), dimnames=list(1:7, 15:59, c('male', 'female'), y_eval))
  
  mort_nat_pop = hivpop * array(matrix(mort_nat, nrow=7, ncol=66 * 2 * y_span, byrow=TRUE), dim=c(7, 66, 2, y_span))
  mort_hiv_pop = hivpop * array(mort_hiv, dim=c(7, 66, 2, y_span))
  
  for (k in 1:y_eval_span) {
    k_rel = k
    k_abs = k + y_eval_init - y_init
    a_cur = 15:59+(k_rel-1)-14
    size_cohort[,,,k_rel]  = hivpop[,a_cur,,k_abs]
    mort_density[,,,k_rel] = mort_nat_pop[,a_cur,,k_abs] + mort_hiv_pop[,a_cur,,k_abs]
  }
  
  return(list(
    prop_surv = prop_surv,                         # matrix: row=age at infection, col=year(15-24, 25-34, 35-44, 45+) x (1997:2017)
    size_cohort = apply(size_cohort, 2:4, sum),    # array: i1=age at infection, i2=sex, i3=year (15:59 x 1:2 x 1997:2017)
    mort_density = apply(mort_density, 2:4, sum))) # array: i1=age at infection, i2=sex, i3=year (15:59 x 1:2 x 1997:2017)
}


eval.cd4.dist = function(fp, param, year) {
  du = 0.1 ## Simulation time step
  
  ## Reaggregate natural history parameters for EPP-ASM age groups
  fp$cd4_initdist = param$dist[,c(1,1,1,2,2,3,3,4,4),]
  fp$cd4_mort     = param$mort[,c(1,1,1,2,2,3,3,4,4),]
  fp$cd4_prog     = param$prog[,c(1,1,1,2,2,3,3,4,4),]
  
  fp$median_cd4init = rep(0, length(fp$median_cd4init)) # Disable the median CD4 at ART initiation for treatment allocation
  fp$med_cd4init_input = rep(0L, length(fp$med_cd4init_input))
  fp$scale_cd4_mort = 1L                                # Enable reductions in mortality proportional to ART coverage
  fp$art_alloc_method = 3L                              # Enable ART allocation weighted between eligibility & expected mortality
  fp$art_alloc_mxweight = param$art
  
  mod = first90::simmod(fp)
  
  y_ind = year - fp$ss$proj_start + 1
  
  ## age groups in attr(mod, "hivpop") are:
  ## 1=15-16, 2=17-19, 3=20-24
  ## 4=25-29, 5=30-34, 6=35-39
  ## 7=40-44, 8=45-49, 9=50-80
  hivpop = attr(mod, "hivpop")                                         # alias the untreated PLHIV population
  cd4age = apply(hivpop[,1:8,,y_ind], c(1,2), sum)                     # aggregate by sex in the target year, drop ages 50+
  cd4age = rbind(cd4age[1:2,], colSums(cd4age[3:4,]), cd4age[5:7,])    # merge 200-249 & 250-349 categories
  cd4age = cbind(rowSums(cd4age[,1:2]), cd4age[,3:8], rowSums(cd4age)) # merge ages 15-16 & 17-19, include 15-49
  dist = sweep(cd4age, 2, colSums(cd4age), FUN='/')                    # normalize to get CD4 distributions by age
  
  dimnames(dist) = list(c(">500","350-500", "200-349", "100-199", "50-99", "<50"),
                        c(sprintf("%d-%d", seq(15,45,5), seq(19,49,5)), '15-49'))
  
  ## Return the distribution by ages 15-19, 20-24, ..., 45-49 and 15-49
  return(dist)
}

eval.art.dist = function(fp, param, years) {
  du = 0.1 ## Simulation time step
  
  ## Reaggregate natural history parameters for EPP-ASM age groups
  fp$cd4_initdist = param$dist[,c(1,1,1,2,2,3,3,4,4),]
  fp$cd4_mort     = param$mort[,c(1,1,1,2,2,3,3,4,4),]
  fp$cd4_prog     = param$prog[,c(1,1,1,2,2,3,3,4,4),]
  
  fp$median_cd4init = rep(0, length(fp$median_cd4init)) # Disable the median CD4 at ART initiation for treatment allocation
  fp$med_cd4init_input = rep(0L, length(fp$med_cd4init_input))
  fp$scale_cd4_mort = 1L                                # Enable reductions in mortality proportional to ART coverage
  fp$art_alloc_method = 3L                              # Enable ART allocation weighted between eligibility & expected mortality
  fp$art_alloc_mxweight = param$art

  mod = first90::simmod(fp)
  sim.years = fp$ss$proj_start:(fp$ss$proj_start + fp$ss$PROJ_YEARS - 1)
  y.ind = which(sim.years %in% years)

  artpop = attr(mod, "artpop")
  cd4art = apply(artpop[1,,,,y.ind], c(1,3,4), sum) # cd4art[h,s,t] is the number of PLHIV on ART for <6 months by CD4 h, sex s, and year t

  dist = list(m = rbind(cd4art[1:2,1,], colSums(cd4art[3:4,1,]), cd4art[5:7,1,]), # merge 200-249 and 250-349 categories
              f = rbind(cd4art[1:2,2,], colSums(cd4art[3:4,2,]), cd4art[5:7,2,]))

  rval = lapply(dist, function(p) {
    prop = sweep(p, 2, colSums(p), FUN='/')
    colnames(prop) = sim.years[y.ind]
    rownames(prop) = c(">500","350-500", "200-349", "100-199", "50-99", "<50")
    return(prop[,])
  })
  
  return(rval)
}
