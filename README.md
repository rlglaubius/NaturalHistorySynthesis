# NaturalHistorySynthesis
Multi-parameter synthesis of evidence on HIV natural history

## Description
The code in this repository was used to perform a Bayesian multi-parameter evidence synthesis to estimate rates of HIV disease progression and mortality. This analysis used data on [immune status among HIV seroconverters soon after HIV acquisition](https://journals.lww.com/aidsonline/Fulltext/2017/05150/Joint_estimation_of_CD4__cell_progression_and.3.aspx), [mortality rates among people living with HIV](https://academic.oup.com/jid/article/197/3/398/2908655), [survival times after HIV seroconversion](https://journals.lww.com/aidsonline/Fulltext/2007/11006/Time_from_HIV_seroconversion_to_death__a.8.aspx), [immune status among household survey respondents found living with untreated HIV infection](https://phia.icap.columbia.edu/), and [immune status of people living with HIV when they started antiretroviral treatment](https://academic.oup.com/cid/article/66/6/893/4823847). Please note that several of these datasets are confidential and are not included in this repository, except for datasets extracted from published study reports (Table 2 [here](https://academic.oup.com/jid/article/197/3/398/2908655) and Table 2 [here](https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(00)02061-4/ppt)).

A manuscript describing this work is under review. Please check back for details on that publication once it is available, including acknowledgements of the many people who contributed to this work.

## Requirements
This analysis is implemented in R and uses several packages available from CRAN to perform its analyses. These are listed in the R code files.

## How to use this code
This repository includes two key entry point scripts, nhfitter.R and visualize-fit.R. Other code files provide supporting functionality required by these entry point scripts. To run the analysis, point your R session working directory at this repository, then run

```
> fit.model = TRUE
> source("nhfitter.R")
> save.image("fit-ws.RData") # optional
```

The parameter estimation code is time-consuming to run, so in the snippet above we have set a flag `fit.model` to tell the `nhfitter.R` to perform model fitting when run. The next line will load all required functions and data, compile some required C++ code (`cohort-sim.cpp`), then run the parameter estimation algorithm, then run `visualize-fit.R` to analyze the resulting parameter estimates. If you source `nhfitter.R` with `fit.model` undefined or FALSE, the code will stop just before running the parameter estimation algorithm. The third line (`save.image("fit-ws.RData")`) saves the resulting R environment in case you would like to perform additional analyses later.

The second entry point, `visualize-fit.R`, is also time-consuming to run because it performs numerous simulations of HIV cohorts and epidemics to evaluate the implications of parameter estimates. If you have saved your workspace (e.g., by running `save.image("fit-ws.RData")` above), you may wish to run just the visualization code:

```
> load("fit-ws.RData")
> source("visualize-fit.R")
```

## Repository contents
### Entry points
- nhfitter.R main entry point, prepares and runs natural history parameter estimation
- visualize-fit.R secondary entry point, analyzes parameter estimates

### Utilities
- cohort-sim.cpp HIV disease progression model implementation used to simulate cohorts of people after HIV acquisition
- nhfitter-utils.R key utilities needed to transform and evaluate natural history parameter values
- prepare-data-cart.R prepares ART initiator data for analysis
- prepare-data-dist.R prepares data on immune status after HIV acquisition
- prepare-data-mort.R prepares data on mortality rates by CD4 category
- prepare-data-phia.R prepares data on immune status among PHIA respondents living with untreated HIV
- prepare-data-surv.R prepares data on survival after HIV seroconversion
- prepare-spec-info.R prepares [Spectrum](https://avenirhealth.org/software-spectrum.php) files used to simulate national HIV epidemics
- write-params.R saves parameter estimates in .csv and Excel formats

### Folders
- Data includes publicly-available data included in this analysis
- Utils provides additional utilities needed for analysis and visualization
