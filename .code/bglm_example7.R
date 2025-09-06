## -----------------------------------------------------------------------------
#| label: setup
#| include: false

knitr::opts_chunk$set(cache.lazy = FALSE,
                      tidy = "styler")
options(tinytex.engine = "xelatex")


## -----------------------------------------------------------------------------
#| label: libraries
#| output: false
#| eval: true
#| warning: false
#| message: false
#| cache: false

library(tidyverse)     #for data wrangling etc
library(rstanarm)      #for fitting models in STAN
library(cmdstanr)      #for cmdstan
library(brms)          #for fitting models in STAN
library(standist)      #for exploring distributions
library(coda)          #for diagnostics
library(bayesplot)     #for diagnostics
library(ggmcmc)        #for MCMC diagnostics
library(DHARMa)        #for residual diagnostics
library(rstan)         #for interfacing with STAN
library(emmeans)       #for marginal means etc
library(broom)         #for tidying outputs
library(tidybayes)     #for more tidying outputs
library(HDInterval)    #for HPD intervals
library(ggeffects)     #for partial plots
library(broom.mixed)   #for summarising models
library(posterior)     #for posterior draws
library(ggeffects)     #for partial effects plots
library(patchwork)     #for multi-panel figures
library(bayestestR)    #for ROPE
library(see)           #for some plots
library(easystats)     #framework for stats, modelling and visualisation
library(geoR)     #framework for stats, modelling and visualisation
library(modelsummary)  #for data and model summaries
theme_set(theme_grey()) #put the default ggplot theme back
source('helperFunctions.R')


## ----readData, results='markdown', eval=TRUE----------------------------------
freitas <- read_csv("../data/freitas.csv", trim_ws = TRUE)


## -----------------------------------------------------------------------------
#| label: examinData
glimpse(freitas)


## -----------------------------------------------------------------------------
## Explore the first 6 rows of the data
head(freitas)


## -----------------------------------------------------------------------------
str(freitas)


## -----------------------------------------------------------------------------
freitas |> datawizard::data_codebook()


## -----------------------------------------------------------------------------
freitas |> modelsummary::datasummary_skim()


## ----dataProcessing, results='markdown', eval=TRUE, mhidden=TRUE--------------
freitas <- freitas |>
  filter(FISH == "Cod_9044") |>
  droplevels() |>
  mutate(DAY = as.numeric(as.factor(DATE)),
         DEPTH = DEPTH_MEAN_DAY)


## ----eda0, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=7, fig.height=5----
freitas |> ggplot(aes(y = DEPTH)) +
  geom_boxplot()
## freitas |> ggplot(aes(y = Dist)) +
##   geom_boxplot()
freitas |> ggplot(aes(y = DEPTH)) +
  geom_boxplot() +
  scale_y_log10()
freitas |> ggplot(aes(y = TEMPERATURE)) +
  geom_boxplot() +
  scale_y_log10()
freitas |> ggplot(aes(y = TEMPERATURE)) +
  geom_boxplot()


## ----eda1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=7, fig.height=5----
## freitas |> ggplot(aes(y = Dist, x = Day)) +
##   geom_point() +
##   geom_line(aes(group = ID))
freitas |> ggplot(aes(y = DEPTH, x = TEMPERATURE)) +
  geom_point() +
  geom_line(aes(colour = DATE))


## ----fitModel2a, results='markdown', eval=TRUE, mhidden=TRUE, paged.print=FALSE, tidy.opts = list(width.cutoff = 80)----
freitas_form <- bf(I(DEPTH) ~ scale(TEMPERATURE),
                family=gaussian(link = 'identity'))
options(width=150)
freitas_form |> get_prior(data = freitas)
options(width=80)


## ----fitModel2h, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
freitas |>
  summarise(
    Int_mu = median(DEPTH),
    Int_sd = mad(DEPTH),
    b_sd = mad(DEPTH))

freitas |>
  mutate(DEPTH = log(DEPTH)) |>
  summarise(
    Int_mu = median(DEPTH),
    Int_sd = mad(DEPTH),
    b_sd = mad(DEPTH))


## ----fitModel2h1, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE-----
priors <- prior(normal(14, 4.5), class = 'Intercept') +
    prior(normal(0, 4.5), class = 'b') +
    prior(student_t(3, 0, 4.5), class = 'sigma')
freitas_form <- bf(DEPTH ~ scale(TEMPERATURE),
                family=gaussian())
## freitas_form <- bf(Dist ~ scale(Day) + scale(TOD),
##                 family=Gamma(link = "log"))
get_prior(freitas_form, data =  freitas)
## priors <- prior(normal(4.2, 0.2), class = 'Intercept') +
##     prior(normal(0, 0.2), class = 'b') +
##   prior(gamma(0.01, 0.01), class = 'shape')
freitas_brm2 <- brm(freitas_form,
                 data = freitas,
                 prior = priors,
                 sample_prior = 'only',
                 iter = 5000,
                 warmup =2500,
                 chains = 3,
                 cores = 3,
                 thin = 10,
                 refresh = 100,
                 seed = 123,
                 control =  list(adapt_delta = 0.99, max_treedepth = 20),
                 backend = "cmdstanr"
                 )



## ----partialPlot2h1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
freitas_brm2 |>
    conditional_effects("TEMPERATURE") |>
    plot(points = TRUE)



## ----fitModel2h1b, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
freitas_brm3 <- update(freitas_brm2,
                       sample_prior = 'yes',
                       chains =  3, cores =  3,
                       refresh = 100)
save(freitas_brm3, file = '../ws/testing/freitas_brm3')


## ----partialPlot2h1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
freitas_brm3 |>
    conditional_effects("TEMPERATURE") |>
  plot(points = TRUE)


freitas_brm3 |> SUYR_prior_and_posterior()


## ----modelValidation2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
freitas_brm3$fit |> stan_trace()
freitas_brm3$fit |> stan_trace(inc_warmup=TRUE)


## ----modelValidation2h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
freitas_brm3$fit |> stan_ac()


## ----modelValidation2i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
freitas_brm3$fit |> stan_rhat()


## ----modelValidation2j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
freitas_brm3$fit |> stan_ess()


## ----modelValidation2k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
freitas_brm3$fit |> stan_dens(separate_chains = TRUE)


## ----modelValidation5a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation5b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
freitas_brm3 |> pp_check(type = 'dens_overlay', ndraws = 100)


## ----modelValidation5c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#freitas_brm3 |> pp_check(type = 'error_scatter_avg')


## ----modelValidation5e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
freitas_brm3 |> pp_check(type='intervals')
## freitas_brm3 |> pp_check(group='DENSITY', type='intervals')


## ----modelValidation5g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(freitas_brm3)


## ----modelValidation6aa, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=10----
freitas_resids <- make_brms_dharma_res(freitas_brm3, integerResponse = FALSE)
wrap_elements(~testUniformity(freitas_resids)) +
  wrap_elements(~plotResiduals(freitas_resids, form = factor(rep(1, nrow(freitas))))) +
  wrap_elements(~plotResiduals(freitas_resids, quantreg = TRUE)) +
  wrap_elements(~testDispersion(freitas_resids))

freitas_resids |> testTemporalAutocorrelation(time = freitas$DAY)

## freitas_resid1 <- freitas_resids |>
##   recalculateResiduals(group = freitas$Day, aggregateBy = mean)
##   ## recalculateResiduals(group = interaction(freitas$Day,  freitas$ID),  aggregateBy = mean)
## freitas_resid1 |> testTemporalAutocorrelation(time=unique(freitas$Day))
autocor_check(freitas, freitas_brm3, variable =  "DAY", n.sim =  250)
residuals(freitas_brm3)[, "Estimate"] |> acf()


## ----fitModel4, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE-------
freitas_form <- bf(DEPTH ~ scale(TEMPERATURE) +
                      ar(time = DAY, p = 1, cov =  TRUE),
                    family =  gaussian())
get_prior(freitas_form, data =  freitas)
priors <- prior(normal(14, 4.5), class = 'Intercept') +
  prior(normal(0, 4.5), class = 'b') +
  prior(student_t(3, 0, 4.5), class = 'sigma') +
  prior(normal(0, 1), class = 'ar', lb = -1, ub = 1)

freitas_brm4 <- brm(freitas_form,
                 data = freitas,
                 prior = priors,
                 sample_prior = 'yes',
                 iter = 5000,
                 warmup =2500,
                 chains = 3,
                 cores = 3,
                 thin = 10,
                 refresh = 1000,
                 seed = 123,
                 control =  list(adapt_delta = 0.99, max_treedepth = 20),
                 backend = "cmdstanr"
                 )

save(freitas_brm4, file = '../ws/testing/freitas_brm4')

## ----partialPlot7h1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
freitas_brm4 |>
    conditional_effects("TEMPERATURE") |>
  plot(points = TRUE) |>
  _[[1]]
## freitas_brm4 |> SUYR_prior_and_posterior()


## ----modelValidation4g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
freitas_brm4$fit |> stan_trace()


## ----modelValidation4h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
freitas_brm4$fit |> stan_ac()


## ----modelValidation4i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
freitas_brm4$fit |> stan_rhat()


## ----modelValidation4j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
freitas_brm4$fit |> stan_ess()


## ----modelValidation4k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
freitas_brm4$fit |> stan_dens(separate_chains = TRUE)


## ----modelValidation6b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
freitas_brm4 |> pp_check(type = 'dens_overlay', ndraws = 100)


## ----modelValidation6ab, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=10----
freitas_resids <- make_brms_dharma_res(freitas_brm4, integerResponse = FALSE)
wrap_elements(~testUniformity(freitas_resids)) +
               wrap_elements(~plotResiduals(freitas_resids, form = factor(rep(1, nrow(freitas))))) +
               wrap_elements(~plotResiduals(freitas_resids, quantreg = TRUE)) +
               wrap_elements(~testDispersion(freitas_resids))


## ----modelValidation6ab2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=4----
freitas_resids |> testTemporalAutocorrelation(time = freitas$DAY)

autocor_check(freitas, freitas_brm3, variable =  "DAY", n.sim =  250)
autocor_check(freitas, freitas_brm4, variable =  "DAY", n.sim =  250)
residuals(freitas_brm3)[, "Estimate"] |> acf()
residuals(freitas_brm4)[, "Estimate"] |> acf()


## -----------------------------------------------------------------------------
#| label: name
#| results: markup
#| eval: true
#| echo: true
#| cache: false

freitas_brm3 |> augment() |>
  ## group_by(ID) |>
  reframe(ACF = as.numeric(acf(.resid, plot =  FALSE)$acf)) |>
  ## group_by(ID) |>
  mutate(N =  1:n()) |>
  ## slice(1:12) |>
  ungroup() |>
  group_by(N) |>
  summarise(ACF = mean(ACF)) |>
  ggplot(aes(y =  ACF, x =  N)) +
  geom_hline(yintercept =  0) +
  geom_segment(aes(yend =  0, xend =  N))
freitas_brm4 |> augment() |>
  ## group_by(ID) |>
  reframe(ACF = as.numeric(acf(.resid, plot =  FALSE)$acf)) |>
  ## group_by(ID) |>
  mutate(N =  1:n()) |>
  ## slice(1:12) |>
  ungroup() |>
  group_by(N) |>
  summarise(ACF = mean(ACF)) |>
  ggplot(aes(y =  ACF, x =  N)) +
  geom_hline(yintercept =  0) +
  geom_segment(aes(yend =  0, xend =  N))


