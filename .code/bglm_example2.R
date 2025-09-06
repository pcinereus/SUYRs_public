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
library(coda)          #for diagnostics
library(bayesplot)     #for diagnostics
library(ggmcmc)        #for MCMC diagnostics
library(DHARMa)        #for residual diagnostics
library(rstan)         #for interfacing with STAN
library(emmeans)       #for marginal means etc
library(broom)         #for tidying outputs
library(tidybayes)     #for more tidying outputs
library(ggeffects)     #for partial plots
library(broom.mixed)   #for summarising models
library(ggeffects)     #for partial effects plots
library(bayestestR)    #for ROPE
library(see)           #for some plots
library(easystats)     #for the easystats ecosystem
library(INLA)          #for approximate Bayes
library(INLAutils)     #for additional INLA outputs
library(patchwork)     #for multiple plots
library(modelsummary)  #for data and model summaries
theme_set(theme_grey()) #put the default ggplot theme back
source("helperFunctions.R")


## -----------------------------------------------------------------------------
#| label: readData
#| output: true
#| eval: true
polis <- read_csv("../data/polis.csv", trim_ws = TRUE)


## -----------------------------------------------------------------------------
#| label: examinData
glimpse(polis)


## -----------------------------------------------------------------------------
## Explore the first 6 rows of the data
head(polis)


## -----------------------------------------------------------------------------
str(polis)


## -----------------------------------------------------------------------------
polis |> datawizard::data_codebook()


## -----------------------------------------------------------------------------
polis |> modelsummary::datasummary_skim()


## ----EDA, results='markdown', eval=TRUE, mhidden=TRUE-------------------------
ggplot(polis, aes(y=PA, x=RATIO))+
  geom_point()
ggplot(polis, aes(y=PA, x=RATIO))+
  geom_point()+
  geom_smooth(method='glm', formula=y~x,
              method.args=list(family='binomial'))


## ----lm, results='markdown', eval=TRUE, mhidden=TRUE--------------------------
summary(glm(PA ~ RATIO, data = polis, family = binomial()))
summary(glm(PA ~ center(RATIO), data = polis, family = binomial()))


## ----fitModel1a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
polis.rstanarm = stan_glm(PA ~ RATIO, data=polis,
                          family=binomial(),
                         iter = 5000, warmup = 1000,
                         chains = 3, thin = 5, refresh = 0)


## ----fitModel1b, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
prior_summary(polis.rstanarm)


## ----fitModel1c, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
mean(polis$PA)


## ----fitModel1d, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
2.5 * 1/sd(polis$RATIO)


## ----fitModel1f, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
polis.rstanarm1 <- update(polis.rstanarm,  prior_PD=TRUE)

## ----fitModel1g, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
ggemmeans(polis.rstanarm1,  ~RATIO) |> plot(show_data=TRUE)


## ----fitModel1h, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
polis.rstanarm2= stan_glm(PA ~ RATIO, data=polis,
                          family=binomial(),
                          prior_intercept = normal(0, 2.5, autoscale=FALSE),
                          prior = normal(0, 0.1, autoscale=FALSE),
                          prior_PD=TRUE,
                          iter = 5000, warmup = 1000,
                          chains = 3, thin = 5, refresh = 0
                          )


## ----fitModel1i, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
ggemmeans(polis.rstanarm2,  ~RATIO) |>
  plot(show_data=TRUE)


## ----fitModel1j, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
polis.rstanarm3= update(polis.rstanarm2,  prior_PD=FALSE)


## ----modelFit1k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
posterior_vs_prior(polis.rstanarm3, color_by='vs', group_by=TRUE,
                   facet_args=list(scales='free_y'))


## ----modelFit1l, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggemmeans(polis.rstanarm3,  ~RATIO) |> plot(show_data=TRUE)
ggemmeans(polis.rstanarm3, terms = "RATIO[0:63]") |>
    plot(show_data=TRUE,
         show_residuals = TRUE,
         jitter = FALSE)
#OR
polis.rstanarm3 |> ggpredict(~RATIO) |> plot(show_data=TRUE)


## ----fitModel2a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
form <- bf(PA|trials(1) ~ RATIO, family = binomial())
#OR
form <- bf(PA ~ RATIO, family = bernoulli())
polis_brm <- brm(form,
                data = polis,
                iter = 5000,
                warmup = 1000,
                chains = 3, cores = 3,
                thin = 5,
                refresh = 0,
                backend = 'cmdstanr')


## ----fitModel2b, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE, paged.print=FALSE,tidy.opts = list(width.cutoff = 80), echo=2----
options(width=100)
polis_brm |> prior_summary()
options(width=80)


## ----fitModel2d, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
form <- bf(PA|trials(1) ~ RATIO, family = binomial())
#OR
polis_form <- bf(PA ~ RATIO, family = bernoulli())
priors <- prior(normal(0, 2.5), class = 'Intercept') +
    prior(normal(0, 1), class = 'b')
polis_brm1 = brm(polis_form,
                 data = polis,
                 prior = priors,
                 sample_prior = 'only',
                 iter = 5000,
                 warmup = 1000,
                 chains = 3, cores = 3,
                 thin = 5,
                 refresh = 0,
                 backend = 'cmdstanr')


## ----fitModel2e, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
polis_brm1 |> ggemmeans(~RATIO) |> plot(show_data=TRUE)
polis_brm1 |> conditional_effects() |>  plot(points=TRUE)


## ----normal2h, results='markdown', eval=TRUE----------------------------------
standist::visualize("normal(0, 2.5)", xlim = c(-10, 10))


## ----normal2h2, results='markdown', eval=TRUE---------------------------------
dat <- data.frame(sigma = c(2.5, 2, 1.5, 1))
dat <- dat |>
  group_by(sigma) |>
  reframe(r = rnorm(10000, 0, sigma),
          p = plogis(r))
ggplot(dat, aes(x =  p)) +
  geom_density(aes(fill = factor(sigma)), alpha =  0.3)


## ----fitModel2h, results='markdown', eval=TRUE, mhidden=TRUE------------------
polis_form <- bf(PA|trials(1) ~ RATIO, family = binomial())
#OR
polis_form <- bf(PA ~ RATIO, family = bernoulli())
priors <- prior(normal(0, 1),  class = 'Intercept') +
    prior(normal(0, 0.1), class = 'b')

polis_brm2 <- brm(polis_form,
                 data = polis,
                 prior = priors,
                 sample_prior = 'only',
                 iter = 5000,
                 warmup = 1000,
                 chains = 3, cores = 3,
                 thin = 5,
                 refresh = 0,
                 backend = 'cmdstanr')


## ----fitModel2i, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
polis_brm2 |> conditional_effects() |>
  plot(points = TRUE)

ggemmeans(polis_brm2,  ~RATIO) |>
  plot(show_data = TRUE)



## ----fitModel2j, results='markdown', eval=TRUE, mhidden=TRUE------------------
polis_brm3 <- polis_brm2 |> update(sample_prior = 'yes', refresh = 0)


## ----fitModel2j1, results='markdown', eval=TRUE, echo = FALSE, mhidden=TRUE----
save(polis_brm3, file = '../ws/testing/polis_brm3.RData')


## ----fitModel2j2, results='markdown', eval=TRUE, mhidden=TRUE-----------------
polis_brm3 |>
  conditional_effects() |>
  plot(points = TRUE)

## To see the fitted trend on the scale of the link function
## E.g. to see that it is a straight line
polis_brm3 |> conditional_effects(method = "posterior_linpred")

ggemmeans(polis_brm3,  ~RATIO) |>
  plot(show_data = TRUE)
ggemmeans(polis_brm3,  terms = "RATIO[0:63]") |>
    plot(show_data = TRUE)



## ----posterior2k, results='markdown', eval=TRUE-------------------------------
polis_brm3 |> get_variables()
## polis_brm3 |> hypothesis('Intercept=0', class='b') |> plot()
polis_brm3 |> hypothesis('RATIO=0') |> plot()


## ----fitModel2k, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
polis_brm3 |> get_variables()
polis_brm3 |> SUYR_prior_and_posterior()


## ----fitModel2l, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
polis_brm3 |> standata()
polis_brm3 |> stancode()


## ----INLApackages, results='markdown', eval=TRUE------------------------------
library(INLA)


## ----fitModel3a, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
polis.inla <- inla(PA ~ RATIO,
                  data = polis,
                  family = 'gaussian',
                  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE)
                  )


## ----fitModel3a2, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE----
polis.inla |> names()


## ----fitModel3b0, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE, paged.print=FALSE,tidy.opts = list(width.cutoff = 80), echo=TRUE----
inla.priors.used(polis.inla)


## ----fitModel3b3, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE, paged.print=FALSE,tidy.opts = list(width.cutoff = 80)----
standist::visualize("gamma(1, 0.00005)", xlim=c(-100,100000))


## ----fitModel3b1, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE, paged.print=FALSE,tidy.opts = list(width.cutoff = 80)----
standist::visualize("normal(0, 31)", xlim=c(-100,100))


## ----fitModel3b9, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE, paged.print=FALSE,tidy.opts = list(width.cutoff = 80)----
mean(polis$PA)


## ----fitModel3b8, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE, paged.print=FALSE,tidy.opts = list(width.cutoff = 80)----
standist::visualize("normal(0, 1)", xlim=c(-2,2))


## ----fitModel3b4, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE, paged.print=FALSE,tidy.opts = list(width.cutoff = 80)----
standist::visualize("normal(0, 2)", xlim=c(-3,3))


## ----fitModel3b6, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE----
polis.inla1 <- inla(PA ~ scale(RATIO, scale=FALSE),
                  data = polis,
                  Ntrials = 1,
                  family = 'binomial',
                  control.fixed = list(
                      mean.intercept = 0,
                      prec.intercept = 0.01,
                      mean = 0,
                      prec = 0.25),
                  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE)
                  )


## ----modelValidation1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_mcmc()


## ----modelValidation1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(polis.rstanarm3, plotfun='mcmc_trace')


## ----modelValidation1c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(polis.rstanarm3, 'acf_bar')


## ----modelValidation1d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(polis.rstanarm3, 'rhat_hist')


## ----modelValidation1e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(polis.rstanarm3, 'neff_hist')


## ----Validation1f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(polis.rstanarm3, 'combo')
plot(polis.rstanarm3, 'violin')


## ----modelValidation1g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_trace(polis.rstanarm3)


## ----modelValidation1h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_ac(polis.rstanarm3)


## ----modelValidation1i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_rhat(polis.rstanarm3)


## ----modelValidation1j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_ess(polis.rstanarm3)


## ----modelValidation1k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_dens(polis.rstanarm3, separate_chains = TRUE)


## ----modelValidation1l, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
polis.ggs <- ggs(polis.rstanarm3)
ggs_traceplot(polis.ggs)


## ----modelValidation1m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_autocorrelation(polis.ggs)


## ----modelValidation1n, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_Rhat(polis.ggs)


## ----modelValidation1o, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_effective(polis.ggs)


## ----modelValidation1p, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_crosscorrelation(polis.ggs)


## ----modelValidation1q, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_grb(polis.ggs)


## ----modelValidation2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_mcmc()


## ----modelValidation2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
polis_brm3 |> mcmc_plot(type = 'trace')


## ----modelValidation2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
polis_brm3 |> mcmc_plot(type = 'acf_bar')


## ----modelValidation2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
polis_brm3 |> mcmc_plot(type = 'rhat_hist')


## ----modelValidation2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
polis_brm3 |> mcmc_plot(type = 'neff_hist')


## ----modelValidation2f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
polis_brm3 |> mcmc_plot(type = 'combo')
polis_brm3 |> mcmc_plot(type = 'violin')


## ----modelValidation2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
polis_brm3$fit |> stan_trace()


## ----modelValidation2h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
polis_brm3$fit |> stan_ac()


## ----modelValidation2i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
polis_brm3$fit |> stan_rhat()


## ----modelValidation2j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
polis_brm3$fit |> stan_ess()


## ----modelValidation2k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_dens(polis_brm3$fit, separate_chains = TRUE)


## ----modelValidation2l, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=7----
polis.ggs <- polis_brm3 |> ggs(inc_warmup = FALSE, burnin = FALSE)
polis.ggs |> ggs_traceplot()


## ----modelValidation2m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=7----
polis.ggs |> ggs_autocorrelation()


## ----modelValidation2n, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
polis.ggs |> ggs_Rhat()


## ----modelValidation2o, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
polis.ggs |> ggs_effective()


## ----modelValidation2p, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
polis.ggs |> ggs_crosscorrelation()


## ----modelValidation2q, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
polis.ggs |> ggs_grb()


## ----modelValidation3a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation3b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(polis.rstanarm3,  plotfun='dens_overlay')


## ----modelValidation3c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#pp_check(polis.rstanarm3, plotfun='error_scatter_avg')


## ----modelValidation3d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(polis.rstanarm3, x=polis$RATIO, plotfun='error_scatter_avg_vs_x')


## ----modelValidation3e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(polis.rstanarm3, x=polis$RATIO, plotfun='intervals')


## ----modelValidation3f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(polis.rstanarm3, x=polis$RATIO, plotfun='ribbon')


## ----modelValidation3g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(polis.rstanarm3)


## ----modelValidation4a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
preds <- posterior_predict(polis.rstanarm3,  ndraws=250,  summary=FALSE)
polis.resids <- createDHARMa(simulatedResponse = t(preds),
                            observedResponse = polis$PA,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse = TRUE)
plot(polis.resids)


## ----modelValidation5a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation5b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
polis_brm3 |> pp_check(type = 'dens_overlay', ndraws=100)


## ----modelValidation5c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## polis_brm3 |> pp_check(type = 'error_scatter_avg')


## ----modelValidation5d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
polis_brm3 |> pp_check(x = 'RATIO', type = 'error_scatter_avg_vs_x')


## ----modelValidation5e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
polis_brm3 |> pp_check(x = 'RATIO', type = 'intervals')


## ----modelValidation5f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
polis_brm3 |> pp_check(x = 'RATIO', type = 'ribbon')


## ----modelValidation5g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(polis_brm3)


## ----modelValidation6aa, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=10----
polis.resids <- make_brms_dharma_res(polis_brm3, integerResponse = FALSE)
wrap_elements(~testUniformity(polis.resids)) +
               wrap_elements(~plotResiduals(polis.resids, form = factor(rep(1, nrow(polis))))) +
               wrap_elements(~plotResiduals(polis.resids, quantreg = FALSE)) +
               wrap_elements(~testDispersion(polis.resids))



## ----modelValidation6a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
preds <- polis_brm3 |> posterior_predict(ndraws = 250,  summary = FALSE)
polis.resids <- createDHARMa(simulatedResponse = t(preds),
                            observedResponse = polis$PA,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse = TRUE)
polis.resids |> plot()



## ----partialPlot1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis.rstanarm3 |> ggpredict() |> plot(show_data=TRUE, jitter = FALSE)


## ----partialPlot1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis.rstanarm3 |> ggemmeans(~RATIO) |> plot(show_data=TRUE)


## ----partialPlot1c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis.rstanarm3 |> epred_draws(newdata=polis) |>
 median_hdci() |>
 ggplot(aes(x=RATIO, y=.epred)) +
 geom_ribbon(aes(ymin=.lower, ymax=.upper), fill='blue', alpha=0.3) +
 geom_line()


## ----partialPlot2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis_brm3 |>
    conditional_effects() |>
    plot(points = TRUE)
polis_brm3 |>
    conditional_effects(spaghetti = TRUE,ndraws = 500) |>
    plot(points = TRUE)


## ----partialPlot2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis_brm3 |>
    ggpredict() |>
    plot(show_data = TRUE, jitter = FALSE)


## ----partialPlot2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis_brm3 |>
    ggemmeans(~RATIO) |>
    plot(show_data = TRUE)


## ----partialPlot2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis_brm3 |> epred_draws(newdata = polis) |>
 median_hdci() |>
 ggplot(aes(x = RATIO, y = .epred)) +
 geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = 'blue', alpha = 0.3) +
 geom_line()

partial.obs <- polis |>
    mutate(fit = fitted(polis_brm3, newdata = polis)[,'Estimate'],
           resid = resid(polis_brm3)[,'Estimate'],
           Obs = fit + resid)

polis_brm3 |>
    epred_draws(newdata = polis) |>
    median_hdci() |>
    ggplot(aes(x = RATIO, y = .epred)) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = 'blue', alpha = 0.3) +
    geom_point(data = partial.obs, aes(y = Obs)) +
    geom_line()


## ----summariseModel1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
summary(polis.rstanarm3)


## ----summariseModel1a1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=FALSE----
polis.sum <- summary(polis.rstanarm3)


## ----summariseModel1dd, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis.rstanarm3$stanfit |>
    summarise_draws(median,
                    HDInterval::hdi,
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel1d2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis.rstanarm3$stanfit |>
    summarise_draws(median,
                    ~HDInterval::hdi(.x, credMass = 0.9),
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tidyMCMC(polis.rstanarm3$stanfit, estimate.method='median',  conf.int=TRUE,
         conf.method='HPDinterval',  rhat=TRUE, ess=TRUE)


## ----summariseModel1b1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
polis.tidy <- tidyMCMC(polis.rstanarm3$stanfit, estimate.method='median',  conf.int=TRUE,  conf.method='HPDinterval',  rhat=TRUE, ess=TRUE)


## ----summariseModel1m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis.rstanarm3$stanfit |> as_draws_df()

## summarised
polis.rstanarm3$stanfit |>
    as_draws_df() |>
    summarise_draws("median",
                    ~ HDInterval::hdi(.x),
                    "rhat",
                    "ess_bulk")
## summarised on fractional scale
polis.rstanarm3$stanfit |>
    as_draws_df() |>
    mutate(across(everything(), exp)) |>
    summarise_draws("median",
                    ~ HDInterval::hdi(.x),
                    "rhat",
                    "ess_bulk")


## ----summariseModel1c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis.draw <- polis.rstanarm3 |> gather_draws(`(Intercept)`, RATIO)
## OR via regex
polis.draw <- polis.rstanarm3 |> gather_draws(`.Intercept.*|RATIO.*`,  regex=TRUE)
polis.draw


## ----summariseModel1c1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis.draw |> median_hdci()


## ----summariseModel1c3, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
polis.gather <- polis.rstanarm3 |> gather_draws(`(Intercept)`,RATIO) |>
  median_hdci()


## ----summariseModel1c4, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
polis.rstanarm3 |>
  gather_draws(`(Intercept)`, RATIO) |>
  ggplot() +
  stat_halfeye(aes(x=.value,  y=.variable)) +
  facet_wrap(~.variable, scales='free')


## ----summariseModel1c5, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
polis.rstanarm3 |>
  gather_draws(`(Intercept)`, RATIO) |>
  group_by(.variable) |>
  mutate(.value=exp(.value)) |>
  median_hdci()


## ----summariseModel1j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis.rstanarm3 |> plot(plotfun='mcmc_intervals')


## ----summariseModel1d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis.rstanarm3 |> tidy_draws()


## ----summariseModel1e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis.rstanarm3 |> spread_draws(`(Intercept)`, RATIO)
# OR via regex
polis.rstanarm3 |> spread_draws(`.Intercept.*|RATIO.*`,  regex=TRUE)

## summarised
polis.rstanarm3 |>
    spread_draws(`(Intercept)`, RATIO) |>
    summarise_draws("median",
                    ~ HDInterval::hdi(.x),
                    "rhat",
                    "ess_bulk")

## summarised on fractional scale
polis.rstanarm3 |>
    spread_draws(`(Intercept)`, RATIO) |>
    mutate(across(everything(), exp)) |>
    summarise_draws("median",
                    ~ HDInterval::hdi(.x),
                    "rhat",
                    "ess_bulk")


## ----summariseModel1f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis.rstanarm3 |> posterior_samples() |> as_tibble()


## ----summariseModel1g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis.rstanarm3 |> bayes_R2() |> median_hdci()


## ----summariseModel2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis_brm3 |> summary()


## ----summariseModel2a1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=FALSE----
polis.sum <- polis_brm3 |> summary()


## ----summariseModel2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis_brm3$fit |>
    tidyMCMC(
      estimate.method = "median", conf.int = TRUE,
      conf.method = "HPDinterval", rhat = TRUE, ess = TRUE
    )

## ----summariseModel2b1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
polis.tidy <- polis_brm3$fit |> tidyMCMC(estimate.method = 'median',  conf.int = TRUE,  conf.method = 'HPDinterval',  rhat = TRUE, ess = TRUE)


## ----summariseModel2dd, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis_brm3 |>
    summarise_draws(median,
                    HDInterval::hdi,
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel2d2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis_brm3 |>
    summarise_draws(median,
                    ~HDInterval::hdi(.x, credMass = 0.9),
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel2m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis_brm3 |> as_draws_df()
## summarised
polis_brm3 |>
  as_draws_df() |>
  summarise_draws(
    "median",
    ~ HDInterval::hdi(.x),
    "rhat",
    "ess_bulk"
  )

## summarised on fractional scale
polis_brm3 |>
    as_draws_df() |>
    dplyr::select(starts_with("b_")) |>
    mutate(across(everything(), exp)) |>
    summarise_draws("median",
                    ~ HDInterval::hdi(.x),
                    "rhat",
                    "ess_bulk")


## ----summariseModel2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis.draw <- polis_brm3 |> gather_draws(b_Intercept, b_RATIO)
## OR via regex
polis.draw <- polis_brm3 |> gather_draws(`b_.*`,  regex=TRUE)
polis.draw


## ----summariseModel2c1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis.draw |> median_hdci()


## ----summariseModel2c3, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
polis.gather <- polis_brm3 |>
    gather_draws(b_Intercept, b_RATIO) |>
    median_hdci()


## ----summariseModel2c4, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
polis_brm3 |>
  gather_draws(b_Intercept, b_RATIO) |>
  ggplot() +
  stat_halfeye(aes(x = .value,  y = .variable)) +
  facet_wrap(~.variable, scales = 'free')

polis.draw |>
    ggplot() +
    stat_halfeye(aes(x = .value,  y = .variable,
                     fill = stat(ggdist::cut_cdf_qi(cdf,
                               .width = c(0.5, 0.8, 0.95),
                               labels = scales::percent_format())))) +
    scale_fill_brewer('Interval', direction = -1, na.translate = FALSE) +
    facet_wrap(~.variable, scales = 'free') +
    theme_bw()
## polis.draw |>
##     ggplot() +
##     stat_halfeye(aes(x = .value,  y = .variable,
##                      fill = stat(ggdist::cut_cdf_qi(cdf,
##                                .width = c(0.5, 0.8, 0.95),
##                                labels = scales::percent_format())))) +
##     scale_fill_brewer('Interval', direction = -1, na.translate = FALSE) +
##     theme_bw()


## ----summariseModel2c5, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
polis_brm3 |>
  gather_draws(b_Intercept, b_RATIO) |>
  group_by(.variable) |>
  mutate(.value = exp(.value)) |>
  median_hdci()


## ----summariseModel2j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis_brm3 |> mcmc_plot(type = 'intervals')


## ----summariseModel2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis_brm3 |> tidy_draws()


## ----summariseModel2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis_brm3 |> spread_draws(b_Intercept, b_RATIO)
# OR via regex
polis_brm3 |> spread_draws(`b_.*`,  regex=TRUE)

## summarised
polis_brm3 |>
    as_draws_df() |>
    dplyr::select(starts_with("b_")) |>
    summarise_draws("median",
                    ~ HDInterval::hdi(.x),
                    "rhat",
                    "ess_bulk")
## summarised on fractional scale
polis_brm3 |>
    as_draws_df() |>
    dplyr::select(starts_with("b_")) |>
    mutate(across(everything(), exp)) |>
    summarise_draws("median",
                    ~ HDInterval::hdi(.x),
                    "rhat",
                    "ess_bulk")


polis_brm3 |>
    tidy_draws() |>
    exp() |>
    summarise_draws(median,HDInterval::hdi, rhat, ess_bulk, ess_tail) |>
    filter(variable %in% c('b_Intercept', 'b_RATIO'))


polis_brm3 |>
    tidy_draws() |>
    exp() |>
    dplyr::select(starts_with("b_")) |>
    summarise_draws(median,HDInterval::hdi, rhat, ess_bulk, ess_tail)


## ----summariseModel2f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis_brm3 |> posterior_samples() |> as_tibble()


## ----summariseModel2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis_brm3 |>
    bayes_R2()
## OR as median and hdci
polis_brm3 |>
    bayes_R2(summary = FALSE) |>
    median_hdci()


## -----------------------------------------------------------------------------
#| label: modelsummary
#| results: markup
#| eval: true
#| echo: true
#| cache: false
polis_brm3 |> modelsummary(
  statistic = c("conf.low", "conf.high"),
  shape = term ~ statistic
)

polis_brm3 |> modelsummary(
  statistic = c("conf.low", "conf.high"),
  shape = term ~ statistic,
  exponentiate = TRUE
)

modelsummary(list("Raw" = polis_brm3, "Exponentiated" = polis_brm3),
  statistic = c("conf.low", "conf.high"),
  shape = term ~ statistic,
  exponentiate = c(FALSE, TRUE)
)


## -----------------------------------------------------------------------------
#| label: modelsummary_plot
#| results: markup
#| eval: true
#| echo: true
#| cache: false
polis_brm3 |> modelplot()
polis_brm3 |> modelplot(exponentiate = TRUE)


## ----LD501a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis.rstanarm3 |> tidy_draws() |>
  mutate(LD50 = -1*`(Intercept)`/RATIO) |>
  pull(LD50) |>
  median_hdci()


## ----LD502a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
polis_brm3 |> tidy_draws() |>
  mutate(ED50 = -1*b_Intercept/b_RATIO) |>
  pull(ED50) |>
  median_hdci()


## ----figureModel1a, results='markdown', eval=TRUE, mhidden=TRUE---------------
## Using emmeans
polis.grid = with(polis, list(RATIO = seq(min(RATIO), max(RATIO), len=100)))

newdata = emmeans(polis.rstanarm3, ~RATIO, at=polis.grid, type='response') |> as.data.frame()
head(newdata)

ggplot(newdata, aes(y=prob, x=RATIO)) +
    geom_point(data=polis, aes(y=PA)) +
    geom_line() +
    geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD), fill='blue', alpha=0.3) +
    scale_y_continuous('PA') +
    scale_x_continuous('RATIO') +
    theme_classic()

## spaghetti plot
newdata = emmeans(polis.rstanarm3, ~RATIO, at=polis.grid) |>
    gather_emmeans_draws() |>
    mutate(.value = plogis(.value))
newdata |> head()
ggplot(newdata,  aes(y=.value,  x=RATIO)) +
  geom_line(aes(group=.draw),  alpha=0.01) +
  geom_point(data=polis,  aes(y=PA))


## ----figureModel1b, results='markdown', eval=TRUE, mhidden=TRUE---------------
## Using emmeans
polis.grid = with(polis, list(RATIO = seq(min(RATIO), max(RATIO), len=100)))

newdata = emmeans(polis.rstanarm3, ~RATIO, at=polis.grid, type='response') |> as.data.frame()
head(newdata)

ggplot(newdata, aes(y=prob, x=RATIO)) +
    geom_point(data=polis, aes(y=PA)) +
    geom_line() +
    geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD), fill='blue', alpha=0.3) +
    scale_y_continuous(expression(Presence/absence~of~italic(Uta)~lizards)) +
    scale_x_continuous(expression(Island~perimeter:Area~ratio)) +
    theme_classic()

## spaghetti plot
newdata = emmeans(polis.rstanarm3, ~RATIO, at=polis.grid) |>
    gather_emmeans_draws() |>
    mutate(.value = plogis(.value))
newdata |> head()
ggplot(newdata,  aes(y=.value,  x=RATIO)) +
    geom_line(aes(group=.draw),  alpha=0.01) +
    geom_point(data=polis,  aes(y=PA)) +
    scale_y_continuous(expression(Presence/absence~of~italic(Uta)~lizards)) +
    scale_x_continuous(expression(Island~perimeter:Area~ratio)) +
    theme_classic()


## ----figureModel2a, results='markdown', eval=TRUE, mhidden=TRUE---------------
## Using emmeans
polis_grid <- with(polis, list(RATIO = modelr::seq_range(RATIO, n = 100)))

newdata <- polis_brm3 |>
    emmeans(~RATIO, at = polis_grid, type = 'response') |>
    as.data.frame()
head(newdata)

## Using raw data for points
newdata |>
    ggplot(aes(y = response, x = RATIO)) +
    geom_point(data = polis, aes(y = PA)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD), fill = 'blue', alpha = 0.3) +
    scale_y_continuous('PA') +
    scale_x_continuous('RATIO') +
    theme_classic()


## Using partial residuals for points

partial.obs <- polis |>
    bind_cols(Pred = predict(polis_brm3)[,'Estimate'],
              Resid = residuals(polis_brm3)[,'Estimate']) |>
    mutate(
        Obs = round(Pred + Resid, 0)
    )

newdata |>
    ggplot(aes(y = response, x = RATIO)) +
    geom_point(data = partial.obs, aes(y = Obs)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD), fill = 'blue', alpha = 0.3) +
    scale_y_continuous('PA') +
    scale_x_continuous('RATIO') +
    theme_classic()

## spaghetti plot
newdata = emmeans(polis_brm3, ~RATIO, at=polis_grid) |>
    gather_emmeans_draws() |>
    mutate(.value = plogis(.value))
newdata |> head()
ggplot(newdata,  aes(y=.value,  x=RATIO)) +
  geom_line(aes(group=.draw),  alpha=0.01) +
  geom_point(data=polis,  aes(y=PA))



## ----figureModel2b, results='markdown', eval=TRUE, mhidden=TRUE---------------
## Using emmeans
polis_grid <- with(polis, list(RATIO = seq(min(RATIO), max(RATIO), len=100)))

newdata <- polis_brm3 |>
  emmeans(~RATIO, at=polis_grid, type='response') |>
  as.data.frame()
head(newdata)

ggplot(newdata, aes(y=response, x=RATIO)) +
    geom_point(data=polis, aes(y=PA)) +
    geom_line() +
    geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD), fill='blue', alpha=0.3) +
    scale_y_continuous(expression(Presence/absence~of~italic(Uta)~lizards)) +
    scale_x_continuous(expression(Island~perimeter:Area~ratio)) +
    theme_classic()

## spaghetti plot
newdata <- polis_brm3 |>
  emmeans(~RATIO, at = polis_grid) |>
    gather_emmeans_draws() |>
    mutate(.value = plogis(.value))
newdata |> head()
ggplot(newdata,  aes(y=.value,  x=RATIO)) +
    geom_line(aes(group=.draw),  alpha=0.01) +
    geom_point(data=polis,  aes(y=PA)) +
    scale_y_continuous(expression(Presence/absence~of~italic(Uta)~lizards)) +
    scale_x_continuous(expression(Island~perimeter:Area~ratio)) +
    theme_classic()

