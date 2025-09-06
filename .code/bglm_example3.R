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


## ----readData, results='markdown', eval=TRUE----------------------------------
peake <- read_csv("../data/peakquinn.csv", trim_ws = TRUE)


## -----------------------------------------------------------------------------
#| label: examinData
glimpse(peake)


## -----------------------------------------------------------------------------
## Explore the first 6 rows of the data
head(peake)


## -----------------------------------------------------------------------------
str(peake)


## -----------------------------------------------------------------------------
peake |> datawizard::data_codebook()


## -----------------------------------------------------------------------------
peake |> modelsummary::datasummary_skim()


## ----EDA, results='markdown', eval=TRUE, mhidded=TRUE, warning=FALSE, message=FALSE----
ggplot(peake, aes(y=INDIV, x=AREA)) +
  geom_point()+
  geom_smooth()


## ----EDA1, results='markdown', eval=TRUE, mhidded=TRUE, warning=FALSE, message=FALSE----
ggplot(peake, aes(y=INDIV)) + geom_boxplot()

ggplot(peake, aes(y=AREA)) + geom_boxplot()


## ----EDA2, results='markdown', eval=TRUE, mhidded=TRUE, warning=FALSE, message=FALSE----
ggplot(peake, aes(y=INDIV, x=AREA)) +
  geom_point()+
  geom_smooth() +
  scale_y_log10() +
  scale_x_log10()


## ----lm, results='markdown', eval=TRUE, mhidden=TRUE--------------------------
summary(glm(INDIV ~ log(AREA), data=peake, family=poisson()))
summary(MASS::glm.nb(INDIV ~ log(AREA), data=peake))


## ----fitModel1a, results='markdown', eval=TRUE, mhidden=TRUE------------------
peake_rstanarm = stan_glm(INDIV ~ log(AREA), data=peake,
                          family=poisson(),
                          iter = 5000, warmup = 1000,
                          chains = 3, thin = 5, refresh = 0)


## ----fitModel1b, results='markdown', eval=TRUE, mhidden=TRUE------------------
prior_summary(peake_rstanarm)


## ----fitModel1d, results='markdown', eval=TRUE, mhidden=TRUE------------------
2.5/sd(log(peake$AREA))


## ----fitModel1f, results='markdown', eval=TRUE, mhidden=TRUE------------------
peake_rstanarm1 <- update(peake_rstanarm,  prior_PD=TRUE)

## ----fitModel1g, results='markdown', eval=TRUE, mhidden=TRUE------------------
ggemmeans(peake_rstanarm1,  ~AREA) |> plot(show_data=TRUE)
ggemmeans(peake_rstanarm1,  ~AREA) |> plot(jitter = FALSE, show_data = TRUE) + scale_y_log10()


## ----fitModel1h, results='markdown', eval=TRUE, mhidden=TRUE------------------
peake_rstanarm2 <- stan_glm(INDIV ~ log(AREA), data=peake,
                          family=poisson(),
                          prior_intercept = normal(6, 2.8, autoscale=FALSE),
                          prior = normal(0, 2.3, autoscale=FALSE),
                          prior_PD=TRUE,
                          iter = 5000, warmup = 1000,
                          chains = 3, thin = 5, refresh = 0
                          )


## ----fitModel1i, results='markdown', eval=TRUE, mhidden=TRUE------------------
ggemmeans(peake_rstanarm2,  ~AREA) |>
    plot(show_data=TRUE) +
    scale_y_log10()


## ----fitModel1j, results='markdown', eval=TRUE, mhidden=TRUE------------------
peake_rstanarm3 <- update(peake_rstanarm2,  prior_PD=FALSE)


## ----modelFit1k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
posterior_vs_prior(peake_rstanarm3, color_by='vs', group_by=TRUE,
                   facet_args=list(scales='free_y'))


## ----modelFit1l, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggemmeans(peake_rstanarm3,  ~AREA) |> plot(show_data=TRUE)


## ----fitModel2a, results='markdown', eval=TRUE, mhidden=TRUE------------------
peake_form <- bf(INDIV ~ scale(log(AREA)), family=poisson())
peake_brm <- brm(peake_form,
                data=peake,
                iter = 5000, warmup = 1000,
                chains = 3, cores = 3,
                thin = 5, refresh = 0,
                backend = "cmdstanr")


## ----fitModel2b, results='markdown', eval=TRUE, mhidden=TRUE, paged.print=FALSE,tidy.opts = list(width.cutoff = 80), echo=2----
options(width=100)
prior_summary(peake_brm)
options(width=80)


## ----fitModel2c, results='markdown', eval=TRUE, mhidden=TRUE, paged.print=FALSE,tidy.opts = list(width.cutoff = 80), echo=2----
median(log(peake$INDIV))
mad(log(peake$INDIV))


## ----fitModel2d, results='markdown', eval=TRUE, mhidden=TRUE------------------
peake_form <- bf(INDIV ~ scale(log(AREA)), family=poisson())
peake_brm1 <- brm(peake_form,
                 data=peake,
                 prior=c(
                   prior(normal(0, 2.8), class='b')),
                 sample_prior = 'only',
                 iter = 5000, warmup = 1000,
                 chains = 3, cores = 3,
                 thin = 5, refresh = 0,
                 backend = "cmdstanr")


## ----fitModel2e, results='markdown', eval=TRUE, mhidden=TRUE------------------
peake_brm1 |>
  conditional_effects() |>
  plot(points=TRUE) |>
  _[[1]] + scale_y_log10()
#OR
plot(conditional_effects(peake_brm1), points=TRUE, plot = FALSE)[[1]] + scale_y_log10()

ggemmeans(peake_brm1,  ~AREA) |> plot(show_data=TRUE)
ggemmeans(peake_brm1,  ~AREA) |> plot(show_data=TRUE) + scale_y_log10()


## ----fitModel2h, results='markdown', eval=TRUE, mhidden=TRUE------------------
peake_form <- bf(INDIV ~ scale(log(AREA)), family=poisson())
priors <- prior(normal(6, 1.1), class = "Intercept") +
  prior(normal(0, 1.1), class = "b")
peake_brm2 <- brm(peake_form,
                 data=peake,
                 prior=priors,
                 sample_prior = 'only',
                 iter = 5000, warmup = 1000,
                 chains = 3, cores = 3,
                 thin = 5, refresh = 0,
                 backend = "cmdstanr")


## ----fitModel2i, results='markdown', eval=TRUE, mhidden=TRUE------------------
peake_brm2 |>
  conditional_effects() |>
  plot(points=TRUE) |>
  _[[1]] + scale_y_log10()
ggemmeans(peake_brm2,  ~AREA) |>
    plot(show_data=TRUE) +
    scale_y_log10()


## ----fitModel2j, results='markdown', eval=TRUE, mhidden=TRUE------------------
peake_brm3 <- update(peake_brm2,  sample_prior=TRUE, refresh=0)


## ----fitModel2j2, results='markdown', eval=TRUE, echo = FALSE, mhidden=TRUE----
save(peake_brm3, file = '../ws/testing/peake_brm3.RData')


## ----fitModel2k, results='markdown', eval=TRUE, mhidden=TRUE------------------
peake_brm3 |> get_variables()
peake_brm3 |>
  hypothesis("scalelogAREA=0") |>
  plot()
peake_brm3 |> SUYR_prior_and_posterior()


## ----fitModel2l, results='markdown', eval=TRUE, mhidden=TRUE------------------
standata(peake_brm3)
stancode(peake_brm3)


## ----modelValidation1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_mcmc()


## ----modelValidation1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(peake_rstanarm3, plotfun='mcmc_trace')


## ----modelValidation1c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(peake_rstanarm3, 'acf_bar')


## ----modelValidation1d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(peake_rstanarm3, 'rhat_hist')


## ----modelValidation1e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(peake_rstanarm3, 'neff_hist')


## ----Validation1f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(peake_rstanarm3, 'combo')
plot(peake_rstanarm3, 'violin')


## ----modelValidation1g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_trace(peake_rstanarm3)


## ----modelValidation1h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_ac(peake_rstanarm3)


## ----modelValidation1i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_rhat(peake_rstanarm3)


## ----modelValidation1j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_ess(peake_rstanarm3)


## ----modelValidation1k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_dens(peake_rstanarm3, separate_chains = TRUE)


## ----modelValidation1l, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
peake_ggs <- ggs(peake_rstanarm3)
ggs_traceplot(peake_ggs)


## ----modelValidation1m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_autocorrelation(peake_ggs)


## ----modelValidation1n, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_Rhat(peake_ggs)


## ----modelValidation1o, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_effective(peake_ggs)


## ----modelValidation1p, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_crosscorrelation(peake_ggs)


## ----modelValidation1q, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_grb(peake_ggs)


## ----modelValidation2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_mcmc()


## ----modelValidation2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
mcmc_plot(peake_brm3, type='trace')


## ----modelValidation2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
mcmc_plot(peake_brm3, type='acf_bar')


## ----modelValidation2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
mcmc_plot(peake_brm3, type='rhat_hist')


## ----modelValidation2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
mcmc_plot(peake_brm3, type='neff_hist')


## ----modelValidation2f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
mcmc_plot(peake_brm3, type='combo')
mcmc_plot(peake_brm3, type='violin')


## ----modelValidation2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_trace(peake_brm3$fit)


## ----modelValidation2h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_ac(peake_brm3$fit)


## ----modelValidation2i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_rhat(peake_brm3$fit)


## ----modelValidation2j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_ess(peake_brm3$fit)


## ----modelValidation2k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_dens(peake_brm3$fit, separate_chains = TRUE)


## ----modelValidation2l, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=7----
peake_ggs <- ggs(peake_brm3)
ggs_traceplot(peake_ggs)


## ----modelValidation2m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=7----
ggs_autocorrelation(peake_ggs)


## ----modelValidation2n, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_Rhat(peake_ggs)


## ----modelValidation2o, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_effective(peake_ggs)


## ----modelValidation2p, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_crosscorrelation(peake_ggs)


## ----modelValidation2q, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_grb(peake_ggs)


## ----modelValidation3a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation3b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(peake_rstanarm3,  plotfun='dens_overlay')


## ----modelValidation3c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#pp_check(peake_rstanarm3, plotfun='error_scatter_avg')


## ----modelValidation3d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(peake_rstanarm3, x=peake$AREA, plotfun='error_scatter_avg_vs_x')


## ----modelValidation3e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(peake_rstanarm3, x=peake$AREA, plotfun='intervals')


## ----modelValidation3f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(peake_rstanarm3, x=peake$AREA, plotfun='ribbon')


## ----modelValidation3g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(peake_rstanarm3)


## ----modelValidation4a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
preds <- posterior_predict(peake_rstanarm3,  ndraws=250,  summary=FALSE)
peake_resids <- createDHARMa(simulatedResponse = t(preds),
                            observedResponse = peake$INDIV,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse = TRUE)
plot(peake_resids)
peake_resids |> testDispersion()


## ----modelValidation5a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation5b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(peake_brm3,  type='dens_overlay', ndraws = 100)


## ----modelValidation5c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#pp_check(peake_brm3, type='error_scatter_avg')


## ----modelValidation5d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(peake_brm3, x='AREA', type='error_scatter_avg_vs_x')


## ----modelValidation5e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(peake_brm3, x='AREA', type='intervals')


## ----modelValidation5f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(peake_brm3, x='AREA', type='ribbon')


## ----modelValidation5g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(peake_brm3)


## ----modelValidation6aa, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=10----
peake_resids <- make_brms_dharma_res(peake_brm3, integerResponse = FALSE)
wrap_elements(~testUniformity(peake_resids)) +
               wrap_elements(~plotResiduals(peake_resids, form = factor(rep(1, nrow(peake))))) +
               wrap_elements(~plotResiduals(peake_resids, quantreg = FALSE)) +
               wrap_elements(~testDispersion(peake_resids))



## ----modelValidation6a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
preds <- posterior_predict(peake_brm3,  ndraws=250,  summary=FALSE)
peake_resids <- createDHARMa(simulatedResponse = t(preds),
                            observedResponse = peake$INDIV,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse = TRUE)
plot(peake_resids)
peake_resids |> testDispersion()


## ----fitNBModel1a, results='markdown', eval=TRUE, mhidden=TRUE----------------
peake_rstanarm4 <- stan_glm(INDIV ~ log(AREA), data=peake,
                          family=neg_binomial_2(),
                          prior_intercept = normal(6, 2.8, autoscale=FALSE),
                          prior = normal(0, 2.3, autoscale=FALSE),
                          prior_aux = rstanarm::exponential(rate=1, autoscale=FALSE),
                          iter = 5000, warmup = 1000,
                          chains = 3, thin = 5, refresh = 0
                          )


## ----fitNBModel1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
posterior_vs_prior(peake_rstanarm4, color_by='vs', group_by=TRUE,
                   facet_args=list(scales='free_y'))


## ----fitNBModel1c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## ggemmeans(peake_rstanarm4,  ~AREA) |> plot(show_data=TRUE)
ggpredict(peake_rstanarm4,  ~AREA) |> plot(show_data=TRUE)


## ----fitNBModel1d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(peake_rstanarm4, plotfun='mcmc_trace')
plot(peake_rstanarm4, 'acf_bar')
plot(peake_rstanarm4, 'rhat_hist')
plot(peake_rstanarm4, 'neff_hist')
pp_check(peake_rstanarm4, x=peake$AREA, plotfun='dens_overlay')
pp_check(peake_rstanarm4, x=peake$AREA, plotfun='error_scatter_avg_vs_x')
pp_check(peake_rstanarm4, x=peake$AREA, plotfun='intervals')


## ----fitNBModel1e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
preds <- posterior_predict(peake_rstanarm4,  ndraws=250,  summary=FALSE)
peake_resids <- createDHARMa(simulatedResponse = t(preds),
                            observedResponse = peake$INDIV,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse = TRUE)
plot(peake_resids)


## ----fitNBModel1f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
(peake_rstanarm3.loo <- loo(peake_rstanarm3))


## ----fitNBModel1g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
(peake_rstanarm4.loo <- loo(peake_rstanarm4))


## ----fitNBModel1h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
loo_compare(peake_rstanarm3.loo, peake_rstanarm4.loo)


## ----fitNBModel2a2, results='markdown', eval=TRUE, mhidden=TRUE---------------
standist::visualize("gamma(0.01, 0.01)", "gamma(2, 1)", "inv_gamma(0.4, 0.3)")


## ----fitNBModel2a, results='markdown', eval=TRUE, mhidden=TRUE----------------
peake_form <- bf(INDIV ~ scale(log(AREA)), family=negbinomial())
get_prior(peake_form, data = peake)
priors <-
  prior(normal(6, 1.1),  class='Intercept') +
  prior(normal(0, 1.1), class='b') +
  prior(gamma(0.01, 0.01), class='shape')
  ## prior(gamma(2, 1), class='shape')
peake_brm4 <- brm(peake_form,
                 data=peake,
                 prior = priors,
                 sample_prior=TRUE,
                 iter = 5000, warmup = 1000,
                 chains = 3, cores = 3,
                 thin = 5, refresh = 0,
                 backend = "cmdstan")



## ----fitNBModel2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
peake_brm4 |> SUYR_prior_and_posterior()


## ----fitNBModel2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
peake_brm4 |>
  conditional_effects() |>
  plot(points=TRUE) |>
  _[[1]] + scale_y_log10()
peake_brm4 |>
  conditional_effects() |>
  plot(points=TRUE) |>
  _[[1]] + scale_y_log10() + scale_x_log10()
## ggemmeans(peake_brm4,  ~AREA) |> plot(show_data=TRUE)
ggpredict(peake_brm4,  ~AREA) |> plot(show_data=TRUE)


## ----fitNBModel2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
peake_brm4$fit |> stan_trace()
peake_brm4$fit |> stan_ac()
peake_brm4$fit |> stan_rhat()
peake_brm4$fit |> stan_ess()
pp_check(peake_brm4, x='AREA',type='dens_overlay', ndraws = 100)
pp_check(peake_brm4, x='AREA', type='intervals')


## ----fitNBModel2e2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=10----
peake_resids <- make_brms_dharma_res(peake_brm4, integerResponse = FALSE)
wrap_elements(~testUniformity(peake_resids)) +
               wrap_elements(~plotResiduals(peake_resids, form = factor(rep(1, nrow(peake))))) +
               wrap_elements(~plotResiduals(peake_resids, quantreg = FALSE)) +
               wrap_elements(~testDispersion(peake_resids))



## ----fitNBModel2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
preds <- posterior_predict(peake_brm4,  ndraws=250,  summary=FALSE)
peake_resids <- createDHARMa(simulatedResponse = t(preds),
                            observedResponse = peake$INDIV,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse = TRUE)
plot(peake_resids)


## ----fitNBModel2f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
(peake_brm3.loo <- loo(peake_brm3))


## ----fitNBModel2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
(peake_brm4.loo <- loo(peake_brm4))


## ----fitNBModel2h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
loo_compare(peake_brm3.loo, peake_brm4.loo)


## ----partialPlot1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_rstanarm4 |> ggpredict() |> plot(show_data=TRUE)


## ----partialPlot1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
## Note, does not seem to backtransform...
peake_rstanarm4 |> ggemmeans(~AREA, type = "fixed", back.transform = TRUE) |> plot()


## ----partialPlot1c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_rstanarm4 |> epred_draws(newdata=peake) |>
  median_hdci() |>
  ggplot(aes(x=AREA, y=.epred)) +
  geom_ribbon(aes(ymin=.lower, ymax=.upper), fill='blue', alpha=0.3) +
  geom_line() +
  geom_point(data=peake,  aes(y=INDIV,  x=AREA))


## ----partialPlot2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_brm4 |> conditional_effects() |> plot(points = TRUE)
peake_brm4 |> conditional_effects(spaghetti=TRUE,ndraws=200) |> plot(points = TRUE) + theme_classic()
ce <- peake_brm4 |> conditional_effects(spaghetti=TRUE,ndraws=200)
plot(ce, points = TRUE, plot = FALSE)[[1]] + theme_classic()


## ----partialPlot2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_brm4 |> ggpredict() |> plot(show_data=TRUE)


## ----partialPlot2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_brm4 |> ggemmeans(~AREA) |> plot(show_data=TRUE)


## ----partialPlot2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_brm4 |> epred_draws(newdata=peake) |>
  median_hdci() |>
  ggplot(aes(x=AREA, y=.epred)) +
  geom_ribbon(aes(ymin=.lower, ymax=.upper), fill='blue', alpha=0.3) +
  geom_line() +
  geom_point(data=peake,  aes(y=INDIV,  x=AREA))


## ----summariseModel1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
summary(peake_rstanarm4)


## ----summariseModel1a1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=FALSE----
peake_sum <- summary(peake_rstanarm4)


## ----summariseModel1dd, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_rstanarm4$stanfit |>
    summarise_draws(median,
                    HDInterval::hdi,
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel1d2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_rstanarm4$stanfit |>
    summarise_draws(median,
                    ~HDInterval::hdi(.x, credMass = 0.9),
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel1d3, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_rstanarm4$stanfit |>
    ## tidy_draws() |>
    ## exp() |>
    summarise_draws(
        ~ median(exp(.x)),
        ~HDInterval::hdi(exp(.x)),
        rhat, length, ess_bulk, ess_tail)


## ----summariseModel1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tidyMCMC(peake_rstanarm4$stanfit, estimate.method='median',  conf.int=TRUE,  conf.method='HPDinterval',  rhat=TRUE, ess=TRUE)

## ----summariseModel1b1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
peake_tidy <- tidyMCMC(peake_rstanarm4$stanfit, estimate.method='median',  conf.int=TRUE,  conf.method='HPDinterval',  rhat=TRUE, ess=TRUE)


## ----summariseModel1m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_rstanarm4$stanfit |> as_draws_df()

## summarised
peake_rstanarm4$stanfit |>
    as_draws_df() |>
    exp() |>
    summarise_draws(median,
                    HDInterval::hdi,
                    rhat,
                    ess_bulk, ess_tail)
## summarised on fractional scale
peake_rstanarm4$stanfit |>
    as_draws_df() |>
    dplyr::select(matches('Intercept|AREA')) |>
    summarise_draws(median,
                    HDInterval::hdi,
                    rhat,
                    length,
                    ess_bulk, ess_tail)
## OR
names <- peake_rstanarm4 |>
    formula() |>
    model.matrix(peake) |>
    colnames()

peake_rstanarm4$stanfit |>
    as_draws_df() |>
    dplyr::select(any_of(names)) |>
    mutate(across(everything(), exp)) |>
    summarise_draws(median,
                    HDInterval::hdi,
                    rhat,
                    length,
                    ess_bulk, ess_tail)


## ----summariseModel1c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_rstanarm4 |> get_variables()
peake_draw <- peake_rstanarm4 |> gather_draws(`.Intercept.*|.*AREA.*`,  regex=TRUE)
peake_draw


## ----summariseModel1c1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_draw |> median_hdci()


## ----summariseModel1c3, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
peake_gather <- peake_rstanarm4 |> gather_draws(`.Intercept.*|.*AREA.*`,  regex=TRUE) |>
  median_hdci()


## ----summariseModel1c4, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
peake_rstanarm4 |>
  gather_draws(`.Intercept.*|.*AREA.*`, regex=TRUE) |>
  ggplot() +
  stat_halfeye(aes(x=.value,  y=.variable)) +
  facet_wrap(~.variable, scales='free')


## ----summariseModel1c5, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
peake_rstanarm4 |>
  gather_draws(`.Intercept.*|.*AREA.*`, regex=TRUE) |>
  group_by(.variable) |>
  mutate(.value=exp(.value)) |>
  median_hdci()


## ----summariseModel1j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_rstanarm4 |> plot(plotfun='mcmc_intervals')


## ----summariseModel1d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_rstanarm4 |> tidy_draws()


## ----summariseModel1e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_rstanarm4 |> spread_draws(`.Intercept.*|.*AREA.*`,  regex=TRUE)

## summarised
peake_rstanarm4 |>
    spread_draws(`.Intercept.*|.*AREA.*`,  regex=TRUE) |>
    summarise_draws("median",
                    ~ HDInterval::hdi(.x),
                    "rhat",
                    "ess_bulk")

## summarised on fractional scale
peake_rstanarm4 |>
    spread_draws(`.Intercept.*|.*AREA.*`,  regex=TRUE) |>
    mutate(across(everything(), exp)) |>
    summarise_draws("median",
                    ~ HDInterval::hdi(.x),
                    "rhat",
                    "ess_bulk")


## ----summariseModel1f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_rstanarm4 |> posterior_samples() |> as_tibble()


## ----summariseModel1g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
## peake_rstanarm4 |> bayes_R2() |> median_hdci


## ----summariseModel2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
summary(peake_brm4)


## ----summariseModel2a1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=FALSE----
peake_sum <- summary(peake_brm4)


## ----summariseModel2dd, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_brm4 |>
    summarise_draws(median,
                    HDInterval::hdi,
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel2d2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_brm4 |>
    summarise_draws(median,
                    ~HDInterval::hdi(.x, credMass = 0.9),
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tidyMCMC(peake_brm4$fit, estimate.method='median',  conf.int=TRUE,  conf.method='HPDinterval',  rhat=TRUE, ess=TRUE)

## ----summariseModel2b1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
peake_tidy <- tidyMCMC(peake_brm4$fit, estimate.method='median',  conf.int=TRUE,  conf.method='HPDinterval',  rhat=TRUE, ess=TRUE)


## ----summariseModel2m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_brm4 |> as_draws_df()
## summarised
peake_brm4 |>
  as_draws_df() |>
  summarise_draws(
    median,
    ~ HDInterval::hdi(.x),
    length,
    rhat,
    ess_bulk, ess_tail
  )

## summarised on fractional scale
peake_brm4 |>
    as_draws_df() |>
    dplyr::select(starts_with("b_")) |>
    mutate(across(everything(), exp)) |>
    ## exp() |>
    summarise_draws(median,
                    HDInterval::hdi,
                    rhat,
                    length,
                    ess_bulk, ess_tail)


## To get the slope on a scale that is based on a unit change in the
## predictor (rather than a unit change in scaled log-transformed
## area)
## Note, type= "response" is would be ignored
peake_brm4 |>
  emtrends(~1, var = "AREA")

## As if we had not scaled the predictor (but did log-transform it)
peake_brm4 |>
  emtrends(~0, var = "log(AREA)") |>
  tidy_draws() |>
  mutate(.value = exp(`1 overall`)) |>
  summarise_draws(median,
                  HDInterval::hdi)

## peake_brm4 |>
##   emmeans(~AREA, type = "response", at = list(AREA = c(1, 2))) |>
##   pairs(reverse = TRUE)

## peake_brm4 |>
##   emmeans(~AREA, type = "response", at = list(AREA = c(10, 20))) |>
##   pairs(reverse = TRUE)


## ----summariseModel2d1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_brm4 |> tidy_draws()

peake_brm4 |> get_variables()
peake_brm4$fit |>
    tidy_draws() |>
    dplyr::select(matches('^b_.*'))  |>
    exp() |>
    summarise_draws(median,
                    HDInterval::hdi,
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_brm4 |>
    tidy_draws() |>
    exp() |>
    gather_variables() |>
    median_hdci()


## ----summariseModel2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_brm4 |> get_variables()
peake_draw <- peake_brm4 |> gather_draws(`b_.*`,  regex=TRUE)
peake_draw


## ----summariseModel2c1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_draw |> median_hdci()


## ----summariseModel2c11, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_brm4 |>
    gather_draws(`b_.*`,  regex=TRUE) |>
    mutate(.value = exp(.value)) |>
    summarise_draws(median,
                    HDInterval::hdi,
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel2c3, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
peake_gather <- peake_brm4 |> gather_draws(`b_.*`,  regex=TRUE) |>
  median_hdci()


## ----summariseModel2c4, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
peake_brm4 |>
  gather_draws(`b_.*`, regex=TRUE) |>
  ggplot() +
  stat_halfeye(aes(x=.value,  y=.variable)) +
  facet_wrap(~.variable, scales='free')


## ----summariseModel2c5, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
peake_brm4 |>
  gather_draws(`.Intercept.*|.*AREA.*`, regex=TRUE) |>
  group_by(.variable) |>
  mutate(.value=exp(.value)) |>
  median_hdci()


## ----summariseModel2j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_brm4 |> plot(plotfun='mcmc_intervals')


## ----summariseModel2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_brm4 |> spread_draws(`b_.*`,  regex=TRUE)
## summarised
peake_brm4 |>
    as_draws_df() |>
    dplyr::select(starts_with("b_")) |>
    summarise_draws("median",
                    ~ HDInterval::hdi(.x),
                    "rhat",
                    "ess_bulk")
## summarised on fractional scale
peake_brm4 |>
    as_draws_df() |>
    dplyr::select(starts_with("b_")) |>
    mutate(across(everything(), exp)) |>
    summarise_draws("median",
                    ~ HDInterval::hdi(.x),
                    "rhat",
                    "ess_bulk")


## ----summariseModel2f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_brm4 |> posterior_samples() |> as_tibble()


## ----summariseModel2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_brm4 |> bayes_R2(summary=FALSE) |> median_hdci()


## -----------------------------------------------------------------------------
#| label: modelsummary
#| results: markup
#| eval: true
#| echo: true
#| cache: false
peake_brm4 |> modelsummary(
  statistic = c("conf.low", "conf.high"),
  shape = term ~ statistic
)

peake_brm4 |> modelsummary(
  statistic = c("conf.low", "conf.high"),
  shape = term ~ statistic,
  exponentiate = TRUE
)

modelsummary(list("Raw" = peake_brm4, "Exponentiated" = peake_brm4),
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
peake_brm4 |> modelplot()
peake_brm4 |> modelplot(exponentiate = TRUE)


## ----Probability1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_rstanarm4 |> as.data.frame() |> rename(lAREA=`log(AREA)`) |> hypothesis('lAREA>0')


## ----Probability1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,paged.print=FALSE----
peake_rstanarm4 |> tidy_draws() |> summarise(P=mean(`log(AREA)`>0))


## ----Probability1c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=1:5----
newdata = list(AREA=c(5000, 10000))
## fractional scale
peake_rstanarm4 |> emmeans(~AREA,  at=newdata, type = 'response') |> pairs(reverse = TRUE)
## absolute scale
peake_rstanarm4 |> emmeans(~AREA,  at=newdata) |> regrid() |> pairs(reverse = TRUE)
peake_mcmc <- peake_rstanarm4 |> emmeans(~AREA,  at=newdata) |> pairs(reverse = TRUE) |> as.data.frame()


## ----Probability1d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_mcmc <- peake_rstanarm4 |> emmeans(~AREA,  at=newdata) |>
    regrid() |>
    tidy_draws() |>
    rename_with(~str_replace(., 'AREA ', 'p')) |>
    mutate(Eff=p10000 - p5000,
           PEff=100*Eff/p5000)
peake_mcmc |> head()


## ----Probability1e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_mcmc |> tidyMCMC(estimate.method='median',
                       conf.int=TRUE, conf.method='HPDinterval')


## ----Probability1f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_mcmc |> median_hdci(PEff)


## ----Probability1g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, paged.print=FALSE----
peake_mcmc |> summarise(P=mean(PEff>50))


## ----Probability1h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_mcmc |> hypothesis('PEff>50')


## ----Probability1i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
## Note, the P value is on absolute difference
newdata = list(AREA=c(5000, 10000))
peake_rstanarm4 |>
    emmeans(~AREA,  at=newdata) |>
    regrid() |>
    pairs(reverse = TRUE) |>
    tidy_draws() |>
    summarise(across(contains('contrast'),
                     list(P = ~ mean(.>50),
                          HDCI = ~ median_hdci(.)),
                     .names = c('{.fn}')
                     )) |>
    tidyr::unpack(HDCI)


## ----Probability1j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
newdata = list(AREA=c(5000, 10000))
## Simple
peake_rstanarm4 |>
    emmeans(~AREA,  at=newdata) |>
    pairs(reverse = TRUE) |>
    regrid()

## More advanced (both P and percent change)
peake_mcmc <- peake_rstanarm4 |>
    emmeans(~AREA,  at=newdata) |>
    pairs(reverse = TRUE) |>
    regrid() |>
    tidy_draws() |>
    mutate(across(contains('contrast'), ~ 100*(. - 1)))

peake_mcmc |>
    summarise(across(contains('contrast'),
                     list(P = ~ mean(.>50),
                          HDCI = ~ median_hdci(.)),
                     .names = c('{.fn}')
                     )) |>
    tidyr::unpack(HDCI)



## ----Probability1k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_rstanarm4 |>
    linpred_draws(newdata = as.data.frame(newdata)) |>
    mutate(.linpred = exp(.linpred)) |>
    ungroup() |>
    group_by(.draw) |>
    summarise(Eff = diff(.linpred),
              PEff = 100*Eff/.linpred[1]) |>
    ungroup() |>
    mutate(P = mean(PEff>50)) |>
    pivot_longer(cols = -.draw) |>
    group_by(name) |>
    median_hdci()

##OR
peake_rstanarm4 |>
    epred_draws(newdata = as.data.frame(newdata)) |>
    ungroup() |>
    group_by(.draw) |>
    summarise(Eff = diff(.epred),
              PEff = 100*Eff/.epred[1]) |>
    ungroup() |>
    mutate(P = mean(PEff>50)) |>
    pivot_longer(cols = -.draw) |>
    group_by(name) |>
    median_hdci()

##OR for prediction of individual values
peake_rstanarm4 |>
    predicted_draws(newdata = as.data.frame(newdata)) |>
    ungroup() |>
    group_by(.draw) |>
    summarise(Eff = diff(.prediction),
              PEff = 100*Eff/.prediction[1]) |>
    ungroup() |>
    mutate(P = mean(PEff>50)) |>
    pivot_longer(cols = -.draw) |>
    group_by(name) |>
    median_hdci()


## ----Probability2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_brm4 |> hypothesis('scalelogAREA>0')


## ----Probability2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,paged.print=FALSE----
peake_brm4 |> tidy_draws() |> summarise(P=mean(b_scalelogAREA>0))


## ----Probability2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=1:5----
newdata = list(AREA=c(5000, 10000))
## fractional scale
peake_brm4 |> emmeans(~AREA,  at=newdata, type = 'response') |> pairs(reverse = TRUE)
## absolute scale
peake_brm4 |> emmeans(~AREA,  at=newdata) |> regrid() |> pairs(reverse = TRUE)
peake_mcmc <- peake_brm4 |> emmeans(~AREA,  at=newdata) |> pairs(reverse = TRUE) |> as.data.frame()


## ----Probability2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_mcmc <- peake_brm4 |> emmeans(~AREA,  at=newdata) |>
    regrid() |>
    tidy_draws() |>
    rename_with(~str_replace(., 'AREA ', 'p')) |>
    mutate(Eff=p10000 - p5000,
           PEff=100*Eff/p5000)
peake_mcmc |> head()


## ----Probability2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_mcmc |> tidyMCMC(estimate.method='median',
                       conf.int=TRUE, conf.method='HPDinterval')


## ----Probability2f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_mcmc |> median_hdci(PEff)


## ----Probability2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, paged.print=FALSE----
peake_mcmc |> summarise(P=mean(PEff>50))


## ----Probability2h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_mcmc |> hypothesis('PEff>50')


## ----Probability2i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
## Note, the P value is on absolute difference
newdata = list(AREA=c(5000, 10000))
peake_brm4 |>
    emmeans(~AREA,  at=newdata) |>
    regrid() |>
    pairs(reverse = TRUE) |>
    tidy_draws() |>
    summarise(across(contains('contrast'),
                     list(P = ~ mean(.>50),
                          HDCI = ~ median_hdci(.)),
                     .names = c('{.fn}')
                     )) |>
    tidyr::unpack(HDCI)


## ----Probability2j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
newdata = list(AREA=c(5000, 10000))
## Simple
peake_brm4 |>
    emmeans(~AREA,  at=newdata) |>
    pairs(reverse = TRUE) |>
    regrid()

## More advanced (both P and percent change)
peake_mcmc <- peake_brm4 |>
    emmeans(~AREA,  at=newdata) |>
    pairs(reverse = TRUE) |>
    regrid() |>
    tidy_draws() |>
    mutate(across(contains('contrast'), ~ 100*(. - 1)))

peake_mcmc |>
    summarise(across(contains('contrast'),
                     list(P = ~ mean(.>50),
                          HDCI = ~ median_hdci(.)),
                     .names = c('{.fn}')
                     )) |>
    tidyr::unpack(HDCI)



## ----Probability2k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
peake_brm4 |>
    linpred_draws(newdata = as.data.frame(newdata)) |>
    mutate(.linpred = exp(.linpred)) |>
    ungroup() |>
    group_by(.draw) |>
    summarise(Eff = diff(.linpred),
              PEff = 100*Eff/.linpred[1]) |>
    ungroup() |>
    mutate(P = mean(PEff>50)) |>
    pivot_longer(cols = -.draw) |>
    group_by(name) |>
    median_hdci()

##OR
peake_brm4 |>
    epred_draws(newdata = as.data.frame(newdata)) |>
    ungroup() |>
    group_by(.draw) |>
    summarise(Eff = diff(.epred),
              PEff = 100*Eff/.epred[1]) |>
    ungroup() |>
    mutate(P = mean(PEff>50)) |>
    pivot_longer(cols = -.draw) |>
    group_by(name) |>
    median_hdci()

##OR for prediction of individual values
peake_brm4 |>
    predicted_draws(newdata = as.data.frame(newdata)) |>
    ungroup() |>
    group_by(.draw) |>
    summarise(Eff = diff(.prediction),
              PEff = 100*Eff/.prediction[1]) |>
    ungroup() |>
    mutate(P = mean(PEff>50)) |>
    pivot_longer(cols = -.draw) |>
    group_by(name) |>
    median_hdci()


## ----figureModel1a, results='markdown', eval=TRUE, mhidden=TRUE---------------
## Using emmeans
peake_grid = with(peake, list(AREA = seq(min(AREA), max(AREA), len=100)))

newdata = emmeans(peake_rstanarm4, ~AREA, at=peake_grid, type='response') |> as.data.frame()
head(newdata)

ggplot(newdata, aes(y=prob, x=AREA)) +
    geom_point(data=peake, aes(y=INDIV)) +
    geom_line() +
    geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD), fill='blue', alpha=0.3) +
    scale_y_continuous('Individuals') +
    scale_x_continuous('Mussel clump area') +
    theme_classic()

## spaghetti plot
newdata = emmeans(peake_rstanarm4, ~AREA, at=peake_grid) |>
    regrid() |>
    gather_emmeans_draws()
newdata |> head()
ggplot(newdata,  aes(y=.value,  x=AREA)) +
    geom_line(aes(group=.draw), colour = 'orange', alpha=0.01) +
    geom_point(data=peake,  aes(y=INDIV)) +
    scale_y_continuous('Number of Individuals') +
    scale_x_continuous(expression(Mussel~clump~area~(per~mm^2))) +
    theme_classic()


## ----figureModel2a, results='markdown', eval=TRUE, mhidden=TRUE---------------
## Using emmeans
peake_grid = with(peake, list(AREA = seq(min(AREA), max(AREA), len=100)))

newdata = emmeans(peake_brm4, ~AREA, at=peake_grid, type='response') |> as.data.frame()
head(newdata)

ggplot(newdata, aes(y=prob, x=AREA)) +
    geom_point(data=peake, aes(y=INDIV)) +
    geom_line() +
    geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD), fill='blue', alpha=0.3) +
    scale_y_continuous('Individuals') +
    scale_x_continuous('Mussel clump area') +
    theme_classic()

## spaghetti plot
newdata = emmeans(peake_brm4, ~AREA, at=peake_grid) |>
    regrid() |>
    gather_emmeans_draws()
newdata |> head()
ggplot(newdata,  aes(y=.value,  x=AREA)) +
    geom_line(aes(group=.draw), colour = 'orange', alpha=0.01) +
    geom_point(data=peake,  aes(y=INDIV)) +
    scale_y_continuous('Number of Individuals') +
    scale_x_continuous(expression(Mussel~clump~area~(per~mm^2))) +
    theme_classic()

