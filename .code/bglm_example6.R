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
library(HDInterval)    #for HPD intervals
library(posterior)     #for posterior draws
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
library(ggridges)      #for ridge plots
library(patchwork)     #for multiple plots
library(modelsummary)  #for data and model summaries
theme_set(theme_grey()) #put the default ggplot theme back
source("helperFunctions.R")


## -----------------------------------------------------------------------------
quinn <- read_csv("../data/quinn.csv", trim_ws = TRUE)


## -----------------------------------------------------------------------------
#| label: examinData
glimpse(quinn)


## -----------------------------------------------------------------------------
## Explore the first 6 rows of the data
head(quinn)


## -----------------------------------------------------------------------------
str(quinn)


## -----------------------------------------------------------------------------
quinn |> datawizard::data_codebook()


## -----------------------------------------------------------------------------
quinn |> modelsummary::datasummary_skim()
quinn |>
    dplyr::select(-SQRTRECRUITS) |>
    modelsummary::datasummary_skim(
      type = "numeric",
      by = c("SEASON", "DENSITY")
    )


## ----dataprep, results='markdown', eval=TRUE----------------------------------
## A recent bug has infiltrated emmeans such that SEASON seems
## to be some sort of keyword.  So in order to prevent downstream issues
## when declaring season as a factor, we will also alter its name...
quinn <- quinn |>
  mutate(fSEASON = factor(SEASON,
                         levels = c("Spring", "Summer", "Autumn", "Winter")),
                         DENSITY = factor(DENSITY))


## ----EDA, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=6, mhidden=TRUE----
quinn |> head()
quinn |> ggplot(aes(y=RECRUITS, x=fSEASON, fill=DENSITY)) + geom_boxplot()


## ----fitModel1a, results='markdown', eval=TRUE, mhidden=TRUE------------------
quinn_rstanarmP <- stan_glm(RECRUITS~fSEASON*DENSITY,
                            data = quinn,
                            family = poisson(link = 'log'),
                            refresh = 0,
                            chains = 3, iter = 5000, thin = 5, warmup = 2000)


## ----fitModel1b, results='markdown', eval=TRUE, mhidden=TRUE------------------
quinn_rstanarmP |> prior_summary()


## ----fitModel1d, results='markdown', eval=TRUE, mhidden=TRUE------------------
2.5/apply(model.matrix(~fSEASON*DENSITY, quinn)[,-1], 2, sd)


## ----fitModel1f, results='markdown', eval=TRUE, mhidden=TRUE------------------
quinn_rstanarm1 <- stan_glm(RECRUITS~fSEASON*DENSITY, data = quinn,
                            family = poisson(link = 'log'),
                            prior_PD = TRUE,
                            refresh = 0,
                            chains = 3, iter = 5000, thin = 5, warmup = 2000)


## ----fitModel1g1, results='markdown', eval=TRUE, mhidden=TRUE-----------------
quinn_rstanarm1 |>
    ggpredict(~fSEASON+DENSITY) |>
    plot(show_data = TRUE)
quinn_rstanarm1 |>
    ggpredict(~fSEASON+DENSITY) |>
    plot(show_data = TRUE) |>
    wrap_plots() &
    scale_y_log10()
## although, since there are zeros...
quinn_rstanarm1 |>
    ggpredict(~fSEASON+DENSITY) |>
    plot(show_data = TRUE, jitter = FALSE) |>
    wrap_plots() &
    scale_y_continuous(trans = scales::pseudo_log_trans())


## ----fitModel1g2, results='markdown', eval=TRUE, mhidden=TRUE-----------------
quinn_rstanarm1 |>
    ggemmeans(~fSEASON+DENSITY) |>
    plot(show_data=TRUE) |>
    plot(show_data = TRUE) |>
    wrap_plots() &
    scale_y_log10()
## although, since there are zeros...
quinn_rstanarm1 |>
    ggemmeans(~fSEASON+DENSITY) |>
    plot(show_data = TRUE, jitter = FALSE) |>
    wrap_plots() &
    scale_y_continuous(trans = scales::pseudo_log_trans())


## ----fitModel1h, results='markdown', eval=TRUE, mhidden=TRUE------------------
quinn |> group_by(fSEASON, DENSITY) |>
    summarise(Mean = log(mean(RECRUITS)),
              SD = log(sd(RECRUITS)))
log(sd(quinn$RECRUITS))/
    apply(model.matrix(~fSEASON*DENSITY, data = quinn), 2, sd)
quinn_rstanarm2 <- stan_glm(RECRUITS~fSEASON*DENSITY, data = quinn,
                            family = poisson(link = 'log'),
                            prior_intercept = normal(2.3, 2, autoscale = FALSE),
                            prior = normal(0, 10, autoscale = FALSE),
                            prior_PD = TRUE,
                            refresh = 0,
                            chains = 3, iter = 5000, thin = 5, warmup = 2000)


## ----fitModel1i1, results='markdown', eval=TRUE, mhidden=TRUE-----------------
quinn_rstanarm2 |>
    ggpredict(~fSEASON+DENSITY) |>
    plot(show_data = TRUE)
quinn_rstanarm2 |>
    ggpredict(~fSEASON+DENSITY) |>
    plot(show_data = TRUE) |>
    wrap_plots() &
    scale_y_log10()
## although, since there are zeros...
quinn_rstanarm2 |>
    ggpredict(~fSEASON+DENSITY) |>
    plot(show_data = TRUE, jitter = FALSE) |>
    wrap_plots() &
    scale_y_continuous(trans = scales::pseudo_log_trans())


## ----fitModel1i2, results='markdown', eval=TRUE, mhidden=TRUE-----------------
quinn_rstanarm2 |>
    ggemmeans(~fSEASON+DENSITY) |>
    plot(show_data=TRUE) |>
    plot(show_data = TRUE) |>
    wrap_plots() &
    scale_y_log10()
## although, since there are zeros...
quinn_rstanarm2 |>
    ggemmeans(~fSEASON+DENSITY) |>
    plot(show_data = TRUE, jitter = FALSE) |>
    wrap_plots() &
    scale_y_continuous(trans = scales::pseudo_log_trans())


## ----fitModel1j, results='markdown', eval=TRUE, mhidden=TRUE------------------
quinn_rstanarm3 <- quinn_rstanarm2 |> update(prior_PD = FALSE)


## ----modelFit1k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
quinn_rstanarm3 |> posterior_vs_prior(color_by = 'vs', group_by = TRUE,
                   facet_args = list(scales = 'free_y'))


## ----modelFit1l1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
quinn_rstanarm3 |>
    ggpredict(~fSEASON+DENSITY) |>
    plot(show_data = TRUE)


## ----modelFit1l2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
quinn_rstanarm3 |>
    ggemmeans(~fSEASON+DENSITY) |>
    plot(show_data=TRUE)


## ----fitModel2a, results='markdown', eval=TRUE, mhidden=TRUE------------------
quinn_form <- bf(RECRUITS ~ fSEASON*DENSITY,  family = poisson(link = 'log'))


## ----fitModel2a1, results='markdown', eval=TRUE, mhidden=TRUE, echo=2---------
options(width=100)
get_prior(quinn_form,  data = quinn)
options(width=80)


## ----fitModel2d1, results='markdown', eval=TRUE-------------------------------
quinn |>
  group_by(fSEASON, DENSITY) |>
  summarise(Mean = mean(log(RECRUITS+1)),
    Median = median(log(RECRUITS+1)),
    SD = sd(log(RECRUITS+1)),
    MAD = mad(log(RECRUITS+1))
    )


## ----fitModel2h, results='markdown', eval=TRUE, mhidden=TRUE------------------
quinn_form <- bf(RECRUITS ~ fSEASON*DENSITY,  family = poisson(link = 'log'))
priors <- prior(normal(2.5, 0.2), class = 'Intercept') +
    prior(normal(0, 2.5), class = 'b')
quinn_brm2 <- brm(quinn_form,
                  data = quinn,
                  prior = priors,
                  sample_prior = "only",
                  refresh = 0,
                  chains = 3,
                  iter = 5000,
                  thin = 5,
                  warmup = 2500,
                  backend = 'cmdstanr')


## ----fitModel2i1, results='markdown', eval=TRUE, mhidden=TRUE-----------------
quinn_brm2 |>
    ggpredict(~fSEASON+DENSITY) |>
    plot(show_data=TRUE)
quinn_brm2 |>
    ggpredict(~fSEASON+DENSITY) |>
    plot(show_data=TRUE) +
    scale_y_log10()
## Or since there are zeros
quinn_brm2 |>
    ggpredict(~fSEASON+DENSITY) |>
    plot(show_data=TRUE) +
    scale_y_continuous(trans = scales::pseudo_log_trans())


## ----fitModel2i2, results='markdown', eval=TRUE, mhidden=TRUE-----------------
quinn_brm2 |>
    ggemmeans(~fSEASON+DENSITY) |>
    plot(show_data=TRUE)
quinn_brm2 |>
    ggemmeans(~fSEASON+DENSITY) |>
    plot(show_data=TRUE) +
    scale_y_log10()
## Or since there are zeros
quinn_brm2 |>
    ggemmeans(~fSEASON+DENSITY) |>
    plot(show_data=TRUE) +
    scale_y_continuous(trans = scales::pseudo_log_trans())


## ----fitModel2i3, results='markdown', eval=TRUE, mhidden=TRUE-----------------
quinn_brm2 |>
    conditional_effects('fSEASON:DENSITY') |>
    plot(points=TRUE)
quinn_brm2 |>
    conditional_effects('fSEASON:DENSITY') |>
    plot(points=TRUE) |>
    wrap_plots() &
    scale_y_continuous(trans = scales::pseudo_log_trans())


## ----fitModel2j, results='markdown', eval=TRUE, mhidden=TRUE------------------
quinn_brmP <- quinn_brm2 |> update(sample_prior = 'yes', refresh = 0)


## ----fitModel2k1, results='markdown', eval=TRUE, mhidden=TRUE-----------------
quinn_brmP |> get_variables()
quinn_brmP |> hypothesis('fSEASONSummer<0') |> plot()
quinn_brmP |> hypothesis('DENSITYLow<0') |> plot()
quinn_brmP |> hypothesis('fSEASONSummer:DENSITYLow<0') |> plot()


## ----fitModel2k2, results='markdown', out.width = 600, fig.width = 8, fig.height = 4, eval=TRUE, mhidden=TRUE----
quinn_brmP |> SUYR_prior_and_posterior()
quinn_brmP |>
  posterior_samples() |>
  dplyr::select(-`lp__`) |>
  pivot_longer(everything(), names_to = 'key') |>
  mutate(Type = ifelse(str_detect(key, 'prior'), 'Prior', 'b'),
         Class = ifelse(str_detect(key, 'Intercept'),  'Intercept',
               ifelse(str_detect(key, 'b'),  'b', 'sigma')),
         Par = str_replace(key, 'b_', '')) |>
  ggplot(aes(x = Type,  y = value, color = Par)) +
  stat_pointinterval(position = position_dodge())+
  facet_wrap(~Class,  scales = 'free')


## ----fitModel2l, results='markdown', eval=TRUE, mhidden=TRUE------------------
quinn_brmP |> standata()
quinn_brmP |> stancode()


## ----modelValidation1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_mcmc()


## ----modelValidation1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(quinn_rstanarm3, plotfun='mcmc_trace')


## ----modelValidation1c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(quinn_rstanarm3, 'acf_bar')


## ----modelValidation1d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(quinn_rstanarm3, 'rhat_hist')


## ----modelValidation1e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(quinn_rstanarm3, 'neff_hist')


## ----Validation1f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(quinn_rstanarm3, 'combo')
plot(quinn_rstanarm3, 'violin')


## ----modelValidation1g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_trace(quinn_rstanarm3)


## ----modelValidation1h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_ac(quinn_rstanarm3)


## ----modelValidation1i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_rhat(quinn_rstanarm3)


## ----modelValidation1j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_ess(quinn_rstanarm3)


## ----modelValidation1k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_dens(quinn_rstanarm3, separate_chains = TRUE)


## ----modelValidation1l, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
quinn_ggs <- ggs(quinn_rstanarm3, burnin = FALSE, inc_warmup = FALSE)
ggs_traceplot(quinn_ggs)


## ----modelValidation1m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_autocorrelation(quinn_ggs)


## ----modelValidation1n, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_Rhat(quinn_ggs)


## ----modelValidation1o, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_effective(quinn_ggs)


## ----modelValidation1p, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_crosscorrelation(quinn_ggs)


## ----modelValidation1q, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_grb(quinn_ggs)


## ----modelValidation2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
quinn_brmP$fit |> stan_trace()
quinn_brmP$fit |> stan_trace(inc_warmup=TRUE)


## ----modelValidation2h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
quinn_brmP$fit |> stan_ac()


## ----modelValidation2i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
quinn_brmP$fit |> stan_rhat()


## ----modelValidation2j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
quinn_brmP$fit |> stan_ess()


## ----modelValidation2k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
quinn_brmP$fit |> stan_dens(separate_chains = TRUE)


## ----modelValidation3a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation3b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(quinn_rstanarm3,  plotfun='dens_overlay')


## ----modelValidation3c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#pp_check(quinn_rstanarm3, plotfun='error_scatter_avg')


## ----modelValidation3d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(quinn_rstanarm3, x=as.numeric(quinn$fSEASON), plotfun='error_scatter_avg_vs_x')
pp_check(quinn_rstanarm3, x=as.numeric(quinn$DENSITY), plotfun='error_scatter_avg_vs_x')


## ----modelValidation3e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(quinn_rstanarm3, x=as.numeric(quinn$fSEASON), plotfun='intervals')


## ----modelValidation3g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(quinn_rstanarm3)


## ----modelValidation4a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
preds <- posterior_predict(quinn_rstanarm3,  nsamples=250,  summary=FALSE)
quinn_resids <- createDHARMa(simulatedResponse = t(preds),
                            observedResponse = quinn$RECRUITS,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse = TRUE)
plot(quinn_resids)


## ----modelValidation5a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation5b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
quinn_brmP |> pp_check(type = 'dens_overlay', ndraws = 100)


## ----modelValidation5c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#quinn_brmP |> pp_check(type = 'error_scatter_avg')


## ----modelValidation5e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
quinn_brmP |> pp_check(type='intervals')
## quinn_brmP |> pp_check(group='DENSITY', type='intervals')


## ----modelValidation5g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(quinn_brmP)


## ----modelValidation6aa, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=10----
quinn_resids <- make_brms_dharma_res(quinn_brmP, integerResponse = TRUE)
wrap_elements(~testUniformity(quinn_resids)) +
               wrap_elements(~plotResiduals(quinn_resids, form = factor(rep(1, nrow(quinn))))) +
               wrap_elements(~plotResiduals(quinn_resids, quantreg = TRUE)) +
               wrap_elements(~testDispersion(quinn_resids))



## ----modelValidation6a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
preds <- quinn_brmP |> posterior_predict(nsamples = 250,  summary = FALSE)
quinn_resids <- createDHARMa(simulatedResponse = t(preds),
                            observedResponse = quinn$RECRUITS,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse = TRUE)
quinn_resids |> plot()

quinn_resids |> testDispersion()
quinn_resids |> testZeroInflation()


## ----fitModel3a, results='markdown', eval=TRUE, mhidden=TRUE------------------
quinn_rstanarmNB <- stan_glm(RECRUITS~fSEASON*DENSITY, data = quinn,
                            family = neg_binomial_2(link = 'log'),
                            prior_intercept = normal(2.3, 2, autoscale = FALSE),
                            prior = normal(0, 10, autoscale = FALSE),
                            prior_aux = rstanarm::exponential(rate = 1, autoscale = FALSE),
                            prior_PD = FALSE,
                            refresh = 0,
                            chains = 3, iter = 5000, thin = 5, warmup = 2000)


## ----fitModel3b, results='markdown', eval=TRUE, mhidden=TRUE------------------
posterior_vs_prior(quinn_rstanarmNB, color_by='vs', group_by=TRUE,
                   facet_args=list(scales='free_y'))


## ----fitModel3b2, results='markdown', eval=TRUE, mhidden=TRUE-----------------
quinn_rstanarmNB |> ggpredict(~fSEASON+DENSITY) |> plot(show_data = TRUE)


## ----fitModel3b3, results='markdown', eval=TRUE, mhidden=TRUE-----------------
quinn_rstanarmNB |>
    ggemmeans(~fSEASON+DENSITY, back.transform = TRUE) |>
    plot(show_data=TRUE)


## ----fitModel3b4, results='markdown', eval=TRUE, mhidden=TRUE-----------------
quinn_rstanarmNB |> plot('mcmc_trace')
quinn_rstanarmNB |> plot('mcmc_acf_bar')
quinn_rstanarmNB |> plot('mcmc_rhat_hist')
quinn_rstanarmNB |> plot('mcmc_neff_hist')


## ----fitModel3b5, results='markdown', eval=TRUE, mhidden=TRUE, fig.width = 8, fig.height = 4, out.width = 600----
preds <- posterior_predict(quinn_rstanarmNB,  nsamples=250,  summary=FALSE)
quinn_resids <- createDHARMa(simulatedResponse = t(preds),
                            observedResponse = quinn$RECRUITS,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse=TRUE)
plot(quinn_resids)
quinn_resids |> testDispersion()


## ----fitModel3c, results='markdown', eval=TRUE, mhidden=TRUE------------------
(loo.P = loo(quinn_rstanarmP))
(loo.NB = loo(quinn_rstanarmNB))
loo_compare(loo.P, loo.NB)


## ----fitModel4a, results='markdown', eval=TRUE, mhidden=TRUE------------------
quinn_form <- bf(RECRUITS ~ fSEASON*DENSITY,  family = negbinomial(link = 'log'))
get_prior(quinn_form,  data = quinn)

priors <- prior(normal(2.5, 0.2), class = 'Intercept') +
    prior(normal(0, 2.5), class = 'b') +
    prior(gamma(0.01, 0.01), class = "shape")
quinn_brmsNB <- brm(quinn_form,
                    data = quinn,
                    prior = priors,
                    refresh = 0,
                    chains = 3,
                    iter = 5000,
                    thin = 5,
                    warmup = 2500,
                    backend = "cmdstanr")


## ----fitModel4a1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width = 8, fig.height = 4, out.width = 600----
quinn_resids <- make_brms_dharma_res(quinn_brmsNB, integerResponse = TRUE)
wrap_elements(~testUniformity(quinn_resids)) +
               ## wrap_elements(~plotResiduals(quinn_resids, form = factor(rep(1, nrow(quinn))))) +
               wrap_elements(~plotResiduals(quinn_resids, quantreg = TRUE)) +
               wrap_elements(~testDispersion(quinn_resids))



## ----partialPlot1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_rstanarmNB |> ggpredict(~fSEASON+DENSITY) |> plot(show_data=TRUE)


## ----partialPlot1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_rstanarmNB |>
    ggemmeans(~fSEASON|DENSITY,  type='fixed', back.transform = TRUE) |>
    plot(show_data=TRUE)


## ----partialPlot1c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_rstanarmNB |>
    fitted_draws(newdata=quinn) |>
    median_hdci() |>
    ggplot(aes(x=fSEASON, colour=DENSITY, y=.value)) +
    geom_pointrange(aes(ymin=.lower, ymax=.upper), position = position_dodge(width=0.2)) +
    geom_line(position = position_dodge(width=0.2)) +
    geom_point(data=quinn,  aes(y=RECRUITS,  x=fSEASON, colour = DENSITY), position = position_dodge(width=0.2))


## ----partialPlot2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_brmsNB |> conditional_effects("fSEASON:DENSITY") |> plot(points = TRUE)


## ----partialPlot2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_brmsNB |> ggpredict(~fSEASON+DENSITY) |> plot(show_data = TRUE)


## ----partialPlot2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_brmsNB |> ggemmeans(~fSEASON|DENSITY) |> plot(show_data = TRUE)


## ----partialPlot2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_brmsNB |>
    fitted_draws(newdata=quinn) |>
    median_hdci() |>
    ggplot(aes(x=fSEASON, colour=DENSITY, y=.value)) +
    geom_pointrange(aes(ymin=.lower, ymax=.upper), position = position_dodge(width=0.2)) +
    geom_line(position = position_dodge(width=0.2)) +
    geom_point(data=quinn,  aes(y=RECRUITS,  x=fSEASON, colour = DENSITY), position = position_dodge(width=0.2))


## ----summariseModel1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_rstanarmNB |> summary()


## ----summariseModel1a1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=FALSE----
quinn_sum <- summary(quinn_rstanarmNB)


## ----summariseModel1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tidyMCMC(quinn_rstanarmNB$stanfit, estimate.method='median',  conf.int=TRUE,  conf.method='HPDinterval',  rhat=TRUE, ess=TRUE)

## ----summariseModel1b1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
quinn_tidy <- tidyMCMC(quinn_rstanarmNB$stanfit, estimate.method='median',  conf.int=TRUE,  conf.method='HPDinterval',  rhat=TRUE, ess=TRUE)


## ----summariseModel1dd, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_rstanarmNB$stanfit |>
    summarise_draws(median,
                    HDInterval::hdi,
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel1d2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_rstanarmNB$stanfit |>
    summarise_draws(median,
                    ~HDInterval::hdi(.x, credMass = 0.9),
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel1d3, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_rstanarmNB$stanfit |>
    summarise_draws(
        ~ median(exp(.x)),
        ~HDInterval::hdi(exp(.x)),
        rhat, length, ess_bulk, ess_tail)


## ----summariseModel1m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_rstanarmNB$stanfit |> as_draws_df()
quinn_rstanarmNB$stanfit |>
    as_draws_df() |>
    summarise_draws(
        median,
        ~ HDInterval::hdi(.x),
        rhat,
        ess_bulk
    )

quinn_rstanarmNB$stanfit |>
    as_draws_df() |>
    exp() |>
    summarise_draws(
        median,
        ~ HDInterval::hdi(.x),
        rhat,
        ess_bulk
    )


## ----summariseModel1c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_rstanarmNB |> get_variables()
quinn_draw <- quinn_rstanarmNB |> gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`,  regex=TRUE)
quinn_draw

exceedP <- function(x, Val = 0) mean(x>Val)

quinn_rstanarmNB |>
    gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`,  regex=TRUE) |>
    mutate(.value = exp(.value)) |>
    summarise_draws(median,
                    HDInterval::hdi,
                    rhat,
                    length,
                    ess_bulk,
                    ess_tail,
                    ~ exceedP(.x, 1))


## ----summariseModel1c1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_draw |> median_hdci()


## ----summariseModel1c8, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_draw |> median_hdci(exp(.value))


## ----summariseModel1c3, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
quinn_gather <- quinn_rstanarmNB |> gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`,  regex=TRUE) |>
  median_hdci()


## ----summariseModel1c4, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
quinn_rstanarmNB |>
  gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
  ggplot() +
  stat_halfeye(aes(x=.value,  y=.variable)) +
  facet_wrap(~.variable, scales='free')
quinn_rstanarmNB |>
  gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
  ggplot() +
    geom_vline(xintercept=0, linetype='dashed') +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    theme_classic()


## ----summariseModel1c5, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
quinn_rstanarmNB |>
  gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
  group_by(.variable) |>
  mutate(.value=exp(.value)) |>
  median_hdci()
quinn_rstanarmNB |>
  gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
  mutate(.value=exp(.value)) |>
  ggplot() +
    geom_vline(xintercept=1, linetype='dashed') +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    scale_x_continuous('', trans = scales::log2_trans(), breaks=unique(as.vector(2^(0:4 %o% c(-1,1))))) +
    theme_classic()


## ----summariseModel1j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_rstanarmNB |> plot(plotfun='mcmc_intervals')


## ----summariseModel2d5, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
## Link scale
quinn_rstanarmNB |>
    gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
    ggplot() +
    stat_slab(aes(x = .value, y = .variable,
                  fill = stat(ggdist::cut_cdf_qi(cdf,
                           .width = c(0.5, 0.8, 0.95),
                           labels = scales::percent_format())
                           )), color='black') +
    geom_vline(xintercept=0, linetype='dashed') +
    scale_fill_brewer('Interval', direction = -1, na.translate = FALSE)
## Fractional scale
quinn_rstanarmNB |>
    gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
    mutate(.value=exp(.value)) |>
    ggplot() +
    stat_slab(aes(x = .value, y = .variable,
                  fill = stat(ggdist::cut_cdf_qi(cdf,
                           .width = c(0.5, 0.8, 0.95),
                           labels = scales::percent_format())
                           )), color='black') +
    geom_vline(xintercept=1, linetype='dashed') +
    scale_fill_brewer('Interval', direction = -1, na.translate = FALSE) +
    scale_x_continuous(trans = scales::log2_trans())


## ----summariseModel1c44, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
quinn_rstanarmNB |>
  gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
  ggplot() +
  stat_halfeye(aes(x=.value,  y=.variable)) +
  facet_wrap(~.variable, scales='free')

quinn_rstanarmNB |>
  gather_draws(`.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
  ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    geom_vline(xintercept = 0, linetype = 'dashed')

quinn_rstanarmNB |>
  gather_draws(`.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
  ggplot() +
    stat_halfeye(aes(x=exp(.value),  y=.variable)) +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    scale_x_continuous(trans = scales::log2_trans())


## ----summariseModel1c7, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
quinn_rstanarmNB |>
  gather_draws(`.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
  ggplot() +
    geom_density_ridges(aes(x=.value, y = .variable), alpha=0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed')
##Or on a fractional scale
quinn_rstanarmNB |>
  gather_draws(`.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
  ggplot() +
    geom_density_ridges_gradient(aes(x=exp(.value),
                                     y = .variable,
                                     fill = stat(x)),
                                 alpha=0.4, colour = 'white',
                                 quantile_lines = TRUE,
                                 quantiles = c(0.025, 0.975)) +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    scale_x_continuous(trans = scales::log2_trans()) +
    scale_fill_viridis_c(option = "C")


## ----summariseModel1d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_rstanarmNB |> tidy_draws()


## ----summariseModel1e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_rstanarmNB |> spread_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`,  regex=TRUE)


## ----summariseModel1f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_rstanarmNB |> posterior_samples() |> as_tibble()


## ----summariseModel1g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
#quinn_rstanarmNB |> bayes_R2() |> median_hdci


## ----summariseModel2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_brmsNB |> summary()


## ----summariseModel2a1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=FALSE----
quinn_sum <- quinn_brmsNB |> summary()
quinn_sum <- quinn_sum$fixed


## ----summariseModel2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_brmsNB$fit |>
    tidyMCMC(estimate.method = 'median',
             conf.int = TRUE,
             conf.method = 'HPDinterval',
             rhat = TRUE,
             ess = TRUE)

## ----summariseModel2b1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
quinn_tidy <- tidyMCMC(quinn_brmsNB$fit, estimate.method='median',  conf.int=TRUE,  conf.method='HPDinterval',  rhat=TRUE, ess=TRUE)


## ----summariseModel2m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_brmsNB |> as_draws_df()
quinn_brmsNB |>
  as_draws_df() |>
  summarise_draws(
    median,
    ~ HDInterval::hdi(.x),
    rhat,
    ess_bulk, ess_tail
  )

quinn_brmsNB |>
    as_draws_df() |>
    exp() |>
    summarise_draws(
        median,
        HDInterval::hdi,
        rhat,
        length,
        ess_bulk, ess_tail
    )


## ----summariseModel2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_brmsNB |> get_variables()
quinn_draw <- quinn_brmsNB |>
    gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`,  regex = TRUE)
quinn_draw

quinn_brmsNB |>
    gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`,  regex = TRUE) |>
    mutate(.value = exp(.value)) |>
    summarise_draws(median,
                    ~HDInterval::hdi(.x, credMass = 0.95),
                    rhat,
                    length,
                    ess_bulk, ess_tail)


exceedP <- function(x, Val = 0) mean(x>Val)
quinn_brmsNB |>
    tidy_draws() |>
    exp() |>
    dplyr::select(starts_with("b_")) |>
    summarise_draws(median,
                    ~HDInterval::hdi(.x, credMass = 0.9),
                    rhat,
                    ess_bulk, ess_tail,
                    ~exceedP(.x, 1))


## ----summariseModel2c1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_draw |> median_hdci()


## ----summariseModel2c5, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
quinn_draw |>
  median_hdci(exp(.value))


## ----summariseModel2c3, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
quinn_gather <- quinn_brmsNB |> gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`,  regex = TRUE) |>
  median_hdci()


## ----summariseModel2c45, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
quinn_brmsNB |>
    gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
    ggplot() +
    geom_vline(xintercept=0, linetype='dashed') +
    stat_slab(aes(x = .value, y = .variable,
                  fill = stat(ggdist::cut_cdf_qi(cdf,
                           .width = c(0.5, 0.8, 0.95),
                           labels = scales::percent_format())
                           )), color='black') +
    scale_fill_brewer('Interval', direction = -1, na.translate = FALSE)

quinn_brmsNB |>
  gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
  ggplot() +
    geom_vline(xintercept=0, linetype='dashed') +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    theme_classic()


## ----summariseModel2c55, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
quinn_brmsNB |>
  gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
  group_by(.variable) |>
  mutate(.value=exp(.value)) |>
  median_hdci()

quinn_brmsNB |>
    gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
    mutate(.value=exp(.value)) |>
    ggplot() +
    geom_vline(xintercept=1, linetype='dashed') +
    stat_slab(aes(x = .value, y = .variable,
                  fill = stat(ggdist::cut_cdf_qi(cdf,
                           .width = c(0.5, 0.8, 0.95),
                           labels = scales::percent_format())
                           )), color='black') +
    scale_x_continuous('', trans = scales::log2_trans(), breaks=unique(as.vector(2^(0:4 %o% c(-1,1))))) +
    scale_fill_brewer('Interval', direction = -1, na.translate = FALSE) +
    theme_classic()


## ----summariseModel2j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_brmsNB$fit |> plot(type='intervals')


## ----summariseModel2d55, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
## Link scale
quinn_brmsNB |>
    gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
    ggplot() +
    stat_slab(aes(x = .value, y = .variable,
                  fill = stat(ggdist::cut_cdf_qi(cdf,
                           .width = c(0.5, 0.8, 0.95),
                           labels = scales::percent_format())
                           )), color='black') +
    geom_vline(xintercept=0, linetype='dashed') +
    scale_fill_brewer('Interval', direction = -1, na.translate = FALSE)
## Fractional scale
quinn_brmsNB |>
    gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
    mutate(.value=exp(.value)) |>
    ggplot() +
    stat_slab(aes(x = .value, y = .variable,
                  fill = stat(ggdist::cut_cdf_qi(cdf,
                           .width = c(0.5, 0.8, 0.95),
                           labels = scales::percent_format())
                           )), color='black') +
    geom_vline(xintercept=1, linetype='dashed') +
    scale_fill_brewer('Interval', direction = -1, na.translate = FALSE) +
    scale_x_continuous(trans = scales::log2_trans())


## ----summariseModel2c4, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
quinn_brmsNB |>
    gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
    ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    facet_wrap(~.variable, scales='free')

quinn_brmsNB |>
    gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
    ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    geom_vline(xintercept = 0, linetype = 'dashed')

quinn_brmsNB |>
    gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
    ggplot() +
    stat_halfeye(aes(x=exp(.value),  y=.variable)) +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    scale_x_continuous(trans = scales::log2_trans())


## ----summariseModel2c7, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
quinn_brmsNB |>
    gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
    ggplot() +
    geom_density_ridges(aes(x=.value, y = .variable), alpha=0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed')
##Or on a fractional scale
quinn_brmsNB |>
    gather_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`, regex=TRUE) |>
    ggplot() +
    geom_density_ridges_gradient(aes(x=exp(.value),
                                     y = .variable,
                                     fill = stat(x)),
                                 alpha=0.4, colour = 'white',
                                 quantile_lines = TRUE,
                                 quantiles = c(0.025, 0.975)) +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    scale_x_continuous(trans = scales::log2_trans()) +
    scale_fill_viridis_c(option = "C")


## ----summariseModel2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_brmsNB |> tidy_draws()


## ----summariseModel2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_brmsNB |> spread_draws(`.Intercept.*|.*fSEASON.*|.*DENSITY.*`,  regex=TRUE)


## ----summariseModel2f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_brmsNB |> posterior_samples() |> as_tibble()


## ----summariseModel2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
quinn_brmsNB |> bayes_R2(summary=FALSE) |> median_hdci()


## ----summariseModel2k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
0.1 * log(sd(quinn$RECRUITS))
quinn_brmsNB |> rope(range = c(-0.28, 0.28))
rope(quinn_brmsNB, range = c(-0.28, 0.28)) |> plot()

## Or based on fractional scale
quinn_brmsNB |>
    as_draws_df('^b_fSEASON.*|^b_DENSITY.*', regex = TRUE) |>
    exp() |>
    ## equivalence_test(range = c(0.755, 1.32))
    rope(range = c(0.755, 1.32))
quinn_mcmc <-
    quinn_brmsNB |>
    as_draws_df('^b_fSEASON.*|^b_DENSITY.*', regex = TRUE) |>
    exp()
quinn_mcmc |>
    rope(range = c(0.755, 1.32))
## note, the following is not quit correct, it does not get the CI correct
quinn_mcmc |>
    rope(range = c(0.755, 1.32)) |>
    plot(data = quinn_brmsNB)


quinn_mcmc |>
    equivalence_test(range = c(0.755, 1.32))


## -----------------------------------------------------------------------------
#| label: modelsummary
#| results: markup
#| eval: true
#| echo: true
#| cache: false
quinn_brmsNB |> modelsummary(
  statistic = c("conf.low", "conf.high"),
  shape = term ~ statistic,
  exponentiate = TRUE
)


## -----------------------------------------------------------------------------
#| label: modelsummary_plot
#| results: markup
#| eval: true
#| echo: true
#| cache: false
quinn_brmsNB |> modelplot(coef_omit = "shape", exponentiate = TRUE)


## ----mainEffects1a, results='markdown', eval=TRUE, mhidden=TRUE---------------
## fold scale
quinn_rstanarmNB |>
    emmeans(~DENSITY|fSEASON, type='response') |>
    pairs()
## absolute response scale
quinn_rstanarmNB |>
    emmeans(~DENSITY|fSEASON, type='link') |>
    regrid() |>
    pairs()


## ----mainEffects1b, results='markdown', eval=TRUE, mhidden=TRUE---------------
quinn_em <- quinn_rstanarmNB |>
    emmeans(~DENSITY|fSEASON, type='link') |>
    pairs() |>
    gather_emmeans_draws() |>
    mutate(Fit=exp(.value))
head(quinn_em)

g2 <- quinn_em |>
  group_by(contrast, fSEASON) |>
  median_hdci() |>
  ggplot() +
  geom_vline(xintercept=1, linetype='dashed') +
  geom_pointrange(aes(x=Fit, y=fSEASON, xmin=Fit.lower, xmax=Fit.upper)) +
  scale_x_continuous('Effect size (High/Low)', trans = scales::log2_trans(), breaks=unique(as.vector(2^(0:4 %o% c(-1,1))))) +
  theme_classic()
g2

ggplot(quinn_em, aes(x=Fit)) +
    geom_histogram() +
    geom_vline(xintercept = 1, linetype='dashed') +
    scale_x_continuous('Effect size (High/Low)', trans = scales::log2_trans(), breaks=unique(as.vector(2^(0:4 %o% c(-1,1))))) +
    facet_wrap(fSEASON~contrast, scales='free')
quinn_em |> group_by(contrast, fSEASON) |> median_hdci(Fit)
# Probability of effect
quinn_em |> group_by(contrast,fSEASON) |> summarize(P=mean(Fit>1))
##Probability of effect greater than 10%
quinn_em |> group_by(contrast,fSEASON) |> summarize(P=mean(Fit>1.1))

## Using summarise_draws
quinn_rstanarmNB |>
    emmeans(~DENSITY|fSEASON, type='link') |>
    pairs() |>
    tidy_draws() |>
    exp() |>
    summarise_draws(median,
                    HDInterval::hdi,
                    P = ~ mean(.x > 1)
                    )


## ----mainEffects1c, results='markdown', eval=TRUE, mhidden=TRUE---------------
newdata <- with(quinn, expand.grid(fSEASON = levels(fSEASON),
                                  DENSITY = levels(DENSITY)))
Xmat<- model.matrix(~fSEASON*DENSITY, data = newdata)
as.matrix(quinn_rstanarmNB) |> head()
## coefs <- as.matrix(quinn_rstanarmNB)
coefs <- as.matrix(as.data.frame(quinn_rstanarmNB) |>
                  dplyr:::select(-reciprocal_dispersion)) |>
    as.matrix()
fit <- exp(coefs %*% t(Xmat))
newdata <- newdata |>
    cbind(tidyMCMC(fit, conf.int = TRUE, conf.method = 'HPDinterval'))
head(newdata)

ggplot(newdata, aes(y = estimate, x = fSEASON, fill = DENSITY)) +
    geom_blank() +
    geom_line(aes(x=as.numeric(fSEASON), ymin=conf.low, ymax=conf.high, linetype=DENSITY),
              position = position_dodge(0.2))+
    geom_pointrange(aes(ymin=conf.low, ymax=conf.high), shape=21,
                    position = position_dodge(0.2))

#Compare high and low in each season
#via contrasts
newdata <- with(quinn, expand.grid(fSEASON = levels(fSEASON),
                                   DENSITY = levels(DENSITY)))
## factor differences
Xmat<- model.matrix(~fSEASON*DENSITY, data=newdata)
Xmat.high <- Xmat[newdata$DENSITY=="High",]
Xmat.low <- Xmat[newdata$DENSITY=="Low",]
Xmat.density <- Xmat.high-Xmat.low
rownames(Xmat.density) <- levels(quinn$fSEASON)
coefs = as.matrix(as.data.frame(quinn_rstanarmNB) |> dplyr:::select(-reciprocal_dispersion))
fit = exp(coefs %*% t(Xmat.density))
tidyMCMC(fit, conf.int=TRUE, conf.method='HPDinterval')
## or absolute
fit.high = coefs %*% t(Xmat.high)
fit.low = coefs %*% t(Xmat.low)
fit = exp(fit.high) - exp(fit.low)
#fit = exp(fit.high - fit.low)
tidyMCMC(fit, conf.int=TRUE, conf.method='HPDinterval')


## ----mainEffects2a, results='markdown', eval=TRUE, mhidden=TRUE---------------
quinn_brmsNB |>
    emmeans(~DENSITY|fSEASON, type='response') |>
    pairs()
## absolute response scale
quinn_brmsNB |>
    emmeans(~DENSITY|fSEASON, type='link') |>
    regrid() |>
    pairs()


## ----mainEffects2b, results='markdown', eval=TRUE, mhidden=TRUE---------------
quinn_em <- quinn_brmsNB |>
    emmeans(~DENSITY|fSEASON, type='link') |>
    pairs() |>
    gather_emmeans_draws() |>
    mutate(Fit=exp(.value))
head(quinn_em)

g2 <- quinn_em |>
  group_by(contrast, fSEASON) |>
  median_hdci() |>
  ggplot() +
  geom_vline(xintercept=1, linetype='dashed') +
  geom_pointrange(aes(x=Fit, y=fSEASON, xmin=Fit.lower, xmax=Fit.upper)) +
  scale_x_continuous('Effect size (High/Low)', trans = scales::log2_trans(), breaks=unique(as.vector(2^(0:4 %o% c(-1,1))))) +
  theme_classic()
g2

ggplot(quinn_em, aes(x=Fit)) +
    geom_histogram() +
    geom_vline(xintercept = 1, linetype='dashed') +
    scale_x_continuous('Effect size (High/Low)', trans = scales::log2_trans(), breaks=unique(as.vector(2^(0:4 %o% c(-1,1))))) +
    facet_wrap(fSEASON~contrast, scales='free')
quinn_em |> group_by(contrast, fSEASON) |> median_hdci(Fit)
# Probability of effect
quinn_em |> group_by(contrast,fSEASON) |> summarize(P=mean(Fit>1))
##Probability of effect greater than 10%
quinn_em |> group_by(contrast,fSEASON) |> summarize(P=mean(Fit>1.1))



## ----summaryFig1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=4----
newdata <- quinn_rstanarmNB |>
    emmeans(~fSEASON|DENSITY, type='response') |>
    as.data.frame()
head(newdata)
g1 <- ggplot(newdata, aes(y=prob, x=fSEASON, color=DENSITY)) +
    geom_pointrange(aes(ymin=lower.HPD, ymax=upper.HPD),
                    position=position_dodge(width=0.2)) +
    theme_classic()
g1 + g2


## ----summaryFig2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=4----
newdata <- quinn_brmsNB %>%
    emmeans(~fSEASON|DENSITY, type='response') |>
    as.data.frame()
head(newdata)
g1 <- ggplot(newdata, aes(y=prob, x=fSEASON, color=DENSITY)) +
    geom_pointrange(aes(ymin=lower.HPD, ymax=upper.HPD),
                    position=position_dodge(width=0.2)) +
    theme_classic()
g1 + g2


## ----fitModel.brms, results='markdown', eval=TRUE, mhidden=TRUE---------------
quinn <- quinn |>
  group_by(fSEASON, DENSITY) |>
  mutate(Obs = factor(1:n()))

quinn_form <- bf(RECRUITS ~ fSEASON*DENSITY + (1|Obs),  family = poisson(link = 'log'))
get_prior(quinn_form,  data = quinn)

quinn_brmsU <- brm(quinn_form,
                   data = quinn,
                   refresh = 0,
                   chains = 3,
                   iter = 5000,
                   thin = 5,
                   warmup = 2000)

preds <- posterior_predict(quinn_brmsU,  nsamples=250,  summary=FALSE)
quinn_resids <- createDHARMa(simulatedResponse = t(preds),
                            observedResponse = quinn$RECRUITS,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse = TRUE)
plot(quinn_resids)
newdata = emmeans(quinn_brmsU, ~fSEASON|DENSITY, type='response') |> as.data.frame()
newdata
ggplot(newdata, aes(y=rate, x=fSEASON, color=DENSITY)) +
    geom_pointrange(aes(ymin=lower.HPD, ymax=upper.HPD),
                    position=position_dodge(width=0.2))


