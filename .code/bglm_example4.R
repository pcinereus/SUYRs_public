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
library(patchwork)     #for multiple plots
library(modelsummary)  #for data and model summaries
library(car)           #for scatterplot matrices
library(ggridges)      #for ridge plots
theme_set(theme_grey()) #put the default ggplot theme back
source("helperFunctions.R")


## ----readData, results='markdown', eval=TRUE----------------------------------
loyn <- read_csv('../data/loyn.csv', trim_ws=TRUE)


## -----------------------------------------------------------------------------
#| label: examinData
#| dependson: readData

loyn |> glimpse()


## -----------------------------------------------------------------------------
#| label: strData
#| dependson: readData
loyn |> str()


## -----------------------------------------------------------------------------
#| label: headData
#| dependson: readData
## Explore the first 6 rows of the data
loyn |> head()


## -----------------------------------------------------------------------------
#| label: easyData
#| dependson: readData
loyn |> datawizard::data_codebook()


## -----------------------------------------------------------------------------
loyn |> modelsummary::datasummary_skim(categorical = TRUE)


## ----processData, results='markdown', eval=TRUE-------------------------------
loyn <- loyn |> mutate(fGRAZE = factor(GRAZE))


## ----EDA1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=7, fig.height=7----
scatterplotMatrix(~ABUND+DIST+LDIST+AREA+GRAZE+ALT+YR.ISOL, data = loyn,
                  diagonal = list(method = 'boxplot'))


## ----EDA1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=7, fig.height=7----
scatterplotMatrix(~ABUND+log(DIST)+log(LDIST)+log(AREA)+GRAZE+ALT+YR.ISOL, data = loyn,
                  diagonal = list(method = 'boxplot'))


## ----lm, results='markdown', eval=TRUE, mhidden=TRUE--------------------------
loyn_glm <- glm(ABUND~scale(log(DIST))+scale(log(LDIST))+scale(log(AREA))+
                    fGRAZE + scale(ALT) + scale(YR.ISOL),
                data = loyn,
                family = gaussian(link='log'))
loyn_glm |> summary()


## ----fitModel1a, results='markdown', eval=TRUE, mhidden=TRUE------------------
loyn_rstanarm = stan_glm(ABUND ~ scale(log(DIST))+
                              scale(log(LDIST))+
                              scale(log(AREA))+
                              fGRAZE+
                              scale(ALT)+
                              scale(YR.ISOL),
                         data=loyn,
                         family=gaussian(link='log'),
                         iter = 5000, warmup = 2500,
                         chains = 3, thin = 5, refresh = 0)


## ----fitModel1a1, results='markdown', eval=TRUE, mhidden=TRUE-----------------
loyn_rstanarm = stan_glm(ABUND ~ scale(log(DIST))+
                              scale(log(LDIST))+
                              scale(log(AREA))+
                              fGRAZE+
                              scale(ALT)+
                              scale(YR.ISOL),
                         data=loyn,
                         family=gaussian(link='log'),
                         iter = 5000, warmup = 2500,
                         chains = 3, thin = 5, refresh = 0,
                         adapt_delta = 0.99)


## ----fitModel1b, results='markdown', eval=TRUE, mhidden=TRUE------------------
prior_summary(loyn_rstanarm)


## ----fitModel1f, results='markdown', eval=TRUE, mhidden=TRUE------------------
loyn_rstanarm1 <- update(loyn_rstanarm,  prior_PD=TRUE)

## ----fitModel1g, results='markdown', eval=TRUE, mhidden=TRUE------------------
ggemmeans(loyn_rstanarm1,  ~AREA) |> plot(show_data=TRUE) + scale_y_log10()
ggpredict(loyn_rstanarm1) |>
    plot(show_data = TRUE) |>
    wrap_plots() &
    scale_y_log10()


## ----fitModel1h, results='markdown', eval=TRUE, mhidden=TRUE------------------
loyn_rstanarm2 <- stan_glm(ABUND ~ scale(log(DIST))+
                             scale(log(LDIST))+
                              scale(log(AREA))+
                              fGRAZE+
                              scale(ALT)+
                              scale(YR.ISOL), data=loyn,
                          family=gaussian(link='log'),
                          prior_intercept = normal(3, 1,  autoscale=FALSE),
                          prior = normal(0, 1, autoscale=FALSE),
                          prior_aux = cauchy(0, 2),
                          prior_PD=TRUE,
                          iter = 5000, thin=5,
                          chains = 3, warmup=2500,
                          refresh=0)


## ----fitModel1i, results='markdown', eval=TRUE, mhidden=TRUE------------------
ggemmeans(loyn_rstanarm2,  ~AREA) |>
  plot(show_data=TRUE) + scale_y_log10()
ggpredict(loyn_rstanarm2) |>
    plot(show_data = TRUE) |>
    wrap_plots() &
    scale_y_log10()


## ----fitModel1j, results='markdown', eval=TRUE, mhidden=TRUE, dependson='fitModel1h'----
loyn_rstanarm3= update(loyn_rstanarm2,  prior_PD=FALSE)


## ----modelFit1k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
posterior_vs_prior(loyn_rstanarm3, color_by='vs', group_by=TRUE,
                   facet_args=list(scales='free_y'))


## ----modelFit1l, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#ggemmeans(loyn_rstanarm3,  ~AREA) |> plot(show_data=TRUE) + scale_y_log10()
ggpredict(loyn_rstanarm3,  terms = "AREA[0:1000]") |>
  plot(jitter = FALSE, show_data=TRUE) +
  scale_x_log10() +
  scale_x_log10()
## ggpredict(loyn_rstanarm3,  terms = "AREA[0:1000]") |>
##     plot(jitter = FALSE, residuals=TRUE, log.y = TRUE) + scale_x_log10()
ggpredict(loyn_rstanarm3) |>
    plot(show_data = TRUE) |>
    wrap_plots() &
    scale_y_log10()


## ----fitModel2a, results='markdown', eval=TRUE, mhidden=TRUE------------------
loyn_form <- bf(ABUND ~ scale(log(DIST))+
                       scale(log(LDIST))+
                       scale(log(AREA))+
                       fGRAZE+
                       scale(ALT)+
                       scale(YR.ISOL),
                   family = gaussian(link = 'log'))
loyn_brm <- brm(loyn_form,
                data = loyn,
                iter = 5000,
                warmup = 2500,
                chains = 3, cores = 3,
                thin = 5,
                refresh = 0,
                backend = "cmdstanr")


## ----fitModel2b, results='markdown', eval=TRUE, mhidden=TRUE, paged.print=FALSE,tidy.opts = list(width.cutoff = 80), echo=2----
options(width=100)
loyn_brm |> prior_summary()
options(width=80)


## ----fitModel2d, results='markdown', eval=TRUE, mhidden=TRUE------------------
loyn |> group_by(fGRAZE) |>
  summarise(median(log(ABUND)),
            mad(log(ABUND))
            )
priors <- prior(normal(0, 2.5), class = 'b')
loyn_form <- bf(ABUND ~ scale(log(DIST))+
                    scale(log(LDIST))+
                    scale(log(AREA))+
                    fGRAZE+
                    scale(ALT)+
                    scale(YR.ISOL),
                family = gaussian(link = 'log'))
loyn_brm1 <- brm(loyn_form,
                 data = loyn,
                 prior = priors,
                 sample_prior = 'only',
                 iter = 5000,
                 warmup = 2500,
                 chains = 3,
                 thin = 5,
                 refresh = 0,
                 backend = "cmdstanr")


## ----fitModel2e, results='markdown', eval=TRUE, mhidden=TRUE------------------
## Individual plots - the following seems to be broken??
##loyn_brm1 |> ggemmeans(~AREA) |> plot(show_data = TRUE) + scale_y_log10()
loyn_brm1 |>
    ggemmeans(~AREA) |>
    plot(show_data = TRUE) + scale_y_log10()

## ----fitModel2e2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width = 8, fig.height = 8----
## All effects
loyn_brm1 |>
    conditional_effects() |>
    plot(points = TRUE, ask = FALSE, plot = FALSE) |>
    wrap_plots() &
    scale_y_log10()

loyn_brm1 |>
    ggpredict() |>
    plot(show_data=TRUE, facet=TRUE) +
    scale_y_log10()


## ----fitModel2h, results='markdown', eval=TRUE, mhidden=TRUE------------------
mod.mat <- model.matrix(as.formula(loyn_form), data = loyn)
mad(log(loyn$ABUND))/
    apply(mod.mat, 2, mad)

loyn |> group_by(fGRAZE) |>
  summarise(median = median(log(ABUND)),
            mad = mad(log(ABUND))
            )
priors <- prior(normal(3.4, 0.1),  class = 'Intercept') +
    prior(normal(0, 2), class = 'b') +
    prior(student_t(3, 0, 1.5), class = 'sigma')
loyn_form <- bf(ABUND ~ scale(log(DIST))+
                     scale(log(LDIST))+
                     scale(log(AREA))+
                     fGRAZE+
                     scale(ALT)+
                     scale(YR.ISOL),
                   family = gaussian(link = 'log'))
loyn_brm2 <- brm(loyn_form,
                data = loyn,
                prior = priors,
                sample_prior = 'only',
                iter = 5000,
                warmup = 2500,
                chains = 3, cores = 3,
                thin = 5,
                refresh = 0,
                backend = "cmdstanr")


## ----fitModel2i, results='markdown', eval=TRUE, mhidden=TRUE------------------
loyn_brm2 |> ggemmeans(~DIST) |>
    plot(show_data = TRUE) +
    scale_y_log10()

## ----fitModel2i2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width = 8, fig.height = 8----
loyn_brm2 |>
    conditional_effects() |>
    plot(points = TRUE, ask = FALSE, plot = FALSE) |>
    ## lapply(function(x) x + scale_y_log10()) |>
    wrap_plots() &
    scale_y_log10()

loyn_brm2 |>
    ggpredict() |>
    plot(show_data=TRUE, facet=TRUE) +
    scale_y_log10()


## ----fitModel2j, results='markdown', eval=TRUE, mhidden=TRUE------------------
loyn_brm3 <- update(loyn_brm2,  sample_prior = 'yes', refresh = 0)


## ----fitModel2j2, results='markdown', eval=TRUE, echo = FALSE, mhidden=TRUE----
save(loyn_brm3, file = '../ws/testing/loyn_brm3.RData')


## ----fitModel2k, results='markdown', eval=TRUE, mhidden=TRUE------------------
loyn_brm3 |> get_variables()
## loyn_brm3 |> hypothesis('Intercept = 0', class = 'b') |> plot
## loyn_brm3 |> hypothesis('Intercept = 0', class = 'prior') |> plot
loyn_brm3 |> hypothesis('scalelogDIST = 0') |> plot()
loyn_brm3 |> hypothesis('scalelogAREA = 0') |> plot()
loyn_brm3 |> hypothesis('sigma = 0', class = '') |> plot()


## ----fitModel2k2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width = 8, fig.height = 4----
loyn_brm3 |> SUYR_prior_and_posterior()


## ----fitModel2l, results='markdown', eval=TRUE, mhidden=TRUE------------------
loyn_brm3 |> standata()
loyn_brm3 |> stancode()


## ----modelValidation1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
available_mcmc()


## ----modelValidation1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
plot(loyn_rstanarm3, plotfun='mcmc_trace')


## ----modelValidation1c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
plot(loyn_rstanarm3, 'acf_bar')


## ----modelValidation1d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(loyn_rstanarm3, 'rhat_hist')


## ----modelValidation1e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(loyn_rstanarm3, 'neff_hist')


## ----Validation1f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
plot(loyn_rstanarm3, 'combo')
plot(loyn_rstanarm3, 'violin')


## ----modelValidation1g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_trace(loyn_rstanarm3)


## ----modelValidation1h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_ac(loyn_rstanarm3)


## ----modelValidation1i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_rhat(loyn_rstanarm3)


## ----modelValidation1j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_ess(loyn_rstanarm3)


## ----modelValidation1k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_dens(loyn_rstanarm3, separate_chains = TRUE)


## ----modelValidation1l, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=12----
loyn_ggs <- ggs(loyn_rstanarm3)
ggs_traceplot(loyn_ggs)


## ----modelValidation1m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=12----
ggs_autocorrelation(loyn_ggs)


## ----modelValidation1n, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_Rhat(loyn_ggs)


## ----modelValidation1o, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_effective(loyn_ggs)


## ----modelValidation1p, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_crosscorrelation(loyn_ggs)


## ----modelValidation1q, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_grb(loyn_ggs)


## ----modelValidation2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
available_mcmc()


## ----modelValidation2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
loyn_brm3 |> mcmc_plot(type = 'trace')


## ----modelValidation2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
loyn_brm3 |> mcmc_plot(type = 'acf_bar')


## ----modelValidation2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
loyn_brm3 |> mcmc_plot(type = 'rhat_hist')


## ----modelValidation2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
loyn_brm3 |> mcmc_plot(type = 'neff_hist')


## ----Validation2f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
loyn_brm3 |> mcmc_plot(type = 'combo')
loyn_brm3 |> mcmc_plot(type = 'violin')


## ----modelValidation2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
loyn_brm3$fit |> stan_trace()


## ----modelValidation2h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
loyn_brm3$fit |> stan_ac()


## ----modelValidation2i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
loyn_brm3$fit |> stan_rhat()


## ----modelValidation2j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
loyn_brm3$fit |> stan_ess()


## ----modelValidation2k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
loyn_brm3$fit |> stan_dens(separate_chains = TRUE)


## ----modelValidation2l, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=12----
loyn_ggs <- loyn_brm3 |> ggs(inc_warmup = FALSE, burnin = FALSE)
loyn_ggs |> ggs_traceplot()


## ----modelValidation2m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=12----
loyn_ggs |> ggs_autocorrelation()


## ----modelValidation2n, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
loyn_ggs |> ggs_Rhat()


## ----modelValidation2o, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
loyn_ggs |> ggs_effective()


## ----modelValidation2p, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
loyn_ggs |> ggs_crosscorrelation()


## ----modelValidation2q, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
loyn_ggs |> ggs_grb()


## ----modelValidation3a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation3b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(loyn_rstanarm3,  plotfun='dens_overlay')


## ----modelValidation3c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#pp_check(loyn_rstanarm3, plotfun='error_scatter_avg')


## ----modelValidation3d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(loyn_rstanarm3, x=loyn$AREA, plotfun='error_scatter_avg_vs_x')


## ----modelValidation3e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(loyn_rstanarm3, x=loyn$AREA, plotfun='intervals')


## ----modelValidation3f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(loyn_rstanarm3, x=loyn$AREA, plotfun='ribbon')


## ----modelValidation3g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(loyn_rstanarm3)


## ----modelValidation4a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
preds <- posterior_predict(loyn_rstanarm3,  ndraws=250,  summary=FALSE)
loyn_resids <- createDHARMa(simulatedResponse = t(preds),
                            observedResponse = loyn$ABUND,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse = FALSE)
plot(loyn_resids)


## ----modelValidation5a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation5b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
loyn_brm3 |> pp_check(type = 'dens_overlay', ndraws = 100)


## ----modelValidation5c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#loyn_brm3 |> pp_check(type = 'error_scatter_avg')


## ----modelValidation5d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
loyn_brm3 |> pp_check(x = 'AREA', type = 'error_scatter_avg_vs_x')


## ----modelValidation5e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
loyn_brm3 |> pp_check(x = 'AREA', type = 'intervals')


## ----modelValidation5f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
loyn_brm3 |> pp_check(x = 'AREA', type = 'ribbon')


## ----modelValidation5g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(loyn_brm3)


## ----modelValidation6a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
preds <- loyn_brm3 |> posterior_predict(ndraws = 250,  summary = FALSE)
loyn_resids <- createDHARMa(simulatedResponse = t(preds),
                            observedResponse = loyn$ABUND,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse = FALSE)
loyn_resids |> plot()


## ----modelValidation6aa, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=10----
loyn_resids <- make_brms_dharma_res(loyn_brm3, integerResponse = FALSE)
wrap_elements(~testUniformity(loyn_resids)) +
               wrap_elements(~plotResiduals(loyn_resids, form = factor(rep(1, nrow(loyn))))) +
               wrap_elements(~plotResiduals(loyn_resids, quantreg = TRUE)) +
               wrap_elements(~testDispersion(loyn_resids))



## ----partialPlot1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=7, fig.height=7----
loyn_rstanarm3 |> ggpredict() |> plot(show_data=TRUE, facet=TRUE)


## ----partialPlot1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
#loyn_rstanarm3 |> ggemmeans(~AREA,  type='fixed') |> plot(show_data=TRUE) + scale_y_log10()


## ----partialPlot1c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
loyn_rstanarm3 |> fitted_draws(newdata=loyn) |>
  median_hdci() |>
  ggplot(aes(x=AREA, y=.value)) +
  geom_ribbon(aes(ymin=.lower, ymax=.upper), fill='blue', alpha=0.3) +
  geom_line() +
  geom_point(data=loyn,  aes(y=ABUND,  x=AREA)) +
  scale_y_log10()


## ----partialPlot2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
loyn_brm3 |>
  conditional_effects() |>
  plot(ask = FALSE, points = TRUE, plot = FALSE) |>
  wrap_plots()

loyn_brm3 |>
    conditional_effects() |>
    plot(ask = FALSE, points = TRUE, plot = FALSE) |>
    wrap_plots() &
    scale_y_log10()

g <- loyn_brm3 |>
  conditional_effects() |>
  plot(ask = FALSE, points = TRUE, plot = FALSE)
library(patchwork)
length(g)
(g[[1]] + scale_x_log10()) +
    (g[[2]] + scale_x_log10()) +
    (g[[3]] + scale_x_log10()) +
    g[[4]] +
    g[[5]] +
    g[[6]]


## ----partialPlot2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
loyn_brm3 |>
    ggpredict() |>
    plot(show_data=TRUE) |>
    wrap_plots()
loyn_brm3 |>
    ggpredict() |>
    plot(show_data=TRUE) |>
    wrap_plots() &
    scale_y_log10()


## ----partialPlot2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8, cache = TRUE----
loyn_brm3 |>
    ggemmeans("AREA[0:1000]") |>
    plot(show_data=TRUE) |>
    wrap_plots() &
    scale_y_log10() &
    scale_x_log10()


## ----partialPlot2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
loyn_brm3 |>
    epred_draws(newdata = loyn) |>
    median_hdci(.epred) |>
    ggplot(aes(x = AREA, y = .epred, colour = fGRAZE, fill = fGRAZE)) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), colour = NA, alpha = 0.3) +
    geom_line() +
    geom_point(data = loyn,  aes(y = ABUND,  x = AREA)) +
    scale_y_log10() +
    scale_x_log10()


## ----summariseModel1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
summary(loyn_rstanarm3)


## ----summariseModel1a1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=FALSE----
loyn_sum <- summary(loyn_rstanarm3)


## ----summariseModel1dd, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
loyn_rstanarm3$stanfit |>
    summarise_draws(median,
                    HDInterval::hdi,
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel1d2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
loyn_rstanarm3$stanfit |>
    summarise_draws(median,
                    ~HDInterval::hdi(.x, credMass = 0.9),
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel1d3, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
loyn_rstanarm3$stanfit |>
    summarise_draws(
        ~ median(exp(.x)),
        ~HDInterval::hdi(exp(.x)),
        rhat, length, ess_bulk, ess_tail)


## ----summariseModel1e1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_rstanarm3 |> tidy_draws()


## ----summariseModel1i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
loyn_rstanarm3$stanfit |> as_draws_df()
loyn_rstanarm3$stanfit |>
  as_draws_df() |>
  summarise_draws(
    "median",
    ~ HDInterval::hdi(.x),
    "rhat",
    "ess_bulk"
  )


## ----summariseModel1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tidyMCMC(loyn_rstanarm3$stanfit, estimate.method='median',  conf.int=TRUE,
         conf.method='HPDinterval',  rhat=TRUE, ess=TRUE)

## ----summariseModel1b1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
loyn_tidy <- tidyMCMC(loyn_rstanarm3$stanfit, estimate.method='median',  conf.int=TRUE,  conf.method='HPDinterval',  rhat=TRUE, ess=TRUE)


## ----summariseModel1c1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_rstanarm3 |> get_variables()


## ----summariseModel1c2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_draw <- loyn_rstanarm3 |> gather_draws(`.Intercept.*|.*AREA.*|.*DIST.*|.*GRAZE.*|.*ALT.*|.*YR.*`,  regex=TRUE)
loyn_draw


## ----summariseModel1d1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_rstanarm3 |> plot(plotfun='mcmc_intervals')


## ----summariseModel1c4, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_rstanarm3 |>
  gather_draws(`.*AREA.*|.*DIST.*|.*GRAZE.*|.*ALT.*|.*YR.*`, regex=TRUE) |>
  ggplot() +
  stat_halfeye(aes(x=.value,  y=.variable)) +
  facet_wrap(~.variable, scales='free')

loyn_rstanarm3 |>
  gather_draws(`.*AREA.*|.*DIST.*|.*GRAZE.*|.*ALT.*|.*YR.*`, regex=TRUE) |>
  ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    geom_vline(xintercept = 1, linetype = 'dashed')


## ----summariseModel1c7, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_rstanarm3 |>
  gather_draws(`.*AREA.*|.*DIST.*|.*GRAZE.*|.*ALT.*|.*YR.*`, regex=TRUE) |>
  ggplot() +
    geom_density_ridges(aes(x=.value, y = .variable), alpha=0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed')
##Or on a fractional scale
loyn_rstanarm3 |>
  gather_draws(`.*AREA.*|.*DIST.*|.*GRAZE.*|.*ALT.*|.*YR.*`, regex=TRUE) |>
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


## ----summariseModel1f1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_rstanarm3 |> spread_draws(`.Intercept.*|.*DIST.*|.*AREA.*|.*GRAZE.*|.*ALT.*|.*YR.*`,  regex=TRUE)


## ----summariseModel1g1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_rstanarm3 |> posterior_samples() |> as.tibble()


## ----summariseModel1h1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_rstanarm3 |> bayes_R2() |> median_hdci()


## ----summariseModel2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
loyn_brm3 |> summary()


## ----summariseModel2a1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=FALSE----
loyn_sum <- summary(loyn_brm3)


## ----summariseModel2e1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_brm3 |> tidy_draws()
loyn_brm3 |>
    tidy_draws() |>
    dplyr::select(starts_with("b_")) |>
    exp() |>
    summarise_draws(median,
                    HDInterval::hdi,
                    rhat,
                    length,
                    ess_bulk, ess_tail)


loyn_brm3 |>
    spread_draws(`^b_.*|sigma`, regex = TRUE) |>
    exp() |>
    summarise_draws(median,
                    HDInterval::hdi,
                    rhat,
                    length,
                    ess_bulk, ess_tail)

## we can also attempt to provide the slopes back on the scale of the unscaled predictors
loyn_brm3 |>
  emtrends(~1, var = "log(AREA)") |>
  tidy_draws() |>
  mutate(.value = exp(`1 overall`)) |>
  summarise_draws(median,
                  HDInterval::hdi)


## ----summariseModel2i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
loyn_brm3 |> as_draws_df()
loyn_brm3 |>
  as_draws_df() |>
  dplyr::select(matches("^b_.*|^sigma$")) |>
  mutate(across(everything(), exp)) |>
  summarise_draws(
    "median",
    ~ HDInterval::hdi(.x),
    "rhat",
    "ess_bulk"
  )



## ----summariseModel2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
loyn_brm3$fit |>
    tidyMCMC(estimate.method = 'median',
             conf.int = TRUE,  conf.method = 'HPDinterval',
             rhat = TRUE, ess = TRUE)

## ----summariseModel2b1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
loyn_tidy <- tidyMCMC(loyn_brm3$fit, estimate.method='median',  conf.int=TRUE,  conf.method='HPDinterval',  rhat=TRUE, ess=TRUE)


## ----summariseModel2c1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_brm3 |> get_variables()


## ----summariseModel2c2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_draw <- loyn_brm3 |> gather_draws(`^b_.*`,  regex = TRUE)
loyn_draw


## ----summariseModel2d1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_brm3 |> mcmc_plot(type = 'intervals')


## ----summariseModel2c4, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_brm3 |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(.variable != 'b_Intercept') |>
    ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    facet_wrap(~.variable, scales='free')

loyn_brm3 |>
  gather_draws(`^b_.*`, regex=TRUE) |>
  mutate(.value = exp(.value)) |>
  filter(.variable != 'b_Intercept') |>
  ggplot() +
  stat_halfeye(aes(x=.value,  y=.variable)) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  scale_x_continuous(trans = scales::log2_trans()) +
  theme_classic()


## ----summariseModel2c7, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_brm3 |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(.variable != 'b_Intercept') |>
    ggplot() +
    geom_density_ridges(aes(x=.value, y = .variable), alpha=0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed')
##Or on a fractional scale
loyn_brm3 |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(.variable != 'b_Intercept') |>
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


## ----summariseModel2f1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_brm3 |> spread_draws(`^b_.*`,  regex = TRUE)


## ----summariseModel2g1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_brm3 |> posterior_samples() |> as_tibble()


## ----summariseModel2h1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_brm3 |> bayes_R2(summary = FALSE) |> median_hdci()


## ----ROP, results='markdown', eval=TRUE---------------------------------------
0.1 * sd(log(loyn$ABUND))
loyn_brm3 |> bayestestR::rope_range()
loyn_brm3 |> bayestestR::rope(range = c(-0.09, 0.09))
loyn_brm3 |> bayestestR::rope(range = c(-0.09, 0.09)) |>
    plot(data = loyn_brm3)


## ----furtherModel1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_rstanarm4a <- update(loyn_rstanarm3,  .~scale(log(DIST))*scale(log(LDIST)),
                          diagnostic_file = file.path(tempdir(), "dfa.csv"))
loyn_rstanarm4b <- update(loyn_rstanarm3,  .~scale(log(AREA)) * fGRAZE,
                          diagnostic_file = file.path(tempdir(), "dfb.csv"))
loyn_rstanarm4c <- update(loyn_rstanarm3,  .~scale(log(AREA)) * fGRAZE * scale(YR.ISOL),
                          diagnostic_file = file.path(tempdir(), "dfc.csv"))
loyn_rstanarm4d <- update(loyn_rstanarm3,  .~scale(ALT),
                          diagnostic_file = file.path(tempdir(), "dfd.csv"))
loyn_rstanarm4e <- update(loyn_rstanarm3,  .~1,
                          diagnostic_file = file.path(tempdir(), "dfe.csv"))
loo_compare(loo(loyn_rstanarm4a),
            loo(loyn_rstanarm4e)
            )
loo_compare(loo(loyn_rstanarm4b),
            loo(loyn_rstanarm4e)
            )
loo_compare(loo(loyn_rstanarm4c),
            loo(loyn_rstanarm4e)
            )
loo_compare(loo(loyn_rstanarm4d),
            loo(loyn_rstanarm4e)
            )


## ----furtherModel1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
bayes_factor(bridge_sampler(loyn_rstanarm4a),
             bridge_sampler(loyn_rstanarm4e))
bayes_factor(bridge_sampler(loyn_rstanarm4b),
             bridge_sampler(loyn_rstanarm4e))
bayes_factor(bridge_sampler(loyn_rstanarm4c),
             bridge_sampler(loyn_rstanarm4e))
bayes_factor(bridge_sampler(loyn_rstanarm4d),
             bridge_sampler(loyn_rstanarm4e))


## ----furtherModel2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_brm4a <- update(loyn_brm3,  .~scale(log(DIST))*scale(log(LDIST)),
                     save_pars = save_pars(all = TRUE), refresh = 0)
loyn_brm4b <- update(loyn_brm3,  .~scale(log(AREA)) * fGRAZE,
                     save_pars = save_pars(all = TRUE), refresh = 0)
loyn_brm4c <- update(loyn_brm3,  .~scale(log(AREA)) * fGRAZE * scale(YR.ISOL),
                     save_pars = save_pars(all = TRUE), refresh = 0)
loyn_brm4d <- update(loyn_brm3,  .~scale(ALT),
                     save_pars = save_pars(all = TRUE), refresh = 0)
loyn_brm4e <- update(loyn_brm3,  .~1,
                     save_pars = save_pars(all = TRUE), refresh = 0)
waic(loyn_brm4a)
loo(loyn_brm4a)
loo(loyn_brm4e)
loo_compare(loo(loyn_brm4a),
            loo(loyn_brm4e)
            )
loo_compare(loo(loyn_brm4b),
            loo(loyn_brm4e)
            )
loo_compare(loo(loyn_brm4b, moment_match = TRUE),
            loo(loyn_brm4e)
            )
loo_compare(loo(loyn_brm4c, moment_match = TRUE),
            loo(loyn_brm4e)
            )
loo_compare(loo(loyn_brm4d),
            loo(loyn_brm4e)
            )


## ----furtherModel2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE, cache = TRUE----
bayes_factor(loyn_brm4a,
            loyn_brm4e)
#OR
bayes_factor(loyn_brm4e,
            loyn_brm4a)

bayes_factor(loyn_brm4b,
             loyn_brm4e)
#OR
bayes_factor(loyn_brm4e,
             loyn_brm4b)

bayes_factor(loyn_brm4c,
             loyn_brm4e)
#OR
bayes_factor(loyn_brm4e,
             loyn_brm4c)

bayes_factor(loyn_brm4d,
             loyn_brm4e)
#OR
bayes_factor(loyn_brm4e,
             loyn_brm4d)


## ----furtherModel2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
loyn_list <- with(loyn, list(AREA = c(min(AREA), mean(AREA), max(AREA))))

loyn_brm4b |>
    emmeans(~fGRAZE|AREA, at = loyn_list, type = "response") |>
    pairs(reverse = FALSE)

newdata <- loyn_brm4b |>
    emmeans(~fGRAZE|AREA, at = loyn_list, type = 'response') |>
    pairs() |>
    as.data.frame()

head(newdata)


newdata <- loyn_brm4b |>
    emmeans(~fGRAZE|AREA, at = loyn_list, type = 'response') |>
    pairs() |>
    gather_emmeans_draws()

newdata |> median_hdci() |>
    ggplot() +
    geom_hline(yintercept = 1, linetype = 'dashed') +
    geom_pointrange(aes(y = .value, ymin = .lower, ymax = .upper, x = contrast)) +
    facet_wrap(~AREA) +
    coord_flip()

loyn_brm4b |>
    emmeans(~fGRAZE|AREA, at = loyn_list, type = 'response') |>
    gather_emmeans_draws()
newdata.p <- newdata |> summarise(P = mean(.value>1))
g <- newdata |>
    ggplot() +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    stat_slab(aes(x  =  .value, y = contrast,
                  fill = stat(ggdist::cut_cdf_qi(cdf,
                           .width = c(0.5, 0.8, 0.95),
                           labels = scales::percent_format())
                           )), color = 'black') +
    scale_fill_brewer('Interval', direction = -1, na.translate = FALSE) +
    facet_grid(~round(AREA,1))


## ----furtherModel2c1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
g + geom_text(data = newdata.p, aes(y = contrast, x = 1, label = paste('P = ',round(P,3))), hjust = -0.2, position = position_nudge(y = 0.5))


## ----summaryFigure1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=TRUE----
loyn_list <- with(loyn, list(AREA = modelr::seq_range(AREA, n=100)))

newdata <- emmeans(loyn_brm3, ~AREA|fGRAZE, at = loyn_list, type='response') |>
    as.data.frame()
head(newdata)

ggplot(newdata, aes(y=response, x=AREA)) +
  geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD, fill=fGRAZE), alpha=0.3) +
  geom_line(aes(color=fGRAZE)) +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10()

spaghetti = emmeans(loyn_rstanarm3, ~AREA|fGRAZE, at = loyn_list, type='response') |>
  gather_emmeans_draws() |>
  mutate(Fit=exp(.value))
wch = sample(1:max(spaghetti$.draw), 100,replace=FALSE)
spaghetti = spaghetti |>
  filter(.draw %in% wch)
ggplot(newdata) +
  geom_line(data=spaghetti, aes(y=Fit, x=AREA, color=fGRAZE,
                                group=interaction(fGRAZE,.draw)), alpha=0.05) +
  geom_line(aes(y=response, x=AREA, color=fGRAZE)) +
  theme_bw() +
  scale_x_log10() + scale_y_log10()

## or honouring the data range
loyn_nd <- loyn |>
    group_by(fGRAZE) |>
    tidyr::expand(AREA = modelr::seq_range(AREA, n=100),
           DIST = mean(DIST), LDIST = mean(LDIST), ALT = mean(ALT), YR.ISOL = mean(YR.ISOL))
loyn_rstanarm3 |>
    epred_draws(newdata = loyn_nd, value = '.value') |>
    median_hdci() |>
    ggplot(aes(x = AREA, y = .value, colour = fGRAZE, fill = fGRAZE)) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), colour = NA, alpha = 0.3) +
    geom_line() +
    geom_point(data = loyn, aes(y = ABUND, x = AREA)) +
    scale_y_log10() +
    scale_x_log10() +
    theme_bw()



## ----summaryFigure1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=TRUE----
loyn_list = with(loyn, list(AREA = seq(min(AREA), max(AREA), len=100)))

newdata = emmeans(loyn_rstanarm4b, ~AREA|fGRAZE, at = loyn_list, type='response') |>
    as.data.frame()
head(newdata)

ggplot(newdata, aes(y=response, x=AREA)) +
  geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD, fill=fGRAZE), alpha=0.3) +
  geom_line(aes(color=fGRAZE)) +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10()

spaghetti = emmeans(loyn_rstanarm4b, ~AREA|fGRAZE, at = loyn_list, type='response') |>
  gather_emmeans_draws() |> mutate(Fit=exp(.value))
wch = sample(1:max(spaghetti$.draw), 100,replace=FALSE)
spaghetti = spaghetti |> filter(.draw %in% wch)
ggplot(newdata) +
  geom_line(data=spaghetti, aes(y=Fit, x=AREA, color=fGRAZE,
                                group=interaction(fGRAZE,.draw)), alpha=0.05) +
  geom_line(aes(y=response, x=AREA, color=fGRAZE)) +
  theme_bw() +
  scale_x_log10() + scale_y_log10()

## or honouring the data range
loyn_nd <- loyn |>
    group_by(fGRAZE) |>
    tidyr::expand(AREA = modelr::seq_range(AREA, n=100))
loyn_rstanarm4b |>
    epred_draws(newdata = loyn_nd, value = '.value') |>
    median_hdci() |>
    ggplot(aes(x = AREA, y = .value, colour = fGRAZE, fill = fGRAZE)) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), colour = NA, alpha = 0.3) +
    geom_line() +
    geom_point(data = loyn, aes(y = ABUND, x = AREA)) +
    scale_y_log10() +
    scale_x_log10() +
    theme_bw()


## ----summaryFigure2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=TRUE----
loyn_list <- with(loyn, list(AREA = modelr::seq_range(AREA, n=100)))

newdata <- emmeans(loyn_brm3, ~AREA|fGRAZE, at = loyn_list, type='response') |>
    as.data.frame()
head(newdata)

ggplot(newdata, aes(y = response, x = AREA)) +
  geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD, fill = fGRAZE), alpha = 0.3) +
  geom_line(aes(color = fGRAZE)) +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10()

spaghetti = emmeans(loyn_brm3, ~AREA|fGRAZE, at = loyn_list, type = 'response') |>
  gather_emmeans_draws() |> mutate(Fit = exp(.value))
wch <- sample(1:max(spaghetti$.draw), 100,replace = FALSE)
spaghetti <- spaghetti |> filter(.draw %in% wch)
ggplot(newdata) +
  geom_line(data = spaghetti, aes(y = Fit, x = AREA, color = fGRAZE,
                                group = interaction(fGRAZE,.draw)), alpha = 0.1) +
  geom_line(aes(y = response, x = AREA, color = fGRAZE)) +
  theme_bw() +
  scale_x_log10() + scale_y_log10()

# or honouring the data range
loyn_nd <- loyn |>
    group_by(fGRAZE) |>
    tidyr::expand(AREA = modelr::seq_range(AREA, n=100),
           DIST = mean(DIST), LDIST = mean(LDIST), ALT = mean(ALT), YR.ISOL = mean(YR.ISOL))
loyn_brm3 |>
    epred_draws(newdata = loyn_nd, value = '.value') |>
    median_hdci() |>
    ggplot(aes(x = AREA, y = .value, colour = fGRAZE, fill = fGRAZE)) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), colour = NA, alpha = 0.3) +
    geom_line() +
    geom_point(data = loyn, aes(y = ABUND, x = AREA)) +
    scale_y_log10() +
    scale_x_log10() +
    theme_bw()



## ----summaryFigure2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=TRUE----
loyn_list <- with(loyn, list(AREA = modelr::seq_range(AREA, n = 100)))

newdata <- emmeans(loyn_brm4b, ~AREA|fGRAZE, at = loyn_list, type='response') |>
    as.data.frame()
head(newdata)

ggplot(newdata, aes(y = response, x = AREA)) +
  geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD, fill = fGRAZE), alpha = 0.3) +
  geom_line(aes(color = fGRAZE)) +
  theme_bw() +
  scale_x_log10() +
  scale_y_log10()

spaghetti = emmeans(loyn_brm4b, ~AREA|fGRAZE, at = loyn_list, type = 'response') |>
  gather_emmeans_draws() |> mutate(Fit = exp(.value))
wch <- sample(1:max(spaghetti$.draw), 100,replace = FALSE)
spaghetti <- spaghetti |> filter(.draw %in% wch)
ggplot(newdata) +
  geom_line(data = spaghetti, aes(y = Fit, x = AREA, color = fGRAZE,
                                group = interaction(fGRAZE,.draw)), alpha = 0.1) +
  geom_line(aes(y = response, x = AREA, color = fGRAZE)) +
  theme_bw() +
  scale_x_log10() + scale_y_log10()

## or honouring the data range
loyn_nd <- loyn |>
    group_by(fGRAZE) |>
    tidyr::expand(AREA = modelr::seq_range(AREA, n=100))
loyn_brm4b |>
    epred_draws(newdata = loyn_nd, value = '.value') |>
    median_hdci() |>
    ggplot(aes(x = AREA, y = .value, colour = fGRAZE, fill = fGRAZE)) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), colour = NA, alpha = 0.3) +
    geom_line() +
    geom_point(data = loyn, aes(y = ABUND, x = AREA)) +
    scale_y_log10() +
    scale_x_log10() +
    theme_bw()

