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
library(ggridges)      #for ridge plots
library(patchwork)     #for multiple plots
library(modelsummary)  #for data and model summaries
theme_set(theme_grey()) #put the default ggplot theme back
source("helperFunctions.R")


## -----------------------------------------------------------------------------
#| label: readData

day <- read_csv('../data/day.csv', trim_ws = TRUE)


## -----------------------------------------------------------------------------
#| label: examinData
glimpse(day)


## -----------------------------------------------------------------------------
## Explore the first 6 rows of the data
head(day)


## -----------------------------------------------------------------------------
str(day)


## -----------------------------------------------------------------------------
day |> datawizard::data_codebook()


## -----------------------------------------------------------------------------
day |> modelsummary::datasummary_skim()
day |> modelsummary::datasummary_skim(by = "TREAT")


## ----prepare, results='markdown', eval=TRUE, mhidden=TRUE---------------------
day <- day |> mutate(TREAT = factor(TREAT))


## ----EDA, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=6, mhidden=TRUE----
day |> ggplot(aes(y = BARNACLE, x = TREAT)) +
    geom_boxplot()+
    geom_point(color = 'red')
day |> ggplot(aes(y = BARNACLE, x = TREAT)) +
    geom_violin()+
    geom_point(color = 'red')


## ----fitModel1a, results='markdown', eval=TRUE, mhidden=TRUE------------------
day_rstanarm = stan_glm(BARNACLE ~ TREAT, data=day,
                        family=poisson(link='log'),
                         iter = 5000, warmup = 1000,
                         chains = 3, thin = 5, refresh = 0)


## ----fitModel1b, results='markdown', eval=TRUE, mhidden=TRUE------------------
prior_summary(day_rstanarm)


## ----fitModel1d, results='markdown', eval=TRUE, mhidden=TRUE------------------
2.5/sd(model.matrix(~TREAT, day)[, 2])


## ----fitModel1f, results='markdown', eval=TRUE, mhidden=TRUE------------------
day_rstanarm1 <- update(day_rstanarm,  prior_PD=TRUE)


## ----fitModel1g, results='markdown', eval=TRUE, mhidden=TRUE------------------
day_rstanarm1 |>
    ggpredict() |>
    plot(show_data = TRUE)

ggpredict(day_rstanarm1) |>
    plot(show_data=TRUE) |>
    wrap_plots() &
    scale_y_log10()


## ----echo = FALSE, eval = FALSE-----------------------------------------------
# day |>
#     group_by(TREAT) |>
#     summarise(mean(log(BARNACLE)),
#               sd(log(BARNACLE)))
# 
# model.matrix(~TREAT, data = day)
# apply(model.matrix(~TREAT, data = day), 2, sd)
# sd(log(day$BARNACLE))/apply(model.matrix(~TREAT, data = day), 2, sd)
# 
# ## day |>
# ##     group_by(TREAT) |>
# ##     summarise(BARNACLE = sample(log(BARNACLE), 1000, replace = TRUE),
# ##               R = 1:1000) |>
# ##     ungroup() |>
# ##     group_by(R) |>
# ##     summarise(D = diff(BARNACLE), Comp = 1:3) |>
# ##     ungroup() |>
# ##     group_by(Comp) |>
# ##     summarise(sd(D))


## ----fitModel1h, results='markdown', eval=TRUE, mhidden=TRUE------------------
day_rstanarm2 <- stan_glm(BARNACLE ~ TREAT, data=day,
                        family=poisson(link='log'),
                         prior_intercept = normal(3, 1, autoscale=FALSE),
                         prior = normal(0, 1, autoscale=FALSE),
                         prior_PD=TRUE,
                         iter = 5000, warmup = 1000,
                         chains = 3, thin = 5, refresh = 0
                         )


## ----fitModel1i, results='markdown', eval=TRUE, mhidden=TRUE------------------
ggpredict(day_rstanarm2) |>
  plot(show_data=TRUE) +
  scale_y_log10()


## ----fitModel1j, results='markdown', eval=TRUE, mhidden=TRUE------------------
day_rstanarm3= update(day_rstanarm2,  prior_PD=FALSE)


## ----modelFit1k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
posterior_vs_prior(day_rstanarm3, color_by='vs', group_by=TRUE,
                   facet_args=list(scales='free_y'))


## ----modelFit1l, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggpredict(day_rstanarm3) |> plot(show_data=TRUE)


## ----fitModel2a, results='markdown', eval=TRUE, mhidden=TRUE------------------
day_form <- bf(BARNACLE ~ TREAT,
               family = poisson(link = 'log'))
day_brm <- brm(day_form,
               data = day,
               iter = 5000,
               warmup = 2500,
               chains = 3, cores = 3,
               thin = 5,
               refresh = 0,
               backend = "cmdstan")


## ----fitModel2b, results='markdown', eval=TRUE, mhidden=TRUE, paged.print=FALSE,tidy.opts = list(width.cutoff = 80), echo=2----
options(width=100)
day_brm |> prior_summary()
options(width=80)


## ----fitModel2d, results='markdown', eval=TRUE, mhidden=TRUE------------------
model.matrix(~TREAT, data = day)
apply(model.matrix(~TREAT, data = day), 2, sd)
sd(log(day$BARNACLE))/apply(model.matrix(~TREAT, data = day), 2, sd)

day_form <- bf(BARNACLE ~ TREAT,
               family = poisson(link = 'log'))
day |>
    group_by(TREAT) |>
    summarise(median(log(BARNACLE)),
              mad(log(BARNACLE)))

priors <- prior(normal(0, 1),  class = 'b')
day_brm1 <- brm(day_form,
               data = day,
               prior = priors,
               sample_prior = 'only',
               iter = 5000,
               warmup = 2500,
               chains = 3, cores = 3,
               thin = 5,
               refresh = 0,
               backend = "cmdstan")


## ----fitModel2e, results='markdown', eval=TRUE, mhidden=TRUE------------------
day_brm1 |> ggpredict() |> plot(show_data=TRUE)
day_brm1 |>
    ggpredict() |>
    plot(show_data=TRUE) |>
    wrap_plots() &
    scale_y_log10()


## ----fitModel2e2, results='markdown', eval=TRUE, mhidden=TRUE-----------------
day_brm1 |>
    ggemmeans(~TREAT) |>
    plot(show_data=TRUE) +
    scale_y_log10()


## ----fitModel2e3, results='markdown', eval=TRUE, mhidden=TRUE-----------------
day_brm1 |>
    conditional_effects() |>
    plot(points=TRUE)
day_brm1 |>
    conditional_effects() |>
    plot(points=TRUE) |>
    wrap_plots() &
    scale_y_log10()


## ----echo = FALSE, eval = FALSE-----------------------------------------------
# day |>
#     group_by(TREAT) |>
#     summarise(median(log(BARNACLE)),
#               mad(log(BARNACLE)))
# 
# ## model.matrix(~TREAT, data = day)
# ## apply(model.matrix(~TREAT, data = day), 2, sd)
# ## sd(log(day$BARNACLE))/apply(model.matrix(~TREAT, data = day), 2, sd)
# 
# ## day |>
# ##     group_by(TREAT) |>
# ##     summarise(BARNACLE = sample(log(BARNACLE), 1000, replace = TRUE),
# ##               R = 1:1000) |>
# ##     ungroup() |>
# ##     group_by(R) |>
# ##     summarise(D = diff(BARNACLE), Comp = 1:3) |>
# ##     ungroup() |>
# ##     group_by(Comp) |>
# ##     summarise(sd(D))


## ----fitModel2h, results='markdown', eval=TRUE, mhidden=TRUE------------------
priors <- prior(normal(3.1, 0.3),  class = 'Intercept') +
    prior(normal(0, 1), class = 'b')
day_form <- bf(BARNACLE ~ TREAT,
                   family = poisson(link = 'log'))
day_brm2 <- brm(day_form,
                data = day,
                prior = priors,
                sample_prior = 'only',
                iter = 5000,
                warmup = 1000,
                chains = 3, cores = 3,
                thin = 5,
                refresh = 0,
                seed = 1,
                backend = "cmdstan")


## ----fitModel2i, results='markdown', eval=TRUE, mhidden=TRUE------------------
day_brm2 |> ggpredict() |> plot(show_data=TRUE)
day_brm2 |>
    ggpredict() |>
    plot(show_data=TRUE) |>
    wrap_plots() &
    scale_y_log10()


## ----fitModel2i2, results='markdown', eval=TRUE, mhidden=TRUE-----------------
day_brm2 |> ggemmeans(~TREAT) |> plot(show_data=TRUE)
day_brm2 |>
    ggemmeans(~TREAT) |>
    plot(show_data=TRUE) |>
    wrap_plots() &
    scale_y_log10()


## ----fitModel2i3, results='markdown', eval=TRUE, mhidden=TRUE-----------------
day_brm2 |> conditional_effects() |>  plot(points=TRUE)
day_brm2 |>
    conditional_effects() |>
    plot(points=TRUE) |>
    wrap_plots() &
    scale_y_log10()


## ----fitModel2j, results='markdown', eval=TRUE, mhidden=TRUE------------------
day_brm3 <- day_brm2 |> update(sample_prior = 'yes', refresh = 0)


## ----fitModel2k, results='markdown', eval=TRUE, mhidden=TRUE------------------
day_brm3 |> get_variables()
day_brm3 |> hypothesis('TREATALG2<0') |> plot()
day_brm3 |> hypothesis('TREATNB<0') |> plot()
day_brm3 |> hypothesis('TREATS<0') |> plot()
day_brm3 |> SUYR_prior_and_posterior()


## ----fitModel2l, results='markdown', eval=TRUE, mhidden=TRUE------------------
day_brm3 |> standata()
day_brm3 |> stancode()


## ----modelValidation1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_mcmc()


## ----modelValidation1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(day_rstanarm3, plotfun='mcmc_trace')


## ----modelValidation1c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(day_rstanarm3, 'acf_bar')


## ----modelValidation1d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(day_rstanarm3, 'rhat_hist')


## ----modelValidation1e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(day_rstanarm3, 'neff_hist')


## ----Validation1f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
plot(day_rstanarm3, 'combo')
plot(day_rstanarm3, 'violin')


## ----modelValidation1g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_trace(day_rstanarm3)


## ----modelValidation1h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_ac(day_rstanarm3)


## ----modelValidation1i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_rhat(day_rstanarm3)


## ----modelValidation1j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_ess(day_rstanarm3)


## ----modelValidation1k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
stan_dens(day_rstanarm3, separate_chains = TRUE)


## ----modelValidation1l, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=12----
day_ggs <- ggs(day_rstanarm3)
ggs_traceplot(day_ggs)


## ----modelValidation1m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=12----
ggs_autocorrelation(day_ggs)


## ----modelValidation1n, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_Rhat(day_ggs)


## ----modelValidation1o, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_effective(day_ggs)


## ----modelValidation1p, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_crosscorrelation(day_ggs)


## ----modelValidation1q, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_grb(day_ggs)


## ----modelValidation2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
day_brm3$fit |> stan_trace()
day_brm3$fit |> stan_trace(inc_warmup=TRUE)


## ----modelValidation2h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
day_brm3$fit |> stan_ac()


## ----modelValidation2i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
day_brm3$fit |> stan_rhat()


## ----modelValidation2j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
day_brm3$fit |> stan_ess()


## ----modelValidation2k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
day_brm3$fit |> stan_dens(separate_chains = TRUE)


## ----modelValidation3a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation3b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(day_rstanarm3,  plotfun='dens_overlay')


## ----modelValidation3c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#pp_check(day_rstanarm3, plotfun='error_scatter_avg')


## ----modelValidation3d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(day_rstanarm3, x=as.numeric(day$TREAT), plotfun='error_scatter_avg_vs_x')


## ----modelValidation3e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pp_check(day_rstanarm3, x=as.numeric(day$TREAT), plotfun='intervals')


## ----modelValidation3g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(day_rstanarm3)


## ----modelValidation4a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
preds <- posterior_predict(day_rstanarm3,  ndraws=250,  summary=FALSE)
day_resids <- createDHARMa(simulatedResponse = t(preds),
                            observedResponse = day$BARNACLE,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse = TRUE)
plot(day_resids)


## ----modelValidation5a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation5b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
day_brm3 |> pp_check(type = 'dens_overlay', ndraws=200)


## ----modelValidation5c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#day_brm3 |> pp_check(type = 'error_scatter_avg')


## ----modelValidation5e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
day_brm3 |> pp_check(group='TREAT', type='intervals')


## ----modelValidation5g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(day_brm3)


## ----modelValidation6a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
preds <- day_brm3 |> posterior_predict(ndraws = 250,  summary = FALSE)
day_resids <- createDHARMa(simulatedResponse = t(preds),
                            observedResponse = day$BARNACLE,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse = TRUE)
day_resids |> plot()

## ----modelValidation6aa, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=10----
day_resids <- make_brms_dharma_res(day_brm3, integerResponse = TRUE)
wrap_elements(~testUniformity(day_resids)) +
               wrap_elements(~plotResiduals(day_resids, form = factor(rep(1, nrow(day))))) +
               wrap_elements(~plotResiduals(day_resids, quantreg = TRUE)) +
               wrap_elements(~testDispersion(day_resids))



## ----partialPlot1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_rstanarm3 |> ggpredict() |> plot(show_data=TRUE)


## ----partialPlot1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_rstanarm3 |> ggemmeans(~TREAT,  type='fixed') |> plot(show_data=TRUE)


## ----partialPlot1c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_rstanarm3 |>
    epred_draws(newdata=day) |>
  median_hdci() |>
  ggplot(aes(x=TREAT, y=.epred)) +
  geom_pointrange(aes(ymin=.lower, ymax=.upper)) +
  geom_line() +
  geom_point(data=day,  aes(y=BARNACLE,  x=TREAT))


## ----partialPlot2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_brm3 |> conditional_effects() |> plot(points = TRUE)


## ----partialPlot2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_brm3 |> ggpredict() |> plot(show_data = TRUE)


## ----partialPlot2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_brm3 |> ggemmeans(~TREAT) |> plot(show_data = TRUE)


## ----partialPlot2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_brm3 |>
    epred_draws(newdata = day) |>
    median_hdci() |>
    ggplot(aes(x = TREAT, y = .epred)) +
    geom_pointrange(aes(ymin = .lower, ymax = .upper)) +
    geom_line() +
    geom_point(data = day,  aes(y = BARNACLE,  x = TREAT))


## ----summariseModel1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
summary(day_rstanarm3)


## ----summariseModel1a1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=FALSE----
day_sum <- summary(day_rstanarm3)


## ----summariseModel1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tidyMCMC(day_rstanarm3$stanfit, estimate.method='median',  conf.int=TRUE,  conf.method='HPDinterval',  rhat=TRUE, ess=TRUE)

## ----summariseModel1b1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
day_tidy <- tidyMCMC(day_rstanarm3$stanfit, estimate.method='median',  conf.int=TRUE,  conf.method='HPDinterval',  rhat=TRUE, ess=TRUE)


## ----summariseModel1dd, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_rstanarm3$stanfit |>
    summarise_draws(median,
                    HDInterval::hdi,
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel1d2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_rstanarm3$stanfit |>
    summarise_draws(median,
                    ~HDInterval::hdi(.x, credMass = 0.9),
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel1d3, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_rstanarm3$stanfit |>
    summarise_draws(
        ~ median(exp(.x)),
        ~HDInterval::hdi(exp(.x)),
        rhat, length, ess_bulk, ess_tail)


## ----summariseModel1m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_rstanarm3$stanfit |> as_draws_df()
day_rstanarm3$stanfit |>
  as_draws_df() |>
  summarise_draws(
    median,
    ~ HDInterval::hdi(.x),
    rhat,
    length,
    ess_bulk, ess_tail
  )

day_rstanarm3$stanfit |>
    as_draws_df() |>
    exp() |>
    summarise_draws(
        median,
        ~ HDInterval::hdi(.x),
        rhat,
        length,
        ess_bulk, ess_tail
    )


## ----summariseModel1c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_rstanarm3 |> get_variables()
day_draw <- day_rstanarm3 |> gather_draws(`.Intercept.*|.*TREAT.*`,  regex=TRUE)
day_draw

exceedP <- function(x, Val = 0) mean(x>Val)

day_rstanarm3 |>
    gather_draws(`.Intercept.*|.*TREAT.*`,  regex=TRUE) |>
    mutate(.value = exp(.value)) |>
    summarise_draws(median,
                    HDInterval::hdi,
                    rhat,
                    length,
                    ess_bulk,
                    ess_tail,
                    ~ exceedP(.x, 1))


## ----summariseModel1c1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_draw |> median_hdci(.value)


## ----summariseModel1c8, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_draw |> median_hdci(exp(.value))


## ----summariseModel1c3, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
day_gather <- day_rstanarm3 |> gather_draws(`.Intercept.*|.*TREAT.*`,  regex=TRUE) |>
  median_hdci()


## ----summariseModel1j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_rstanarm3 |> plot(plotfun='mcmc_intervals')


## ----summariseModel1c4, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
day_rstanarm3 |>
  gather_draws(`.Intercept.*|.*TREAT.*`, regex=TRUE) |>
  ggplot() +
  stat_halfeye(aes(x=.value,  y=.variable)) +
  facet_wrap(~.variable, scales='free')

day_rstanarm3 |>
  gather_draws(`.*TREAT.*`, regex=TRUE) |>
  ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    geom_vline(xintercept = 0, linetype = 'dashed')

day_rstanarm3 |>
  gather_draws(`.*TREAT.*`, regex=TRUE) |>
  ggplot() +
    stat_halfeye(aes(x=exp(.value),  y=.variable)) +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    scale_x_continuous(trans = scales::log2_trans())


## ----summariseModel1c7, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
day_rstanarm3 |>
  gather_draws(`.*TREAT.*`, regex=TRUE) |>
  ggplot() +
    geom_density_ridges(aes(x=.value, y = .variable), alpha=0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed')
##Or on a fractional scale
day_rstanarm3 |>
  gather_draws(`.*TREAT.*`, regex=TRUE) |>
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
day_rstanarm3 |> tidy_draws()


## ----summariseModel1e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_rstanarm3 |> spread_draws(`.Intercept.*|.*TREAT.*`,  regex=TRUE)


## ----summariseModel1f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_rstanarm3 |> posterior_samples() |> as_tibble()


## ----summariseModel1g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
#day_rstanarm3 |> bayes_R2() |> median_hdci


## ----summariseModel2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_brm3 |> summary()


## ----summariseModel2a1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=FALSE----
day_sum <- summary(day_brm3)


## ----summariseModel2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_brm3$fit |> tidyMCMC(estimate.method = 'median',
                          conf.int = TRUE,
                          conf.method = 'HPDinterval',
                          rhat = TRUE,
                          ess = TRUE)

## ----summariseModel2b1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
day_tidy <- tidyMCMC(day_brm3$fit, estimate.method='median',  conf.int=TRUE,  conf.method='HPDinterval',  rhat=TRUE, ess=TRUE)


## ----summariseModel2m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_brm3 |> as_draws_df()
day_brm3 |>
  as_draws_df() |>
  summarise_draws(
    median,
    HDInterval::hdi,
    rhat,
    length,
    ess_bulk, ess_tail
  )

day_brm3 |>
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
day_brm3 |> get_variables()
day_draw <- day_brm3 |>
    gather_draws(`.Intercept.*|.*TREAT.*`,  regex = TRUE)
day_draw

day_brm3 |>
    gather_draws(`.Intercept.*|.*TREAT.*`,  regex = TRUE) |>
    mutate(.value = exp(.value)) |>
    summarise_draws(median,
                    ~HDInterval::hdi(.x, credMass = 0.95),
                    rhat,
                    length,
                    ess_bulk, ess_tail)


exceedP <- function(x, Val = 0) mean(x>Val)
day_brm3 |>
    tidy_draws() |>
    exp() |>
    dplyr::select(starts_with("b_")) |>
    summarise_draws(median,
                    ~HDInterval::hdi(.x, credMass = 0.9),
                    rhat,
                    ess_bulk, ess_tail,
                    ~exceedP(.x, 1))


## ----summariseModel2c1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_draw |> median_hdci(.value)


## ----summariseModel2c5, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
day_draw |>
  median_hdci(exp(.value))


## ----summariseModel2c3, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
day_gather <- day_brm3 |> gather_draws(`.Intercept.*|.*TREAT.*`,  regex = TRUE) |>
  median_hdci()


## ----summariseModel2c4a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
day_brm3 |>
  gather_draws(`.Intercept.*|.*TREAT.*`, regex=TRUE) |>
  ggplot() +
  stat_halfeye(aes(x=.value,  y=.variable)) +
  facet_wrap(~.variable, scales='free')


## ----summariseModel2j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_brm3$fit |> plot(type='intervals')


## ----summariseModel2d5, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
## Link scale
day_brm3 |>
    gather_draws(`.Intercept.*|.*TREAT.*`, regex=TRUE) |>
    ggplot() +
    stat_slab(aes(x = .value, y = .variable,
                  fill = stat(ggdist::cut_cdf_qi(cdf,
                           .width = c(0.5, 0.8, 0.95),
                           labels = scales::percent_format())
                           )), color='black') +
    geom_vline(xintercept=0, linetype='dashed') +
    scale_fill_brewer('Interval', direction = -1, na.translate = FALSE)
## Fractional scale
day_brm3 |>
    gather_draws(`.Intercept.*|.*TREAT.*`, regex=TRUE) |>
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
day_brm3 |>
  gather_draws(`.Intercept.*|.*TREAT.*`, regex=TRUE) |>
  ggplot() +
  stat_halfeye(aes(x=.value,  y=.variable)) +
  facet_wrap(~.variable, scales='free')

day_brm3 |>
  gather_draws(`.*TREAT.*`, regex=TRUE) |>
  ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    geom_vline(xintercept = 0, linetype = 'dashed')

day_brm3 |>
  gather_draws(`.*TREAT.*`, regex=TRUE) |>
  ggplot() +
    stat_halfeye(aes(x=exp(.value),  y=.variable)) +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    scale_x_continuous(trans = scales::log2_trans())


## ----summariseModel2c7, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
day_brm3 |>
  gather_draws(`.*TREAT.*`, regex=TRUE) |>
  ggplot() +
    geom_density_ridges(aes(x=.value, y = .variable), alpha=0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed')
##Or on a fractional scale
day_brm3 |>
  gather_draws(`.*TREAT.*`, regex=TRUE) |>
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
day_brm3 |> tidy_draws()


## ----summariseModel2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_brm3 |> spread_draws(`.Intercept.*|.*TREAT.*`,  regex=TRUE)


## ----summariseModel2f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_brm3 |> posterior_samples() |> as_tibble()


## ----summariseModel2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_brm3 |> bayes_R2(summary=FALSE) |> median_hdci()


## ----summariseModel2k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
0.1 * sd(log(day$BARNACLE))
day_brm3 |> rope(range = c(-0.04, 0.04))
rope(day_brm3, range = c(-0.04, 0.04)) |> plot()

## Or based on fractional scale
day_brm3 |>
    as_draws_df('^b_TREAT.*', regex = TRUE) |>
    exp() |>
    ## equivalence_test(range = c(0.9, 1.1))
    rope(range = c(0.9, 1.1))
day_mcmc <-
    day_brm3 |> as_draws_df('^b_TREAT.*', regex = TRUE) |>
    exp()
day_mcmc |>
    rope(range = c(0.9, 1.1))
day_mcmc |>
    rope(range = c(0.9, 1.1)) |>
    plot(day_mcmc)

day_mcmc |>
    equivalence_test(range = c(0.9, 1.1))


## ----Probability1a, results='markdown', echo=1,eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_rstanarm3 |> emmeans(pairwise ~TREAT, type='response')
day_pairwise <- (day_rstanarm3 |>
  emmeans(pairwise ~TREAT, type='response') |>
  confint())$contrasts |>
  as.data.frame()


## ----Probability1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_em = emmeans(day_rstanarm3, pairwise~TREAT, type='link')$contrasts |>
    gather_emmeans_draws() |>
    mutate(Fit=exp(.value))
day_em |> head()
day_em |> group_by(contrast) |>
    ggplot(aes(x=Fit)) +
    geom_histogram() +
    geom_vline(xintercept=1, color='red') +
    facet_wrap(~contrast, scales='free')
day_em |> group_by(contrast) |>
    ggplot(aes(x=Fit)) +
    geom_density_ridges_gradient(aes(y = contrast, fill = stat(x)),
                                 alpha = 0.4, color = 'white',
                                 quantile_lines = TRUE,
                                 quantiles = c(0.025, 0.975)) +
    geom_vline(xintercept=1, linetype = 'dashed') +
    scale_fill_viridis_c(option = "C") +
    scale_x_continuous(trans = scales::log2_trans())

day_em |>
  ggplot(aes(x = Fit)) +
    geom_density_ridges_gradient(aes(y = contrast,
                                     fill = factor(stat(x>0))),
                                 alpha=0.4, colour = 'white',
                                 quantile_lines = TRUE,
                                 quantiles = c(0.025, 0.975)) +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    scale_x_continuous(trans = scales::log2_trans()) +
    scale_fill_viridis_d()

day_em |> group_by(contrast) |> median_hdi()
# Probability of effect
day_em |> group_by(contrast) |> summarize(P=mean(Fit>1))
day_em |> group_by(contrast) |> summarize(P=mean(Fit<1))
##Probability of effect greater than 10%
day_em |> group_by(contrast) |> summarize(P=mean(Fit>1.1))


## ----Probability2a, results='markdown', echo=TRUE,eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_brm3 |>
    emmeans(~TREAT, type = 'response') |>
    pairs()
#OR
day_pairwise <- day_brm3 |>
    emmeans(~TREAT, type = 'response') |>
    pairs() |>
    as.data.frame()
##OR
day_brm3 |>
    emmeans(~TREAT) |>
    pairs() |>
    tidy_draws() |>
    exp() |>
    summarise_draws(median)
##OR
day_brm3 |>
    emmeans(~TREAT) |>
    pairs() |>
    gather_emmeans_draws() |>
    mutate(.ratio = exp(.value)) |>
    ## median_hdci(.ratio)
    summarise(
        median_hdci(.ratio),
        Pl = mean(.ratio < 1),
        Pg = mean(.ratio > 1),
        ROPE = rope(.ratio, range = c(0.9, 1.1))$ROPE_Percentage)
    ## summarise(across(c(.value, .ratio), c(median, HDInterval::hdi)))
    ## summarise(across(c(.value, .ratio), c(median_hdci)))

day_mcmc <-
    day_brm3 |>
    emmeans(~TREAT) |>
    pairs() |>
    tidy_draws() |>
    dplyr::select(-.chain, -.iteration, -.draw) |>
    exp()
rope(day_mcmc, range = c(0.9, 1.1)) |> plot()


## ----Probability2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
day_brm3 |> emmeans(~TREAT, type = 'response') |>
    pairs()
day_em <- day_brm3 |>
    emmeans(~TREAT, type = 'link') |>
    pairs() |>
    gather_emmeans_draws() |>
    mutate(Fit = exp(.value))

day_em |>
  ggplot(aes(x = Fit)) +
  geom_density_ridges_gradient(aes(y = contrast,
                                   fill = factor(stat(x < 1.1 & x > 0.9))),
                               alpha=0.4, colour = 'white',
                               quantile_lines = TRUE,
                               quantiles = c(0.025, 0.975),
                               show.legend = FALSE) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  scale_x_continuous("Effect size (fold scale)", trans = scales::log2_trans()) +
  scale_y_discrete("") +
  scale_fill_viridis_d() +
  theme_classic()

day_em |> group_by(contrast) |>
    ggplot(aes(x=Fit)) +
    geom_density_ridges_gradient(aes(y = contrast, fill = stat(x)),
                                 alpha = 0.4, color = 'white',
                                 quantile_lines = TRUE,
                                 quantiles = c(0.025, 0.975)) +
    geom_vline(xintercept=1, linetype = 'dashed') +
    scale_fill_viridis_c(option = "C") +
    scale_x_continuous(trans = scales::log2_trans())

day_em |>
  ggplot(aes(x = Fit)) +
    geom_density_ridges_gradient(aes(y = contrast,
                                     fill = factor(stat(x>0))),
                                 alpha=0.4, colour = 'white',
                                 quantile_lines = TRUE,
                                 quantiles = c(0.025, 0.975)) +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    scale_x_continuous(trans = scales::log2_trans()) +
    scale_fill_viridis_d() +
    scale_x_continuous(trans = scales::log2_trans())

day_em |> head()
day_em |>
    group_by(contrast) |>
    ggplot(aes(x = Fit)) +
    ## geom_histogram() +
    geom_halfeyeh() +
    geom_vline(xintercept = 1, color = 'red') +
    facet_wrap(~contrast, scales = 'free')
day_em |>
    group_by(contrast) |>
    median_hdi()
# Probability of effect
day_em |> group_by(contrast) |> summarize(P = mean(Fit>1))
##Probability of effect greater than 10%
day_em |> group_by(contrast) |> summarize(P = mean(Fit>1.1))

## Effect size on absolute scale
day_em <- day_brm3 |>
    emmeans(~TREAT, type = 'link') |>
    regrid() |>
    pairs() |>
    gather_emmeans_draws() |>
    median_hdci(.value)

day_brm3 |>
    emmeans(~TREAT, type = 'link') |>
    regrid() |>
    pairs() |>
    gather_emmeans_draws() |>
    group_by(contrast) |>
    ggplot(aes(x = .value)) +
    geom_halfeyeh() +
    geom_vline(xintercept = 0, color = 'red') +
    facet_wrap(~contrast, scales = 'free')


day_brm3 |>
    emmeans(~TREAT, type = 'response') |>
    pairs() |>
    gather_emmeans_draws() |>
    mutate(.value = exp(.value)) |>
    ggplot(aes(x = .value)) +
    geom_density_ridges_gradient(aes(y = contrast,
                                     fill = factor(stat(x>0))),
                                 alpha=0.4, colour = 'white',
                                 quantile_lines = TRUE,
                                 quantiles = c(0.025, 0.975)) +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    scale_x_continuous(trans = scales::log2_trans()) +
    scale_fill_viridis_d()



## ----Probability1c, results='markdown', echo=TRUE, eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
##Planned contrasts
cmat<-cbind('Alg2_Alg1'=c(-1,1,0,0),
              'NB_S'=c(0,0,1,-1),
             'Alg_Bare'=c(0.5,0.5,-0.5,-0.5),
             'Alg_NB'=c(0.5,0.5,-1,0))
# On the link scale
emmeans(day_rstanarm3, ~TREAT, contr=list(TREAT=cmat), type='link')
# On the response scale
emmeans(day_rstanarm3, ~TREAT, contr=list(TREAT=cmat), type='response')

day_em = emmeans(day_rstanarm3, ~TREAT, contr=list(TREAT=cmat), type='link')$contrasts |>
      gather_emmeans_draws() |> mutate(Fit=exp(.value))
day_em |> group_by(contrast) |> mean_hdi()
# Probability of effect
day_em |> group_by(contrast) |> summarize(P=mean(Fit>1))
##Probability of effect greater than 10%
day_em |> group_by(contrast) |> summarize(P=mean(Fit>1.5))


day_sum <- day_em |>
  group_by(contrast) |>
  median_hdci(.width=c(0.8, 0.95))
day_sum
ggplot(day_sum) +
  geom_hline(yintercept=1, linetype='dashed') +
  geom_pointrange(aes(x=contrast, y=Fit, ymin=Fit.lower, ymax=Fit.upper, size=factor(.width)),
                  show.legend = FALSE) +
  scale_size_manual(values=c(1, 0.5)) +
  coord_flip()

g1 <- ggplot(day_sum) +
  geom_hline(yintercept=1) +
  geom_pointrange(aes(x=contrast, y=Fit, ymin=Fit.lower, ymax=Fit.upper, size=factor(.width)), show.legend = FALSE) +
  scale_size_manual(values=c(1, 0.5)) +
  scale_y_continuous(trans=scales::log2_trans(),  breaks=c(0.5, 1, 2, 4)) +
  coord_flip() +
  theme_classic()
g1


## ----Probability2c, results='markdown', echo=TRUE, eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
##Planned contrasts
cmat<-cbind('Alg2_Alg1' = c(-1,1,0,0),
            'NB_S' = c(0,0,1,-1),
            'Alg_Bare' = c(0.5,0.5,-0.5,-0.5),
            'Alg_NB' = c(0.5,0.5,-1,0))
# On the link scale
day_brm3 |> emmeans(~TREAT, type = 'link') |>
    contrast(method = list(TREAT = cmat))
# On the response scale
day_brm3 |>
    emmeans(~TREAT, type = 'response') |>
    contrast(method = list(TREAT = cmat))

day_em <- day_brm3 |>
    emmeans(~TREAT, type = 'link') |>
    contrast(method = list(TREAT = cmat)) |>
    gather_emmeans_draws() |>
    mutate(Fit = exp(.value))
day_em |> median_hdi(Fit)
# Probability of effect
day_em |> summarize(P = mean(Fit>1))
##Probability of effect greater than 50%
day_em |> summarize(P = mean(Fit>1.5))

day_sum <- day_em |>
  group_by(contrast) |>
  median_hdci(.value, .width = c(0.8, 0.95))
day_sum
g1 <- ggplot(day_sum) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_pointrange(aes(y = contrast, x = .value, xmin = .lower, xmax = .upper,
                      size = factor(.width)),
                  show.legend  =  FALSE) +
  scale_size_manual(values = c(1, 0.5))

day_sum <- day_em |>
  group_by(contrast) |>
  median_hdci(Fit, .width = c(0.8, 0.95))
day_sum
g1 <- ggplot(day_sum) +
  geom_vline(xintercept = 1, linetype='dashed') +
  geom_pointrange(aes(y = contrast, x = Fit, xmin = .lower, xmax = .upper, size = factor(.width)), show.legend = FALSE) +
  scale_size_manual(values = c(1, 0.5)) +
  scale_x_continuous(trans = scales::log2_trans(),  breaks = c(0.5, 0.8,1,1.2, 1.5, 2, 4)) +
  theme_classic()
g1

g1a <-
    day_em |>
    ggplot() +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    #geom_vline(xintercept = 1.5, alpha=0.3, linetype = 'dashed') +
    stat_slab(aes(x = Fit, y = contrast,
                  fill = stat(ggdist::cut_cdf_qi(cdf,
                            .width = c(0.5, 0.8, 0.95),
                            labels = scales::percent_format())
                            )), color = 'black') +
    scale_fill_brewer('Interval', direction  =  -1, na.translate = FALSE) +
    scale_x_continuous('Effect', trans = scales::log2_trans(),
                       breaks = c(0.5, 0.8,1,1.2, 1.5, 2, 4)) +
    scale_y_discrete('', breaks=c('TREAT.NB_S', 'TREAT.Alg2_Alg1', 'TREAT.Alg_NB', 'TREAT.Alg_Bare'),
                     labels = c('Nat. Bare vs Scapped', 'Algae 1 vs 2', 'Algae vs Nat. Bare', 'Algae vs Bare')) +
    theme_classic()
g1 + g1a

library(ggridges)
day_em |>
  ggplot() +
    geom_density_ridges(aes(x=Fit, y = contrast), alpha=0.4) +
    geom_vline(xintercept = 1, linetype = 'dashed')
day_em |>
  ggplot() +
    geom_density_ridges_gradient(aes(x=Fit,
                                     y = contrast,
                                     fill = stat(x)),
                                 alpha=0.4, colour = 'white',
                                 quantile_lines = TRUE,
                                 quantiles = c(0.025, 0.975)) +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    scale_x_continuous(trans = scales::log2_trans()) +
    scale_fill_viridis_c(option = "C")
day_em |>
  ggplot() +
    geom_density_ridges_gradient(aes(x=100*(Fit-1),
                                     y = contrast,
                                     fill = stat(x)),
                                 alpha=0.4, colour = 'white',
                                 quantile_lines = TRUE,
                                 quantiles = c(0.025, 0.975)) +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    scale_x_continuous('Percentage change') +
    scale_fill_viridis_c(option = "C")
##Or on a fractional scale
day_brm3 |>
  gather_draws(`.*TREAT.*`, regex=TRUE) |>
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

## Perhaps nicer
g1b <-
  day_brm3 |>
    emmeans(~TREAT, type = 'link') |>
    contrast(method = list(TREAT = cmat)) |>
    gather_emmeans_draws() |>
    mutate(Fit = exp(.value)) |>
  ggplot() +
  geom_density_ridges_gradient(aes(x=Fit,
    y = contrast,
    fill = factor(stat(x < 1.1 & x > 0.9))),
    alpha=0.4, colour = 'white',
    quantile_lines = TRUE,
    quantiles = c(0.025, 0.975)) +
  geom_vline(xintercept = 1, linetype = 'dashed') +
  scale_x_continuous("Effect size (fold scale)", trans = scales::log2_trans()) +
  scale_fill_viridis_d("Inside ROPE") +
  theme_classic()

g1 + g1b


## ----summaryFig1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=3----
newdata = emmeans(day_rstanarm3, ~TREAT, type='response') |>
    as.data.frame()
newdata
## A quick version
g2 <- ggplot(newdata, aes(y=rate, x=TREAT)) +
    geom_pointrange(aes(ymin=lower.HPD, ymax=upper.HPD)) +
    theme_classic()

g2 + g1


## ----summaryFig2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=3----
newdata <- day_brm3 |>
    emmeans(~TREAT, type='response') |>
    as.data.frame()
newdata
## A quick version
g2 <- ggplot(newdata, aes(y=rate, x=TREAT)) +
    geom_pointrange(aes(ymin=lower.HPD, ymax=upper.HPD)) +
    scale_y_continuous('Number of newly recruited barnacles') +
    scale_x_discrete('', breaks=c('ALG1', 'ALG2', 'NB', 'S'),
                     labels = c('Algae 1', 'Algae 2', 'Nat. Bare', 'Scraped')) +
    theme_classic()

g2 + g1


## ----summaryFig2a2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=3----
(g2 + ggtitle('a)')) + (g1a + ggtitle('b)'))
(g2 + ggtitle('a)')) + (g1b + ggtitle('b)'))

