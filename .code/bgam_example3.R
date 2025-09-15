## ----setup, include=FALSE, warnings=FALSE, message=FALSE----------------------
knitr::opts_chunk$set(cache.lazy = FALSE, tidy='styler')


## -----------------------------------------------------------------------------
#| label: libraries
#| output: false
#| eval: true
#| warning: false
#| message: false
#| cache: false

library(tidyverse)  #for data wrangling etc
library(cmdstanr)   #for cmdstan
library(brms)       #for fitting models in STAN
library(standist)   #for exploring distributions
library(coda)       #for diagnostics
library(bayesplot)  #for diagnostics
library(ggmcmc)     #for MCMC diagnostics
library(DHARMa)     #for residual diagnostics
library(rstan)      #for interfacing with STAN
library(emmeans)    #for marginal means etc
library(broom)      #for tidying outputs
library(tidybayes)  #for more tidying outputs
library(HDInterval) #for HPD intervals
library(ggeffects)  #for partial plots
library(broom.mixed)#for summarising models
library(posterior)  #for posterior draws
library(ggeffects)  #for partial effects plots
library(patchwork)  #for multi-panel figures
library(bayestestR) #for ROPE
library(see)        #for some plots
library(easystats)     #framework for stats, modelling and visualisation
library(mgcv)
library(gratia)
theme_set(theme_grey()) #put the default ggplot theme back
source('helperFunctions.R')


## -----------------------------------------------------------------------------
#| label: readData
#| cache: false
wq <- read_csv("../data/aims_wq.csv", trim_ws = TRUE)


## -----------------------------------------------------------------------------
#| label: prepareData
#| cache: false
wq <- wq |>
    mutate(
        Site = factor(Site),
        Region = factor(Region),
        Subregion = factor(Subregion),
        Season = factor(Season)
    )


## ----EDA1a, results='markdown', eval=TRUE, mhidden=TRUE-----------------------
wq |> ggplot(aes(y=NOx, x=Date)) + geom_point()


## -----------------------------------------------------------------------------
#| label: EDA1b
#| eval: true
#| mhidden: true
#| fig-width: 15
#| fig-height: 12
#| out-width: 800px
#| warning: false
#| message: false
ggplot(wq, aes(y = NOx, x = Date)) +
    geom_point(aes(colour = Season)) +
    geom_smooth() +
    facet_wrap(~Region + Site, scales = 'free_y')


## -----------------------------------------------------------------------------
#| label: EDA1c
#| eval: true
#| mhidden: true
#| fig-width: 15
#| fig-height: 12
#| out-width: 800px
#| warning: false
#| message: false
ggplot(wq, aes(y=NOx, x=Date)) +
  geom_point() +
  geom_smooth() +
    scale_y_log10() +
  scale_y_continuous(trans=scales::pseudo_log_trans()) +
  facet_wrap(~Site,  scales='free_y')



## -----------------------------------------------------------------------------
#| label: prepareData2
#| eval: true
#| mhidden: true
#| cache: false
wq <- wq |>
    mutate(dt_num = decimal_date(Date)) |>
    mutate(Mnth = month(Date)) |>
    mutate(NOx_flag = ifelse(NOx == 0.01, "left", "none"))


## -----------------------------------------------------------------------------
#| label: model1a
#| eval: true
#| cache: false
#| mhidden: true
wq_sub <- wq |>
  filter(Site == "Double Island") |>
  droplevels()

wq_sub |>
  ggplot(aes(y = NOx, x = Date)) +
  geom_point(aes(colour = Season)) +
  geom_smooth()
  scale_y_log10()

wq_sub |>
  filter(!is.na(NOx)) |>
  ggplot(aes(y=NOx, x=Date)) +
  geom_point() +
  geom_smooth(method = 'gam', formula = y ~ s(x),
    method.args = list(family = Gamma(link = "log")))

wq_sub |>
    pull(NOx) |>
    summary()


## ----fitModel1a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
wq.form <- bf(NOx | cens(NOx_flag) ~ s(dt_num), family = Gamma(link = "log"))
get_prior(wq.form, data =  wq_sub)

wq.brm <- brm(wq.form,
                 data = wq_sub,
                 prior = prior(normal(0, 2.5), class = 'b'),
                 sample_prior = 'only',
                 iter = 5000,
                 warmup = 1000,
                 chains = 3,
                 thin = 5,
                 backend = 'cmdstan',
                 refresh = 0)


## ----fitModel1d2, results='markdown', eval=TRUE, mhidden=TRUE-----------------
wq.brm |> conditional_effects() |>  plot(points=TRUE) |>
  scale_y_log10()


## ----fitModel2a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
wq.form <- bf(NOx | cens(NOx_flag) ~ s(scale(dt_num)), family = Gamma(link = "log"))
wq_sub |> summarise(
    median(log(NOx)), mad(log(NOx)),
    mad(log(NOx)) / mad(scale(dt_num))
)
get_prior(wq.form, data = wq_sub)
priors <- prior(normal(-3.1, 1),  class='Intercept') +
  prior(normal(0, 10), class='b') +
  prior(gamma(0.01, 0.01),  class='shape') +
  prior(student_t(3, 0, 10), class =  "sds")

wq.brm2 <- brm(wq.form,
                 data = wq_sub,
                 prior = priors,
                 sample_prior = 'only',
                 iter = 5000,
                 warmup = 1000,
                 thin = 5,
                 chains = 3, cores =  3,
                 backend = 'cmdstan',
                 seed =  123,
                 control =  list(adapt_delta = 0.99, max_treedepth = 20),
                 refresh = 0)


## ----fitModel2b, results='markdown', eval=TRUE, mhidden=TRUE------------------
wq.brm2 |> conditional_effects() |>  plot(points=TRUE)
wq.brm2 |> conditional_effects() |>  plot(points=TRUE) |>
  _[[1]] +
  scale_y_log10()

wq.brm2 |> conditional_effects() |>
  plot(points=TRUE, point_args = list(colour = "red", size = 4)) |>
  _[[1]] +
  scale_y_log10() +
  geom_point(data = wq_sub, inherit.aes = FALSE,
             aes(y = NOx, x = dt_num))


## ----fitModel2c, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
wq.brm3 <- update(wq.brm2, sample_prior = "yes", cores =  3, refresh =  0)


## ----fitModel2d, results='markdown', eval=TRUE, mhidden=TRUE------------------
wq.brm3 |> conditional_effects() |>  plot(points=TRUE) |>
  _[[1]] +
  scale_y_log10()

wq.brm3 |> conditional_effects(spaghetti = TRUE, ndraws =  250) |>  plot(points=TRUE) |>
  _[[1]] +
  scale_y_log10() +
  geom_point(data = wq_sub, inherit.aes = FALSE,
             aes(y = NOx, x = dt_num))

wq.brm3 |> conditional_effects() |>
  plot(points=TRUE, point_args = list(colour = "red", size = 4)) |>
  _[[1]] +
  scale_y_log10() +
  geom_point(data = wq_sub, inherit.aes = FALSE,
             aes(y = NOx, x = dt_num))

wq.brm3 |> conditional_effects() |>
  plot(points=TRUE, point_args = list(colour = "red", size = 4)) |>
  _[[1]] +
  geom_point(data = wq_sub, inherit.aes = FALSE,
             aes(y = NOx, x = dt_num))


## semi-manual plot
wq.brm3 |>
    conditional_effects(spaghetti = TRUE, ndraws = 250) |>
  as_draws_df() |>
  mutate(across(c(dt_num, estimate__, lower__, upper__), as.numeric)) |>
  mutate(Date =  lubridate::date_decimal(dt_num)) |>
  ggplot(aes(x = as.Date(Date), y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "orange", alpha = 0.3) +
  geom_line() +
  scale_y_continuous("Estimated NOx (µM)") +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  theme_bw()


# semi-manual spaghetti plot
newdata <- with(wq_sub, tibble(dt_num = seq(min(dt_num), max(dt_num), length = 100))) |>
  mutate(Date = as.Date(lubridate::date_decimal(dt_num)))
wq.brm3 |>
  add_epred_draws(newdata = newdata) |>
  filter(.draw < 250) |>
  ggplot(aes(x = Date, y = .epred, group = .draw)) +
  geom_line(alpha = 0.1) +
  stat_summary(
    aes(group = 1),
    fun = median,
    geom = "line",
    color = "orange",
    linewidth = 1
  ) +
  geom_point(data = wq_sub, inherit.aes = FALSE,
             aes(x = Date, y = NOx), alpha = 0.5) +
  scale_x_date("", date_breaks = "2 years", date_labels = "%Y") +
  ## scale_y_log10() +
  scale_y_continuous("Estimated NOx (µM)") +
  theme_bw()



## ----fitModel2ab, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE-----
wq.form <- bf(NOx | cens(NOx_flag) ~ gp(dt_num), family = Gamma(link = "log"))
wq_sub |> summarise(
    median(log(NOx)), mad(log(NOx)),
    mad(log(NOx)) / mad(scale(dt_num))
)
get_prior(wq.form, data = wq_sub)
priors <- prior(normal(-3.1, 1), class = "Intercept") +
  prior(gamma(0.01, 0.01), class = "shape") +
  prior(inv_gamma(1.5, 0.05), class = "lscale", coef = "gpdt_num") +
  prior(student_t(3, 0, 1), class =  "sdgp")

wq.brm2b <- brm(wq.form,
                 data = wq_sub,
                 prior = priors,
                 sample_prior = 'only',
                 iter = 5000,
                 warmup = 1000,
                 thin = 5,
                 chains = 3, cores =  3,
                 backend = 'cmdstan',
                 seed =  123,
                 control =  list(adapt_delta = 0.99, max_treedepth = 20),
                 refresh = 0)


## ----fitModel2bb, results='markdown', eval=FALSE, mhidden=TRUE----------------
# wq.brm2b |> conditional_effects() |>  plot(points=TRUE)
# wq.brm2b |> conditional_effects() |>  plot(points=TRUE) |>
#   _[[1]] +
#   scale_y_log10()


## ----fitModel2cb, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE-----
wq.brm3b <- update(wq.brm2b, sample_prior = "yes", cores =  3, refresh =  0)


## ----fitModel2db, results='markdown', eval=TRUE, mhidden=TRUE-----------------
wq.brm3b |> conditional_effects() |>  plot(points=TRUE) |>
  _[[1]] +
  scale_y_log10()
wq.brm3b |> conditional_effects(spaghetti = TRUE, ndraws =  250) |>  plot(points=TRUE) |>
  _[[1]] +
  scale_y_log10()


## ----fitModel2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width = 8, fig.height = 4, error = TRUE----
try({
wq.brm3 |> get_variables()
wq.brm3b |> hypothesis('sdgp_gpdt_num = 0', class = '') |> plot()
wq.brm3 |> hypothesis('bs_sscaledt_num_1 = 0', class = '') |> plot()
#wq.brm3b |> SUYR_prior_and_posterior()
})


## ----modelValidation2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
wq.brm3$fit |> stan_trace()


## ----modelValidation2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
wq.brm3$fit |> stan_ac()


## ----modelValidation2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
wq.brm3$fit |> stan_rhat()

## ----modelValidation2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
wq.brm3$fit |> stan_ess()


## ----modelValidation3a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
wq.brm3 |> pp_check( type='dens_overlay', ndraws=100)


## ----modelValidation3b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=10----

wq.resids <- make_brms_dharma_res(wq.brm3, integerResponse = FALSE)
wrap_elements(~testUniformity(wq.resids)) +
  wrap_elements(~plotResiduals(wq.resids, form = factor(rep(1, nrow(wq_sub))))) +
  wrap_elements(~plotResiduals(wq.resids, quantreg = TRUE)) +
  wrap_elements(~testDispersion(wq.resids))


## ----modelValidation3ab, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
wq.brm3b |> pp_check( type='dens_overlay', ndraws=100)


## ----modelValidation3bb, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=10----

wq.resids <- make_brms_dharma_res(wq.brm3b, integerResponse = FALSE)
wrap_elements(~testUniformity(wq.resids)) +
  wrap_elements(~plotResiduals(wq.resids, form = factor(rep(1, nrow(wq_sub))))) +
  wrap_elements(~plotResiduals(wq.resids, quantreg = TRUE)) +
  wrap_elements(~testDispersion(wq.resids))


## ----partialPlot1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
wq.brm3 |> conditional_effects() |> plot(points = TRUE)
wq.brm3 |>
    conditional_effects(spaghetti = TRUE, ndraws = 250) |>
    plot(points = TRUE)


newdata <- with(wq_sub, tibble(dt_num = seq(min(dt_num), max(dt_num), length = 100))) |>
  mutate(Date = as.Date(lubridate::date_decimal(dt_num)))
# semi-manual spaghetti plot
wq.brm3 |>
  add_epred_draws(newdata = newdata) |>
  filter(.draw < 250) |>
  ggplot(aes(x = Date, y = .epred, group = .draw)) +
  geom_line(alpha = 0.1) +
  stat_summary(
    aes(group = 1),
    fun = median,
    geom = "line",
    color = "orange",
    linewidth = 1
  ) +
  geom_point(data = wq_sub, inherit.aes = FALSE,
             aes(x = Date, y = NOx), alpha = 0.5) +
  scale_x_date("", date_breaks = "2 years", date_labels = "%Y") +
  ## scale_y_log10() +
  scale_y_log10("Estimated NOx (µM)") +
  theme_bw()


## ----summariseModel1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
wq.brm3 |> summary()


## ----summariseModel1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
wq.brm3 |> as.data.frame() |>
  dplyr::select(matches("^b_.*|^bs.*|^sds.*|^shape$|^s_s.*")) |>
  summarise_draws(median,
    HDInterval::hdi)


## ----additionalfitModel1a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
wq.form <- bf(NOx | cens(NOx_flag) ~ s(scale(dt_num), by = Season), family = Gamma(link = "log"))
get_prior(wq.form, data = wq_sub)
priors <- prior(normal(-3.1, 1),  class='Intercept') +
  prior(normal(0, 10), class='b') +
  prior(gamma(0.01, 0.01),  class='shape') +
  prior(student_t(3, 0, 10), class =  "sds")

wq.brm4 <- brm(wq.form,
                 data = wq_sub,
                 prior = priors,
                 sample_prior = 'only',
                 iter = 5000,
                 warmup = 1000,
                 chains = 3, cores =  3,
                 thin = 5,
                 backend = 'cmdstan',
                 seed =  123,
                 control =  list(adapt_delta =  0.99, max_treedepth = 20),
                 refresh = 0)


## ----additionalfitModel1b, results='markdown', eval=TRUE, mhidden=TRUE--------
wq.brm4 |> conditional_effects(effects = "dt_num:Season") |>  plot(points=TRUE)
wq.brm4 |> conditional_effects(effects = "dt_num:Season") |>  plot(points=TRUE) |>
  _[[1]] +
  scale_y_log10()


## ----additionalfitModel1c, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
wq.brm5 <- update(wq.brm4, sample_prior = "yes", cores =  3, refresh =  0)


## ----additionalfitModel1d, results='markdown', eval=TRUE, mhidden=TRUE--------
wq.brm5 |> conditional_effects(effects = "dt_num:Season") |>  plot(points=TRUE)
wq.brm5 |> conditional_effects(effects = "dt_num:Season") |>  plot(points=TRUE) |>
  _[[1]] +
  scale_y_log10()
wq.brm5 |> conditional_effects(effects = "dt_num:Season", spaghetti = TRUE, ndraws =  300) |>  plot(points=TRUE) |>
  _[[1]] +
  scale_y_log10()


## ----additionalfitModel1e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width = 8, fig.height = 4, error = TRUE----
try({
wq.brm5 |> get_variables()
wq.brm5 |> hypothesis('bs_sdt_num:SeasonDry_1 = 0', class = '') |> plot()
wq.brm5 |> hypothesis('sds_sdt_numSeasonDry_1 = 0', class = '') |> plot()
wq.brm5 |> SUYR_prior_and_posterior()
})


## ----additionalmodelValidation1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=3----
wq.resids <- make_brms_dharma_res(wq.brm5, integerResponse = FALSE)
wrap_elements(~testUniformity(wq.resids)) +
  ## wrap_elements(~plotResiduals(wq.resids, form = factor(rep(1, nrow(wq_sub))))) +
  wrap_elements(~plotResiduals(wq.resids, quantreg = TRUE)) +
  wrap_elements(~testDispersion(wq.resids))

## testTemporalAutocorrelation(wq.resids, time =  wq_sub$dt_num)
## resids1 <- recalculateResiduals(wq.resids, group =  wq_sub$dt_num, aggregateBy =  mean)
## testTemporalAutocorrelation(resids1, time =  wq_sub$dt_num)


## ----additionalfitModel3a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
wq.form <- bf(NOx | cens(NOx_flag) ~ s(scale(dt_num)) + s(scale(Mnth), bs = 'cc', k = 6),
  family = Gamma(link = "log"))
priors <- prior(normal(-3.1, 1),  class='Intercept') +
  prior(normal(0, 10), class='b') +
  prior(gamma(0.01, 0.01),  class='shape') +
  prior(student_t(3, 0, 10), class =  "sds")

wq.brm6 <- brm(wq.form,
                 data = wq_sub,
                 prior = priors,
                 knots = list(Mnth = seq(1, 12, len = 6)),
                 sample_prior = 'only',
                 iter = 5000,
                 warmup = 1000,
                 chains = 3, cores =  3,
                 thin = 5,
                 backend = 'cmdstan',
                 seed =  123,
                 control =  list(adapt_delta =  0.99, max_treedepth = 20),
                 refresh = 0)


## ----additionalfitModel3b, results='markdown', eval=TRUE, mhidden=TRUE--------
wq.brm6 |> conditional_effects(effects = "dt_num:Mnth") |>  plot(points=TRUE)
wq.brm6 |> conditional_effects(effects = "dt_num:Mnth") |>  plot(points=TRUE) |>
  _[[1]] +
  scale_y_log10()
wq.brm6 |> conditional_effects(effects = "Mnth") |>  plot(points=TRUE)


## ----additionalfitModel2b, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
wq.brm7 <- update(wq.brm6, sample_prior = "yes", cores =  3, refresh =  0)


## ----additionalfitModel2c, results='markdown', eval=TRUE, mhidden=TRUE--------
wq.brm7 |> conditional_effects(effects = "dt_num") |>  plot(points=TRUE)
wq.brm7 |> conditional_effects(effects = "Mnth") |>  plot(points=TRUE)


## ----additionalfitModel2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width = 8, fig.height = 4, error = TRUE----
try({
wq.brm7 |> get_variables()
wq.brm7 |> hypothesis('bs_sdt_num_1 = 0', class = '') |> plot()
wq.brm7 |> hypothesis('sds_sdt_num_1 = 0', class = '') |> plot()
wq.brm7 |> hypothesis('sds_sMnth_1 = 0', class = '') |> plot()
wq.brm7 |> SUYR_prior_and_posterior()
})


## ----additionalmodelValidation2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=10----
wq.resids <- make_brms_dharma_res(wq.brm7, integerResponse = FALSE)
wrap_elements(~testUniformity(wq.resids)) +
  ## wrap_elements(~plotResiduals(wq.resids, form = factor(rep(1, nrow(wq_sub))))) +
  wrap_elements(~plotResiduals(wq.resids, quantreg = TRUE)) +
  wrap_elements(~testDispersion(wq.resids))


## ----process3, results='markdown', eval=TRUE, mhidden=TRUE--------------------
wq_sub <- wq |>
  mutate(dt_num =  lubridate::decimal_date(Date)) |>
  group_by(reef.alias) |>
  mutate(Min=min(dt_num)) |>
  ungroup() |>
  filter(Min<2012, Region != 'Fitzroy', reef.alias != 'Daydream') |>
  droplevels()

## reef=wq |>
##   group_by(reef.alias) |>
##   dplyr:::summarise(Min=min(dt_num)) |>
##   filter(Min<2012) |>
##   pull(reef.alias)
## reef
## wq2=wq |> filter(reef.alias %in% reef) |> droplevels
ggplot(wq_sub, aes(y=NOx,x=dt_num)) +
    geom_point() +
    facet_wrap(~reef.alias, scales = 'free_y') +
    geom_smooth() +
    scale_y_log10()




## ----EDA2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=15, fig.height=12----
ggplot(wq_sub, aes(y=NOx,x=dt_num)) +
    geom_point() +
    facet_wrap(~reef.alias, scales='free_y')
## Some reefs dont have the full time series


## ----additionalfitModel8a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
wq.form <- bf(NOx | cens(NOx_flag) ~ s(scale(dt_num)) + (1|reef.alias),
  family = Gamma(link = "log"))
get_prior(wq.form, data = wq_sub)
priors <- prior(normal(-3.1, 1),  class='Intercept') +
  prior(normal(0, 10), class='b') +
  prior(gamma(0.01, 0.01),  class='shape') +
  prior(student_t(3, 0, 10), class =  "sds") +
  prior(student_t(3, 0, 1), class =  "sd")

wq.brm8 <- brm(wq.form,
                 data = wq_sub,
                 prior = priors,
                 sample_prior = 'only',
                 iter = 5000,
                 warmup = 1000,
                 chains = 3, cores =  3,
                 thin = 5,
                 backend = 'cmdstan',
                 seed =  123,
                 control =  list(adapt_delta =  0.99, max_treedepth = 20),
                 refresh = 0)


## ----additionalfitModel8b, results='markdown', eval=TRUE, mhidden=TRUE--------
wq.brm8 |> conditional_effects(effects = "dt_num") |>  plot(points=TRUE)
wq.brm8 |> conditional_effects(effects = "dt_num") |>  plot(points=TRUE) |>
  _[[1]] +
  scale_y_log10()


## ----additionalfitModel8c, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
wq.brm9 <- update(wq.brm8, sample_prior = "yes", cores =  3, refresh =  100)


## ----additionalfitModel8d, results='markdown', eval=TRUE, mhidden=TRUE--------
wq.brm9 |> conditional_effects(effects = "dt_num") |>  plot(points=TRUE)
wq.brm9 |> conditional_effects(effects = "dt_num")

wq.brm9 |>
    conditional_effects(effects = "dt_num") |>
    plot(points = TRUE) |>
  _[[1]] +
  geom_point(data = wq_sub, inherit.aes = FALSE,
             aes(y = NOx, x = dt_num), alpha = 0.5)


## ----additionalfitModel8e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width = 8, fig.height = 4, error = TRUE----
try({
wq.brm9 |> get_variables()
wq.brm9 |> hypothesis('bs_sdt_num_1 = 0', class = '') |> plot()
wq.brm9 |> hypothesis('sds_sdt_num_1 = 0', class = '') |> plot()
wq.brm9 |> SUYR_prior_and_posterior()
})


## ----additionalmodelValidation8a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=10----
wq.resids <- make_brms_dharma_res(wq.brm9, integerResponse = FALSE)
wrap_elements(~testUniformity(wq.resids)) +
  wrap_elements(~plotResiduals(wq.resids, quantreg = TRUE)) +
  wrap_elements(~testDispersion(wq.resids)) +
  plot_layout(nrow =  2)


## ----additionalfitModel10a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
wq.form <- bf(NOx | cens(NOx_flag) ~ s(scale(dt_num)) +
                s(scale(Mnth), bs = "cc", k = 6) + (1|reef.alias),
  family = Gamma(link = "log"))
get_prior(wq.form, data = wq_sub)
priors <- prior(normal(-3.1, 1),  class='Intercept') +
  prior(normal(0, 2), class='b') +
  prior(gamma(0.01, 0.01),  class='shape') +
  prior(student_t(3, 0, 10), class =  "sds") +
  prior(student_t(3, 0, 0.5), class =  "sd")

wq.brm10 <- brm(wq.form,
                 data = wq_sub,
                 knots = list(Mnth = seq(1, 12,len = 6)),
                 prior = priors,
                 sample_prior = 'only',
                 iter = 5000,
                 warmup = 1000,
                 chains = 3, cores =  3,
                 thin = 5,
                 backend = 'cmdstan',
                 seed =  123,
                 control =  list(adapt_delta =  0.99, max_treedepth = 20),
                 refresh = 0)


## ----additionalfitModel10b, results='markdown', eval=TRUE, mhidden=TRUE-------
wq.brm10 |> conditional_effects(effects = "dt_num") |>  plot(points=TRUE)
wq.brm10 |> conditional_effects(effects = "dt_num") |>  plot(points=TRUE) |>
  _[[1]] +
  scale_y_log10()


## ----additionalfitModel10c, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
wq.brm11 <- update(wq.brm10, sample_prior = "yes", cores =  3, refresh =  100)


## ----additionalfitModel10d, results='markdown', eval=TRUE, mhidden=TRUE-------
wq.brm11 |> conditional_effects(effects = "dt_num") |>  plot(points=TRUE)
wq.brm11 |> conditional_effects(effects = "dt_num", spaghetti = TRUE, ndraws =  250) |>  plot()
wq.brm11 |> conditional_effects(effects = "Mnth") |>  plot(points=TRUE)
wq.brm11 |> conditional_effects(effects = "Mnth") |>  plot()


## ----additionalfitModel10e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width = 8, fig.height = 4, error = TRUE----
try({
wq.brm11 |> get_variables()
wq.brm11 |> hypothesis('bs_sdt_num_1 = 0', class = '') |> plot()
wq.brm11 |> hypothesis('sds_sdt_num_1 = 0', class = '') |> plot()
wq.brm11 |> SUYR_prior_and_posterior()
})


## ----additionalmodelValidation10a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=10----
wq.resids <- make_brms_dharma_res(wq.brm11, integerResponse = FALSE)
wrap_elements(~testUniformity(wq.resids)) +
  wrap_elements(~plotResiduals(wq.resids, quantreg = TRUE)) +
  wrap_elements(~testDispersion(wq.resids)) +
  plot_layout(nrow =  2)


## ----Derivatives1a, results='markdown', eval=TRUE, mhidden=TRUE---------------


## if based on highest value
newdata <- with(wq_sub, data.frame(dt_num = seq(min(2010), max(2017), length = 1000),
  reef.alias = NA, Mnth = 5))
wq.peak <-
    add_epred_draws(
        object = wq.brm11, newdata = newdata,
        re_formula = NA,
        ndraws = 1000
    ) |>
    ungroup() |>
    group_by(.draw) |>
  summarise(x = dt_num[which.max(.epred)],
            Nox = .epred[which.max(.epred)]) |>
    ungroup() |>
  pivot_longer(cols = c(x, Nox), names_to = "variable", values_to = "value") |>
    group_by(variable) |>
  median_hdci(.width = 0.95) |>
  as.data.frame()

## based on gradient of 0
newdata <- with(wq_sub, data.frame(dt_num = seq(min(2012), max(2015), length = 1000),
  reef.alias = NA, Mnth = 5))
wq.peak <-
  add_epred_draws(object = wq.brm11, newdata = newdata,
    re_formula = NA ,
    ndraws = 1000) |>
  ungroup() |>
  group_by(.draw) |>
   #summarise(x = x[which.max(.epred)]) |>
  mutate(diff = .epred - lag(.epred)) |>
  summarise(dt_num = dt_num[which.min(abs(diff))]) |>
  median_hdci(dt_num, .width = 0.95)
wq.peak

## ## lets plot this
## data_gam.preds <-
##   data_gam.brm11 |>
##   add_epred_draws(newdata = newdata, object = _) |>
##   ungroup() |>
##   dplyr::select(-.row, -.chain, -.iteration) |>
##   group_by(x) |>
##   summarise_draws(median, HDInterval::hdi) |>
##   ungroup() |>
##   mutate(Flag = between(x, data_gam.peak$.lower, data_gam.peak$.upper),
##     Grp = data.table::rleid(Flag)
##     )
## data_gam.preds |> head()

## ggplot(data_gam.preds, aes(y = median, x = x)) +
##   geom_line(aes(colour = Flag, group = Grp))+
##   geom_ribbon(aes(ymin = lower, ymax = upper, fill = Flag, group = Grp), alpha = 0.2)


## ----additionalfitModel12a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
wq.form <- bf(NOx | cens(NOx_flag) ~ s(scale(dt_num), by =  Region) +
                s(scale(Mnth), bs = "cc", k = 6, by = Region) +
                (1|reef.alias),
  family = Gamma(link = "log"))
get_prior(wq.form, data = wq_sub)
priors <- prior(normal(-3.1, 1),  class='Intercept') +
  prior(normal(0, 2), class='b') +
  prior(gamma(0.01, 0.01),  class='shape') +
  prior(student_t(3, 0, 10), class =  "sds") +
  prior(student_t(3, 0, 1), class =  "sd")

wq.brm14 <- brm(wq.form,
                 data = wq_sub,
                 knots = list(Mnth = seq(1, 12,len = 6)),
                 prior = priors,
                 sample_prior = 'yes',
                 iter = 5000,
                 warmup = 1000,
                 chains = 3, cores =  3,
                 thin = 5,
                 backend = 'cmdstan',
                 seed =  123,
                 control =  list(adapt_delta =  0.99, max_treedepth = 20),
                 refresh = 0)


## ----additionalfitModel12d, results='markdown', eval=TRUE, mhidden=TRUE-------
wq.brm14 |> conditional_effects(effects = "dt_num:Region") |>  plot(points=TRUE)
wq.brm14 |> conditional_effects(effects = "dt_num:Region", spaghetti = TRUE, ndraws =  250) |>  plot()
wq.brm14 |> conditional_effects(effects = "Mnth:Region") |>  plot(points=TRUE)
wq.brm14 |> conditional_effects(effects = "Mnth:Region", spaghetti =  TRUE, draws =  250) |>  plot()


## ----Derivatives2a, results='markdown', eval=FALSE, mhidden=TRUE--------------
# 
# 
# ## if based on highest value
# newdata <- with(wq_sub, data.frame(dt_num = seq(min(2010), max(2017), length = 1000),
#   reef.alias = NA, Mnth = 5), Region = "Townsville")
# wq.peak <-
#     add_epred_draws(
#         object = wq.brm14, newdata = newdata,
#         re_formula = NA,
#         ndraws = 1000
#     ) |>
#     ungroup() |>
#     group_by(.draw, Region) |>
#   summarise(x = dt_num[which.max(.epred)],
#             Nox = .epred[which.max(.epred)]) |>
#     ungroup() |>
#   pivot_longer(cols = c(x, Nox), names_to = "variable", values_to = "value") |>
#     group_by(Region, variable) |>
#   median_hdci(.width = 0.95) |>
#   as.data.frame()
# 
# ## based on gradient of 0
# newdata <- with(wq_sub, data.frame(dt_num = seq(min(2012), max(2015), length = 1000),
#   reef.alias = NA, Mnth = 5), Region = "Townsville")
# wq.peak <-
#   add_epred_draws(object = wq.brm14, newdata = newdata,
#     re_formula = NA ,
#     ndraws = 1000) |>
#   ungroup() |>
#   group_by(.draw) |>
#    #summarise(x = x[which.max(.epred)]) |>
#   mutate(diff = .epred - lag(.epred)) |>
#   summarise(dt_num = dt_num[which.min(abs(diff))]) |>
#   median_hdci(dt_num, .width = 0.95)
# wq.peak
# 
# ## ## lets plot this
# ## data_gam.preds <-
# ##   data_gam.brm11 |>
# ##   add_epred_draws(newdata = newdata, object = _) |>
# ##   ungroup() |>
# ##   dplyr::select(-.row, -.chain, -.iteration) |>
# ##   group_by(x) |>
# ##   summarise_draws(median, HDInterval::hdi) |>
# ##   ungroup() |>
# ##   mutate(Flag = between(x, data_gam.peak$.lower, data_gam.peak$.upper),
# ##     Grp = data.table::rleid(Flag)
# ##     )
# ## data_gam.preds |> head()
# 
# ## ggplot(data_gam.preds, aes(y = median, x = x)) +
# ##   geom_line(aes(colour = Flag, group = Grp))+
# ##   geom_ribbon(aes(ymin = lower, ymax = upper, fill = Flag, group = Grp), alpha = 0.2)

