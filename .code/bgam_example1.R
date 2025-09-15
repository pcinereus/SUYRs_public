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


## ----readData, results='markdown', eval=TRUE----------------------------------
data_gam <- read_csv("../data/data_gam.csv", trim_ws = TRUE)


## -----------------------------------------------------------------------------
#| label: examinData
#| dependson: readData

data_gam |> glimpse()


## -----------------------------------------------------------------------------
#| label: headData
#| dependson: readData
## Explore the first 6 rows of the data
data_gam |> head()


## -----------------------------------------------------------------------------
#| label: strData
#| dependson: readData
data_gam |> str()


## -----------------------------------------------------------------------------
#| label: easyData
#| dependson: readData
data_gam |> datawizard::data_codebook()


## ----EDA1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=4, fig.height=4----
ggplot(data_gam, aes(y=y, x=x))+
    geom_point()+
    geom_line()


## ----EDA1b, results='markdown', eval=TRUE, mhidden=TRUE,warning=FALSE,message=FALSE, fig.width=4, fig.height=4----
ggplot(data_gam, aes(y=y, x=x))+
    geom_point()+
    geom_smooth()


## ----EDA1c, results='markdown', eval=TRUE, mhidden=TRUE,warning=FALSE,message=FALSE, fig.width=4, fig.height=4----
ggplot(data_gam, aes(y=y, x=x))+
    geom_point()+
    geom_smooth(method='lm')


## ----EDA1d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=4, fig.height=4----
ggplot(data_gam, aes(y=y, x=x))+
    geom_point()+
    geom_smooth(method='gam', formula=y~s(x,k=3))


## ----fitModel1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=4, fig.height=4----
data.frame(smoothCon(s(x, k=3),  data=data_gam)[[1]]$X) %>%
  bind_cols(data_gam)


## ----fitModel2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=4, fig.height=4----
basis(s(x, k=3),  data=data_gam) %>% draw()
basis(s(x, k=3, bs='cr'),  data=data_gam) %>% draw()


## ----fitModel3, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=4, fig.height=4----
newdata <-
    data.frame(smoothCon(s(x, k=3),  data=data_gam)[[1]]$X) %>%
    bind_cols(data_gam)
ggplot(newdata,  aes(x=x)) +
    geom_line(aes(y=X1)) +
    geom_line(aes(y=X2)) +
    geom_line(aes(y=X3))


## ----fitModel1a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE, error = TRUE----
try({
data_gam.form <- bf(y ~ s(x), family = gaussian())
data_gam.brm <- brm(data_gam.form,
               data = data_gam,
               iter = 5000,
               warmup = 1000,
               chains = 3,
               thin = 5,
               refresh = 0,
               backend = 'cmdstan')
})


## ----fitModel1b, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
data_gam.form <- bf(y ~ s(x, k = 3), family = gaussian())
get_prior(data_gam.form, data =  data_gam)
data_gam.brm <- brm(data_gam.form,
               data = data_gam,
               iter = 5000,
               warmup = 1000,
               chains = 3, cores = 3,
               thin = 5,
               refresh = 0,
               backend = 'cmdstan')

## ----fitModel1c, results='markdown', eval=TRUE, mhidden=TRUE, paged.print=FALSE,tidy.opts = list(width.cutoff = 80), echo=2----
options(width=100)
prior_summary(data_gam.brm)
options(width=80)


## ----fitModel1d, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
data_gam.brm <- brm(data_gam.form,
                 data = data_gam,
                 prior = prior(normal(0, 2.5), class = 'b'),
                 sample_prior = 'only',
                 iter = 5000,
                 warmup = 1000,
                 chains = 3,
                 thin = 5,
                 backend = 'cmdstan',
                 refresh = 0)


## ----fitModel1d2, results='markdown', eval=TRUE, mhidden=TRUE-----------------
data_gam.brm |> conditional_effects() |>  plot(points=TRUE)


## ----fitModel2a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
data_gam |> summarise(median(y), mad(y))
priors <- prior(normal(5, 1.5),  class='Intercept') +
  prior(normal(0, 1.5), class='b') +
  prior(student_t(3, 0, 1.5),  class='sigma') +
  prior(student_t(3, 0, 10), class =  "sds")

data_gam.form <- bf(y ~ s(x, k = 3))
data_gam.brm2 <- brm(data_gam.form,
                 data = data_gam,
                 prior = priors,
                 sample_prior = 'only',
                 iter = 5000,
                 warmup = 1000,
                 chains = 3, cores =  3,
                 thin = 5,
                 backend = 'cmdstan',
                 control =  list(adapt_delta =  0.99),
                 refresh = 0)


## ----fitModel2b, results='markdown', eval=TRUE, mhidden=TRUE------------------
data_gam.brm2 |> conditional_effects() |>  plot(points=TRUE)
data_gam.brm2 |> conditional_smooths() |>  plot()


## ----fitModel2c, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
data_gam.brm3 <- update(data_gam.brm2, sample_prior = "yes", cores =  3, refresh =  0)


## ----fitModel2d, results='markdown', eval=TRUE, mhidden=TRUE------------------
data_gam.brm3 |> conditional_effects() |>  plot(points=TRUE)
data_gam.brm3 |> conditional_effects(spaghetti = TRUE, ndraws =  200) |>  plot(points=TRUE)


## ----fitModel2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width = 8, fig.height = 4, error = TRUE----
try({
data_gam.brm3 |> get_variables()
data_gam.brm3 |> hypothesis('bs_sx_1 = 0', class = '') |> plot()
data_gam.brm3 |> hypothesis('sds_sx_1 = 0', class = '') |> plot()
data_gam.brm3 |> hypothesis('sigma = 0', class = '') |> plot()
data_gam.brm3 |> SUYR_prior_and_posterior()
})


## ----modelValidation2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
data_gam.brm3$fit |> stan_trace()


## ----modelValidation2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
data_gam.brm3$fit |> stan_ac()


## ----modelValidation2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
data_gam.brm3$fit |> stan_rhat()

## ----modelValidation2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
data_gam.brm3$fit |> stan_ess()


## ----modelValidation3a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
data_gam.brm3 |> pp_check( type='dens_overlay', ndraws=100)


## ----modelValidation3b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=10----

data_gam.resids <- make_brms_dharma_res(data_gam.brm3, integerResponse = FALSE)
wrap_elements(~testUniformity(data_gam.resids)) +
               wrap_elements(~plotResiduals(data_gam.resids, form = factor(rep(1, nrow(data_gam))))) +
               wrap_elements(~plotResiduals(data_gam.resids, quantreg = FALSE)) +
               wrap_elements(~testDispersion(data_gam.resids))
testDispersion(data_gam.resids)


## ----partialPlot1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
data_gam.brm3 |> conditional_effects() |> plot(points = TRUE)
data_gam.brm3 |>
    conditional_effects(spaghetti = TRUE, ndraws = 250) |>
    plot(points = TRUE)


## ----summariseModel1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
data_gam.brm3 |> summary()


## ----summariseModel1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
data_gam.brm3 |> get_variables()
data_gam.brm3 |>
  as_draws_df() |>
  dplyr::select(matches("^b_.*|^bs.*|^sds.*|^sigma$|^s_s.*")) |>
  summarise_draws(median,
    HDInterval::hdi,
    Pl = ~mean(.x < 0),
    Pg = ~mean(.x > 0)
    )


## ----aaa, results='markdown', eval=TRUE---------------------------------------
newdata <- with(data_gam, data.frame(x = c(min(x), 9)))
add_epred_draws(object = data_gam.brm3, newdata = newdata,
  ndraws =  2400) |>
  ungroup() |>
  group_by(.draw) |>
  summarise(Diff =  diff(.epred)) |>
  summarise(median_hdci(Diff),
    Pl = mean(Diff < 0),
    Pg =  mean(Diff > 0))



## ----Derivatives1a, results='markdown', eval=TRUE, mhidden=TRUE---------------
newdata <- with(data_gam, data.frame(x = seq(min(x), max(x), length = 1000)))
data_gam.peak <-
  add_epred_draws(object = data_gam.brm3, newdata = newdata, ndraws = 1000) |>
  ungroup() |>
  group_by(.draw) |>
   #summarise(x = x[which.max(.epred)]) |>
  mutate(diff = .epred - lag(.epred)) |>
  summarise(x = x[which.min(abs(diff))]) |>
  median_hdci(x, .width = 0.95)
data_gam.peak

## lets plot this
data_gam.preds <-
  data_gam.brm3 |>
  add_epred_draws(newdata = newdata, object = _) |>
  ungroup() |>
  dplyr::select(-.row, -.chain, -.iteration) |>
  group_by(x) |>
  summarise_draws(median, HDInterval::hdi) |>
  ungroup() |>
  mutate(Flag = between(x, data_gam.peak$.lower, data_gam.peak$.upper),
    Grp = data.table::rleid(Flag)
    )
data_gam.preds |> head()

ggplot(data_gam.preds, aes(y = median, x = x)) +
  geom_line(aes(colour = Flag, group = Grp))+
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = Flag, group = Grp), alpha = 0.2)


## ----Derivatives1b, results='markdown', eval=TRUE, mhidden=TRUE---------------
data_gam.brm3 |>
  estimate_relation(keep_iterations = TRUE, length = 1000) |>
  estimate_smooth(x = "x")


## ----Derivatives1c, results='markdown', eval=TRUE, mhidden=TRUE---------------
newdata <- with(data_gam, data.frame(x = seq(min(x), max(x), length = 1000)))
data_gam.brm3 |>
  add_epred_draws(newdata = newdata, object = _) |>
  ungroup() |>
  group_by(.draw) |>
  mutate(diff = .epred - lag(.epred)) |>
  summarise(
    maxGrad = max(abs(diff), na.rm = TRUE),
    x = x[which.max(diff)]) |>
  summarise_draws(median, HDInterval::hdi)



## ----Derivatives1d, results='markdown', eval=TRUE, mhidden=TRUE---------------
newdata <- with(data_gam, data.frame(x = seq(min(x), max(x), length = 1000)))
data_gam.brm3 |>
  add_epred_draws(newdata = newdata, object = _) |>
  filter(x > 3, x <13) |>
  ungroup() |>
  group_by(.draw) |>
  mutate(diff = .epred - lag(.epred),
    diff2 = diff - lag(diff)) |>
  summarise(
    maxChange = max(abs(diff2), na.rm = TRUE),
    x = x[which.max(diff)]) |>
  summarise_draws(median, HDInterval::hdi)

