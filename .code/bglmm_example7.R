## -----------------------------------------------------------------------------
#| label: setup
#| include: false
#| cache: false

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

library(car)       #for regression diagnostics
library(broom)     #for tidy output
library(ggfortify) #for model diagnostics
library(knitr)     #for kable
library(ggeffects)
library(emmeans)   #for estimating marginal means
library(MASS)      #for glm.nb
library(tidyverse) #for data wrangling
library(brms)
library(tidybayes)
library(broom.mixed)
library(rstan)
library(cmdstanr)
library(patchwork)
library(DHARMa)
library(easystats)
library(modelsummary)
source("helperFunctions.R")


## -----------------------------------------------------------------------------
#| label: readData
mullens <- read_csv('../data/mullens.csv', trim_ws=TRUE)


## -----------------------------------------------------------------------------
#| label: examinData
glimpse(mullens)


## -----------------------------------------------------------------------------
## Explore the first 6 rows of the data
head(mullens)


## -----------------------------------------------------------------------------
str(mullens)


## -----------------------------------------------------------------------------
mullens |> datawizard::data_codebook()


## -----------------------------------------------------------------------------
mullens |> modelsummary::datasummary_skim()
mullens |> modelsummary::datasummary_skim(by = c("BREATH", "O2LEVEL"))


## ----dataPreparation, results='markdown', eval=TRUE, mhidden=FALSE------------
mullens <- mullens |>
  mutate(BREATH = factor(BREATH),
         TOAD = factor(TOAD),
         pBUC = FREQBUC/100,
         pzBUC = ifelse(pBUC == 0,0.01, pBUC))


## ----eda1a, results='markdown', eval=TRUE, mhidden=TRUE-----------------------
ggplot(mullens,aes(y=FREQBUC, x=factor(O2LEVEL), color=BREATH)) +
    geom_boxplot()


## ----eda1b, results='markdown', eval=TRUE, mhidden=TRUE-----------------------
ggplot(mullens,aes(y=pzBUC, x=O2LEVEL, color=BREATH)) +
    geom_smooth() + geom_point()



## ----eda1c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=7, fig.height=7----
ggplot(mullens,aes(y=pzBUC, x=O2LEVEL, color=BREATH)) +
  geom_smooth() + geom_point() +
  facet_wrap(~BREATH+TOAD, scales='free')
  #facet_grid(TOAD~BREATH)


## ----fitModel2a, results='markdown', eval=TRUE, mhidden=TRUE------------------
mullens.form <- bf(pBUC ~ BREATH*poly(O2LEVEL,3) + (1|TOAD),
  family=zero_one_inflated_beta())
get_prior(mullens.form, data =  mullens)



## ----fitModel3a, results='markdown', eval=TRUE, mhidden=TRUE------------------
mullens |> group_by(BREATH) |>
    summarise(logit(median(pBUC)),
              logit(mad(pBUC)))
median(logit(mullens$pBUC))
mad(logit(mullens$pBUC))

sd(logit(mullens$pBUC))/apply(model.matrix(~BREATH*poly(O2LEVEL,3), data = mullens), 2, sd)

standist::visualize('gamma(0.01, 0.01)', 'gamma(2,1)', 'gamma(1,1)', xlim=c(0,10))
standist::visualize('beta(1,1)', xlim=c(0,1))
priors <- prior(normal(-2, 2), class='Intercept') +
    prior(normal(0,1), class='b') +
    ## prior(normal(0,1), class='b', coef = "BREATHlung") +
    #prior(gamma(2,1), class='sd') +
    prior(student_t(3, 0, 2), class='sd') +
    prior(gamma(0.01, 0.01), class='phi') +
    ## prior(gamma(2, 1), class='phi') +
    prior(beta(1,1), class='zoi') +
    prior(beta(1,1), class='coi')

mullens.form <- bf(pBUC ~ BREATH*poly(O2LEVEL,3) + (1|TOAD),
  family=zero_one_inflated_beta())

mullens.brm <- brm(mullens.form,
  data=mullens,
  prior = priors,
  sample_prior = 'only',
  iter=5000,
  warmup=2500,
  thin=5,
  chains=3, cores=3,
  refresh =  0,
  seed = 123,
  control=list(adapt_delta=0.99, max_treedepth = 20),
  backend =  'cmdstan')


## ----fitModel3b, results='markdown', eval=TRUE, mhidden=TRUE------------------
mullens.brm |> conditional_effects(effects = "O2LEVEL:BREATH") |> plot(points =  TRUE)


## ----fitModel3c, results='markdown', eval=TRUE, mhidden=TRUE------------------
mullens.brm2 <- update(mullens.brm, sample_prior = "yes", refresh =  0,
  cores = 3, seed =  123)
priors <- prior(normal(-2, 2), class='Intercept') +
    prior(normal(0,2), class='b') +
    prior(student_t(3, 0, 2), class='sd') +
    prior(gamma(0.01, 0.01), class='phi') +
    prior(beta(1,1), class='zoi') +
    prior(beta(1,1), class='coi')
mullens.brm3 <- update(mullens.brm2, prior =  priors, refresh =  0,
  cores = 3, seed =  123)


## ----fitModel3d, results='markdown', eval=TRUE, mhidden=TRUE------------------
mullens.brm2 |> conditional_effects(effects = "O2LEVEL:BREATH") |> plot(points =  TRUE)
mullens.brm2 |> ggpredict(~O2LEVEL|BREATH) |> plot(show_data = TRUE, jitter =  FALSE)
mullens.brm3 |> conditional_effects(effects = "O2LEVEL:BREATH") |> plot(points =  TRUE)


## ----fitModel3e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width = 10, fig.height=8----
mullens.brm2 |> SUYR_prior_and_posterior()
mullens.brm3 |> SUYR_prior_and_posterior()


## ----modelValidation1a, results='markdown', eval=FALSE, mhidden=TRUE----------
# mullens.brm3$fit |> stan_trace()
# mullens.brm3$fit |> stan_ac()
# mullens.brm3$fit |> stan_rhat()
# mullens.brm3$fit |> stan_ess()


## ----modelValidation2a, results='markdown', eval=FALSE, mhidden=TRUE----------
# mullens.brm3 |> pp_check(type = 'dens_overlay', ndraws=200)
# 
# mullens.brm3 |> pp_check(group="BREATH", type='intervals_grouped')
# mullens.brm3 |> pp_check(group="O2LEVEL", type='intervals_grouped')


## ----modelValidation2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=10----
mullens.resids <- make_brms_dharma_res(mullens.brm3, integerResponse = FALSE)
wrap_elements(~testUniformity(mullens.resids)) +
               wrap_elements(~plotResiduals(mullens.resids, form = factor(rep(1, nrow(mullens))))) +
               wrap_elements(~plotResiduals(mullens.resids, quantreg = TRUE)) +
               wrap_elements(~testDispersion(mullens.resids))


## ----partialPlot2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mullens.brm3 |> conditional_effects(effects = "O2LEVEL:BREATH") |> plot(points =  TRUE)


## ----partialPlot2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mullens.brm3 |> ggpredict(~O2LEVEL|BREATH) |> plot(show_data = TRUE, jitter =  FALSE)


## ----partialPlot2f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mullens.brm3 |> ggemmeans(~O2LEVEL|BREATH) |> plot(show_data = TRUE, jitter =  FALSE)

newdata <- with(mullens, list(
  O2LEVEL = seq(min(O2LEVEL), max(O2LEVEL), len =  100),
  BREATH =  levels(BREATH)))

mullens.em <-
  mullens.brm3 |>
  emmeans(~O2LEVEL|BREATH, at = newdata, type = "response") |>
  as.data.frame()
mullens.em |>
  ggplot(aes(y = response, x = O2LEVEL, colour = BREATH, fill = BREATH)) +
  geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD),
    alpha =  0.2, colour = NA) +
  geom_line()


## ----summariseModel2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mullens.brm2 |> summary()


## ----summariseModel2m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mullens.brm3 |>
  as_draws_df() |>
  mutate(across(everything(), exp)) |>
  dplyr::select(matches("^b_.*|^sd_.*")) |>
  summarise_draws(
    median,
    HDInterval::hdi,
    rhat,
    length,
    ess_bulk, ess_tail,
    Pl = ~ mean(.x < 1),
    Pg = ~ mean(.x > 1)
  ) |>
  knitr::kable()

mean(mullens$O2LEVEL)
0.163/(1+0.163)

mullens.brm3 |>
  as_draws_df() |>
  dplyr::select(matches("^b_.*")) |>
  exp() |>
  summarise_draws(
    median,
    HDInterval::hdi,
    rhat,
    length,
    ess_bulk, ess_tail,
    Pl = ~ mean(.x < 1),
    Pg = ~ mean(.x > 1)
  )
mullens.brm2 |> performance::r2()
mullens.brm2 |> performance::r2_posterior() |> as.data.frame() |>
  median_hdci()


## ----furtherInvestigations1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mullens.brm3 |> emtrends(specs = 'BREATH',  var = 'O2LEVEL',  max.degree = 3)
mullens.brm3 |> emtrends(specs = 'BREATH',  var = 'O2LEVEL',  max.degree = 3) |> gather_emmeans_draws() |>
  summarise(median_hdci(.value),
    Pl = mean(.value < 0),
    Pg =  mean(.value > 0)
    )


## ----furtherInvestigations1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mullens.brm3 |> emmeans(~BREATH|O2LEVEL, at =  list(O2LEVEL = c(0, 21, 50))) |>
  pairs() |>
  gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  summarise(median_hdci(.value),
    Pl = mean(.value < 1),
    Pg =  mean(.value > 1)
    )


mullens.grid <- with(mullens,
  list(O2LEVEL = seq(min(O2LEVEL), max(O2LEVEL), len = 1000),
    BREATH = levels(BREATH)))
## mullens.grid <- with(mullens,
##   list(O2LEVEL = modelr::seq_range(O2LEVEL, n = 1000),
##     BREATH = levels(BREATH)))

## At what point does the evidence for the curves being different stop
mullens.brm3 |> emmeans(~BREATH|O2LEVEL, at =  mullens.grid) |>
  pairs() |>
  gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  summarise(median_hdci(.value),
    Pl = mean(.value < 1),
    Pg =  mean(.value > 1)
  ) |>
  filter(Pg < 0.90) |>
  slice(1:3)

## Find the peak
newdata <- mullens.brm3 |>
  emmeans(~O2LEVEL|BREATH,  at = mullens.grid) |>
  as.data.frame()
newdata |> group_by(BREATH) |>
  summarise(value = O2LEVEL[which.max(emmean)])
## With uncertainty
mullens.brm3 |>
  emmeans(~O2LEVEL|BREATH,  at = mullens.grid) |>
  gather_emmeans_draws() |>
  group_by(.draw, BREATH) |>
  summarise(value = O2LEVEL[which.max(.value)]) |>
  ungroup() |>
  group_by(BREATH) |>
  summarise(median_hdci(value))


# Or via Derivatives
## mullens.brm3 |> estimate_slopes(trend =  "O2LEVEL",
##   at = c("BREATH='lung'","O2LEVEL"))

## mullens.brm2 |>
##   emtrends(~O2LEVEL|BREATH, var = "O2LEVEL",
##     cov.red = \(x) seq(min(x), max(x), length = 10)
##   ) |>
##   summary() |>
##   ggplot(aes(y = O2LEVEL.trend, x = O2LEVEL, colour =  BREATH)) +
##   geom_hline(yintercept = 0, linetype = "dashed") +
##   geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD, fill = BREATH), alpha = 0.2) +
##   geom_line()


## ----summaryFigure1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mullens.grid <- with(mullens,
   list(BREATH = levels(BREATH),
     O2LEVEL = modelr::seq_range(O2LEVEL, n=100)
     )
)
newdata <- mullens.brm2 |>
    emmeans(~O2LEVEL|BREATH, at = mullens.grid, type = 'response') |>
    as.data.frame()
head(newdata)
ggplot() +
    geom_ribbon(data = newdata,
                aes(ymin = lower.HPD, ymax = upper.HPD,
                    x = O2LEVEL, fill = BREATH), alpha = 0.3)+
    geom_line(data = newdata,
              aes(y = response, x = O2LEVEL, color = BREATH)) +
    scale_y_continuous('Buccal breathing rate', labels = function(x) 100*x) +
    theme_classic()



## ----old, results='markdown', eval=FALSE, mhidden=TRUE------------------------
# ##prior_summary(mullens.brm)
# 
# 
# ## pars <- mullens.brm |> get_variables()
# ## wch <- grepl('^b_.*|^sd_.*|phi', pars, perl=TRUE)
# 
# ## g <- vector('list', length=sum(wch)-1)
# ## names(g) <- pars[wch][-1]
# ## for (i in pars[wch]) {
# ##     print(i)
# ##     if (i == 'b_Intercept') next
# ##     p <- mullens.brm |> hypothesis(paste0(i,'=0'), class='') |> plot()
# ##     g[[i]] <- p[[1]]
# ## }
# ## patchwork::wrap_plots(g)
# 
# stan_trace(mullens.brm$fit, pars = pars[wch])
# stan_ac(mullens.brm$fit, pars = pars[wch])
# stan_rhat(mullens.brm$fit, pars = pars[wch])
# stan_rhat(mullens.brm$fit)
# stan_ess(mullens.brm$fit)
# 
# 
# preds <- posterior_predict(mullens.brm,  nsamples=250,  summary=FALSE)
# mullens.resids <- createDHARMa(simulatedResponse = t(preds),
#                             observedResponse = mullens$pBUC,
#                             fittedPredictedResponse = apply(preds, 2, median),
#                             integerResponse = FALSE)
# plot(mullens.resids)
# testDispersion(mullens.resids)

