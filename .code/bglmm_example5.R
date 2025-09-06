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
library(tidyverse) #for data wrangling
library(car)       #for regression diagnostics
library(broom)     #for tidy output
library(ggfortify) #for model diagnostics
library(knitr)     #for kable
library(emmeans)   #for estimating marginal means
library(MASS)      #for glm.nb
library(brms)
library(broom.mixed)
library(tidybayes)
library(bayesplot)
library(standist)   #for visualizing distributions
library(rstanarm)
library(cmdstanr)   #for cmdstan
library(ggeffects)
library(rstan)
library(DHARMa)
library(ggridges)
library(easystats)     #framework for stats, modelling and visualisation
library(patchwork)
library(modelsummary)
library(geoR)
source('helperFunctions.R')


## ----readData, results='markdown', eval=TRUE----------------------------------
owls <- read_csv("../data/owls.csv", trim_ws = TRUE)


## -----------------------------------------------------------------------------
#| label: examinData
glimpse(owls)


## -----------------------------------------------------------------------------
## Explore the first 6 rows of the data
head(owls)


## -----------------------------------------------------------------------------
str(owls)


## -----------------------------------------------------------------------------
owls |> datawizard::data_codebook()


## -----------------------------------------------------------------------------
owls |> modelsummary::datasummary_skim()
owls |> modelsummary::datasummary_skim(by = c("SexParent", "FoodTreatment"))


## ----dataProcessing, results='markdown', eval=TRUE, mhidden=TRUE--------------
## Amount of Sibling negotiation (vocalizations when parents are absent)
## Foot treatment (deprived or satiated
## Sex of parent
## Arrival time of parent
## Nest as random
## Brood size offset
owls <- owls |> mutate(Nest = factor(Nest),
                       FoodTreatment = factor(FoodTreatment),
                       SexParent = factor(SexParent),
                       NCalls = SiblingNegotiation)


## ----eda1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=7, fig.height=5----
ggplot(data = owls, aes(y = NCalls, x = FoodTreatment,  color=SexParent)) +
  geom_violin()
ggplot(data = owls, aes(y = NCalls, x = FoodTreatment,  color=SexParent)) +
  geom_violin() +
  geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width=0.9))
ggplot(data = owls, aes(y = NCalls, x = FoodTreatment,  color=SexParent)) +
  geom_violin() +
  geom_point(position=position_jitterdodge(jitter.height=0,  dodge.width=1))+
  scale_y_continuous(trans=scales::pseudo_log_trans())


## ----eda2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=10, fig.height=10----
ggplot(data=owls) +
  geom_point(aes(y=NCalls,  x=FoodTreatment,  color=SexParent),  position=position_dodge(0.5)) +
  facet_wrap(~Nest)
ggplot(data=owls) +
  geom_violin(aes(y = NCalls, x = FoodTreatment,  color=SexParent)) +
  geom_point(aes(y=NCalls,  x=FoodTreatment,  color=SexParent),  position=position_dodge(0.5)) +
  facet_wrap(~Nest)


## ----eda3, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=10, fig.height=5----
ggplot(data = owls,aes(y = NCalls, x = BroodSize, color=SexParent)) +
  geom_point() +
  geom_smooth(method='lm') +
  facet_grid(~FoodTreatment) +
  scale_y_continuous(trans=scales::pseudo_log_trans()) +
  scale_x_log10()


## ----fitModel1a, results='markdown', eval=TRUE, mhidden=TRUE------------------
owls_rstanP <- stan_glmer(NCalls ~ FoodTreatment*SexParent +
                           offset(log(BroodSize)) + (1|Nest),
                          data = owls,
                          family = poisson(link = 'log'),
                          refresh = 0,
                          iter = 5000,
                          warmup = 2000,
                          thin = 10,
                          chains = 3,
                          cores = 3)


## ----fitModel1b, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
owls_rstanP |> prior_summary()


## ----fitModel1c, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
2.5/sd(owls$NCalls)


## ----fitModel1d, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
model.matrix(~FoodTreatment*SexParent+offset(log(BroodSize)), data = owls) |>
    apply(2, sd) * 1/2.5


## ----fitModel1f, results='markdown', eval=TRUE, mhidden=TRUE------------------
owls_rstanarmP1 <- update(owls_rstanP,  prior_PD=TRUE)


## ----fitModel1g, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
owls_rstanarmP1 |>
    ggpredict(~FoodTreatment*SexParent) |>
    plot(show_data=TRUE, jitter=c(0.25,0)) +
    scale_y_continuous('', trans=scales::pseudo_log_trans())


## ----fitModel1h, results='markdown', eval=TRUE, mhidden=TRUE------------------
owls_rstanarmP2 <- stan_glmer(NCalls ~ FoodTreatment*SexParent +
                               offset(log(BroodSize)) + (1|Nest),
                           data = owls,
                           family = poisson(link = 'log'),
                           prior_intercept = normal(0, 1.5, autoscale = FALSE),
                           prior = normal(0, c(2.2,2.2,2.5), autoscale = FALSE),
                           prior_aux = exponential(1),
                           prior_covariance = decov(1, 1, 1, 1),
                           refresh = 0,
                           iter = 5000,
                           prior_PD = TRUE,
                           warmup = 2000,
                           thin = 10,
                           chains = 3,
                           cores = 3)


## ----fitModel1i, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
owls_rstanarmP2 |>
    ggpredict(~FoodTreatment*SexParent) |>
    plot(show_data = TRUE, jitter = c(0.25, 0)) +
    scale_y_continuous('', trans=scales::pseudo_log_trans())


## ----fitModel1j, results='markdown', eval=TRUE, mhidden=TRUE, dependson='fitModel1h'----
owls_rstanarmP3 <- update(owls_rstanarmP2,  prior_PD=FALSE)


## ----modelFit1k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
posterior_vs_prior(owls_rstanarmP3, color_by='vs', group_by=TRUE,
                   facet_args=list(scales='free_y'))


## -----------------------------------------------------------------------------
#| label: test
#| results: markup
#| eval: false
#| echo: false
#| cache: false

# owls_form <- bf(NCalls ~
#                     offset(log(BroodSize)) + (1|Nest),
#                 family=poisson(link='log'))
# options(width = 150)
# get_prior(owls_form, data = owls)
# ## owls |>
# ##     group_by(FoodTreatment, SexParent) |>
# ##     summarise(Mean = log(median(NCalls/BroodSize)),
# ##               SD = log(sd(NCalls/BroodSize)),
# ##               MAD = log(mad(NCalls/BroodSize)))
# ## priors <- prior(normal(0.2, 0.6), class = 'Intercept') +
# ##     prior(student_t(3, 0, 0.6), class = 'sd')
# 
# ## owls_brm2 <- brm(owls_form,
# ##                  data = owls,
# ##                  prior = priors,
# ##                  sample_prior = 'only',
# ##                  iter = 5000,
# ##                  warmup =2500,
# ##                  chains = 3,
# ##                  cores = 3,
# ##                  thin = 10,
# ##                  refresh = 0,
# ##                  seed = 123,
# ##                  control =  list(adapt_delta = 0.99),
# ##                  backend = "cmdstanr"
# ##                  )
# 
# ## owls_brm2 |>
# ##     conditional_effects("FoodTreatment:SexParent") |>
# ##   plot(points = TRUE) |>
# ##   _[[1]] +
# ##   scale_y_continuous(trans = scales::pseudo_log_trans())
# 
# ## owls_form <- bf(NCalls ~
# ##                     offset(log(BroodSize)) +
# ##                     (1|Nest),
# ##                 zi ~ 1,
# ##                 family=zero_inflated_poisson(link='log'))
# ## priors <- prior(normal(0.2, 0.6), class = "Intercept") +
# ##   prior(student_t(3, 0, 0.6), class = "sd") +
# ## prior(logistic(0,1), class = "Intercept", dpar = "zi")
# 
# ## owls_brm3 <- update(owls_brm2,
# ##                        sample_prior = 'yes',
# ##                        iter =  5000,
# ##                        warmup =  2500,
# ##                        control = list(adapt_delta = 0.99, max_treedepth = 20),
# ##                        backend =  "cmdstanr",
# ##                        refresh = 0)
# 
# ## owls_brm3 |> SUYR_prior_and_posterior()
# 
# 
# ## owls_brm3 |> pp_check(type = 'dens_overlay', ndraws = 100)
# ## owls_brm3 |> pp_check(type = 'dens_overlay', ndraws = 100) +
# ##   scale_x_log10()
# 
# ## owls_resids <- make_brms_dharma_res(owls_brm3, integerResponse = TRUE)
# ## wrap_elements(~testUniformity(owls_resids)) +
# ##                wrap_elements(~plotResiduals(owls_resids, form = factor(rep(1, nrow(owls))))) +
# ##                wrap_elements(~plotResiduals(owls_resids, quantreg = TRUE)) +
# ##                wrap_elements(~testDispersion(owls_resids))
# ## preds <- owls_brm3 |> posterior_predict(nsamples = 250,  summary = FALSE)
# 
# ## owls_resids |> testZeroInflation()
# ## owls_resids |> testDispersion()
# ## owls_resids |> testUniformity()
# ## owls_resids |> testQuantiles()
# ## owls_resids |> testResiduals()
# 
# 
# 
# ## owls |>
# ##     group_by(FoodTreatment, SexParent) |>
# ##     summarise(Mean = log(median(NCalls/BroodSize)),
# ##               SD = log(sd(NCalls/BroodSize)),
# ##               MAD = log(mad(NCalls/BroodSize)))
# ## owls_form <- bf(NCalls ~
# ##                     offset(log(BroodSize)) +
# ##                     (FoodTreatment*SexParent|Nest),
# ##                 zi ~ FoodTreatment*SexParent,
# ##                 family=zero_inflated_poisson(link='log'))
# ## priors <- prior(normal(0.2, 0.6), class = "Intercept") +
# ##   prior(student_t(3, 0, 0.6), class = "sd") +
# ## prior(logistic(0,1), class = "Intercept", dpar = "zi") +
# ##     prior(normal(0,1), class='b', dpar='zi')
# 
# ## owls_brm3 <- brm(owls_form,
# ##                  data=owls,
# ##                  prior = priors,
# ##                  sample_prior = 'yes',
# ##                  iter=5000,
# ##                  warmup=2500,
# ##                  thin=10,
# ##                  chains=3,
# ##                  refresh=0,
# ##                  backend = "cmdstanr",
# ##                  cores=3)
# ## owls_brm3 |> SUYR_prior_and_posterior()
# 
# 
# ## owls_brm3 |> pp_check(type = 'dens_overlay', ndraws = 100)
# ## owls_brm3 |> pp_check(type = 'dens_overlay', ndraws = 100) +
# ##   scale_x_log10()
# 
# ## owls_resids <- make_brms_dharma_res(owls_brm3, integerResponse = TRUE)
# ## wrap_elements(~testUniformity(owls_resids)) +
# ##                wrap_elements(~plotResiduals(owls_resids, form = factor(rep(1, nrow(owls))))) +
# ##                wrap_elements(~plotResiduals(owls_resids, quantreg = TRUE)) +
# ##                wrap_elements(~testDispersion(owls_resids))
# ## preds <- owls_brm3 |> posterior_predict(nsamples = 250,  summary = FALSE)
# 
# ## owls_resids |> testZeroInflation()
# ## owls_resids |> testDispersion()
# ## owls_resids |> testUniformity()
# ## owls_resids |> testQuantiles()
# ## owls_resids |> testResiduals()
# 
# ## priors <- prior(normal(0.4,0.7), class = 'Intercept') +
# ##     ## prior(normal(0, 2.2), class = 'b', coef = 'FoodTreatmentSatiated') +
# ##     ## prior(normal(0, 2.2), class = 'b', coef = 'SexParentMale') +
# ##     prior(normal(0, 2), class = 'b') +
# ##     prior(student_t(3, 0, 0.7), class = 'sd') +
# ##     prior(lkj_corr_cholesky(1), class = 'cor') +
# ##     prior(logistic(0,1), class='Intercept', dpar='zi') +
# ##     prior(normal(0,1), class='b', dpar='zi')
# ## owls_form <- bf(NCalls ~ FoodTreatment*SexParent +
# ##                     offset(log(BroodSize)) +
# ##                     (FoodTreatment*SexParent|Nest),
# ##                 zi ~ FoodTreatment*SexParent,
# ##                 family=zero_inflated_poisson(link='log'))
# 
# ## owls_brm6 <- brm(owls_form,
# ##                  data=owls,
# ##                  prior = priors,
# ##                  sample_prior = 'yes',
# ##                  iter=5000,
# ##                  warmup=2500,
# ##                  thin=10,
# ##                  chains=3,
# ##                  refresh=0,
# ##                  backend = "cmdstanr",
# ##                  cores=3)


## ----fitModel2a, results='markdown', eval=TRUE, mhidden=TRUE, paged.print=FALSE, tidy.opts = list(width.cutoff = 80)----
owls_form <- bf(NCalls ~ FoodTreatment*SexParent +
                    offset(log(BroodSize)) + (1|Nest),
                family=poisson(link='log'))
options(width=150)
owls_form |> get_prior(data = owls)
options(width=80)


## ----fitModel2h_1a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
owls_form <- bf(NCalls ~ +
                    offset(log(BroodSize)) + (1|Nest),
                family=poisson(link='log'))
options(width=150)
owls_form |> get_prior(data = owls)
options(width=80)
owls |>
    summarise(Mean = log(median(NCalls/BroodSize)),
              SD = log(sd(NCalls/BroodSize)),
              MAD = log(mad(NCalls/BroodSize)))
priors <- prior(normal(0.2, 0.6), class = 'Intercept') +
    prior(student_t(3, 0, 0.6), class = 'sd')
owls_brm1a <- brm(owls_form,
                 data = owls,
                 prior = priors,
                 sample_prior = 'only',
                 iter = 5000,
                 warmup =2500,
                 chains = 3,
                 cores = 3,
                 thin = 10,
                 refresh = 0,
                 seed = 123,
                 control =  list(adapt_delta = 0.99),
                 backend = "cmdstanr"
                 )
owls_brm1a |>
  emmeans(~1, type = "response") |>
  as.data.frame() |>
  ggplot(aes(x = 1, y = rate)) +
  geom_pointrange(aes(ymin = lower.HPD, ymax = upper.HPD)) +
  geom_point(data = owls, aes(y = NCalls/BroodSize),
    position = position_jitter(width = 0.2))


## ----fitModel2h_1b, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
owls_brm1b <- update(owls_brm1a,
                       sample_prior = 'yes',
                       iter =  5000,
                       warmup =  2500,
                       control = list(adapt_delta = 0.99, max_treedepth = 20),
                       backend =  "cmdstanr",
                       refresh = 0)
owls_brm1b |>
  emmeans(~1, type = "response") |>
  as.data.frame() |>
  ggplot(aes(x = 1, y = rate)) +
  geom_pointrange(aes(ymin = lower.HPD, ymax = upper.HPD)) +
  geom_point(data = owls, aes(y = NCalls/BroodSize),
    position = position_jitter(width = 0.2))

## ----posterior2h_1c, results='markdown', eval=TRUE, fig.width = 7, fig.height = 5----
owls_brm1b |> SUYR_prior_and_posterior()


## ----fitModel2h, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
owls |>
    group_by(FoodTreatment, SexParent) |>
    summarise(Mean = log(median(NCalls/BroodSize)),
              SD = log(sd(NCalls/BroodSize)),
              MAD = log(mad(NCalls/BroodSize)))
standist::visualize("normal(0.4, 0.5)", xlim=c(0,20))
standist::visualize("student_t(3, 0, 2.5)",
                    "cauchy(0,1)",
                    xlim=c(-10,25))


## ----modelValidation1c_a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_mcmc()


## ----modelValidation1c_b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pars <- owls_brm1b |> get_variables()
pars <- pars |> str_extract('^b.Intercept|^b_FoodTreatment.*|^b_SexParent.*|[sS]igma|^sd.*') |>
    na.omit()
pars
owls_brm1b |> mcmc_plot(type='trace', variable = pars)
##OR
owls_brm1b |> mcmc_plot(type='trace',
                        variable = '^b.Intercept|^b_FoodTreatment.*|^b_SexParent.*|[sS]igma|^sd.*',
                        regex = TRUE)



## ----modelValidation1c_c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm1b |> mcmc_plot(type='acf_bar', variable = pars)
##OR
owls_brm1b |> mcmc_plot(type='acf_bar',
                        variable = '^b.Intercept|^b_FoodTreatment.*|^b_SexParent.*|[sS]igma|^sd.*',
                        regex = TRUE)


## ----modelValidation1c_d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm1b |> mcmc_plot(type='rhat_hist')


## ----modelValidation1c_e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm1b |> mcmc_plot(type='neff_hist')


## ----modelValidation1c_f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm1b |> mcmc_plot(type='combo', pars = pars)
owls_brm1b |> mcmc_plot(type='violin', pars = pars)


## ----modelValidation1c_g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm1b |> get_variables()
pars <- owls_brm1b |> get_variables()
pars <- str_extract(pars, '^b_.*|^sigma$|^sd.*') |> na.omit()

owls_brm1b$fit |>
    stan_trace(pars = pars)


## ----modelValidation1c_h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm1b$fit |>
    stan_ac(pars = pars)


## ----modelValidation1c_i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm1b$fit |> stan_rhat()


## ----modelValidation1c_j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm1b$fit |> stan_ess()


## ----modelValidation1c_k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm1b$fit |>
    stan_dens(separate_chains = TRUE, pars = pars)


## ----modelValidation1c_l, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=7----
## owls_ggs <- owls_brm3 |> ggs(burnin = FALSE, inc_warmup = FALSE)
## owls_ggs |> ggs_traceplot()


## ----modelValidation1c_m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=7----
## ggs_autocorrelation(owls_ggs)


## ----modelValidation1c_n, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## ggs_Rhat(owls_ggs)


## ----modelValidation1c_o, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## ggs_effective(owls_ggs)


## ----modelValidation1c_p, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## ggs_crosscorrelation(owls_ggs)


## ----modelValidation1c_q, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## ggs_grb(owls_ggs)


## ----modelValidation1c_r, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation1c_s, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm1b |> pp_check(type = 'dens_overlay', ndraws = 100)


## ----modelValidation1c_t, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## owls_brm1b |> pp_check(type = 'error_scatter_avg')


## ----modelValidation1c_u, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm1b |> pp_check(group = 'Nest', type = 'intervals')


## ----modelValidation1c_v, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(owls_brm2)


## ----modelValidation1c_w, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
owls_resids <- make_brms_dharma_res(owls_brm1b, integerResponse = TRUE)
wrap_elements(~testUniformity(owls_resids)) +
               wrap_elements(~plotResiduals(owls_resids, form = factor(rep(1, nrow(owls))))) +
               wrap_elements(~plotResiduals(owls_resids, quantreg = TRUE)) +
               wrap_elements(~testDispersion(owls_resids))


## ----validation1c_x, results='markdown', eval=TRUE, error=TRUE,mhidden=TRUE, fig.width=7, fig.height=5, cache=FALSE, message=FALSE, warning=FALSE----
try({
owls_resids |> testZeroInflation()
})


## ----validation1c_y, results='markdown', eval=TRUE, error=TRUE,mhidden=TRUE, fig.width=7, fig.height=5, cache=TRUE, message=FALSE, warning=FALSE----
try({
owls_resids |> testTemporalAutocorrelation(time=owls$ArrivalTime)
## owls_resid1 <- owls_resids |> recalculateResiduals(group=interaction(owls$ArrivalTime,  owls$Nest),  aggregateBy = mean)
## owls_resid1 <- owls_resids |> recalculateResiduals(group=interaction(owls$ArrivalTime,  owls$Nest),  aggregateBy = sum)
## resids1 <- owls_resids |> recalculateResiduals(group = interaction(owls$ArrivalTime), aggregateBy = sum)
resids1 <- owls_resids |> recalculateResiduals(group = interaction(owls$ArrivalTime), aggregateBy = mean)
resids1 |> testTemporalAutocorrelation(time=unique(owls$ArrivalTime))

autocor_check(owls, owls_brm1b,
              variable =  "ArrivalTime",
              grouping = "Nest",
              n.sim =  250)
})


## ----fitModel2h_2a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
owls_form <- bf(NCalls ~ +
                  offset(log(BroodSize)) + (1|Nest),
  zi = ~ 1,
  family=zero_inflated_poisson(link='log'))
options(width=150)
owls_form |> get_prior(data = owls)
options(width=80)
owls |>
    summarise(Mean = log(median(NCalls/BroodSize)),
              SD = log(sd(NCalls/BroodSize)),
              MAD = log(mad(NCalls/BroodSize)))
priors <- prior(normal(0.2, 0.6), class = 'Intercept') +
  prior(student_t(3, 0, 0.6), class = 'sd') +
  prior(normal(0, 1), class = 'Intercept', dpar = 'zi')
owls_brm2a <- brm(owls_form,
                 data = owls,
                 prior = priors,
                 sample_prior = 'only',
                 iter = 5000,
                 warmup =2500,
                 chains = 3,
                 cores = 3,
                 thin = 10,
                 refresh = 0,
                 seed = 123,
                 control =  list(adapt_delta = 0.99),
                 backend = "cmdstanr"
                 )
owls_brm2a |>
  emmeans(~1, type = "response") |>
  as.data.frame() |>
  ggplot(aes(x = 1, y = rate)) +
  geom_pointrange(aes(ymin = lower.HPD, ymax = upper.HPD)) +
  geom_point(data = owls, aes(y = NCalls/BroodSize),
    position = position_jitter(width = 0.2))


## ----fitModel2h_2b, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
owls_brm2b <- update(owls_brm2a,
                       sample_prior = 'yes',
                       iter =  5000,
                       warmup =  2500,
                       control = list(adapt_delta = 0.99, max_treedepth = 20),
                       backend =  "cmdstanr",
                       refresh = 0)
owls_brm2b |>
  emmeans(~1, type = "response") |>
  as.data.frame() |>
  ggplot(aes(x = 1, y = rate)) +
  geom_pointrange(aes(ymin = lower.HPD, ymax = upper.HPD)) +
  geom_point(data = owls, aes(y = NCalls/BroodSize),
    position = position_jitter(width = 0.2))

## ----posterior2h_2c, results='markdown', eval=TRUE, fig.width = 7, fig.height = 5----
owls_brm2b |> SUYR_prior_and_posterior()


## ----modelValidation2c_g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm2b |> get_variables()
pars <- owls_brm2b |> get_variables()
pars <- str_extract(pars, '^b_.*|^sigma$|^sd.*') |> na.omit()

owls_brm2b$fit |>
    stan_trace(pars = pars)

## ----modelValidation2c_h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm2b$fit |>
    stan_ac(pars = pars)

## ----modelValidation2c_i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm2b$fit |> stan_rhat()

## ----modelValidation2c_j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm2b$fit |> stan_ess()

## ----modelValidation2c_k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm2b$fit |>
    stan_dens(separate_chains = TRUE, pars = pars)


## ----modelValidation2c_s, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm2b |> pp_check(type = 'dens_overlay', ndraws = 100)


## ----modelValidation2c_t, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## owls_brm2b |> pp_check(type = 'error_scatter_avg')


## ----modelValidation2c_u, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm2b |> pp_check(group = 'Nest', type = 'intervals')


## ----modelValidation2c_w, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
owls_resids <- make_brms_dharma_res(owls_brm2b, integerResponse = TRUE)
wrap_elements(~testUniformity(owls_resids)) +
               wrap_elements(~plotResiduals(owls_resids, form = factor(rep(1, nrow(owls))))) +
               wrap_elements(~plotResiduals(owls_resids, quantreg = TRUE)) +
               wrap_elements(~testDispersion(owls_resids))


## ----fitModel2h_3a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
owls_form <- bf(NCalls ~ +
                  offset(log(BroodSize)) + (1|Nest),
  zi = ~ 1,
  family=zero_inflated_negbinomial(link='log'))
options(width=150)
owls_form |> get_prior(data = owls)
options(width=80)
owls |>
    summarise(Mean = log(median(NCalls/BroodSize)),
              SD = log(sd(NCalls/BroodSize)),
              MAD = log(mad(NCalls/BroodSize)))
priors <- prior(normal(0.2, 0.6), class = 'Intercept') +
  prior(student_t(3, 0, 0.6), class = 'sd') +
  prior(normal(0, 1), class = 'Intercept', dpar = 'zi') +
  prior(gamma(0.01, 0.01), class = "shape")
owls_brm3a <- brm(owls_form,
                 data = owls,
                 prior = priors,
                 sample_prior = 'only',
                 iter = 5000,
                 warmup =2500,
                 chains = 3,
                 cores = 3,
                 thin = 10,
                 refresh = 0,
                 seed = 123,
                 control =  list(adapt_delta = 0.99),
                 backend = "cmdstanr"
                 )


## ----fitModel2h_3b, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
owls_brm3b <- update(owls_brm3a,
                       sample_prior = 'yes',
                       iter =  5000,
                       warmup =  2500,
                       control = list(adapt_delta = 0.99, max_treedepth = 20),
                       backend =  "cmdstanr",
                       refresh = 0)
owls_brm3b |>
  emmeans(~1, type = "response") |>
  as.data.frame() |>
  ggplot(aes(x = 1, y = prob)) +
  geom_pointrange(aes(ymin = lower.HPD, ymax = upper.HPD)) +
  geom_point(data = owls, aes(y = NCalls/BroodSize),
    position = position_jitter(width = 0.2))

## ----posterior2h_3c, results='markdown', eval=TRUE, fig.width = 7, fig.height = 5----
owls_brm3b |> SUYR_prior_and_posterior()


## ----modelValidation3c_g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm3b |> get_variables()
pars <- owls_brm3b |> get_variables()
pars <- str_extract(pars, '^b_.*|^sigma$|^sd.*') |> na.omit()

owls_brm3b$fit |>
    stan_trace(pars = pars)

## ----modelValidation3c_h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm3b$fit |>
    stan_ac(pars = pars)

## ----modelValidation3c_i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm3b$fit |> stan_rhat()

## ----modelValidation3c_j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm3b$fit |> stan_ess()

## ----modelValidation3c_k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm3b$fit |>
    stan_dens(separate_chains = TRUE, pars = pars)


## ----modelValidation3c_s, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm3b |> pp_check(type = 'dens_overlay', ndraws = 100)


## ----modelValidation3c_t, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## owls_brm3b |> pp_check(type = 'error_scatter_avg')


## ----modelValidation3c_u, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm3b |> pp_check(group = 'Nest', type = 'intervals')


## ----modelValidation3c_w, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
owls_resids <- make_brms_dharma_res(owls_brm3b, integerResponse = TRUE)
wrap_elements(~testUniformity(owls_resids)) +
               wrap_elements(~plotResiduals(owls_resids, form = factor(rep(1, nrow(owls))))) +
               wrap_elements(~plotResiduals(owls_resids, quantreg = TRUE)) +
               wrap_elements(~testDispersion(owls_resids))


## ----fitModel2h_4a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
owls_form <- bf(NCalls ~ FoodTreatment*SexParent +
                  offset(log(BroodSize)) + (1|Nest),
  zi = ~ 1,
  family=zero_inflated_negbinomial(link='log'))
options(width=150)
owls_form |> get_prior(data = owls)
options(width=80)
owls |>
  group_by(FoodTreatment, SexParent) |>
    summarise(Mean = log(median(NCalls/BroodSize)),
              SD = log(sd(NCalls/BroodSize)),
              MAD = log(mad(NCalls/BroodSize)))
priors <- prior(normal(0.4, 0.7), class = 'Intercept') +
  prior(normal(0, 1), class = "b") +
  prior(student_t(3, 0, 0.7), class = 'sd') +
  prior(normal(0, 1), class = 'Intercept', dpar = 'zi') +
  prior(gamma(0.01, 0.01), class = "shape")
owls_brm4a <- brm(owls_form,
                 data = owls,
                 prior = priors,
                 sample_prior = 'only',
                 iter = 5000,
                 warmup =2500,
                 chains = 3,
                 cores = 3,
                 thin = 10,
                 refresh = 0,
                 seed = 123,
                 control =  list(adapt_delta = 0.99),
                 backend = "cmdstanr"
                 )


## ----fitModel2h_4b, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
owls_brm4b <- update(owls_brm4a,
                       sample_prior = 'yes',
                       iter =  5000,
                       warmup =  2500,
                       control = list(adapt_delta = 0.99, max_treedepth = 20),
                       backend =  "cmdstanr",
                       refresh = 0)
owls_brm4b |>
  conditional_effects("FoodTreatment:SexParent") |>
  plot(points = TRUE)

## ----posterior2h_4c, results='markdown', eval=TRUE, fig.width = 7, fig.height = 5----
owls_brm4b |> SUYR_prior_and_posterior()


## ----modelValidation4c_g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm4b |> get_variables()
pars <- owls_brm4b |> get_variables()
pars <- str_extract(pars, '^b_.*|^sigma$|^sd.*') |> na.omit()

owls_brm4b$fit |>
    stan_trace(pars = pars)

## ----modelValidation4c_h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm4b$fit |>
    stan_ac(pars = pars)

## ----modelValidation4c_i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm4b$fit |> stan_rhat()

## ----modelValidation4c_j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm4b$fit |> stan_ess()

## ----modelValidation4c_k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm4b$fit |>
    stan_dens(separate_chains = TRUE, pars = pars)


## ----modelValidation4c_s, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm4b |> pp_check(type = 'dens_overlay', ndraws = 100)


## ----modelValidation4c_t, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## owls_brm4b |> pp_check(type = 'error_scatter_avg')


## ----modelValidation4c_u, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm4b |> pp_check(group = 'Nest', type = 'intervals')


## ----modelValidation4c_w, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
owls_resids <- make_brms_dharma_res(owls_brm4b, integerResponse = TRUE)
wrap_elements(~testUniformity(owls_resids)) +
               wrap_elements(~plotResiduals(owls_resids, form = factor(rep(1, nrow(owls))))) +
               wrap_elements(~plotResiduals(owls_resids, quantreg = TRUE)) +
               wrap_elements(~testDispersion(owls_resids))


## ----fitModel2h_5a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
owls_form <- bf(NCalls ~ FoodTreatment*SexParent +
                  offset(log(BroodSize)) +
                  (FoodTreatment*SexParent|Nest),
  zi = ~ FoodTreatment*SexParent,
  family=zero_inflated_negbinomial(link='log'))
options(width=150)
owls_form |> get_prior(data = owls)
options(width=80)
owls |>
  group_by(FoodTreatment, SexParent) |>
    summarise(Mean = log(median(NCalls/BroodSize)),
              SD = log(sd(NCalls/BroodSize)),
              MAD = log(mad(NCalls/BroodSize)))
priors <- prior(normal(0.4, 0.7), class = 'Intercept') +
  prior(normal(0, 1), class = "b") +
  prior(student_t(3, 0, 0.7), class = 'sd') +
  prior(normal(0, 1), class = 'Intercept', dpar = 'zi') +
  prior(normal(0, 1), class = "b", dpar = 'zi') +
  prior(gamma(0.01, 0.01), class = "shape")
owls_brm5a <- brm(owls_form,
                 data = owls,
                 prior = priors,
                 sample_prior = 'only',
                 iter = 5000,
                 warmup =2500,
                 chains = 3,
                 cores = 3,
                 thin = 10,
                 refresh = 0,
                 seed = 123,
                 control =  list(adapt_delta = 0.99),
                 backend = "cmdstanr"
                 )


## ----fitModel2h_5b, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
owls_brm5b <- update(owls_brm5a,
                       sample_prior = 'yes',
                       iter =  5000,
                       warmup =  2500,
                       control = list(adapt_delta = 0.99, max_treedepth = 20),
                       backend =  "cmdstanr",
                       refresh = 0)
owls_brm5b |>
  conditional_effects("FoodTreatment:SexParent") |>
  plot(points = TRUE)

## ----posterior2h_5c, results='markdown', eval=TRUE, fig.width = 7, fig.height = 5----
## owls_brm5b |> SUYR_prior_and_posterior()


## ----modelValidation5c_g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm5b |> get_variables()
pars <- owls_brm5b |> get_variables()
pars <- str_extract(pars, '^b_.*|^sigma$|^sd.*') |> na.omit()

owls_brm5b$fit |>
    stan_trace(pars = pars)

## ----modelValidation5c_h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm5b$fit |>
    stan_ac(pars = pars)

## ----modelValidation5c_i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm5b$fit |> stan_rhat()

## ----modelValidation5c_j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm5b$fit |> stan_ess()

## ----modelValidation5c_k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm5b$fit |>
    stan_dens(separate_chains = TRUE, pars = pars)


## ----modelValidation5c_s, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm5b |> pp_check(type = 'dens_overlay', ndraws = 100)


## ----modelValidation5c_t, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## owls_brm5b |> pp_check(type = 'error_scatter_avg')


## ----modelValidation5c_u, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
owls_brm5b |> pp_check(group = 'Nest', type = 'intervals')


## ----modelValidation5c_w, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
owls_resids <- make_brms_dharma_res(owls_brm5b, integerResponse = TRUE)
wrap_elements(~testUniformity(owls_resids)) +
               wrap_elements(~plotResiduals(owls_resids, form = factor(rep(1, nrow(owls))))) +
               wrap_elements(~plotResiduals(owls_resids, quantreg = TRUE)) +
               wrap_elements(~testDispersion(owls_resids))


## ----modelValidation6a_1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
(loo_3 <- loo::loo(owls_brm3b))
(loo_4 <- loo::loo(owls_brm4b))
(loo_5 <- loo::loo(owls_brm5b))
loo::loo_compare(loo_3, loo_4, loo_5)



## ----partialPlot2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
owls_brm5b |>
    conditional_effects("FoodTreatment:SexParent") |>
    plot(points = TRUE, jitter_width = 0.25)
## If we want to present the conditional effects on the scale of calls
## per chick, then we need to condition on a broodsize of 1
owls_brm5b |>
    conditional_effects("FoodTreatment:SexParent", conditions = list(BroodSize = 1))


## ----partialPlot2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
owls_brm5b |>
    ggpredict(~FoodTreatment*SexParent) |>
    plot(show_data = TRUE, jitter = c(0.25, 0))


## ----partialPlot2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
off <- owls |> summarize(Mean=mean(BroodSize))
owls_brm5b |>
    ggemmeans(~FoodTreatment*SexParent, offset = log(off$Mean)) |>
    plot()


## ----summariseModel2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
owls_brm5b |> summary()


## ----summariseModel2a1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=FALSE----
owls_sum <- summary(owls_brm5b)


## ----summariseModel2bm, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
owls_brm5b |> as_draws_df()
owls_brm5b |>
  as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(
    median,
    HDInterval::hdi,
    Pl = ~mean(.x < 1),
    Pg = ~mean(.x > 1),
    rhat,
    length,
    ess_bulk,
    ess_tail
  )
0.205/(1+0.205)
(0.205*2.26)/(1+(0.205*2.26))
(0.205*0.449)/(1+(0.205*0.449))


## ----summariseModel2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
owls_brm5b$fit |>
    tidyMCMC(estimate.method = 'median',
             conf.int = TRUE,  conf.method = 'HPDinterval',
             rhat = TRUE, ess = TRUE)

## ----summariseModel2b1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
owls_tidy <- tidyMCMC(owls_brm5b$fit, estimate.method='median',
                         conf.int=TRUE,  conf.method='HPDinterval',
                         rhat=TRUE, ess=TRUE)


## ----summariseModel2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
owls_brm5b |> get_variables()
owls_draw <- owls_brm5b |>
    gather_draws(`b.Intercept.*|b_FoodTreatment.*|b_SexParent.*`,  regex=TRUE)
owls_draw


## ----summariseModel2c1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
owls_draw |> median_hdci()
## On a fractional scale
owls_draw |>
    mutate(.value = exp(.value)) |>
    median_hdci()


## ----summariseModel2c3, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
owls_gather <- owls_brm5b |>
    gather_draws(`b_Intercept.*|b_FoodTreatment.*|b_SexParent.*`,  regex=TRUE) |>
    mutate(.value = exp(.value)) |>
    median_hdci()


## ----summariseModel2c4, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
owls_brm5b |>
    gather_draws(`b_Intercept.*|b_FoodTreatment.*|b_SexParent.*`, regex=TRUE) |>
    ## mutate(.value = exp(.value)) |>
    ggplot() +
    geom_vline(xintercept=0, linetype='dashed') +
    stat_slab(aes(x = .value, y = .variable,
                  fill = stat(ggdist::cut_cdf_qi(cdf,
                                                 .width = c(0.5, 0.8, 0.95),
                                                 labels = scales::percent_format())
                              )), color='black') +
    scale_fill_brewer('Interval', direction = -1, na.translate = FALSE)

owls_brm5b |>
    gather_draws(`.Intercept.*|b_FoodTreatment.*|b_SexParent.*`, regex=TRUE) |>
    ggplot() +
    geom_vline(xintercept = 0, linetype='dashed') +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    theme_classic()


## ----summariseModel2j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
owls_brm5b$fit |> plot(type='intervals')


## ----summariseModel2ka, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
owls_brm5b |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(.variable != 'b_Intercept') |>
    ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    facet_wrap(~.variable, scales='free') +
    theme(axis.text.y = element_blank())

owls_brm5b |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(.variable != 'b_Intercept') |>
    ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    geom_vline(xintercept = 0, linetype = 'dashed')


## ----summariseModel2c7, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
owls_brm5b |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(str_detect(.variable, 'b_.*Intercept', negate = TRUE)) |>
    ggplot() +
    geom_density_ridges(aes(x=.value, y = .variable), alpha=0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed')
##Or in colour
owls_brm5b |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(str_detect(.variable, 'b_.*Intercept', negate = TRUE)) |>
    ggplot() +
    geom_density_ridges_gradient(aes(x=(.value),
                                     y = .variable,
                                     fill = stat(x)),
                                 alpha=0.4, colour = 'white',
                                 quantile_lines = TRUE,
                                 quantiles = c(0.025, 0.975)) +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    scale_x_continuous() +
    scale_fill_viridis_c(option = "C")

## Fractional scale
owls_brm5b |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(str_detect(.variable, 'b_.*Intercept', negate = TRUE)) |>
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
owls_brm5b |> tidy_draws()


## ----summariseModel2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
owls_brm5b |> spread_draws(`.*Intercept.*|b_FoodTreatment.*|b_SexParent.*`,  regex=TRUE)


## ----summariseModel2f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
owls_brm5b |> posterior_samples() |> as_tibble()


## ----summariseModel2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
owls_brm5b |>
    bayes_R2(re.form = NA, summary=FALSE) |>
    median_hdci()
owls_brm5b |>
    bayes_R2(re.form = ~(1|Nest), summary=FALSE) |>
    median_hdci()
owls_brm5b |>
    bayes_R2(re.form = ~(FoodTreatment*SexParent|Nest), summary=FALSE) |>
    median_hdci()


## -----------------------------------------------------------------------------
#| label: modelsummary
#| results: markup
#| eval: false
#| echo: true
#| cache: false
# owls_brm5b |> modelsummary(
#   statistic = c("conf.low", "conf.high"),
#   shape = term ~ statistic,
#   exponentiate = TRUE
# )


## -----------------------------------------------------------------------------
#| label: modelsummary_plot
#| results: markup
#| eval: false
#| echo: true
#| cache: false
# owls_brm5b |> modelplot(exponentiate = TRUE)


## ----postHoc1a, results='markdown', eval=TRUE,mhidden=TRUE--------------------
## ## The following should work, but there is a bug and therefore it does not
## ## (although it has been reported - so may get fixed at some point).
## ## The offset seems to get handled incorrectly
## newdata <- owls_brm5b |>
##     emmeans(~FoodTreatment|SexParent, offset=0, type='response') |>
##     as.data.frame()
newdata <- owls_brm5b |>
    emmeans(~FoodTreatment|SexParent, type='response') |>
    as.data.frame()
head(newdata)
## As an alternative, we can do the following...
newdata <- owls_brm5b |>
    emmeans(~FoodTreatment|SexParent,
          at = list(BroodSize = 1), type='response') |>
    as.data.frame()
head(newdata)
ggplot(newdata) +
    geom_pointrange(aes(y=prob,  x=FoodTreatment,  color=SexParent,
                        ymin=lower.HPD,  ymax=upper.HPD),
                    position=position_dodge(width=0.2)) +
    theme_classic()


## ----fitmodelrstanarm, results='markdown', eval=FALSE,mhidden=TRUE, echo=FALSE----
# starling.rstan |> get_variables()
# plot(starling.rstan,  'mcmc_trace', regex_pars = '^.intercept|^situation|^month|[ss]igma')
# plot(starling.rstan,  'mcmc_acf_bar', regex_pars = '^.intercept|^situation|^month|[ss]igma')
# plot(starling.rstan,  'mcmc_rhat_hist', regex_pars = '^.intercept|^situation|^month|[ss]igma')
# plot(starling.rstan,  'mcmc_neff_hist', regex_pars = '^.intercept|^situation|^month|[ss]igma')
# 
# 
# preds <- posterior_predict(starling.rstan,  nsamples=250,  summary=false)
# starling.resids <- createdharma(simulatedresponse = t(preds),
#                             observedresponse = starling$mass,
#                             fittedpredictedresponse = apply(preds, 2, median))
# plot(starling.resids)
# 
# 
# starling.rstan1 = stan_glmer(mass ~ month*situation+(month|bird),data=starling,
#                             iter=5000, warmup=2000, thin=5, chains=3, refresh=0)
# starling.rstan1 = stan_glmer(mass ~ month*situation+(month|bird),data=starling,
#                              iter=5000, warmup=2000, thin=5, chains=3, refresh=0,
#                              adapt_delta = 0.99)
# #pairs(starling.rstan1,  pars=c('(intercept)', 'monthnov'))
# starling.rstan1 |> get_variables()
# pairs(starling.rstan1,  regex_pars=c('situation', 'sigma'))
# prior_summary(starling.rstan1)
# 
# plot(starling.rstan1,  'mcmc_trace', regex_pars = '^.intercept|^situation|^month|[ss]igma')
# plot(starling.rstan1,  'mcmc_acf_bar', regex_pars = '^.intercept|^situation|^month|[ss]igma')
# plot(starling.rstan1,  'mcmc_rhat_hist', regex_pars = '^.intercept|^situation|^month|[ss]igma')
# plot(starling.rstan1,  'mcmc_neff_hist', regex_pars = '^.intercept|^situation|^month|[ss]igma')
# 
# starling.rstan1 = stan_glmer(mass ~ month*situation+(month|bird),data=starling,
#                              iter=10000, warmup=5000, thin=15, chains=3, refresh=0,
#                              adapt_delta = 0.99)
# 
# plot(starling.rstan1,  'mcmc_trace', regex_pars = '^.intercept|^situation|^month|[ss]igma')
# plot(starling.rstan1,  'mcmc_acf_bar', regex_pars = '^.intercept|^situation|^month|[ss]igma')
# plot(starling.rstan1,  'mcmc_rhat_hist', regex_pars = '^.intercept|^situation|^month|[ss]igma')
# plot(starling.rstan1,  'mcmc_neff_hist', regex_pars = '^.intercept|^situation|^month|[ss]igma')
# preds <- posterior_predict(starling.rstan1,  nsamples=250,  summary=FALSE)
# starling.resids <- createDHARMa(simulatedResponse = t(preds),
#                             observedResponse = starling$MASS,
#                             fittedPredictedResponse = apply(preds, 2, median))
# plot(starling.resids)
# 
# (l.1 <- loo(starling.rstan))
# (l.2 <- loo(starling.rstan1))
# loo_compare(l.1, l.2)
# 
# as.matrix(starling.rstan) |> colnames()
# posterior_vs_prior(starling.rstan1, color_by='vs', group_by=TRUE, regex_pars=c('^MONTH','^SITUATION','^[sS]igma'),
#                    facet_args=list(scales='free_y'))
# 
# 
# g=ggpredict(starling.rstan1) |> plot()
# do.call('grid.arrange',  g)
# ggemmeans(starling.rstan1, ~SITUATION|MONTH) |> plot()
# 
# summary(starling.rstan1)
# 
# nms <- starling.rstan1 |> get_variables()
# nms
# wch <- grep('^.Intercept|^MONTH|^SITUATION|[sS]igma', nms)
# tidyMCMC(starling.rstan1$stanfit,conf.int=TRUE, conf.method='HPDinterval',
#          rhat=TRUE, ess=TRUE, pars=nms[wch])
# 
# emmeans(starling.rstan1, pairwise~MONTH|SITUATION)
# starling.em = emmeans(starling.rstan1, ~MONTH|SITUATION) |>
#     gather_emmeans_draws() |> spread(key=MONTH, value=.value) |>
#     mutate(Eff=Jan-Nov,
#            PEff=100*(Jan-Nov)/Nov)
# starling.em |> head()
# 
# starling.em |> ungroup() |>
#   dplyr::select(SITUATION,Eff,PEff) |>
#   group_by(SITUATION) |>
#   median_hdi()
# 
# starling.em |> ungroup() |>
#   dplyr::select(SITUATION,Eff,PEff) |>
#   group_by(SITUATION) |>
#   summarize(Prob=mean(PEff>10))
# 
# bayes_R2(starling.rstan1, re.form=NA) |> median_hdi()
# bayes_R2(starling.rstan1, re.form=~(1|BIRD)) |> median_hdi()
# bayes_R2(starling.rstan1, re.form=~(MONTH|BIRD)) |> median_hdi()
# 
# newdata = emmeans(starling.rstan1, ~MONTH|SITUATION) |> as.data.frame()
# head(newdata)
# ggplot(newdata, aes(y=emmean, x=SITUATION)) +
#     geom_pointrange(aes(ymin=lower.HPD, ymax=upper.HPD, fill=MONTH),
#                     position=position_dodge(width=0.3), shape=21)


## ----fitModel, results='markdown', eval=FALSE, mhidden=TRUE, echo=FALSE-------
# owls_rstanP <- stan_glmer(NCalls ~ FoodTreatment+SexParent +
#                            offset(log(BroodSize)) + (1|Nest),
#                          dat=owls,  family=poisson(link='log'), refresh=0,
#                          iter=5000,  warmup=2000,  thin=10,  chains=3, cores=3)
# 
# 
# owls_rstanP <- stan_glmer(NCalls ~ FoodTreatment+scale(ArrivalTime) +
#                            offset(log(BroodSize)) + (1|Nest),
#                          dat=owls,  family=poisson(link='log'), refresh=0,
#                          iter=5000,  warmup=2000,  thin=10,  chains=3, cores=3)
# 
# owls_rstanP %>% get_variables()
# plot(owls_rstanP,  'mcmc_trace', regex_pars='^.Intercept|Food|Arrival|[sS]igma')
# plot(owls_rstanP,  'mcmc_acf_bar', regex_pars='^.Intercept|Food|Arrival|[sS]igma')
# plot(owls_rstanP,  'mcmc_rhat_hist', regex_pars='^.Intercept|Food|Arrival|[sS]igma')
# plot(owls_rstanP,  'mcmc_neff_hist', regex_pars='^.Intercept|Food|Arrival|[sS]igma')
# 
# 
# preds <- posterior_predict(owls_rstanP,  nsamples=250,  summary=FALSE)
# owls_resids <- createDHARMa(simulatedResponse = t(preds),
#                             observedResponse = owls$NCalls,
#                             fittedPredictedResponse = apply(preds, 2, median))
# plot(owls_resids)
# 
# 
# owls_rstanNB <- stan_glmer(NCalls ~ FoodTreatment+SexParent +
#                            offset(log(BroodSize)) + (1|Nest),
#                          dat=owls,  family=neg_binomial_2(link='log'), refresh=0,
#                          iter=5000,  warmup=2000,  thin=10,  chains=3, cores=3)
# 
# owls_rstanNB <- stan_glmer(NCalls ~ FoodTreatment+scale(ArrivalTime) +
#                            offset(log(BroodSize)) + (1|Nest),
#                          dat=owls,  family=neg_binomial_2(link='log'), refresh=0,
#                          iter=5000,  warmup=2000,  thin=10,  chains=3, cores=3)
# 
# owls_rstanNB %>% get_variables()
# plot(owls_rstanNB,  'mcmc_trace', regex_pars='^.Intercept|Food|Arrival|[sS]igma')
# plot(owls_rstanNB,  'mcmc_acf_bar', regex_pars='^.Intercept|Food|Arrival|[sS]igma')
# plot(owls_rstanNB,  'mcmc_rhat_hist', regex_pars='^.Intercept|Food|Arrival|[sS]igma')
# plot(owls_rstanNB,  'mcmc_neff_hist', regex_pars='^.Intercept|Food|Arrival|[sS]igma')
# 
# 
# preds <- posterior_predict(owls_rstanNB,  nsamples=250,  summary=FALSE)
# owls_resids <- createDHARMa(simulatedResponse = t(preds),
#                             observedResponse = owls$NCalls,
#                             fittedPredictedResponse = apply(preds, 2, median))
# plot(owls_resids)
# testZeroInflation(owls_resids)
# 
# testTemporalAutocorrelation(owls_resids,  time=owls$ArrivalTime)
# owls_resids1 <- recalculateResiduals(owls_resids,  group=interaction(owls$ArrivalTime,  owls$Nest),  aggregateBy = mean)
# testTemporalAutocorrelation(owls_resids1,  time=unique(owls$ArrivalTime))
# 
# ##Cant use zero inflation with glmer
# 
# owls_rstan1 <- stan_glmer(NCalls ~ FoodTreatment+scale(ArrivalTime) +
#                            offset(log(BroodSize)) + (FoodTreatment*scale(ArrivalTime)|Nest),
#                          dat=owls,  family=neg_binomial_2(link='log'), refresh=0,
#                          iter=5000,  warmup=2000,  thin=10,  chains=3, cores=3)
# 
# owls_rstan1 %>% get_variables()
# plot(owls_rstan1,  'mcmc_trace', regex_pars='^.Intercept|^Food|^scale.Arrival|[sS]igma')
# plot(owls_rstan1,  'mcmc_acf_bar', regex_pars='^.Intercept|^Food|^scale.Arrival|[sS]igma')
# plot(owls_rstan1,  'mcmc_rhat_hist', regex_pars='^.Intercept|^Food|^scale.Arrival|[sS]igma')
# plot(owls_rstan1,  'mcmc_neff_hist', regex_pars='^.Intercept|^Food|^scale.Arrival|[sS]igma')
# 
# preds <- posterior_predict(owls_rstan1,  nsamples=250,  summary=FALSE)
# owls_resids <- createDHARMa(simulatedResponse = t(preds),
#                             observedResponse = owls$NCalls,
#                             fittedPredictedResponse = apply(preds, 2, median))
# plot(owls_resids)
# 


## ----fitModel.brms, results='markdown', eval=FALSE, mhidden=TRUE, echo=FALSE----
# owls_form <- bf(NCalls ~ FoodTreatment*SexParent +
#                     offset(log(BroodSize)) + (1|Nest),
#                 family=poisson(link='log'))
# owls_brm1 <- brm(owls_form,
#                   data=owls,
#                   sample_prior = 'yes',
#                   iter=5000,  warmup=2000,
#                   thin=10,  chains=3, cores=3,
#                   refresh=0,
#                   )
# prior_summary(owls_brm1)
# 
# owls %>%
#     group_by(FoodTreatment, SexParent) %>%
#     summarise(log(median(NCalls)),
#               log(mad(NCalls)))
# priors <- prior(normal(1.8, 5), class='Intercept') +
#     prior(normal(0, 2), class='b') +
#     prior(gamma(2,1), class='sd')
# owls_form <- bf(NCalls ~ FoodTreatment*SexParent +
#                     offset(log(BroodSize)) + (1|Nest),
#                 family=poisson(link='log'))
# owls_brm1 <- brm(owls_form,
#                  data=owls,
#                  prior = priors,
#                  sample_prior = 'yes',
#                  iter=5000,  warmup=2000,
#                  thin=10,  chains=3, cores=3,
#                  refresh=0,
#                  )
# 
# owls_form <- bf(NCalls ~ FoodTreatment*SexParent +
#                     offset(log(BroodSize)) + (FoodTreatment*SexParent|Nest),
#                 family=poisson(link='log'))
# owls_brm2 <- brm(owls_form,
#                  data=owls,
#                  prior = priors,
#                  sample_prior = 'yes',
#                  iter=5000,  warmup=2500,
#                  thin=10,  chains=3, cores=3,
#                  refresh=0,
#                  )
# 
# (l.1 <- loo(owls_brm1))
# (l.2 <- loo(owls_brm2))
# loo_compare(l.1, l.2)
# 
# 
# owls_brm2 %>% get_variables()
# owls_brm2 %>% hypothesis('FoodTreatmentSatiated=0') %>% plot()
# owls_brm2 %>% hypothesis('b_FoodTreatmentSatiated=0', class='') %>% plot()
# 
# pars <- owls_brm2 %>% get_variables()
# wch <- grepl('^b_.*|^sd_.*', pars, perl=TRUE)
# wch <- grepl('^[bsd]{1,2}_.*', pars, perl=TRUE)
# 
# 
# g <- vector('list', length=sum(wch))
# names(g) <- pars[wch]
# 
# for (i in pars[wch]) {
#     print(i)
# }
# for (i in pars[wch]) {
#     print(i)
#     if (i == 'b_Intercept') next
#     p <- owls_brm2 %>% hypothesis(paste0(i,'=0'), class='') %>% plot()
#     g[[i]] <- p[[1]]
# }
# patchwork::wrap_plots(g)
# 
# p <- owls_brm2 %>% hypothesis('FoodTreatmentSatiated=0') %>% plot()
# 
# g <- vector('list', length=sum(wch)-1)
# names(g) <- pars[wch][-1]
# for (i in pars[wch]) {
#     print(i)
#     if (i == 'b_Intercept') next
#     p <- owls_brm2 %>% hypothesis(paste0(i,'=0'), class='') %>% plot()
#     g[[i]] <- p[[1]]
# }
# patchwork::wrap_plots(g)
# 
# stan_trace(owls_brm2$fit, pars = pars[wch])
# stan_ac(owls_brm2$fit, pars = pars[wch])
# stan_rhat(owls_brm2$fit, pars = pars[wch])
# stan_rhat(owls_brm2$fit)
# stan_ess(owls_brm2$fit)
# 
# ##mcmc_plot(owls_brmsP,  type='trace')
# ##mcmc_plot(owls_brmsP,  type='acf_bar')
# ##mcmc_plot(owls_brmsP,  type='rhat_hist')
# ##mcmc_plot(owls_brmsP,  type='neff_hist')
# ##mcmc_plot(owls_brmsP,  type='trace', regex_pars='^b.Intercept|Food|Arrival|sd')
# ##mcmc_plot(owls_brmsP,  type='acf_bar', regex_pars='^b.Intercept|Food|Arrival|sd')
# ##mcmc_plot(owls_brmsP,  type='rhat_hist', regex_pars='^b.Intercept|Food|Arrival|sd')
# ##mcmc_plot(owls_brmsP,  type='neff_hist', regex_pars='^b.Intercept|Food|Arrival|sd')
# 
# preds <- posterior_predict(owls_brm2,  nsamples=250,  summary=FALSE)
# owls_resids <- createDHARMa(simulatedResponse = t(preds),
#                             observedResponse = owls$NCalls,
#                             fittedPredictedResponse = apply(preds, 2, median),
#                             integerResponse = TRUE)
# plot(owls_resids)
# testZeroInflation(owls_resids)
# 
# 
# 
# ##owls_form <- bf(NCalls ~ FoodTreatment*SexParent +
# ##                    offset(log(BroodSize)) + (FoodTreatment*SexParent|Nest),
# ##                family=negbinomial(link='log'))
# ##owls_brmsNB <- brm(owls_form,  data=owls, refresh=0,
# ##                  iter=5000,  warmup=2000,  thin=10,  chains=3, cores=3)
# 
# ##mcmc_plot(owls_brmsNB,  type='trace')
# ##mcmc_plot(owls_brmsNB,  type='acf_bar')
# ##mcmc_plot(owls_brmsNB,  type='rhat_hist')
# ##mcmc_plot(owls_brmsNB,  type='neff_hist')
# ##preds <- posterior_predict(owls_brmsNB,  nsamples=250,  summary=FALSE)
# ##owls_resids <- createDHARMa(simulatedResponse = t(preds),
# ##                            observedResponse = owls$NCalls,
# ##                            fittedPredictedResponse = apply(preds, 2, median),
# ##                            integerResponse = TRUE)
# ##plot(owls_resids)
# ##testZeroInflation(owls_resids)
# 
# priors <- prior(normal(1.8, 5), class='Intercept') +
#     prior(normal(0, 5), class='b') +
#     prior(gamma(2,1), class='sd') +
#     prior(logistic(0,1), class='Intercept', dpar='zi') +
#     prior(normal(0,1), class='b', dpar='zi')
# owls_form <- bf(NCalls ~ FoodTreatment*SexParent +
#                     offset(log(BroodSize)) +
#                     (FoodTreatment*SexParent|Nest),
#                 zi ~ FoodTreatment+SexParent,
#                 family=zero_inflated_poisson(link='log'))
#                 ##family=zero_inflated_negbinomial(link='log'))
# 
# owls_brm3 <- brm(owls_form,
#                  data=owls,
#                  prior = priors,
#                  sample_prior = 'yes',
#                  iter=5000,  warmup=2500,
#                  thin=10,  chains=3, cores=3,
#                  refresh=0)
# prior_summary(owls_brm3)
# 
# pars <- owls_brm3 %>% get_variables()
# pars
# wch <- grepl('^b_.*|^sd_.*', pars, perl=TRUE)
# 
# g <- vector('list', length=sum(wch))
# names(g) <- pars[wch]
# for (i in pars[wch]) {
#     print(i)
#     p <- owls_brm3 %>% hypothesis(paste0(i,'=0'), class='') %>% plot()
#     g[[i]] <- p[[1]]
# }
# patchwork::wrap_plots(g)
# 
# g <- vector('list', length=sum(wch)-1)
# names(g) <- pars[wch][-1]
# for (i in pars[wch]) {
#     print(i)
#     if (i == 'b_Intercept') next
#     p <- owls_brm3 %>% hypothesis(paste0(i,'=0'), class='') %>% plot()
#     g[[i]] <- p[[1]]
# }
# patchwork::wrap_plots(g[[2:5]])
# 
# stan_trace(owls_brm3$fit, pars = pars[wch])
# stan_ac(owls_brm3$fit, pars = pars[wch])
# stan_rhat(owls_brm3$fit, pars = pars[wch])
# stan_rhat(owls_brm3$fit)
# stan_ess(owls_brm3$fit)
# 
# 
# 
# preds <- posterior_predict(owls_brm3,  nsamples=250,  summary=FALSE)
# owls_resids <- createDHARMa(simulatedResponse = t(preds),
#                             observedResponse = owls$NCalls,
#                             fittedPredictedResponse = apply(preds, 2, median),
#                             integerResponse = TRUE)
# plot(owls_resids)
# testZeroInflation(owls_resids)
# testDispersion(owls_resids)
# 
# priors <- prior(normal(1.8, 5), class='Intercept') +
#     prior(normal(0, 5), class='b') +
#     prior(gamma(2,1), class='sd') +
#     prior(logistic(0,1), class='Intercept', dpar='zi') +
#     prior(normal(0,1), class='b', dpar='zi') +
#     prior(gamma(0.01, 0.01), class='shape')
# owls_form <- bf(NCalls ~ FoodTreatment*SexParent +
#                     offset(log(BroodSize)) +
#                     (FoodTreatment*SexParent|Nest),
#                 zi ~ FoodTreatment+SexParent,
#                 family=zero_inflated_negbinomial(link='log'))
# 
# owls_brm4 <- brm(owls_form,
#                  data=owls,
#                  prior = priors,
#                  sample_prior = 'yes',
#                  iter=5000,  warmup=2500,
#                  thin=10,  chains=3, cores=3,
#                  refresh=0)
# 
# pars <- owls_brm4 %>% get_variables()
# pars
# wch <- grepl('^b_.*|^sd_.*', pars, perl=TRUE)
# 
# g <- vector('list', length=sum(wch))
# names(g) <- pars[wch]
# for (i in pars[wch]) {
#     print(i)
#     p <- owls_brm4 %>% hypothesis(paste0(i,'=0'), class='') %>% plot()
#     g[[i]] <- p[[1]]
# }
# patchwork::wrap_plots(g)
# 
# wch <- grepl('^b_.*|^sd_.*|.*shape.*', pars, perl=TRUE)
# stan_trace(owls_brm4$fit, pars = pars[wch])
# stan_ac(owls_brm4$fit, pars = pars[wch])
# stan_rhat(owls_brm4$fit, pars = pars[wch])
# stan_rhat(owls_brm4$fit)
# stan_ess(owls_brm4$fit)
# 
# preds <- posterior_predict(owls_brm4,  nsamples=250,  summary=FALSE)
# owls_resids <- createDHARMa(simulatedResponse = t(preds),
#                             observedResponse = owls$NCalls,
#                             fittedPredictedResponse = apply(preds, 2, median),
#                             integerResponse = TRUE)
# plot(owls_resids)
# testZeroInflation(owls_resids)
# testDispersion(owls_resids)
# 
# testTemporalAutocorrelation(owls_resids,  time=owls$ArrivalTime)
# owls_resids1 <- recalculateResiduals(owls_resids,  group=interaction(owls$ArrivalTime,  owls$Nest),  aggregateBy = mean)
# testTemporalAutocorrelation(owls_resids1,  time=unique(owls$ArrivalTime))
# 
# g <- owls_brm4 %>%
#     conditional_effects() %>%
#     plot(ask=FALSE, plot=FALSE, points=TRUE)
# patchwork::wrap_plots(g)
# 
# summary(owls_brm4)
# 
# tidyMCMC(owls_brmsZINB2$fit,  estimate.method='median',
#          conf.int=TRUE,  conf.method='HPDinterval',
#          rhat=TRUE,  ess=TRUE)
# 
# owls_brm4 %>% bayes_R2(re.form = NA, summary=FALSE) %>% median_hdci()
# owls_brm4 %>% bayes_R2(re.form = ~(1|Nest), summary=FALSE) %>% median_hdci()
# owls_brm4 %>% bayes_R2(re.form = ~(FoodTreatment*SexParent|Nest), summary=FALSE) %>% median_hdci()
# 
# owls_brm4 %>%
#     emmeans(~FoodTreatment, type='response')
# owls_brm4 %>%
#     emmeans(~FoodTreatment, offset=0, type='response')
# 
# newdata <- owls_brm4 %>%
#     emmeans(~FoodTreatment, offset=0, type='response') %>%
#     as.data.frame
# head(newdata)
# ggplot(newdata) +
#     geom_pointrange(aes(y=prob,  x=FoodTreatment,  color=SexParent,
#                         ymin=lower.HPD,  ymax=upper.HPD),
#                     position=position_dodge(width=0.2))
# ggplot(newdata, aes(y=prob, x=as.numeric(SexParent)) +
#     geom_line(aes(color=FoodTreatment)) +
#     geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD, fill=FoodTreatment), alpha=0.2)
# 

