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
source('helperFunctions.R')


## ----readData, results='markdown', eval=TRUE----------------------------------
mckeon <- read_csv("../data/mckeon.csv", trim_ws = TRUE)


## -----------------------------------------------------------------------------
#| label: examinData
glimpse(mckeon)


## -----------------------------------------------------------------------------
## Explore the first 6 rows of the data
head(mckeon)


## -----------------------------------------------------------------------------
str(mckeon)


## -----------------------------------------------------------------------------
mckeon |> datawizard::data_codebook()


## ----processData, results='markdown', eval=TRUE, mhidden=TRUE-----------------
mckeon <- mckeon |>
  mutate(BLOCK = factor(BLOCK),
         SYMBIONT = factor(SYMBIONT, levels = c('none', 'crabs', 'shrimp', 'both')))


## ----eda1a, results='markdown', eval=TRUE, mhidden=TRUE-----------------------
ggplot(mckeon, aes(y=PREDATION, x=SYMBIONT)) +
    geom_point(position=position_jitter(width=0.2, height=0))+
    facet_wrap(~BLOCK)


## ----fitModel1a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
mckeon_rstanarm <- stan_glmer(PREDATION ~ SYMBIONT + (1|BLOCK),
                           data = mckeon,
                           family = binomial(link = 'logit'),
                           iter = 5000,
                           warmup = 2000,
                           chains = 3,
                           thin = 5,
                           refresh = 0,
                           cores = 3)


## ----fitModel1b, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
mckeon_rstanarm %>% prior_summary()


## ----fitModel1d, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
model.matrix(~SYMBIONT, data=mckeon) %>%
    apply(2,sd) %>%
    (function(x) 2.5/x)


## ----fitModel1f, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
mckeon_rstanarm1 <- update(mckeon_rstanarm,  prior_PD=TRUE)


## ----fitModel1g, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
mckeon_rstanarm1 %>% ggpredict() %>% plot(show_data=TRUE, jitter=c(0.5,0))


## ----fitModel1g2, results='markdown', eval=TRUE, mhidden=TRUE, cache = FALSE----
mckeon_rstanarm1 %>% ggemmeans(~SYMBIONT) %>% plot(show_data=TRUE, jitter=c(0.5,0))


## ----fitModel1g3, results='markdown', eval=TRUE, mhidden=TRUE, cache = FALSE----
mckeon_rstanarm1 %>% emmeans(~SYMBIONT, type = 'link') %>%
    as.data.frame() %>%
    ggplot(aes(y = emmean, x = SYMBIONT)) +
    geom_hline(yintercept = c(5,-5), linetype = 'dashed') +
    geom_pointrange(aes(ymin = lower.HPD, ymax = upper.HPD)) +
    geom_point(data = mckeon, aes(y = PREDATION),
               position = position_jitter(width=0.2, height = 0),
               alpha=0.4, color = 'red')


## ----fitModel1h, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
mckeon_rstanarm2 <- stan_glmer(PREDATION ~ SYMBIONT + (1|BLOCK),
                                data = mckeon,
                                family = binomial(link='logit'),
                                prior_intercept = normal(0, 2.5, autoscale = FALSE),
                                prior = normal(0, 6, autoscale = FALSE),
                                prior_covariance = decov(1, 1, 1, 1),
                                prior_PD = TRUE,
                                iter = 5000,
                                warmup = 1000,
                                chains = 3,
                                thin = 5,
                                refresh = 0
                                )


## ----fitModel1i, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
mckeon_rstanarm2 %>%
    ggpredict(~SYMBIONT) %>%
    plot(show_data = TRUE, jitter = c(0.5, 0))
mckeon_rstanarm2 %>% ggemmeans(~SYMBIONT) %>% plot(show_data=TRUE, jitter=c(0.5,0))


## ----fitModel1j, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE, dependson='fitModel1h'----
mckeon_rstanarm3 <- update(mckeon_rstanarm2,  prior_PD=FALSE)


## ----modelFit1k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
posterior_vs_prior(mckeon_rstanarm3, color_by='vs', group_by=TRUE,
                   facet_args=list(scales='free_y'))


## ----modelFit1l, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggemmeans(mckeon_rstanarm3,  ~SYMBIONT) %>% plot(show_data=TRUE)
ggpredict(mckeon_rstanarm3,  ~SYMBIONT) %>% plot(show_data=TRUE)


## ----fitModel2a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE, paged.print=FALSE, tidy.opts = list(width.cutoff = 80)----
mckeon_form <- bf(PREDATION | trials(1) ~ SYMBIONT + (1|BLOCK),
                  family=binomial(link='logit'))
#OR
mckeon_form <- bf(PREDATION ~ SYMBIONT + (1|BLOCK),
                  family=bernoulli(link='logit'))
options(width=150)
mckeon_form |> get_prior(data = mckeon)
options(width=80)


## ----fitModel2h, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
mckeon_form <- bf(PREDATION | trials(1) ~ SYMBIONT + (1|BLOCK),
                  family=binomial(link='logit'))
#OR
mckeon_form <- bf(PREDATION ~ SYMBIONT + (1|BLOCK),
                  family=bernoulli(link='logit'))
options(width=150)
mckeon_form |> get_prior(data = mckeon)
options(width=80)

mckeon |> group_by(SYMBIONT) |>
  summarise(Mean = logit(mean(PREDATION)),
    Median = logit(median(PREDATION)),
    sd = logit(abs(sd(PREDATION))))
2.5/model.matrix(~SYMBIONT, data=mckeon) |> apply(2,sd)


mckeon |>
    group_by(SYMBIONT) |>
    summarise(
        mean_response = mean(PREDATION),
        mad_response = mad(PREDATION),
        sd_response = sd(PREDATION),
    ) |>
    mutate(
        mean_logit = logit(mean_response),
        # Delta method approximation
      sd_logit = sd_response / (mean_response * (1 - mean_response))
    )

standist::visualize("student_t(3, 0, 3.5)",
                    "gamma(2,0.5)",
                    "cauchy(0,2)",
                    xlim=c(-10,25))
standist::visualize("student_t(3, 0, 3.5)",
                    "gamma(2,0.5)",
                    "cauchy(0,2)",
                    xlim=c(-10,25))


## ----fitModel2h1, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE-----
mckeon_form <- bf(PREDATION | trials(1) ~ SYMBIONT + (1|BLOCK),
                  family=binomial(link='logit'))
#OR
mckeon_form <- bf(PREDATION ~ SYMBIONT + (1|BLOCK),
                  family=bernoulli(link='logit'))
get_prior(mckeon_form, data = mckeon)
priors <-
    prior(normal(0, 2.5), class = 'Intercept') +
    prior(normal(0, 3), class = 'b') +
    prior(student_t(3, 0, 1.5), class = 'sd')
priors <-
    prior(normal(0, 2), class = 'Intercept') +
    prior(normal(0, 3), class = 'b') +
    prior(student_t(3, 0, 1.5), class = 'sd')

mckeon_brm2 <- brm(mckeon_form,
                  data = mckeon,
                  prior = priors,
                  sample_prior = 'only',
                  iter = 5000,
                  warmup = 1000,
                  chains = 3, cores = 3,
                  thin = 5,
                  control = list(adapt_delta = 0.99, max_treedepth =  20),
                  refresh = 0,
                  backend = "cmdstan"
                  )


## ----fitModel10, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
mckeon_brm2 |>
    conditional_effects() |>
    plot(points = TRUE)

mckeon_brm2 |>
    ggpredict(~SYMBIONT) |>
    plot(show_data = TRUE)

mckeon_brm2 |>
    ggemmeans(~SYMBIONT) |>
    plot(show_data = TRUE)


## ----fitModel1i2, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE----
mckeon_brm2 |> emmeans(~SYMBIONT, type = 'link') |>
    as.data.frame() |>
    ggplot(aes(y = emmean, x = SYMBIONT)) +
    geom_hline(yintercept = c(5,-5), linetype = 'dashed') +
    geom_pointrange(aes(ymin = lower.HPD, ymax = upper.HPD)) +
    geom_point(data = mckeon, aes(y = PREDATION),
               position = position_jitter(width=0.2, height = 0),
               alpha=0.4, color = 'red')


## ----fitModel1j1, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE, dependson='fitModel1h'----
mckeon_brm3 <- update(mckeon_brm2,
                      sample_prior = 'yes',
                      cores = 3,
                      refresh = 0)
save(mckeon_brm3, file = '../ws/testing/mckeon_brm3')


## ----partialPlot2h1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mckeon_brm3 |>
    conditional_effects() |>
    plot(points = TRUE)

mckeon_brm3 |>
    ggpredict(~SYMBIONT) |>
    plot(show_data = TRUE)

mckeon_brm3 |>
    ggemmeans(~SYMBIONT) |>
    plot(show_data = TRUE)

mckeon_brm3 |> emmeans(~SYMBIONT, type = 'link') |>
    as.data.frame() |>
    ggplot(aes(y = emmean, x = SYMBIONT)) +
    geom_hline(yintercept = c(5,-5), linetype = 'dashed') +
    geom_pointrange(aes(ymin = lower.HPD, ymax = upper.HPD)) +
    geom_point(data = mckeon, aes(y = PREDATION),
               position = position_jitter(width=0.2, height = 0),
               alpha=0.4, color = 'red')


## ----posterior2h2, results='markdown', eval=TRUE------------------------------
mckeon_brm3 |> get_variables()
mckeon_brm3 |> hypothesis('SYMBIONTcrabs=0') |> plot()
mckeon_brm3 |> hypothesis('SYMBIONTshrimp=0') |> plot()


## ----posterior2h2a, results='markdown', eval=TRUE, fig.width = 7, fig.height = 5----
mckeon_brm3 |> SUYR_prior_and_posterior()


## ----fitModel2h3, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE-----
mckeon_form <- bf(PREDATION | trials(1) ~ SYMBIONT + (SYMBIONT|BLOCK),
                  family=binomial(link='logit'))
#OR
mckeon_form <- bf(PREDATION ~ SYMBIONT + (SYMBIONT|BLOCK),
                  family=bernoulli(link='logit'))
get_prior(mckeon_form, mckeon)
## As there are not many observations, the following might be too ambitious without
## stronger priors
priors <-
    prior(normal(0, 2), class = 'Intercept') +
    prior(normal(0, 3), class = 'b') +
    ## prior(student_t(3,0,1.5), class = 'sd', coef =  "Intercept", group =  "BLOCK") +
    prior(student_t(3, 0, 1.5), class = 'sd') +
    prior(lkj_corr_cholesky(1), class = 'cor')

mckeon_brm4 <- brm(mckeon_form,
                  data = mckeon,
                  prior = priors,
                  sample_prior = 'yes',
                  iter = 5000,
                  warmup = 1000,
                  chains = 3, cores = 3,
                  thin = 5,
                  refresh = 0,
                  control = list(adapt_delta = 0.99, max_treedepth = 20),
                  backend = 'cmdstan'
                  )

save(mckeon_brm4, file = '../ws/testing/mckeon_brm4')


## ----posterior2k, results='markdown', eval=TRUE-------------------------------
mckeon_brm4 |> get_variables()
mckeon_brm4 |> hypothesis('SYMBIONTcrabs=0') |> plot()
mckeon_brm4 |> hypothesis('SYMBIONTshrimp=0') |> plot()


## ----posterior2k1, results='markdown', eval=TRUE, fig.width = 7, fig.height = 5----
mckeon_brm4 |> SUYR_prior_and_posterior()
mckeon_brm3 |> SUYR_prior_and_posterior()
mckeon_brm4 |>
    conditional_effects() |>
    plot(points = TRUE)
mckeon_brm3 |>
    conditional_effects() |>
    plot(points = TRUE)
mckeon_brm4 |> summary()


## ----posterior2k2, results='markdown', eval=TRUE, fig.width=10, fig.height=4----
mckeon_brm4 |>
  posterior_samples() |>
  dplyr::select(-`lp__`) |>
  pivot_longer(everything(), names_to = 'key') |>
  filter(!str_detect(key, '^r')) |>
  mutate(Type = ifelse(str_detect(key, 'prior'), 'Prior', 'Posterior'),
         ## Class = ifelse(str_detect(key, 'Intercept'),  'Intercept',
         ##         ifelse(str_detect(key, 'b'),  'b', 'sigma')),
         Class = case_when(
               str_detect(key, '(^b|^prior).*Intercept$') ~ 'Intercept',
               str_detect(key, 'b_SYMBIONT.*|prior_b_SYMBIONT.*') &
               !str_detect(key, '.*:.*') ~ 'SYMBIONT',
               str_detect(key, 'sd') ~ 'sd',
               str_detect(key, '^cor|prior_cor') ~ 'cor',
             str_detect(key, 'sigma') ~ 'sigma'
             ),
         Par = str_replace(key, 'b_', '')) |>
  ggplot(aes(x = Type,  y = value, color = Par)) +
  stat_pointinterval(position = position_dodge())+
  facet_wrap(~Class,  scales = 'free')


## ----fitModel2h3a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
(l.1 <- mckeon_brm3 |> loo())
(l.2 <- mckeon_brm4 |> loo())
loo_compare(l.1, l.2)


## ----posterior2k20, results='markdown', eval=TRUE, fig.width=10, fig.height=4----
mckeon_brm4 %>%
  posterior_samples %>%
  dplyr::select(-`lp__`) %>%
  pivot_longer(everything(), names_to = 'key') %>%
  filter(!str_detect(key, '^r')) %>%
  mutate(Type = ifelse(str_detect(key, 'prior'), 'Prior', 'Posterior'),
         Class = case_when(
             str_detect(key, '(^b|^prior).*Intercept$') ~ 'Intercept',
             str_detect(key, 'b_SYMBIONT.*|prior_b') ~ 'TREATMENT',
             str_detect(key, 'sd') ~ 'sd',
             str_detect(key, '^cor|prior_cor') ~ 'cor',
             str_detect(key, 'sigma') ~ 'sigma'),
         Par = str_replace(key, 'b_', '')) %>%
  ggplot(aes(x = Type,  y = value, color = Par)) +
  stat_pointinterval(position = position_dodge(), show.legend = FALSE)+
  facet_wrap(~Class,  scales = 'free')



## ----modelValidation2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_mcmc()


## ----modelValidation2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
pars <- mckeon_brm4 |> get_variables()
pars <- pars |> str_extract('^b.Intercept|^b_SYMBIONT.*|[sS]igma|^sd.*') |>
    na.omit()
pars
mckeon_brm4 |> mcmc_plot(type='trace', variables = pars)
#OR
mckeon_brm4 |> mcmc_plot(type='trace',
                          variable = '^b.Intercept|^b_SYMBIONT.*|[sS]igma|^sd.*',
                          regex = TRUE)



## ----modelValidation2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
mckeon_brm4 |> mcmc_plot(type='acf_bar', variable = pars)
##OR
mckeon_brm4 |> mcmc_plot(type='acf_bar',
                          variable = '^b.Intercept|^b_SYMBIONT.*|[sS]igma|^sd.*',
                          regex = TRUE)



## ----modelValidation2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
mckeon_brm4 |> mcmc_plot(type='rhat_hist')


## ----modelValidation2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
mckeon_brm4 |> mcmc_plot(type='neff_hist')


## ----modelValidation2f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
mckeon_brm4 |> mcmc_plot(type='combo', variable = pars)
mckeon_brm4 |> mcmc_plot(type='violin', variable = pars)


## ----modelValidation2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
mckeon_brm4 |> get_variables()
pars <- mckeon_brm4 |> get_variables()
pars <- str_extract(pars, '^b_.*|^sigma$|^sd.*') |> na.omit()

mckeon_brm4$fit |>
    stan_trace(pars = pars)
## mckeon_brm3$fit |> stan_trace()


## ----modelValidation2h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
mckeon_brm4$fit |>
    stan_ac(pars = pars)
## mckeon_brm3$fit |> stan_ac()


## ----modelValidation2i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
mckeon_brm4$fit |> stan_rhat()


## ----modelValidation2j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
mckeon_brm4$fit |> stan_ess()


## ----modelValidation2k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
mckeon_brm4$fit |>
    stan_dens(separate_chains = TRUE, pars = pars)


## ----modelValidation2l, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=7----
## mckeon_ggs <- mckeon_brm3 %>% ggs(burnin = FALSE, inc_warmup = FALSE)
## mckeon_ggs %>% ggs_traceplot()


## ----modelValidation2m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=7----
## ggs_autocorrelation(mckeon_ggs)


## ----modelValidation2n, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## ggs_Rhat(mckeon_ggs)


## ----modelValidation2o, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## ggs_effective(mckeon_ggs)


## ----modelValidation2p, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## ggs_crosscorrelation(mckeon_ggs)


## ----modelValidation2q, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## ggs_grb(mckeon_ggs)


## ----modelValidation5a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation5b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
mckeon_brm4 |> pp_check(type = 'dens_overlay', nsamples = 250)


## ----modelValidation5c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## mckeon_brm4 |> pp_check(type = 'error_scatter_avg')


## ----modelValidation5e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
mckeon_brm4 |> pp_check(group = 'BLOCK', type = 'intervals')


## ----modelValidation5g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(mckeon_brm2)


## ----modelValidation6aa, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
mckeon_resids <- make_brms_dharma_res(mckeon_brm4, integerResponse = FALSE)
wrap_elements(~testUniformity(mckeon_resids)) +
               wrap_elements(~plotResiduals(mckeon_resids, form = factor(rep(1, nrow(mckeon))))) +
               wrap_elements(~plotResiduals(mckeon_resids, quantreg = FALSE)) +
               wrap_elements(~testDispersion(mckeon_resids))



## ----partialPlot2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mckeon_brm4 |>
    conditional_effects() |>
    plot(points = TRUE)


## ----partialPlot2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mckeon_brm4 |>
    ggpredict() |>
    plot(show_data = TRUE, jitter=c(0.5,0))


## ----partialPlot2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mckeon_brm4 |>
    ggemmeans(~SYMBIONT) |>
    plot(show_data = TRUE, jitter=c(0.5,0))


## ----partialPlot2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
## Partial residuals in binomial models are too confusing for the average viewer
## as they will yeild values that are not exactly 0 or 1 and this seems wrong.
## Partial.obs <- mckeon_brm3$data %>%
##     mutate(Pred = predict(mckeon_brm3, re.form=NA)[,'Estimate'],
##            Resid = resid(mckeon_brm3)[,'Estimate'],
##            Obs = Pred + Resid)

mckeon_brm4 |>
  epred_draws(newdata = mckeon, re_formula = NA) |>
  median_hdci() |>
  ggplot(aes(x = SYMBIONT, y = .epred)) +
  geom_pointrange(aes(ymin = .lower, ymax = .upper)) +
  geom_line() +
  geom_point(data = mckeon,  aes(y = PREDATION,  x = SYMBIONT),
             alpha=0.2,
             position = position_jitter(width= 0.2, height = 0))


## ----summariseModel2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mckeon_brm4 |> summary()


## ----summariseModel2a1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=FALSE----
mckeon_sum <- summary(mckeon_brm4)


## ----summariseModel2bm, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
mckeon_brm4 |> as_draws_df()
mckeon_brm4 |>
  as_draws_df() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(
    median,
    HDInterval::hdi,
    Pl = ~mean(.x < 1),
    Pg = ~mean(.x > 1),
    rhat,
    ess_bulk,
    ess_tail
  ) |>
  knitr::kable()

mckeon_brm4 |>
    as_draws_df() |>
    exp() |>
    dplyr::select(matches("^b_.*|^sd_.*")) |>
  summarise_draws(
    median,
    HDInterval::hdi,
    Pl = ~mean(.x < 1),
    Pg = ~mean(.x > 1),
    rhat,
    ess_bulk,
    ess_tail
  ) |>
  knitr::kable()


## ----summariseModel2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mckeon_brm4$fit |>
    tidyMCMC(estimate.method = 'median',
             conf.int = TRUE,  conf.method = 'HPDinterval',
             rhat = TRUE, ess = TRUE)

## ----summariseModel2b1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
mckeon_tidy <- tidyMCMC(mckeon_brm4$fit, estimate.method='median',
                         conf.int=TRUE,  conf.method='HPDinterval',
                         rhat=TRUE, ess=TRUE)


## ----summariseModel2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mckeon_brm4 |> get_variables()
mckeon_draw <- mckeon_brm4 |>
    gather_draws(`b.Intercept.*|b_SYMBIONT.*`,  regex=TRUE)
mckeon_draw


## ----summariseModel2c1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mckeon_draw |> median_hdci()
## On a odd ratio scale
mckeon_draw |>
    mutate(.value = exp(.value)) |>
    median_hdci()


## ----summariseModel2c3, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
mckeon_gather <- mckeon_brm4 |>
    gather_draws(`b_Intercept.*|b_SYMBIONT.*`, regex = TRUE) |>
    median_hdci()


## ----summariseModel2c4, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
mckeon_brm4 |>
    gather_draws(`b_Intercept.*|b_SYMBIONT.*`, regex=TRUE) |>
    ggplot() +
    geom_vline(xintercept=0, linetype='dashed') +
    stat_slab(aes(x = .value, y = .variable,
                  fill = stat(ggdist::cut_cdf_qi(cdf,
                                                 .width = c(0.5, 0.8, 0.95),
                                                 labels = scales::percent_format())
                              )), color='black') +
    scale_fill_brewer('Interval', direction = -1, na.translate = FALSE)

mckeon_brm4 |>
    gather_draws(`.Intercept.*|b_SYMBIONT.*`, regex=TRUE) |>
    ggplot() +
    geom_vline(xintercept = 0, linetype='dashed') +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    theme_classic()


## ----summariseModel2j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mckeon_brm4$fit |> plot(type='intervals')


## ----summariseModel2ka, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
mckeon_brm4 |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(.variable != 'b_Intercept') |>
    ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    facet_wrap(~.variable, scales='free')

mckeon_brm4 |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(.variable != 'b_Intercept') |>
    ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    geom_vline(xintercept = 0, linetype = 'dashed')


## ----summariseModel2c7, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
mckeon_brm4 |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(.variable != 'b_Intercept') |>
    ggplot() +
    geom_density_ridges(aes(x=.value, y = .variable), alpha=0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed')
##Or in colour
mckeon_brm4 |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(.variable != 'b_Intercept') |>
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
mckeon_brm4 |>
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


## ----summariseModel2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mckeon_brm4 |> tidy_draws()


## ----summariseModel2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mckeon_brm4 |> spread_draws(`.*Intercept.*|b_SYMBIONT.*`,  regex=TRUE)


## ----summariseModel2f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mckeon_brm4 |> posterior_samples() |> as_tibble()


## ----summariseModel2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
mckeon_brm4 |>
    bayes_R2(re.form = NA, summary = FALSE) |>
    median_hdci()
mckeon_brm4 |>
    bayes_R2(re.form = ~ (1 | BLOCK), summary = FALSE) |>
    median_hdci()
## if we had random intercept/slope
mckeon_brm4 |>
    bayes_R2(re.form = ~ (SYMBIONT | BLOCK), summary = FALSE) |>
    median_hdci()


## -----------------------------------------------------------------------------
#| label: modelsummary
#| results: markup
#| eval: true
#| echo: true
#| cache: false
mckeon_brm4 |> modelsummary(
  statistic = c("conf.low", "conf.high"),
  shape = term ~ statistic,
  exponentiate = FALSE
)


## -----------------------------------------------------------------------------
#| label: modelsummary_plot
#| results: markup
#| eval: true
#| echo: true
#| cache: false
mckeon_brm4 |> modelplot(exponentiate = FALSE)


## ----posthoc1a, results='markdown', eval=TRUE, mhidden=TRUE-------------------
mckeon_brm4 |>
    emmeans(~SYMBIONT, type='response') |>
    pairs()

# via tidy_draws
mckeon_em <-
    mckeon_brm4 |>
    emmeans(~SYMBIONT, type = "link") |>
    pairs() |>
    tidy_draws() |>
    mutate(across(starts_with("contrast"), exp)) |>
    summarise_draws(median,
                    HDInterval::hdi,
                    Pl = ~mean(.x < 1),
                    Pg = ~mean(.x > 1),
                    ROPE = ~ mean(.x < 1.1 & .x > 0.9)
                    )
mckeon_em
# or via gather_emmeans_draws()
mckeon_em <-
    mckeon_brm4 |>
    emmeans(~SYMBIONT, type = "link") |>
    pairs() |>
    gather_emmeans_draws() |>
  mutate(.value = exp(.value)) |>
  summarise(median_hdci(.value),
    Pl = mean(.value < 1),
    Pg = mean(.value > 1)
    )
mckeon_em

## On a probability scale
mckeon_em <-
    mckeon_brm4 |>
    emmeans(~SYMBIONT, type = "link") |>
    regrid() |>
    pairs() |>
    tidy_draws() |>
    summarise_draws(median,
                    HDInterval::hdi,
                    Pl = ~mean(.x < 0),
                    Pg = ~mean(.x > 0)
                    )
mckeon_em

# or via gather_emmeans_draws()
## mckeon_em <- mckeon_brm4 |>
##     emmeans(~SYMBIONT, type='link') |>
##     pairs() |>
##     gather_emmeans_draws() |>
##     mutate(Eff=exp(.value),
##            PEff=100*(Eff-1))#,
##              # Prob = plogis(.value))
## mckeon_em |> head()
## mckeon_em |>
##     group_by(contrast) |>
##     dplyr::select(contrast, Eff) |>
##     median_hdi()
## mckeon_em |>
##   group_by(contrast) |>
##   summarize(Prob=sum(Eff>1)/n())

## On a probability scale
mckeon_em <- mckeon_brm4 |>
    emmeans(~SYMBIONT, type='link') |>
    regrid() |>
    pairs() |>
    gather_emmeans_draws() |>
    mutate(Eff=.value)#,
mckeon_em |> head()
mckeon_em |>
    group_by(contrast) |>
    dplyr::select(contrast, Eff) |>
    median_hdi()

## Cell means
mckeon_em = emmeans(mckeon_brm4, ~SYMBIONT, type='link') |>
      gather_emmeans_draws()
mckeon_em |> mutate(P=plogis(.value)) |> median_hdci(P)


## ----posteriors1a, results='markdown', eval=TRUE, mhidden=TRUE----------------
cmat <- cbind(
    "Crab vs shrimp" =c(0,1,-1,0),
    "Both vs One"=c(0,-1/2,-1/2,1),
    "Any vs None"=c(-1, 1/3, 1/3, 1/3)
)

mckeon_em <- mckeon_brm4 |>
    emmeans(~SYMBIONT, type = 'response') |>
    contrast(method = list(cmat))

mckeon_em <-
    mckeon_brm4 |>
    emmeans(~SYMBIONT, type = "link") |>
    contrast(method = list(cmat)) |>
    tidy_draws() |>
    mutate(across(starts_with("contrast"), exp)) |>
    summarise_draws(median,
                    HDInterval::hdi,
                    Pl = ~mean(.x < 1),
                    Pg = ~mean(.x > 1)
                    )
mckeon_em
## or via gather_emmeans_draws
mckeon_em <-
  mckeon_brm4 |>
    emmeans(~SYMBIONT, type='link') |>
    contrast(method=list(cmat)) |>
    gather_emmeans_draws() |>
  mutate(Fit=exp(.value)) |>
  summarise(median_hdci(Fit),
    Pl =  mean(Fit < 0),
    Pg =  mean(Fit > 0)
    )
mckeon_em

newdata <- emmeans(mckeon_brm4, ~SYMBIONT, type = "response") |> as.data.frame()
head(newdata)
ggplot(newdata, aes(y=response, x=SYMBIONT)) +
    geom_pointrange(aes(ymin=lower.HPD, ymax=upper.HPD))


## ----R21a, results='markdown', eval=TRUE, mhidden=TRUE------------------------
mckeon_brm4 |> bayes_R2(re.form=NA)
mckeon_brm4 |> bayes_R2(re.form=NA, summary=FALSE) |> median_hdci()
mckeon_brm4 |> bayes_R2(re.form=~(1|BLOCK), summary=FALSE) |> median_hdci()
## for random intercept/slope model
mckeon_brm4 |> bayes_R2(re.form=~(SYMBIONT|BLOCK), summary=FALSE) |> median_hdci()

