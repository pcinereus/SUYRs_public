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

library(tidyverse)  #for data wrangling etc
library(rstanarm)   #for fitting models in STAN
library(cmdstanr)   #for cmdstan
library(brms)       #for fitting models in STAN
library(standist)   #for exploring distributions
library(HDInterval) #for HPD intervals
library(posterior)  #for posterior draws
library(coda)       #for diagnostics
library(bayesplot)  #for diagnostics
library(ggmcmc)     #for MCMC diagnostics
library(rstan)      #for interfacing with STAN
library(emmeans)    #for marginal means etc
library(broom)      #for tidying outputs
library(DHARMa)     #for residual diagnostics
library(tidybayes)  #for more tidying outputs
library(ggeffects)  #for partial plots
library(broom.mixed)#for tidying MCMC outputs
library(patchwork)  #for multiple plots
library(ggridges)   #for ridge plots
library(bayestestR) #for ROPE
library(see)        #for some plots
library(easystats)     #framework for stats, modelling and visualisation
library(modelsummary)
source('helperFunctions.R')


## \tikzstyle{HandLabel} = [font={\fontspec[Scale=1.1]{xkcd}}]
## \tikzstyle{Messy} = [decorate,decoration={random steps,segment length=3pt, amplitude=0.5pt}]
## \tikzset{%
## every node/.style={%
## draw=black,
## inner sep=1mm,
## outer sep=0,
## Messy, HandLabel,
## minimum size=2.5cm,
## minimum height=8mm,
## align=center,
## anchor=north,
## },
## Rnd/.style={%
## draw=black!90,
## fill=black!30,
## },
## Trt/.style={%
## %rounded corners,
## %Messy,
## draw=black,
## fill=none,
## %top color=blue!10,
## %bottom color=blue!30
## },
## Latent/.style={%
## %rounded corners,
## %Messy,
## draw=black!40,
## text=black!40,
## fill=none,
## %top color=blue!10,
## %bottom color=blue!30
## },
## Th/.style={%
## %rounded corners,
## draw=black!90
## },
## Control/.style={%
## rounded corners,
## draw=green!90,
## top color=green!10,
## bottom color=green!30,
## },
## Comment/.style={%
## draw=none,
## inner sep=0mm,
## outer sep=0mm,
## minimum height=5mm,
## align=right
## },
## }
## 
## \forestset{myst/.style={%
## for tree={%
## parent anchor=south,
## child anchor=north,
## l sep=1cm,
## s sep=0.5cm,
## edge path={\noexpand\path[\forestoption{edge},-{latex}]
## (!u.parent anchor) |- ($(!u.parent anchor)!.5!(.child anchor)$) -| (.child anchor)
## \forestoption{edge label};}
## }
## }
## }
## 
## \begin{forest} myst,
## [,phantom, s=1cm
## [FishID.1, Rnd, name=Random
## [{Low Salinity}, Trt, name=Trial
## [SMR, Trt, name=SMR]
## ]
## [{High Salinity}, Trt
## [SMR, Trt]
## ]
## [{Hypoxia}, Trt
## [SMR, Trt]
## ]
## ]
## [FishID.2, Rnd
## [{Low Salinity}, Trt
## [SMR, Trt]
## ]
## [{High Salinity}, Trt
## [SMR, Trt]
## ]
## [{Hypoxia}, Trt
## [SMR, Trt]
## ]
## ]
## [..., Comment]
## [FishID.n, Rnd
## [{Low Salinity}, Trt
## [SMR, Trt]
## ]
## [{High Salinity}, Trt
## [SMR, Trt]
## ]
## [{Hypoxia}, Trt
## [SMR, Trt]
## ]
## ]
## ]
## \node[left=1cm of Trial, Comment] (lTrial) {TRIAL};
## \node[left=1cm of SMR, Comment] (lSMR) {SMR};
## \node[Comment] at (lTrial |- Random.west) {FISHID};
## \end{forest}
## 

## ----readData, results='markdown', eval=TRUE----------------------------------
norin <- read_csv("../data/norin.csv", trim_ws = TRUE)


## -----------------------------------------------------------------------------
#| label: examinData
glimpse(norin)


## -----------------------------------------------------------------------------
## Explore the first 6 rows of the data
head(norin)


## -----------------------------------------------------------------------------
str(norin)


## -----------------------------------------------------------------------------
norin |> datawizard::data_codebook()


## -----------------------------------------------------------------------------
norin |> modelsummary::datasummary_skim()
norin |> modelsummary::datasummary_skim(by = "TRIAL")


## ----eda1, results='markdown', eval=TRUE, mhidden=TRUE------------------------
norin <- norin |>
  mutate(FISHID = factor(FISHID),
         TRIAL = factor(TRIAL))


## ----eda2, results='markdown', eval=TRUE, mhidden=TRUE------------------------
ggplot(norin, aes(y=CHANGE, x=TRIAL)) + geom_boxplot()


## ----eda3, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=7, fig.height=4----
ggplot(norin, aes(y=CHANGE, x=SMR_contr, shape=TRIAL, color=TRIAL)) +
    geom_smooth(method='lm') + geom_point()
ggplot(norin, aes(y=CHANGE, x=SMR_contr, shape=TRIAL, color=TRIAL)) +
  geom_smooth() + geom_point()


## ----eda4, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=7, fig.height=4----
ggplot(norin, aes(y=CHANGE, x=as.numeric(FISHID), color=TRIAL)) +
    geom_point() + geom_line()


## ----eda5, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=7, fig.height=4----
ggplot(norin, aes(y=CHANGE, x=MASS, color=TRIAL)) +
  geom_point() +
  geom_smooth(method='lm')


## ----fitModel1a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
norin_rstanarm <- stan_glmer(CHANGE ~  (1|FISHID)+TRIAL*scale(SMR_contr, scale = FALSE)+scale(MASS, scale = FALSE),
                             data = norin,
                             family = gaussian(),
                             iter = 5000,
                             warmup = 2000,
                             chains = 3,
                             thin = 5,
                             refresh = 0)


## ----fitModel1b, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
norin_rstanarm |> prior_summary()


## ----fitModel1c, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
2.5*sd(norin$CHANGE)


## ----fitModel1d, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
2.5*sd(norin$CHANGE)/apply(model.matrix(~TRIAL*SMR_contr+MASS, norin), 2, sd)


## ----fitModel1e, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
1/sd(norin$CHANGE)


## ----fitModel1f, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
norin_rstanarm1 <- update(norin_rstanarm,  prior_PD=TRUE)


## ----fitModel1g, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
norin_rstanarm1 |>
  ggpredict(~SMR_contr*TRIAL) |>
  plot(show_data=TRUE)
#OR
norin_rstanarm1 |>
  ggemmeans(~SMR_contr*TRIAL) |>
  plot(show_data=TRUE)


## ----fitModel1h, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
norin_rstanarm2 <- stan_glmer(CHANGE ~  (1|FISHID)+TRIAL*scale(SMR_contr, scale = FALSE)+offset(MASS),
                                data = norin,
                                family = gaussian(),
                                prior_intercept = normal(17, 35, autoscale = FALSE),
                                prior = normal(0, 70, autoscale = FALSE),
                                prior_aux=rstanarm::exponential(0.03, autoscale = FALSE),
                                prior_covariance = decov(1, 1, 1, 1),
                                prior_PD = TRUE,
                                iter = 5000,
                                warmup = 1000,
                                chains = 3,
                                thin = 5,
                                refresh = 0
                                )


## ----fitModel1i, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
norin_rstanarm2 |>
    ggpredict(~SMR_contr * TRIAL) |>
    plot(show_data = TRUE)


## ----fitModel1j, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE, dependson='fitModel1h'----
norin_rstanarm3 <- update(norin_rstanarm2,  prior_PD=FALSE)


## ----fitModel1j2, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE, dependson='fitModel1h'----
norin_rstanarm4 <- stan_glmer(CHANGE ~  (TRIAL|FISHID)+TRIAL*SMR_contr + offset(MASS),
                                data = norin,
                                family = gaussian(),
                                prior_intercept = normal(17, 35, autoscale = FALSE),
                                prior = normal(0, 70, autoscale = FALSE),
                                prior_aux=rstanarm::exponential(0.03, autoscale = FALSE),
                                prior_covariance = decov(1, 1, 1, 1),
                                iter = 5000,
                                warmup = 1000,
                                chains = 3,
                                thin = 5,
                                refresh = 0
                                )
preds <- norin_rstanarm4 |> posterior_predict(ndraws = 250,  summary = FALSE)
norin_resids <- createDHARMa(simulatedResponse = t(preds),
                            observedResponse = norin$CHANGE,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse = FALSE)
plot(norin_resids, quantreg = FALSE)

## Clearly an issue here!


## ----modelFit1k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
get_variables(norin_rstanarm3)
posterior_vs_prior(norin_rstanarm3, color_by='vs', group_by=TRUE,
                   facet_args=list(scales='free_y'),
                   regex_pars = "^.Intercept|TRIAL|SMR_contr|MASS|sigma|Sigma")


## ----modelFit1l, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
norin_rstanarm3 |>
    ggpredict(~SMR_contr * TRIAL) |>
    plot(show_data = TRUE)
##OR
norin_rstanarm3 |>
    ggemmeans(~SMR_contr * TRIAL) |>
    plot(show_data = TRUE)


## ----fitModel2a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE, paged.print=FALSE, tidy.opts = list(width.cutoff = 80), echo=c(-4,-6)----
norin_form <- bf(CHANGE ~  (1|FISHID)+TRIAL*SMR_contr+offset(MASS),
                   family = gaussian()
                   )
options(width=100)
norin_form |> get_prior(data=norin)
options(width=80)
norin_brm <- brm(norin_form,
                 data=norin,
                 iter = 5000,
                 warmup = 1000,
                 chains = 3, cores = 3,
                 thin = 5,
                 refresh = 0,
                 backend = 'cmdstanr')


## ----fitModel2h, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE, fig.width = 10, fig.height = 7----
norin |>
    group_by(TRIAL) |>
    summarise(median(CHANGE),
              mad(CHANGE))

norin |>
    group_by(TRIAL, FISHID) |>
    summarise(median = median(CHANGE),
              MAD = mad(CHANGE)) |>
    ungroup(FISHID) |>
    summarise(sd(median))

sd(norin$CHANGE)/apply(model.matrix(~TRIAL*scale(SMR_contr)+scale(MASS), norin)[, -1], 2, sd)
## mad(norin$CHANGE)/apply(model.matrix(~TRIAL*scale(SMR_contr)+scale(MASS), norin)[, -1], 2, mad)

standist::visualize("normal(53, 25)", xlim=c(-10,100))
standist::visualize("normal(0, 60)", xlim=c(-200,200))
standist::visualize("gamma(2, 1)", "gamma(35, 1)",
                    "student_t(3,0, 40)",
                    "cauchy(0, 5.8)",
                    xlim=c(-10,100))


## ----fitModel2h1, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE-----
norin_form <- bf(CHANGE ~ (1|FISHID) + TRIAL*scale(SMR_contr) + scale(MASS),
                     family = gaussian()
                   )
get_prior(norin_form, data = norin)

priors <- prior(normal(53, 25), class = 'Intercept') +
    ## prior(normal(0, 60), class = 'b', coef = 'TRIALHypoxia') +
    ## prior(normal(0, 70), class = 'b', coef = 'TRIALLowSalinity') +
    ## prior(normal(0, 54), class = 'b', coef = 'SMR_contr') +
    ## prior(normal(0, 4), class = 'b', coef = 'MASS') +
    prior(normal(0, 60), class = 'b') +
    prior(student_t(3,0,40), class = 'sd') +
    prior(student_t(3,0,40), class = 'sigma')
    ## prior(gamma(35, 1), class = 'sigma') +
    ## prior(cauchy(0, 5.8), class = 'sd')
norin_brm2 <- brm(norin_form,
                  data = norin,
                  prior = priors,
                  sample_prior = 'only',
                  iter = 5000,
                  warmup = 1000,
                  chains = 3, cores = 3,
                  thin = 5,
                  refresh = 0,
                  backend = "cmdstanr"
                  )


## ----partialPlot2h1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
norin_brm2 |>
  conditional_effects(~"SMR_contr:TRIAL") |>
  plot(points = TRUE)
norin_brm2 |>
    ggpredict(~SMR_contr*TRIAL) |>
    plot(show_data = TRUE)


## ----fitModel2h1b, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
norin_brm3 <- update(norin_brm2,
                     sample_prior = 'yes',
                     control = list(adapt_delta = 0.99),
                     refresh = 0,
                     cores = 3)
save(norin_brm3, file = '../ws/testing/norin_brm3')


## ----partialPlot2h1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
norin_brm3 |>
  conditional_effects(~"SMR_contr:TRIAL") |>
  plot(points = TRUE)
norin_brm3 |>
    ggpredict(~SMR_contr*TRIAL) |>
    plot(show_data = TRUE)


## ----posterior2h2, results='markdown', eval=TRUE------------------------------
norin_brm3 |> get_variables()
norin_brm3 |> hypothesis('TRIALHypoxia=0') |> plot()
norin_brm3 |> hypothesis('scaleSMR_contr=0') |> plot()


## ----posterior2h2a, results='markdown', eval=TRUE, fig.width = 7, fig.height = 5----
norin_brm3 |> SUYR_prior_and_posterior()


## ----fitModel2h3, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE-----
priors <- prior(normal(53, 35), class = 'Intercept') +
    ## prior(normal(0, 70), class = 'b', coef = 'TRIALHypoxia') +
    ## prior(normal(0, 70), class = 'b', coef = 'TRIALLowSalinity') +
    ## prior(normal(0, 54), class = 'b', coef = 'SMR_contr') +
    ## prior(normal(0, 60), class = 'b', coef = 'MASS') +
    prior(normal(0, 60), class = 'b') +
    prior(student_t(3,0,40), class = 'sigma') +
    prior(student_t(3, 0, 40), class = 'sd') +
    prior(lkj_corr_cholesky(1), class = 'cor')
norin_form <- bf(CHANGE ~ (TRIAL|FISHID) + TRIAL*scale(SMR_contr) + scale(MASS),
                 ## sigma ~ TRIAL*SMR_contr + offset(MASS) + (1|FISHID),
                 family = gaussian()
                 )
norin_brm4 <-  brm(norin_form,
                  data = norin,
                  prior = priors,
                  sample_prior = 'yes',
                  iter = 5000,
                  warmup = 1000,
                  chains = 3, cores = 3,
                  thin = 10,
                  refresh = 0,
                  control = list(adapt_delta=0.99),
                  backend = "cmdstanr"
                  )
save(norin_brm4, file = '../ws/testing/norin_brm4')


## ----posterior2k, results='markdown', eval=TRUE-------------------------------
norin_brm4 |> get_variables()
norin_brm4 |> hypothesis('TRIALHypoxia=0') |> plot()
norin_brm4 |> hypothesis('scaleSMR_contr=0') |> plot()


## ----posterior2k1, results='markdown', eval=TRUE, fig.width = 7, fig.height = 5----
norin_brm4 |> SUYR_prior_and_posterior()


## ----posterior2k2, results='markdown', eval=TRUE, fig.width=10, fig.height=4----
norin_brm4 |>
  posterior_samples() |>
  dplyr::select(-`lp__`) |>
  pivot_longer(everything(), names_to = 'key') |>
  filter(!str_detect(key, '^r')) |>
  mutate(Type = ifelse(str_detect(key, 'prior'), 'Prior', 'Posterior'),
         Class = case_when(
             str_detect(key, '(^b|^prior).*Intercept$') ~ 'Intercept',
             str_detect(key, 'b_TRIAL.*|prior_b_TRIAL.*') & !str_detect(key, '.*\\:.*') ~ 'TRIAL',
             str_detect(key, 'b_scaleSMR_contr|prior_b_SMR_contr') ~ 'SMR',
             str_detect(key, 'b_scaleMASS|prior_b_MASS') ~ 'MASS',
             str_detect(key, '.*\\:.*|prior_b_.*\\:.*') ~ 'Interaction',
             str_detect(key, 'sd') ~ 'sd',
             str_detect(key, '^cor|prior_cor') ~ 'cor',
             str_detect(key, 'sigma') ~ 'sigma'),
         Par = str_replace(key, 'b_', '')) |>
  ggplot(aes(x = Type,  y = value, color = Par)) +
  stat_pointinterval(position = position_dodge())+
  facet_wrap(~Class,  scales = 'free')



## ----fitModel2h3a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
(l.1 <- norin_brm3 |> loo())
(l.2 <- norin_brm4 |> loo())
loo_compare(l.1, l.2)


## ----modelValidation2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_mcmc()


## ----modelValidation2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
norin_brm4 |> mcmc_plot(type='trace')


## ----modelValidation2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
norin_brm4 |> mcmc_plot(type='acf_bar')


## ----modelValidation2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
norin_brm4 |> mcmc_plot(type='rhat_hist')


## ----modelValidation2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
norin_brm4 |> mcmc_plot(type='neff_hist')


## ----modelValidation2f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
norin_brm4 |> mcmc_plot(type='combo')
norin_brm4 |> mcmc_plot(type='violin')


## ----modelValidation2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
norin_brm4 |> get_variables()
pars <- norin_brm4 |> get_variables()
pars <- str_extract(pars, '^b_.*|^sigma$|^sd.*') |> na.omit()

norin_brm4$fit |>
    stan_trace(pars = pars)


## ----modelValidation2h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
norin_brm4$fit |>
    stan_ac(pars = pars)


## ----modelValidation2i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
norin_brm4$fit |> stan_rhat()


## ----modelValidation2j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
norin_brm4$fit |> stan_ess()


## ----modelValidation2k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
norin_brm4$fit |>
    stan_dens(separate_chains = TRUE, pars = pars)


## ----modelValidation2l, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=7----
norin_ggs <- norin_brm4 |> ggs(burnin = FALSE, inc_warmup = FALSE)
norin_ggs |> ggs_traceplot()


## ----modelValidation2m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=7----
ggs_autocorrelation(norin_ggs)


## ----modelValidation2n, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_Rhat(norin_ggs, scaling = 1.01)


## ----modelValidation2o, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_effective(norin_ggs)


## ----modelValidation2p, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_crosscorrelation(norin_ggs)


## ----modelValidation2q, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_grb(norin_ggs)


## ----modelValidation5a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation5b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
norin_brm4 |> pp_check(type = 'dens_overlay', ndraws = 100)


## ----modelValidation5c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## norin_brm4 |> pp_check(type = 'error_scatter_avg')


## ----modelValidation5e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
norin_brm4 |> pp_check(group = 'TRIAL', type = 'intervals')
norin_brm3 |> pp_check(group = 'TRIAL', type = 'intervals_grouped')
norin_brm3 |> pp_check(group = 'TRIAL', type = 'violin_grouped')


## ----modelValidation5g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(norin_brm2)


## ----modelValidation6aa, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=10----
norin_resids <- make_brms_dharma_res(norin_brm4, integerResponse = FALSE)
wrap_elements(~testUniformity(norin_resids)) +
               wrap_elements(~plotResiduals(norin_resids, form = factor(rep(1, nrow(norin))))) +
               wrap_elements(~plotResiduals(norin_resids)) +
               wrap_elements(~testDispersion(norin_resids))


## ----partialPlot2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
norin_brm4 |>
    conditional_effects() |>
    plot(ask = FALSE, points = TRUE, plot = FALSE) |>
    wrap_plots()


## ----partialPlot2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
norin_brm4 |>
    ggpredict(terms = c("SMR_contr", "TRIAL")) |>
    plot(show_data = TRUE)


## ----partialPlot2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
norin_brm4 |>
    ggemmeans(~SMR_contr*TRIAL) |>
    plot(show_data = TRUE)


## ----summariseModel2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
norin_brm4 |> summary()


## ----summariseModel2a1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=FALSE----
norin_sum <- summary(norin_brm4)


## ----summariseModel2i2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
norin_brm4 |>
  summarise_draws(
    median,
    ~ HDInterval::hdi(.x),
    rhat,
    ess_bulk,
    ess_tail
  )

## or if you want to exclude some parameters
norin_brm4 |>
  summarise_draws(
    median,
    ~ HDInterval::hdi(.x),
    rhat,
    ess_bulk,
    ess_tail
  ) |>
  filter(str_detect(variable, 'prior|^r_|^lp__', negate = TRUE))


## ----summariseModel2i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
norin_brm4 |> as_draws_df()
norin_brm4 |>
  as_draws_df() |>
  summarise_draws(
    median,
    ~ HDInterval::hdi(.x),
    rhat,
    ess_bulk,
    ess_tail
  )
## or if you want to exclude some parameters
norin_brm4 |>
  as_draws_df() |>
  summarise_draws(
    median,
    ~ HDInterval::hdi(.x),
    rhat,
    ess_bulk,
    ess_tail
  ) |>
  filter(str_detect(variable, 'prior|^r_|^lp__', negate = TRUE))


## ----summariseModel2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
norin_brm4$fit |>
    tidyMCMC(estimate.method = 'median',
             conf.int = TRUE,  conf.method = 'HPDinterval',
             rhat = TRUE, ess = TRUE)

## ----summariseModel2b1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
norin_tidy <- tidyMCMC(norin_brm4$fit, estimate.method='median',
                         conf.int=TRUE,  conf.method='HPDinterval',
                         rhat=TRUE, ess=TRUE)


## ----summariseModel2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
norin_brm4 |> get_variables()
norin_draw <- norin_brm4 |>
    gather_draws(`^b.Intercept$|b_.*|sd_.*|sigma`,  regex=TRUE)
norin_draw


## ----summariseModel2c1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
norin_draw |> median_hdci()


## ----summariseModel2c3, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
norin_gather <- norin_brm4 |>
    gather_draws(`b_Intercept|b_TREAT.*|sd_.*|sigma`,  regex=TRUE) %>%
  median_hdci()


## ----summariseModel2c4, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
norin_brm4 |>
    gather_draws(`b_Intercept|b_.*`, regex=TRUE) |>
    ggplot() +
    geom_vline(xintercept=0, linetype='dashed') +
    stat_slab(aes(x = .value, y = .variable,
                  fill = stat(ggdist::cut_cdf_qi(cdf,
                           .width = c(0.5, 0.8, 0.95),
                           labels = scales::percent_format())
                           )), color='black') +
    scale_fill_brewer('Interval', direction = -1, na.translate = FALSE)

norin_brm4 |>
    gather_draws(`b_Intercept|b_.*`, regex=TRUE) |>
    ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    facet_wrap(~.variable, scales='free')

norin_brm4 |>
    gather_draws(`b_Intercept|b_.*`, regex=TRUE) |>
    ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    theme_classic()


## ----summariseModel2j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
norin_brm4$fit |> plot(type='intervals')


## ----summariseModel2ka, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
norin_brm4 |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(.variable != 'b_Intercept') |>
    ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    facet_wrap(~.variable, scales='free')

norin_brm4 |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(.variable != 'b_Intercept') |>
    ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    geom_vline(xintercept = 0, linetype = 'dashed')


## ----summariseModel2c7, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
norin_brm4 |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(.variable != 'b_Intercept') |>
    ggplot() +
    geom_density_ridges(aes(x=.value, y = .variable), alpha=0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed')
##Or in colour
norin_brm4 |>
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
norin_brm4 |> tidy_draws()


## ----summariseModel2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
norin_brm4 |> spread_draws(`.*Intercept.*|^b_.*`,  regex=TRUE)


## ----summariseModel2f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
norin_brm4 |> posterior_samples() |> as_tibble()


## ----summariseModel2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
norin_brm4 |>
    bayes_R2(re.form = NA, summary=FALSE) |>
    median_hdci()
norin_brm4 |>
    bayes_R2(re.form = ~(1|FISHID), summary=FALSE) |>
    median_hdci()
norin_brm4 |>
    bayes_R2(re.form = ~(TRIAL|FISHID), summary=FALSE) |>
    median_hdci()


## ----summariseModel2k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
0.1 * sd(norin$CHANGE)
norin_brm4 |> rope(range = c(-3.38, 3.38))
rope(norin_brm4, range = c(-3.38, 3.38)) |> plot()



## -----------------------------------------------------------------------------
#| label: modelsummary
#| results: markup
#| eval: true
#| echo: true
#| cache: false
norin_brm4 |> modelsummary(
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
norin_brm4 |> modelplot(exponentiate = FALSE)


## ----posteriors1a, results='markdown', eval=TRUE, echo=TRUE, mhidden=TRUE-----
norin_brm4 |> emtrends(~TRIAL, var='SMR_contr')
norin_brm4 |> emtrends(~TRIAL, var='SMR_contr') |> pairs()
norin_brm4 |> emtrends(~TRIAL, var='SMR_contr') |> pairs() |>
  tidy_draws() |>
  summarise_draws(median, HDInterval::hdi,
                  Pl = ~ mean(.x < 0),
                  Pg = ~ mean(.x > 0))
norin_emt <- norin_brm4 |> emtrends(~TRIAL, var='SMR_contr') |> pairs() |> as.data.frame()


## ----posteriors1b, results='markdown', eval=TRUE, mhidden=TRUE----------------
norin_grid <- with(norin,  list(SMR_contr=Hmisc::smean.sdl(SMR_contr)))
norin_grid
norin_brm4 |> emmeans(~TRIAL|SMR_contr,  at=norin_grid) |> pairs()

norin_brm4 |>
    emmeans(~TRIAL|SMR_contr,  at=norin_grid) |>
    pairs() |>
    gather_emmeans_draws() |>
    median_hdci()
norin_brm4 |>
    emmeans(~TRIAL|SMR_contr,  at=norin_grid) |>
    pairs() |>
    tidy_draws() |>
    summarise_draws(median)


## ----posteriors1b2, results='markdown', eval=TRUE, mhidden=TRUE---------------
norin_em <- norin_brm4 |>
    emmeans(~TRIAL|SMR_contr, at=norin_grid) |>
    pairs() |>
    gather_emmeans_draws() |>
    mutate(Fit=.value)
norin_em
norin_em |> group_by(contrast) |> median_hdci(Fit)
norin_em |> group_by(contrast, SMR_contr) |> median_hdci(Fit)
## norin_em |>
##     group_by(contrast) |>
##     summarize(P=sum(Fit>0)/n())
norin_em |>
    group_by(contrast, SMR_contr) |>
    summarise(P=mean(Fit>0),
              P2 = 1 - P)


## ----summaryFigure1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=5----
norin_grid <- with(norin, list(SMR_contr = modelr::seq_range(SMR_contr, n = 100)))
newdata <- norin_brm4 |> emmeans(~SMR_contr|TRIAL, at = norin_grid) |>
    as.data.frame()
head(newdata)

ggplot(data = newdata, aes(y = emmean, x = SMR_contr)) +
  geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD, fill = TRIAL), alpha = 0.3) +
  geom_line(aes(, color = TRIAL)) +
  theme_classic() +
  theme(legend.position = c(0.99, 0.99),
        legend.justification = c(1, 1))

## The .fixed values are the predicted values without random effects
obs <- norin_brm4 |>
  augment() |>
  mutate(PartialObs=.fitted + .resid)

ggplot(data=newdata, aes(y=emmean, x=SMR_contr)) +
  geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD, fill=TRIAL), alpha=0.3) +
  geom_line(aes(, color=TRIAL)) +
  geom_point(data=obs,  aes(y=PartialObs,  color=TRIAL)) +
  ## geom_point(data=norin,  aes(y=CHANGE), color='gray') +
  theme_classic()


## ----summaryFigure1a2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=9, fig.height=4----
g1 <- ggplot(data=newdata, aes(y=emmean, x=SMR_contr)) +
  geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD, fill=TRIAL), alpha=0.3) +
  geom_line(aes(, color=TRIAL)) +
  geom_point(data=obs,  aes(y=PartialObs,  color=TRIAL)) +
  theme_classic()
g2 <- ggplot(data=newdata, aes(y=emmean, x=SMR_contr)) +
  geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD, fill=TRIAL), alpha=0.3) +
  geom_line(aes(, color=TRIAL)) +
  geom_point(data=norin,  aes(y=CHANGE,  color=TRIAL)) +
  theme_classic()
g1 + g2


## ----fitModels, results='markdown', eval=FALSE, mhidden=TRUE------------------
# norin = norin |> mutate(FISHID=factor(FISHID),
#                          TRIAL=factor(TRIAL))
# 
# ggplot(norin, aes(y=CHANGE, x=TRIAL)) + geom_boxplot()
# ggplot(norin, aes(y=CHANGE, x=SMR_contr, shape=TRIAL, color=TRIAL)) +
#     geom_smooth(method='lm') + geom_point()
# #ggplot(norin, aes(y=CHANGE, x=MASS, shape=TRIAL, color=TRIAL)) +
# #geom_smooth(method='lm') + geom_point()
# ggplot(norin, aes(y=CHANGE, x=as.numeric(FISHID), color=TRIAL)) +
#     geom_point() + geom_line()
# 
# #ggplot(norin, aes(y=MASS, x=TRIAL)) + geom_boxplot()
# ggplot(norin, aes(y=CHANGE, x=MASS, color=TRIAL)) + geom_point() + geom_smooth(method='lm')
# 
# norin_rstanarm = stan_glmer(CHANGE ~ (1|FISHID)+TRIAL*SMR_contr+MASS, data=norin,
#                             prior_PD=TRUE,
#                          iter=5000, warmup=2000, chains=3, thin=5, refresh=0)
# prior_summary(norin_rstanarm)
# 
# posterior_vs_prior(norin_rstanarm, color_by='vs', group_by=TRUE,
#                    facet_args=list(scales='free_y'), pars=c('(Intercept)'))
# ggpredict(norin_rstanarm, ~TRIAL*SMR_contr) |> plot(show_data=TRUE)
# 
# norin_rstanarm |> get_variables()
# plot(norin_rstanarm,  'mcmc_trace', regex_pars='^.Intercept|TRIAL|SMR|MASS|[sS]igma')
# plot(norin_rstanarm,  'mcmc_acf_bar', regex_pars='^.Intercept|TRIAL|SMR|MASS|[sS]igma')
# plot(norin_rstanarm,  'mcmc_rhat_hist', regex_pars='^.Intercept|TRIAL|SMR|MASS|[sS]igma')
# plot(norin_rstanarm,  'mcmc_neff_hist', regex_pars='^.Intercept|TRIAL|SMR|MASS|[sS]igma')
# 
# #norin_rstan1 = stan_glmer(CHANGE ~ (TRIAL|FISHID)+TRIAL*SMR_contr+MASS, data=norin,
# #                          iter=5000, warmup=2000, chains=3, thin=5, refresh=0, cores=3)
# norin_rstanarm1 = stan_glmer(CHANGE ~ (SMR_contr|FISHID) + TRIAL*SMR_contr+MASS, data=norin,
#                              prior_PD=FALSE,
#                           iter=5000, warmup=2000, chains=3, thin=5, refresh=0, cores=3)
# norin_rstanarm1 = update(norin_rstanarm1,  prior_PD=FALSE)
# 
# 
# 
# norin_rstanarm2 = stan_glmer(CHANGE ~ (TRIAL*SMR_contr|FISHID) + TRIAL*SMR_contr+MASS, data=norin,
#                              prior_PD=FALSE,
#                           iter=5000, warmup=2000, chains=3, thin=5, refresh=0, cores=3)
# 
# posterior_vs_prior(norin_rstanarm1, color_by='vs', group_by=TRUE,
#                    facet_args=list(scales='free_y'), pars=c('(Intercept)'))
# 
# ggpredict(norin_rstanarm1, ~TRIAL*SMR_contr) |> plot(show_data=TRUE)
# 
# norin_rstanarm1 |> get_variables()
# plot(norin_rstanarm1,  'mcmc_trace', regex_pars='^.Intercept|TRIAL|^SMR|MASS|[sS]igma')
# plot(norin_rstanarm1,  'mcmc_acf_bar', regex_pars='^.Intercept|TRIAL|^SMR|MASS|[sS]igma')
# plot(norin_rstanarm1,  'mcmc_rhat_hist', regex_pars='^.Intercept|TRIAL|^SMR|MASS|[sS]igma')
# plot(norin_rstanarm1,  'mcmc_neff_hist', regex_pars='^.Intercept|TRIAL|^SMR|MASS|[sS]igma')
# 
# 
# (l.1 <- loo(norin_rstanarm))
# (l.2 <- loo(norin_rstanarm1))
# loo_compare(l.1,  l.2)
# 
# 
# preds <- posterior_predict(norin_rstanarm,  nsamples=250,  summary=FALSE)
# norin_resids <- createDHARMa(simulatedResponse = t(preds),
#                             observedResponse = norin$CHANGE,
#                             fittedPredictedResponse = apply(preds, 2, median))
# plot(norin_resids)
# 
# 
# g=ggpredict(norin_rstanarm) |> plot()
# do.call('grid.arrange', g)
# 
# #ggemmeans(norin_rstan, ~TRIAL)
# 
# summary(norin_rstanarm)
# nms <- norin_rstanarm1 |> get_variables()
# wch <- grep('^.Intercept|TRIAL|^SMR|[sS]igma', nms)
# tidyMCMC(norin_rstanarm$stanfit,conf.int=TRUE, conf.method='HPDinterval',
#          rhat=TRUE, ess=TRUE, pars=nms[wch], estimate.method='median')
# 
# tidyMCMC(norin_rstanarm1$stanfit,conf.int=TRUE, conf.method='HPDinterval',
#          rhat=TRUE, ess=TRUE, pars=nms[wch], estimate.method='median')
# 
# 
# norin_grid = with(norin, list(SMR_contr=seq(min(SMR_contr),max(SMR_contr), len=100)))
# newdata = emmeans(norin_rstanarm, ~TRIAL|SMR_contr, at=norin_grid) |> as.data.frame()
# head(newdata)
# ggplot(newdata, aes(y=emmean, x=SMR_contr, color=TRIAL)) +
#     geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD, fill=TRIAL), alpha=0.3,color=NA) +
#     geom_line()
# 
# norin_grid = with(norin, list(SMR_contr=c(min(SMR_contr),mean(SMR_contr),max(SMR_contr))))
# 
# emmeans(norin_rstan, pairwise~TRIAL|SMR_contr, at=norin_grid)
# 
# norin_em = emmeans(norin_rstan, pairwise~TRIAL|SMR_contr, at=norin_grid)$contrast |>
#               gather_emmeans_draws() |>
#               mutate(Fit=.value)
# norin_em
# norin_em |> group_by(contrast) |> median_hdci(Fit)
# norin_em |> group_by(contrast, SMR_contr) |> median_hdci(Fit)
# ## norin_em |>
# ##     group_by(contrast) |>
# ##     summarize(P=sum(Fit>0)/n())
# norin_em |>
#     group_by(contrast, SMR_contr) |>
#     summarize(P=mean(Fit>0))
# 
# 
# bayes_R2(norin_rstanarm, re.form=NA) |> median_hdi()
# bayes_R2(norin_rstanarm, re.form=~(1|FISHID)) |> median_hdi()
# #bayes_R2(norin_rstan1, re.form=~(SMR_contr|FISHID)) |> median_hdi
# 


## ----fitModels.brms, results='markdown', eval=FALSE, mhidden=TRUE-------------
# norin = norin |> mutate(FISHID=factor(FISHID),
#                          TRIAL=factor(TRIAL))
# 
# ggplot(norin, aes(y=CHANGE, x=TRIAL)) + geom_boxplot()
# ggplot(norin, aes(y=CHANGE, x=SMR_contr, shape=TRIAL, color=TRIAL)) +
#     geom_smooth(method='lm') + geom_point()
# ggplot(norin, aes(y=CHANGE, x=MASS, shape=TRIAL, color=TRIAL)) +
#     geom_smooth(method='lm') + geom_point()
# ggplot(norin, aes(y=CHANGE, x=as.numeric(FISHID), color=TRIAL)) +
#     geom_point() + geom_line()
# 
# ##ggplot(norin, aes(y=MASS, x=TRIAL)) + geom_boxplot()
# ##ggplot(norin, aes(y=CHANGE, x=MASS, color=TRIAL)) + geom_point() + geom_smooth(method='lm')
# 
# norin |> group_by(TRIAL) |>
#     summarise(median(CHANGE),
#               mad(CHANGE))
# priors <- prior(normal(50, 20), class='Intercept') +
#     prior(normal(0, 10), class='b') +
#     prior(gamma(2,1), class='sigma') +
#     prior(gamma(2,1), class='sd')
# 
# norin_form <- bf(CHANGE ~ (1|FISHID)+TRIAL*SMR_contr+MASS,
#                  family=gaussian)
# 
# norin_brm1 = brm(norin_form,
#                  data=norin,
#                  prior = priors,
#                  sample_prior = 'yes',
#                  iter=5000, warmup=2000,
#                  chains=3, thin=5, refresh=0)
# 
# norin_brm1 |> get_variables()
# pars <- norin_brm1 |> get_variables()
# wch <- grepl('^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd', pars, perl=TRUE)
# 
# stan_trace(norin_brm1$fit, pars=pars[wch])
# stan_ac(norin_brm1$fit, pars=pars[wch])
# stan_rhat(norin_brm1$fit, pars=pars[wch])
# stan_ess(norin_brm1$fit, pars=pars[wch])
# 
# ##mcmc_plot(norin_brms,  type='trace',
# ##          regex_pars='^b.*|sigma|^sd')
# ##mcmc_plot(norin_brms,  type='trace',
# ##          regex_pars='^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd')
# ##mcmc_plot(norin_brms,  type='acf_bar',
# ##          regex_pars='^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd')
# ##mcmc_plot(norin_brms,  type='rhat_hist',
# ##          regex_pars='^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd')
# ##mcmc_plot(norin_brms,  type='neff_hist',
# ##          regex_pars='^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd')
# 
# preds <- posterior_predict(norin_brm1,  nsamples=250,  summary=FALSE)
# norin_resids <- createDHARMa(simulatedResponse = t(preds),
#                             observedResponse = norin$CHANGE,
#                             fittedPredictedResponse = apply(preds, 2, median),
#                             integerResponse =FALSE)
# plot(norin_resids)
# #norin_rstan1 = stan_glmer(CHANGE ~ (TRIAL|FISHID)+TRIAL*SMR_contr+MASS, data=norin,
# #                          iter=5000, warmup=2000, chains=3, thin=5, refresh=0, cores=3)
# norin_form <- bf(CHANGE ~ (TRIAL|FISHID) + TRIAL*SMR_contr+MASS,
#                  family=gaussian)
# norin_brm2 = brm(norin_form, data=norin,
#                  prior = priors,
#                  sample_prior = 'yes',
#                  iter=5000, warmup=2000,
#                  chains=3, thin=5, refresh=0, cores=3,
#                  control=list(adapt_delta=0.99))
# 
# norin_brm2 |> get_variables()
# 
# pars <- norin_brm2 |> get_variables()
# ## wch <- grepl('^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd', pars, perl=TRUE)
# wch <- grepl('^b_.*|[sS]igma|^sd_.*', pars, perl=TRUE)
# 
# stan_trace(norin_brm2$fit, pars=pars[wch])
# stan_ac(norin_brm2$fit, pars=pars[wch])
# stan_rhat(norin_brm2$fit)#, pars=pars[wch])
# stan_ess(norin_brm2$fit)#, pars=pars[wch])
# ##mcmc_plot(norin_brm2,  type='trace',
# ##          regex_pars='^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd')
# ##mcmc_plot(norin_brm2,  type='trace',
# ##          regex_pars='^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd')
# ##mcmc_plot(norin_brm2,  type='acf_bar',
# ##          regex_pars='^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd')
# ##mcmc_plot(norin_brm2,  type='rhat_hist',
# ##          regex_pars='^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd')
# ##mcmc_plot(norin_brm2,  type='neff_hist',
# ##          regex_pars='^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd')
# 
# (l.1 <- loo(norin_brm1))
# (l.2 <- loo(norin_brm2))
# loo_compare(l.1,  l.2)
# 
# 
# preds <- posterior_predict(norin_brm2,  nsamples=250,  summary=FALSE)
# norin_resids <- createDHARMa(simulatedResponse = t(preds),
#                             observedResponse = norin$CHANGE,
#                             fittedPredictedResponse = apply(preds, 2, median))
# plot(norin_resids)
# 
# g <- norin_brm2 |>
#     conditional_effects() |>
#     plot(points=TRUE, ask=FALSE)
# library(patchwork)
# g[[1]] + g[[2]] + g[[3]] + g[[4]]
# 
# 
# ##g=ggpredict(norin_brms1) |> plot
# ##library(patchwork)
# ##g[[1]] + g[[2]] + g[[3]]
# 
# ##do.call('grid.arrange', g)
# 
# ggemmeans(norin_brm2, ~TRIAL) |> plot()
# 
# summary(norin_brm2)
# 
# tidyMCMC(norin_brm2$fit,conf.int=TRUE, conf.method='HPDinterval',
#          rhat=TRUE, ess=TRUE, estimate.method='median') |>
#   slice(1:11)
# 
# pars <- norin_brm2 |> get_variables()
# wch <- grep('^b.Intercept|TRIAL|^b.*SMR|[sS]igma|^sd', pars)
# tidyMCMC(norin_brms1$fit,conf.int=TRUE, conf.method='HPDinterval',
#          rhat=TRUE, ess=TRUE, pars=pars[wch], estimate.method='median')
# 
# bayes_R2(norin_brm2, re.form=NA,  summary=FALSE) |>
#     median_hdci()
# bayes_R2(norin_brm2, re.form=~(1|FISHID), summary=FALSE) |>
#     median_hdci()
# bayes_R2(norin_brm2, re.form=~(TRIAL|FISHID), summary=FALSE) |>
#     median_hdci()
# 
# emmeans(norin_brm2, pairwise~TRIAL)
# 
# 
# norin_em <- norin_brm2 |>
#     emmeans(~TRIAL) |>
#     pairs() |>
#     gather_emmeans_draws() |>
#     mutate(Fit=.value)
# 
# norin_em |>
#   group_by(contrast) |>
#   median_hdi()
# 
# norin_em |>
#     ggplot() +
#     geom_vline(xintercept=0, linetype='dashed') +
#     stat_slab(aes(x=.value, y=contrast,
#                   fill = stat(ggdist::cut_cdf_qi(cdf,
#                             .width = c(0.5, 0.8, 0.95),
#                             labels = scales::percent_format())
#                             )), color='black') +
#     scale_fill_brewer('Interval', direction = -1, na.translate = FALSE) +
#     theme_bw()
# 
# norin_em |>
#     group_by(contrast) |>
#   summarize(P=mean(Fit>0))
# 
# 
# norin_grid <- with(norin, list(SMR_contr=c(min(SMR_contr),
#                                            mean(SMR_contr),
#                                            max(SMR_contr))))
# 
# norin_em <- norin_brm2 |>
#     emmeans(~TRIAL|SMR_contr, at=norin_grid) |>
#     pairs() |>
#     gather_emmeans_draws()
# 
# norin_em |> head()
# norin_em |>
#     group_by(contrast, SMR_contr) |>
#     median_hdi()
# 
# norin_em |>
#     group_by(contrast, SMR_contr) |>
#     summarize(P=mean(.value>0))
# 
# norin_grid <- with(norin, list(SMR_contr=modelr::seq_range(SMR_contr, n=100)))
# newdata <- norin_brm2 |>
#     emmeans(~SMR_contr|TRIAL, at=norin_grid) |>
#     as.data.frame()
# head(newdata)
# partial.obs <- norin |>
#     mutate(Pred = predict(norin_brm2, re.form = NA, summary=TRUE)[,'Estimate'],
#            Resid = resid(norin_brm2)[,'Estimate'],
#            Obs = Pred + Resid)
# ggplot(newdata, aes(y=emmean, x=SMR_contr, color=TRIAL)) +
#     geom_point(data=partial.obs, aes(y=Obs)) +
#     ##geom_point(data=partial.obs, aes(y=CHANGE), shape=2) +
#     geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD, fill=TRIAL), alpha=0.3,color=NA) +
#     geom_line()

