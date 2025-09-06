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

library(tidyverse)  #for data wrangling etc
library(rstanarm)   #for fitting models in STAN
library(cmdstanr)   #for cmdstan
library(brms)       #for fitting models in STAN
library(standist)   #for exploring distributions
library(HDInterval) #for HPD intervals
library(posterior)  #for posterior draws
library(coda)       #for diagnostics
library(bayesplot)  #for diagnostics
library(ggmcmc)     #for diagnostics
library(rstan)      #for interfacing with STAN
library(DHARMa)     #for residual diagnostics
library(emmeans)    #for marginal means etc
library(broom)      #for tidying outputs
library(broom.mixed) #for tidying MCMC outputs
library(tidybayes)  #for more tidying outputs
library(ggeffects)  #for partial plots
library(patchwork)  #for multiple figures
library(bayestestR) #for ROPE
library(see)        #for some plots
library(ggridges)   #for ridge plots
library(easystats)     #framework for stats, modelling and visualisation
library(modelsummary)
source('helperFunctions.R')


## -----------------------------------------------------------------------------
#| label: readData
tobacco <- read_csv("../data/tobacco.csv", trim_ws = TRUE)


## -----------------------------------------------------------------------------
#| label: examinData
glimpse(tobacco)


## -----------------------------------------------------------------------------
## Explore the first 6 rows of the data
head(tobacco)


## -----------------------------------------------------------------------------
str(tobacco)


## -----------------------------------------------------------------------------
tobacco |> datawizard::data_codebook()


## -----------------------------------------------------------------------------
tobacco |> modelsummary::datasummary_skim()
tobacco |> modelsummary::datasummary_skim(by = "TREATMENT")


## ----processData, results='markdown', eval=TRUE-------------------------------
tobacco <- tobacco |> mutate(LEAF = factor(LEAF),
                             TREATMENT = factor(TREATMENT))
tobacco |> head()


## ----tobaccoEDA2, results='markdown', eval=TRUE, mhidden=TRUE-----------------
ggplot(tobacco,  aes(y = NUMBER,  x = TREATMENT)) +
  geom_boxplot()


## ----tobaccoEDA3, results='markdown', eval=TRUE, mhidden=TRUE-----------------
ggplot(tobacco,  aes(y = NUMBER,  x = as.numeric(LEAF))) +
  geom_line(aes(linetype = TREATMENT))

## If we want to retain the original LEAF labels
ggplot(tobacco,  aes(y = NUMBER,  x = as.numeric(LEAF))) +
  geom_blank(aes(x = LEAF)) +
  geom_line(aes(linetype = TREATMENT))


## ----tobaccoEDA4, results='markdown', eval=TRUE, mhidden=TRUE-----------------
ggplot(tobacco,  aes(y = NUMBER,  x = TREATMENT,  group = LEAF)) +
  geom_point() +
  geom_line(aes(x = as.numeric(TREATMENT)))


## ----fitModel1a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
tobacco_rstanarm <- stan_glmer(NUMBER ~ (1|LEAF) + TREATMENT,
                               data = tobacco,
                               family = gaussian(),
                               iter = 5000,
                               warmup = 2000,
                               chains = 3,
                               thin = 5,
                               refresh = 0)


## ----fitModel1b, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
tobacco_rstanarm |> prior_summary()


## ----fitModel1c, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
2.5*sd(tobacco$NUMBER)


## ----fitModel1d, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
2.5*sd(tobacco$NUMBER)/apply(model.matrix(~TREATMENT, tobacco), 2, sd)


## ----fitModel1e, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
1/sd(tobacco$NUMBER)


## ----fitModel1f, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
tobacco_rstanarm1 <- update(tobacco_rstanarm,  prior_PD=TRUE)


## ----fitModel1g, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
ggpredict(tobacco_rstanarm1) |> plot(show_data = TRUE)


## ----fitModel1h, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE------
tobacco_rstanarm2 <- stan_glmer(NUMBER ~ (1|LEAF) + TREATMENT,
                                data = tobacco,
                                family = gaussian(),
                                prior_intercept = normal(35, 7, autoscale = FALSE),
                                prior = normal(0, 13, autoscale = FALSE),
                                prior_aux=rstanarm::exponential(0.15, autoscale = FALSE),
                                prior_covariance = decov(1, 1, 1, 1),
                                prior_PD = TRUE,
                                iter = 5000,
                                warmup = 1000,
                                chains = 3,
                                thin = 5,
                                refresh = 0
                                )


## ----fitModel1i, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE-----
tobacco_rstanarm2 |>
    ggpredict() |>
    plot(show_data = TRUE)


## ----fitModel1j, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE, dependson='fitModel1h'----
tobacco_rstanarm3 <- update(tobacco_rstanarm2,  prior_PD=FALSE)


## ----modelFit1k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
posterior_vs_prior(tobacco_rstanarm3, color_by='vs', group_by=TRUE,
                   facet_args=list(scales='free_y'))


## ----modelFit1l, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggemmeans(tobacco_rstanarm3,  ~TREATMENT) |> plot(show_data=TRUE)
ggpredict(tobacco_rstanarm3,  ~TREATMENT) |> plot(show_data=TRUE)


## ----fitModel2a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE, paged.print=FALSE, tidy.opts = list(width.cutoff = 80), echo=c(-4,-6)----
tobacco_form <- bf(NUMBER ~ (1|LEAF) + TREATMENT,
                   family = gaussian()
                   )
options(width=100)
tobacco_form |> get_prior(data=tobacco)
options(width=80)
## tobacco_brm <- brm(tobacco_form,
##                   data=tobacco,
##                   iter = 5000,
##                   warmup = 1000,
##                   chains = 3,
##                   thin = 5,
##                   refresh = 0)


## ----fitModel2h, results='markdown', eval=TRUE, mhidden=TRUE, cache=FALSE, fig.width = 10, fig.height = 7----
tobacco |>
    group_by(TREATMENT) |>
    summarise(median(NUMBER),
              mad(NUMBER))
standist::visualize("normal(35,3)", xlim=c(-10,100))
standist::visualize("normal(0, 10)", xlim=c(-20,20))
standist::visualize(
            "student_t(3,0,10)",
            xlim=c(-30,50))


## ----fitModel2h1, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE-----
priors <- prior(normal(35, 3), class = 'Intercept') +
    prior(normal(0, 10), class = 'b') +
    prior(student_t(3, 0, 5), class = 'sigma') +
    prior(student_t(3, 0, 5), class = 'sd')
tobacco_form <- bf(NUMBER ~ (1|LEAF) + TREATMENT,
                     family = gaussian()
                   )
tobacco_brm2 <- brm(tobacco_form,
                  data = tobacco,
                  prior = priors,
                  sample_prior = 'only',
                  iter = 5000,
                  warmup = 2500,
                  chains = 3, cores = 3,
                  thin = 5,
                  refresh = 0,
                  backend = "cmdstanr"
                  )


## ----fitModel2h1a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE, fig.width = 8, fig.height = 3----
tobacco_np <- nuts_params(tobacco_brm2)
tobacco_mcmc <- as.array(tobacco_brm2)
mcmc_parcoord(x = tobacco_mcmc, np = tobacco_np) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
tobacco_brm2 |> mcmc_parcoord(np = tobacco_np) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
tobacco_brm2 |> mcmc_parcoord(regex_pars = "^b.*|^r.*") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## ----partialPlot2h1a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tobacco_brm2 |>
  conditional_effects() |>
  plot(points = TRUE)
tobacco_brm2 |>
    ggpredict() |>
    plot(show_data = TRUE)


## ----fitModel2h1b, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
tobacco_brm3 <- update(tobacco_brm2,
                       sample_prior = 'yes',
                       control = list(adapt_delta = 0.99, max_treedepth =  20),
                       refresh = 0, cores = 3)
save(tobacco_brm3, file = '../ws/testing/tobacco_brm3')


## ----partialPlot2h1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tobacco_brm3 |>
  conditional_effects() |>
  plot(points = TRUE)
tobacco_brm3 |>
    ggpredict() |>
    plot(show_data = TRUE)


## ----posterior2h2, results='markdown', eval=TRUE------------------------------
tobacco_brm3 |> get_variables()
tobacco_brm3 |> hypothesis('TREATMENTWeak=0') |> plot()


## ----posterior2h2a, results='markdown', eval=TRUE, fig.width = 7, fig.height = 5----
tobacco_brm3 |> SUYR_prior_and_posterior()


## ----fitModel2h3, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE-----
priors <- prior(normal(35, 3), class = 'Intercept') +
    prior(normal(0, 10), class = 'b') +
    prior(student_t(3, 0, 5), class = 'sigma') +
    prior(student_t(3, 0, 5), class = 'sd') +
    prior(lkj_corr_cholesky(1), class = 'cor')
tobacco_form <- bf(NUMBER ~ (TREATMENT|LEAF) + TREATMENT,
                     family = gaussian()
                   )

tobacco_brm4 <-  brm(tobacco_form,
                  data = tobacco,
                  prior = priors,
                  sample_prior = 'yes',
                  iter = 5000,
                  warmup = 1000,
                  chains = 3,
                  thin = 5,
                  refresh = 0,
                  control = list(adapt_delta=0.99),
                  backend = "cmdstanr"
                  )
save(tobacco_brm4, file = '../ws/testing/tobacco_brm4')


## ----posterior2k, results='markdown', eval=TRUE-------------------------------
tobacco_brm4 |> get_variables()
tobacco_brm4 |> hypothesis('TREATMENTWeak=0') |> plot()


## ----posterior2k1, results='markdown', eval=TRUE, fig.width = 7, fig.height = 5----
tobacco_brm4 |> SUYR_prior_and_posterior()


## ----posterior2k2, results='markdown', eval=TRUE, fig.width=10, fig.height=4----
tobacco_brm4 |>
  posterior_samples() |>
  dplyr::select(-`lp__`) |>
  pivot_longer(everything(), names_to = 'key') |>
  filter(!str_detect(key, '^r')) |>
  mutate(Type = ifelse(str_detect(key, 'prior'), 'Prior', 'Posterior'),
         ## Class = ifelse(str_detect(key, 'Intercept'),  'Intercept',
         ##         ifelse(str_detect(key, 'b'),  'b', 'sigma')),
         Class = case_when(
             str_detect(key, '(^b|^prior).*Intercept$') ~ 'Intercept',
             str_detect(key, 'b_TREATMENT|prior_b') ~ 'TREATMENT',
             str_detect(key, 'sd') ~ 'sd',
             str_detect(key, '^cor|prior_cor') ~ 'cor',
             str_detect(key, 'sigma') ~ 'sigma'),
         Par = str_replace(key, 'b_', '')) |>
  ggplot(aes(x = Type,  y = value, color = Par)) +
  stat_pointinterval(position = position_dodge())+
  facet_wrap(~Class,  scales = 'free')



## ----fitModel2h3a, results='markdown', eval=TRUE, mhidden=TRUE, cache=TRUE----
(l.1 <- tobacco_brm3 |> loo())
(l.2 <- tobacco_brm4 |> loo())
loo_compare(l.1, l.2)


## ----modelValidation2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_mcmc()


## ----modelValidation2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
tobacco_brm3 |> mcmc_plot(type='trace')


## ----modelValidation2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
tobacco_brm3 |> mcmc_plot(type='acf_bar')


## ----modelValidation2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
tobacco_brm3 |> mcmc_plot(type='rhat_hist')


## ----modelValidation2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
tobacco_brm2 |> mcmc_plot(type='neff_hist')


## ----modelValidation2f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
tobacco_brm3 |> mcmc_plot(type='combo')
tobacco_brm3 |> mcmc_plot(type='violin')


## ----modelValidation2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
tobacco_brm3 |> get_variables()
pars <- tobacco_brm3 |> get_variables()
pars <- str_extract(pars, '^b_.*|^sigma$|^sd.*') |> na.omit()

tobacco_brm3$fit |> stan_trace(pars = pars)


## ----modelValidation2h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
tobacco_brm3$fit |>
    stan_ac(pars = pars)


## ----modelValidation2i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
tobacco_brm3$fit |> stan_rhat()


## ----modelValidation2j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
tobacco_brm2$fit |> stan_ess()


## ----modelValidation2k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
tobacco_brm3$fit |>
    stan_dens(separate_chains = TRUE, pars = pars)


## ----modelValidation2l, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=7----
tobacco_ggs <- tobacco_brm3 |> ggs(burnin = FALSE, inc_warmup = FALSE)
tobacco_ggs |> ggs_traceplot()


## ----modelValidation2m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=7----
ggs_autocorrelation(tobacco_ggs)


## ----modelValidation2n, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_Rhat(tobacco_ggs)


## ----modelValidation2o, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_effective(tobacco_ggs)


## ----modelValidation2p, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_crosscorrelation(tobacco_ggs)


## ----modelValidation2q, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
ggs_grb(tobacco_ggs)


## ----modelValidation5a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation5b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
tobacco_brm3 |> pp_check(type = 'dens_overlay', ndraws = 100)


## ----modelValidation5c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#tobacco_brm3 |> pp_check(type = 'error_scatter_avg')


## ----modelValidation5e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
tobacco_brm3 |> pp_check(group = 'TREATMENT', type = 'intervals')
tobacco_brm3 |> pp_check(group = 'TREATMENT', type = 'intervals_grouped')
tobacco_brm3 |> pp_check(group = 'TREATMENT', type = 'violin_grouped')


## ----modelValidation5g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(tobacco_brm2)


## ----modelValidation6a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
preds <- tobacco_brm4 |> posterior_predict(ndraws = 250,  summary = FALSE)
tobacco_resids <- createDHARMa(simulatedResponse = t(preds),
                            observedResponse = tobacco$NUMBER,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse = FALSE)
plot(tobacco_resids, quantreg = FALSE)


## ----modelValidation6aa, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=10----
tobacco_resids <- make_brms_dharma_res(tobacco_brm3, integerResponse = FALSE)
wrap_elements(~testUniformity(tobacco_resids)) +
               wrap_elements(~plotResiduals(tobacco_resids, form = factor(rep(1, nrow(tobacco))))) +
               wrap_elements(~plotResiduals(tobacco_resids, quantreg = FALSE)) +
               wrap_elements(~testDispersion(tobacco_resids))


## ----partialPlot2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tobacco_brm3 |>
    conditional_effects() |>
    plot(points = TRUE)


## ----partialPlot2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tobacco_brm3 |>
    ggpredict() |>
    plot(show_data = TRUE)


## ----partialPlot2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tobacco_brm3 |>
    ggemmeans(~TREATMENT) |>
    plot(show_data = TRUE) +
    geom_point(data = tobacco, aes(y = NUMBER, x = as.numeric(TREATMENT)))


## ----partialPlot2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
Partial.obs <- tobacco_brm3$data |>
    mutate(Pred = predict(tobacco_brm3)[,'Estimate'],
           Resid = resid(tobacco_brm3)[,'Estimate'],
           Obs = Pred + Resid)

tobacco_brm3 |>
    fitted_draws(newdata = tobacco) |>
    median_hdci() |>
    ggplot(aes(x = TREATMENT, y = .value)) +
    geom_pointrange(aes(ymin = .lower, ymax = .upper)) +
    geom_line() +
    geom_point(data = Partial.obs,  aes(y = Obs,  x = TREATMENT), color = 'red',
               position = position_nudge(x = 0.1)) +
    geom_point(data = tobacco,  aes(y = NUMBER,  x = TREATMENT),
               position = position_nudge(x = 0.05))

tobacco_brm3 |>
    epred_draws(newdata = tobacco) |>
    ggplot() +
    geom_violin(data = tobacco, aes(y = NUMBER, x = TREATMENT), fill = 'blue', alpha = 0.2) +
    geom_point(data = tobacco, aes(y = NUMBER, x = TREATMENT),
               position = position_jitter(width = 0.1, height = 0)) +
    geom_violin(aes(y = .epred, x = TREATMENT), fill = 'orange', alpha = 0.2) +
    theme_bw()


## ----summariseModel2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tobacco_brm3 |> summary()


## ----summariseModel2a1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=FALSE----
tobacco_sum <- summary(tobacco_brm3)


## ----summariseModel2i2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tobacco_brm3 |>
  summarise_draws(
    median,
    ~ HDInterval::hdi(.x),
    N =  ~ length(.x),
    Pl =  ~mean(.x < 0),
    Pg =  ~mean(.x > 0),
    rhat,
    ess_bulk,
    ess_tail
  )

## or if you want to exclude some parameters
tobacco_brm3 |>
  summarise_draws(
    median,
    ~ HDInterval::hdi(.x),
    rhat,
    ess_bulk,
    ess_tail
  ) |>
  filter(str_detect(variable, 'prior|^r_|^lp__', negate = TRUE))


## ----summariseModel2i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tobacco_brm3 |> as_draws_df()
tobacco_brm3 |>
  as_draws_df() |>
  dplyr::select(matches("^b_.*|^sigma$|^sd_.*")) |>
  summarise_draws(
    median,
    ~ HDInterval::hdi(.x),
    Pg = ~ mean(.x > 0),
    Pl = ~ mean(.x < 0),
    rhat,
    ess_bulk,
    ess_tail
  )
## or if you want to exclude some parameters
tobacco_brm3 |>
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
tobacco_brm3$fit |>
    tidyMCMC(estimate.method = 'median',
             conf.int = TRUE,  conf.method = 'HPDinterval',
             rhat = TRUE, ess = TRUE)

## ----summariseModel2b1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
tobacco_tidy <- tidyMCMC(tobacco_brm3$fit, estimate.method='median',
                         conf.int=TRUE,  conf.method='HPDinterval',
                         rhat=TRUE, ess=TRUE)


## ----summariseModel2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tobacco_brm3 |> get_variables()
tobacco_draw <- tobacco_brm3 |>
    gather_draws(`b.Intercept.*|b_TREAT.*|sd_.*|sigma`,  regex=TRUE)
tobacco_draw


## ----summariseModel2c1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tobacco_draw |> median_hdci()


## ----summariseModel2c3, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
tobacco_gather <- tobacco_brm3 |>
    gather_draws(`b_Intercept.*|b_TREAT.*|sd_.*|sigma`,  regex=TRUE) |>
  median_hdci()


## ----summariseModel2c4, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
tobacco_brm3 |>
    gather_draws(`b_Intercept.*|b_TREAT.*`, regex=TRUE) |>
    ggplot() +
    geom_vline(xintercept=0, linetype='dashed') +
    stat_slab(aes(x = .value, y = .variable,
                  fill = stat(ggdist::cut_cdf_qi(cdf,
                           .width = c(0.5, 0.8, 0.95),
                           labels = scales::percent_format())
                           )), color='black') +
    scale_fill_brewer('Interval', direction = -1, na.translate = FALSE)

tobacco_brm3 |>
  gather_draws(`.Intercept.*|.*TREAT.*`, regex=TRUE) |>
  ggplot() +
  stat_halfeye(aes(x=.value,  y=.variable)) +
  facet_wrap(~.variable, scales='free')

tobacco_brm3 |>
  gather_draws(`.Intercept.*|.*TREAT.*`, regex=TRUE) |>
  ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    theme_classic()


## ----summariseModel2j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tobacco_brm3$fit |> plot(type='intervals')


## ----summariseModel2ka, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
tobacco_brm3 |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(.variable != 'b_Intercept') |>
    ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    facet_wrap(~.variable, scales='free')

tobacco_brm3 |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(.variable != 'b_Intercept') |>
    ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    geom_vline(xintercept = 0, linetype = 'dashed')


## ----summariseModel2c7, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
tobacco_brm3 |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(.variable != 'b_Intercept') |>
    ggplot() +
    geom_density_ridges(aes(x=.value, y = .variable), alpha=0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed')
##Or in colour
tobacco_brm3 |>
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
tobacco_brm3 |> tidy_draws()


## ----summariseModel2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tobacco_brm3 |> spread_draws(`.*Intercept.*|.*TREAT.*`,  regex=TRUE)


## ----summariseModel2f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tobacco_brm3 |> posterior_samples() |> as_tibble()


## ----summariseModel2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tobacco_brm3 |>
    bayes_R2(re.form = NA, summary=FALSE) |>
    median_hdci()
tobacco_brm3 |>
    bayes_R2(re.form = ~(1|LEAF), summary=FALSE) |>
    median_hdci()
tobacco_brm3 |>
     bayes_R2(re.form = ~(TREATMENT|LEAF), summary=FALSE) |>
     median_hdci()


## ----summariseModel2k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
0.1 * sd(tobacco$NUMBER)
tobacco_brm3 |> rope(range = c(-0.65, 0.65))
rope(tobacco_brm3, range = c(-0.65, 0.65)) |> plot()

## Or based on fractional scale
tobacco_brm3 |> emmeans(~TREATMENT) |>
    gather_emmeans_draws() |>
    group_by(.draw) |>
    arrange(desc(TREATMENT)) |>
    summarise(Diff = 100*(exp(diff(log(.value))) -1)) |>
    rope(range = c(-10,10))



## -----------------------------------------------------------------------------
#| label: modelsummary
#| results: markup
#| eval: true
#| echo: true
#| cache: false
tobacco_brm3 |> modelsummary(
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
tobacco_brm3 |> modelplot(exponentiate = TRUE)


## ----predictions2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
tobacco_brm3 |> emmeans(~TREATMENT) |>
    gather_emmeans_draws() |>
    group_by(.draw) |>
    arrange(desc(TREATMENT)) |>
    summarise(Diff = 100*(exp(diff(log(.value))) -1)) |>
    summarise_draws(
        median,
        ~ HDInterval::hdi(.x),
        rhat,
        ess_bulk,
        ess_tail
        )

## Or via gather and pivot
newdata <- tobacco_brm3 |>
    emmeans(~TREATMENT) |>
    gather_emmeans_draws() |>
    pivot_wider(names_from=TREATMENT,values_from=.value) |>
    mutate(Eff = Strong - Weak,
           PEff = 100*Eff/Weak)
newdata |> median_hdci(PEff)
newdata |> summarise(P = mean(PEff>0))
newdata |> summarise(P = mean(PEff>20))
newdata |>
    dplyr::select(-.chain, -.iteration) |>
    hypothesis('PEff>20')

newdata <- tobacco_brm3 |> emmeans(~TREATMENT) |> as.data.frame()
head(newdata)
ggplot(newdata, aes(y=emmean, x=TREATMENT)) +
    geom_pointrange(aes(ymin=lower.HPD, ymax=upper.HPD)) +
    theme_bw()

tobacco_brm3 |>
    emmeans(~TREATMENT) |>
    gather_emmeans_draws() |>
    ggplot() +
    geom_density_ridges(aes(x = .value, y = TREATMENT), alpha = 0.5, fill = 'orange') +
    scale_x_continuous("Average number of lesions") +
    theme_bw()



## ----predictions2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
newdat <- tobacco |> tidyr::expand(TREATMENT)
newdata <- tobacco_brm3 |>
    brms::posterior_epred(newdat, re_formula = NA) |>
    as.data.frame() |>
    rename_with(~as.character(newdat$TREATMENT)) |>
    mutate(Eff = Strong - Weak,
           PEff = 100*Eff/Weak)
head(newdata)
newdata |> median_hdci(PEff)
newdata |> summarise(P = mean(PEff>0))
newdata |> summarise(P = mean(PEff>20))


## ----predictions2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
newdata <- tobacco_brm3 |>
    emmeans(~TREATMENT) |>
    pairs() |>
    gather_emmeans_draws()
newdata |> median_hdci()

## OR on percentage scale
newdata <- tobacco_brm3 |>
    emmeans(~TREATMENT) |>
    regrid(trans = 'log') |>
    pairs() |>
    regrid() |>
    gather_emmeans_draws() |>
    mutate(.value = (.value - 1) * 100)
newdata |> median_hdci()
newdata |> summarise(P = mean(.value>0))
newdata |> summarise(P = mean(.value>20))

