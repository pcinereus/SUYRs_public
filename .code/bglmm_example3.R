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
## minimum size=1cm,
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
## [Tree, Trt, name=Situation
## [{tree1}, Rnd, name=Bird
## [Nov, Trt, name=Month]
## [Jan, Trt]
## ]
## [{tree2}, Rnd
## [Nov, Trt]
## [Jan, Trt]
## ]
## [{tree ...}, Rnd
## [Nov, Trt]
## [Jan, Trt]
## ]
## ]
## [Nest-box, Trt, name=NB
## [nest-box1, Rnd [Nov, Trt][Jan, Trt]]
## [nest-box2, Rnd [Nov, Trt][Jan, Trt]]
## [{nest-box ...}, Rnd [Nov, Trt][Jan, Trt]]
## ]
## [Other, Trt, name=Ot
## [other1, Rnd [Nov, Trt][Jan, Trt]]
## [other2, Rnd [Nov, Trt][Jan, Trt]]
## [{other ...}, Rnd [Nov, Trt][Jan, Trt]]
## ]
## ]
## \node[left=1cm of Month.west, Comment,anchor=east] (lMonth) {Month};
## \node[Comment,anchor=east] at (lMonth.east |- Bird.west) (lBird) {Bird};
## \node[Comment,anchor=east] at (lMonth.east |- Situation.west) (lSituation) {Situation};
## \node [Comment] at ($(NB) !0.5! (Ot)$) {....};
## \end{forest}
## 
## 

## ----readData, results='markdown', eval=TRUE----------------------------------
starling <- read_csv('../data/starling_full.csv', trim_ws = TRUE)


## -----------------------------------------------------------------------------
#| label: examinData
glimpse(starling)


## -----------------------------------------------------------------------------
## Explore the first 6 rows of the data
head(starling)


## -----------------------------------------------------------------------------
str(starling)


## -----------------------------------------------------------------------------
starling |> datawizard::data_codebook()


## -----------------------------------------------------------------------------
starling |> modelsummary::datasummary_skim()
starling |> modelsummary::datasummary_skim(by = c("SITUATION", "MONTH"))


## ----processData, results='markdown', eval=TRUE-------------------------------
starling <- starling |>
    mutate(BIRD = factor(BIRD),
           SITUATION = factor(SITUATION),
           MONTH = factor(MONTH, levels=c('Nov', 'Jan')))


## ----eda1, results='markdown', eval=TRUE--------------------------------------
ggplot(starling, aes(y=MASS, x=MONTH)) +
    geom_boxplot() +
    facet_grid(~SITUATION)
ggplot(starling, aes(y=MASS, x=SITUATION, fill = MONTH)) +
    geom_boxplot()
## Better still
ggplot(starling, aes(y=MASS, x=MONTH, group=BIRD)) +
    geom_point() +
    geom_line() +
    facet_grid(~SITUATION)


## ----fitModel1a, results='markdown', eval=TRUE, cache=TRUE--------------------
starling_rstanarm <- stan_glmer(MASS ~ MONTH*SITUATION+(1|BIRD),
                                data = starling,
                                family = gaussian(),
                                iter = 5000,
                                warmup = 2000,
                                chains = 3,
                                thin = 5,
                                refresh = 0)


## ----fitModel1b, results='markdown', eval=TRUE, cache=FALSE-------------------
starling_rstanarm %>% prior_summary()


## ----fitModel1c, results='markdown', eval=TRUE, cache=FALSE-------------------
2.5*sd(starling$MASS)


## ----fitModel1d, results='markdown', eval=TRUE, cache=FALSE-------------------
2.5*sd(starling$MASS)/apply(model.matrix(~MONTH*SITUATION, starling)[,-1], 2, sd)


## ----fitModel1e, results='markdown', eval=TRUE, cache=TRUE--------------------
1/sd(starling$MASS)


## ----fitModel1f, results='markdown', eval=TRUE, cache=TRUE--------------------
starling_rstanarm1 <- update(starling_rstanarm,  prior_PD=TRUE)


## ----fitModel1g, results='markdown', eval=TRUE, cache=FALSE-------------------
starling_rstanarm1 %>%
    ggpredict(~MONTH*SITUATION) %>%
    plot(show_data=TRUE)


## ----fitModel1h, results='markdown', eval=TRUE, cache=TRUE--------------------
starling_rstanarm2 <- stan_glmer(MASS ~ MONTH*SITUATION+(1|BIRD),
                                data = starling,
                                family = gaussian(),
                                prior_intercept = normal(84, 17, autoscale = FALSE),
                                prior = normal(0, c(33, 39, 39, 39, 50, 50, 50), autoscale = FALSE),
                                prior_aux=rstanarm::exponential(0.15, autoscale = FALSE),
                                prior_covariance = decov(1, 1, 1, 1),
                                prior_PD = TRUE,
                                iter = 5000,
                                warmup = 1000,
                                chains = 3,
                                thin = 5,
                                refresh = 0
                                )


## ----fitModel1i, results='markdown', eval=TRUE, cache=FALSE-------------------
starling_rstanarm2 %>%
    ggpredict(~SITUATION*MONTH) %>%
    plot(show_data = TRUE)


## ----fitModel1j, results='markdown', eval=TRUE, cache=TRUE, dependson='fitModel1h'----
starling_rstanarm3 <- update(starling_rstanarm2,  prior_PD=FALSE)


## ----modelFit1k, results='markdown', eval=TRUE, fig.width=8, fig.height=8-----
posterior_vs_prior(starling_rstanarm3, color_by='vs', group_by=TRUE,
                   facet_args=list(scales='free_y'))


## ----modelFit1l, results='markdown', eval=TRUE, fig.width=6, fig.height=4-----
ggemmeans(starling_rstanarm3,  ~SITUATION*MONTH) %>% plot(show_data=TRUE)
ggpredict(starling_rstanarm3,  ~SITUATION*MONTH) %>% plot(show_data=TRUE)


## ----fitModel2a, results='markdown', eval=TRUE, cache=TRUE, paged.print=FALSE, tidy.opts = list(width.cutoff = 80)----
starling_form <- bf(MASS ~ MONTH*SITUATION+(MONTH + MONTH:SITUATION|BIRD),
  family = gaussian()
)
options(width=150)
starling_form |> get_prior(data = starling)
options(width=80)


## ----fitModel2h, results='markdown', eval=TRUE, cache=FALSE-------------------
starling_form <- bf(MASS ~ MONTH*SITUATION+(1|BIRD),
                     family = gaussian()
                   )
get_prior(starling_form, data = starling)

starling |>
  group_by(SITUATION, MONTH) |>
  summarise(
    mean(MASS),
    median(MASS),
    sd(MASS),
    mad(MASS))

standist::visualize("normal(80, 2.5)", xlim=c(50,100))
standist::visualize("student_t(3, 0, 6)",
                    xlim=c(-10,25))


## ----fitModel2h1, results='markdown', eval=TRUE, cache=TRUE-------------------
priors <- prior(normal(80, 2.5), class = 'Intercept') +
    ## prior(normal(0, 13), class = 'b', coef = "MONTHJan") +
    ## prior(normal(0, 15), class = 'b', coef = "SITUATIONnestMbox") +
    ## prior(normal(0, 15), class = 'b', coef = "SITUATIONother") +
    ## prior(normal(0, 15), class = 'b', coef = "SITUATIONtree") +
    prior(normal(0, 10), class = 'b') +
    ## prior(gamma(6.5,1), class = 'sigma') +
    prior(student_t(3, 0, 6), class = 'sigma') +
    prior(student_t(3, 0, 6), class = 'sd')
starling_brm2 <- brm(starling_form,
                  data = starling,
                  prior = priors,
                  sample_prior = 'only',
                  iter = 5000,
                  warmup = 1000,
                  chains = 3, cores = 3,
                  thin = 5,
                  refresh = 0,
                  backend = "cmdstanr"
                  )


## ----partialPlot2h1a, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
starling_brm2 |>
    conditional_effects("SITUATION:MONTH") |>
    plot(points = TRUE)
starling_brm2 |>
    ggpredict(~SITUATION*MONTH) |>
    plot(show_data = TRUE)


## ----fitModel2h1b, results='markdown', eval=TRUE, cache=TRUE------------------
starling_brm3 <- update(starling_brm2,
                       sample_prior = 'yes',
                       control = list(adapt_delta = 0.99),
                       refresh = 0)
save(starling_brm3, file = '../ws/testing/starling_brm3')


## ----partialPlot2h1b, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
starling_brm3 |>
    conditional_effects("SITUATION:MONTH") |>
    plot(points = TRUE)
starling_brm3 |>
    ggpredict(~SITUATION*MONTH) |>
    plot(show_data = TRUE)


## ----posterior2h2, results='markdown', eval=TRUE------------------------------
starling_brm3 |> get_variables()
starling_brm3 |> hypothesis('MONTHJan=0') |> plot()
starling_brm3 |> hypothesis('SITUATIONnestMbox=0') |> plot()


## ----posterior2h2a, results='markdown', eval=TRUE, fig.width = 7, fig.height = 5----
starling_brm3 |> SUYR_prior_and_posterior()


## ----fitModel2h3, results='markdown', eval=TRUE, cache=TRUE-------------------
# I am going to narrow the priors so as to stabalise the model
priors <- prior(normal(80, 2), class = 'Intercept') +
    #prior(normal(0, 13), class = 'b', coef = "MONTHJan") +
    #prior(normal(0, 15), class = 'b', coef = "SITUATIONnestMbox") +
    #prior(normal(0, 15), class = 'b', coef = "SITUATIONother") +
    #prior(normal(0, 15), class = 'b', coef = "SITUATIONtree") +
    prior(normal(0, 10), class = 'b') +
    ## prior(gamma(6.5,1), class = 'sigma') +
    prior(student_t(3, 0, 3), class = 'sigma') +
    prior(student_t(3, 0, 3), class = 'sd') +
    prior(lkj_corr_cholesky(1), class = 'cor')

starling_form <- bf(MASS ~ MONTH*SITUATION+(MONTH|BIRD),
                     family = gaussian()
                   )
starling_brm3a <- brm(starling_form,
                  data = starling,
                  prior = priors,
                  sample_prior = 'only',
                  iter = 5000,
                  warmup = 2500,
                  chains = 3, cores =  3,
                  thin = 5,
                  refresh = 0,
                  control = list(adapt_delta = 0.99),
                  backend = 'cmdstanr'
                  )
starling_brm4 <- brm(starling_form,
                  data = starling,
                  prior = priors,
                  sample_prior = 'yes',
                  iter = 10000,  # note, these changes are in response to issues
                  warmup = 5000,
                  chains = 3, cores = 3,
                  thin = 10,     #some autocorrelation noted
                  refresh = 0,
                  control = list(adapt_delta = 0.99, max_treedepth = 20)
                  )
save(starling_brm4, file = '../ws/testing/starling_brm4')


## ----posterior2k, results='markdown', eval=TRUE-------------------------------
starling_brm4 |> get_variables()
starling_brm4 |> hypothesis('MONTHJan=0') |> plot()
starling_brm4 |> hypothesis('SITUATIONnestMbox=0') |> plot()


## ----posterior2k1, results='markdown', eval=TRUE, fig.width = 7, fig.height = 5----
starling_brm4 |> SUYR_prior_and_posterior()


## ----posterior2k2, results='markdown', eval=TRUE, fig.width=10, fig.height=4----
starling_brm4 |>
  posterior_samples() |>
  dplyr::select(-`lp__`) |>
  pivot_longer(everything(), names_to = 'key') |>
  filter(!str_detect(key, '^r')) |>
  mutate(Type = ifelse(str_detect(key, 'prior'), 'Prior', 'Posterior'),
         ## Class = ifelse(str_detect(key, 'Intercept'),  'Intercept',
         ##         ifelse(str_detect(key, 'b'),  'b', 'sigma')),
         Class = case_when(
               str_detect(key, '(^b|^prior).*Intercept$') ~ 'Intercept',
               str_detect(key, 'b_SITUATION.*|prior_b_SITUATION.*') &
               !str_detect(key, '.*:.*') ~ 'SITUATION',
               str_detect(key, 'b_MONTH.*|prior_b_MONTH.*') &
               !str_detect(key, '.*\\:.*') ~ 'MONTH',
               str_detect(key, '.*\\:.*|prior_b_.*\\:.*') ~ 'Interaction',
               str_detect(key, 'sd') ~ 'sd',
               str_detect(key, '^cor|prior_cor') ~ 'cor',
             str_detect(key, 'sigma') ~ 'sigma'
             ),
         Par = str_replace(key, 'b_', '')) |>
  ggplot(aes(x = Type,  y = value, color = Par)) +
  stat_pointinterval(position = position_dodge())+
  facet_wrap(~Class,  scales = 'free')


## ----fitModel2h3a, results='markdown', eval=TRUE, cache=TRUE------------------
(l.1 <- starling_brm3 |> rstan::loo())
(l.2 <- starling_brm4 |> rstan::loo())
loo_compare(l.1, l.2)


## ----modelValidation2a, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
available_mcmc()


## ----modelValidation2b, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
pars <- starling_brm4 |> get_variables()
pars <- pars |>
    str_extract('^b_.*|[sS]igma|^sd.*') |>
    ## str_extract('^b.Intercept|^b_SITUTATION.*|^b_MONTH.*|[sS]igma|^sd.*') |>
    na.omit()
pars
starling_brm4 |> mcmc_plot(type='trace', variable = pars)


## ----modelValidation2c, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
starling_brm4 |> mcmc_plot(type='acf_bar', variable = pars)


## ----modelValidation2d, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
starling_brm4 |> mcmc_plot(type='rhat_hist')


## ----modelValidation2e, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
starling_brm4 |> mcmc_plot(type='neff_hist')


## ----modelValidation2f, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
starling_brm4 |> mcmc_plot(type='combo', variable = pars)
starling_brm4 |> mcmc_plot(type='violin', variable = pars)


## ----modelValidation2g, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
starling_brm4 |> get_variables()
pars <- starling_brm4 |> get_variables()
pars <- str_extract(pars, '^b_.*|^sigma$|^sd.*') |> na.omit()

starling_brm4$fit |>
    stan_trace(pars = pars)


## ----modelValidation2h, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
starling_brm4$fit |>
    stan_ac(pars = pars)


## ----modelValidation2i, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
starling_brm4$fit |> stan_rhat()


## ----modelValidation2j, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
starling_brm4$fit |> stan_ess()


## ----modelValidation2k, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
starling_brm4$fit |>
    stan_dens(separate_chains = TRUE, pars = pars)


## ----modelValidation2l, results='markdown', eval=TRUE, fig.width=6, fig.height=7----
## starling_ggs <- starling_brm3 %>% ggs(burnin = FALSE, inc_warmup = FALSE)
## starling_ggs %>% ggs_traceplot()


## ----modelValidation2m, results='markdown', eval=TRUE, fig.width=6, fig.height=7----
## ggs_autocorrelation(starling_ggs)


## ----modelValidation2n, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
## ggs_Rhat(starling_ggs)


## ----modelValidation2o, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
## ggs_effective(starling_ggs)


## ----modelValidation2p, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
## ggs_crosscorrelation(starling_ggs)


## ----modelValidation2q, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
## ggs_grb(starling_ggs)


## ----modelValidation5a, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation5b, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
starling_brm4 |> pp_check(type = 'dens_overlay', nsamples = 100)


## ----modelValidation5c, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
## starling_brm4 |> pp_check(type = 'error_scatter_avg')


## ----modelValidation5e, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
starling_brm4 |> pp_check(group = 'BIRD', type = 'intervals')


## ----modelValidation5g, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(starling_brm2)


## ----modelValidation6aa, results='markdown', eval=TRUE, fig.width=8, fig.height=10----
starling_resids <- make_brms_dharma_res(starling_brm4, integerResponse = FALSE)
wrap_elements(~testUniformity(starling_resids)) +
               wrap_elements(~plotResiduals(starling_resids, form = factor(rep(1, nrow(starling))))) +
               wrap_elements(~plotResiduals(starling_resids, quantreg = TRUE)) +
               wrap_elements(~testDispersion(starling_resids))


## ----partialPlot2d, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
starling_brm4 |>
    conditional_effects("SITUATION:MONTH") |>
    plot(points = TRUE)


## ----partialPlot2a, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
starling_brm4 |>
    ggpredict(~SITUATION*MONTH) |>
    plot(show_data = TRUE)


## ----partialPlot2b, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
starling_brm4 |>
    ggemmeans(~SITUATION*MONTH) |>
    plot(show_data = TRUE)


## ----partialPlot2c, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
Partial.obs <- starling_brm4$data |>
    mutate(Pred = predict(starling_brm4)[,'Estimate'],
           Resid = resid(starling_brm4)[,'Estimate'],
           Obs = Pred + Resid)

starling_brm4 |>
    fitted_draws(newdata = starling, re_formula = NA) |>
    median_hdci() |>
    ggplot(aes(x = SITUATION, y = .value, color = MONTH)) +
    geom_pointrange(aes(ymin = .lower, ymax = .upper)) +
    geom_line() +
    geom_point(data = Partial.obs,  aes(y = Obs,  x = SITUATION, color = MONTH),
               position = position_nudge(x = 0.1)) +
    geom_point(data = starling,  aes(y = MASS,  x = SITUATION, color = MONTH), alpha=0.2,
               position = position_nudge(x = 0.05))


## ----summariseModel2a, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
starling_brm4 |> summary()


## ----summariseModel2a1, results='markdown', eval=TRUE, fig.width=8, fig.height=5, echo=FALSE----
starling_sum <- summary(starling_brm4)


## ----summariseModel2bm, results='markdown', eval=TRUE, fig.width=8, fig.height=5,echo=FALSE----
starling_brm4 |> as_draws_df()
starling_brm4 |>
  as_draws_df() |>
  summarise_draws(
    median,
    HDInterval::hdi,
    Pl = ~mean(.x < 0),
    Pg = ~mean(.x > 0),
    rhat,
    ess_bulk
  )


## ----summariseModel2b, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
starling_brm4$fit |>
    tidyMCMC(estimate.method = 'median',
             conf.int = TRUE,  conf.method = 'HPDinterval',
             rhat = TRUE, ess = TRUE)

## ----summariseModel2b1, results='markdown', eval=TRUE, fig.width=8, fig.height=5,echo=FALSE----
starling_tidy <- tidyMCMC(starling_brm4$fit, estimate.method='median',
                         conf.int=TRUE,  conf.method='HPDinterval',
                         rhat=TRUE, ess=TRUE)


## ----summariseModel2c, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
starling_brm4 |> get_variables()
starling_draw <- starling_brm4 |>
    gather_draws(`b.Intercept.*|b_SITUATION.*|b_MONTH.*`,  regex=TRUE)
starling_draw


## ----summariseModel2c1, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
starling_draw |> median_hdci()


## ----summariseModel2c3, results='markdown', eval=TRUE, fig.width=8, fig.height=5,echo=FALSE----
starling_gather <- starling_brm4 |>
    gather_draws(`b_Intercept.*|b_SITUATION.*|b_MONTH.*`,  regex=TRUE) |>
    median_hdci()


## ----summariseModel2c4, results='markdown', eval=TRUE, fig.width=8, fig.height=5,echo=TRUE----
starling_brm3 |>
    gather_draws(`b_Intercept.*|b_SITUATION.*|b_MONTH.*`, regex=TRUE) |>
    ggplot() +
    geom_vline(xintercept=0, linetype='dashed') +
    stat_slab(aes(x = .value, y = .variable,
                  fill = stat(ggdist::cut_cdf_qi(cdf,
                                                 .width = c(0.5, 0.8, 0.95),
                                                 labels = scales::percent_format())
                              )), color='black') +
    scale_fill_brewer('Interval', direction = -1, na.translate = FALSE)

starling_brm3 |>
    gather_draws(`.Intercept.*|b_SITUATION.*|b_MONTH.*`, regex=TRUE) |>
    ggplot() +
    geom_vline(xintercept = 0, linetype='dashed') +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    theme_classic()


## ----summariseModel2j, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
starling_brm3$fit %>% plot(type='intervals')


## ----summariseModel2ka, results='markdown', eval=TRUE, fig.width=8, fig.height=5,echo=TRUE----
starling_brm4 |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(.variable != 'b_Intercept') |>
    ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    facet_wrap(~.variable, scales='free')

starling_brm4 |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(.variable != 'b_Intercept') |>
    ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    geom_vline(xintercept = 0, linetype = 'dashed')


## ----summariseModel2c7, results='markdown', eval=TRUE, fig.width=8, fig.height=5,echo=TRUE----
starling_brm4 |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(.variable != 'b_Intercept') |>
    ggplot() +
    geom_density_ridges(aes(x=.value, y = .variable), alpha=0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed')
##Or in colour
starling_brm4 |>
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
starling_brm4 |>
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


## ----summariseModel2d, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
starling_brm4 |> tidy_draws()


## ----summariseModel2e, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
starling_brm4 |> spread_draws(`.*Intercept.*|b_SITUATION.*|b_MONTH.*`,  regex=TRUE)


## ----summariseModel2f, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
starling_brm4 |> posterior_samples() |> as_tibble()


## ----summariseModel2g, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
starling_brm4 |>
    bayes_R2(re.form = NA, summary=FALSE) |>
    median_hdci()
starling_brm4 |>
    bayes_R2(re.form = ~(1|BIRD), summary=FALSE) |>
    median_hdci()
starling_brm4 |>
    bayes_R2(re.form = ~(MONTH|BIRD), summary=FALSE) |>
    median_hdci()
starling_brm4 |>
    bayes_R2(re.form = ~(MONTH + MONTH:SITUATION|BIRD), summary=FALSE) |>
    median_hdci()
starling_brm4 |>
    bayes_R2(re.form = ~(MONTH*SITUATION|BIRD), summary=FALSE) |>
    median_hdci()


## -----------------------------------------------------------------------------
#| label: modelsummary
#| results: markup
#| eval: true
#| echo: true
#| cache: false
starling_brm4 |> modelsummary(
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
starling_brm4 |> modelplot(exponentiate = FALSE)


## ----postHoc1a, results='markdown', eval=TRUE, echo=1-------------------------
starling_brm4 |>
    emmeans(~SITUATION) |>
    pairs()

starling_brm4 |>
    emmeans(~MONTH) |>
    pairs()

starling_brm4 |>
    emmeans(~MONTH|SITUATION) |>
    pairs()

## To get percentage change
starling_brm4 |>
    emmeans(~MONTH|SITUATION) |>
    regrid(transform = 'log') |>
    pairs(reverse = TRUE) |>
    regrid() |>
    gather_emmeans_draws() |>
    mutate(.value = 100*(.value - 1)) |>
    median_hdci()

## OR if you want both absolute and percentage change
starling_em <- starling_brm4 |>
    emmeans(~MONTH|SITUATION) |>
    gather_emmeans_draws() |>
    spread(key=MONTH, value=.value) |>
    mutate(Eff=Jan-Nov,
           PEff=100*(Jan-Nov)/Nov)
starling_em |> head()

starling_em |>
    ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 10, linetype = "dashed", color='red') +
    stat_halfeye(aes(x = PEff, y = SITUATION)) +
    theme_classic()

starling_em |> median_hdci(PEff)

starling_em |>
    summarize(
        Prob=sum(PEff>0)/n(),
        `Prob>10`=sum(PEff>10)/n())


## ----postHoc2a, results='markdown', eval=TRUE, echo=1-------------------------
levels(starling$SITUATION)
cmat <- cbind(Comp1=c(0.5, 0.5, -1, 0),
              Comp2=c(1, -0.5, -0.5, 0))
cmat <- cbind("Nat vs Art"=c(-1/3, -1/3, -1/3, 1),
              "Tree vs I/NB" =  c(-1/2, -1/2, 0, 1),
              "Tree vs NB" =  c(0, -1, 0, 1))

starling_brm4 |>
    emmeans(~SITUATION) |>
    contrast(list(SITUATION=cmat))

starling_brm4 |>
  emmeans(~SITUATION) |>
  contrast(list(SITUATION=cmat)) |>
  gather_emmeans_draws()

starling_brm4 |>
    emmeans(~SITUATION|MONTH) |>
    contrast(list(SITUATION=cmat))

starling_em <- starling_brm4 |>
    emmeans(~SITUATION|MONTH) |>
    contrast(list(SITUATION=cmat)) |>
    gather_emmeans_draws() |>
    spread(key=MONTH, value=.value) |>
    mutate(Eff=Jan-Nov,
           PEff=100*(Jan-Nov)/Nov)
starling_em |> head()


## ----fitModel, results='markdown', eval=FALSE, echo=FALSE---------------------
# starling_rstan |> get_variables()
# plot(starling_rstan,  'mcmc_trace', regex_pars = '^.Intercept|^SITUATION|^MONTH|[sS]igma')
# plot(starling_rstan,  'mcmc_acf_bar', regex_pars = '^.Intercept|^SITUATION|^MONTH|[sS]igma')
# plot(starling_rstan,  'mcmc_rhat_hist', regex_pars = '^.Intercept|^SITUATION|^MONTH|[sS]igma')
# plot(starling_rstan,  'mcmc_neff_hist', regex_pars = '^.Intercept|^SITUATION|^MONTH|[sS]igma')
# 
# 
# preds <- posterior_predict(starling_rstan,  nsamples=250,  summary=FALSE)
# starling_resids <- createDHARMa(simulatedResponse = t(preds),
#                             observedResponse = starling$MASS,
#                             fittedPredictedResponse = apply(preds, 2, median))
# plot(starling_resids)
# 
# 
# starling_rstan1 = stan_glmer(MASS ~ MONTH*SITUATION+(MONTH|BIRD),data=starling,
#                             iter=5000, warmup=2000, thin=5, chains=3, refresh=0)
# starling_rstan1 = stan_glmer(MASS ~ MONTH*SITUATION+(MONTH|BIRD),data=starling,
#                              iter=5000, warmup=2000, thin=5, chains=3, refresh=0,
#                              adapt_delta = 0.99)
# #pairs(starling_rstan1,  pars=c('(Intercept)', 'MONTHNov'))
# starling_rstan1 %>% get_variables()
# pairs(starling_rstan1,  regex_pars=c('SITUATION', 'sigma'))
# prior_summary(starling_rstan1)
# 
# plot(starling_rstan1,  'mcmc_trace', regex_pars = '^.Intercept|^SITUATION|^MONTH|[sS]igma')
# plot(starling_rstan1,  'mcmc_acf_bar', regex_pars = '^.Intercept|^SITUATION|^MONTH|[sS]igma')
# plot(starling_rstan1,  'mcmc_rhat_hist', regex_pars = '^.Intercept|^SITUATION|^MONTH|[sS]igma')
# plot(starling_rstan1,  'mcmc_neff_hist', regex_pars = '^.Intercept|^SITUATION|^MONTH|[sS]igma')
# 
# starling_rstan1 = stan_glmer(MASS ~ MONTH*SITUATION+(MONTH|BIRD),data=starling,
#                              iter=10000, warmup=5000, thin=15, chains=3, refresh=0,
#                              adapt_delta = 0.99)
# 
# plot(starling_rstan1,  'mcmc_trace', regex_pars = '^.Intercept|^SITUATION|^MONTH|[sS]igma')
# plot(starling_rstan1,  'mcmc_acf_bar', regex_pars = '^.Intercept|^SITUATION|^MONTH|[sS]igma')
# plot(starling_rstan1,  'mcmc_rhat_hist', regex_pars = '^.Intercept|^SITUATION|^MONTH|[sS]igma')
# plot(starling_rstan1,  'mcmc_neff_hist', regex_pars = '^.Intercept|^SITUATION|^MONTH|[sS]igma')
# preds <- posterior_predict(starling_rstan1,  nsamples=250,  summary=FALSE)
# starling_resids <- createDHARMa(simulatedResponse = t(preds),
#                             observedResponse = starling$MASS,
#                             fittedPredictedResponse = apply(preds, 2, median))
# plot(starling_resids)
# 
# (l.1 <- loo(starling_rstan))
# (l.2 <- loo(starling_rstan1))
# loo_compare(l.1, l.2)
# 
# as.matrix(starling_rstan) %>% colnames
# posterior_vs_prior(starling_rstan1, color_by='vs', group_by=TRUE, regex_pars=c('^MONTH','^SITUATION','^[sS]igma'),
#                    facet_args=list(scales='free_y'))
# 
# 
# g=ggpredict(starling_rstan1) %>% plot
# do.call('grid.arrange',  g)
# ggemmeans(starling_rstan1, ~SITUATION|MONTH) %>% plot
# 
# summary(starling_rstan1)
# 
# nms <- starling_rstan1 %>% get_variables()
# nms
# wch <- grep('^.Intercept|^MONTH|^SITUATION|[sS]igma', nms)
# tidyMCMC(starling_rstan1$stanfit,conf.int=TRUE, conf.method='HPDinterval',
#          rhat=TRUE, ess=TRUE, pars=nms[wch])
# 
# emmeans(starling_rstan1, pairwise~MONTH|SITUATION)
# starling_em = emmeans(starling_rstan1, ~MONTH|SITUATION) %>%
#     gather_emmeans_draws() %>% spread(key=MONTH, value=.value) %>%
#     mutate(Eff=Jan-Nov,
#            PEff=100*(Jan-Nov)/Nov)
# starling_em %>% head
# 
# starling_em %>% ungroup %>%
#     dplyr::select(SITUATION,Eff,PEff) %>% group_by(SITUATION) %>% median_hdi
# 
# starling_em %>% ungroup %>%
#     dplyr::select(SITUATION,Eff,PEff) %>% group_by(SITUATION) %>%
#     summarize(Prob=sum(PEff>10)/n())
# 
# bayes_R2(starling_rstan1, re.form=NA) %>% median_hdi
# bayes_R2(starling_rstan1, re.form=~(1|BIRD)) %>% median_hdi
# bayes_R2(starling_rstan1, re.form=~(MONTH|BIRD)) %>% median_hdi
# 
# newdata = emmeans(starling_rstan1, ~MONTH|SITUATION) %>% as.data.frame
# head(newdata)
# ggplot(newdata, aes(y=emmean, x=SITUATION)) +
#     geom_pointrange(aes(ymin=lower.HPD, ymax=upper.HPD, fill=MONTH),
#                     position=position_dodge(width=0.3), shape=21)

