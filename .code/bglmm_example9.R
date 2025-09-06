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
library(effects)   #for partial effects plots
library(emmeans)   #for estimating marginal means
library(MASS)      #for glm.nb
library(tidyverse) #for data wrangling
library(brms)
library(tidybayes)
library(bayesplot)
library(broom.mixed)
library(rstan)
library(patchwork)
library(modelsummary)
library(DHARMa)
source('helperFunctions.R')


## ----readDataP, results='markdown', eval=FALSE--------------------------------
# hughes <- read_csv('../data/hughes_full.csv', trim_ws=TRUE)
# glimpse(hughes)


## ----processDataP, results='markdown', eval=FALSE-----------------------------
# hughes <- hughes |>
#     mutate(fYear=factor(Year),
#            Score=ifelse(Score==5,4,Score),
#            oScore = factor(Score, ordered=TRUE),
#            nScore = as.numeric(factor(Score, ordered=TRUE)),
#            SectorThree=factor(SectorThree, levels=c('North','Central','South')),
#            fReef=factor(ReefID),
#            nReef=as.numeric(fReef))
# # now make a version that is just 2016
# hughes <- hughes |> filter(fYear==2016) |>
#     dplyr::select(REEF=ReefID, HABITAT=Habitat, SECTOR=SectorThree, SCORE=Score)
# write_csv(hughes, file='../data/hughes.csv')
# ## hughes.colors = c('#FFFFFF', rev(heat.colors(length(levels(hughes$oScore))))[-1])


## ----readData, results='markdown', eval=TRUE----------------------------------
hughes = read_csv('../data/hughes.csv', trim_ws=TRUE)
glimpse(hughes)


## -----------------------------------------------------------------------------
#| label: examinData
glimpse(hughes)


## -----------------------------------------------------------------------------
## Explore the first 6 rows of the data
head(hughes)


## -----------------------------------------------------------------------------
str(hughes)


## -----------------------------------------------------------------------------
hughes |> datawizard::data_codebook()


## -----------------------------------------------------------------------------
hughes |> modelsummary::datasummary_skim(by = "HABITAT")


## ----processData, results='markdown', eval=TRUE, mhidden=TRUE-----------------
hughes <- hughes |>
    mutate(oSCORE = factor(SCORE, ordered = TRUE),
           HABITAT = factor(HABITAT),
           SECTOR = factor(SECTOR, levels = c("North", "Central", "South")),
           REEF = factor(REEF))
## hughes.colors = c("#FFFFFF", rev(heat.colors(length(levels(hughes$oSCORE))))[-1])


## ----EDA2, results='markdown', eval=TRUE, fig.width=10, fig.height=10, mhidden=TRUE----
# Scatterplot
hughes |>
    ggplot(aes(y = oSCORE, x = HABITAT)) +
    geom_point(position = position_jitter()) +
    facet_wrap(~SECTOR)


## ----EDA2a, results='markdown', eval=TRUE, fig.width=10, fig.height=10, mhidden=TRUE----
hughes |>
    group_by(SECTOR, REEF, HABITAT) |>
    summarise(SCORE = mean(SCORE)) |>
    ungroup() |>
    ggplot(aes(y = SCORE, x = as.numeric(HABITAT), group = REEF)) +
    geom_blank(aes(x = HABITAT)) +
    geom_line() +
    facet_grid(~SECTOR)


## ----EDA2b, results='markdown', eval=TRUE, fig.width=10, fig.height=10, mhidden=TRUE----
hughes |>
    group_by(SECTOR, HABITAT, oSCORE) |>
    summarise(n = n()) |>
    ungroup() |>
    group_by(SECTOR, HABITAT) |>
    mutate(prop = n/sum(n)) |>
    mutate(oSCORE = factor(oSCORE, levels = rev(levels(oSCORE)))) ->
    hughes.sum

## hughes.sum <- hughes |>
##     count(SECTOR,HABITAT,oSCORE) |>
##     group_by(SECTOR, HABITAT) |>
##     mutate(prop=prop.table(n),
##            oSCORE=factor(oSCORE, levels=rev(levels(oSCORE))))

hughes.sum |> head()

ggplot(data=hughes.sum, aes(y=prop, x=HABITAT)) +
    geom_bar(stat='Identity', aes(fill=oSCORE), color='black') +
    facet_grid(~SECTOR) +
    ## scale_fill_manual('Bleaching score', values=rev(hughes.colors) ) +
    scale_fill_manual('Bleaching score', values=c(heat.colors(5)[-5], '#FFFFFF') ) +
    scale_y_continuous('Proportion of Reef', expand=c(0,0))+
    theme_bw() +
    theme(panel.spacing.y=unit(10,'pt'))


## ----fitModel2a, results='markdown', eval=TRUE, cache=TRUE, paged.print=FALSE, tidy.opts = list(width.cutoff = 80)----
hughes.form <- bf(oSCORE ~ HABITAT*SECTOR + (1|REEF),
  family = cumulative(link = "logit", threshold = "flexible")
)
options(width=150)
hughes.form %>% get_prior(data = hughes)
options(width=80)


## ----fitModel2h, results='markdown', eval=TRUE, cache=FALSE-------------------


## ----fitModel2h1, results='markdown', eval=TRUE, cache=TRUE-------------------
priors <- prior(normal(0, 1), class = 'Intercept') +
    prior(normal(0, 1), class = 'b') +
    prior(student_t(3, 0, 1), class = 'sd')
hughes.form <- bf(oSCORE ~ HABITAT*SECTOR + (1|REEF),
  family = cumulative(link = "logit", threshold = "flexible")
)
hughes.brm2 <- brm(hughes.form,
                  data = hughes,
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
## hughes.brm2 %>%
##     ggpredict(~HABITAT*SECTOR) |>
##     plot(show_data = TRUE)

hughes.brm2 |>
  conditional_effects("HABITAT",
    conditions = make_conditions(hughes.brm2, "SECTOR"),
    categorical = TRUE
)


## ----fitModel2h1b, results='markdown', eval=TRUE, cache=TRUE------------------
hughes.brm3 <- update(hughes.brm2,
                       sample_prior = 'yes',
                       control = list(adapt_delta = 0.99),
                       refresh = 100)
## save(hughes.brm3, file = '../ws/testing/hughes.brm3')
## load(file = '../ws/testing/hughes.brm3')


## ----partialPlot2h1b, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
hughes.brm3 |>
    conditional_effects("HABITAT", categorical = TRUE) |>
    plot(points = TRUE)
hughes.brm3 |>
    conditional_effects("SECTOR", categorical = TRUE) |>
    plot(points = TRUE)

hughes.brm3 |>
  conditional_effects("HABITAT",
    conditions = make_conditions(hughes.brm3, "SECTOR"),
    categorical = TRUE
)

a <-
hughes.brm3 |>
  conditional_effects("HABITAT",
    conditions = make_conditions(hughes.brm3, "SECTOR"),
    categorical = TRUE,
)


## ----posterior2h2, results='markdown', eval=TRUE------------------------------
hughes.brm3 %>% get_variables()
hughes.brm3 %>% hypothesis('HABITATF=0') %>% plot
hughes.brm3 %>% hypothesis('SECTORCentral=0') %>% plot


## ----posterior2h2a, results='markdown', eval=TRUE, fig.width = 7, fig.height = 5----
# hughes.brm3 %>% SUYR_prior_and_posterior()


## ----modelValidation2a, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
available_mcmc()


## ----modelValidation2b, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
pars <- hughes.brm3 %>% get_variables()
pars <- hughes.brm3 |>
    get_variables() |>
    str_subset("^b_.*|[sS]igma|^sd.*")
pars
hughes.brm3 |> mcmc_plot(type='trace', variable = pars)


## ----modelValidation2c, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
hughes.brm3 |> mcmc_plot(type='acf_bar', variable = pars)


## ----modelValidation2d, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
hughes.brm3 |> mcmc_plot(type='rhat_hist')


## ----modelValidation2e, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
hughes.brm2 |> mcmc_plot(type='neff_hist')


## ----modelValidation2f, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
hughes.brm3 |> mcmc_plot(type='combo', variable = pars)
hughes.brm3 |> mcmc_plot(type='violin', variable = pars)


## ----modelValidation2g, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
hughes.brm3 |> get_variables()
pars <- hughes.brm3 |>
    get_variables() |>
    str_subset("^b_.*|[sS]igma|^sd.*")
pars
hughes.brm3$fit |>
    stan_trace(pars = pars)


## ----modelValidation2h, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
hughes.brm3$fit |>
    stan_ac(pars = pars)


## ----modelValidation2i, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
hughes.brm3$fit |> stan_rhat()


## ----modelValidation2j, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
hughes.brm3$fit |> stan_ess()


## ----modelValidation2k, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
hughes.brm3$fit |>
    stan_dens(separate_chains = TRUE, pars = pars)


## ----modelValidation5a, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation5b, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
hughes.brm3 |> pp_check(type = 'dens_overlay', ndraws = 100)


## ----modelValidation5c, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
## hughes.brm3 %>% pp_check(type = 'error_scatter_avg')


## ----modelValidation5e, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
hughes.brm3 %>% pp_check(group = 'REEF', type = 'intervals')


## ----modelValidation5g, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(hughes.brm2)


## ----modelValidation6aa, results='markdown', eval=TRUE, fig.width=8, fig.height=10----
hughes.resids <- make_brms_dharma_res(hughes.brm3, integerResponse = FALSE)
wrap_elements(~testUniformity(hughes.resids)) +
               wrap_elements(~plotResiduals(hughes.resids, form = factor(rep(1, nrow(hughes))))) +
               wrap_elements(~plotResiduals(hughes.resids, quantreg = TRUE)) +
               wrap_elements(~plotResiduals(hughes.resids, quantreg = TRUE, asFactor = FALSE)) +
               wrap_elements(~testDispersion(hughes.resids))


## ----summariseModel2a, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
hughes.brm3 |>
  conditional_effects("HABITAT",
    conditions = make_conditions(hughes.brm3, "SECTOR"),
    categorical = TRUE
)

hughes.brm3 %>% summary()


## ----summariseModel2a1, results='markdown', eval=TRUE, fig.width=8, fig.height=5, echo=FALSE----
hughes.sum <- summary(hughes.brm3)


## ----summariseModel2bm, results='markdown', eval=TRUE, fig.width=8, fig.height=5,echo=TRUE----
hughes.brm3 |> as_draws_df()
hughes.brm3 |>
    as_draws_df() |>
    exp() |>
    ## dplyr::select(matches("^b_Intercept.*|^b_HABITAT.*|^b_SECTOR.*")) |>
    dplyr::select(matches("^b_*")) |>
  summarise_draws(
    median,
    HDInterval::hdi,
    Pl = ~mean(.x < 1),
    Pg = ~mean(.x > 1),
    rhat,
    ess_bulk
  )


## ----summariseModel2g, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
hughes.brm3 |>
    bayes_R2(re.form = NA, summary=FALSE) |>
    median_hdci()
hughes.brm3 |>
    bayes_R2(re.form = ~(1|REEF), summary=FALSE) |>
    median_hdci()


## ----postHoc1a, results='markdown', eval=TRUE, echo=TRUE----------------------
newdata <- with(hughes, expand.grid(
    HABITAT = levels(HABITAT),
    SECTOR = levels(SECTOR)
))
newdata
predict(hughes.brm3, newdata = newdata, re_formula = NA)

sum(0:4 * (c(0.0000000000, 0.0041666667, 0.04666667, 0.35000000, 0.5991666667)))

add_epred_draws(hughes.brm3,
    newdata = newdata,
    re_formula = NA
) |>
    ## filter(.draw == 1)
  mutate(fit = as.numeric(as.character(.category)) * .epred) |>
  group_by(HABITAT, SECTOR, .draw) |>
  summarise(fit = sum(fit)) |>
  summarise_draws(
    median,
    HDInterval::hdi
  ) |>
  arrange(SECTOR, HABITAT)

hughes.cellmeans <-
    add_epred_draws(hughes.brm3, newdata = newdata, re_formula = NA) |>
    mutate(fit = as.numeric(as.character(.category)) * .epred) |>
    group_by(HABITAT, SECTOR, .draw) |>
    summarise(fit = sum(fit)) |>
    summarise_draws(
        median,
        HDInterval::hdi
    ) |>
  arrange(SECTOR, HABITAT)

fig1 <-
hughes.cellmeans |>
  ggplot(aes(y = median, x = HABITAT)) +
    geom_hline(yintercept=1, linetype='dashed', size=0.1) +
    geom_hline(yintercept=2, linetype='dashed', size=0.1) +
    geom_hline(yintercept=3, linetype='dashed', size=0.1) +
    geom_pointrange(aes(ymin = lower, ymax = upper)) +
    facet_grid(~SECTOR) +
    scale_y_continuous('Bleaching score', breaks=(0:4), labels=0:4, limits=c(0,4),expand=c(0,0)) +
    theme_bw() +
    theme(panel.spacing.y=unit(10,'pt'))



## Alternative
## marginaleffects::avg_predictions(hughes.brm3, newdata = newdata, c("SECTOR", "HABITAT"), re_formula = NA) |>
##     marginaleffects::posterior_draws() |>
##   dplyr::select(SECTOR, HABITAT, group, estimate) |>
##     group_by(SECTOR, HABITAT, group) |>
##     summarise_draws(median)

## marginaleffects::avg_predictions(hughes.brm3, newdata = newdata, c("SECTOR", "HABITAT"), re_formula = NA) |>
##     marginaleffects::posterior_draws() |>
##     ggplot(aes(x = estimate, y = HABITAT, fill = group)) +
##   stat_halfeye() +
##   facet_grid(~SECTOR)
##     head()


## ----postHoc1b, results='markdown', eval=TRUE, echo=TRUE, width = 10, height = 3----
newdata <- with(hughes, expand.grid(
    HABITAT = levels(HABITAT),
    SECTOR = levels(SECTOR)
))

hughes.eff <-
add_epred_draws(hughes.brm3,
  newdata = newdata,
  re_formula = NA
) |>
  mutate(fit = as.numeric(as.character(.category)) * .epred) |>
  group_by(HABITAT, SECTOR, .draw) |>
  summarise(fit = log(sum(fit))) |>
  ## compare_levels(var = fit, by = HABITAT, comparison = "pairwise") |>
  tidybayes::compare_levels(
    var = fit, by = HABITAT,
    ## comparison = emmeans_comparison("tukey", reverse = TRUE)) |>
    comparison = emmeans_comparison("revpairwise")
    ## comparison = emmeans_comparison("pairwise")
    ## comparison = "pairwise"
  ) |>
  mutate(fit = exp(fit)) |>
  group_by(SECTOR, HABITAT) |>
  summarise_draws(
    median,
    HDInterval::hdi
  )

fig2 <-
hughes.eff |>
  ggplot(aes(x = median, y = HABITAT)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  geom_pointrange(aes(xmin = lower, xmax = upper)) +
  ## facet_grid(~SECTOR, scales = "free") +
  facet_grid(~SECTOR) +
  scale_x_continuous("Effect size (percentage change in bleaching category)",
    trans = scales::log2_trans(),
    breaks = seq(0.25, 2.5, by = 0.25),
    ## breaks = scales::trans_breaks("log10", function(x) exp(x)),
    ## breaks = scales::trans_breaks("log2", function(x) 2^x),
    ## labels = scales::trans_format("log2", function(x) x),
    ## labels = function(x) x
    labels = function(x) (x - 1) * 100
  )+
  theme_bw()

fig1 + fig2

## Alternative with slab

hughes.eff <-
add_epred_draws(hughes.brm3,
  newdata = newdata,
  re_formula = NA
) |>
  mutate(fit = as.numeric(as.character(.category)) * .epred) |>
  group_by(HABITAT, SECTOR, .draw) |>
  summarise(fit = log(sum(fit))) |>
  ## compare_levels(var = fit, by = HABITAT, comparison = "pairwise") |>
  tidybayes::compare_levels(
    var = fit, by = HABITAT,
    ## comparison = emmeans_comparison("tukey", reverse = TRUE)) |>
    comparison = emmeans_comparison("revpairwise")
    ## comparison = emmeans_comparison("pairwise")
    ## comparison = "pairwise"
  ) |>
  mutate(fit = exp(fit))

fig2 <-
    hughes.eff |>
    ggplot(aes(x = fit, y = HABITAT)) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    stat_halfeye(slab_fill = "orange", slab_alpha = 0.5, normalize = "panels") +
  facet_grid(~SECTOR, scales = "free") +
  ## facet_grid(~SECTOR) +
  scale_x_continuous("Effect size (percentage change in bleaching category)",
    trans = scales::log2_trans(),
    breaks = seq(0.25, 2.5, by = 0.25),
    ## breaks = scales::trans_breaks("log10", function(x) exp(x)),
    ## breaks = scales::trans_breaks("log2", function(x) 2^x),
    ## labels = scales::trans_format("log2", function(x) x),
    ## labels = function(x) x
    labels = function(x) (x - 1) * 100
  )+
  theme_bw()

fig1 + fig2


## -----------------------------------------------------------------------------
#| label: temp
#| results: markup
#| eval: false
#| echo: false
#| cache: false
# summary(hughes.brm3)
# brms::inv_logit_scaled(-0.33)
# brms::inv_logit_scaled(-1.34)
# 
# newdata <- with(hughes, expand.grid(
#     HABITAT = levels(HABITAT),
#     SECTOR = levels(SECTOR)
# ))
# ## Cellmeans
# add_epred_draws(hughes.brm3, newdata = newdata, re_formula = NA) |>
#   mutate(fit = as.numeric(as.character(.category)) * .epred) |>
#   group_by(HABITAT, SECTOR, .draw) |>
#   summarise(fit = sum(fit)) |>
#   summarise_draws(
#     median,
#     HDInterval::hdi
#   ) |>
#   arrange(SECTOR, HABITAT) |>
#   ggplot(aes(y = median, x = HABITAT)) +
#     geom_hline(yintercept=1, linetype='dashed', size=0.1) +
#     geom_hline(yintercept=2, linetype='dashed', size=0.1) +
#     geom_hline(yintercept=3, linetype='dashed', size=0.1) +
#     geom_pointrange(aes(ymin = lower, ymax = upper)) +
#     facet_grid(~SECTOR) +
#     scale_y_continuous('Bleaching score', breaks=(0:4), labels=0:4, limits=c(0,4),expand=c(0,0)) +
#     theme_bw() +
#     theme(panel.spacing.y=unit(10,'pt'))
# 
#     ## Effects - habitats (absolute)
# add_epred_draws(hughes.brm3, newdata = newdata, re_formula = NA) |>
#   mutate(fit = as.numeric(as.character(.category)) * .epred) |>
#   group_by(HABITAT, SECTOR, .draw) |>
#   summarise(fit = sum(fit)) |>
#   ## compare_levels(var = fit, by = HABITAT, comparison = "pairwise") |>
#   tidybayes::compare_levels(var = fit, by = HABITAT,
#     ## comparison = emmeans_comparison("tukey", reverse = TRUE)) |>
#     comparison = emmeans_comparison("revpairwise")) |>
#   group_by(SECTOR, HABITAT) |>
#   summarise_draws(
#     median,
#     HDInterval::hdi
#   ) |>
#   ggplot(aes(x = median, y = HABITAT)) +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   geom_pointrange(aes(xmin = lower, xmax = upper)) +
#   facet_grid(~SECTOR) +
#   scale_x_continuous('Effect size (absolute change in bleaching category)')+
#   theme_bw()
# 
# ## Effects - habitats (fractional change)
# add_epred_draws(hughes.brm3, newdata = newdata, re_formula = NA) |>
#     mutate(fit = as.numeric(as.character(.category)) * .epred) |>
#     group_by(HABITAT, SECTOR, .draw) |>
#     summarise(fit = log(sum(fit))) |>
#     ## compare_levels(var = fit, by = HABITAT, comparison = "pairwise") |>
#     tidybayes::compare_levels(
#         var = fit, by = HABITAT,
#         ## comparison = emmeans_comparison("tukey", reverse = TRUE)) |>
#         comparison = emmeans_comparison("revpairwise")
#     ) |>
#     mutate(fit = exp(fit)) |>
#     group_by(SECTOR, HABITAT) |>
#     summarise_draws(
#         median,
#         HDInterval::hdi
#     ) |>
#     ggplot(aes(x = median, y = HABITAT)) +
#     geom_vline(xintercept = 1, linetype = "dashed") +
#     geom_pointrange(aes(xmin = lower, xmax = upper)) +
#     facet_grid(~SECTOR, scales = "free") +
#     scale_x_continuous("Effect size (percentage change in bleaching category)",
#         trans = scales::log2_trans(),
#         breaks = scales::trans_breaks("log2", function(x) 2^x),
#         ## labels = scales::trans_format("log2", function(x) x)
#         ## labels = function(x) x
#         labels = function(x) (x - 1) * 100
#   )+
#   theme_bw()
# 
# 
# 
# ## More manual
# add_epred_draws(hughes.brm3, newdata = newdata, re_formula = NA) |>
#     mutate(fit = as.numeric(as.character(.category)) * .epred) |>
#     group_by(HABITAT, SECTOR, .draw) |>
#     summarise(fit = sum(fit)) |>
#     group_by(SECTOR, .draw) |>
#   reframe(
#     contrast = colnames(emmeans:::pairwise.emmc(HABITAT)),
#     fit1 = t(as.vector(fit) %*% as.matrix(emmeans:::pairwise.emmc(HABITAT)))
#     ) |>
#     group_by(SECTOR, contrast) |>
#   summarise_draws(
#     median,
#     HDInterval::hdi
#     )
# 
# 
# sum(0:4 * posterior_epred(hughes.brm3, newdata = newdata, re_formula = NA)[,4,] |> colMeans())
# 
# 
# 
# 
# 
# 
# 
# 
# 
# hughes.brm3 |> emmeans(~HABITAT | SECTOR, mode = "latent")
# emmeans(hughes.clmm, ~oSCORE|HABITAT|SECTOR, mode='prob')
# emmeans(hughes.clmm, ~HABITAT|SECTOR)
# 
# newdata <- with(hughes, expand.grid(
#     HABITAT = levels(HABITAT),
#     SECTOR = levels(SECTOR),
#     REEF = NA
# ) )
# hughes.brm3 |>
#   epred_draws(newdata) |>
#   filter(HABITAT == "C", SECTOR == "North") |>
#   dplyr::select(HABITAT, SECTOR, .draw, .epred) |>
#     ungroup() |>
#   summarise(mean(.epred))
# hughes.brm3 |>
#     linpred_draws(newdata) |>
#     filter(HABITAT == "C", SECTOR == "North") |>
#     dplyr::select(HABITAT, SECTOR, .draw, .linpred) |>
#     ungroup() |>
#   summarise(mean(.linpred))
# 
# hughes.brm3 |>
#   predicted_draws(newdata) |>
#   filter(HABITAT == "C", SECTOR == "North") |>
#   dplyr::select(HABITAT, SECTOR, .draw, .prediction)
# 
# newdata <- cbind(newdata, predict(hughes.brm3, newdata = newdata))
# newdata
# apply(sweep(newdata[, 4:8], 2, 0:4, "*"), 1, sum)
# rowSums(1:5 * newdata[, 4:8])
# 1:5 * newdata[, 4:8]
# 
# 
# coefs <- hughes.brm3 |>
#     as_draws_df() |>
#     dplyr::select(matches("^b_HABITAT.*|^b_SECTOR.*")) |>
#     as.matrix()
# threshs <- hughes.brm3 |>
#     as_draws_df() |>
#     dplyr::select(matches("^b_Intercept.*")) |>
#     as.matrix()
# newdata <- with(hughes, expand.grid(
#     HABITAT = levels(HABITAT),
#     SECTOR = levels(SECTOR),
#     REEF = NA
# ) )
# Xmat = model.matrix(~HABITAT*SECTOR, data=newdata)[,-1]
# fit = coefs %*% t(Xmat)
# fit=sapply(1:4, function(i) threshs[,i] - fit, simplify='array')
# fit = binomial()$linkinv(fit)
# fit = aperm(fit, c(1,3,2))
# fit = sapply(1:dim(fit)[3], function(i) cbind(fit[,1,i], fit[,-1,i]-fit[,-length(fit[1,,i]),i]), simplify='array')
# out1 = sapply(1:dim(fit)[3], function(i) cbind(fit[,,i], 1-rowSums(fit[,,i])), simplify='array')
# 
# out2=sweep(out1, 2, 1:5, '*')
# out3 = apply(out2, 3, rowSums)
# 
# out3 = out3-1
# library(broom)
# newdata = newdata %>% cbind(tidyMCMC(as.mcmc(out3), conf.int=TRUE, conf.method='HPDinterval'))
# ScoreBoundaries = data.frame(Score=factor(0:4), ymin=c(0:4), ymax=c(1:5))
# ggplot(newdata) +
#     geom_blank(aes(y=estimate, x=HABITAT)) +
#     geom_hline(yintercept=1, linetype='dashed', size=0.1) +
#     geom_hline(yintercept=2, linetype='dashed', size=0.1) +
#     geom_hline(yintercept=3, linetype='dashed', size=0.1) +
#     #geom_rect(data=ScoreBoundaries, aes(ymin=ymin, ymax=ymax, xmin=-Inf, xmax=Inf, fill=Score), alpha=0.2) +
#     geom_pointrange(aes(y=estimate, x=HABITAT, ymin=conf.low, ymax=conf.high)) +
#     facet_grid(~SECTOR) +
#     scale_y_continuous('Bleaching score', breaks=(0:4), labels=0:4, limits=c(0,4),expand=c(0,0)) +
#     theme_bw() +
#     theme(panel.spacing.y=unit(10,'pt'))
# ###################
# n = rep(1, length(levels(hughes$HABITAT)))
# names(n) <- levels(hughes$HABITAT)
# tuk.cont = multcomp::contrMat(n,'Tukey')
# newdata=with(hughes, expand.grid(HABITAT=levels(HABITAT),
#                                     SECTOR=levels(SECTOR)))
# Xmat = model.matrix(~-1+HABITAT*SECTOR, data=newdata)
# Xmat = diag(ncol(Xmat))
# newdata = newdata %>% cbind(Xmat)
# Xmat = newdata %>% group_by(SECTOR) %>%
#     do({
#         x=.
#         xx = x[,-1:-3] %>% as.matrix
#         data.frame(Contrasts = rownames(tuk.cont), tuk.cont %*% xx)
#     })
# xs = Xmat[,1:3]
# Xmat = Xmat[,-1:-3] %>% as.matrix
# out4=out3 %*% t(Xmat)
# newdata = tidyMCMC(as.mcmc(out4), conf.int=TRUE, conf.method='HPDinterval') %>% bind_cols(xs)
# ggplot(newdata) +
#     geom_hline(yintercept=0) +
#     geom_pointrange(aes(y=estimate, x=Contrasts, ymin=conf.low, ymax=conf.high, color=SECTOR),
#                     position=position_dodge(width=0.5)) +
#     facet_grid(SECTOR~fYear) +
#     coord_flip() +
#     scale_y_continuous('Effect size')+
#     theme_bw()
# 
# 
# 
# hughes.clmm=ordinal::clmm(oSCORE ~ HABITAT*SECTOR+(1|REEF), data=hughes)
# hughes.clmm1=ordinal::clmm(oSCORE ~ HABITAT*SECTOR+(HABITAT|REEF), data=hughes)
# 
# hughes.clmm1 %>% ggemmeans(~HABITAT|SECTOR) %>% plot
# 
# 
# summary(hughes.clmm)
# summary(hughes.clmm1)
# 
# emmeans(hughes.clmm, ~oSCORE|HABITAT|SECTOR, mode='prob')
# ## emmeans(hughes.clmm1, ~oSCORE|HABITAT+SECTOR, mode='prob')
# emmeans(hughes.clmm, ~HABITAT|SECTOR, mode='mean.class')
# emmeans(hughes.clmm1, ~HABITAT|SECTOR, mode='mean.class')
# emmeans(hughes.clmm, ~HABITAT|SECTOR, mode='mean.class') %>% pairs()
# ## emmeans(hughes.clmm1, ~HABITAT|SECTOR, mode='mean.class') %>% pairs()

