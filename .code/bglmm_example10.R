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

library(broom)     #for tidy output
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
# ltmp <- read_csv(file='../data/ltmp_full.csv')
# glimpse(ltmp)


## ----processDataP, results='markdown', eval=FALSE-----------------------------
# coral_fams <- ltmp |>
#     filter(HC>0) |>
#     pull(FAMILY_2021) |>
#     unique()
# ltmp <- ltmp |>
#   filter(FAMILY_2021 %in% coral_fams) |>
#   dplyr::rename(P_CODE = P_CODE.y) |>
#   distinct() |>
#   group_by(P_CODE, AIMS_REEF_NAME, REPORT_YEAR, SITE_NO, TRANSECT_NO) |>
#   summarise(
#     HC = sum(HC, na.rm = TRUE),
#     n.points = sum(n.points),
#     total.points = unique(total.points)
#   ) |>
#   ungroup() |>
#   dplyr::select(AIMS_REEF_NAME, REPORT_YEAR, SITE_NO, TRANSECT_NO, HC, n.points, total.points) |>
#   filter(!is.na(AIMS_REEF_NAME)) |>
#   droplevels()
# write_csv(ltmp, file='../data/ltmp.csv')
# rm("ltmp")
# gc()


## ----readData, results='markdown', eval=TRUE----------------------------------
ltmp <- read_csv('../data/ltmp.csv', trim_ws=TRUE)
glimpse(ltmp)


## -----------------------------------------------------------------------------
#| label: examinData
glimpse(ltmp)


## -----------------------------------------------------------------------------
## Explore the first 6 rows of the data
head(ltmp)


## -----------------------------------------------------------------------------
str(ltmp)


## -----------------------------------------------------------------------------
ltmp |> datawizard::data_codebook()


## -----------------------------------------------------------------------------
ltmp |> modelsummary::datasummary_skim(by = "SITE_NO")


## ----processData, results='markdown', eval=TRUE, mhidden=TRUE-----------------
ltmp_sub <- ltmp |>
    filter(AIMS_REEF_NAME == "Agincourt Reef No.1") |>
    droplevels() |>
    mutate(
        REEF_SITE = factor(paste(AIMS_REEF_NAME, SITE_NO)),
        REEF_SITE_TRANSECT = factor(paste(REEF_SITE, TRANSECT_NO))
    )


## ----EDA2b, results='markdown', eval=TRUE, fig.width=12, fig.height=5, mhidden=TRUE----

ltmp_sub |>
    ggplot(aes(y = HC, x = REPORT_YEAR)) +
  geom_point()

ltmp_sub |>
    ggplot(aes(y = HC, x = REPORT_YEAR)) +
    geom_line(aes(group = REEF_SITE_TRANSECT, colour = REEF_SITE), alpha = 1) +
  geom_point()

give.n <- function(val, ypos){
  return(data.frame(y = ypos, label = round(mean(val), 1)))
  }
ltmp_sub |>
    ggplot(aes(y = HC, x = factor(REPORT_YEAR))) +
    geom_violin(fill = "orange") +
    geom_line(aes(group = REEF_SITE_TRANSECT), alpha = 0.2) +
    geom_point() +
    stat_summary(
      geom = "text",
        aes(y = total.points, ymax =HC),
      fun.data = give.n,
        fun.args = list(ypos = 0)
    )
    ##   fun = mean, aes(y = total.points)
    ## )


## ----processData2, results='markdown', eval=TRUE, mhidden=TRUE----------------
ltmp_sub <- ltmp_sub |>
    mutate(
        fREPORT_YEAR = factor(REPORT_YEAR, levels = rev(sort(unique(REPORT_YEAR))))
    )


## ----fitModel2a, results='markdown', eval=TRUE, cache=TRUE, paged.print=FALSE, tidy.opts = list(width.cutoff = 80)----
ltmp_sub.form <- bf(
    n.points | trials(total.points) ~ fREPORT_YEAR +
        (1 | REEF_SITE) + (1 | REEF_SITE_TRANSECT),
    family = binomial(link = "logit")
)
options(width=150)
ltmp_sub.form %>% get_prior(data = ltmp_sub)
options(width=80)


## ----fitModel2h, results='markdown', eval=TRUE, cache=FALSE-------------------
ltmp_sub |>
    group_by(fREPORT_YEAR) |>
    summarise(
        Median = median(qlogis(n.points / total.points)),
        MAD = mad(qlogis(n.points / total.points)),
        N = mean(total.points)
    )


## ----fitModel2h1, results='markdown', eval=TRUE, cache=TRUE-------------------
priors <- prior(normal(-0.7, 0.3), class = 'Intercept') +
    prior(normal(0, 2), class = 'b') +
    prior(student_t(3, 0, 0.5), class = 'sd')
ltmp_sub.form <- bf(
    n.points | trials(total.points) ~ fREPORT_YEAR +
        (1 | REEF_SITE) + (1 | REEF_SITE_TRANSECT),
    family = binomial(link = "logit")
)
ltmp_sub.brm2 <- brm(ltmp_sub.form,
                  data = ltmp_sub,
                  prior = priors,
                  sample_prior = 'only',
                  iter = 10000,
                  warmup = 5000,
                  chains = 3, cores = 3,
                  thin = 10,
                  refresh = 0,
                  control = list(adapt_delta = 0.99, max_treedepth = 20),
                  backend = "cmdstanr"
                  )



## ----partialPlot2h1a, results='markdown', eval=TRUE, fig.width=8, fig.height=5----

ltmp_sub.brm2 |>
  conditional_effects() |>
  plot(points = TRUE)
ltmp_sub.brm2 |>
  conditional_effects(conditions = data.frame(total.points = 200)) |>
  plot(points = TRUE)

ltmp_sub.brm2 |>
  conditional_effects(conditions = data.frame(total.points = 200)) |>
  plot()


## ----fitModel2h1b, results='markdown', eval=TRUE, cache=TRUE------------------
ltmp_sub.brm3 <- update(ltmp_sub.brm2,
  sample_prior = 'yes',
  chains = 3, cores = 3,
  control = list(adapt_delta = 0.99, max_treedepth = 20),
  thin = 10, iter = 10000, warmup = 5000,
  refresh = 100)
save(ltmp_sub.brm3, file = '../ws/testing/ltmp_sub.brm3')


## ----partialPlot2h1b, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
ltmp_sub.brm3 |>
  conditional_effects(conditions = data.frame(total.points = 200)) |>
  plot(points = TRUE)


## ----posterior2h2, results='markdown', eval=TRUE------------------------------
ltmp_sub.brm3 |> get_variables()
ltmp_sub.brm3 |> hypothesis('fREPORT_YEAR2021=0') %>% plot


## ----posterior2h2a, results='markdown', eval=TRUE, fig.width = 7, fig.height = 5----
## ltmp_sub.brm3 %>% SUYR_prior_and_posterior()


## ----modelValidation2a, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
available_mcmc()


## ----modelValidation2b, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
pars <- ltmp_sub.brm3 |>
    get_variables() |>
    str_subset("^b_.*|^sd.*")
pars
ltmp_sub.brm3 |> mcmc_plot(type='trace', variable = pars)


## ----modelValidation2c, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
ltmp_sub.brm3 |> mcmc_plot(type='acf_bar', variable = pars)


## ----modelValidation2d, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
ltmp_sub.brm3 |> mcmc_plot(type='rhat_hist')


## ----modelValidation2e, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
ltmp_sub.brm2 |> mcmc_plot(type='neff_hist')


## ----modelValidation2f, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
ltmp_sub.brm3 |> mcmc_plot(type='combo', variable = pars)
ltmp_sub.brm3 |> mcmc_plot(type='violin', variable = pars)


## ----modelValidation2g, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
pars <- ltmp_sub.brm3 |>
    get_variables() |>
    str_subset("^b_.*|^sd.*")
pars
ltmp_sub.brm3$fit |>
    stan_trace(pars = pars)


## ----modelValidation2h, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
ltmp_sub.brm3$fit |>
    stan_ac(pars = pars)


## ----modelValidation2i, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
ltmp_sub.brm3$fit |> stan_rhat()


## ----modelValidation2j, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
ltmp_sub.brm3$fit |> stan_ess()


## ----modelValidation2k, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
ltmp_sub.brm3$fit |>
    stan_dens(separate_chains = TRUE, pars = pars)


## ----modelValidation5a, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation5b, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
ltmp_sub.brm3 |> pp_check(type = 'dens_overlay', ndraws = 100)


## ----modelValidation5c, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
## ltmp_sub.brm3 %>% pp_check(type = 'error_scatter_avg')


## ----modelValidation5e, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
ltmp_sub.brm3 %>% pp_check(group = 'REEF', type = 'intervals')


## ----modelValidation5g, results='markdown', eval=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(ltmp_sub.brm2)


## ----modelValidation6aa, results='markdown', eval=TRUE, fig.width=8, fig.height=10----

ltmp_sub.resids <- make_brms_dharma_res(ltmp_sub.brm3, integerResponse = TRUE)
wrap_elements(~testUniformity(ltmp_sub.resids)) +
               wrap_elements(~plotResiduals(ltmp_sub.resids, form = factor(rep(1, nrow(ltmp_sub))))) +
               wrap_elements(~plotResiduals(ltmp_sub.resids, quantreg = TRUE)) +
               wrap_elements(~testDispersion(ltmp_sub.resids))


## ----modelValidation6aaa, results='markdown', eval=TRUE, fig.width=8, fig.height=10----
resids2 <- recalculateResiduals(ltmp_sub.resids, group = unique(ltmp_sub$fREPORT_YEAR))
testTemporalAutocorrelation(resids2, time = unique(ltmp_sub$fREPORT_YEAR))

library(geoR)
autocor_check(ltmp_sub, ltmp_sub.brm3, variable =  "fREPORT_YEAR", n.sim =  250)




## ----summariseModel2a, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
ltmp_sub.brm3 %>% summary()


## ----summariseModel2a1, results='markdown', eval=TRUE, fig.width=8, fig.height=5, echo=FALSE----
ltmp_sub.sum <- summary(ltmp_sub.brm3)


## ----summariseModel2bm, results='markdown', eval=TRUE, fig.width=8, fig.height=5,echo=TRUE----
ltmp_sub.brm3 |> as_draws_df()
ltmp_sub.brm3 |>
    as_draws_df() |>
    dplyr::select(matches("^b_.*|^sd_.*")) |>
    mutate(across(matches("b_.*"), plogis)) |>
  summarise_draws(
    median,
    HDInterval::hdi,
    Pl = ~mean(.x < 1),
    Pg = ~mean(.x > 1),
    rhat,
    ess_bulk
  )


## ----summariseModel2g, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
ltmp_sub.brm3 |>
    bayes_R2(re.form = NA, summary=FALSE) |>
    median_hdci()
ltmp_sub.brm3 |>
    bayes_R2(re.form = ~(1|REEF_SITE), summary=FALSE) |>
    median_hdci()
ltmp_sub.brm3 |>
    bayes_R2(re.form = ~(1|REEF_SITE) + (1| REEF_SITE_TRANSECT), summary=FALSE) |>
    median_hdci()


## ----postHoc1a, results='markdown', eval=TRUE, echo=TRUE----------------------
ltmp_sub.brm3 |>
    emmeans(~fREPORT_YEAR, type = "response") |>
    as.data.frame() |>
    mutate(
        fREPORT_YEAR = factor(fREPORT_YEAR, levels = rev(levels(fREPORT_YEAR))),
        REPORT_YEAR = as.numeric(as.character(fREPORT_YEAR))
    ) |>
    ggplot(aes(y = prob, x = REPORT_YEAR)) +
    geom_ribbon(aes(ymin = lower.HPD, ymax = upper.HPD), fill = "orange", alpha = 0.3) +
    geom_line() +
    geom_point() +
    theme_classic() +
    scale_y_continuous("Live hard coral cover", label = scales::percent_format())



## ----postHoc1b, results='markdown', eval=TRUE, echo=TRUE, width = 10, height = 3----
newdata <- list(fREPORT_YEAR = c(2011, 2012))
ltmp_sub.brm3 |>
  emmeans(~fREPORT_YEAR,
    at = newdata,
    type = "response"
    )

ltmp_sub.brm3 |>
  emmeans(~fREPORT_YEAR,
    at = newdata,
    type = "response"
  ) |>
  pairs(reverse = TRUE)


newdata <- list(fREPORT_YEAR = c(2000:2005, 2012:2013))
cmat <- cbind(c(rep(1 / 6, 6), rep(-1 / 2, 2)))
ltmp_sub.brm3 |>
    emmeans(~fREPORT_YEAR,
        at = newdata,
        type = "response"
    ) |>
    regrid() |>
  contrast(method = list(fREPORT_YEAR = cmat))



ltmp_sub.brm3 |>
    emmeans(~fREPORT_YEAR,
        at = list(fREPORT_YEAR = c(2011, 2012)),
        type = "link"
    ) |>
    pairs(reverse = TRUE) |>
    tidy_draws() |>
    exp() |>
    summarise_draws(
      median,
      HDInterval::hdi,
      Pl = ~ mean(.x < 1),
      Pg = ~ mean(.x > 1)
    )


cmat <- cbind(c(rep(-1 / 7, 7), 1))
ltmp_sub.brm3 |>
    emmeans(~fREPORT_YEAR,
        at = list(fREPORT_YEAR = c(1999:2005, 2012)),
        type = "response"
    ) |>
    contrast(method = list(fREPORT_YEAR = cmat))


ltmp_sub.brm3 |>
    emmeans(~fREPORT_YEAR,
        at = list(fREPORT_YEAR = c(1999:2005, 2012)),
        type = "response"
    ) |>
    regrid() |>
    contrast(method = list(fREPORT_YEAR = cmat))


## -----------------------------------------------------------------------------
#| label: as a gam
#| results: markup
#| eval: true
#| echo: true
#| cache: false
ltmp.form <- bf(
    n.points | trials(total.points) ~ s(REPORT_YEAR) +
        (1 | REEF_SITE) + (1 | REEF_SITE_TRANSECT),
    family = binomial(link = "log")
)
get_prior(ltmp.form, data = ltmp_sub)
ltmp_sub |>
    summarise(
        Median = median(qlogis(HC/100)),
        MAD = mad(qlogis(HC/100)),
        N = mean(total.points)
    )
priors <- prior(normal(-1.1, 0.8), class = "Intercept") +
    prior(normal(0, 1), class = 'b') +
    prior(student_t(3, 0, 1), class = 'sd') +
    prior(student_t(3, 0, 10), class = 'sds')
ltmp_sub.brm4 <- brm(ltmp.form,
                  data = ltmp_sub,
                  prior = priors,
                  sample_prior = 'yes',
                  iter = 5000,
                  warmup = 2000,
                  chains = 3, cores = 3,
                  thin = 10,
                  refresh = 100,
                  control = list(adapt_delta = 0.99, max_treedepth = 20),
                  backend = "cmdstanr"
                  )
## save(ltmp_sub.brm4, file = '../ws/testing/ltmp_sub.brm4')
ltmp_sub.brm4 |>
  conditional_effects(conditions = data.frame(total.points = 200)) |>
  plot(points = TRUE)



data.frame(mgcv::smoothCon(s(REPORT_YEAR), data = ltmp_sub)[[1]]$X) |>
    head()
## data.frame(mgcv::smoothCon(s(REPORT_YEAR, k=3),  data=ltmp_sub)[[1]]$REPORT_YEAR) |>
##   bind_cols(ltmp_sub)

gratia::basis(s(REPORT_YEAR),  data = ltmp_sub) |> gratia::draw()
gratia::basis(s(REPORT_YEAR, bs='cr'),  data = ltmp_sub) |> gratia::draw()
gratia::basis(s(REPORT_YEAR, bs='cr', k = 3),  data = ltmp_sub) |> gratia::draw()


ltmp_glmm <- ltmp_sub.brm3 |>
    emmeans(~fREPORT_YEAR, type = "response") |>
    as.data.frame() |>
    mutate(
        fREPORT_YEAR = factor(fREPORT_YEAR, levels = rev(levels(fREPORT_YEAR))),
        REPORT_YEAR = as.numeric(as.character(fREPORT_YEAR))
    )

ltmp_sub.brm4 |>
  conditional_effects() |>
  plot() |>
  _[[1]] +
  geom_line(data = ltmp_glmm, inherit.aes = FALSE,
    aes(y = prob, x = REPORT_YEAR), colour = "blue") +
  geom_ribbon(data = ltmp_glmm, inherit.aes = FALSE,
    aes(y = prob, x = REPORT_YEAR,
      min = lower.HPD, ymax = upper.HPD), fill = "orange", alpha = 0.3) +
  theme_classic() +
  scale_y_continuous("Live hard coral cover", label = scales::percent_format())

