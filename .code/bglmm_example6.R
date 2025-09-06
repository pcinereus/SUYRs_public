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


## \tikzstyle{Messy} = [decorate,decoration={random steps,segment length=3pt, amplitude=0.5pt}]
## \tikzstyle{HandTitle} = [font={\fontspec[Scale=2.1]{CabinSketch}}]
## \tikzstyle{HandBoxTitle} = [font={\fontspec[Scale=1.5]{Complete in Him}}]
## \tikzstyle{HandLabel} = [font={\fontspec[Scale=1.1]{Hannahs Messy Handwriting}}]
## \tikzstyle{Plot} = [rectangle,draw,Messy,fill=white,HandLabel, minimum height=2em]
## \tikzstyle{Plate} = [circle,draw,Messy,fill=white,HandLabel, minimum height=5.5cm]
## \tikzstyle{Dist} = [rectangle,draw,Messy]
## 
## \pgfdeclarelayer{Plates}
## \pgfdeclarelayer{Dists}
## \pgfsetlayers{Plates,Dists,main}
## 
## \newcommand{\mybox}[2][]{
## \node[Plate,fill=blue!2] (Plate1) {#1};
## \path  (Plate1.north) +(0,-0.3) node [HandLabel] (Dist4Title) {\textbf{Dist 4}};
## % \node[Plot1, right of=Plot1, anchor=east,node distance=2.5cm] (Plot2) {#1};
## \node[Plate,fill=blue!5,minimum height=4.5cm] (Dist3) {#1};
## \path  (Dist3.north) +(0,-0.3) node [HandLabel] (Dist3Title) {\textbf{Dist 3}};
## \node[Plate,fill=blue!10,minimum height=3.5cm] (Dist2) {#1};
## \path  (Dist2.north) +(0,-0.3) node [HandLabel] (Dist2Title) {\textbf{Dist 2}};
## \node[Plate,fill=blue!20,minimum height=2.5cm] (Dist1) {#1};
## \path  (Dist1.north) +(0,-0.3) node [HandLabel] (Dist1Title) {\textbf{Dist 1}};
## \node[Plate,fill=blue!30, minimum height=1.5cm] (Core) {#1};
## % \draw (0,0) circle [Messy,radius=2.25];
## \begin{pgfonlayer}{Plates}
## \path (Plate1.west |- Plate1.north) +(-0.2,+0.5) node (S1nw) {};
## \path (Plate1.east |- Plate1.south) +(+0.2,-0.2) node (S1se) {};
## %% title
## \path  ($ (S1nw.west |- S1nw.north) !0.5! (S1se.east |- S1nw.north)$) +(0,-0) node [HandBoxTitle] (Plate1Title) {\textbf{Plates #2}};
## \path  ($ (S1nw.west |- S1nw.north) !0.5! (S1se.east |- S1nw.north)$) +(0,-0.35) node [HandBoxTitle] (Plate1Title) {\textbf{#1}};
## 
## \end{pgfonlayer}
## }
## 
## 
## \begin{tikzpicture} \path node (Plates1) {
## \begin{tikzpicture}
## \mybox[Control]{1,2,3,4,5}
## \end{tikzpicture}
## };
## \path (Plates1.east) +(3,0) node (Plates2) {
## \begin{tikzpicture}
## \mybox[Week 1]{6,7,8,9,10}
## \end{tikzpicture}
## };
## \path (Plates2.east) +(3,0) node (Plates3) {
## \begin{tikzpicture}
## \mybox[Week 2]{11,12,13,14,15}
## \end{tikzpicture}
## };
## 
## \end{tikzpicture}

## -----------------------------------------------------------------------------
#| label: readData
copper <- read_csv("../data/copper.csv", trim_ws = TRUE)


## -----------------------------------------------------------------------------
#| label: examinData
glimpse(copper)


## -----------------------------------------------------------------------------
## Explore the first 6 rows of the data
head(copper)


## -----------------------------------------------------------------------------
str(copper)


## -----------------------------------------------------------------------------
copper |> datawizard::data_codebook()


## -----------------------------------------------------------------------------
copper |> modelsummary::datasummary_skim()
copper |> modelsummary::datasummary_skim(by = c("COPPER", "DIST"))


## ----dataProcessing, results='markdown', eval=TRUE, mhidden=TRUE--------------
copper <- copper |> mutate(
    COPPER = factor(COPPER),
    PLATE = factor(PLATE),
    DIST = factor(DIST),
    AREA = 4,
    COUNT = WORMS * AREA)


## ----eda1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=7, fig.height=5----
ggplot(copper, aes(y=WORMS, x=DIST, fill=COPPER)) +
  geom_boxplot()

ggplot(copper, aes(y=WORMS, x=DIST, fill=COPPER)) +
    geom_boxplot() +
    scale_y_continuous(trans = scales::pseudo_log_trans())

ggplot(copper, aes(y = WORMS, x = DIST, colour = COPPER)) +
    geom_point(aes(x = as.numeric(DIST))) +
    geom_line(aes(x = as.numeric(DIST), group = PLATE)) +
    scale_y_continuous(trans = scales::pseudo_log_trans())


## ----eda2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=7, fig.height=5----
copper |> filter(WORMS > 0) |>
    summarise(min(WORMS)/2)


## ----fitModel2a, results='markdown', eval=TRUE, mhidden=TRUE, paged.print=FALSE, tidy.opts = list(width.cutoff = 80)----
copper.form <- bf(WORMS ~ COPPER * DIST + (1|PLATE),
                family=Gamma(link='log'))
options(width=150)
copper.form |> get_prior(data = copper)
options(width=80)


## ----fitModel2h, results='markdown', eval=TRUE, mhidden=TRUE------------------
copper |>
    group_by(COPPER, DIST) |>
    summarise(log(median(WORMS)),
              log(mad(WORMS)))
standist::visualize("normal(3,0.45)", xlim=c(0,20))
standist::visualize("student_t(3, 0, 2.5)",
                    "cauchy(0,2)",
                    xlim=c(-10,25))


## ----fitModel2h0, results='markdown', eval=TRUE, mhidden=TRUE-----------------
priors <- prior(normal(2.5, 0.6), class = 'Intercept') +
    prior(normal(0, 1.5), class = 'b') +
    prior(student_t(3, 0, 1), class = 'sd')

copper.form <- bf(COUNT~ offset(log(AREA)) + COPPER * DIST + (1|PLATE),
                family = poisson(link = "log"))
copper.brm2 <- brm(copper.form,
                 data = copper,
                 prior = priors,
                 sample_prior = 'only',
                 iter = 5000,
                 warmup =2500,
                 chains = 3,
                 cores = 3,
                 thin = 10,
                 refresh = 0,
                 seed = 123,
                 control = list(adapt_delta = 0.99)
                 )



## ----partialPlot2h0a, results='markdown', warning=FALSE, message=FALSE, eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
copper.brm2 |>
    conditional_effects("COPPER:DIST") |>
    plot(points = TRUE)
copper.brm2 |>
    ggpredict(~COPPER*DIST) |>
    plot(show_data = TRUE)


## ----fitModel2h0ab, results='markdown', warning=FALSE, eval=TRUE, mhidden=TRUE----
copper.brm3 <- update(copper.brm2,
                       sample_prior = 'yes',
                       control = list(adapt_delta = 0.99),
                       refresh = 0)
save(copper.brm3, file = '../ws/testing/copper.brm3')


## ----partialPlot2h00a, results='markdown', warning=FALSE, message=FALSE, eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
copper.brm3 |>
    conditional_effects("COPPER:DIST") |>
    plot(points = TRUE)
copper.brm3 |>
    ggpredict(~COPPER*DIST) |>
    plot(show_data = TRUE)


## ----fitModel2h1, results='markdown', eval=TRUE, mhidden=TRUE-----------------
priors <- prior(normal(3, 0.45), class = 'Intercept') +
    prior(normal(0, 0.9), class = 'b', coef = 'COPPERWeek1') +
    prior(normal(0, 0.9), class = 'b', coef = 'COPPERWeek2') +
    prior(normal(0, 1), class = 'b', coef = 'DIST2') +
    prior(normal(0, 1), class = 'b', coef = 'DIST3') +
    prior(normal(0, 1), class = 'b', coef = 'DIST4') +
    prior(normal(0, 1.5), class = 'b') +
    prior(cauchy(0,2), class = 'sd') +
    prior(gamma(2,1), class = "sigma")

copper.form <- bf(I(WORMS + 0.125)~ COPPER * DIST + (1|PLATE),
                family=lognormal())
copper.brm2a <- brm(copper.form,
                 data = copper,
                 prior = priors,
                 sample_prior = 'only',
                 iter = 5000,
                 warmup =2500,
                 chains = 3,
                 cores = 3,
                 thin = 10,
                 refresh = 0,
                 seed = 123,
                 control = list(adapt_delta = 0.99)
                 )



## ----partialPlot2h1a, results='markdown', warning=FALSE, message=FALSE, eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
copper.brm2a |>
    ggpredict(~COPPER*DIST) |>
    plot(show_data = TRUE)


## ----fitModel2h1ab, results='markdown', warning=FALSE, eval=TRUE, mhidden=TRUE----
copper.brm3a <- update(copper.brm2a,
                       sample_prior = 'yes',
                       control = list(adapt_delta = 0.99),
                       refresh = 0)
save(copper.brm3a, file = '../ws/testing/copper.brm3a')


## ----fitModel2h1b, results='markdown', eval=TRUE, mhidden=TRUE----------------
priors <- prior(normal(3, 0.45), class = 'Intercept') +
    prior(normal(0, 0.9), class = 'b', coef = 'COPPERWeek1') +
    prior(normal(0, 0.9), class = 'b', coef = 'COPPERWeek2') +
    prior(normal(0, 1), class = 'b', coef = 'DIST2') +
    prior(normal(0, 1), class = 'b', coef = 'DIST3') +
    prior(normal(0, 1), class = 'b', coef = 'DIST4') +
    prior(normal(0, 1.5), class = 'b') +
    prior(cauchy(0,1), class = 'sd') +
    prior(gamma(0.01, 0.01), class = "shape")

copper.form <- bf(I(WORMS + 0.125)~ COPPER * DIST + ( DIST |PLATE),
                family=Gamma(link='log'))
copper.brm2b <- brm(copper.form,
                 data = copper,
                 prior = priors,
                 sample_prior = 'only',
                 iter = 5000,
                 warmup =2500,
                 chains = 3,
                 cores = 3,
                 thin = 10,
                 refresh = 0,
                 seed = 123,
                 control = list(adapt_delta = 0.99)
                 )



## ----fitModel2h2b, results='markdown', eval=TRUE, mhidden=TRUE----------------
copper.brm3b <- update(copper.brm2b,
                       sample_prior = 'yes',
                       control = list(adapt_delta = 0.99),
                       refresh = 0)
save(copper.brm3b, file = '../ws/testing/copper.brm3b')


## ----partialPlot2h2b, results='markdown', warning=FALSE, message=FALSE, eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
copper.brm3b |>
    ggpredict(~COPPER*DIST) |>
    plot(show_data = TRUE)
copper.brm3b |>
    ggpredict(~DIST*COPPER) |>
    plot(show_data = TRUE)


## ----posterior2h0, results='markdown', eval=TRUE------------------------------
copper.brm3 |> get_variables()
copper.brm3 |> hypothesis('COPPERWeek1=0') |> plot()
copper.brm3 |> hypothesis('DIST2=0') |> plot()


## ----posterior2h0a, results='markdown', eval=TRUE, fig.width=10, fig.height=10----
copper.brm3 |> SUYR_prior_and_posterior()


## ----posterior2h2, results='markdown', eval=TRUE------------------------------
copper.brm3a |> get_variables()
copper.brm3a |> hypothesis('COPPERWeek1=0') |> plot()
copper.brm3a |> hypothesis('DIST2=0') |> plot()


## ----posterior2h2a, results='markdown', eval=TRUE, fig.width=10, fig.height=10----
copper.brm3a |> SUYR_prior_and_posterior()


## ----posterior2i2, results='markdown', eval=TRUE------------------------------
copper.brm3b |> get_variables()
copper.brm3b |> hypothesis('COPPERWeek1=0') |> plot()
copper.brm3b |> hypothesis('DIST2=0') |> plot()


## ----posterior2i2a, results='markdown', eval=TRUE, fig.width=10, fig.height=10----
copper.brm3b |> SUYR_prior_and_posterior()


## ----modelValidation2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
available_mcmc()


## ----modelValidation2b0, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=10, fig.height=10----
pars <- copper.brm3 |> get_variables()
pars <- pars |> str_extract('^b.Intercept|^b_COPPER.*|^b_DIST.*|[sS]igma|^sd.*') |>
    na.omit()
pars
copper.brm3 |> mcmc_plot(type='trace', variable = pars)
##OR
copper.brm3 |> mcmc_plot(type='trace',
                        variable = '^b.Intercept|^b_COPPER.*|^b_DIST.*|[sS]igma|^sd.*',
                        regex = TRUE)



## ----modelValidation2c0, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=10, fig.height=4----
copper.brm3 |> mcmc_plot(type='acf_bar', variable = pars)
##OR
copper.brm3 |> mcmc_plot(type='acf_bar',
                        variable = '^b.Intercept|^b_COPPER.*|^b_DIST.*|[sS]igma|^sd.*',
                        regex = TRUE)


## ----modelValidation2d0, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3 |> mcmc_plot(type='rhat_hist')


## ----modelValidation2e0, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3 |> mcmc_plot(type='neff_hist')


## ----modelValidation2f0, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3 |> mcmc_plot(type='combo', pars = pars)
copper.brm3 |> mcmc_plot(type='violin', pars = pars)


## ----modelValidation2b2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=10, fig.height=10----
pars <- copper.brm3a |> get_variables()
pars <- pars |> str_extract('^b.Intercept|^b_COPPER.*|^b_DIST.*|[sS]igma|^sd.*') |>
    na.omit()
pars
copper.brm3a |> mcmc_plot(type='trace', variable = pars)
##OR
copper.brm3a |> mcmc_plot(type='trace',
                        variable = '^b.Intercept|^b_COPPER.*|^b_DIST.*|[sS]igma|^sd.*',
                        regex = TRUE)



## ----modelValidation2c2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=10, fig.height=4----
copper.brm3a |> mcmc_plot(type='acf_bar', variable = pars)
##OR
copper.brm3a |> mcmc_plot(type='acf_bar',
                        variable = '^b.Intercept|^b_COPPER.*|^b_DIST.*|[sS]igma|^sd.*',
                        regex = TRUE)


## ----modelValidation2d2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3a |> mcmc_plot(type='rhat_hist')


## ----modelValidation2e2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3a |> mcmc_plot(type='neff_hist')


## ----modelValidation2f2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3a |> mcmc_plot(type='combo', pars = pars)
copper.brm3a |> mcmc_plot(type='violin', pars = pars)


## ----modelValidation2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=10, fig.height=10----
pars <- copper.brm3b |> get_variables()
pars <- pars |> str_extract('^b.Intercept|^b_COPPER.*|^b_DIST.*|[sS]igma|^sd.*') |>
    na.omit()
pars
copper.brm3b |> mcmc_plot(type='trace', variable = pars)
##OR
copper.brm3b |> mcmc_plot(type='trace',
                        variable = '^b.Intercept|^b_COPPER.*|^b_DIST.*|[sS]igma|^sd.*',
                        regex = TRUE)



## ----modelValidation2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=10, fig.height=4----
copper.brm3b |> mcmc_plot(type='acf_bar', variable = pars)
##OR
copper.brm3b |> mcmc_plot(type='acf_bar',
                        variable = '^b.Intercept|^b_COPPER.*|^b_DIST.*|[sS]igma|^sd.*',
                        regex = TRUE)


## ----modelValidation2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3b |> mcmc_plot(type='rhat_hist')


## ----modelValidation2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3b |> mcmc_plot(type='neff_hist')


## ----modelValidation2f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3b |> mcmc_plot(type='combo', pars = pars)
copper.brm3b |> mcmc_plot(type='violin', pars = pars)


## ----modelValidation2g0, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=10, fig.height=10----
copper.brm3 |> get_variables()
pars <- copper.brm3 |> get_variables()
pars <- str_extract(pars, '^b_.*|^sigma$|^sd.*') |> na.omit()

copper.brm3$fit |>
    stan_trace(pars = pars)


## ----modelValidation2h0, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=10, fig.height=10----
copper.brm3$fit |>
    stan_ac(pars = pars)


## ----modelValidation2i0, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3$fit |> stan_rhat()


## ----modelValidation2j0, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3$fit |> stan_ess()


## ----modelValidation2k0, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=10, fig.height=10----
copper.brm3$fit |>
    stan_dens(separate_chains = TRUE, pars = pars)


## ----modelValidation2g2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=10, fig.height=10----
copper.brm3a |> get_variables()
pars <- copper.brm3a |> get_variables()
pars <- str_extract(pars, '^b_.*|^sigma$|^sd.*') |> na.omit()

copper.brm3a$fit |>
    stan_trace(pars = pars)


## ----modelValidation2h2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=10, fig.height=10----
copper.brm3a$fit |>
    stan_ac(pars = pars)


## ----modelValidation2i2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3a$fit |> stan_rhat()


## ----modelValidation2j2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3a$fit |> stan_ess()


## ----modelValidation2k2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=10, fig.height=10----
copper.brm3a$fit |>
    stan_dens(separate_chains = TRUE, pars = pars)


## ----modelValidation2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=10, fig.height=10----
copper.brm3b |> get_variables()
pars <- copper.brm3b |> get_variables()
pars <- str_extract(pars, '^b_.*|^sigma$|^sd.*') |> na.omit()

copper.brm3b$fit |>
    stan_trace(pars = pars)


## ----modelValidation2h, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=10, fig.height=10----
copper.brm3b$fit |>
    stan_ac(pars = pars)


## ----modelValidation2i, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3b$fit |> stan_rhat()


## ----modelValidation2j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3b$fit |> stan_ess()


## ----modelValidation2k, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=10, fig.height=10----
copper.brm3b$fit |>
    stan_dens(separate_chains = TRUE, pars = pars)


## ----modelValidation2l, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=7----
## owls.ggs <- copper.brm3 |> ggs(burnin = FALSE, inc_warmup = FALSE)
## copper.ggs |> ggs_traceplot()


## ----modelValidation2m, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=7----
## ggs_autocorrelation(copper.ggs)


## ----modelValidation2n, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## ggs_Rhat(copper.ggs)


## ----modelValidation2o, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## ggs_effective(copper.ggs)


## ----modelValidation2p, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## ggs_crosscorrelation(owls.ggs)


## ----modelValidation2q, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## ggs_grb(owls.ggs)


## ----modelValidation5a0, results='markdown', eval=FALSE, mhidden=TRUE, fig.width=6, fig.height=4----
# available_ppc()


## ----modelValidation5b0, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3 |> pp_check(type = 'dens_overlay', ndraws = 100)


## ----modelValidation5c0, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## copper.brm3 |> pp_check(type = 'error_scatter_avg')


## ----modelValidation5e0, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3 |> pp_check(group = 'Nest', type = 'intervals')


## ----modelValidation5g0, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(copper.brm2)


## ----modelValidation5a2, results='markdown', eval=FALSE, mhidden=TRUE, fig.width=6, fig.height=4----
# available_ppc()


## ----modelValidation5b2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3a |> pp_check(type = 'dens_overlay', ndraws = 100)


## ----modelValidation5c2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## copper.brm3a |> pp_check(type = 'error_scatter_avg')


## ----modelValidation5e2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3a |> pp_check(group = 'Nest', type = 'intervals')


## ----modelValidation5g2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(copper.brm2)


## ----modelValidation5a, results='markdown', eval=FALSE, mhidden=TRUE, fig.width=6, fig.height=4----
# available_ppc()


## ----modelValidation5b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3b |> pp_check(type = 'dens_overlay', ndraws = 100)


## ----modelValidation5c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
## copper.brm3b |> pp_check(type = 'error_scatter_avg')


## ----modelValidation5e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
copper.brm3b |> pp_check(group = 'Nest', type = 'intervals')


## ----modelValidation5g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(copper.brm2)


## ----modelValidation6a0, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=10, fig.height=10----
copper.resids <- make_brms_dharma_res(copper.brm3, integerResponse = TRUE)
wrap_elements(~testUniformity(copper.resids)) +
               wrap_elements(~plotResiduals(copper.resids, form = factor(rep(1, nrow(copper))))) +
               wrap_elements(~plotResiduals(copper.resids, quantreg = TRUE)) +
               wrap_elements(~testDispersion(copper.resids))


## ----validation2g0, results='markdown', eval=TRUE, error=TRUE,mhidden=TRUE, fig.width=7, fig.height=5, cache=FALSE, message=FALSE, warning=FALSE----
try({
copper.resids |> testZeroInflation()
})


## ----modelValidation6a2, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=10, fig.height=10----
## preds <- copper.brm3a |> posterior_predict(nsamples = 250,  summary = FALSE)
## copper.resids <- createDHARMa(simulatedResponse = t(preds),
##                             observedResponse = copper$NCalls,
##                             fittedPredictedResponse = apply(preds, 2, median),
##                             integerResponse = TRUE)
## plot(copper.resids)

copper.resids <- make_brms_dharma_res(copper.brm3a, integerResponse = TRUE)
wrap_elements(~testUniformity(copper.resids)) +
               wrap_elements(~plotResiduals(copper.resids, form = factor(rep(1, nrow(copper))))) +
               wrap_elements(~plotResiduals(copper.resids, quantreg = TRUE)) +
               wrap_elements(~testDispersion(copper.resids))


## ----validation2g2, results='markdown', eval=TRUE, error=TRUE,mhidden=TRUE, fig.width=7, fig.height=5, cache=FALSE, message=FALSE, warning=FALSE----
try({
copper.resids |> testZeroInflation()
})


## ----modelValidation6a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=10, fig.height=10----
## preds <- copper.brm3b |> posterior_predict(nsamples = 250,  summary = FALSE)
## copper.resids <- createDHARMa(simulatedResponse = t(preds),
##                             observedResponse = copper$NCalls,
##                             fittedPredictedResponse = apply(preds, 2, median),
##                             integerResponse = TRUE)
## plot(copper.resids)

copper.resids <- make_brms_dharma_res(copper.brm3b, integerResponse = TRUE)
wrap_elements(~testUniformity(copper.resids)) +
               wrap_elements(~plotResiduals(copper.resids, form = factor(rep(1, nrow(copper))))) +
               wrap_elements(~plotResiduals(copper.resids, quantreg = TRUE)) +
               wrap_elements(~testDispersion(copper.resids))


## ----validation2g, results='markdown', eval=TRUE, error=TRUE,mhidden=TRUE, fig.width=7, fig.height=5, cache=FALSE, message=FALSE, warning=FALSE----
try({
copper.resids |> testZeroInflation()
})


## ----partialPlot2d, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
copper.brm3 |>
    conditional_effects("COPPER:DIST") |>
    plot(points = TRUE, jitter_width = 0.25)


## ----partialPlot2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
copper.brm3 |>
    ggpredict(~COPPER*DIST) |>
    plot(show_data = TRUE, jitter = c(0.25, 0))


## ----partialPlot2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
copper.brm3 |>
    ggemmeans(~COPPER*DIST) |>
    plot()


## ----summariseModel2a, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
copper.brm3 |> summary()


## ----summariseModel2a1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5, echo=FALSE----
copper.sum <- summary(copper.brm3)


## ----summariseModel2bm, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
copper.brm3 |> as_draws_df()
copper.brm3 |>
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


## ----summariseModel2b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
copper.brm3$fit |>
    tidyMCMC(estimate.method = 'median',
             conf.int = TRUE,  conf.method = 'HPDinterval',
             rhat = TRUE, ess = TRUE)

## ----summariseModel2b1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
copper.tidy <- tidyMCMC(copper.brm3$fit, estimate.method='median',
                         conf.int=TRUE,  conf.method='HPDinterval',
                         rhat=TRUE, ess=TRUE)


## ----summariseModel2c, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
copper.brm3 |> get_variables()
copper.draw <- copper.brm3 |>
    gather_draws(`b.Intercept.*|b_COPPER.*|b_DIST.*`,  regex=TRUE)
copper.draw


## ----summariseModel2c1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
copper.draw |> median_hdci()
## On a fractional scale
copper.draw |>
    mutate(.value = exp(.value)) |>
    median_hdci()


## ----summariseModel2c3, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
copper.gather <- copper.brm3 |>
    gather_draws(`b_Intercept.*|b_COPPER.*|b_DIST.*`,  regex=TRUE) |>
    mutate(.value = exp(.value)) |>
    median_hdci()


## ----summariseModel2c4, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
copper.brm3 |>
    gather_draws(`b_Intercept.*|b_COPPER.*|b_DIST.*`, regex=TRUE) |>
    ## mutate(.value = exp(.value)) |>
    ggplot() +
    geom_vline(xintercept=0, linetype='dashed') +
    stat_slab(aes(x = .value, y = .variable,
                  fill = stat(ggdist::cut_cdf_qi(cdf,
                                                 .width = c(0.5, 0.8, 0.95),
                                                 labels = scales::percent_format())
                              )), color='black') +
    scale_fill_brewer('Interval', direction = -1, na.translate = FALSE)

copper.brm3 |>
    gather_draws(`.Intercept.*|b_COPPER.*|b_DIST.*`, regex=TRUE) |>
    ggplot() +
    geom_vline(xintercept = 0, linetype='dashed') +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    theme_classic()


## ----summariseModel2j, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
copper.brm3$fit |> plot(type='intervals')


## ----summariseModel2ka, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
copper.brm3 |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(.variable != 'b_Intercept') |>
    ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    facet_wrap(~.variable, scales='free') +
    theme(axis.text.y = element_blank())

copper.brm3 |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(.variable != 'b_Intercept') |>
    ggplot() +
    stat_halfeye(aes(x=.value,  y=.variable)) +
    geom_vline(xintercept = 0, linetype = 'dashed')


## ----summariseModel2c7, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
copper.brm3 |>
    gather_draws(`^b_.*`, regex=TRUE) |>
    filter(str_detect(.variable, 'b_.*Intercept', negate = TRUE)) |>
    ggplot() +
    geom_density_ridges(aes(x=.value, y = .variable), alpha=0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed')
##Or in colour
copper.brm3 |>
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
copper.brm3 |>
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
copper.brm3 |> tidy_draws()


## ----summariseModel2e, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
copper.brm3 |> spread_draws(`.*Intercept.*|b_COPPER.*|b_DIST.*`,  regex=TRUE)


## ----summariseModel2f, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
copper.brm3 |> posterior_samples() |> as_tibble()


## ----summariseModel2g, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=5----
copper.brm3 |>
    bayes_R2(re.form = NA, summary=FALSE) |>
    median_hdci()
copper.brm3 |>
    bayes_R2(re.form = ~(1|PLATE), summary=FALSE) |>
    median_hdci()
copper.brm3 |>
    bayes_R2(re.form = ~(COPPER*DIST|PLATE), summary=FALSE) |>
    median_hdci()


## -----------------------------------------------------------------------------
#| label: modelsummary
#| results: markup
#| eval: false
#| echo: true
#| cache: false
# copper.brm3 |> modelsummary(
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
# copper.brm3 |> modelplot(exponentiate = TRUE)


## ----postHoc1a, results='markdown', eval=TRUE,mhidden=TRUE--------------------
newdata <- copper.brm3 |>
    emmeans(~COPPER|DIST, type='response') |>
    as.data.frame()
head(newdata)
ggplot(newdata) +
    geom_pointrange(aes(y=rate,  x=COPPER,  color=DIST,
                        ymin=lower.HPD,  ymax=upper.HPD),
                    position=position_dodge(width=0.2)) +
    theme_classic()




## ----postHoc2, results='markdown', eval=TRUE,mhidden=TRUE---------------------
copper.brm3 |>
  emmeans(~COPPER|DIST, type='response') |>
  pairs()
copper.brm3 |>
  emmeans(~COPPER|DIST, type='response') |>
  pairs() |>
  tidy_draws() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median, HDInterval::hdi,
                  Pl = ~ mean(.x < 1), Pg = ~ mean(.x > 1))


copper.brm3 |>
  emmeans(~DIST | COPPER, type='response') |>
  pairs()
copper.brm3 |>
  emmeans(~DIST | COPPER, type='response') |>
  pairs() |>
  tidy_draws() |>
  mutate(across(everything(), exp)) |>
  summarise_draws(median, HDInterval::hdi,
                  Pl = ~ mean(.x < 1), Pg = ~ mean(.x > 1))



