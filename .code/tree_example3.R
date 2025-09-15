## -----------------------------------------------------------------------------
#| label: setup
#| include: false

knitr::opts_chunk$set(cache = TRUE, cache.lazy = FALSE, tidy='styler')


## -----------------------------------------------------------------------------
#| label: libraries
#| output: false
#| eval: true
#| warning: false
#| message: false
#| cache: false

library(gbm)         #for gradient boosted models
## library(gbm3)         #for gradient boosted models
library(car)
library(pdp)
library(ggfortify)
library(randomForest)
library(tidyverse)
library(patchwork)


## ----readData, results='markdown', eval=TRUE----------------------------------
newman <- read_csv('../data/newman.csv', trim_ws=TRUE)
glimpse(newman)


## -----------------------------------------------------------------------------
#| label: examinData
#| dependson: readData

newman |> glimpse()


## -----------------------------------------------------------------------------
#| label: headData
#| dependson: readData
## Explore the first 6 rows of the data
newman |> head()


## -----------------------------------------------------------------------------
#| label: strData
#| dependson: readData
newman |> str()


## -----------------------------------------------------------------------------
#| label: easyData
#| dependson: readData
newman |> datawizard::data_codebook()


## -----------------------------------------------------------------------------
#| label: process data
#| results: markup
#| eval: true
#| echo: true
#| cache: false
names(newman) <- str_replace_all(names(newman), " ", "_")
newman <- newman |>
  mutate(
    fCountry = factor(Country),
    nCountry = as.numeric(fCountry),
    fReef_complexity = factor(Reef_complexity)
  )


## ----EDA, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=8, fig.height=8----
#| dependson: readData
scatterplotMatrix(~total_rich + Reef_complexity + No_tall_corals +
                    No_corals + Sponge_max_height + Octocoral_max_height +
                    Slope_angle + nCountry, data = newman,
                    diagonal = list(method='boxplot'))


## ----fitModel1, results='hide', eval=TRUE, mhidden=TRUE, cache=TRUE-----------
newman_gbm <- gbm(fish_rich ~ fReef_complexity +
                   No_tall_corals +
                   No_corals +
                   Sponge_max_height +
                   Octocoral_max_height +
                   Slope_angle +
                  fCountry,
  data = newman,
  distribution = "poisson",
  n.trees = 10000,
  var.monotone = c(0, 1, 1, 1, 0, 1, 0),
  interaction.depth = 5,
  shrinkage = 0.001,
  bag.fraction = 0.5,
  train.fraction = 1,
  n.minobsinnode = 10,
  verbose = FALSE,
  cv.folds = 10,
  n.cores = 1
)


## ----fitModel2, results='markdown', eval=TRUE, mhidden=TRUE-------------------
(best.iter <- gbm::gbm.perf(newman_gbm, method='OOB'))
(best.iter <- gbm::gbm.perf(newman_gbm, method = 'cv'))
#OR
#(best_iter <- gbm3::gbmt_performance(newman_gbm, method = "cv"))
#plot(best_iter)


## ----relativeInfluence1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=7, fig.height=10----
#| dependson: fitModel3
summary(newman_gbm, n.trees = best.iter)
## summary(newman_gbm, num_trees = best.iter)
100 / 7



## ----fitModel2a, results='hide', eval=TRUE, mhidden=TRUE, cache=TRUE----------
newman_gbm <- gbm(fish_rich ~ fReef_complexity +
                   No_tall_corals +
                   No_corals +
                   Sponge_max_height +
                   Octocoral_max_height +
                   Slope_angle,
  data = newman,
  distribution = "poisson",
  n.trees = 10000,
  var.monotone = c(0, 1, 1, 1, 0, 1),
  interaction.depth = 5,
  shrinkage = 0.001,
  bag.fraction = 0.5,
  train.fraction = 1,
  n.minobsinnode = 10,
  verbose = FALSE,
  cv.folds = 10,
  n.cores = 1
)


## ----fitModel2b, results='markdown', eval=TRUE, mhidden=TRUE------------------
(best.iter <- gbm::gbm.perf(newman_gbm, method='OOB'))
(best.iter <- gbm::gbm.perf(newman_gbm, method = 'cv'))
#OR
#(best_iter <- gbm3::gbmt_performance(newman_gbm, method = "cv"))
#plot(best_iter)


## ----relativeInfluence1b, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=7, fig.height=10----
#| dependson: fitModel3
summary(newman_gbm, n.trees = best.iter)
100 / 6



## ----partialEffects1, results='markdown', eval=TRUE, mhidden=TRUE, fig.width=4, fig.height=4----
#| dependson: fitModel3

attr(newman_gbm$Terms,"term.labels")

interact.gbm(newman_gbm, newman, c(1,5), n.tree=best.iter)

attr(newman_gbm$Terms,"term.labels")
plot(newman_gbm, 5, n.tree = best.iter)
plot(newman_gbm, 1, n.tree = best.iter, log = 'x')
## plot(newman_gbm, 7, num.tree = best.iter)
## plot(newman_gbm, 1, num.tree = best.iter, log = 'x')
newman_gbm |>
    pdp::partial(pred.var='fReef_complexity',
                 n.trees = best.iter,
                 type = "regression",
                 recursive = FALSE,
                 ice = TRUE,
                 inv.link = exp) |>
  autoplot()
## newman_gbm |>
##     pdp::partial(pred.var='fCountry',
##                  n.trees = best.iter,
##                  type = "regression",
##                  recursive = FALSE,
##                  ice = TRUE,
##                  inv.link = exp) |>
##   autoplot()

newman_gbm |>
    pdp::partial(pred.var='Octocoral_max_height',
                 n.trees = best.iter,
                 type = "regression",
                 recursive = FALSE,
                 ice = TRUE,
                 inv.link = exp) |>
  autoplot()


## ----bootstrapping, results='markdown', eval=TRUE, mhidden=TRUE, fig.width = 7, fig.height = 5----

## define a function that resamples the data
resample_data <- function(dat) {
  dat <- dat |> sample_n(size = n(), replace = TRUE)
  return(dat)
}

fit_tree <- function(dat) {
  mod_gbm <- gbm(fish_rich ~ fReef_complexity +
                      No_tall_corals +
                      No_corals +
                      Sponge_max_height +
                      Octocoral_max_height +
                      Slope_angle,
                      ## fCountry,
                      data = dat,
                      distribution = "poisson",
                      n.trees = 10000,
                      var.monotone = c(0, 1, 1, 1, 0, 1),
                      interaction.depth = 5,
                      shrinkage = 0.001,
                      bag.fraction = 0.5,
                      train.fraction = 1,
                      n.minobsinnode = 10,
                      verbose = FALSE,
                      cv.folds = 10,
                      n.cores = 1
  )
  return(mod_gbm)
}

## newdata <- with(newman,
##   expand.grid(
##     Octocoral_max_height = seq(min(Octocoral_max_height), max(Octocoral_max_height), len = 100),
##     fCountry = levels(fCountry),
##     fReef_complexity =  NA,
##     No_tall_corals = NA,
##     No_corals = NA,
##     Sponge_max_height =  NA,
##     Slope_angle = NA)
## ) |>
##   mutate(fCountry = factor(fCountry))

newdata <- with(
  newman,
  expand.grid(
    fReef_complexity = levels(fReef_complexity),
    Octocoral_max_height = seq(min(Octocoral_max_height), max(Octocoral_max_height), len = 100),
    ## fCountry = levels(fCountry),
    No_tall_corals = NA,
    No_corals = NA,
    Sponge_max_height = NA,
    Slope_angle = NA
  )
) |>
  mutate(
    ## fCountry = factor(fCountry),
    fReef_complexity = factor(fReef_complexity))

nBoot <- 10
pred.list <- vector('list', nBoot)
rel.inf.list <- vector('list', nBoot)

for (i in 1:nBoot) {
  print(paste0('Boot number: ', i))
  ## Resample the data
  dat <- resample_data(newman)
  ## Fit the tree
  mod <- fit_tree(dat)
  ## Determine the best number of trees
  best.iter <- gbm.perf(mod, method = 'cv')
  ## predict based on shell weight
  fit <- predict(mod, newdata = newdata, n.trees = best.iter)
  pred.list[[i]] <- data.frame(newdata, Boot = i, Fit = fit)
  ## relative influence
  rel.inf.list[[i]] <- summary(mod, n.trees = best.iter)
}

newman.fit <- do.call('rbind', pred.list)
## newman.fit <- newman.fit |>
##     group_by(Octocoral_max_height, fCountry) |>
##     ggdist::median_hdci(Fit)
newman.fit <- newman.fit |>
    group_by(fReef_complexity, Octocoral_max_height) |>
    tidybayes::median_hdci(Fit)

## g1 <-
##   newman.fit |> ggplot(aes(y=Fit, x = Octocoral_max_height,
##   fill = fCountry, colour = fCountry)) +
##   geom_ribbon(aes(ymin=.lower, ymax=.upper), alpha=0.3, color=NA) +
##   geom_line() +
##   scale_fill_viridis_d() +
##   scale_colour_viridis_d() +
##   scale_x_log10() +
##   theme_classic()
g1 <-
  newman.fit |>
  ggplot(aes(y=Fit, x = Octocoral_max_height, color = fReef_complexity)) +
  ## fill = fCountry, colour = fCountry)) +
  geom_ribbon(aes(ymin =.lower, ymax =.upper, fill = fReef_complexity),
    colour = NA, alpha = 0.3) +
  geom_line() +
  scale_colour_viridis_d("Reef\ncomplexity") +
  scale_fill_viridis_d("Reef\ncomplexity") +
  theme_classic()

rel.inf<- do.call('rbind', rel.inf.list)
rel.inf <- rel.inf |>
    group_by(var) |>
    tidybayes::median_hdci(rel.inf)

g2 <-
  rel.inf |>
  arrange(rel.inf) |>
  mutate(var =  factor(var, levels = unique(var))) |>
  ggplot(aes(y=var, x=rel.inf)) +
    geom_vline(xintercept=12.5, linetype='dashed') +
    geom_pointrange(aes(xmin=.lower, xmax=.upper)) +
    theme_classic()

g2 + patchwork::inset_element(g1, left=0.4, bottom=0.01, right=1, top=0.4)
g1 + patchwork::inset_element(g2, left=0.01, bottom=0.01, right=0.5, top=0.35)


## -----------------------------------------------------------------------------
library(randomForest)
newman.rf <- randomForest(
  fish_rich ~ fReef_complexity +
                      No_tall_corals +
                      No_corals +
                      Sponge_max_height +
                      Octocoral_max_height +
                      Slope_angle,
                      ## fCountry,
                       data=newman, importance=TRUE,
                       ntree=1000)
newman.imp <- randomForest::importance(newman.rf)


## Rank by either:
## *MSE (mean decrease in accuracy)
## For each tree, calculate OOB prediction error.
## This also done after permuting predictors.
## Then average diff of prediction errors for each tree
## *NodePurity (mean decrease in node impurity)
## Measure of the total decline of impurity due to each
## predictor averaged over trees
100*newman.imp/sum(newman.imp)
varImpPlot(newman.rf)
## use brute force
newman.rf |>
    pdp::partial("No_tall_corals") |>
    autoplot()

