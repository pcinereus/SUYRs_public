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
library(car)
library(dismo)
library(pdp)
library(ggfortify)
library(randomForest)
library(tidyverse)
library(gridExtra)
library(patchwork)
library(easystats)


## ----readData, results='markdown', eval=TRUE----------------------------------
leathwick <- read_csv('../data/leathwick.csv', trim_ws=TRUE)
glimpse(leathwick)


## ----readData1, results='markdown', eval=TRUE---------------------------------
leathwick_test <- read_csv('../data/leathwick_test.csv', trim_ws=TRUE)
glimpse(leathwick_test)


## ----preparation, results='markdown', eval=TRUE, hidden=FALSE-----------------
leathwick <- leathwick |>
  mutate(Method=factor(Method),
    LocSed = as.numeric(LocSed)) |>
  ## LocSed=factor(LocSed)) |>
  as.data.frame()


## ----preparation1, results='markdown', eval=TRUE, hidden=FALSE----------------
leathwick_test <- leathwick_test |>
  mutate(Method=factor(Method),
    LocSed=as.numeric(LocSed)) |>
  as.data.frame()


## ----EDA, results='markdown', eval=TRUE, fig.width=15, fig.height=15----------
## scatterplotMatrix(~Angaus+SegSumT+SegTSeas+SegLowFlow+DSDist+DSMaxSlope+DSDam+
##                     USAvgT+USRainDays+USSlope+USNative+Method+LocSed,  data=leathwick,
##                     diagonal=list(method='boxplot'))
scatterplotMatrix(~Angaus+SegSumT+SegTSeas+SegLowFlow+DSDist+DSMaxSlope+DSDam+
                    USAvgT+USRainDays+USSlope+USNative+LocSed,
  data=leathwick,
  diagonal=list(method='boxplot'))


## ----fitModel1, results='markdown', eval=TRUE, cache=TRUE---------------------
set.seed(123)
## leathwick.gbm <- gbm(Angaus ~ SegSumT+SegTSeas+SegLowFlow+DSDist+DSMaxSlope+DSDam+
##                        USAvgT+USRainDays+USSlope+USNative+Method+LocSed,
##   data=leathwick,
##   distribution='bernoulli',
##   var.monotone=c(1,1,0,-1,-1,0,1,-1,-1,-1,0,-1),
##   n.trees=10000,
##   interaction.depth=5,
##   bag.fraction=0.5,
##   shrinkage=0.001,
##   train.fraction=1,
##   cv.folds=3)
leathwick.gbm <- gbm(Angaus ~ SegSumT+SegTSeas+SegLowFlow+DSDist+DSMaxSlope+DSDam+
                       USAvgT+USRainDays+USSlope+USNative+LocSed,
  data=leathwick,
  distribution = 'bernoulli',
  var.monotone = c(1,1,0,-1,-1,0,1,-1,-1,-1,-1),
  n.trees = 10000,
  interaction.depth = 5,
  bag.fraction = 0.5,
  shrinkage = 0.001,
  train.fraction = 1,
  cv.folds = 3,
  n.cores = 1  #should not need this, it seems to be an wsl issue
)


## ----fitModel2, results='markdown', eval=TRUE, cache=FALSE--------------------
(best.iter = gbm.perf(leathwick.gbm,method='OOB'))
(best.iter = gbm.perf(leathwick.gbm,method='cv'))


## ----fitModel3, results='markdown', eval=TRUE, cache=FALSE, fig.width=10, fig.height=10----
leathwick.gbm |> summary(n.trees = best.iter)


## ----partialPlots1, results='markdown', eval=TRUE, cache=TRUE, fig.width=10, fig.height=10----
leathwick.gbm |>
  pdp::partial(pred.var='SegSumT',
    n.trees=best.iter,
    inv.link=plogis,
    recursive=FALSE,
    ice = TRUE,
    type='regression') |>
  autoplot() +
  ylim(0,1)


## ----partialPlots2, results='markdown', eval=TRUE, cache=TRUE, fig.width=10, fig.height=10----
## Treat simplar to posteriors
fits <-
  leathwick.gbm |>
  pdp::partial(pred.var='SegSumT',
    n.trees=best.iter,
    inv.link=plogis,
    recursive=FALSE,
    ice = TRUE,
    type='regression')
fits |>
  group_by(SegSumT) |>
  tidybayes::median_hdci(yhat) |>
  ggplot(aes(y = yhat, x = SegSumT)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) +
  geom_line() +
  theme_classic()


## ----partialPlots3, results='markdown', eval=TRUE, cache=TRUE, fig.width=10, fig.height=10----
nms <- colnames(leathwick)
partial_list <- vector('list', 11)
names(partial_list) <- nms[c(3:12, 14)]
for (nm in names(partial_list)) {
  print(nm)
  partial_list[[nm]] <- leathwick.gbm |> pdp::partial(pred.var=nm,
    n.trees=best.iter,
    inv.link=plogis,
    recursive=FALSE,
    ice = TRUE,
    type='regression') |>
    autoplot() + ylim(0, 1)
}


wrap_plots(partial_list) & theme_bw()
## do.call('grid.arrange', partial_list)


## ----Accuracy1, results='markdown', eval=TRUE, cache=FALSE, fig.width=10, fig.height=10----
leathwick_test |>
  bind_cols(Pred = predict(leathwick.gbm,newdata=leathwick_test,
    n.tree=best.iter, type='response')) |>
  ggplot() +
  geom_boxplot(aes(y=Pred,  x=as.factor(Angaus_obs))) +
  geom_point(aes(y=Pred, x=as.factor(Angaus_obs)),
    position=position_jitter(width=0.05))


#preds <- predict.gbm(leathwick.gbm, newdata=leathwick_test,
#                     n.trees=best.iter,  type='response')
preds <- leathwick_test |>
  bind_cols(Pred = predict(leathwick.gbm,newdata=leathwick_test,
    n.tree=best.iter, type='response'))
pres <- preds |> filter(Angaus_obs==1) |> pull(Pred)
abs <- preds |> filter(Angaus_obs==0) |> pull(Pred)
e <- dismo::evaluate(p=pres,  a=abs)
e


## ----name, results='markdown', eval=TRUE, fig.width=15, fig.height=15---------
data(Anguilla_grids)
leathwick.grid = Anguilla_grids
glimpse(leathwick.grid)
glimpse(Anguilla_grids)

plot(Anguilla_grids)

## Method <- factor('electric', levels=levels(leathwick$Method))
## Method = as.data.frame(Method)

fit <- predict(Anguilla_grids, leathwick.gbm,
               n.trees=best.iter,  type='response')
## fit <- predict(leathwick.grid, leathwick.gbm,  const=Method,
##                n.trees=best.iter,  type='response')
#fit <- mask(fit,  raster(leathwick.grid, 1))
library(stars)
fit= stars::st_as_stars(fit)


ggplot() +
  geom_stars(data=fit) +
  scale_fill_gradient(low='red', high='green', 'Probability\nof occurrance', na.value=NA) +
  coord_sf(expand=FALSE) +
  theme_bw()

