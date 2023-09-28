
nBoot <- 10
loyn.pred <- with(loyn,
                     expand.grid(lAREA = modelr::seq_range(log(AREA), n = 100),
                                fGRAZE = levels(fGRAZE),
                                DIST = NA,
                                LDIST = NA,
                                ALT = NA,
                                YR.ISOL = NA)
) |>
  mutate(fGRAZE = factor(fGRAZE),
         AREA = exp(lAREA)) |>
  dplyr::select(AREA, lAREA, fGRAZE, DIST, LDIST, ALT, YR.ISOL)

loyn.list <- vector('list', nBoot) 
##loyn.list
loyn.sum <- vector('list', nBoot) 
for (i in 1:nBoot) {
    print(paste0('Boot number: ', i))
    ## Create random set
    loyn.rnd <- loyn %>%
        sample_n(size = n(), replace=TRUE)
    ## Fit the trees
    loyn.gbm = gbm(ABUND ~ AREA + fGRAZE + DIST + LDIST + ALT + YR.ISOL,
                      data=loyn.rnd,
                      distribution='gaussian',
                      var.monotone=c(1,0,1,1,1,1),
                      n.trees=5000,
                      interaction.depth=5,
                      bag.fraction=0.5,
                      shrinkage=0.001,
                      train.fraction=1,
                      n.minobsinnode = 2,
                      cv.folds=3)
    ## Determine the best number of trees
    (best.iter = gbm.perf(loyn.gbm,method='cv'))
    ## predict based on shell weight
    fit <- predict(loyn.gbm, newdata = loyn.pred, n.trees = best.iter) 
    loyn.list[[i]] <- data.frame(loyn.pred, Boot = i, Fit = fit)
    ## relative influence
    loyn.sum[[i]] <- summary(loyn.gbm, n.trees = best.iter)
}
loyn.fit <- do.call('rbind', loyn.list)
loyn.fit <- loyn.fit %>%
    group_by(AREA, fGRAZE) %>%
    ## summarise(Median = median(Fit),
    ##           Lower = quantile(Fit, p=0.025),
    ##           Upper = quantile(Fit, p=0.975))
    ggdist::median_hdci(Fit)       
g1 <- loyn.fit %>% ggplot(aes(y=Fit, x=AREA, fill=fGRAZE, color=fGRAZE)) +
    geom_ribbon(aes(ymin=.lower, ymax=.upper), alpha=0.3, color=NA) +
    geom_line() +
    scale_fill_viridis_d() +
  scale_colour_viridis_d() +
  scale_x_log10() +
    theme_classic()

loyn.inf <- do.call('rbind', loyn.sum)
loyn.inf <- loyn.inf %>%
    group_by(var) %>%
    ggdist::median_hdci(rel.inf)       

g2 <- loyn.inf %>%
  arrange(rel.inf) |> 
  mutate(var =  factor(var, levels = unique(var))) |> 
  ggplot(aes(y=var, x=rel.inf)) +
    geom_vline(xintercept=12.5, linetype='dashed') +
    geom_pointrange(aes(xmin=.lower, xmax=.upper)) +
    theme_classic()

g2 + patchwork::inset_element(g1, left=0.3, bottom=0.01, right=1, top=0.7)
g1 + patchwork::inset_element(g2, left=0.5, bottom=0.01, right=1, top=0.4)
