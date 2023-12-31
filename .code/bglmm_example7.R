## ----setup, include=FALSE, warnings=FALSE, message=FALSE----------------------
knitr::opts_chunk$set(echo = TRUE)


## ----libraries, results='markdown', eval=TRUE, message=FALSE, warning=FALSE----
library(car)       #for regression diagnostics
library(broom)     #for tidy output
library(ggfortify) #for model diagnostics
library(sjPlot)    #for outputs
library(knitr)     #for kable
library(effects)   #for partial effects plots
library(emmeans)   #for estimating marginal means
library(MASS)      #for glm.nb
library(MuMIn)     #for AICc
library(tidyverse) #for data wrangling
library(brms)
library(tidybayes)
library(broom.mixed)
library(rstan)
library(patchwork)
library(DHARMa)


## ----readData, results='markdown', eval=TRUE----------------------------------
mullens = read_csv('../data/mullens.csv', trim_ws=TRUE)
glimpse(mullens)


## ----dataPreparation, results='markdown', eval=TRUE, hidden=FALSE-------------
mullens = mullens %>%
  mutate(BREATH=factor(BREATH),
         TOAD=factor(TOAD),
         pBUC=FREQBUC/100,
         pzBUC=ifelse(pBUC==0,0.01,pBUC))


## ----eda1a, results='markdown', eval=TRUE, hidden=FALSE-----------------------
ggplot(mullens,aes(y=FREQBUC, x=factor(O2LEVEL), color=BREATH)) +
    geom_boxplot()


## ----eda1b, results='markdown', eval=TRUE, hidden=FALSE-----------------------
ggplot(mullens,aes(y=pzBUC, x=O2LEVEL, color=BREATH)) +
    geom_smooth() + geom_point()



## ----eda1c, results='markdown', eval=TRUE, hidden=FALSE, fig.width=7, fig.height=7----
ggplot(mullens,aes(y=pzBUC, x=O2LEVEL, color=BREATH)) +
  geom_smooth() + geom_point() +
  facet_wrap(~BREATH+TOAD, scales='free')
  #facet_grid(TOAD~BREATH)


## ----fitModel3a, results='markdown', eval=FALSE, hidden=FALSE-----------------
## mullens %>% group_by(BREATH) %>%
##     summarise(logit(median(pBUC)),
##               logit(mad(pBUC)))
## 
## standist::visualize('gamma(0.01, 0.01)', xlim=c(0,1))
## standist::visualize('beta(1,1)', xlim=c(0,1))
## priors <- prior(normal(0, 5), class='Intercept') +
##     prior(normal(0,5), class='b') +
##     prior(gamma(2,1), class='sd') +
##     prior(gamma(0.01, 0.01), class='phi') +
##     prior(beta(1,1), class='zoi')
## 
## mullens.form <- bf(pBUC ~ BREATH*poly(O2LEVEL,3) + (1|TOAD),
##                    family=zero_one_inflated_beta())
## 
## mullens.brm <- brm(mullens.form,
##                    data=mullens,
##                    prior = priors,
##                    sample_prior = 'yes',
##                    iter=5000, warmup=2500,
##                    thin=5, chains=3, cores=3)
## ##prior_summary(mullens.brm)
## 
## 
## pars <- mullens.brm %>% get_variables()
## wch <- grepl('^b_.*|^sd_.*|phi', pars, perl=TRUE)
## 
## g <- vector('list', length=sum(wch)-1)
## names(g) <- pars[wch][-1]
## for (i in pars[wch]) {
##     print(i)
##     if (i == 'b_Intercept') next
##     p <- mullens.brm %>% hypothesis(paste0(i,'=0'), class='') %>% plot()
##     g[[i]] <- p[[1]]
## }
## patchwork::wrap_plots(g)
## 
## stan_trace(mullens.brm$fit, pars = pars[wch])
## stan_ac(mullens.brm$fit, pars = pars[wch])
## stan_rhat(mullens.brm$fit, pars = pars[wch])
## stan_rhat(mullens.brm$fit)
## stan_ess(mullens.brm$fit)
## 
## 
## preds <- posterior_predict(mullens.brm,  nsamples=250,  summary=FALSE)
## mullens.resids <- createDHARMa(simulatedResponse = t(preds),
##                             observedResponse = mullens$pBUC,
##                             fittedPredictedResponse = apply(preds, 2, median),
##                             integerResponse = FALSE)
## plot(mullens.resids)
## testDispersion(mullens.resids)
## 
## 
## priors <- prior(normal(-1.9, 2), class='Intercept') +
##     prior(normal(0,2), class='b') +
##     prior(cauchy(0,5), class='sd') +
##     prior(normal(0, 5), class='Intercept', dpar='phi') +
##     prior(logistic(0,1), class='Intercept', dpar='coi') +
##     prior(logistic(0,1), class='Intercept', dpar='zoi') +
##     prior(cauchy(0,5), class='sd', dpar='coi') +
##     prior(cauchy(0,5), class='sd', dpar='phi') +
##     prior(cauchy(0,5), class='sd', dpar='zoi')
## 
## prior(gamma(2,1), class='sd', dpar='coi') +
##     prior(gamma(2,1), class='sd', dpar='phi') +
##     prior(gamma(2,1), class='sd', dpar='zoi')
## 
## 
## mullens.form <- bf(pBUC ~ BREATH*poly(O2LEVEL,3) + (poly(O2LEVEL,3)|TOAD),
##                    phi ~ (1|TOAD),
##                    coi ~ (1|TOAD),
##                    zoi ~ (1|TOAD),#BREATH*poly(O2LEVEL,3),
##                    family=zero_one_inflated_beta())
## 
## mullens.brm <- brm(mullens.form,
##                    data=mullens,
##                    prior = priors,
##                    sample_prior = 'yes',
##                    iter=5000, warmup=2500,
##                    thin=5, chains=3, cores=3,
##                    control=list(adapt_delta=0.99))
## 
## prior_summary(mullens.brm)
## mullens.brm %>% conditional_effects() %>% plot(points=TRUE)
## 
## g <- mullens.brm %>%
##     conditional_effects() %>%
##     plot(points=TRUE, ask=FALSE, plot=FALSE)
## patchwork::wrap_plots(g)
## 
## pars <- mullens.brm %>% get_variables()
## pars
## wch <- grepl('^b_.*|^sd_.*', pars, perl=TRUE)
## 
## g <- vector('list', length=sum(wch))
## names(g) <- pars[wch]
## for (i in pars[wch]) {
##     print(i)
##     p <- mullens.brm %>% hypothesis(paste0(i,'=0'), class='') %>% plot()
##     g[[i]] <- p[[1]]
## }
## patchwork::wrap_plots(g)
## 
## stan_trace(mullens.brm$fit, pars = pars[wch])
## stan_ac(mullens.brm$fit, pars = pars[wch])
## stan_rhat(mullens.brm$fit, pars = pars[wch])
## stan_rhat(mullens.brm$fit)
## stan_ess(mullens.brm$fit)
## 
## 
## summary(mullens.brm)
## 
## 
## mullens.brm %>% emtrends(specs='BREATH', var='O2LEVEL', max.degree=3) %>%
##     gather_emmeans_draws() %>%
##     median_hdci()
## 
## mullens.grid <- with(mullens, list(O2LEVEL = modelr::seq_range(O2LEVEL, n=10000)))
## mullens.em <- mullens.brm %>%
##     emmeans(~O2LEVEL|BREATH, at=mullens.grid, type='response') %>%
##     gather_emmeans_draws() %>%
##     mutate(Fit=exp(.value))
## 
## mullens.em %>%
##     median_hdci(Fit) %>%
##     ggplot(aes(y=Fit, x=O2LEVEL, color=BREATH)) +
##     geom_line() +
##     geom_ribbon(aes(ymin=.lower, ymax=.upper, fill=BREATH), color=NA, alpha=0.3)
## 
## 
## mullens.em %>%
##     ungroup() %>%
##     group_by(BREATH, .draw) %>%
##     summarise(Max=O2LEVEL[which.max(.value)]) %>%
##     ungroup() %>%
##     group_by(BREATH) %>%
##     median_hdci(Max)


## ----fitModel3aa, results='markdown', eval=FALSE, hidden=FALSE----------------
## mullens.glmmTMB2a = glmmTMB(pzBUC ~ BREATH*poly(O2LEVEL, 3) + (1|TOAD), data=mullens,
##                             family=beta_family(link = "logit"),  REML=FALSE)
## mullens.glmmTMB2b = update(mullens.glmmTMB2a,  .~.-BREATH:poly(O2LEVEL,  3))
## AICc(mullens.glmmTMB2a,  mullens.glmmTMB2b)


## ----fitModel3b, results='markdown', eval=FALSE, hidden=FALSE-----------------
## mullens.glmmTMB2a = update(mullens.glmmTMB2a,  REML=TRUE)
## mullens.glmmTMB2b = update(mullens.glmmTMB2a,  ~ . - (1|TOAD) + (poly(O2LEVEL, 3)|TOAD))
## AICc(mullens.glmmTMB2a,  mullens.glmmTMB2b)


## ----modelValidation1a, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## plot_model(mullens.lmer1a, type='diag')[-2] %>% plot_grid


## ----modelValidation1b, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## performance::check_model(mullens.lmer1a)


## ----modelValidation1c, results='markdown', eval=FALSE, hidden=FALSE, fig.width=8, fig.height=4----
## mullens.resid <- simulateResiduals(mullens.lmer1a,  plot=TRUE)


## ----modelValidation2a, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## plot_model(mullens.glmmTMB1a, type='diag')[-2] %>% plot_grid


## ----modelValidation2b, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## performance::check_model(mullens.glmmTMB1a)


## ----modelValidation2c, results='markdown', eval=FALSE, hidden=FALSE, fig.width=8, fig.height=4----
## mullens.resid <- simulateResiduals(mullens.glmmTMB1a,  plot=TRUE)


## ----modelValidation3a, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## plot_model(mullens.glmmTMB2b, type='diag')[-2] %>% plot_grid


## ----modelValidation3b, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## performance::check_model(mullens.glmmTMB2b)


## ----modelValidation3c, results='markdown', eval=FALSE, hidden=FALSE, fig.width=8, fig.height=4----
## mullens.resid <- simulateResiduals(mullens.glmmTMB2a,  plot=TRUE)
## mullens.resid <- simulateResiduals(mullens.glmmTMB2b,  plot=TRUE)


## ----modelValidation3cdd, results='markdown', eval=FALSE, hidden=FALSE, fig.width=8, fig.height=4----
## ggemmeans(mullens.glmmTMB2b, ~O2LEVEL|BREATH) %>% plot(add.data=TRUE)


## ----summary1a, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## summary(mullens.lmer1a)


## ----summary1a1, results='markdown', eval=FALSE, hidden=FALSE, echo=FALSE-----
## mullens.tidy <- tidy(mullens.lmer1a)


## ----summary1b, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## tidy(mullens.lmer1a,  conf.int=TRUE)
## tidy(mullens.lmer1a,  conf.int=TRUE) %>% kable


## ----summary1c, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## # warning this is only appropriate for html output
## sjPlot::tab_model(mullens.lmer1a, show.se=TRUE, show.aic=TRUE)


## ----summary1d, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## # with Satterthwaite degrees of freedom calculations
## anova(mullens.lmer1a)
## #OR with Keyward-Roger degrees of freedom calculations
## anova(mullens.lmer1a,  ddf='Kenward-Roger')


## ----summary2a, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## summary(mullens.glmmTMB1a)


## ----summary2a1, results='markdown', eval=FALSE, hidden=FALSE, echo=FALSE-----
## mullens.tidy <- tidy(mullens.glmmTMB1a)


## ----summary2b, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## tidy(mullens.glmmTMB1a,  conf.int=TRUE)
## tidy(mullens.glmmTMB1a,  conf.int=TRUE) %>% kable


## ----summary2c, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## # warning this is only appropriate for html output
## sjPlot::tab_model(mullens.glmmTMB1a, show.se=TRUE, show.aic=TRUE)


## ----summary3a, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## summary(mullens.glmmTMB2b)


## ----summary3a1, results='markdown', eval=FALSE, hidden=FALSE, echo=FALSE-----
## mullens.tidy <- tidy(mullens.glmmTMB2b)


## ----summary3b, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## tidy(mullens.glmmTMB2b,  conf.int=TRUE)
## tidy(mullens.glmmTMB2b,  conf.int=TRUE, exponentiate=TRUE)
## tidy(mullens.glmmTMB2b,  conf.int=TRUE, exponentiate=TRUE) %>% kable


## ----summary3c, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## # warning this is only appropriate for html output
## sjPlot::tab_model(mullens.glmmTMB2b, show.se=TRUE, show.aic=TRUE)


## ----predictions1a, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## emtrends(mullens.lmer1a,  specs='BREATH',  var='O2LEVEL',  max.degree=3,  infer=c(TRUE, TRUE))
## 
## newdata <- with(mullens,  list(O2LEVEL=seq(min(O2LEVEL),  max(O2LEVEL),  len=100),
##                                BREATH=levels(BREATH)))
## mullens.grid <- emmeans(mullens.lmer1a,  ~O2LEVEL|BREATH,  at=newdata) %>% as.data.frame
## mullens.grid %>% group_by(BREATH) %>%
##   summarise(value = O2LEVEL[which.max(emmean)])
## emtrends(mullens.glmmTMB2b,  specs='BREATH',  var='O2LEVEL',
##          max.degree=3,  infer=c(TRUE, TRUE))
## 
## newdata <- with(mullens,  list(O2LEVEL=seq(min(O2LEVEL),  max(O2LEVEL),  len=100),
##                                BREATH=levels(BREATH)))
## mullens.grid <- emmeans(mullens.lmer1a,  ~O2LEVEL|BREATH,  at=newdata) %>% as.data.frame
## mullens.grid %>% group_by(BREATH) %>%
##   summarise(value = O2LEVEL[which.max(emmean)])


## ----predictions2a, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## emtrends(mullens.glmmTMB1a,  specs='BREATH',  var='O2LEVEL',  max.degree=3,  infer=c(TRUE, TRUE))
## 
## newdata <- with(mullens,  list(O2LEVEL=modelr::seq_range(O2LEVEL, n=1000),
##                                BREATH=levels(BREATH)))
## mullens.grid <- emmeans(mullens.glmmTMB1a,  ~O2LEVEL|BREATH,  at=newdata) %>% as.data.frame
## mullens.grid %>% group_by(BREATH) %>%
##   summarise(value = O2LEVEL[which.max(emmean)])


## ----predictions3a, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## emtrends(mullens.glmmTMB2b,  specs='BREATH',  var='O2LEVEL',  max.degree=3,  infer=c(TRUE, TRUE))
## 
## ## newdata <- with(mullens,  list(O2LEVEL=seq(min(O2LEVEL),  max(O2LEVEL),  len=100),
## ##                                BREATH=levels(BREATH)))
## mullens.grid <- with(mullens,  list(O2LEVEL=modelr::seq_range(O2LEVEL, n=1000),
##                                BREATH=levels(BREATH)))
## newdata <- emmeans(mullens.glmmTMB2b,  ~O2LEVEL|BREATH,  at=mullens.grid) %>% as.data.frame
## newdata %>% group_by(BREATH) %>%
##   summarise(value = O2LEVEL[which.max(emmean)])
## ## r.squaredGLMM(mullens.glmmTMB2a)
## performance::r2_nakagawa(mullens.glmmTMB2a)


## ----summaryFigures1a, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## mullens.grid = with(mullens,
##    list(BREATH=levels(BREATH),
##      O2LEVEL=seq(min(O2LEVEL), max(O2LEVEL), len=100)
##    )
## )
## newdata = emmeans(mullens.lmer1a, ~O2LEVEL|BREATH,
##                   at=mullens.grid, type='response') %>% as.data.frame %>%
##           mutate(emmean = emmean^2,
##                  lower.CL=lower.CL^2,
##                  upper.CL=upper.CL^2)
## head(newdata)
## 
## ggplot() +
##     geom_ribbon(data=newdata,
##                 aes(ymin=lower.CL,ymax=upper.CL,
##                     x=O2LEVEL, fill=BREATH), alpha=0.3)+
##     geom_line(data=newdata,
##               aes(y=emmean, x=O2LEVEL, color=BREATH)) +
##     theme_classic()


## ----summaryFigures2a, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## mullens.grid = with(mullens,
##    list(BREATH=levels(BREATH),
##      O2LEVEL=seq(min(O2LEVEL), max(O2LEVEL), len=100)
##    )
## )
## newdata = emmeans(mullens.glmmTMB1a, ~O2LEVEL|BREATH,
##                   at=mullens.grid, type='response') %>% as.data.frame %>%
##           mutate(emmean = emmean^2,
##                  lower.CL=lower.CL^2,
##                  upper.CL=upper.CL^2)
## head(newdata)
## 
## ggplot() +
##     geom_ribbon(data=newdata,
##                 aes(ymin=lower.CL,ymax=upper.CL,
##                     x=O2LEVEL, fill=BREATH), alpha=0.3)+
##     geom_line(data=newdata,
##               aes(y=emmean, x=O2LEVEL, color=BREATH)) +
##     theme_classic()


## ----summaryFigures3a, results='markdown', eval=FALSE, hidden=FALSE, fig.width=7, fig.height=7----
## mullens.grid = with(mullens,
##    list(BREATH=levels(BREATH),
##      O2LEVEL=modelr::seq_range(O2LEVEL, n=1000)
##    )
## )
## newdata = emmeans(mullens.glmmTMB2b, ~O2LEVEL|BREATH,
##                   at=mullens.grid, type='response') %>% as.data.frame
## head(newdata)
## 
## ggplot() +
##     geom_ribbon(data=newdata,
##                 aes(ymin=lower.CL,ymax=upper.CL,
##                     x=O2LEVEL, fill=BREATH), alpha=0.3)+
##     geom_line(data=newdata,
##               aes(y=response, x=O2LEVEL, color=BREATH)) +
##     theme_classic()
## 
## mullens <- mullens %>%
##     mutate(resid=resid(mullens.glmmTMB2b),
##            fitted=fitted(mullens.glmmTMB2b),
##            obs=fitted+resid
##            )
## 
## ggplot() +
##     geom_point(data=mullens, aes(y=obs, x=O2LEVEL, color=BREATH)) +
##     geom_ribbon(data=newdata,
##                 aes(ymin=lower.CL,ymax=upper.CL,
##                     x=O2LEVEL, fill=BREATH), alpha=0.3)+
##     geom_line(data=newdata,
##               aes(y=response, x=O2LEVEL, color=BREATH)) +
##     theme_classic()


## ----name, results='markdown', eval=FALSE, hidden=FALSE-----------------------
## mullens.glmmTMB = glmmTMB(pzBUC ~ BREATH*poly(O2LEVEL, 3) + (1|TOAD), data=mullens,
##                           family=beta_family(link = "logit"))
## 
## #mullens.glmmTMB1 = glmmTMB(pzBUC ~ BREATH*poly(scale(O2LEVEL), 3) +
## #                               (poly(scale(O2LEVEL), 3)|TOAD), data=mullens,
## #                          family=beta_family(link = "logit"))
## 
## mullens.glmmTMB1 = glmmTMB(pzBUC ~ BREATH*poly(O2LEVEL, 3) + (O2LEVEL|TOAD), data=mullens,
##                            family=beta_family(link = "logit"))
## AIC(mullens.glmmTMB, mullens.glmmTMB1)
## 
## 
## ## library(lme4)
## ## mullens.glmer = glmer(pBUC ~ BREATH*poly(O2LEVEL, 3) + (1|TOAD),
## ##                       data=mullens, family=Beta(link='logit'))
## ## mullens.glmer1 = glmer(pBUC ~ BREATH*poly(O2LEVEL, 3) + (BREATH|TOAD),
## ##                       data=mullens, family=binomial(link='logit'))
## ## AIC(mullens.glmer,mullens.glmer1)
## ##                                         # Try this trick instead
## ## mullens.glmer = glmer(cbind(FREQBUC,100-FREQBUC) ~ BREATH*poly(O2LEVEL, 3) + (1|TOAD),
## ##                       data=mullens, family=binomial(link='logit'))
## 
## ## mullens.glmer1 = glmer(cbind(FREQBUC,100-FREQBUC) ~ BREATH+poly(O2LEVEL, 3) + (1|TOAD),
## ##                       data=mullens, family=binomial(link='logit'))
## ## AICc(mullens.glmer, mullens.glmer1)
## 
## #ggplot() +
## #    geom_point(data=NULL, aes(y=resid(mullens.glmmTMB), x=fitted(mullens.glmmTMB)))
## 
## 
## #plot(mullens.glmmTMB)
## 
## 
## plot_model(mullens.glmmTMB1, type='diag')
## performance::check_model(mullens.glmmTMB1)
## mullens.resid = simulateResiduals(mullens.glmmTMB1,  plot=TRUE)
## 
## plot(allEffects(mullens.glmmTMB1))
## plot(allEffects(mullens.glmmTMB1),  multiline=TRUE,  ci.style='bands')
## plot_model(mullens.glmmTMB1, type='eff', terms=c('O2LEVEL', 'BREATH'))
## 
## 
## summary(mullens.glmmTMB1)
## tidy(mullens.glmmTMB1, conf.int=TRUE,  exponentiate=TRUE)
## 
## #model.matrix(~BREATH*poly(O2LEVEL,3), mullens)
## emtrends(mullens.glmmTMB1,  ~BREATH,  var='O2LEVEL',  max.degree=3)
## 
## emmeans(mullens.glmmTMB1, ~O2LEVEL|BREATH) #%>% as.data.frame
## 
## emmeans(mullens.glmmTMB1, ~O2LEVEL|BREATH,type='response')
## #contrast(emmeans(mullens.glmmTMB, ~O2LEVEL|BREATH,type='response'),interaction='poly')
## 
## mullens.grid = with(mullens,
##    list(BREATH=levels(BREATH),
##      O2LEVEL=seq(min(O2LEVEL), max(O2LEVEL), len=100)
##    )
## )
## newdata = emmeans(mullens.glmmTMB1, ~O2LEVEL|BREATH,
##         at=mullens.grid, type='response') %>% as.data.frame
## head(newdata)
## 
## ggplot() +
##     geom_ribbon(data=newdata,
##                 aes(ymin=lower.CL,ymax=upper.CL,
##                     x=O2LEVEL, fill=BREATH), alpha=0.3)+
##     geom_line(data=newdata,
##               aes(y=response, x=O2LEVEL, color=BREATH)) +
##     theme_classic()
## 
## r.squaredGLMM(mullens.glmmTMB1)
## performance::r2_nakagawa(mullens.glmmTMB1)
## 
## 

