## ----setup, include=FALSE, warnings=FALSE, message=FALSE----------------------
knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE,cache.lazy = FALSE, tidy='styler')


## ----libraries, results='markdown', eval=TRUE, message=FALSE, warning=FALSE----
library(car)       #for regression diagnostics
library(broom)     #for tidy output
library(ggfortify) #for model diagnostics
library(sjPlot)    #for outputs
library(knitr)     #for kable
library(effects)   #for partial effects plots
library(emmeans)   #for estimating marginal means
library(ggeffects) #for plotting marginal means
library(MASS)      #for glm.nb
library(MuMIn)     #for AICc
library(tidyverse) #for data wrangling
library(modelr)    #for auxillary modelling functions
library(DHARMa)    #for residual diagnostics plots
library(performance) #for residuals diagnostics
library(see)         #for plotting residuals
library(patchwork) #grid of plots
library(scales)    #for more scales


## ----readData, results='markdown', eval=TRUE----------------------------------
quinn = read_csv('../data/quinn.csv', trim_ws=TRUE)
glimpse(quinn)
summary(quinn)


## ----dataprep, results='markdown', eval=TRUE, hidden=TRUE---------------------
quinn = quinn %>% mutate(SEASON = factor(SEASON, levels=c('Spring', 'Summer', 'Autumn', 'Winter')),
                         DENSITY = factor(DENSITY))


## ----eda, hidden=TRUE---------------------------------------------------------
ggplot(quinn, aes(y=RECRUITS, x=SEASON, fill=DENSITY)) +
     geom_boxplot()


## ----eda1, hidden=TRUE--------------------------------------------------------
ggplot(quinn, aes(y=RECRUITS, x=SEASON, fill=DENSITY)) +
  geom_boxplot() +
  scale_y_log10()


## ----fitModel1a, hidden=TRUE--------------------------------------------------
quinn.glmG <- glm(log(RECRUITS+1) ~ DENSITY*SEASON, data=quinn, family=gaussian)


## ----fitModel1b, hidden=TRUE--------------------------------------------------
quinn.glm <- glm(RECRUITS ~ DENSITY*SEASON, data=quinn,
                  family=poisson(link='log'))


## ----ValidateModel1a, results='markdown', eval=TRUE, fig.width=7, fig.height=7, hidden=TRUE----
quinn.glm %>% autoplot(which=1:6)


## ----ValidateModel1b, results='markdown', eval=TRUE, fig.width=7, fig.height=7, hidden=TRUE----
quinn.glm %>% performance::check_model()
quinn.glm %>% performance::check_overdispersion()
quinn.glm %>% performance::check_zeroinflation()


## ----ValidateModel1c, results='markdown', eval=TRUE, fig.width=8, fig.height=5, hidden=TRUE----
quinn.resid <- quinn.glm %>% simulateResiduals(plot=TRUE)
quinn.resid %>% testResiduals()
quinn.resid %>% testDispersion()
quinn.resid %>% testZeroInflation()
#testTemporalAutocorrelation(quinn.glm1)


## ----modelValidation1d, results='markdown', eval=TRUE, hidden=TRUE------------
## goodness of fit
1-pchisq(quinn.glm$deviance, df=quinn.glm$df.residual)
## any evidence of overdispersion
quinn.glm$deviance/quinn.glm$df.residual


## ----modelValidation3, results='markdown', eval=TRUE, hidden=TRUE-------------
quinn %>% group_by(SEASON, DENSITY) %>%
  summarise(Zeros = sum(RECRUITS==0), 
            Prop = Zeros/n(),
            Mean = mean(RECRUITS))
x <- rpois(100000,lambda = 2.67)
tab.1 <- table(x == 0)
tab.1/sum(tab.1)

## OR,  over the entire data
## is this due to excessive zeros (zero-inflation)
tab <- table(quinn$RECRUITS == 0)
tab/sum(tab)
## 5% is not many.. so it cant be zero-inflated
## how many 0's would we expect from a poisson distribution with a mean similar to our mean
mean(quinn$RECRUITS)
x <- rpois(100000, lambda = mean(quinn$RECRUITS))
tab.1 <- table(x == 0)
tab.1/sum(tab.1)


## ----fitModel2, results='markdown', eval=TRUE, hidden=TRUE--------------------
quinn.nb <- MASS::glm.nb(RECRUITS ~ DENSITY*SEASON, data=quinn)


## ----modelValidation4a, results='markdown', eval=TRUE, fig.width=7, fig.height=7, hidden=TRUE----
quinn.nb %>% autoplot(which=1:6)


## ----modelValidation4b, results='markdown', eval=TRUE, fig.width=7, fig.height=7, hidden=TRUE----
quinn.nb %>% performance::check_model()


## ----modelValidation4c, results='markdown', eval=TRUE, fig.width=7, fig.height=7, hidden=TRUE----
quinn.resid <- quinn.nb %>% simulateResiduals(plot=TRUE)
quinn.resid %>% testResiduals()
quinn.resid %>% testDispersion()
quinn.resid %>% testZeroInflation()


## ----modelValidation5, results='markdown', eval=TRUE, hidden=TRUE-------------
## goodness of fit
1-pchisq(quinn.nb$deviance, df=quinn.nb$df.residual)
## any evidence of overdispersion
quinn.nb$deviance/quinn.nb$df.residual


## ----modelValidation6, results='markdown', eval=TRUE, hidden=TRUE-------------
AICc(quinn.glm, quinn.nb)


## ----partialplots1a, results='markdown', eval=TRUE, hidden=TRUE---------------
quinn.nb %>% plot_model(type='eff',  terms=c('SEASON', 'DENSITY'))


## ----partialplots1b, results='markdown', eval=TRUE, hidden=TRUE---------------
quinn.nb %>% allEffects() %>% plot(multiline=TRUE, ci.style='bar')
quinn.nb %>% allEffects() %>% plot(multiline=TRUE, ci.style='bar', type='link')


## ----partialplots1c, results='markdown', eval=TRUE, hidden=TRUE---------------
quinn.nb %>% ggpredict(c('SEASON', 'DENSITY')) %>% plot()


## ----partialplots1d, results='markdown', eval=TRUE, hidden=TRUE---------------
quinn.nb %>% ggemmeans(~SEASON*DENSITY) %>% plot()


## ----summarys1a, results='markdown', eval=TRUE, hidden=TRUE-------------------
quinn.nb %>% summary()


## ----summarys1b, results='markdown', eval=TRUE, hidden=TRUE-------------------
quinn.nb %>% tidy(conf.int=TRUE)


## ----summarys1c, results='markdown', eval=TRUE, hidden=TRUE-------------------
quinn.nb %>% tidy(conf.int=TRUE, exponentiate=TRUE)


## ----summarys1d, results='markdown', eval=TRUE, hidden=TRUE-------------------
quinn.glm %>% tidy(conf.int=TRUE, exponentiate=TRUE)


## ----mainEffects1a, results='markdown', eval=TRUE, hidden=TRUE----------------
quinn.nb %>% emmeans(~DENSITY|SEASON) %>% pairs() %>% summary(infer=TRUE)


## ----mainEffects1a1, results='markdown', eval=TRUE, echo=FALSE, hidden=TRUE----
eff <- quinn.nb %>%
    emmeans(~DENSITY|SEASON, type='link') %>%
    pairs() %>%
    as.data.frame()


## ----mainEffects1b, results='markdown', eval=TRUE, hidden=TRUE----------------
quinn.nb %>% emmeans(~DENSITY|SEASON, type='response') %>% pairs()


## ----mainEffects1b2, results='markdown', eval=TRUE, hidden=TRUE---------------
quinn.nb %>% emmeans(~DENSITY|SEASON, type='response') %>% regrid() %>% pairs()


## ----mainEffects1c, results='markdown', eval=TRUE, hidden=TRUE----------------
quinn.nb %>% emmeans(~DENSITY|SEASON,  type='response') %>% pairs() %>% summary(infer=TRUE)


## ----mainEffects1d, results='markdown', eval=TRUE, hidden=TRUE----------------
eff <- quinn.nb %>%
    emmeans(~DENSITY|SEASON,  type='response') %>%
    pairs() %>%
    summary(infer=TRUE) %>%
    as.data.frame

eff %>%
    ggplot(aes(y=ratio, x=SEASON)) +
    geom_pointrange(aes(ymin=asymp.LCL, ymax=asymp.UCL)) +
    geom_text(aes(label=sprintf("p=%.2f", p.value)), nudge_x=0.25) +
    geom_hline(yintercept=1, linetype='dashed') +
    scale_x_discrete(name='') +
    scale_y_continuous(name='Density effect (High vs Low)', trans=scales::log2_trans(),
                       breaks=scales::breaks_log(base=2)) +
    coord_flip(ylim=c(0.25, 4)) +
    theme_classic()


## ----summaryFig, results='markdown', eval=TRUE, hidden=TRUE-------------------
newdata <- emmeans(quinn.nb, ~DENSITY|SEASON, type = 'response') %>% as.data.frame
head(newdata)
ggplot(newdata, aes(y = response, x = SEASON, fill = DENSITY)) +
    geom_pointrange(aes(ymin = asymp.LCL, ymax = asymp.UCL, shape = DENSITY),
                    position = position_dodge(width = 0.2)) +
    theme_classic() +
    theme(axis.title.x = element_blank(),
          legend.position = c(0.01,1), legend.justification = c(0,1)) +
    annotate(geom = 'text', x = 'Summer', y = 70, label = '*', size = 7) +
    scale_shape_manual(values = c(21, 22))


## ----zeroinflate, results='markdown', eval=TRUE, hidden=TRUE, fig.width=10, fig.height=5----
library(pscl)
library(glmmTMB)
quinn.zip <- zeroinfl(RECRUITS ~ DENSITY*SEASON | 1, data=quinn,  dist='poisson')
quinn.zip <- glmmTMB(RECRUITS ~ DENSITY*SEASON, zi=~1, data=quinn,  family=poisson())
quinn.resid <- quinn.zip %>% simulateResiduals(plot=TRUE)
## The following does not work due to a lack of pearson residuals for glmmTMB
## quinn.zip %>% performance::check_overdispersion()

quinn.zip %>% summary()
#tidy(quinn.zip,  conf.int=TRUE, exponentiate = TRUE)
exp(-3.0037)

## quinn.zip1 <- zeroinfl(RECRUITS ~ DENSITY*SEASON | SEASON, data=quinn,  dist='poisson')
quinn.zip1 <- glmmTMB(RECRUITS ~ DENSITY*SEASON, zi=~SEASON, data=quinn,  family=poisson())
quinn.resid <- quinn.zip1 %>% simulateResiduals(plot=TRUE)

quinn.zip1 %>% summary()
plogis(-20.468)
plogis(-20.468 + 19.214)


## quinn.zinb <- zeroinfl(RECRUITS ~ DENSITY*SEASON | 1, data=quinn,  dist='negbin')
quinn.zinb <- glmmTMB(RECRUITS ~ DENSITY*SEASON, zi=~SEASON, data=quinn,  family=nbinom2())
quinn.resid <- quinn.zinb %>% simulateResiduals(plot=TRUE)
AICc(quinn.zip,  quinn.zip1, quinn.zinb)

quinn.zinb %>% summary()

