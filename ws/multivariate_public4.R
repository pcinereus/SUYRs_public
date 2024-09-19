## ---- hierachical
brink <- read_csv(file='../data/brink.csv', trim_ws=TRUE)
brink <- brink %>% mutate(WEEK = factor(WEEK),
                          TREATMENT = factor(TREATMENT),
                          DITCH = factor(DITCH))
## Isolate just the invertegrate data
inverts <- brink %>% dplyr::select(-WEEK, -TREATMENT, -DITCH)


inverts.rda <-  rda(wisconsin(inverts^0.25) ~ TREATMENT*WEEK + Condition(DITCH), data=brink)
inverts.rda <-  rda(wisconsin(inverts^0.25) ~ TREATMENT*WEEK + Condition(WEEK), data=brink)
summary(inverts.rda, display=NULL)
inverts.rda %>% autoplot(geom = 'text')
inverts.rda %>% autoplot()
anova(inverts.rda)
anova(inverts.rda, by='terms')


aa = adonis2(inverts~ DITCH+TREATMENT*WEEK, data=brink, strata = brink$DITCH)
aa
## ----end

