## ---- MVABUND spiders
mva = mvabund(spider.abund)
spider.mod <- manyglm(mva~
                          scale(soil.dry)+
                          scale(moss)+
                          scale(herb.layer)+
                          scale(bare.sand),
                      family=poisson(link='log'),
                      data=spider.env)

plot(spider.mod)

spider.mod1 <- manyglm(mva~
                          scale(soil.dry)+
                          scale(moss)+
                          scale(herb.layer)+
                           scale(bare.sand), 
                       family="negative.binomial",
                       data=spider.env)

plot(spider.mod1)
drop1(spider.mod1)
spider.mod1
spider.mod1 |> summary()


anova(spider.mod1)
anova(spider.mod1, test='LR')
anova(spider.mod1, cor.type = 'R')
anova(spider.mod1, cor.type = 'shrink')
anova(spider.mod1, p.uni='adjusted')
summary(spider.mod1, test="LR")
## ----end

