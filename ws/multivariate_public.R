############################################################
##             RDA (Correlations, Euclidean distance)
##            /  \
##Unconstrained   Constrained ---> ordination
##   (PCA)            (RDA)   ---> anova
##
##              CA (Chisq distance)
##             /  \
##Unconstrained   Constrained ---> ordination
## (CA)             (CCA)     ---> anova
##
##             PCoA (any distance)
##             /  \
##Unconstrained   Constrained ---> ordination
##                            ---> anova
##
##             dbRDA (any distance)
##             /  \
##Unconstrained   Constrained ---> ordination
##                            ---> anova
##
##Unconstrained  ---> ordination
##               ---> envfit (overlay enviromental data) (permutation test)
##               ---> lm/glm etc (response or predictor)
#############################################################
##     Dissimilarity
##            --> MDS      ---> ordination
##            --> bioenv   ---> variable importance       (perm test)
##            --> adonis2* ---> anova                     (perm test)
##            --> simper   ---> similarity percentages
##            --- betadisp ---> homogeneity of dispersion (perm test)
#############################################################
##     Model based ordination
##            ---> glmmTMB (via reduced rank / latent variable)
##            ---> gllvm (generalised latent variable models)
#############################################################
##     Model based
##            ---> manyglm ---> anova
##            ---> gllvm (generalized latent variable models)
#############################################################

## ---- libraries
library(tidyverse)
library(vegan)
library(GGally)
library(corrplot)
library(car)
library(ggvegan)
library(ggrepel)
## ----end


## ---- read Spider
spider.abund <- read_csv(file = "../data/spider.abund.csv", trim_ws = TRUE)
spider.env <- read_csv(file = "../data/spider.env.csv", trim_ws = TRUE)
glimpse(spider.abund)
glimpse(spider.env)
## ----end

## ---- EDA spider
spider.abund |>
  cor() |>
  corrplot(type = 'upper',
    diag = FALSE)
## And now with axes arrange according to first princomp
spider.abund |>
  cor() |>
  corrplot(type = 'upper',
    order = 'FPC',
    diag = FALSE)
## ----end

## ---- EDA 2
spider.abund |>
  ggpairs(lower = list(continuous = "smooth"),
    diag = list(continuous = "density"),
    axisLabels = "show")
## ----end

## PCA ------------------------------------------------------

spider.std <- spider.abund |>
  mutate(across(everything(), ~.x^0.25)) |>
  wisconsin()
#OR
spider.std <- (spider.abund^0.25) |>
  wisconsin()
spider.std

spider.rda <- rda(spider.abund, scale=TRUE)
summary(spider.rda, display=NULL)

screeplot(spider.rda)
abline(a=1,b=0)

## skip
scores(spider.rda, choices=1:3, display='sites')
scores(spider.rda, choices=1:3, display='species')

## Quick and nasty ordination plots
biplot(spider.rda, scaling='species')
biplot(spider.rda, scaling='sites')

## Quick and nasty ordination plots
pl <- vegan::ordiplot(spider.rda)
points(pl, "sites", pch=21, col="red", bg="yellow")
text(pl, "sites", col="red", cex=0.9)
text(pl, "species", col="blue", cex=0.9)

autoplot(spider.rda)
autoplot(spider.rda) + theme_bw()
autoplot(spider.rda,geom='text') + theme_bw()


spider.rda.scores <- spider.rda |>
  fortify()
spider.rda.scores

g <-
  ggplot(data = NULL, aes(y=PC2, x=PC1)) +
  geom_hline(yintercept=0, linetype='dotted') +
  geom_vline(xintercept=0, linetype='dotted') +
  geom_point(data=spider.rda.scores |> filter(score=='sites')) +
  geom_text(data=spider.rda.scores |> filter(score=='sites'),
    aes(label=label), hjust=-0.2) +
  geom_segment(data=spider.rda.scores |> filter(score=='species'),
    aes(y=0, x=0, yend=PC2, xend=PC1),
    arrow=arrow(length=unit(0.3,'lines')), color='red') +
  geom_text_repel(data=spider.rda.scores |> filter(score=='species'),
    aes(y=PC2*1.1, x=PC1*1.1, label=label), color='red') +
  theme_bw()
g

## Nice axes titles
eig <- eigenvals(spider.rda)

g <- g +
  scale_y_continuous(paste(names(eig[2]),
    sprintf('(%0.1f%% explained var.)',
    100 * eig[2]/sum(eig))))+
  scale_x_continuous(paste(names(eig[1]),
    sprintf('(%0.1f%% explained var.)',
    100 * eig[1]/sum(eig))))
g

#put a circle
circle.prob <- 0.95
r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(spider.rda$CA$u[,1:2]^2))^(1/4)
theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
circle <- data.frame(PC1 = r * cos(theta), PC2 = r * sin(theta))
g <- g +
  geom_path(
    data = circle,
    aes(y = PC2, x = PC1),
    color = muted("white"), size = 1 / 2, alpha = 1 / 3
  )
g

## ----end


## ---- Envfit
spider.env |>
  cor() |>
  corrplot(type = 'upper',
    order = 'FPC',
    diag = FALSE)

spider.env |>
  ggpairs(lower = list(continuous = "smooth"),
    diag = list(continuous = "density"),
    axislabels = "show")
spider.envfit <- envfit(spider.rda, env = spider.env)
spider.envfit


spider.env.scores <- spider.envfit |>
  fortify() |>
  mutate(Flag = factor(ifelse(sqrt(PC1^2 + PC2^2) > r, 1, 0)))
g <- g +
  geom_segment(data=spider.env.scores,
    aes(y=0, x=0, yend=PC2, xend=PC1, alpha = Flag, show.legend = FALSE),
    arrow=arrow(length=unit(0.3,'lines')), color='blue') +
  geom_text(data=spider.env.scores,
    aes(y=PC2*1.1, x=PC1*1.1, label=label, alpa = Flag),
    color='blue', show.legend = FALSE)
g

## ----end

## ---- lm
pc1 <- spider.rda.scores |> filter(score=='sites') |> pull(PC1)
pc2 <- spider.rda.scores |> filter(score=='sites') |> pull(PC2)

lm(1:nrow(spider.env) ~ soil.dry + bare.sand + fallen.leaves +
     moss + herb.layer + reflection, data =  spider.env) |>
  vif()
lm(1:nrow(spider.env) ~ herb.layer + fallen.leaves + bare.sand + moss, data=spider.env) |>
  vif()
lm(pc1 ~ herb.layer + fallen.leaves + bare.sand + moss, data=spider.env) |>
  summary()
lm(pc2 ~ herb.layer + fallen.leaves + bare.sand + moss, data=spider.env) |>
  summary()

## ----end

## ---- RDA spiders

## ---- RDA
spider.rda <- rda(
  spider.std ~
    scale(herb.layer) +
    scale(fallen.leaves) +
    scale(bare.sand) +
    scale(moss),
  data = spider.env,
  scale = FALSE
)
vif.cca(spider.rda)
summary(spider.rda, display=NULL)
## ----end

## ---- goodness of fit
goodness(spider.rda)
goodness(spider.rda, display = "sites")
inertcomp(spider.rda)
inertcomp(spider.rda, proportional = TRUE)
## ----end

## ---- Anova
anova(spider.rda)
anova(spider.rda, by='axis')
anova(spider.rda, by='margin')
## ----end

## ---- other parameters
coef(spider.rda)
RsquareAdj(spider.rda)
## ----end

## ---- ordination plot
screeplot(spider.rda)
autoplot(spider.rda, geom='text')
## ----end

## ---- CA
spider.std <- spider.abund |>
  mutate(across(everything(), ~.x^0.25)) |>
  wisconsin()

spider.std <- (spider.abund^0.25) |>
  wisconsin()
spider.std
spider.ca <- cca(spider.std, scale=FALSE)
## ----end

## ---- CA summary
summary(spider.ca, display=NULL)
## ----end

## ---- CA ordination plot
screeplot(spider.ca)
sum(eigenvals(spider.ca))/length(eigenvals(spider.ca))
eigenvals(spider.ca)/sum(eigenvals(spider.ca))
plot(spider.ca, scaling='species')

autoplot(spider.ca)
autoplot(spider.ca) + theme_bw()
autoplot(spider.ca, geom='text') + theme_bw()
## ----end

## ---- CA ordination plot pretty
spider.ca.scores <- spider.ca |>
  fortify()
spider.ca.scores |> head()

g <-
  ggplot(data = NULL, aes(y=CA2, x=CA1)) +
  geom_hline(yintercept=0, linetype='dotted') +
  geom_vline(xintercept=0, linetype='dotted') +
  geom_point(data=spider.ca.scores %>% filter(score=='sites')) +
  geom_text(data=spider.ca.scores %>% filter(score=='sites'),
    aes(label=label), hjust=-0.2) +
  geom_segment(data=spider.ca.scores %>% filter(score=='species'),
    aes(y=0, x=0, yend=CA2, xend=CA1),
    arrow=arrow(length=unit(0.3,'lines')), color='red') +
  ## geom_text(data=spider.rda.scores %>% filter(score=='species'),
  ##           aes(y=PC2*1.1, x=PC1*1.1, label=label), color='red') +
  geom_text_repel(data=spider.ca.scores %>% filter(score=='species'),
    aes(y=CA2*1.1, x=CA1*1.1, label=label), color='red') +
  theme_bw()
g

## ----end

## ---- CA envfit
spider.envfit <- envfit(spider.ca, env=spider.env)
spider.envfit
autoplot(spider.envfit)

spider.env.scores <- spider.envfit |> fortify()
g <- g +
  geom_segment(data=spider.env.scores,
    aes(y=0, x=0, yend=CA2, xend=CA1),
    arrow=arrow(length=unit(0.3,'lines')), color='blue') +
  geom_text(data=spider.env.scores,
    aes(y=CA2*1.1, x=CA1*1.1, label=label), color='blue')
g

## ---- PCoA
## principal coordinates analysis
spider.dist <- vegdist(spider.std, method='bray')
spider.capscale <- capscale(spider.dist~1, data=spider.env)
summary(spider.capscale, display=NULL)
plot(spider.capscale)
autoplot(spider.capscale, geom='text')

# Distance based redundancy analysis
spider.capscale <- capscale(spider.dist ~
    scale(herb.layer) +
    scale(fallen.leaves) +
    scale(bare.sand) +
    scale(moss),
  data = spider.env)
summary(spider.capscale, display=NULL)
plot(spider.capscale)

summary(spider.capscale, display=NULL)
anova(spider.capscale)

anova(spider.capscale, by='margin')
screeplot(spider.capscale)
sum(eigenvals(spider.capscale))/length(eigenvals(spider.capscale))
eigenvals(spider.capscale)/sum(eigenvals(spider.capscale))


## ---- MDS macnally
macnally <- read.csv('../public/data/macnally_full.csv',strip.white=TRUE)
head(macnally)
macnally <- macnally |>
  mutate(HABITAT = factor(HABITAT, levels = c(
    "Mixed", "Gipps.Manna",
    "Montane Forest", "Foothills Woodland", "Box-Ironbark", "River Red Gum"
  )))

macnally.mds <- metaMDS(macnally[,-1], k=2,  plot=TRUE)
macnally.mds

macnally.mds$stress
stressplot(macnally.mds)

plot(macnally.mds)

macnally.mds.scores <- macnally.mds |>
  fortify() |>
  full_join(macnally |>
             rownames_to_column(var='label'),
    by =  'label')

g <-
    ggplot(data = NULL, aes(y=NMDS2, x=NMDS1)) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_vline(xintercept=0, linetype='dotted') +
    geom_point(data=macnally.mds.scores %>% filter(score=='sites'),
               aes(color=HABITAT)) +
    geom_text(data=macnally.mds.scores %>% filter(score=='sites'),
              aes(label=label, color=HABITAT), hjust=-0.2, show.legend = FALSE) +
    geom_segment(data=macnally.mds.scores %>% filter(score=='species'),
                 aes(y=0, x=0, yend=NMDS2, xend=NMDS1),
                 arrow=arrow(length=unit(0.3,'lines')), color='red',
      alpha =  0.2) +
    geom_text(data=macnally.mds.scores %>% filter(score=='species'),
      aes(y=NMDS2*1.1, x=NMDS1*1.1, label=label), color='red',
      alpha =  0.2)
g


g1 <-
    ggplot(data = NULL, aes(y=NMDS2, x=NMDS1)) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_vline(xintercept=0, linetype='dotted') +
    geom_point(data=macnally.mds.scores %>% filter(score=='sites'),
               aes(color=HABITAT))
g1

g2 <- g +
  stat_density_2d(
    data = macnally.mds.scores %>% filter(score == "sites"),
    geom = "polygon",
    aes(y = NMDS2, x = NMDS1, fill = HABITAT),
    contour_var = "ndensity",
    breaks = c(0.05, 0.1),
    alpha = 0.3,
    position = "identity",
    show.legend = FALSE
  )
g2

centroids <- macnally.mds.scores |>
  filter(score == "sites") |>
  group_by(HABITAT) |>
  summarise(across(c(NMDS1, NMDS2), list(c = mean)))

macnally.mds.scores <- macnally.mds.scores |>
  full_join(centroids)

macnally.mds.scores.centroids <- macnally.mds.scores |>
  filter(score == "sites") |>
  group_by(HABITAT) |>
  summarise(across(c(NMDS1, NMDS2), list(c = mean)))
macnally.mds.scores <- macnally.mds.scores |>
  full_join(macnally.mds.scores.centroids)

g1 <- g1 +
  stat_density_2d(
    data = macnally.mds.scores %>% filter(score == "sites"),
    geom = "polygon",
    aes(y = NMDS2, x = NMDS1, fill = HABITAT),
    contour_var = "ndensity",
    breaks = c(0.05, 0.1),
    alpha = 0.1,
    position = "identity",
    show.legend = FALSE) +
  geom_segment(data = macnally.mds.scores,
  aes(x = NMDS1_c, xend = NMDS1, y = NMDS2_c, yend = NMDS2, colour = HABITAT)) +
  theme_classic()
g1

Xmat <- model.matrix(~-1+HABITAT, data = macnally)
colnames(Xmat) <-gsub("HABITAT","",colnames(Xmat))
envfit <- envfit(macnally.mds, env=Xmat)
envfit


macnally.env.scores <- envfit |> fortify()
g3 <- g1 +
    geom_segment(data=macnally.env.scores,
                 aes(y=0, x=0, yend=NMDS2, xend=NMDS1),
                 arrow=arrow(length=unit(0.3,'lines')), color='blue') +
    geom_text(data=macnally.env.scores,
              aes(y=NMDS2*1.1, x=NMDS1*1.1, label=label), color='blue')
g3

macnally.dist <- vegdist(macnally[,-1], 'bray')

adonis2(macnally.dist ~ HABITAT, data=macnally)

mm <-  model.matrix(~-1 + HABITAT, data=macnally)
head(mm)
colnames(mm) <-gsub("HABITAT","",colnames(mm))
mm <- data.frame(mm)

macnally.adonis<-adonis2(macnally.dist ~
                           Mixed + Box.Ironbark + Foothills.Woodland + Gipps.Manna +
                           Montane.Forest + River.Red.Gum,
  data=mm,
  by = "terms",
  perm=9999)
print(macnally.adonis)

macnally.adonis<-adonis2(macnally.dist ~ Box.Ironbark + Foothills.Woodland + Gipps.Manna +
                           Montane.Forest + River.Red.Gum,
  data=mm,
  by = "margin",
  perm=9999)
print(macnally.adonis)

macnally.disp <- betadisper(macnally.dist, macnally$HABITAT)
boxplot(macnally.disp)
plot(macnally.disp)
anova(macnally.disp)
permutest(macnally.disp, pairwise = TRUE)
TukeyHSD(macnally.disp)

macnally.std <- wisconsin(macnally[,c(-1)]^0.25)
simper(macnally.std, macnally$HABITAT) |> summary()

## ----end
