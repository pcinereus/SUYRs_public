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
library(mvabund)
library(scales)
library(ggvegan)
library(ggrepel)
library(glmmTMB)
library(gllvm)
library(EcolUtils)
library(rstan)
library(ggforce)
library(concaveman)
## ----end


## PCA ------------------------------------------------------

## ---- PCA - spiders

## ---- readSpider
spider.abund <- read_csv(file = "../public/data/spider.abund.csv", trim_ws = TRUE)
spider.env <- read_csv(file = "../public/data/spider.env.csv", trim_ws = TRUE)
glimpse(spider.abund)
glimpse(spider.env)
## ----end

glimpse(spider.abund)

## ---- EDA spider
spider.abund %>%
  cor %>%
  corrplot(type = 'upper',
    diag = FALSE)
## And now with axes arrange according to first princomp
spider.abund %>%
  cor %>%
  corrplot(type = 'upper',
    order = 'FPC',
    diag = FALSE)
## ----end
## ---- EDA 2
spider.abund %>%
  ggpairs(lower = list(continuous = "smooth"),
    diag = list(continuous = "density"),
    axisLabels = "show")
## - clearly not normal

spider.abund^0.25 %>%
  ggpairs(lower = list(continuous = "smooth"),
    diag = list(continuous = "density"),
    axisLabels = "show")
## - still not normal
## - but at least trends sort of linear
## ----end

## Conclusions:

## ---- PCA prep
spider.std <- spider.abund |> 
  mutate(across(everything(), ~.x^0.25)) |> 
  wisconsin()
spider.std
## ----end

## Run PCA - unconstrained axis rotations
## normally we would not scale, but for illustrative purposes...
## ---- PCA
## run one of the following
spider.rda <- rda(spider.abund, scale=TRUE)
spider.rda <- rda(spider.std, scale=TRUE)
spider.rda <- rda(spider.std, scale=FALSE)
spider.std
## ----end

## ---- PCA summary
summary(spider.rda, display=NULL)

screeplot(spider.rda)
abline(a=1,b=0)
## ----end


## ---- PCA ordinations part 1
scores(spider.rda, choices=1:3, display='sites')
scores(spider.rda, choices=1:3, display='species')

## Quick and nasty ordination plots
biplot(spider.rda, scaling='species')
biplot(spider.rda, scaling='sites')
## ----end

## ---- PCA ordinations part 2
## Quick and nasty ordination plots
pl<-vegan::ordiplot(spider.rda)
points(pl, "sites", pch=21, col="red", bg="yellow")
text(pl, "sites", col="red", cex=0.9)
text(pl, "species", col="blue", cex=0.9)
## line(pl, "sites")
## ----end


## ---- PCA ordinations part 3
## ggvegan provides ways to use ggplot
## library(devtools)
## devtools::install_github("gavinsimpson/ggvegan")
## library(ggvegan)
autoplot(spider.rda)
autoplot(spider.rda) + theme_bw()
autoplot(spider.rda,geom='text') + theme_bw()
## ----end


## ---- PCA ordinations part 4
spider.rda.scores <- spider.rda |> 
  fortify()
spider.rda.scores

ggplot(data = NULL, aes(y=PC2, x=PC1)) +
  geom_hline(yintercept=0, linetype='dotted') +
  geom_vline(xintercept=0, linetype='dotted') +
  geom_point(data=spider.rda.scores %>% filter(Score=='sites')) +
  geom_text(data=spider.rda.scores %>% filter(Score=='sites'),
    aes(label=Label), hjust=-0.2) +
  geom_segment(data=spider.rda.scores %>% filter(Score=='species'),
    aes(y=0, x=0, yend=PC2, xend=PC1),
    arrow=arrow(length=unit(0.3,'lines')), color='red') +
  geom_text(data=spider.rda.scores %>% filter(Score=='species'),
    aes(y=PC2*1.1, x=PC1*1.1, label=Label), color='red')
## ----end

## lets start to build this up
## ---- PCA ordinations part 5
g <-
  ggplot(data = NULL, aes(y=PC2, x=PC1)) +
  geom_hline(yintercept=0, linetype='dotted') +
  geom_vline(xintercept=0, linetype='dotted') +
  geom_point(data=spider.rda.scores %>% filter(Score=='sites')) +
  geom_text(data=spider.rda.scores %>% filter(Score=='sites'),
    aes(label=Label), hjust=-0.2) +
  geom_segment(data=spider.rda.scores %>% filter(Score=='species'),
    aes(y=0, x=0, yend=PC2, xend=PC1),
    arrow=arrow(length=unit(0.3,'lines')), color='red') +
  geom_text_repel(data=spider.rda.scores %>% filter(Score=='species'),
    aes(y=PC2*1.1, x=PC1*1.1, label=Label), color='red') +
  theme_bw()
g
## ----end

## Nice axes titles
## ---- PCA ordinations part 6
eig <- eigenvals(spider.rda)

g <- g +
  scale_y_continuous(paste(names(eig[2]), sprintf('(%0.1f%% explained var.)',
    100 * eig[2]/sum(eig))))+
  scale_x_continuous(paste(names(eig[1]), sprintf('(%0.1f%% explained var.)',
    100 * eig[1]/sum(eig))))

## ----end

#put a circle
## ---- PCA ordinations part 6
circle.prob <- 0.68
r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(spider.rda$CA$u[,1:2]^2))^(1/4)
theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
circle <- data.frame(PC1 = r * cos(theta), PC2 = r * sin(theta))
g <- g + geom_path(data = circle, aes(y=PC2,x=PC1), color = muted('white'), size = 1/2, alpha = 1/3)
g

## ----end

## de-emphasize those that are inside the circle (to prevent cluttering)
## ---- PCA ordinations part 7
spider.rda.scores.sites <- spider.rda.scores |>
  filter(Score == "sites") |>
  mutate(Flag = factor(ifelse(sqrt(PC1^2 + PC2^2) > r, 1, 0)))
spider.rda.scores.species <- spider.rda.scores |>
  filter(Score == "species") |>
  mutate(Flag = factor(ifelse(sqrt(PC1^2 + PC2^2) > r, 1, 0)))

g <-
  ggplot(data = NULL, aes(y=PC2, x=PC1)) +
  geom_hline(yintercept=0, linetype='dotted') +
  geom_vline(xintercept=0, linetype='dotted') +
  geom_point(data=spider.rda.scores.sites, aes(alpha = Flag), show.legend = FALSE) +
  geom_text(data=spider.rda.scores.sites,
    aes(label=Label, alpha = Flag), hjust=-0.2, show.legend = FALSE) +
  geom_segment(data=spider.rda.scores.species,
    aes(y=0, x=0, yend=PC2, xend=PC1, alpha = Flag),
    arrow=arrow(length=unit(0.3,'lines')), color='red', show.legend = FALSE) +
  geom_text_repel(data=spider.rda.scores.species,
    aes(y=PC2*1.1, x=PC1*1.1, label=Label, alpha = Flag),
    color='red', show.legend = FALSE) +
  geom_path(data = circle, aes(y=PC2,x=PC1),
    color = muted('white'), size = 1/2, alpha = 1/3) +
  scale_y_continuous(paste(names(eig[2]), sprintf('(%0.1f%% explained var.)',
    100 * eig[2]/sum(eig))))+
  scale_x_continuous(paste(names(eig[1]), sprintf('(%0.1f%% explained var.)',
    100 * eig[1]/sum(eig)))) +
  theme_bw()
g
## ----end

## ---- Envfit prep
spider.env |>
  cor() |> 
  corrplot(type = 'upper',
    order = 'FPC',
    diag = FALSE)

spider.env |> 
  ggpairs(lower = list(continuous = "smooth"),
    diag = list(continuous = "density"),
    axisLabels = "show")
## ----end

## ---- Envfit prep
spider.envfit <- envfit(spider.rda, env = spider.env)
spider.envfit

## ----end

## ---- Envfit plot
spider.env.scores <- spider.envfit |>
  fortify() |> 
  mutate(Flag = factor(ifelse(sqrt(PC1^2 + PC2^2) > r, 1, 0)))
g <- g + 
  geom_segment(data=spider.env.scores,
    aes(y=0, x=0, yend=PC2, xend=PC1, alpha = Flag, show.legend = FALSE),
    arrow=arrow(length=unit(0.3,'lines')), color='blue') +
  geom_text(data=spider.env.scores,
    aes(y=PC2*1.1, x=PC1*1.1, label=Label, alpa = Flag),
    color='blue', show.legend = FALSE)
g

## ----end

## ---- lm
pc1 <- spider.rda.scores |> filter(Score=='sites') |> pull(PC1)
pc2 <- spider.rda.scores |> filter(Score=='sites') |> pull(PC2)

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

## ----end
## ---- RDA spiders
## ---- RDA
spider.rda <- rda(spider.std ~ 
                    scale(herb.layer)+
                    scale(fallen.leaves) +
                    scale(bare.sand) +
                    scale(moss),
  data=spider.env, scale=FALSE)
## ----end

## ---- RDA summary
summary(spider.rda, display=NULL)

vif.cca(spider.rda)
## ----end

## ---- goodness of fit
goodness(spider.rda)
goodness(spider.rda, display = "sites")
inertcomp(spider.rda)
inertcomp(spider.rda, proportional = TRUE)
## ----end

## ---- Anova
##overall test
anova(spider.rda)
anova(spider.rda, by='axis')
anova(spider.rda, by='margin')

RsquareAdj(spider.rda)

screeplot(spider.rda)
## ----end

## ---- RDA plot
autoplot(spider.rda, geom='text')
## ----end

## ----END

## ---- CA vegetation
## ---- CA get data
data <- read_csv('../public/data/data.csv', trim_ws=TRUE)
head(data) 
enviro <- read_csv('../public/data/enviro.csv', trim_ws=TRUE)
head(enviro)
enviro <- enviro %>% mutate(Substrate=factor(Substrate))
## ----end

## ---- CA std data
data.std <- data %>%
  dplyr::select(-Sites) %>% 
  decostand(method="total",MARGIN=2)
data.std

data.std %>% cor %>%
  corrplot(diag=FALSE)

data.std %>% cor %>%
  corrplot(diag=FALSE, order='FPC')
## ----end

## ---- CA fit
data.ca <- cca(data.std, scale=FALSE)
summary(data.ca, display=NULL)
#anova(data.ca)

screeplot(data.ca)
sum(eigenvals(data.ca))/length(eigenvals(data.ca))
eigenvals(data.ca)/sum(eigenvals(data.ca))
## ----end

## ---- CA plot part 1
autoplot(data.ca)
autoplot(data.ca) + theme_bw()
autoplot(data.ca, geom='text') + theme_bw()

data.ca.scores <- data.ca |> 
  fortify()
data.ca.scores |> head()

## ----end

## ---- CA plot part 2
g <-
  ggplot(data = NULL, aes(y=CA2, x=CA1)) +
  geom_hline(yintercept=0, linetype='dotted') +
  geom_vline(xintercept=0, linetype='dotted') +
  geom_point(data=data.ca.scores %>% filter(Score=='sites')) +
  geom_text(data=data.ca.scores %>% filter(Score=='sites'),
    aes(label=Label), hjust=-0.2) +
  geom_segment(data=data.ca.scores %>% filter(Score=='species'),
    aes(y=0, x=0, yend=CA2, xend=CA1),
    arrow=arrow(length=unit(0.3,'lines')), color='red') +
  ## geom_text(data=data.rda.scores %>% filter(Score=='species'),
  ##           aes(y=PC2*1.1, x=PC1*1.1, label=Label), color='red') +
  geom_text_repel(data=data.ca.scores %>% filter(Score=='species'),
    aes(y=CA2*1.1, x=CA1*1.1, label=Label), color='red') +
  theme_bw()
g
## ----end

## ---- CA envfit
Xmat <- model.matrix(~-1+pH+Slope+Altitude+Substrate, enviro)
data.envfit <- envfit(data.ca, env=Xmat)
data.envfit
## ----end

## ---- CA envfit plot
autoplot(data.envfit)

data.env.scores <- data.envfit |> fortify()
g <- g + 
  geom_segment(data=data.env.scores,
    aes(y=0, x=0, yend=CA2, xend=CA1),
    arrow=arrow(length=unit(0.3,'lines')), color='blue') +
  geom_text(data=data.env.scores,
    aes(y=CA2*1.1, x=CA1*1.1, label=Label), color='blue')
g
## ----end

## ---- CA envfit lm
data.ca.scores <- data.ca %>% fortify()
CA1 <- data.ca.scores %>% filter(Score =='sites') %>% pull(CA1)
CA2 <- data.ca.scores %>% filter(Score =='sites') %>% pull(CA2)
summary(lm(CA1 ~ pH+Slope+Altitude+Substrate, data=enviro))
summary(lm(CA2 ~ pH+Slope+Altitude+Substrate, data=enviro))
## ----end

## ---- CCA

## ---- CCA fit
data.cca <- cca(data.std~pH + Altitude + Substrate + Slope, data=enviro, scale=FALSE)

## ----end

## ---- CCA summary
summary(data.cca, display=NULL)

## ----end

## ---- CCA plot
autoplot(data.cca)
## ----end

## ---- CCA anova
vif.cca(data.cca)
#overall test
anova(data.cca)
anova(data.cca, by='axis')
anova(data.cca, by='margin')
#anova(data.cca, by='margin', scope="pH")

coef(data.cca)

RsquareAdj(data.cca)

screeplot(data.cca)
int <- data.cca$tot.chi/length(data.cca$CA$eig)
abline(h=int)

## ----end
## ----end


## ---- PCoA

## ---- PCoA prep
## principal coordinates analysis
data.std <-
  data %>% dplyr::select(-Sites) %>%
  decostand(method="total",MARGIN=2)
data.std
data.dist = vegdist(data.std, method='bray')
## ----end

## ---- PCoA fit
data.capscale = capscale(data.dist~1, data=enviro)
## ----end

## ---- PCoA summary
summary(data.capscale, display=NULL)
## ----end

## ---- PCoA plot
autoplot(data.capscale, geom='text')
## ----end

## ---- PCoA 2 fit
data.capscale = capscale(data.dist~scale(pH) + scale(Altitude) + Substrate + scale(Slope), data=enviro)
## ----end

## ---- PCoA 2 summary
summary(data.capscale, display=NULL)
anova(data.capscale)

anova(data.capscale, by='margin')
screeplot(data.capscale)
sum(eigenvals(data.capscale))/length(eigenvals(data.capscale))
eigenvals(data.capscale)/sum(eigenvals(data.capscale))
## ----end

## ---- PCoA 2 plot
autoplot(data.capscale)

## ----end

## ---- PCoA 3
data.capscale = capscale(data.dist~pH + Condition(Altitude) + Substrate + Slope, data=enviro)
summary(data.capscale, display=NULL)
plot(data.capscale)
## ----end
## ----end

## ---- MDS macnally
## ---- MDS read data
## MUST READ IN THIS WAY..
macnally <- read.csv('../public/data/macnally_full.csv',strip.white=TRUE)
head(macnally)
macnally[1:5,1:5]

## ----end

## ---- MDS fit
macnally.mds <- metaMDS(macnally[,-1], k=2,  plot=TRUE)
macnally.mds
stressplot(macnally.mds)
## ----end

## ---- MDS plot part 1
macnally.mds.scores <- macnally.mds |> 
  fortify() |> 
  ## mutate(across(c(NMDS1, NMDS2), as.numeric)) |> 
  full_join(macnally |> add_rownames(var='Label'))

g <-
  ggplot(data = NULL, aes(y=NMDS2, x=NMDS1)) +
  geom_hline(yintercept=0, linetype='dotted') +
  geom_vline(xintercept=0, linetype='dotted') +
  geom_point(data=macnally.mds.scores %>% filter(Score=='sites'),
    aes(color=HABITAT)) +
  geom_text(data=macnally.mds.scores %>% filter(Score=='sites'),
    aes(label=Label, color=HABITAT), hjust=-0.2) +
  geom_segment(data=macnally.mds.scores %>% filter(Score=='species'),
    aes(y=0, x=0, yend=NMDS2, xend=NMDS1),
    arrow=arrow(length=unit(0.3,'lines')), color='red',
    alpha =  0.2) +
  geom_text(data=macnally.mds.scores %>% filter(Score=='species'),
    aes(y=NMDS2*1.1, x=NMDS1*1.1, label=Label), color='red',
    alpha =  0.2) 
g
## ----end

# Nice axes titles
## ---- MDS plot part 2
## these are alternatives
g + ggforce::geom_mark_ellipse(data=macnally.mds.scores %>% filter(Score=='sites'),
  aes(y=NMDS2, x=NMDS1, fill=HABITAT), expand=0) 
g + ggforce::geom_mark_hull(data=macnally.mds.scores %>% filter(Score=='sites'),
  aes(y=NMDS2, x=NMDS1, fill=HABITAT), expand=0.01) 
g + ggforce::geom_mark_hull(data=macnally.mds.scores %>% filter(Score=='sites'),
  aes(y=NMDS2, x=NMDS1, fill=HABITAT), expand=0.03, concavity = 10) 

g <- g + ggforce::geom_mark_hull(data=macnally.mds.scores %>% filter(Score=='sites'),
  aes(y=NMDS2, x=NMDS1, fill=HABITAT), expand=0.03, concavity = 10) 
macnally.mds.scores.centroids <- macnally.mds.scores |>
  filter(Score == "sites") |>
  group_by(HABITAT) |>
  summarise(across(c(NMDS1, NMDS2), list(c = mean)))
macnally.mds.scores <- macnally.mds.scores |>
  full_join(macnally.mds.scores.centroids)
g + geom_segment(data = macnally.mds.scores,
  aes(x = NMDS1_c, xend = NMDS1, y = NMDS2_c, yend = NMDS2, colour = HABITAT))

## ----end

# Nice axes titles
## ---- MDS envfit 
Xmat <- model.matrix(~-1+HABITAT, data=macnally)
colnames(Xmat) <-gsub("HABITAT","",colnames(Xmat))
envfit <- envfit(macnally.mds, env=Xmat)
envfit


macnally.env.scores <- envfit %>% fortify()
g <- g + 
  geom_segment(data=macnally.env.scores,
    aes(y=0, x=0, yend=NMDS2, xend=NMDS1),
    arrow=arrow(length=unit(0.3,'lines')), color='blue') +
  geom_text(data=macnally.env.scores,
    aes(y=NMDS2*1.1, x=NMDS1*1.1, label=Label), color='blue')
g

## ----end

## ---- MDS simper 
macnally.std <- wisconsin(macnally[,c(-1)]^0.25)
simper(macnally.std, macnally$HABITAT) |> summary()
## ----end

## ---- MDS betadisp part 1
macnally.disp <- betadisper(macnally.dist, macnally$HABITAT)
boxplot(macnally.disp)
plot(macnally.disp)
anova(macnally.disp)
permutest(macnally.disp, pairwise = TRUE)
TukeyHSD(macnally.disp)
## ----end
## ---- MDS betadisp part 1
macnally.disp <- betadisper(macnally.dist, macnally$HABITAT, type="median",bias.adjust = TRUE)
boxplot(macnally.disp)
plot(macnally.disp)
anova(macnally.disp)
permutest(macnally.disp, pairwise = TRUE)
TukeyHSD(macnally.disp)
## ----end
## ----end


## ---- dune
## ---- dune read data
dune <- read_csv('../public/data/dune.csv', trim_ws=TRUE)
dune <- dune %>% mutate(MANAGEMENT=factor(MANAGEMENT,  levels=c("NM","BF","HF","SF"))) %>%
  as.data.frame()
#dune <- read.csv('../downloads/data/dune.csv')
dune
## ----end

## ---- dune fit
dune.mds <- metaMDS(dune[,-1], k=2)

## ----end

## ---- dune plot
autoplot(dune.mds, geom=c('text'))
## ----end

## ---- dune adonis2 part 1
dune.adonis<-adonis2(dune.dist~MANAGEMENT,  data=dune)
dune.adonis

## ----end

## ---- dune adonis2 part 2
mm <- model.matrix(~MANAGEMENT, data=dune)
head(mm)
colnames(mm) <-gsub("MANAGEMENT","",colnames(mm))
mm <- data.frame(mm)
dune.adonis<-adonis2(dune.dist~BF+HF+SF, data=mm,
  perm=9999)
dune.adonis

## ----end

## ---- dune adonis2 part 3
library(pairwiseAdonis)
pairwise.adonis(dune.dist, dune$MANAGEMENT)

library(EcolUtils)
adonis.pair(dune.dist, dune$MANAGEMENT, nper = 10000)

## ----end

## ---- dune simper
dune.simper=simper(dune[,-1], dune[,1], permutations = 999)
summary(dune.simper)
## ----end

## ---- dune mrpp
dune.mrpp = mrpp(dune.dist, dune[,1], permutations=999)
dune.mrpp
hist(dune.mrpp$boot.deltas)
# Chance corrected within-group agreement = 1-Obs delta / exp delta
dune.meandist = meandist(dune.dist, dune[,1], permutations=999)
dune.meandist
summary(dune.meandist)
plot(dune.meandist)
## ----end

## ---- dune betadisp
dune.disp <- betadisper(dune.dist,  group=dune$MANAGEMENT)
permutest(dune.disp)
permutest(dune.disp, pairwise = TRUE)
boxplot(dune.disp)
plot(dune.disp)
anova(dune.disp)
TukeyHSD(dune.disp)
## ----end

## ----end

## ---- varveg
## ---- varveg read data
vareveg <- read.csv('../public/data/vareveg.csv')
head(vareveg)
vareenv <- read.csv('../public/data/vareenv.csv')
head(vareenv)
## ----end

## ---- varveg bioenv
vareveg.dist <- vegdist(wisconsin(vareveg[,-1]^0.25),'bray')
vareenv.std <- decostand(vareenv[,-1], "standardize")
vareenv.dist <- vegdist(vareenv.std, "euc")
bioenv(vareveg.dist, vareenv.std)
## ----end

## ---- varveg adonis
adonis2(vareveg.dist~Ca+Fe+Mn+Baresoil, data=vareenv.std)
## ----end
## ----end

## what about hierarchical designs?
## ---- hierachical
## ---- hierachical read data
brink <- read_csv(file='../public/data/brink.csv', trim_ws=TRUE)
brink <- brink %>% mutate(WEEK = factor(WEEK),
  TREATMENT = factor(TREATMENT),
  DITCH = factor(DITCH))
## Isolate just the invertegrate data
inverts <- brink %>% dplyr::select(-WEEK, -TREATMENT, -DITCH)
## ----end

## ---- hierachical RDA
inverts.rda <-  rda(wisconsin(inverts^0.25) ~ TREATMENT*WEEK + Condition(DITCH), data=brink)
summary(inverts.rda, display=NULL)
inverts.rda %>% autoplot(geom = 'text')
inverts.rda %>% autoplot()
anova(inverts.rda)
anova(inverts.rda, by='terms')
anova(inverts.rda, by='margin')
## ----end

## ---- hierachical adonis
adonis2(inverts~ DITCH+TREATMENT*WEEK, data=brink, strata = brink$DITCH)
## ----end
## ----end


## ---- MVABUND spiders
## ---- MVABUND fit part 1
mva <- mvabund(spider.abund)
spider.mod <- manyglm(mva~
                        scale(soil.dry)+
                        scale(moss)+
                        scale(herb.layer)+
                        scale(bare.sand),
  family=poisson(link='log'),
  data=spider.env)
## ----end
## ---- MVABUND fit part 2
plot(spider.mod)
## ----end
## ---- MVABUND fit part 3
spider.mod1 <- manyglm(mva~
                         scale(soil.dry)+
                         scale(moss)+
                         scale(herb.layer)+
                         scale(bare.sand), 
  family="negative.binomial",
  data=spider.env)
## ----end
## ---- MVABUND fit part 4
plot(spider.mod1)
spider.mod1
spider.mod1 |> summary()
## ----end
## ---- MVABUND fit part 5
anova(spider.mod)
anova(spider.mod, test='LR')
anova(spider.mod, cor.type = 'R')
anova(spider.mod, cor.type = 'shrink')
anova(spider.mod, p.uni='adjusted')
summary(spider.mod, test="LR")
## ----end
## ----end

## ---- gllvm spiders
fitx <- gllvm(y = spider.abund, X=spider.env,
  formula =  ~ soil.dry + fallen.leaves + moss + herb.layer + bare.sand,
  family = "negative.binomial")
fitx
par(mfrow = c(1,2))
plot(fitx, which = 1:2)
summary(fitx)
coefplot(fitx, mfrow = c(3,2), cex.ylab = 0.8)
crx <- getResidualCor(fitx)
corrplot(crx, diag = FALSE, type = "lower", method = "square", tl.srt = 25)

ordiplot(fitx, biplot = TRUE)
ordiplot(fitx, biplot = TRUE)
abline(h = 0, v = 0, lty=2)

## ----end

## ---- gllvm microbial
data(microbialdata)
X <- microbialdata$Xenv
y <- microbialdata$Y[, order(colMeans(microbialdata$Y > 0), 
  decreasing = TRUE)[21:40]]
fit <- gllvm(y, X, formula = ~ pH + Phosp, family = poisson())
fit$logL
ordiplot(fit)
coefplot(fit)
Site<-data.frame(Site=X$Site)
Xsoils <- cbind(scale(X[, 1:3]),Site)
ftXph <- gllvm(y, Xsoils, formula = ~pH, family = "negative.binomial", 
  row.eff = ~(1|Site), num.lv = 2)
Xenv <- data.frame(X, Region = factor(X$Region),
  Soiltype = factor(X$Soiltype))
ftXi <- gllvm(y, Xenv, formula = ~ SOM + pH + Phosp + Region, 
  family = "negative.binomial", row.eff = ~(1|Site), num.lv = 2,
  sd.errors = FALSE)

ph <- Xenv$pH
rbPal <- colorRampPalette(c('mediumspringgreen', 'blue'))
Colorsph <- rbPal(20)[as.numeric(cut(ph, breaks = 20))]
pchr = NULL
pchr[Xenv$Region == "Kil"] = 1
pchr[Xenv$Region == "NyA"] = 2
pchr[Xenv$Region == "Aus"] = 3
ordiplot(ftXi, main = "Ordination of sites",  
  symbols = TRUE, pch = pchr, s.colors = Colorsph)
legend("topleft", legend = c("Kil", "NyA", "Mayr"), pch = c(1, 2, 3), bty = "n")

ftNULL <- gllvm(y, X = data.frame(Site = X[,5]), 
  family = "negative.binomial", row.eff = ~(1|Site), num.lv = 2,
  sd.errors = FALSE)
1 - getResidualCov(ftXi)$trace/getResidualCov(ftNULL)$trace
## ----end

## ---- boral spiders
#####################
library(boral)

data(spider)
y <- spider$abun
X <- scale(spider$x)
n <- nrow(y)
p <- ncol(y)
example_mcmc_control <- list(n.burnin = 10, n.iteration = 100,
  n.thin = 1)
spiderfit_nb <- boral(y, family = "negative.binomial", 
  lv.control = list(num.lv = 2), row.eff = "fixed", 
  mcmc.control = example_mcmc_control, model.name = "spider.jags"
)

lvsplot(spiderfit_nb, biplot = TRUE)

spiderfit_nb <- boral(y, x = X, family = "negative.binomial", 
  lv.control = list(num.lv = 2), row.eff = "fixed", 
  mcmc.control = example_mcmc_control, model.name = "spider.jags"
)
lvsplot(spiderfit_nb, biplot = TRUE)


## ----end

## ---- stan3
## Prepare data for stan model
Y <- as.matrix(spider.abund)
N <- nrow(Y)
P <- ncol(Y)
D <- 2


## Fit the model with stan
StanDat <- list(Y = Y, N = N, P = P, D = D)
lvm_stan <- rstan::stan("model2.stan", data = StanDat,
  iter = 2000, chain = 3, cores = 3)## ---- stan
lvm_stan

StanDat <- list(N = N, S = P, D = D, Y = as.matrix(Y))
lvm_cmdstan <- cmdstanr::cmdstan_model("model3.stan")
##lvm_cmdstan$print()
lvm_fit <- lvm_cmdstan$sample(data = StanDat,
  seed = 1,
  chains = 3,
  parallel_chains = 3,
  iter_warmup = 1000,
  iter_sampling = 5000,
  adapt_delta =  0.99,
  max_treedepth = 20
)
## lvm_fit$summary() |> as.data.frame() 
## traceplot(as.mcmc(as.matrix(lvm_fit)), c("alpha"))
coords <- lvm_fit$summary() |> filter(str_detect(variable, "LV\\[")) |> select(variable, median)
## coords
coords <- coords |> mutate(site = rep(1:28, each = 2), x = rep(1:2, 28))
coords <- coords |> select(median, site, x) |> pivot_wider(id_cols = site, names_from = x, values_from = median)
loadings <- lvm_fit$summary() |> filter(str_detect(variable, "Lambda\\[")) |> select(variable, median)
## loadings
loadings <- loadings |> mutate(Species = rep(1:12, 2), X = rep(1:2, each = 12))
loadings <- loadings |> select(median, Species, X) |> pivot_wider(id_cols = Species, names_from = X, values_from = median) |>
  mutate(Species = colnames(spider.abund))
ggplot(coords, aes(y = `1`, x = `2`)) +
  geom_text(aes(label =  site)) +
  geom_text(data =  loadings, aes(label =  Species), color =  "red") +
  geom_segment(data = loadings, aes(y = 0, yend = `1`, x = 0, xend = `2`),
    arrow=arrow(length=unit(0.3,'lines')), colour =  "red") 


## ----end


## ---- PCA grasslands
grasslands.spe <- read.delim ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/grasslands-spe.txt', row.names = 1)
grasslands.env <- read.delim ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/data/grasslands-env.txt')
{
    ## ---- EDA 1
    grasslands.spe %>%
        cor %>%
        corrplot(type = 'upper',
                 diag = FALSE)
    ## And now with axes arrange according to first princomp
    grasslands.spe %>%
        cor %>%
        corrplot(type = 'upper',
                 order = 'FPC',
                 diag = FALSE)
    ## ----end
    ## ---- EDA 2
    scatterplotMatrix(x = sqrt(grasslands.spe[,sample(1:171, 20)]))
    grasslands.spe[,sample(1:171, 20)] %>%
        ggpairs(lower = list(continuous = "smooth"),
                diag = list(continuous = "density"),
                axisLabels = "show")
    ## - clearly not normal

    grasslands.spe[,sample(1:171, 20)]^0.25 %>%
        ggpairs(lower = list(continuous = "smooth"),
                diag = list(continuous = "density"),
                axisLabels = "show")
    grasslands.spe1 <- decostand(grasslands.spe^0.25, "hell")
    grasslands.spe1[,sample(1:171, 20)] %>%
        ggpairs(lower = list(continuous = "smooth"),
                diag = list(continuous = "density"),
                axisLabels = "show")
    ## - still not normal
    ## ----end
    ## ---- DCA (detrended correspondance analysis
    grasslands.dca <- decorana(grasslands.spe)
    grasslands.dca
    
    grasslands.dca <- decorana(grasslands.spe1)
    grasslands.dca
    ## ----end
    ## ---- PCA
    grasslands.std <- grasslands.spe %>%
        mutate(across(everything(), function(x) x^0.25)) %>%
        wisconsin()
    grasslands.std
    grasslands.hel <- grasslands.spe %>%
        decostand(method = "hellinger")

    ## Run PCA - unconstrained axis rotations
    ## normally we would not scale, but for illustrative purposes...
    grasslands.rda <- rda(grasslands.spe, scale=TRUE)
    grasslands.rda <- rda(grasslands.std, scale=TRUE)
    grasslands.rda <- rda(grasslands.std, scale=FALSE)
    grasslands.rda <- rda(grasslands.hel, scale=FALSE)
    ## grasslands.rda <- rda(grasslands.std, scale=FALSE)
    summary(grasslands.rda, display=NULL)
    
    screeplot(grasslands.rda)
    abline(a=1,b=0)
    
    scores(grasslands.rda, choices=1:3, display='sites')
    scores(grasslands.rda, choices=1:3, display='species')
    
    ## Quick and nasty ordination plots
    biplot(grasslands.rda, scaling='species')
    biplot(grasslands.rda, scaling='sites')
    
    ## Quick and nasty ordination plots
    pl<-vegan::ordiplot(grasslands.rda)
    points(pl, "sites", pch=21, col="red", bg="yellow")
    text(pl, "sites", col="red", cex=0.9)
    text(pl, "species", col="blue", cex=0.9)
    ## line(pl, "sites")

    
    ## ggvegan provides ways to use ggplot
    ## library(devtools)
    ## devtools::install_github("gavinsimpson/ggvegan")
    ## library(ggvegan)
    autoplot(grasslands.rda)
    autoplot(grasslands.rda) + theme_bw()
    autoplot(grasslands.rda,geom='text') + theme_bw()

    ## Alternatively, we can extract the scores from the PCA/RDA etc
    ## into a dataframe.
    ## - fortify() is a general method used to extract/convert data into
    ##   a form that is suitable for ggplot.
    ## - the ggvegan package defines some fortify methods specifically
    ##   for analyses created with the vegan package

    grasslands.rda.scores <- grasslands.rda %>%
        fortify()
    grasslands.rda.scores

    g <-
        ggplot(data = NULL, aes(y=PC2, x=PC1)) +
        geom_hline(yintercept=0, linetype='dotted') +
        geom_vline(xintercept=0, linetype='dotted') +
        geom_point(data=grasslands.rda.scores %>% filter(Score=='sites')) +
        geom_text(data=grasslands.rda.scores %>% filter(Score=='sites'),
                  aes(label=Label), hjust=-0.2) +
        geom_segment(data=grasslands.rda.scores %>% filter(Score=='species'),
                     aes(y=0, x=0, yend=PC2, xend=PC1),
                     arrow=arrow(length=unit(0.3,'lines')), color='red') +
        ## geom_text(data=grasslands.rda.scores %>% filter(Score=='species'),
        ##           aes(y=PC2*1.1, x=PC1*1.1, label=Label), color='red') +
        geom_text_repel(data=grasslands.rda.scores %>% filter(Score=='species'),
                        aes(y=PC2*1.1, x=PC1*1.1, label=Label), color='red') +
        theme_bw()
    g

                                        # Nice axes titles
    eig <- eigenvals(grasslands.rda)

    paste(names(eig[2]), sprintf('(%0.1f%% explained var.)', 100 * eig[2]/sum(eig)))
    g <- g +
        scale_y_continuous(paste(names(eig[2]), sprintf('(%0.1f%% explained var.)',
                                                        100 * eig[2]/sum(eig))))+
        scale_x_continuous(paste(names(eig[1]), sprintf('(%0.1f%% explained var.)',
                                                        100 * eig[1]/sum(eig))))

                                        #put a circle
    circle.prob <- 0.68
    ## circle.prob <- 0.95
    ## circle.prob <- 0.95
    r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(spider.rda$CA$u[,1:2]^2))^(1/4)
    theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
    circle <- data.frame(PC1 = r * cos(theta), PC2 = r * sin(theta))
    g <- g + geom_path(data = circle, aes(y=PC2,x=PC1), color = muted('white'), size = 1/2, alpha = 1/3)
    g
    
    ## ----end
    ## ---- Envfit
    X <- grasslands.env %>% as.data.frame()
    X %>%
        cor %>%
        corrplot(type = 'upper',
                 order = 'FPC',
                 diag = FALSE)
    X %>%
        as.data.frame() %>%
        ggpairs(lower = list(continuous = "smooth"),
                diag = list(continuous = "density"),
                axisLabels = "show")
    grasslands.env <- envfit(grasslands.rda, env = X, na.rm = TRUE)
    grasslands.env

    
    grasslands.env.scores <- grasslands.env %>% fortify()
    g <- g + 
        geom_segment(data=grasslands.env.scores,
                     aes(y=0, x=0, yend=PC2, xend=PC1),
                     arrow=arrow(length=unit(0.3,'lines')), color='blue') +
        geom_text(data=grasslands.env.scores,
                  aes(y=PC2*1.1, x=PC1*1.1, label=Label), color='blue')
    g

    ## ----end
}


## ----end


## ---- PCA - Vegetation
data <- read_csv('../public/data/data.csv', trim_ws=TRUE)
head(data) 
enviro <- read_csv('../public/data/enviro.csv', trim_ws=TRUE)
head(enviro)

enviro <- enviro %>% mutate(Substrate=factor(Substrate))
## EDA
data[,-1] %>%
    cor %>%
    corrplot(type='upper', diag=FALSE)
data[,-1] %>%
    cor %>%
    corrplot(type='upper', order='FPC', diag=FALSE)

ggpairs(data[,-1], lower = list(continuous = "smooth"),
        diag = list(continuous = "density"),
        axisLabels = "show")

ggpairs(data[,-1]^0.25, lower = list(continuous = "smooth"),
        diag = list(continuous = "density"),
        axisLabels = "show")
#there is evidence of non-normality and non-linearity
#although we could attempt to normalize via sqrt transform, this is
#unlikely to fix all.vars'
#linearity also not good.
#It is likely that these issues are the result of the sampling
#occuring over a larger scale than the natural range of the taxa.
#This will result in some species being left skewed and others being
#right skewed. It will also result in non-linearity and the horseshoe
#effect.
data.dca <- decorana(data[,-1])
data.dca

data.std <- data %>%
    dplyr::select(-Sites) %>%
    mutate(across(everything(), function(x) x^0.25)) %>%
    wisconsin()
data.std

## Run PCA - unconstrained axis rotations
## normally we would not scale, but for illustrative purposes...
data.rda <- rda(data.std, scale=TRUE)
## data.rda <- rda(data.std, scale=FALSE)
summary(data.rda, display=NULL)

screeplot(data.rda)
abline(a=1,b=0)

scores(data.rda, choices=1:3, display='sites')
scores(data.rda, choices=1:3, display='species')

## Quick and nasty ordination plots
biplot(data.rda, scaling='species')
biplot(data.rda, scaling='sites')

## Quick and nasty ordination plots
pl<-ordiplot(data.rda)
points(pl, "sites", pch=21, col="red", bg="yellow")
text(pl, "sites", col="red", cex=0.9)
text(pl, "species", col="blue", cex=0.9)
## line(pl, "sites")


## ggvegan provides ways to use ggplot
## library(devtools)
## devtools::install_github("gavinsimpson/ggvegan")
## library(ggvegan)
autoplot(data.rda)
autoplot(data.rda) + theme_bw()
autoplot(data.rda,geom='text') + theme_bw()

## Alternatively, we can extract the scores from the PCA/RDA etc
## into a dataframe.
## - fortify() is a general method used to extract/convert data into
##   a form that is suitable for ggplot.
## - the ggvegan package defines some fortify methods specifically
##   for analyses created with the vegan package

data.rda.scores <- data.rda %>%
    fortify()
data.rda.scores


#######################################################################
g <-
    ggplot(data = NULL, aes(y=PC2, x=PC1)) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_vline(xintercept=0, linetype='dotted') +
    geom_point(data=data.rda.scores %>% filter(Score=='sites')) +
    geom_text(data=data.rda.scores %>% filter(Score=='sites'),
              aes(label=Label), hjust=-0.2) +
    geom_segment(data=data.rda.scores %>% filter(Score=='species'),
                 aes(y=0, x=0, yend=PC2, xend=PC1),
                 arrow=arrow(length=unit(0.3,'lines')), color='red') +
    ## geom_text(data=data.rda.scores %>% filter(Score=='species'),
    ##           aes(y=PC2*1.1, x=PC1*1.1, label=Label), color='red') +
    geom_text_repel(data=data.rda.scores %>% filter(Score=='species'),
              aes(y=PC2*1.1, x=PC1*1.1, label=Label), color='red') +
    theme_bw()
g

# Nice axes titles
eig <- eigenvals(data.rda)

paste(names(eig[2]), sprintf('(%0.1f%% explained var.)', 100 * eig[2]/sum(eig)))
g <- g +
    scale_y_continuous(paste(names(eig[2]), sprintf('(%0.1f%% explained var.)',
                                                    100 * eig[2]/sum(eig))))+
    scale_x_continuous(paste(names(eig[1]), sprintf('(%0.1f%% explained var.)',
                                                    100 * eig[1]/sum(eig))))

#put a circle
circle.prob <- 0.68
circle.prob <- 0.95
circle.prob <- 0.95
r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(data.rda$CA$u[,1:2]^2))^(1/4)
theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
circle <- data.frame(PC1 = r * cos(theta), PC2 = r * sin(theta))
g <- g + geom_path(data = circle, aes(y=PC2,x=PC1), color = muted('white'), size = 1/2, alpha = 1/3)
g
######################################################################


## data.sites.scores <- data.rda %>%
##     scores(display='sites') %>%
##     as.data.frame() %>%
##     bind_cols(data)

## data.species.scores <- data.rda %>%
##     scores(display='species') %>%
##     as.data.frame() %>%
##     mutate(Species=rownames(.))


## ggplot() +
##   geom_point(data=data.sites.scores,  aes(y=PC2,  x=PC1))

## g <- ggplot() +
##     geom_segment(data=NULL, aes(y=-Inf,x=0,yend=Inf,xend=0),
##                  linetype='dotted')+
##     geom_segment(data=NULL, aes(y=0,x=-Inf,yend=0,xend=Inf),
##                  linetype='dotted')
## g <- g +
##     geom_point(data=data.sites.scores, aes(y=PC2,x=PC1))+
##     geom_text(data=data.sites.scores, aes(y=PC2,x=PC1, label=Sites,hjust=-0.2),
##               show.legend=FALSE) +
##     geom_segment(data=data.species.scores, aes(y=0,x=0,yend=PC2,xend=PC1),
##                  arrow=arrow(length=unit(0.3,'lines')), color='red')+
##     theme_bw()
## g




## OR
#data.species.scores <- data.rda %>%
#    fortify(display='species')
#data.site.scores <- data.rda %>%
#    fortify(display='sites')
#g<-ggplot() +
#    geom_segment(data=NULL, aes(y=-Inf,x=0,yend=Inf,xend=0),
#                 linetype='dotted')+
#    geom_segment(data=NULL, aes(y=0,x=-Inf,yend=0,xend=Inf),
#                linetype='dotted')+
#   geom_point(data=data.species.scores, aes(y=PC2,x=PC1)) +
#   geom_text(data=data.site.scores, aes(y=PC2,x=PC1, label=Label,hjust=-0.2),
#             show.legend=FALSE)+
#    geom_segment(data=data.species.scores, aes(y=0,x=0,yend=PC2,xend=PC1),
#                 arrow=arrow(length=unit(0.3,'lines')), color='red')+
#    theme_bw()

## #put a circle
## circle.prob <- 0.68
## r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(data.rda$CA$u[,1:2]^2))^(1/4)
## theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
## circle <- data.frame(PC1 = r * cos(theta), PC2 = r * sin(theta))
## g <- g + geom_path(data = circle, aes(y=PC2,x=PC1),color = muted('white'), size = 1/2, alpha = 1/3)

## g

## ##Put on some nice species labels
## ## First filter out all those outside the circle
## data.species.scores.sub <- data.species.scores %>%
##   mutate(Length=sqrt(PC1^2 + PC2^2)) %>%
##   filter(Length>r)
## hjust <- ifelse(data.species.scores.sub$PC1>0,0,1)
## vjust <- ifelse(data.species.scores.sub$PC2>0,0,1)
## g <- g + geom_text(data=data.species.scores.sub, aes(y=PC2,x=PC1, label=Species),
##                    hjust=hjust, vjust=vjust, color='red')
## g

## #Put on some nice axes titles
## eig <- eigenvals(data.rda)
## paste(names(eig[2]), sprintf('(%0.1f%% explained var.)', 100 * eig[2]/sum(eig)))
## g <- g + scale_y_continuous(paste(names(eig[2]), sprintf('(%0.1f%% explained var.)', 100 * eig[2]/sum(eig))))+
##     scale_x_continuous(paste(names(eig[1]), sprintf('(%0.1f%% explained var.)', 100 * eig[1]/sum(eig))))
## g




## library(devtools)
## devtools::install_github("gavinsimpson/ggvegan")
## library(ggvegan)
## autoplot(data.rda) + theme_bw()
## autoplot(data.rda, arrows=TRUE, display='species') + theme_bw()
## autoplot(data.rda,geom='text') + theme_bw()

## Envfit
library(car)
vif(lm(1:nrow(enviro)~pH+Slope+Pressure+Altitude+Substrate, data=enviro))
vif(lm(1:nrow(enviro)~pH+Slope+Altitude+Substrate, data=enviro))
lm(1:nrow(enviro)~pH+Slope+Pressure+Altitude+Substrate, data=enviro) %>% vif
#take out either Altitude or Pressure
vif(lm(1:nrow(enviro)~pH+Slope+Altitude+Substrate, data=enviro))

Xmat <- model.matrix(~-1+pH+Slope+Altitude+Substrate, enviro)
Xmat <- model.matrix(~-1+pH+Slope+Altitude+Substrate+Pressure, enviro)
data.env <- envfit(data.rda, env=Xmat)
data.env

########
data.env.scores <- data.env %>% fortify()
g <- g + 
    geom_segment(data=data.env.scores,
                 aes(y=0, x=0, yend=PC2, xend=PC1),
                 arrow=arrow(length=unit(0.3,'lines')), color='blue') +
    geom_text(data=data.env.scores,
              aes(y=PC2*1.1, x=PC1*1.1, label=Label), color='blue')
g
##############

## data.env.scores <- data.env %>%
##     scores(display='vector') %>%
##     as.data.frame %>%
##     mutate(Effect=rownames(.))

## hjust <- ifelse(data.env.scores$PC1>0,0,1)
## vjust <- ifelse(data.env.scores$PC2>0,0,1)
## g<- g + geom_segment(data=data.env.scores, aes(y=0,x=0,yend=PC2,xend=PC1),
##                      arrow=arrow(length=unit(0.3,'lines')), color='blue')+
##     geom_text(data=data.env.scores, aes(y=PC2,x=PC1, label=Effect),
##               hjust=hjust, vjust=vjust, color='blue')
## g


##extract PC's to use as responses
PC1 <- data.rda.scores %>% filter(Score=='sites') %>% pull(PC1)
PC2 <- data.rda.scores %>% filter(Score=='sites') %>% pull(PC2)

lm(PC1 ~ pH+Slope+Altitude+Substrate, data=enviro) %>%
    summary()
lm(PC2 ~ pH+Slope+Altitude+Substrate, data=enviro) %>%
    summary()
## ----end
## ---- RDA
data.rda <- rda(data.std ~ scale(pH)+scale(Slope)+scale(Altitude)+Substrate,
                data=enviro, scale=FALSE)
summary(data.rda, display=NULL)
## Mention conditioning

#data.rda <- rda(data.stnd ~ scale(pH)+scale(Slope)+Condition(scale(Altitude))+Substrate,
#                data=enviro, scale=FALSE)
vif.cca(data.rda)

## Goodness of fit test
## proportion of inertia accounted for by species
## up to chosen axes
## proportion can be assessed by either species or sites
## (depending on argument display)
## THis is the cumulative proportion of variance explained by each axis
goodness(data.rda)
## Proportion of inertia explained by constrained and unconstrained
inertcomp(data.rda)
## or proportionally
inertcomp(data.rda, proportional = TRUE)


#overall test
anova(data.rda)
anova(data.rda, by='axis')
anova(data.rda, by='margin')
#anova(data.rda, by='margin', scope="Altitude")

## see the regression coefficients
coef(data.rda)

RsquareAdj(data.rda)

screeplot(data.rda)

autoplot(data.rda, geom='text')


#######################################################################
data.rda.scores <- data.rda %>% fortify()
data.rda.scores
ggplot(data = NULL, aes(y=PC2, x=PC1)) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_vline(xintercept=0, linetype='dotted') +
    geom_point(data=data.rda.scores %>% filter(Score=='sites')) +
    geom_text(data=data.rda.scores %>% filter(Score=='sites'),
              aes(label=Label), hjust=-0.2) +
      MANAGEMENTjs(label=LabMANAGEMENT data=dunej=-0.2) +
    geom_segment(data=data.rda.scores %>% filter(Score=='species'),
                 aes(y=0, x=0, yend=PC2, xend=PC1),
                 arrow=arrow(length=unit(0.3,'lines')), color='red') +
    geom_text(data=data.rda.scores %>% filter(Score=='species'),
              aes(y=PC2*1.1, x=PC1*1.1, label=Label), color='red') ->
    g
g

# Nice axes titles
eig <- eigenvals(data.rda)
paste(names(eig[2]), sprintf('(%0.1f%% explained var.)', 100 * eig[2]/sum(eig)))
g <- g + scale_y_continuous(paste(names(eig[2]), sprintf('(%0.1f%% explained var.)',
                                                         100 * eig[2]/sum(eig))))+
    scale_x_continuous(paste(names(eig[1]), sprintf('(%0.1f%% explained var.)',
                                                    100 * eig[1]/sum(eig))))

#put a circle
## circle.prob <- 0.68
circle.prob <- 0.95
r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(data.rda$CA$u[,1:2]^2))^(1/4)
theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
circle <- data.frame(PC1 = r * cos(theta), PC2 = r * sin(theta))
g <- g + geom_path(data = circle, aes(y=PC2,x=PC1),color = muted('white'), size = 1/2, alpha = 1/3)
g


g <- g + 
    geom_segment(data=data.rda.scores %>% filter(Score=='biplot'),
                 aes(y=0, x=0, yend=RDA2, xend=RDA1),
                 arrow=arrow(length=unit(0.3,'lines')), color='blue') +
    geom_text(data=data.rda.scores %>% filter(Score=='biplot'),
              aes(y=RDA2*1.1, x=RDA1*1.1, label=Label), color='blue')
g


## ----end
## ---- CA spiders

spider.std <- spider$abund %>%
    decostand(method="total",MARGIN=2)
spider.std

spider.std %>% cor %>%
    corrplot(diag=FALSE, order='FPC')

spider.ca <- cca(spider.std, scale=FALSE)
summary(spider.ca, display=NULL)
#anova(spider.ca)

screeplot(spider.ca)
sum(eigenvals(spider.ca))/length(eigenvals(spider.ca))
eigenvals(spider.ca)/sum(eigenvals(spider.ca))
g <- autoplot(spider.ca, geom='text')
spider.env <- envfit(spider.ca, env=X)
spider.env
g + autoplot(spider.env)



## spider.dca <- decorana(spider.abund)
## summary(spider.dca, display="none")
## ordiplot(spider.dca)
## autoplot(spider.dca,  geom="text")
## ----end
## ---- CA
#convert to frequencies (of row and columns)
# MARGIN=2 indicates columns
data.std <- data %>%
    dplyr::select(-Sites) %>% 
    decostand(method="total",MARGIN=2)
data.std

data.std %>% cor %>%
    corrplot(diag=FALSE)

data.std %>% cor %>%
    corrplot(diag=FALSE, order='FPC')

data.ca <- cca(data.std, scale=FALSE)
summary(data.ca, display=NULL)
#anova(data.ca)

screeplot(data.ca)
sum(eigenvals(data.ca))/length(eigenvals(data.ca))
eigenvals(data.ca)/sum(eigenvals(data.ca))
plot(data.ca, scaling='species')
Xmat <- model.matrix(~-1+pH+Slope+Altitude+Substrate, enviro)
data.env <- envfit(data.ca, env=Xmat)
data.env
autoplot(data.env)

##                                         #extract CA's to use as responses
## data.sites.scores <- as.data.frame(scores(data.ca, display='sites'))
## head(data.sites.scores)
## data.sites.scores <- data.frame(data.sites.scores, data)
## data.species.scores <- as.data.frame(scores(data.ca, display = 'species'))
## head(data.species.scores)
## data.species.scores$Species <- colnames(data[,-1])
## head(data.species.scores)

autoplot(data.ca, geom='text')

data.ca.scores <- data.ca %>% fortify()
CA1 <- data.ca.scores %>% filter(Score =='sites') %>% pull(CA1)
CA2 <- data.ca.scores %>% filter(Score =='sites') %>% pull(CA2)
summary(lm(CA1 ~ pH+Slope+Altitude+Substrate, data=enviro))
summary(lm(CA2 ~ pH+Slope+Altitude+Substrate, data=enviro))

## ----end
## ---- CCA
data.cca <- cca(data.std~pH + Altitude + Substrate + Slope, data=enviro, scale=FALSE)

summary(data.cca, display=NULL)
anova(data.cca)

plot(data.cca)
plot(data.cca, scaling='species')
plot(data.cca, scaling='sites')
autoplot(data.cca)

data.cs <- data.cca$CA$u[,1]
data.rs <- data.cca$CA$v[,1]
data3 <- data.std[order(data.cs),order(data.rs)]
data3 <- data.frame(Sites=rownames(data3),data3)
data4 = data3 %>%
    gather(key=Species, value=value, -Sites) %>%
    mutate(Species=factor(Species, levels=unique(Species)),
           Sites=factor(Sites, levels=unique(Sites))) 

vif.cca(data.cca)
#overall test
anova(data.cca)
anova(data.cca, by='axis')
anova(data.cca, by='margin')
#anova(data.cca, by='margin', scope="pH")

coef(data.cca)

RsquareAdj(data.cca)

screeplot(data.cca)
int <- data.cca$tot.chi/length(data.cca$CA$eig)
abline(h=int)

## ----end


## ---- MDS spiders
## MUST READ IN THIS WAY..

spider.mds <- metaMDS(spider.abund, k=2,  plot=TRUE)
spider.mds

spider.std <- wisconsin(spider.abund^0.25)

#apply(macnally.std[,c(-1,-2)],2,max)
#apply(macnally.std[,c(-1,-2)],2,var, na.rm=TRUE)

#vegdist(macnally[,-1], method='bray')

spider.dist <- vegdist(spider.std,"bray")

spider.mds <- metaMDS(spider.std, k=2, plot=TRUE)
spider.mds <- metaMDS(spider.dist, k=2, plot=TRUE)
spider.mds <- metaMDS(spider.abund, k=2)


wascores(spider.mds$points, spider.abund)


spider.mds$stress

stressplot(spider.mds)

plot(macnally.mds)

## autoplot(macnally.mds)
spider.mds.scores <- spider.mds %>%
        fortify()

spider.env <- spider$x
g <-
    ggplot(data = NULL, aes(y=NMDS2, x=NMDS1)) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_vline(xintercept=0, linetype='dotted') +
    geom_point(data=spider.mds.scores %>%
                   filter(Score=='sites'),
               aes(color=spider.env$soil.dry)) +
    geom_text(data=spider.mds.scores %>%
                  filter(Score=='sites'),
              aes(label=Label,
                  color=spider.env$soil.dry), hjust=-0.2) +
    scale_colour_binned(type = "viridis")
g

spider.env <- spider.env %>%
    mutate(fSoil = cut(soil.dry, breaks = c(0, 1, 2, 4)))
g <-
    ggplot(data = NULL, aes(y=NMDS2, x=NMDS1)) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_vline(xintercept=0, linetype='dotted') +
    geom_point(data=spider.mds.scores %>%
                   filter(Score=='sites'),
               aes(color=spider.env$fSoil)) +
    geom_text(data=spider.mds.scores %>%
                  filter(Score=='sites'),
              aes(label=Label,
                  color=spider.env$fSoil), hjust=-0.2) 
    scale_colour_viridis_c()
g


spider.env <- spider.env %>%
    mutate(fFallen = cut(fallen.leaves, breaks = c(-1, 1, 3, 5)))
g <-
    ggplot(data = NULL, aes(y=NMDS2, x=NMDS1)) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_vline(xintercept=0, linetype='dotted') +
    geom_point(data=spider.mds.scores %>%
                   filter(Score=='sites'),
               aes(color=spider.env$fFallen)) +
    geom_text(data=spider.mds.scores %>%
                  filter(Score=='sites'),
              aes(label=Label,
                  color=spider.env$fFallen), hjust=-0.2) 
    scale_colour_viridis_c()
g


g + ggforce::geom_mark_ellipse(data=spider.mds.scores %>%
                                   filter(Score=='sites'),
                               aes(y=NMDS2, x=NMDS1, fill=spider.env$fFallen),
                               expand=0) 
g + ggforce::geom_mark_hull(data=spider.mds.scores %>%
                                filter(Score=='sites'),
                            aes(y=NMDS2, x=NMDS1, fill=spider.env$fFallen),
                            expand=0, concavity = 10) 

Xmat <- model.matrix(~soil.dry + bare.sand + fallen.leaves + moss +
                         herb.layer + reflection, data=spider.env)
envfit <- envfit(spider.mds, env=Xmat)
envfit <- envfit(spider.mds, env=spider.env)
envfit


spider.env.scores <- envfit %>% fortify()
g <- g + 
    geom_segment(data=spider.env.scores,
                 aes(y=0, x=0, yend=NMDS2, xend=NMDS1),
                 arrow=arrow(length=unit(0.3,'lines')), color='blue') +
    geom_text(data=spider.env.scores,
              aes(y=NMDS2*1.1, x=NMDS1*1.1, label=Label), color='blue')
g


## ----end
## ---- bioenv spider

spider.dist <- vegdist(wisconsin(spider.abund^0.25),"bray")
spider.bioenv <- bioenv(spider.dist,
                      decostand(spider.env,"standardize"))
spider.bioenv
## ----end








## ---- MVABUND

combined.data <- cbind(data, enviro)
names(combined.data)
mva = mvabund(data[,-1])

meanvar.plot(mva)
plot(mva)
X = enviro$Substrate
## enviro = enviro %>% mutate(ph=cut(pH, breaks=c(0,2,4,6,8,10)))

data.mod <- manyglm(mva~scale(pH) + scale(Altitude) + Substrate + scale(Slope),
                    family=poisson(link='log'), data=enviro)

plot(data.mod)

data.mod <- manyglm(mva~scale(pH) + scale(Altitude) + Substrate + scale(Slope),
                    family='negative.binomial', data=enviro)
plot(data.mod)
data.mod
anova(data.mod, test='LR')
##compute analysis of deviance table
## To compute P values, the default is to use PIT (Probability Integral Transform residuals) or "PIT-trap"
## which have been found to five the most reliable Type I error rates.
## All methods are robust to mean/variance relationship as well as the correlation between variables.
## test = 'LR' - likelihood ratio - these are better if the estimated abundances are close to zero for
## either count or binomial data.
anova(data.mod, test='LR')
## LR has good properties, but is only available for cor.type="I"

## The above assumes that the individual responses are not correlated to one another - although the
## permutations do handle this..
## If the number of species is small compared to the number of rows, it is recommended to account for this
## correlation with cor.type="R".
## However, if the number of species are large relative to the number of rows, they can be
## unstable and expensive and low power
anova(data.mod, cor.type = 'R')
## If the number of species are not small relative to the number of rows, then 
## it is recomemnded if possible to specify cor.type="shrink".  Whilst this is much more
## computationally expensive, they are a good comprimise between "I" and "R"
anova(data.mod, cor.type = 'shrink')
## We can also explore the individal univariate tests.
anova(data.mod, p.uni='adjusted')
#the following is better if number of species is much larger than the number of rows.
summary(data.mod, test="LR")
## summary(data.mod, resamp="residual")
## summary(data.mod, resamp="perm.resid")
#summary(data.mod, resamp="monte.carlo", test="wald", nBoot=300) 
#plot(data.mod)


inverts.mva <- mvabund(inverts)
inverts.mglmP <- manyglm(inverts.mva ~ TREATMENT * WEEK, data = brink, family = 'poisson')
plot(inverts.mglmP) 
inverts.mglmNB <- manyglm(inverts.mva ~ TREATMENT * WEEK, data = brink, family = 'negative.binomial')
plot(inverts.mglmNB) 
control <- how(within = Within(type = 'none'),
               Plots(strata = brink$DITCH, type = 'free'),
               nperm = 50)
permutations <- shuffleSet(nrow(inverts.mva), control = control)
inverts.mglmNB2 <- manyglm(inverts.mva ~ TREATMENT + WEEK, 
                                data = brink, family = 'negative.binomial')
inverts_aov <- anova(inverts.mglmNB, inverts.mglmNB2, 
                     bootID = permutations,  
                     p.uni = 'adjusted', test = 'LR') 
inverts_aov 

## Compare to model without any treatment - so test for effect of treatment
inverts.mglmNB3 <- manyglm(inverts.mva ~ WEEK, data = brink, 
                       family = 'negative.binomial')
inverts_aov2 <- anova(inverts.mglmNB, inverts.mglmNB3 , bootID = permutations,  
      p.uni = 'adjusted', test = 'LR') 
inverts_aov2 



mod_pt <- NULL
for (i in levels(brink$WEEK)) {
    brink.sub <- brink %>% filter(WEEK == i)
    inverts.sub <- brink.sub %>% dplyr::select(-TREATMENT, -WEEK, -DITCH) %>%
        mvabund()
    ## model
    ##mod_pt[[i]]$mod <- manyglm(inverts.sub ~ TREATMENT, data = brink.sub)
    mod <- manyglm(inverts.sub ~ TREATMENT, data = brink.sub)
    aov <- anova(mod, nBoot = 100, 
                 p.uni = 'adjusted', test = 'LR', show.time = "none")
    sum <- summary(mod, nBoot = 100, 
                   p.uni = 'adjusted', test = 'LR')
    
    P <- c(community = aov$table[2,4],
           aov$uni.p[2,])
    mod_pt[[i]] <- list(mod = mod, aov=aov, P=P)
}
dd <- do.call('rbind', lapply(mod_pt, function(x) x$P)) %>%
    as.data.frame() %>% 
    rownames_to_column(var = 'WEEK')
dd

## purrr alternative
library(purrr)
d = bind_cols(inverts = inverts.mva, brink %>% dplyr::select(TREATMENT, WEEK, DITCH))
dd <- d %>% group_by(WEEK) %>%
    nest() %>%
    mutate(mod = purrr::map(data, function(x) {
        manyglm(inverts ~ TREATMENT, data=x)
    })) %>% 
    mutate(aov = purrr::map(mod, function(x) {
        anova(x, nBoot=100, p.uni = 'adjusted', test = 'LR', show.time = 'none')
    })) %>%
    mutate(sum = purrr::map(mod, function(x) {
        summary(x, nBoot=100, p.uni = 'adjusted', test = 'LR')
    })) %>%
    mutate(P = purrr::map(aov, function(x) {
        c(Community = x$table[2,4], x$uni.p[2,])
        }))
dd %>% dplyr::select(WEEK, P) %>% unnest_wider(P)

g <- 
    dd %>% mutate(Deviance = purrr::map(aov, function(x) {
        x$uni.test[2,]
    })) %>%
    dplyr::select(WEEK, Deviance) %>% 
    unnest_wider(Deviance) %>%
    pivot_longer(cols=-WEEK) %>%
    ungroup %>%
    mutate(name = forcats::fct_reorder(name, value, 'sum', .desc = TRUE)) %>%
    ggplot(aes(y=value, x=as.numeric(as.character(WEEK)), fill=name)) +
    geom_area() +
    geom_vline(aes(xintercept = 0)) 
g





## data(pyrifos)
## take <- c('binitent', 'olchaeta', 'caenhora', 'cloedipt', 
##           'chaoobsc', 'gammpule', 'libellae', 'agdasphr')
## abu <- pyrifos[ , names(pyrifos) %in% take]
## abu <- round((exp(abu) - 1)/10)
## head(abu)
## env <- data.frame(time = rep(c(-4, -1, 0.1, 1, 2, 4, 8, 12, 15, 19, 24), each = 12),
##                   treatment = rep(c(0.1, 0, 0, 0.9, 0, 44, 6, 0.1, 44, 0.9, 0, 6), 11),
##                   replicate = gl(12, 1, length=132))
## head(env)
## abu <- mvabund(abu)
## mod_pois <- manyglm(abu ~ factor(treatment) * factor(time), 
##                     data = env, family = 'poisson')
## plot(mod_pois) 
## mod_nb <- manyglm(abu ~ factor(treatment) * factor(time), 
##                   data = env, family = 'negative.binomial')

## plot(mod_nb)
## control <- how(within = Within(type = 'none'),
##                plots = Plots(strata = env$replicate, type = 'free'),
##                nperm = 50)
## permutations <- shuffleSet(nrow(abu), control = control)
## mod_nb_nointeraction <- manyglm(abu ~ factor(treatment) + factor(time), 
##                                 data = env, family = 'negative.binomial')
## mod_nb_aov <- anova(mod_nb, mod_nb_nointeraction, 
##                     bootID = permutations,  
##                     p.uni = 'adjusted', test = 'LR') 
## mod_nb_aov

## ## Compare to model without any treatment - so test for effect of treatment
## mod_nb_null <- manyglm(abu ~ factor(time), data = env, 
##                        family = 'negative.binomial')
## mod_treat_aov <- anova(mod_nb, mod_nb_null , bootID = permutations,  
##       p.uni = 'adjusted', test = 'LR') 
## mod_treat_aov

## #Effect of treatment at each time
## mod_pt <- NULL
## for (i in levels(factor(env$time))) {
##   take_abu <- abu[env$time == i, ]
##   take_env <- env[env$time == i, ]
##   # model
##   mod_pt[[i]]$mod <- manyglm(take_abu ~ factor(treatment), data = take_env)
##   mod_pt[[i]]$aov <- anova(mod_pt[[i]]$mod, nBoot = 100, 
##                            p.uni = 'adjusted', test = 'LR', show.time = "none")
##   mod_pt[[i]]$sum <- summary(mod_pt[[i]]$mod, nBoot = 100, 
##                              p.uni = 'adjusted', test = 'LR')
## }
## get_pvals <- function(x){
##   comm <- c(community = x$aov$table[2, 4])
##   spec <- x$aov$uni.p[2, ]
##   c(comm, spec)
## }
## plyr::ldply(mod_pt, get_pvals)

## ## purrr alternative
## d = bind_cols(abu = abu, env)
## dd <- d %>% group_by(time) %>%
##     nest() %>%
##     mutate(mod = map(data, function(x) {
##         manyglm(abu ~ factor(treatment), data=x)
##     })) %>% 
##     mutate(aov = map(mod, function(x) {
##         anova(x, nBoot=100, p.uni = 'adjusted', test = 'LR', show.time = 'none')
##     })) %>%
##     mutate(sum = map(mod, function(x) {
##         summary(x, nBoot=100, p.uni = 'adjusted', test = 'LR')
##     })) %>%
##     mutate(P = map(aov, function(x) {
##         c(Community = x$table[2,4], x$uni.p[2,])
##         }))
## dd %>% dplyr::select(time, P) %>% unnest_wider(P)


## devs <- plyr::ldply(mod_pt, function(x) x$aov$uni.test[2, ])
## plotdf <- reshape::melt(devs, id.vars = '.id')
## ggplot(plotdf, aes(x = as.numeric(as.character(.id)), y = value, fill = variable)) +
##   geom_area(col = 'grey15') + 
##   geom_vline(aes(xintercept = 0)) +
##   theme_bw() +
##   labs(x = 'time', y = 'Deviance (=Effect)')


## ## alternative
## dd %>% mutate(Deviance = map(aov, function(x) {
##     x$uni.test[2,]
## })) %>%
##     dplyr::select(time, Deviance) %>% 
##     unnest_wider(Deviance) %>%
##     pivot_longer(cols=-time) %>%
##     ## mutate(name = forcats::fct_reorder(name, value, max)) %>%
##     ggplot(aes(y=value, x=time, fill=name)) +
##     geom_area()


## ----end


## ---- glmmTMB
data.1 <- data %>%
    pivot_longer(cols = -Sites,
                 names_to = 'Species',
                 values_to = 'Abund')
    
library(glmmTMB)
data.glmmTMB <- glmmTMB(Abund ~ Species + rr(Species + 0|Sites, d = 2),
                          family = nbinom2(),
                          dat = data.1
                          )
summary(data.glmmTMB)
library(DHARMa)
data.resids <- simulateResiduals(data.glmmTMB, plot = TRUE)
testDispersion(data.resids)
data.loadings <- data.glmmTMB$obj$env$report(data.glmmTMB$fit$parfull)$fact_load[[1]] %>%
                                                                           as.data.frame() %>%
                                                                           mutate(Species = colnames(data[,-1]))
fit <-
    ranef(data.glmmTMB)[[1]]$Sites %>%
    mutate(Site = rownames(.))

ggplot(fit, aes(y = SpeciesSp1, x = SpeciesSp10)) +
    geom_text(aes(label = Site)) #+
    geom_text(data = data.loadings, aes(y = V2, x = V1, label = Species), color = 'blue')


data.2 <- data %>%
    cbind(enviro) %>%
    pivot_longer(cols = c(-Sites, -Site, -pH, -Slope, -Pressure, -Altitude, -Substrate),
                 names_to = 'Species',
                 values_to = 'Abund')
    
data.glmmTMB <- glmmTMB(Abund ~ (pH+Slope+Altitude+Substrate) + (0 + pH+Slope+Altitude+Substrate|Species) + (1|Sites),
                        family = nbinom2(),
                        dat = data.2
                        ## control = glmmTMBControl(optimizer = 'optim',
                        ##                          optArgs = 'Nelder-Mead')
                        )
summary(data.glmmTMB)
library(DHARMa)
data.resids <- simulateResiduals(data.glmmTMB, plot = TRUE)
testDispersion(data.resids)
## ----end




## ---- brms
library(brms)

combined.data = cbind(data, enviro)
head(combined.data)
data.form <- bf(mvbind(Sp1,Sp2,Sp3,Sp4,Sp5,Sp6,Sp7,Sp8,Sp9,Sp10) ~ pH+Slope+Altitude+Substrate,
                family=negbinomial(link='log'))
a <- brm(data.form,
         data=combined.data)
summary(a)

data.form <- bf(mvbind(Sp1,Sp2,Sp3,Sp4,Sp5,Sp6,Sp7,Sp8,Sp9,Sp10) ~ Altitude,
                family=gaussian())
b <- brm(data.form,
         data=combined.data)


## ----end


## ---- MDS macnally
## MUST READ IN THIS WAY..
macnally <- read.csv('../public/data/macnally_full.csv',strip.white=TRUE)
head(macnally)
macnally[1:5,1:5]




apply(macnally[,c(-1)],2,mean, na.rm=TRUE)
apply(macnally[,c(-1)],2,max)
apply(macnally[,c(-1)],2,sum)
apply(macnally[,c(-1)],2,var, na.rm=TRUE)


macnally.mds <- metaMDS(macnally[,-1], k=2,  plot=TRUE)
macnally.mds

library(vegan)
macnally.std <- wisconsin(macnally[,c(-1)]^0.25)

#apply(macnally.std[,c(-1,-2)],2,max)
#apply(macnally.std[,c(-1,-2)],2,var, na.rm=TRUE)

#vegdist(macnally[,-1], method='bray')

macnally.dist <- vegdist(macnally.std,"bray")

macnally.mds <- metaMDS(macnally.std, k=2, plot=TRUE)
macnally.mds <- metaMDS(macnally.dist, k=2, plot=TRUE)
macnally.mds <- metaMDS(macnally[,-1], k=2)


wascores(macnally.mds$points, macnally[, -1])


macnally.mds$stress

stressplot(macnally.mds)

plot(macnally.mds)

## autoplot(macnally.mds)
macnally.mds.scores <- macnally.mds %>%
        fortify() %>%
    full_join(macnally %>% add_rownames(var='Label'))

g <-
    ggplot(data = NULL, aes(y=NMDS2, x=NMDS1)) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_vline(xintercept=0, linetype='dotted') +
    geom_point(data=macnally.mds.scores %>% filter(Score=='sites'),
               aes(color=HABITAT)) +
    geom_text(data=macnally.mds.scores %>% filter(Score=='sites'),
              aes(label=Label, color=HABITAT), hjust=-0.2) 
    ## geom_segment(data=macnally.mds.scores %>% filter(Score=='species'),
    ##              aes(y=0, x=0, yend=NMDS2, xend=NMDS1),
    ##              arrow=arrow(length=unit(0.3,'lines')), color='red') +
    ## geom_text(data=macnally.mds.scores %>% filter(Score=='species'),
    ##           aes(y=NMDS2*1.1, x=NMDS1*1.1, label=Label), color='red') 
g

# Nice axes titles
## this is not available for MDS as it is not eigen based...



## species <- wascores(macnally.mds$points,  macnally[, -1]) %>%
##   as.data.frame %>%
##   mutate(Species = row.names(.))

## #g=autoplot(macnally.mds)
## #g

## macnally.sites.scores <- macnally.mds %>%
##     scores(display='sites') %>%
##     as.data.frame() %>%
##     bind_cols(macnally)
## ## If used raw data in metaMDS
## macnally.species.scores <- macnally.mds %>%
##     scores(display='species') %>%
##     as.data.frame() %>%
##     mutate(Species = rownames(.))
## ## If used distance matrix in MDS
## macnally.species.scores <- wascores(macnally.mds$points,  macnally[, -1]) %>%
##     as.data.frame() %>%
##   mutate(Species = rownames(.),
##          NMDS1=MDS1,
##          NMDS2=MDS2)
## head(macnally.species.scores)

## g<-ggplot() +
##     geom_point(data= macnally.sites.scores, aes(y=NMDS2,x=NMDS1, color=HABITAT))+
##                                         #geom_text(data=macnally.sites.scores, aes(y=NMDS2,x=NMDS1, label=SITE,hjust=-0.2,color=HABITAT), show_guide=FALSE)+
##      geom_point(data=macnally.species.scores, aes(y=NMDS2,x=NMDS1), color='grey') +
##     geom_text(data=macnally.species.scores, aes(y=NMDS2,x=NMDS1, label=Species,
##                                                 hjust=-0.2),show.legend=FALSE, color='grey') +
##     geom_segment(data=NULL, aes(y=-Inf,x=0,yend=Inf,xend=0), linetype='dotted')+
##     geom_segment(data=NULL, aes(y=0,x=-Inf,yend=0,xend=Inf), linetype='dotted')
## #    geom_segment(data=macnally.species.scores, aes(y=0,x=0,yend=NMDS2,xend=NMDS1), arrow=arrow(length=unit(0.3,'lines')))
## g


## macnally.hull = macnally.mds.scores %>%
##     filter(Score=='sites') %>%
##   group_by(HABITAT) %>%
##   slice(chull(NMDS1, NMDS2))

## g <- g + geom_polygon(data=macnally.hull, aes(y=NMDS2,x=NMDS1, fill=HABITAT), alpha=0.2)+
##     theme_classic()
## g

g + ggforce::geom_mark_ellipse(data=macnally.mds.scores %>% filter(Score=='sites'),
                      aes(y=NMDS2, x=NMDS1, fill=HABITAT), expand=0) 
g + ggforce::geom_mark_hull(data=macnally.mds.scores %>% filter(Score=='sites'),
                      aes(y=NMDS2, x=NMDS1, fill=HABITAT), expand=0) 
g + ggforce::geom_mark_hull(data=macnally.mds.scores %>% filter(Score=='sites'),
                      aes(y=NMDS2, x=NMDS1, fill=HABITAT), expand=0, concavity = 10) 

    ## geom_segment(data=centroids.long(macnally.sites.scores, grouping=HABITAT))

#ordiplot(macnally.mds, display="sites", type="n")
#text(macnally.mds,lab=rownames(macnally), col=as.numeric(macnally$HABITAT))

## habitat <- model.matrix(~-1+macnally$HABITAT)
## colnames(habitat) <-gsub("macnally\\$HABITAT","",colnames(habitat))
## envfit <- envfit(macnally.mds, env=habitat)
## envfit



## ordiplot(macnally.mds, display="sites", type="n")
## text(macnally.mds,lab=rownames(macnally), col=as.numeric(macnally$HABITAT))
## text(macnally.mds,lab=macnally$HABITAT, col=as.numeric(macnally$HABITAT))

Xmat <- model.matrix(~-1+HABITAT, data=macnally)
colnames(Xmat) <-gsub("HABITAT","",colnames(Xmat))
envfit <- envfit(macnally.mds, env=Xmat)
envfit


data.env.scores <- envfit %>% fortify()
g <- g + 
    geom_segment(data=data.env.scores,
                 aes(y=0, x=0, yend=NMDS2, xend=NMDS1),
                 arrow=arrow(length=unit(0.3,'lines')), color='blue') +
    geom_text(data=data.env.scores,
              aes(y=NMDS2*1.1, x=NMDS1*1.1, label=Label), color='blue')
g


plot(envfit, col="gray")

plot(envfit, col="gray")


vegdist(data.stnd)

macnally.dist <- vegdist(macnally[,-1], 'bray')
## bioenv(macnally.dist, macnally$HABITAT)
##                       decostand(macnally$,"standardize"))

## adonis(macnally.dist ~ HABITAT, data=macnally)
## mm <-  model.matrix(~HABITAT, data=macnally)
## adonis(macnally.dist ~ mm)


## Similarity percentage: decomposition of Bray-Curtis index
## into contributions of each species to average between-group
## dissim.
## pairwise
## displays most important species for each pair.
## Caused by variation in species abundances and only partly
## by differences among groups.
## So driven by most abundant species.
simper(macnally.std, macnally$HABITAT)
## betadisp = PERMDISP2
## analysis of multivariate homogeneity of group dispersions (variances) - test for homogeneity of variance
## to test whether the dispersions (variances) of one or more groups
## are differnt, the distances of group memebers to the group centroid
## are subject to ANOVA.
macnally.disp <- betadisper(macnally.dist, macnally$HABITAT)
boxplot(macnally.disp)
plot(macnally.disp)
anova(macnally.disp)
permutest(macnally.disp, pairwise = TRUE)
TukeyHSD(macnally.disp)
## now with bias correction
## use the spatial median rather than the centroid (default)
## if the median or centroid is calculated from the same data as then used to calculate dispersion,
## dispersion will be biased downwards (less dispersion)
## adjust for small sample bias (default is false)
## can also be biased by unequal sample sizes
## good idea to correct for bias when sample sizes are small or uneven
macnally.disp <- betadisper(macnally.dist, macnally$HABITAT, type="median",bias.adjust = TRUE)
boxplot(macnally.disp)
plot(macnally.disp)
anova(macnally.disp)
permutest(macnally.disp, pairwise = TRUE)
TukeyHSD(macnally.disp)
## ----end

## ---- dune
dune <- read_csv('../public/data/dune.csv', trim_ws=TRUE)
dune <- dune %>% mutate(MANAGEMENT=factor(MANAGEMENT,  levels=c("NM","BF","HF","SF"))) %>%
  as.data.frame()
#dune <- read.csv('../downloads/data/dune.csv')
dune

#species means
apply(dune[,-1],2, mean, na.rm=TRUE)
#species maximums
apply(dune[,-1],2, max)
#species sums
apply(dune[,-1],2, sum, na.rm=TRUE)
#species variance
apply(dune[,-1],2, var, na.rm=TRUE)




library(vegan)
dune.dist <- vegdist(wisconsin(dune[,-1]^0.25), "bray")

dune.mds = metaMDS(dune.dist, k=2)
plot(dune.mds, type="text", display="sites" )

dune.mds$species <- wascores(dune.mds$points, dune[,-1], expand = TRUE)
pl <- vegan::ordiplot(dune.mds, type = "none")
points(pl, "sites", pch=21, col="red", bg="yellow")
text(pl, "species", col="blue", cex=0.9)
ordihull(dune.mds, dune$MANAGEMENT, col=as.numeric(dune$MANAGEMENT))
autoplot(dune.mds, geom=c('text'))

#ordiplot(dune.mds)
#ordihull(dune.mds, dune$MANAGEMENT, col=as.numeric(dune$MANAGEMENT))

dune.adonis<-adonis2(dune.dist~MANAGEMENT,  data=dune)
dune.adonis

#plot(isomap(dune.dist, k=3), new=TRUE)


management <-factor(dune$MANAGEMENT, levels=c("NM","BF","HF","SF"))
mm <- model.matrix(~management)
mm <- model.matrix(~MANAGEMENT, data=dune)
head(mm)
colnames(mm) <-gsub("MANAGEMENT","",colnames(mm))
mm <- data.frame(mm)
dune.adonis<-adonis2(dune.dist~BF+HF+SF, data=mm,
                    perm=9999)
dune.adonis
dune.adonis<-adonis2(dune.dist~BF+HF+SF, data=mm,
                    perm=9999)
dune.adonis

library(pairwiseAdonis)
pairwise.adonis(dune.dist, dune$MANAGEMENT)

library(EcolUtils)
adonis.pair(dune.dist, dune$MANAGEMENT, nper = 10000)

dune.simper=simper(dune[,-1], dune[,1], permutations = 999)
summary(dune.simper)

## Multiple Response Permutation Procedure (MRPP) provides a test of
##      whether there is a significant difference between two or more
##      groups of sampling units. Function ‘meandist’ finds the mean
##      within and between block dissimilarities.
## Compares dissim within and between groups
## delta: weighted mean of the within-group means of pairwise dissim
## among sampling units.

dune.mrpp = mrpp(dune.dist, dune[,1], permutations=999)
dune.mrpp
hist(dune.mrpp$boot.deltas)
# Chance corrected within-group agreement = 1-Obs delta / exp delta
dune.meandist = meandist(dune.dist, dune[,1], permutations=999)
dune.meandist
summary(dune.meandist)
plot(dune.meandist)

#PERMDISP2 - multivariate homogeneity of group dispersions (variances)
dune.disp <- betadisper(dune.dist,  group=dune$MANAGEMENT)
permutest(dune.disp)
permutest(dune.disp, pairwise = TRUE)
boxplot(dune.disp)
plot(dune.disp)
anova(dune.disp)
TukeyHSD(dune.disp)
## ----end


## ---- bioenv data
data.dist <- vegdist(wisconsin(data[,-1]^0.25),"bray")
Xmat = model.matrix(~ -1 + pH + Slope + Altitude + Substrate, data=enviro[,-1])
data.bioenv <- bioenv(data.dist,
                      decostand(Xmat,"standardize"))
data.bioenv
## ----end

## ---- adonis spider
## When covariates are continuous, adonis uses
## Distance-based redundancy analysis (dbRDA)
spider.dbrda <- dbrda(spider.dist~ soil.dry + bare.sand + fallen.leaves + moss,
                        data=spider$x)
#ordiplot(spider.dbrda)
autoplot(spider.dbrda)
autoplot(spider.dbrda,geom='text') + theme_bw()
spider.dbrda
eigenvals(spider.dbrda)
spider.eig <- eigenvals(spider.dbrda)
spider.eig/sum(spider.eig)

anova(spider.dbrda, by = "axis")
anova(spider.dbrda, by = "term")
anova(spider.dbrda, by = "margin")
## The above, really just call permutest
permutest(spider.dbrda)
permutest(spider.dbrda, by = "onedf")
permutest(spider.dbrda, by = "terms")


spider.adonis <- adonis2(spider.dist ~ soil.dry + bare.sand + fallen.leaves + moss,
                        data=spider$x)
spider.adonis
## The above is sequentially
## Note, if we change the order (mainly important for interactions)
spider.adonis <- adonis(spider.dist ~ fallen.leaves + soil.dry + bare.sand + moss,
                        data=spider$x)
spider.adonis
## Adonis2 can do either sequential (as next)
spider.adonis <- adonis2(spider.dist ~ soil.dry + bare.sand + fallen.leaves + moss,
                        data=spider$x)
spider.adonis
## or marginal (next)
spider.adonis <- adonis2(spider.dist ~ soil.dry + bare.sand + fallen.leaves + moss,
                        data=spider$x, by = 'margin')
spider.adonis
## or overall
spider.adonis <- adonis2(spider.dist ~ soil.dry + bare.sand + fallen.leaves + moss,
                        data=spider$x, by = NULL)
spider.adonis

spider.disp <- betadisper(spider.dist, group = spider$x$fallen.leaves)
spider.disp <- betadisper(spider.dist, group = cut(spider$x$fallen.leaves, breaks = 3))
spider.disp

## boxplot(spider.disp)
plot(spider.disp)
anova(spider.disp)
## TukeyHSD(spider.disp)
## permutest(spider.disp, pairwise = TRUE)
## ----end

## ---- adonis data
data.adonis <- adonis(data.dist ~ pH + Slope + Altitude + 
                          Substrate, data=enviro)
## data.adonis <- adonis(data.dist ~ Altitude, data=enviro)
data.adonis

## ----end

## ---- varveg
vareveg <- read.csv('../public/data/vareveg.csv')
head(vareveg)
vareenv <- read.csv('../public/data/vareenv.csv')
head(vareenv)


vareveg.dist <- vegdist(wisconsin(vareveg[,-1]^0.25),'bray')
#vareveg.dist <- vegdist(vareveg.std, "bray")
#environmental variables 
vareenv.std <- decostand(vareenv[,-1], "standardize")
vareenv.dist <- vegdist(vareenv.std, "euc")

#bioenv(vareveg.std, vareenv.std)

bioenv(vareveg.dist, vareenv.std)
#adonis(vareveg.std~Ca+Fe+Mn+Baresoil, data=vareenv)
adonis(vareveg.dist~Ca+Fe+Mn+Baresoil, data=vareenv.std)
adonis2(vareveg.dist~Ca+Fe+Mn+Baresoil, data=vareenv.std)

#meandist(vareveg.dist)
#vareveg.simper = simper(vareveg[,-1], vareenv.std$Ca, permutations=100)
#vareveg.simper
#summary(vareveg.simper)

## ----end


library(codyn)
X <- data %>%
    gather(key=Species, value=value,-Sites)
community_diversity(X, abundance.var = "value",
                    replicate.var = 'Sites')
multivariate_change(X,
                    species.var='Species',
                    abundance.var = "value",
                    time.var = 'Sites')
