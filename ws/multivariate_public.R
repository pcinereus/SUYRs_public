#############################################################
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
##Unconstrained  ---> ordination
##               ---> envfit (overlay enviromental data) (permutation test)
##               ---> lm/glm etc (response or predictor)
#############################################################
##     Dissimilarity
##            --> MDS      ---> ordination
##            --> bioenv   ---> variable importance       (perm test)
##            --> adonis   ---> anova                     (perm test)
##            --> simper   ---> similarity percentages
##            --- betadisp ---> homogeneity of dispersion (perm test)
#############################################################
##     Model based
##            ---> manyglm ---> anova
#############################################################

## ---- libraries
library(tidyverse)
library(vegan)        # for multivariate analyses
library(GGally)       # for scatterplot matrices
library(corrplot)     # for association plots
library(car)
library(mvabund)      # for model-based multivariate analyses
library(scales)
library(ggvegan)      ## ggplot support for vegan
library(ggrepel)      ## for geom_text_repel
## ----end


## PCA ------------------------------------------------------


## ---- readSpider
spider.abund <- read_csv(file = "../data/spider.abund.csv", trim_ws = TRUE)
spider.env <- read_csv(file = "../data/spider.env.csv", trim_ws = TRUE)
glimpse(spider.abund)
glimpse(spider.env)
## ----end

## Exploratory data analysis

## 3D plot - do not put in qmd
library(rgl)
spider.abund.df <- as.data.frame(spider.abund)
plot3d( 
  x=spider.abund.df[,'Pardpull']^0.25,
  y=spider.abund.df[,'Pardnigr']^0.25,
  z=spider.abund.df[,'Trocterr']^0.25, 
  type = 's', 
  radius = 0.1)


spider.env.df <- as.data.frame(spider.env)
plot3d( 
  x=spider.env.df[,'fallen.leaves'],
  y=spider.env.df[,'soil.dry'],
  z=spider.env.df[,'moss'], 
  type = 's', 
  radius = 0.1)

## ---- EDA spider
spider.abund |>
  cor() |>
  corrplot(type = 'upper',
    diag = FALSE)
spider.abund |>
  cor() |>
  corrplot(type = 'upper',
    order = 'FPC',
    diag = FALSE)
## ----end

## From here on, it is up to you to create your chunks etc

### Univariate measures

## diversity and richness
spider.abund.comp <- spider.abund |>
  mutate(Site = factor(1:n())) |>
  pivot_longer(cols = -Site, names_to = "Species", values_to = "Abundance") |>
  group_by(Site) |>
  mutate(diversity = diversity(Abundance, index = "shannon")) |>
  mutate(specnumber = specnumber(Abundance)) |>
  ungroup()
spider.abund.comp |>
  ungroup() |>
  group_by(Site) |>
  summarise(
    diversity = unique(diversity),
    specnumber = unique(specnumber)
  ) |>
  ggplot(aes(x = Site, y = diversity * 3)) +
  geom_blank() +
  geom_line(aes(x = as.numeric(Site))) +
  geom_line(aes(x = as.numeric(Site), y = specnumber), color = "red") +
  scale_y_continuous(
    name = "Diversity",
    label = function(x) x / 3, ,
    sec.axis = sec_axis(~ . * 1,
      name = "Species number"
    )
  ) +
  theme_classic() +
  theme(axis.title.y.right = element_text(color = "red"))

## rarefaction
spider.abund.comp |>
  ungroup() |>
  group_by(Site) |>
  summarise(
    diversity = unique(diversity),
    specnumber = unique(specnumber),
    rarefy = rarefy(Abundance, sample = min(rowSums(spider.abund)))    
  ) |>
  ggplot(aes(x = Site, y = rarefy)) +
  geom_blank() +
  geom_line(aes(x = as.numeric(Site))) +
  geom_line(aes(x = as.numeric(Site), y = specnumber), color = "red") +
  scale_y_continuous(
    name = "Rarefy",
    ## label = function(x) x / 3, ,
    sec.axis = sec_axis(~ . * 1,
      name = "Species number"
    )
  ) +
  theme_classic() +
  theme(axis.title.y.right = element_text(color = "red"))



## More EDA

spider.abund |>
  ggpairs(lower = list(continuous = "smooth"),
    diag = list(continuous = "density"),
    axisLabels = "show")

spider.abund^0.25 |>
  ggpairs(lower = list(continuous = "smooth"),
    diag = list(continuous = "density"),
    axisLabels = "show")

spider.std <- (spider.abund^0.25) |> 
  wisconsin()
spider.std

## PCA
spider.rda <- rda(spider.std, scale=TRUE)
summary(spider.rda, display=NULL)
screeplot(spider.rda)
abline(a=1,b=0)
spider.rda <- rda(spider.std, scale=FALSE)
summary(spider.rda, display=NULL)


autoplot(spider.rda)
autoplot(spider.rda) + theme_bw()
autoplot(spider.rda, geom='text') + theme_bw()


spider.rda.scores <- spider.rda |> 
  fortify()

ggplot(data = NULL, aes(y=PC2, x=PC1)) +
  geom_hline(yintercept=0, linetype='dotted') +
  geom_vline(xintercept=0, linetype='dotted') +
  geom_point(data=spider.rda.scores |> filter(score=='sites')) +
  geom_text(data=spider.rda.scores |> filter(score=='sites'),
    aes(label=label), hjust=-0.2) +
  geom_segment(data=spider.rda.scores |> filter(score=='species'),
    aes(y=0, x=0, yend=PC2, xend=PC1),
    arrow=arrow(length=unit(0.3,'lines')), color='red') +
  geom_text(data=spider.rda.scores |> filter(score=='species'),
    aes(y=PC2*1.1, x=PC1*1.1, label=label), color='red') 
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

## The following is a demonstration to illustrate the use of sprintf
eig <- eigenvals(spider.rda)
sprintf('(%0.1f%% explained var.)', 100 * eig[2]/sum(eig))
paste(names(eig[2]), sprintf('(%0.1f%% explained var.)', 100 * eig[2]/sum(eig)))


g <- g +
  scale_y_continuous(paste(names(eig[2]), sprintf('(%0.1f%% explained var.)',
    100 * eig[2]/sum(eig))))+
  scale_x_continuous(paste(names(eig[1]), sprintf('(%0.1f%% explained var.)',
    100 * eig[1]/sum(eig))))

circle.prob <- 0.68
r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(spider.rda$CA$u[,1:2]^2))^(1/4)
theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
circle <- data.frame(PC1 = r * cos(theta), PC2 = r * sin(theta))
g <- g + geom_path(data = circle, aes(y=PC2,x=PC1), color = muted('white'), size = 1/2, alpha = 1/3)
g


## de-emphasize those that are inside the circle (to prevent cluttering)
spider.rda.scores.sites <- spider.rda.scores |>
  filter(score == "sites") |>
  mutate(Flag = factor(ifelse(sqrt(PC1^2 + PC2^2) > r, 1, 0)))
spider.rda.scores.species <- spider.rda.scores |>
  filter(score == "species") |>
  mutate(Flag = factor(ifelse(sqrt(PC1^2 + PC2^2) > r, 1, 0)))

g <-
  ggplot(data = NULL, aes(y=PC2, x=PC1)) +
  geom_hline(yintercept=0, linetype='dotted') +
  geom_vline(xintercept=0, linetype='dotted') +
  geom_point(data=spider.rda.scores.sites, aes(alpha = Flag), show.legend = FALSE) +
  geom_text(data=spider.rda.scores.sites,
    aes(label=label, alpha = Flag), hjust=-0.2, show.legend = FALSE) +
  geom_segment(data=spider.rda.scores.species,
    aes(y=0, x=0, yend=PC2, xend=PC1, alpha = Flag),
    arrow=arrow(length=unit(0.3,'lines')), color='red', show.legend = FALSE) +
  geom_text_repel(data=spider.rda.scores.species,
    aes(y=PC2*1.1, x=PC1*1.1, label=label, alpha = Flag),
    color='red', show.legend = FALSE) +
  geom_path(data = circle, aes(y=PC2,x=PC1),
    color = muted('white'), size = 1/2, alpha = 1/3) +
  scale_y_continuous(paste(names(eig[2]), sprintf('(%0.1f%% explained var.)',
    100 * eig[2]/sum(eig))))+
  scale_x_continuous(paste(names(eig[1]), sprintf('(%0.1f%% explained var.)',
    100 * eig[1]/sum(eig)))) +
  theme_bw()
g

  
# Envfit
spider.env |>
  cor() |> 
  corrplot(type = 'upper',
    order = 'FPC',
    diag = FALSE)

spider.env |> 
  ggpairs(lower = list(continuous = "smooth"),
    diag = list(continuous = "density"),
    axisLabels = "show")

spider.envfit <- envfit(spider.rda, env = spider.env)
spider.envfit

spider.env.scores <- spider.envfit |> fortify()
g <-
  g + 
  geom_segment(data=spider.env.scores,
    aes(y=0, x=0, yend=PC2, xend=PC1),
    arrow=arrow(length=unit(0.3,'lines')), color='blue') +
  geom_text(data=spider.env.scores,
    aes(y=PC2*1.1, x=PC1*1.1, label=label), color='blue')
l


pc1 <- spider.rda.scores |> filter(score=='sites') |> pull(PC1)
pc2 <- spider.rda.scores |> filter(score=='sites') |> pull(PC2)

lm(1:nrow(spider.env) ~ soil.dry + bare.sand + fallen.leaves +
     moss + herb.layer + reflection, data=spider.env) |>
  vif()
lm(1:nrow(spider.env) ~ herb.layer + fallen.leaves + bare.sand + moss, data=spider.env) |>
  vif()

lm(pc1 ~ herb.layer + fallen.leaves + bare.sand + moss, data=spider.env) %>%
  summary()
lm(pc2 ~ herb.layer + fallen.leaves + bare.sand + moss, data=spider.env) %>%
  summary()

## END PCA ------------------------------------------------------

## RDA ----------------------------------------------------------
spider.rda <- rda(spider.std ~ 
                    scale(herb.layer)+
                    scale(fallen.leaves) +
                    scale(bare.sand) +
                    scale(moss),
  data=spider.env, scale=FALSE)
summary(spider.rda, display=NULL)

vif.cca(spider.rda)


goodness(spider.rda)
goodness(spider.rda, display = "sites")

inertcomp(spider.rda)
inertcomp(spider.rda, proportional = TRUE)


anova(spider.rda)
anova(spider.rda, by='axis')
anova(spider.rda, by='margin')

coef(spider.rda)

RsquareAdj(spider.rda)

screeplot(spider.rda)

autoplot(spider.rda, geom='text')
## END RDA ------------------------------------------------------

## CA ------------------------------------------------------
data <- read_csv('../public/data/data.csv', trim_ws=TRUE)
head(data) 
enviro <- read_csv('../public/data/enviro.csv', trim_ws=TRUE)
head(enviro)
enviro <- enviro |> mutate(Substrate=factor(Substrate))

data.std <- data |>
    dplyr::select(-Sites) |> 
    decostand(method="total",MARGIN=2)
data.std

data.std |>
  cor() |>
    corrplot(diag=FALSE)

data.std |> cor() |>
    corrplot(diag=FALSE, order='FPC')

data.ca <- cca(data.std, scale=FALSE)
summary(data.ca, display=NULL)

screeplot(data.ca)
sum(eigenvals(data.ca))/length(eigenvals(data.ca))
eigenvals(data.ca)/sum(eigenvals(data.ca))

autoplot(data.ca)
autoplot(data.ca) + theme_bw()
autoplot(data.ca, geom='text') + theme_bw()

data.ca.scores <- data.ca |> 
    fortify()
data.ca.scores |> head()


g <-
    ggplot(data = NULL, aes(y=CA2, x=CA1)) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_vline(xintercept=0, linetype='dotted') +
    geom_point(data=data.ca.scores |> filter(score=='sites')) +
    geom_text(data=data.ca.scores |> filter(score=='sites'),
              aes(label=label), hjust=-0.2) +
    geom_segment(data=data.ca.scores |> filter(score=='species'),
                 aes(y=0, x=0, yend=CA2, xend=CA1),
                 arrow=arrow(length=unit(0.3,'lines')), color='red') +
    ## geom_text(data=data.rda.scores |> filter(score=='species'),
    ##           aes(y=PC2*1.1, x=PC1*1.1, label=label), color='red') +
    geom_text_repel(data=data.ca.scores |> filter(score=='species'),
              aes(y=CA2*1.1, x=CA1*1.1, label=label), color='red') +
    theme_bw()
g



Xmat <- model.matrix(~ -1+pH+Slope+Altitude+Substrate, data = enviro)
data.envfit <- envfit(data.ca, env=Xmat)
data.envfit
autoplot(data.envfit)

data.env.scores <- data.envfit |> fortify()
g <- g + 
  geom_segment(data=data.env.scores,
    aes(y=0, x=0, yend=CA2, xend=CA1),
    arrow=arrow(length=unit(0.3,'lines')), color='blue') +
  geom_text(data=data.env.scores,
    aes(y=CA2*1.1, x=CA1*1.1, label=label), color='blue')
g


data.ca.scores <- data.ca %>% fortify()
CA1 <- data.ca.scores %>% filter(Score =='sites') %>% pull(CA1)
CA2 <- data.ca.scores %>% filter(Score =='sites') %>% pull(CA2)
summary(lm(CA1 ~ pH+Slope+Altitude+Substrate, data=enviro))
summary(lm(CA2 ~ pH+Slope+Altitude+Substrate, data=enviro))

## END CA --------------------------------------------------

## CCA ----------------------------------------------------------
data.cca <- cca(data.std~pH + Altitude + Substrate + Slope, data=enviro, scale=FALSE)

summary(data.cca, display=NULL)
anova(data.cca)

autoplot(data.cca)


vif.cca(data.cca)
#overall test
anova(data.cca)
anova(data.cca, by='axis')
anova(data.cca, by='margin')

coef(data.cca)

RsquareAdj(data.cca)

screeplot(data.cca)
## int <- data.cca$tot.chi/length(data.cca$CA$eig)
## abline(h=int)

## END CCA ------------------------------------------------------

## PCoA ----------------------------------------------------------
data.std <-
    data |> dplyr::select(-Sites) |>
    decostand(method="total",MARGIN=2)
data.std

data.dist <- vegdist(data.std, method='bray')
data.capscale <- capscale(data.dist~1, data=enviro)

summary(data.capscale, display=NULL)
autoplot(data.capscale, geom='text')

#Distance based redundancy analysis
data.capscale <- capscale(data.dist~scale(pH) + scale(Altitude) + Substrate + scale(Slope), data=enviro)
summary(data.capscale, display=NULL)
plot(data.capscale)
autoplot(data.capscale, geom='text')

summary(data.capscale, display=NULL)
anova(data.capscale)

anova(data.capscale, by='margin')
screeplot(data.capscale)
sum(eigenvals(data.capscale))/length(eigenvals(data.capscale))
eigenvals(data.capscale)/sum(eigenvals(data.capscale))

## Conditioning on
data.capscale <- capscale(data.dist~pH + Condition(Altitude) + Substrate + Slope, data=enviro)
summary(data.capscale, display=NULL)
plot(data.capscale)
autoplot(data.capscale, geom='text')
## END PCoA ------------------------------------------------------


## MDS ------------------------------------------------------
macnally <- read.csv('../data/macnally_full.csv',strip.white=TRUE)
head(macnally)
macnally[1:5,1:5]

macnally <- macnally |>
  mutate(HABITAT = factor(HABITAT, levels = c(
    "Mixed", "Gipps.Manna",
    "Montane Forest", "Foothills Woodland", "Box-Ironbark", "River Red Gum"
  )))

macnally.mds <- metaMDS(macnally[,-1], k=2,  plot=TRUE)
macnally.mds

macnally.std <- wisconsin(macnally[,c(-1)]^0.25)

macnally.dist <- vegdist(macnally.std,"bray")

macnally.mds <- metaMDS(macnally.std, k=2, plot=TRUE)
macnally.mds <- metaMDS(macnally.dist, k=2, plot=TRUE)
macnally.mds <- metaMDS(macnally[,-1], k=2)

macnally.mds$stress

stressplot(macnally.mds)


macnally.mds.scores <- macnally.mds |> 
  fortify() |> 
  full_join(macnally |>
             rownames_to_column(var='label'),
    by =  'label')

g <-
    ggplot(data = NULL, aes(y=NMDS2, x=NMDS1)) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_vline(xintercept=0, linetype='dotted') +
    geom_point(data=macnally.mds.scores |> filter(score=='sites'),
               aes(color=HABITAT)) +
    geom_text(data=macnally.mds.scores |> filter(score=='sites'),
              aes(label=label, color=HABITAT), hjust=-0.2) +
    geom_segment(data=macnally.mds.scores |> filter(score=='species'),
                 aes(y=0, x=0, yend=NMDS2, xend=NMDS1),
                 arrow=arrow(length=unit(0.3,'lines')), color='red',
      alpha =  0.2) +
    geom_text(data=macnally.mds.scores |> filter(score=='species'),
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

g + ggforce::geom_mark_ellipse(data=macnally.mds.scores |> filter(score=='sites'),
                      aes(y=NMDS2, x=NMDS1, fill=HABITAT), expand=0) 
## For the following you will be asked to install concaveman
## g + ggforce::geom_mark_hull(data=macnally.mds.scores |> filter(score=='sites'),
##                       aes(y=NMDS2, x=NMDS1, fill=HABITAT), expand=0) 
## g + ggforce::geom_mark_hull(data=macnally.mds.scores |> filter(score=='sites'),
##                       aes(y=NMDS2, x=NMDS1, fill=HABITAT), expand=0, concavity = 20) 


macnally.hull = macnally.mds.scores %>%
    filter(score=='sites') %>%
  group_by(HABITAT) %>%
  slice(chull(NMDS1, NMDS2))

g + geom_polygon(data=macnally.hull, aes(y=NMDS2,x=NMDS1, fill=HABITAT),
  alpha=0.6, colour = "black")+
  theme_classic()
g

## spider plots
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
g1 + geom_segment(data = macnally.mds.scores,
  aes(x = NMDS1_c, xend = NMDS1, y = NMDS2_c, yend = NMDS2, colour = HABITAT)) + geom_polygon(data=macnally.hull, aes(y=NMDS2,x=NMDS1, fill=HABITAT),
  alpha=0.3, colour = "black")+
  theme_classic()



Xmat <- model.matrix(~-1+HABITAT, data=macnally)
colnames(Xmat) <-gsub("HABITAT","",colnames(Xmat))
envfit <- envfit(macnally.mds, env=Xmat)
envfit


macnally.env.scores <- envfit |> fortify()
g <- g + 
    geom_segment(data=macnally.env.scores,
                 aes(y=0, x=0, yend=NMDS2, xend=NMDS1),
                 arrow=arrow(length=unit(0.3,'lines')), color='blue') +
    geom_text(data=macnally.env.scores,
              aes(y=NMDS2*1.1, x=NMDS1*1.1, label=label), color='blue')
g


macnally.dist <- vegdist(macnally[,-1], 'bray')
## bioenv(macnally.dist, macnally$HABITAT)
##                       decostand(macnally$,"standardize"))

adonis2(macnally.dist ~ HABITAT, data=macnally)
mm <-  model.matrix(~-1 + HABITAT, data=macnally)
## adonis2(macnally.dist ~ mm)
head(mm)
colnames(mm) <-gsub("HABITAT","",colnames(mm))
mm <- data.frame(mm)
macnally.adonis<-adonis2(macnally.dist ~ Box.Ironbark + Foothills.Woodland + Gipps.Manna +
                       Montane.Forest + River.Red.Gum, data=mm,
                    perm=9999)
macnally.adonis

library(pairwiseAdonis)
pairwise.adonis(macnally.dist, macnally$HABITAT)



macnally.disp <- betadisper(macnally.dist, macnally$HABITAT)
boxplot(macnally.disp)
plot(macnally.disp)
anova(macnally.disp)
permutest(macnally.disp, pairwise = TRUE)
TukeyHSD(macnally.disp)


macnally.disp <- betadisper(macnally.dist, macnally$HABITAT, type="median",bias.adjust = TRUE)
boxplot(macnally.disp)
plot(macnally.disp)
anova(macnally.disp)
permutest(macnally.disp, pairwise = TRUE)
TukeyHSD(macnally.disp)


macnally.std <- wisconsin(macnally[,c(-1)]^0.25)
simper(macnally.std, macnally$HABITAT)
## END MDS ------------------------------------------------------

## Another analysis -------------------------------------
dune <- read_csv('../data/dune.csv', trim_ws=TRUE)
dune <- dune %>% mutate(MANAGEMENT=factor(MANAGEMENT,  levels=c("NM","BF","HF","SF"))) %>%
  as.data.frame()
#dune <- read.csv('../downloads/data/dune.csv')
dune |> head()

dune.dist <- vegdist(wisconsin(dune[,-1]^0.25), "bray")
dune.mds = metaMDS(dune.dist, k=2)
dune.mds = metaMDS(dune[,-1], k=2)

autoplot(dune.mds, geom=c('text'))

dune.adonis<-adonis2(dune.dist~MANAGEMENT,  data=dune)
dune.adonis

mm <- model.matrix(~ MANAGEMENT, data=dune)
head(mm)
colnames(mm) <-gsub("MANAGEMENT","",colnames(mm))
mm <- data.frame(mm)
dune.adonis<-adonis2(dune.dist~BF+HF+SF, data=mm,
                    perm=9999)
dune.adonis

library(pairwiseAdonis)
pairwise.adonis(dune.dist, dune$MANAGEMENT)

## library(EcolUtils)
## adonis.pair(dune.dist, dune$MANAGEMENT, nper = 10000)

dune.simper=simper(dune[,-1], dune[,1], permutations = 999)
summary(dune.simper)


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
## End Another analysis -------------------------------------


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

## ---- glmmTMB spider
dat.spider.1 <- spider.abund %>%
    as.data.frame() %>%
    mutate(Site = factor(1:n())) %>%
    pivot_longer(cols = -Site,
                 names_to = 'Species',
                 values_to = 'Abund')
library(glmmTMB)
spider.glmmTMB <- glmmTMB(Abund ~ 1 + rr(Species + 0|Site, d = 2),
                          family = nbinom2(),
                          dat = dat.spider.1
                          )
spider.loadings <- spider.glmmTMB$obj$env$report(
                     spider.glmmTMB$fit$parfull)$fact_load[[1]] %>%
                      as.data.frame() %>%
                      mutate(Species = colnames(spider.abund))
fit <-
    ranef(spider.glmmTMB)[[1]]$Site %>%
    mutate(Site = 1:n())

ggplot(fit, aes(y = SpeciesAlopcune, x = SpeciesAlopacce)) +
    geom_text(aes(label = Site)) +
    geom_text(data = spider.loadings, aes(y = V2, x = V1, label = Species), color = 'blue')
## ----end

## ---- gllvm spiders
library(gllvm)
fitx <- gllvm(y = spider$abund, X=spider.env, family = "negative.binomial")
fitx
par(mfrow = c(1,2))
plot(fitx, which = 1:2)
summary(fitx)
coefplot(fitx, mfrow = c(3,2), cex.ylab = 0.8)
crx <- getResidualCor(fitx)
corrplot(crx, diag = FALSE, type = "lower", method = "square", tl.srt = 25)

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
anova(data.mod, cor.type = 'R')
anova(data.mod, cor.type = 'shrink')
## We can also explore the individal univariate tests.
anova(data.mod, p.uni='adjusted')
summary(data.mod, test="LR")

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




## ----end

## mvgam
library(mvgam)

spider_abund <- spider.abund |>
  as_tibble() |>
  mutate(time = 1:n()) |>
  pivot_longer(
    cols = c(everything(), -time),
    names_to = "species", values_to = "abund"
  ) |>
  mutate(series = factor(species))

spider_abund |>
  summarise(
    Median = log(median(abund)),
    Mean = log(mean(abund)),
    MAD = log(mad(abund)),
    SD = log(sd(abund))
  )

get_mvgam_priors(
  abund ~ 1,
  trend_model = AR(),
  data = spider_abund,
  use_lv = TRUE,
  n_lv = 2,
  family = poisson()
)

mod_mvgam <- mvgam(abund ~ 1,
  trend_model = AR(),
  data = spider_abund,
  use_lv = TRUE,
  n_lv = 2,
  family = poisson(),
  return_model_data = TRUE,
)

mod_mvgam |> summary()
plot(mod_mvgam, type = 'residuals', series = 1)
## plot(mod_mvgam, type = 'residuals', series = 2)
##  plot(mod_mvgam, type = 'factors')

## site scores
mcmc <- mvgam:::mcmc_chains(mod_mvgam$model_output, "LV")
betas <- apply(mcmc, 2, median)
beta1 <- betas[seq(1, 56, 2)]
beta2 <- betas[seq(2, 56, 2)]
b <- data.frame(
  beta1 = beta1,
  beta2 = beta2, site = 1:28
)
## species scores
scores <- mvgam:::mcmc_chains(mod_mvgam$model_output, "lv_coefs")
b_scores <- apply(scores, 2, median)
beta1 <- b_scores[seq(1, 24, 2)]
beta2 <- b_scores[seq(2, 24, 2)]
scores<- data.frame(
  beta1 = scales::rescale(beta1, from = range(beta1), to = range(b$beta1)),
  beta2 = scales::rescale(beta2, from = range(beta2), to = range(b$beta2)), species = colnames(spider_abund)
)

ggplot(b, aes(x = beta2, y = beta1)) +
  geom_text(aes(label = site)) +
  geom_text(data = scores, aes(x = beta2, y = beta1, label = species), color = "red")

