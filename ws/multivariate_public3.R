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


