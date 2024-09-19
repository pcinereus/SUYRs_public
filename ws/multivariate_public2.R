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

