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

