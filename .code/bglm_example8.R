## -----------------------------------------------------------------------------
#| label: setup
#| include: false
#| cache: false
knitr::opts_chunk$set(cache.lazy = FALSE,
                      tidy = "styler")
options(tinytex.engine = "xelatex")


## -----------------------------------------------------------------------------
#| label: libraries
#| output: false
#| eval: true
#| warning: false
#| message: false
#| cache: false
library(tidyverse)
library(brms)
library(dagitty)
library(ggdag)
library(patchwork)
source('helperFunctions.R')


## -----------------------------------------------------------------------------
#| label: readData
king <- read_csv("../data/king.csv", trim_ws = TRUE)


## -----------------------------------------------------------------------------
#| label: examinData
glimpse(king)


## -----------------------------------------------------------------------------
## Explore the first 6 rows of the data
head(king)


## -----------------------------------------------------------------------------
str(king)


## -----------------------------------------------------------------------------
king |> datawizard::data_codebook()


## -----------------------------------------------------------------------------
#| label: processData
#| eval: false
# king <- king |>
#     mutate(
#       fLight = factor(Light),
#       fDiuron = factor(Diuron),
#       Block = factor(Block)
#     )


## -----------------------------------------------------------------------------
#| label: eda_mod0a
#| results: hide
#| eval: true
#| echo: true
#| cache: true
#| fig-width: 13
#| fig-height: 4
g1 <-
    king |>
    ggplot(aes(y = PI, x = Diuron, colour = factor(Light))) +
    geom_point() +
    geom_smooth(method = "lm") +
    theme_classic()
g2 <-
  king |>
  ggplot(aes(y = celld, x = PI, colour = factor(Light))) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()
g3 <-
  king |>
  ggplot(aes(y = celld, x = Diuron, colour = factor(Light))) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic()
g1 + g2 + g3 + plot_layout(guides = "collect")


## -----------------------------------------------------------------------------
#| label: brms_mod0a
#| eval: true
#| echo: true
#| cache: true
form <- bf(PI ~ scale(Diuron) * factor(Light))
## summarise data to help inform priors
king |>
    group_by(Light) |>
    summarise(PI_mu = median(PI), PI_sd = sd(PI), PI_mad = mad(PI))

## define priors
priors <- prior(normal(1.4, 0.3), class = "Intercept") +
  prior(normal(0, 0.5), class = "b") +
  prior(student_t(3, 0, 0.3), class = "sigma")
## fit the mode
king_brms0a <- brm(
  form,
  prior = priors,
  data = king,
  refresh = 0,
  backend = "cmdstanr"
)


## -----------------------------------------------------------------------------
#| label: brms_mod0a_sum
#| results: markup
#| eval: true
#| echo: true
#| cache: false
summary(king_brms0a)


## -----------------------------------------------------------------------------
#| label: brms_mod0b
#| eval: true
#| echo: true
#| cache: true
form <- bf(celld ~ scale(Diuron) + scale(PI) + Light)
## summarise data to help inform priors
king |>
    group_by(Light) |>
    summarise(celld_mu = median(celld), celld_sd = sd(celld), celld_mad = mad(celld))

## define priors
priors <- prior(normal(150, 50), class = "Intercept") +
  prior(normal(0, 20), class = "b") +
  prior(student_t(3, 0, 50), class = "sigma")
## fit the mode
king_brms0b <- brm(
  form,
  prior = priors,
  data = king |> mutate(Light = factor(Light)),
  refresh = 0,
  backend = "cmdstanr"
)


## -----------------------------------------------------------------------------
#| label: brms_mod0b_sum
#| results: markup
#| eval: true
#| echo: true
#| cache: false
summary(king_brms0b)


## -----------------------------------------------------------------------------
#| label: dag
#| eval: true
#| cache: false
king_dag <- dagify(
    celld ~ PI + Diuron,
    PI ~ Light + Diuron,
    exposure =  "Light",
    outcome = "celld"
)
ggdag(king_dag, text = TRUE, text_size = 2.5) +
  theme_dag_blank()


## -----------------------------------------------------------------------------
#| label: implied conditional independencies
#| eval: true
#| cache: false
dagitty::impliedConditionalIndependencies(king_dag)


## -----------------------------------------------------------------------------
#| label: local tests
#| eval: true
#| cache: false
tests <- localTests(x = king_dag, data = king)
tests


## -----------------------------------------------------------------------------
#| label: ggdag paths 1
#| eval: true
#| cache: false
#| fig-width: 7
#| fig-height: 5
ggdag_paths(king_dag, from = "Light", to = "PI", text_col = "black", shadow = TRUE) +
    theme_dag_blank(panel.border = element_rect(fill = NA))


## -----------------------------------------------------------------------------
#| label: adjustment sets 1a
#| eval: true
#| cache: false
adjustmentSets(king_dag,
  exposure = "Light",
  outcome = "PI",
  type = "minimal",
  effect = "total")


## -----------------------------------------------------------------------------
#| label: adjustment sets 1b
#| eval: true
#| cache: false
adjustmentSets(king_dag,
  exposure = "Light",
  outcome = "PI",
  type = "canonical",
  effect = "total")


## -----------------------------------------------------------------------------
#| label: adjustment sets 1c
#| eval: true
#| cache: false
adjustmentSets(king_dag,
  exposure = "Light",
  outcome = "PI",
  type = "minimal",
  effect = "direct")


## -----------------------------------------------------------------------------
#| label: ggdag paths 2
#| eval: true
#| cache: false
#| fig-width: 7
#| fig-height: 5
ggdag_paths(king_dag, from = "Light", to = "celld", text_col = "black", shadow = TRUE) +
    theme_dag_blank(panel.border = element_rect(fill = NA))


## -----------------------------------------------------------------------------
#| label: adjustment sets 2a
#| eval: true
#| cache: false
adjustmentSets(king_dag,
  exposure = "Light",
  outcome = "celld",
  type = "minimal",
  effect = "total")


## -----------------------------------------------------------------------------
#| label: adjustment sets 2b
#| eval: true
#| cache: false
adjustmentSets(king_dag,
  exposure = "Light",
  outcome = "celld",
  type = "canonical",
  effect = "total")


## -----------------------------------------------------------------------------
#| label: adjustment sets 2c
#| eval: true
#| cache: false
adjustmentSets(king_dag,
  exposure = "Light",
  outcome = "celld",
  type = "minimal",
  effect = "direct")


## -----------------------------------------------------------------------------
#| label: ggdag paths 4
#| eval: true
#| cache: false
#| fig-width: 7
#| fig-height: 5
ggdag_paths(king_dag, from = "PI", to = "celld", text_col = "black", shadow = TRUE) +
    theme_dag_blank(panel.border = element_rect(fill = NA))


## -----------------------------------------------------------------------------
#| label: adjustment sets 4a
#| eval: true
#| cache: false
adjustmentSets(king_dag,
  exposure = "PI",
  outcome = "celld",
  type = "minimal",
  effect = "total")


## -----------------------------------------------------------------------------
#| label: adjustment sets 4b
#| eval: true
#| cache: false
adjustmentSets(king_dag,
  exposure = "PI",
  outcome = "celld",
  type = "canonical",
  effect = "total")


## -----------------------------------------------------------------------------
#| label: adjustment sets 4c
#| eval: true
#| cache: false
adjustmentSets(king_dag,
  exposure = "PI",
  outcome = "celld",
  type = "minimal",
  effect = "direct")


## -----------------------------------------------------------------------------
#| label: ggdag paths 3
#| eval: true
#| cache: false
#| fig-width: 12
#| fig-height: 5
ggdag_paths(king_dag, from = "Diuron", to = "celld", text_col = "black", shadow = TRUE) +
    theme_dag_blank(panel.border = element_rect(fill = NA))


## -----------------------------------------------------------------------------
#| label: adjustment sets 3a
#| eval: true
#| cache: false
adjustmentSets(king_dag,
  exposure = "Diuron",
  outcome = "celld",
  type = "minimal",
  effect = "total")


## -----------------------------------------------------------------------------
#| label: adjustment sets 3b
#| eval: true
#| cache: false
adjustmentSets(king_dag,
  exposure = "Diuron",
  outcome = "celld",
  type = "canonical",
  effect = "total")


## -----------------------------------------------------------------------------
#| label: adjustment sets 3c
#| eval: true
#| cache: false
adjustmentSets(king_dag,
  exposure = "Diuron",
  outcome = "celld",
  type = "minimal",
  effect = "direct")


## -----------------------------------------------------------------------------
#| label: ggdag paths 5
#| eval: true
#| cache: false
#| fig-width: 12
#| fig-height: 5
ggdag_paths(king_dag, from = "Diuron", to = "PI", text_col = "black", shadow = TRUE) +
    theme_dag_blank(panel.border = element_rect(fill = NA))


## -----------------------------------------------------------------------------
#| label: adjustment sets 5a
#| eval: true
#| cache: false
adjustmentSets(king_dag,
  exposure = "Diuron",
  outcome = "PI",
  type = "minimal",
  effect = "total")


## -----------------------------------------------------------------------------
#| label: adjustment sets 5b
#| eval: true
#| cache: false
adjustmentSets(king_dag,
  exposure = "Diuron",
  outcome = "PI",
  type = "canonical",
  effect = "total")


## -----------------------------------------------------------------------------
#| label: adjustment sets 5c
#| eval: true
#| cache: false
adjustmentSets(king_dag,
  exposure = "Diuron",
  outcome = "PI",
  type = "minimal",
  effect = "direct")


## -----------------------------------------------------------------------------
#| label: processData
#| results: hide
#| eval: true
#| echo: true
#| cache: false
king <- king |>
    mutate(
      fLight = factor(Light),
      fDiuron = factor(Diuron),
      Block = factor(Block)
    )


## -----------------------------------------------------------------------------
#| label: eda_mod1a
#| results: hide
#| eval: true
#| echo: true
#| cache: true
king |>
    ggplot(aes(y = PI, x = Diuron, fill = fLight)) +
    geom_boxplot() +
    theme_classic()


## -----------------------------------------------------------------------------
#| label: brms_mod1a
#| eval: true
#| echo: true
#| cache: true
form <- bf(PI ~ fLight * scale(Diuron))
## summarise data to help inform priors
king |>
    group_by(fLight) |>
    summarise(PI_mu = median(PI), PI_sd = sd(PI), PI_mad = mad(PI))

## define priors
priors <- prior(normal(1.4, 0.3), class = "Intercept") +
  prior(normal(0, 0.5), class = "b") +
  prior(student_t(3, 0, 0.3), class = "sigma")
## fit the mode
king_brms1a <- brm(
  form,
  prior = priors,
  data = king,
  refresh = 0,
  backend = "cmdstanr"
)


## -----------------------------------------------------------------------------
#| label: brms_mod1a_sum
#| results: markup
#| eval: true
#| echo: true
#| cache: false
summary(king_brms1a)


## -----------------------------------------------------------------------------
#| label: eda_mod2a
#| results: hide
#| eval: true
#| echo: true
#| cache: true
king |>
    ggplot(aes(y = celld, x = Diuron, fill = fLight)) +
    geom_boxplot() +
    theme_classic()


## -----------------------------------------------------------------------------
#| label: brms_mod2a
#| eval: true
#| echo: true
#| cache: true
form <- bf(celld ~ fLight * scale(Diuron))
## summarise data to help inform priors
king |>
    group_by(fLight) |>
    summarise(celld_mu = median(celld), celld_sd = sd(celld), celld_mad = mad(celld))

## define priors
priors <- prior(normal(150, 50), class = "Intercept") +
  prior(normal(0, 30), class = "b") +
  prior(student_t(3, 0, 50), class = "sigma")
## fit the mode
king_brms2a <- brm(
  form,
  prior = priors,
  data = king,
  refresh = 0,
  backend = "cmdstanr"
)


## -----------------------------------------------------------------------------
#| label: brms_mod2a_sum
#| results: markup
#| eval: true
#| echo: true
#| cache: false
summary(king_brms2a)


## -----------------------------------------------------------------------------
#| label: eda_mod3a
#| results: hide
#| eval: true
#| echo: true
#| cache: true
king |>
    ggplot(aes(y = celld, x = PI, fill = fLight)) +
    geom_boxplot() +
    theme_classic()


## -----------------------------------------------------------------------------
#| label: brms_mod3a
#| eval: true
#| echo: true
#| cache: true
form <- bf(celld ~ scale(PI) * fLight * scale(Diuron))
## summarise data to help inform priors
king |>
    group_by(fLight) |>
    summarise(celld_mu = median(celld), celld_sd = sd(celld), celld_mad = mad(celld))

## define priors
priors <- prior(normal(150, 50), class = "Intercept") +
  prior(normal(0, 30), class = "b") +
  prior(student_t(3, 0, 50), class = "sigma")
## fit the mode
king_brms3a <- brm(
  form,
  prior = priors,
  data = king,
  refresh = 0,
  backend = "cmdstanr"
)


## -----------------------------------------------------------------------------
#| label: brms_mod3a_sum
#| results: markup
#| eval: true
#| echo: true
#| cache: false
summary(king_brms3a)


## -----------------------------------------------------------------------------
#| label: eda_mod4a
#| results: hide
#| eval: true
#| echo: true
#| cache: true
king |>
    ggplot(aes(y = celld, x = Diuron, fill = fLight)) +
    geom_boxplot() +
    theme_classic()


## -----------------------------------------------------------------------------
#| label: brms_mod4a
#| eval: true
#| echo: true
#| cache: true
form <- bf(celld ~ fLight * scale(Diuron))
## summarise data to help inform priors
king |>
    group_by(fLight) |>
    summarise(celld_mu = median(celld), celld_sd = sd(celld), celld_mad = mad(celld))

## define priors
priors <- prior(normal(150, 50), class = "Intercept") +
  prior(normal(0, 30), class = "b") +
  prior(student_t(3, 0, 50), class = "sigma")
## fit the mode
king_brms4a <- brm(
  form,
  prior = priors,
  data = king,
  refresh = 0,
  backend = "cmdstanr"
)


## -----------------------------------------------------------------------------
#| label: brms_mod4a_sum
#| results: markup
#| eval: true
#| echo: true
#| cache: false
summary(king_brms4a)

