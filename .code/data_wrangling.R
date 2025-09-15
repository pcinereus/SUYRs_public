## -----------------------------------------------------------------------------
#| label: setup
#| include: false
#| cache: false

knitr::opts_chunk$set(cache.lazy = FALSE,
                      tidy = "styler")
options(tinytex.engine = "xelatex")


## -----------------------------------------------------------------------------
#| label: libraries
#| results: markup
#| eval: true
#| echo: true
#| cache: false
#| message: false
#| warning: false
library(tidyverse) #for data wrangling


## -----------------------------------------------------------------------------
#| label: getData
#| results: markup
#| eval: true
#| echo: true
#| cache: false
load(file='../data/manipulationDatasets.RData')


## -----------------------------------------------------------------------------
#| label: getData1
#| results: markup
#| eval: true
#| echo: true
#| cache: false
dat.1 |> head()


## -----------------------------------------------------------------------------
#| label: exploringdata
#| results: markup
#| eval: true
#| echo: true
#| cache: false
## replace this with code to explore imported data
## From now on I will not provide many code chunks.
## Instead you are incouraged to create your own chunks
## with discussed code
## In Rstudio, you can create a chunk with Cntr-Alt-I


## -----------------------------------------------------------------------------
#| label: piping
#| results: markup
#| eval: true
#| echo: true
#| cache: false
head(dat.1)
#OR
dat.1 |> head()

