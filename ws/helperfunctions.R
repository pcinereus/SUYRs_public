## @title Extract the prior and posterior distributions of the parameters
## @description This function extracts the prior and posterior distributions of the parameters from a brms model
## @param mod A brms model
## @return A ggplot object
## @examples
## mod <- brm(count ~ zAge + zBase * Trt + (1|patient), data = epilepsy, family = poisson())
## SUYR_prior_and_posterior(mod)
## @export
SUYR_prior_and_posterior <- function(mod) {
    dat <- mod$data
    terms <- attr(dat, 'terms')
    response <- all.vars(update(terms, .~ 1))
    predictors <- all.vars(terms)[-1]
    ## rhs <- mod$formula %>% as.formula %>% brms:::str_rhs()

    f <- mod$formula %>% as.formula %>% update(NULL ~.)
    rhs  <-
        deparse1(f) %>%
        str_remove("~") %>%
        ## paste(f[2],f[3],sep='~') %>%
        str_split("\\+") %>%
        unlist() %>%
        str_trim()
    rhs
    ## ## exclude any terms with a "|" or offset
    ## rhs <-
    ##     rhs[-grep("\\|",rhs)]
    wch.rnd <- rhs[grep("\\|", rhs)]
    if (length(wch.rnd)>0) f <- update(f, paste("~ . -",wch.rnd))
    no.offset <- function(x, preserve = NULL) {
      k <- 0
      proc <- function(x) {
        if (length(x) == 1) return(x)
        if (x[[1]] == as.name("offset") && !((k<<-k+1) %in% preserve)) return(x[[1]])
        replace(x, -1, lapply(x[-1], proc))
      }
      update(proc(x), ~ . - offset)
    }
    f <- no.offset(f)

    Xmat <- model.matrix(f, dat)[,-1] %>%
        as.matrix() %>% 
        colMeans() 
    ## if (length(Xmat)==1) Xmat <- Xmat %>% as.matrix()
    ## Xmat <- dat %>%
    ##     select(any_of(rhs)) %>%
    ##     summarise(across(everything(), mean)) %>%
    ##     as.matrix()

    b <- mod %>%
        as_draws_df() %>%
        dplyr::select(starts_with('b_'),
                      -contains('Intercept')) %>%
        as.matrix() %>%
        suppressWarnings()
    
    scal <- as.vector(Xmat %*% t(b))
    
    ## fixed effects
    ## brms_names <- brms:::change_effects(brmsterms(mod$formula),
    ##                                     data = dat,
    ##                                     pars = variables(mod))
    ## brms_names <- sapply(brms_names, function(x) str_remove(x$fnames, "b_"))
    priors <- mod %>% get_variables() %>% str_subset("^prior_.*") %>% str_subset("lprior", negate = TRUE) 
    pars <- priors |>
      str_subset("Intercept|b_") |>
      str_replace_all("^prior_", "") |> 
      str_replace("Intercept", "b_Intercept")
    ## pars <- mod %>% get_variables() %>%
    ##     str_subset(paste0("^b_(Intercept|",paste0(brms_names,collapse="|"),")"))
    
    ## auxillary
    aux <- priors %>% str_subset("prior_(Intercept|b)", negate = TRUE) %>%
        str_remove("^prior_")
    pars <- c(pars, aux)
    ## random effects
    get_variables(mod)
    if (length(wch.rnd)>0) {
        ## ran.pars <- brms:::change_re(mod$ranef, pars = variables(mod))[[1]]$fnames
        ran.pars <- get_variables(mod) |>
          str_subset("^sd_")
        pars <- c(pars, ran.pars)
    }
    ## variables(mod)

    ## Alternative...?
    vars <- variables(mod)
    priors <- vars %>% str_subset("prior") %>% str_subset("lprior", negate = TRUE)
    all.pars <- priors %>% str_remove("prior_")
    fixed.pars <- vars %>% str_subset("^b_")
    other.pars <- all.pars %>% str_subset("^Intercept$|^b_", negate = TRUE)
    other.pars <- vars %>% str_subset(paste0("^", other.pars, collapse = '|'))
    pars <- c(fixed.pars, other.pars)
    
    ## coefs <- prior_summary(mod)$class %>% unique()
    ## coefs.regex <- paste0("^b_", coefs, collapse = "|")

    ## get_priors(mod) |> filter(str_detect(Parameter, "prior_"))

    brms_names <- fixed.pars |> str_subset("Intercept", negate = TRUE) |> str_replace("^b_", "")
    mod.pp <- mod %>%
        as_draws_df() %>%
        dplyr::select(any_of(c(pars, priors))) %>%
        mutate(b_Intercept = b_Intercept + scal) %>%
        pivot_longer(cols=everything(), names_to='key', values_to='value') %>% 
        mutate(Type = ifelse(str_detect(key, 'prior'), 'Prior', 'Posterior'),
               Parameter = ifelse(Type == 'Prior',
                                  str_remove(key, "^prior_"),
                                  str_remove(key, "^b_")
                                  ),
               ## Parameter = ifelse(Type == 'Posterior',
               ##                     str_remove(Parameter, "__.*"),
               ##                    Parameter),
               ## Class = ifelse(Parameter %in% brms_names, 'b', Parameter),
               Class = ifelse(Parameter %in% brms_names, 'b', Parameter),
               Class = ifelse(Type == 'Posterior', str_remove(Class, "__.*"), Class),
               Class = ifelse(Type == 'Posterior' & Parameter %in% str_remove(priors, "prior_b_"), paste0("b_", Parameter), Class)
          ) %>%
        suppressWarnings()
    
    return(
        ggplot(data = NULL, aes(x=Type,  y=value)) +
        stat_pointinterval(data = mod.pp %>% filter(Type == 'Prior')) +
        stat_pointinterval(data = mod.pp %>% filter(Type != 'Prior' & (Class != 'b' | Parameter == 'b')),
                           aes(colour = Parameter), position = position_dodge()) +
        stat_pointinterval(data = mod.pp %>% filter(Type != 'Prior' & (Class == 'b' & Parameter != 'b')),
                           aes(colour = Parameter), position = position_dodge())+
        facet_wrap(~Class,  scales='free')
        )
}


make_brms_dharma_res <- function(brms_model, seed = 10, ...) {
                                        # equivalent to `simulateResiduals(lme4_model, use.u = FALSE)`
                                        # cores are set to 1 just to ensure reproducibility
    options(mc.cores = 1)
    on.exit(options(mc.cores = parallel::detectCores()))
    response <- brms::standata(brms_model)$Y
    ndraws <- nrow(as_draws_df(brms_model))
    manual_preds_brms <- matrix(0, ndraws, nrow(brms_model$data))
    random_terms <- insight::find_random(
                                 brms_model, split_nested = TRUE, flatten = TRUE
                             )
                                        # for this to have a similar output to `glmmTMB`'s default, we need to
                                        #   create new levels in the hierarchical variables, so then we can
                                        #   use `allow_new_levels = TRUE` and `sample_new_levels = "gaussian"` in
                                        #   `brms::posterior_epred`. This is equivalent to
                                        #   `simulateResiduals(lme4_model, use.u = FALSE)`. See details in
                                        #   `lme4:::simulate.merMod` and `glmmTMB:::simulate.glmmTMB`
    new_data <- brms_model$data |>
        dplyr::mutate(across(
                   all_of(random_terms), \(x)paste0("NEW_", x) |> as.factor()
               ))
    set.seed(seed)
    brms_sims <- brms::posterior_predict(
                           brms_model, re_formula = NULL, newdata = new_data,
                           allow_new_levels = TRUE, sample_new_levels = "gaussian"
                       ) |>
        t()
    fitted_median_brms <- apply(brms_sims, 1, median)
    ## fitted_median_brms <- apply(
    ##     t(brms::posterior_epred(brms_model, ndraws = ndraws, re.form = NA)),
    ##     1,
    ##     mean)
    DHARMa::createDHARMa(
                simulatedResponse = brms_sims,
                observedResponse = response,
                fittedPredictedResponse = fitted_median_brms,
                ...
            )
}
