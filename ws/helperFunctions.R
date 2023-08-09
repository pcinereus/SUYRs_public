## SUYR_prior_and_posterior_old <- function(mod) {
##     dat <- mod$data
##     terms <- attr(dat, 'terms')
##     response <- all.vars(update(terms, .~ 1))
##     predictors <- all.vars(terms)[-1]

##     Xmat <- dat %>%
##         dplyr::select(any_of(predictors)) %>%
##         colMeans() %>%
##         as.matrix()
##     b <- mod %>%
##         as_draws_df() %>%
##         dplyr::select(starts_with('b_'),
##                       -contains('Intercept')) %>%
##         as.matrix()

##     scal <- as.vector(Xmat %*% t(b))

##     priors <- mod %>% get_variables() %>% str_subset("^prior_.*") 
##     pars <- mod %>% get_variables() %>%
##         str_subset(paste0("^b_(Intercept|",paste0(predictors,collapse="|"),")"))
##     aux <- priors %>% str_subset("prior_(Intercept|b)", negate = TRUE) %>%
##         str_remove("^prior_")
##     pars <- c(pars, aux)
##     ## ## nms <- mod %>% get_variables() %>% str_detect(pars)

##     mod.pp <- mod %>%
##         as_draws_df() %>%
##         select(any_of(c(pars, priors))) %>%
##         mutate(b_Intercept = b_Intercept + scal) %>%
##         pivot_longer(cols=everything(), names_to='key', values_to='value') %>% 
##         mutate(Type = ifelse(str_detect(key, 'prior'), 'Prior', 'Posterior'),
##                Parameter = ifelse(Type == 'Prior',
##                                   str_remove(key, "^prior_"),
##                                   str_remove(key, "^b_")
##                                   ),
##                Class = ifelse(Parameter %in% brms_names, 'b', Parameter))

##     return(
##         ggplot(data = NULL, aes(x=Type,  y=value)) +
##         stat_pointinterval(data = mod.pp %>% filter(Class != 'b' | Parameter == 'b')) +
##         stat_pointinterval(data = mod.pp %>% filter(Class == 'b' & Parameter != 'b'),
##                            aes(colour = Parameter))+
##         facet_wrap(~Class,  scales='free')
##         )
## }


## SUYR_prior_and_posterior <- function(mod) {
##     dat <- mod$data
##     terms <- attr(dat, 'terms')
##     response <- all.vars(update(terms, .~ 1))
##     predictors <- all.vars(terms)[-1]
##     rhs <- mod$formula %>% as.formula %>% brms:::str_rhs() #%>%
##         ## str_split('\\+', simplify = TRUE) %>%
##         ## str_trim()

##     Xmat <- model.matrix(as.formula(paste0('~',rhs)), dat)[,-1] %>%
##         colMeans() %>%
##         as.matrix()
    
##     ## Xmat <- dat %>%
##     ##     select(any_of(rhs)) %>%
##     ##     summarise(across(is.numeric, mean)) %>%
##     ##     as.matrix() 

##     b <- mod %>%
##         as_draws_df() %>%
##         dplyr::select(starts_with('b_'),
##                       -contains('Intercept')) %>%
##         as.matrix()
    
##     scal <- as.vector(t(Xmat) %*% t(b))
    
##     brms_names <- brms:::change_effects(brmsterms(mod$formula),
##                                         data = dat,
##                                         pars = variables(mod))
##     brms_names <- sapply(brms_names, function(x) str_remove(x$fnames, "b_"))
##     priors <- mod %>% get_variables() %>% str_subset("^prior_.*") 
##     pars <- mod %>% get_variables() %>%
##         str_subset(paste0("^b_(Intercept|",paste0(brms_names,collapse="|"),")"))
##     aux <- priors %>% str_subset("prior_(Intercept|b)", negate = TRUE) %>%
##         str_remove("^prior_")
##     pars <- c(pars, aux)

##     mod.pp <- mod %>%
##         as_draws_df() %>%
##         select(any_of(c(pars, priors))) %>%
##         mutate(b_Intercept = b_Intercept + scal) %>%
##         pivot_longer(cols=everything(), names_to='key', values_to='value') %>% 
##         mutate(Type = ifelse(str_detect(key, 'prior'), 'Prior', 'Posterior'),
##                Parameter = ifelse(Type == 'Prior',
##                                   str_remove(key, "^prior_"),
##                                   str_remove(key, "^b_")
##                                   ),
##                Class = ifelse(Parameter %in% brms_names, 'b', Parameter))

##     return(
##         ggplot(data = NULL, aes(x=Type,  y=value)) +
##         stat_pointinterval(data = mod.pp %>% filter(Class != 'b' | Parameter == 'b')) +
##         stat_pointinterval(data = mod.pp %>% filter(Class == 'b' & Parameter != 'b'),
##                            aes(colour = Parameter))+
##         facet_wrap(~Class,  scales='free')
##         )
## }

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
    brms_names <- brms:::change_effects(brmsterms(mod$formula),
                                        data = dat,
                                        pars = variables(mod))
    brms_names <- sapply(brms_names, function(x) str_remove(x$fnames, "b_"))
    priors <- mod %>% get_variables() %>% str_subset("^prior_.*") %>% str_subset("lprior", negate = TRUE) 
    pars <- mod %>% get_variables() %>%
        str_subset(paste0("^b_(Intercept|",paste0(brms_names,collapse="|"),")"))
    
    ## auxillary
    aux <- priors %>% str_subset("prior_(Intercept|b)", negate = TRUE) %>%
        str_remove("^prior_")
    pars <- c(pars, aux)
    ## random effects
    if (length(wch.rnd)>0) {
        ran.pars <- brms:::change_re(mod$ranef, pars = variables(mod))[[1]]$fnames
        pars <- c(pars, ran.pars)
    }
    variables(mod)

    vars <- variables(mod)
    priors <- vars %>% str_subset("prior") %>% str_subset("lprior", negate = TRUE)
    all.pars <- priors %>% str_remove("prior_")
    fixed.pars <- vars %>% str_subset("^b_")
    other.pars <- all.pars %>% str_subset("^Intercept$|^b$", negate = TRUE)
    other.pars <- vars %>% str_subset(paste0("^", other.pars, collapse = '|'))
    pars <- c(fixed.pars, other.pars)
    
    ## coefs <- prior_summary(mod)$class %>% unique()
    ## coefs.regex <- paste0("^b_", coefs, collapse = "|")
    
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
               Class = ifelse(Parameter %in% brms_names, 'b', Parameter),
               Class = ifelse(Type == 'Posterior', str_remove(Class, "__.*"), Class)) %>%
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


posterior_predict.inla <- function(object, newdata = NULL, ndraws=250, include_re = TRUE, new_random_levels = FALSE) {
  draws <- inla.posterior.sample(n=ndraws, result = object)
  contents <- object$misc$configs$contents

  ## get the linear predictor
  if (is.null(newdata)) {
    mm <- object$model.matrix
    ## wch <- which(contents$tag=='Predictor')
    ## i_data <- contents$start[wch]:(1+contents$length[wch])
    ## mm <- mm[i_data,]
  } else {
      ## Fixed effects
      contrasts <- object$.args$contrasts
      form <- object$.args$formula
      gf <- INLA:::inla.interpret.formula(form)
      Xmat <- model.matrix(update(gf$fixf, NULL ~ .), newdata, contrasts = contrasts)

      ## Random effects
      if (!is.null(gf$randf)) {
          nms <- unlist(lapply(gf$random.spec, `[[`, "label"))
          ## newdata <- newdata %>% mutate(across(!!nms, ~ paste0('NEW',.x) %>% as.factor()))  
          Zmat <- lapply(nms, function(n) model.matrix(as.formula(paste0("~ 0 + ",n)), newdata))
          Zmat <- do.call('cbind', Zmat)
      }
  }
  ## Fixed effects
  nms <- colnames(Xmat)
  i_lp <- contents$start[contents$tag %in% nms]
  lp <- t(sapply(draws, function(x) x$latent[i_lp]))
  if (nrow(lp)==1) lp <- t(lp)
  b <- tcrossprod(Xmat, lp)
  link <- object$misc$linkfunctions$names
  out <- my_ilink(b, link) 

  ## family theta
  if ((ntheta <- object$misc$configs$ntheta) > 0) {
    theta <- vector('list', ntheta)
    for (t in 1:ntheta) {
      theta[[t]] <- sapply(draws, function(x) 1/sqrt(x$hyperpar[[t]]))
    }
  }

  ## Random effects
  if (gf$n.random>0 & include_re) {
      nms <- unlist(lapply(gf$random.spec, `[[`, "label"))
      i_lp <- contents$start[contents$tag %in% nms]
      lp <- ii_lp <- vector('list', length(nms))
      for (i in 1:length(nms)) {
          ii_lp[[i]] <- i_lp[[i]]:(i_lp[[i]] + (contents$length[contents$tag %in% nms] - 1)[i])
          lp[[i]] <- t(sapply(draws, function(x) x$latent[ii_lp[[i]]]))
          if (new_random_levels) lp[[i]] <- apply(lp[[i]], 2,
                                                  function(x) rnorm(n = length(x), mean = 0, sd = theta[[1+i]]))
      }
      lp <- do.call(`cbind`, lp)
      ## i_lp <- unlist(lapply(1:length(i_lp),
      ##                       function(i)
      ##                           i_lp[i]:(i_lp[i] + (contents$length[contents$tag %in% nms] - 1)[i])))
      
      ## i_lp <- i_lp:(i_lp + contents$length[contents$tag %in% nms] - 1)
      ## lp <- t(sapply(draws, function(x) x$latent[i_lp]))
      ## if (new_random_levels) lp <- apply(lp, 2, function(x) rnorm(n = length(x), mean = 0, sd = theta[[2]]))
      r <- tcrossprod(Zmat, lp)
      link <- object$misc$linkfunctions$names
      out <- my_ilink(b+r, link) 
  }

  fam <- object$.args$family
  N <- nrow(out)
  if (fam == 'gaussian') {
    rdist <- function(N, x, sigma) rnorm(N, x, sigma)
    out <- apply(out, 2, rdist, N=N, sigma=theta[[1]])
  }
  if (fam=='poisson') {
    rdist <- function(N, x) rpois(N, x)
    out <- apply(out, 2, rdist, N=N)
  }
  if (fam=="nbinomial") {
      wch <- grep(
          "size for the nbinomial observations",
          names(draws[[1]]$hyperpar)
      )
      size <- sapply(draws, function(x) x$hyperpar[[wch]])
      ## size <- inla.hyperpar.sample(n = ndraws, result = mod.inla)[, wch]
      ## rdist <- function(N, x, size) MASS::rnegbin(N, mu = x, size) #rnbinom(N, mu = x, prob)
      rdist <- function(N, x, size) rnbinom(N, mu = x, size) #rnbinom(N, mu = x, prob)
      out <- apply(out, 2, rdist, N=N, size=size)
  }
  if (fam=="binomial") {
    ## size <- ifelse(is.null(object$.args$Ntrials), 1, object$.args$Ntrials)
    if (is.null(object$.args$Ntrials)) {
        size <- rep(1, N)
    } else if (length(object$.args$Ntrials) == 1) {
        size <- rep(object$.args$Ntrials, N)
    } else {
        size <- object$.args$Ntrials[1:N]
    }
    rdist <- function(N, x, size) rbinom(N, prob = x, size)
    out <- apply(out, 2, rdist, N=N, size=size)
  }
  if (fam=="beta") {
      wch <- grep(
          "precision parameter for the beta observations",
          names(draws[[1]]$hyperpar)
      )
      phi <- sapply(draws, function(x) x$hyperpar[[wch]])
      rdist <- function(N, x, phi) rbeta(N, shape1 = x*phi, shape2 = (1 - x)*phi)
      ## out <- apply(out, 2, rdist, N=N, phi=phi)
      out <- sapply(seq_len(ncol(out)), function(i)
          rdist(N = N, x = out[,i], phi = phi[i])) %>%
          matrix(ncol = ncol(out), byrow = TRUE)
  }
   if (fam=="betabinomial") {
      wch <- grep(
          "overdispersion for the betabinomial observations",
          names(draws[[1]]$hyperpar)
      )
      phi <-1/sapply(draws, function(x) x$hyperpar[[wch]])
      if (is.null(object$.args$Ntrials)) {
          size <- rep(1, N)
      } else {
          size <- object$.args$Ntrials[1:N]
      }
      ## size <- ifelse(is.null(object$.args$Ntrials), 1, object$.args$Ntrials)
      rdist <- function(N, x, phi, size) {
          rbinom(N,
                 prob = rbeta(N, shape1 = x*phi, shape2 = (1 - x)*phi),
                 size = size)
          }
      out <- sapply(seq_len(ncol(out)), function(i)
          rdist(N = N, x = out[,i], phi = phi[i], size = size)) %>%
          matrix(ncol = ncol(out), byrow = TRUE)
  }
  if (fam=="zeroinflatedpoisson1") {
      wch <- grep(
          "zero-probability parameter for zero-inflated poisson_1",
          names(draws[[1]]$hyperpar)
      )
      phi <- sapply(draws, function(x) x$hyperpar[[wch]])
      rdist <- function(N, x, phi) rbinom(N, size = 1, prob = 1 - phi) * rpois(N, lambda = x)
      out <- apply(out, 2, rdist, N=N, phi=phi)
  }
  
  
  t(out)
}

my_ilink <- function(x, link) {
  switch(link, identity = x, log = exp(x), logm1 = expp1(x), 
         log1p = expm1(x), inverse = 1/x, sqrt = x^2, `1/mu^2` = 1/sqrt(x), 
         tan_half = 2 * atan(x), logit = inv_logit(x), probit = pnorm(x), 
         cauchit = pcauchy(x), cloglog = inv_cloglog(x), probit_approx = pnorm(x), 
        softplus = log1p_exp(x), stop2("Link '", link, "' not supported."))
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


plot0.1 <- theme(axis.title.x = 
                   element_text(
                     margin = margin(t = 10),
                     color = "Black", 
                     size = 15), 
                 axis.title.y = 
                   element_text(
                     margin = margin(r = 20),
                     color = "Black", 
                     size = 15), 
                 axis.text = 
                   element_text(
                     color = "Black",
                     size = 12), 
                 axis.line = 
                   element_line(
                     color = "Black",
                     linewidth = 0.5), 
                 axis.ticks = element_line(color = "Black"),
                 panel.background = element_rect(fill = "grey98"),
                 panel.grid.minor = element_line(color = "White"),
                 plot.title = 
                   element_text(
                     size = 15,
                     face = "bold",
                     hjust = 0.5),
                 plot.subtitle = 
                   element_text(
                     size = 15,
                     hjust = 0.5),
                 strip.text.x = 
                   element_text(
                     size = 15,
                     color = "black", 
                     face = "bold"),
                 strip.background = 
                   element_rect(
                     color = "white",
                     fill = "white", 
                     linewidth = 1, 
                     linetype = "solid"),
                 legend.background = element_blank(), 
                 legend.title = 
                   element_text(
                     colour = "black",
                     size = 15, 
                     face = "bold"), 
                 legend.text = 
                   element_text(
                     colour = "black",
                     size = 12,
                     face = "bold.italic"),
                 legend.key.width = unit(1.5, "lines"), 
                 legend.key = element_rect(fill = NA))

# Same theme without a legend:
plot0 <- plot0.1 + theme(legend.position = "none")

autocor_check <- function(data,  modelTMB, variable, grouping, maxlag, n.sim) {
  
  data <- as.data.frame(data)
  
  if(missing(variable) || is.null(variable) ||  is.na(variable)) 
  {stop("ERROR: Temporal or spatial predictor missing. Please specify one, otherwise skip this step.")
  } else {
    
    if(length(variable) > 1) {   # dealing with spatial dependency here
      
      dists <- dist(data[, variable]) # this computes distance matrix 
      coord <- cbind(data[, variable[1]], data[, variable[2]])
      
      if(missing(maxlag) || is.null(maxlag) ||  is.na(maxlag)) {
        lags <- seq(min(dists), max(dists), length = 10) # set 10 lags at which to calculate semivariance  
      } else {lags <- seq(min(dists), max(maxlag), length = 10)}
      
      
    } else {
      
      if(!missing(grouping) && !is.null(grouping) &&  !is.na(grouping)) {  # dealing with multiple time series here  
        
        if(missing(maxlag) || is.null(maxlag) ||  is.na(maxlag)) {
          maxlag <- max(as.data.frame(data %>%
                                        group_by_at(grouping) %>%
                                        summarise(lag = max(!!sym(variable)) - min(!!sym(variable))))[,2])}
        
        new_time <- data[, variable] + maxlag * 3 * (as.integer(data[, grouping])-1)
        data$new_time <- new_time
        
        dists <- dist(cbind(new_time, rep(0, nrow(data)))) # this computes distance matrix   
        coord <- cbind(data$new_time, rep(0, nrow(data)))
        
        lags <- seq(min(dists), maxlag, length = 10) # set 10 lags at which to calculate semivariance  
        
      } else { # dealing with single time series here
        
        dists <- dist(cbind(data[, variable], rep(0, nrow(data)))) # this computes distance matrix   
        coord <- cbind(data[, variable], rep(0, nrow(data)))
        
        if(missing(maxlag) || is.null(maxlag) ||  is.na(maxlag)){
          lags <- seq(min(dists), max(dists), length = 10) # set 10 lags at which to calculate semivariance  
        } else {
          lags <- seq(min(dists), maxlag, length = 10)}
      }}
    if (class(modelTMB) %in% "glmmTMB") 
      res <- scale(residuals(modelTMB))
    if (class(modelTMB) %in% "brmsfit") 
      res <- scale(rowMeans(residuals(modelTMB)))
    data.Var <- as.data.frame(variog(coords = coord, 
                                     data = res, 
                                     breaks = lags)[c("u","v","n")])
    
    # Here we permute the observed semivariances across the observed pairwise distances (= grey lines in graph)
    data.rand.Var <- replicate(n.sim, variog(coords = coord, 
                                             data = sample(res, replace = F),
                                             breaks = lags)[c("v")])
    
    names(data.rand.Var) <- as.character(c(1:n.sim))
    data.rand.Var <- do.call(cbind.data.frame, data.rand.Var)
    data.rand.Var <- reshape::melt(data.rand.Var)
    data.rand.Var$u <- data.Var$u
    colnames(data.rand.Var) <- c("rep.run","v.rand","u")
    
    ggplot(data = data.Var, 
           aes(x = u, y = v)) + # this extracts distance lags and semi variance
      geom_line(data = data.rand.Var,
                aes(x = u, y = v.rand, 
                    group = rep.run), col = "grey80", size = 1, alpha = 0.5) +
      geom_point(col = "blue", size = 5, alpha = 0.5) +
      geom_smooth(method = "loess", se = F, col = "blue", size = 1) +
      geom_hline(yintercept = 1, col = "red", lty = 2) +
      xlab("Distance between observations") + 
      ylab("Standardised semivariance") +
      scale_y_continuous() +
      scale_x_continuous(breaks = data.Var$u, labels = formatC(round(data.Var$u, digits = 2), 2, format = "f")) +
      annotate("text", x = data.Var$u, y = min(c(data.rand.Var$v.rand,data.Var$v)) - 0.02*max(c(data.rand.Var$v.rand,data.Var$v)), label = c(paste("N =",data.Var$n[1]), data.Var$n[-1]), size = 5) +
      plot0
  }}




## The following is required to fix a bug in emmeans and brms for poly terms
## https://github.com/rvlenth/emmeans/issues/43
recover_data.brmsfit <- function(object, data, ...) {
    bt = brms::parse_bf(formula(object))
    if (class(bt) != "brmsterms")
        stop("This model is currently not supported.")
    trms <- attr(model.frame(bt$dpars$mu$fe, data = object$data), "terms")
    # we don't have a call component so I'll just put in NULL
    emmeans:::recover_data.call(NULL, trms, "na.omit", data = object$data, ...)
}
