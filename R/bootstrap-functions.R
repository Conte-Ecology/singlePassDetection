#' @title Predictions from glmer model
#'
#' @description
#' \code{myPredict} predictions from glmer model for use in apply function
#'
#' @param model merMod object from glmer model fit
#' @param newdata Dataframe of independent variables to predict on
#' 
#' @details
#' blah, blah, something, something
myPredict <- function(model, newdata = data.full, re.form) {
  preds <- predict(model, newdata = newdata, type = "response", re.form = re.form)
  return(preds)
}

#' @title Bootstrap hierarchical data
#'
#' @description
#' \code{groupResample} hierarchical bootstraps of data
#'
#' @param dat Data to bootstrap
#' @param group Grouping factor for data
#' @param replace Logical whether to sample with replacement or not
#' 
#' @details
#' adapted from http://biostat.mc.vanderbilt.edu/wiki/Main/HowToBootstrapCorrelatedData
groupResample <- function(dat, group, replace) {
  # adapted from http://biostat.mc.vanderbilt.edu/wiki/Main/HowToBootstrapCorrelatedData
  # exit early for trivial data
  if(nrow(dat) == 1 || all(replace==FALSE))
    return(dat)
  
  # sample the grouping factor
  cls <- sample(unique(dat[[group[1]]]), replace=replace[1])
  
  # subset on the sampled grouping factors
  sub <- lapply(cls, function(b) subset(dat, dat[[group[1]]]==b))
  
  # sample lower levels of hierarchy (if any)
  if(length(group) > 1)
    sub <- lapply(sub, resample, group=group[-1], replace=replace[-1])
  
  # join and return samples
  do.call(rbind, sub)
}

#' @title Get bootstrapped predictions from functions
#'
#' @description
#' \code{bootSingleVisit} hierarchical bootstraps of data to get predictions
#'
#' @param x Model output from a regression such as lm, glm, glmer, or svabu
#' @param fit.data Data used to fit model x
#' @param new.data Data used for predictions
#' @param nsim Number of bootstrap simulations
#' @param seed Numeric value for reproduction of random boot strap iterations
#' 
#' @details
#' blah, blah, blah
bootSingleVisit <- function(x, fit.data, new.data, nsim, seed = NULL) {
  pred.vals <- matrix(NA, nrow = dim(new.data)[1], ncol = nsim)
  if (!is.null(seed)) 
    set.seed(seed)
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  for(i in 1:nsim) {
    bootData <- groupResample(fit.data, group = "site", replace = TRUE)
    result <- svabu(formula(x)$full, data = bootData) # x$model) #bootData[[i]])
    pred.vals[ , i] <- exp(predict(result, new.data))
    setTxtProgressBar(pb, i)
  }
  close(pb)
  return(pred.vals)
}

#' @title Get bootstrapped predictions from functions
#'
#' @description
#' \code{parallelBoot} hierarchical bootstraps of data to get predictions
#'
#' @param x Model output from a regression such as lm, glm, glmer, or svabu
#' @param fit.data Data used to fit model x
#' @param new.data Data used for predictions
#' @param FUN Function for predictions
#' @param nsim Number of bootstrap simulations
#' @param seed Numeric value for reproduction of random boot strap iterations
#' @param parallel Run the sampler in parallel using "multicore" on Unix machines or "snow"
#' @param ncpus Number of cores to run in parallel
#' @param cl Cluster socket for parallel
#' 
#' @details
#' blah, blah, blah
parallelBoot <- function (x, fit.data, new.data, FUN, nsim = 1, seed = NULL, 
                          parallel = c("no", "multicore", "snow"), 
                          ncpus = getOption("boot.ncpus", 1L), cl = NULL) {
  stopifnot((nsim <- as.integer(nsim[1])) > 0)
  if (missing(parallel)) {
    parallel <- getOption("boot.parallel", "no")
    parallel <- match.arg(parallel)
    have_mc <- have_snow <- FALSE
  }
  if (parallel != "no" && ncpus > 1L) {
    if (parallel == "multicore") {
      have_mc <- .Platform$OS.type != "windows"
    } else if (parallel == "snow") {
      have_snow <- TRUE
    }
    if (!have_mc && !have_snow) {
      ncpus <- 1L
    }
  }
  do_parallel <- (ncpus > 1L && (have_mc || have_snow))
  FUN <- match.fun(FUN)
  if (!is.null(seed)) 
    set.seed(seed)
  else if (!exists(".Random.seed", envir = .GlobalEnv)) 
    runif(1)
  mc <- match.call()
  t0 <- FUN(x, fit.data = fit.data, new.data = new.data, nsim = nsim, seed = seed)
  if (!is.numeric(t0)) 
    stop("parallelBoot currently only handles functions that return numeric vectors")

  fit.data <- fit.data
  new.data <- new.data
  nsim = nsim
  ffun <- local({
    new.data
    fit.data
    FUN
    bootSingleVisit
    x
    do_parallel
  })
  simvec <- seq_len(nsim)
  res <- if (do_parallel) {
    if (have_mc) {
      parallel::mclapply(simvec, ffun, mc.cores = ncpus, fit.data, new.data, nsim, seed)
    }
    else if (have_snow) {
      if (is.null(cl)) {
        cl <- parallel::makePSOCKcluster(rep("localhost", 
                                             ncpus))
        parallel::clusterExport(cl, varlist = getNamespaceExports("parallelBoot"), envir = environment())
        if (RNGkind()[1L] == "L'Ecuyer-CMRG") 
          parallel::clusterSetRNGStream(cl)
        res <- parallel::parLapply(cl, simvec, ffun, fit.data, new.data, nsim, seed)
        parallel::stopCluster(cl)
        res
      }
      else parallel::parLapply(cl, simvec, ffun, fit.data, new.data, nsim, seed)
    }
  }
  else lapply(simvec, ffun, fit.data, new.data, nsim, seed)
  t.star <- do.call(cbind, res)
  rownames(t.star) <- names(t0)
  if ((numFail <- sum(apply(is.na(t.star), 2, all))) > 0) {
    warning("some bootstrap runs failed (", numFail, "/", 
            nsim, ")")
  }
  s <- structure(list(t0 = t0, t = t(t.star), R = nsim, data = old.data, 
                      seed = .Random.seed, statistic = FUN, sim = "parametric", 
                      call = mc), 
                 class = "boot")
  attr(s, "bootFail") <- numFail
  s
}


parallelBoot(yoy.detect.1, fit.data = pass1, new.data = data.full, FUN = bootSingleVisit, nsim = 3, seed = 10101, parallel = "multicore", ncpus = 3)
