#' @title Run Dail-Madsen Generalized N-mixture Model in JAGS
#'
#' @description
#' \code{runJagsModel} samples jags in parallel chunks with rjags jags.samples instead of coda.samples.
#'
#' @param model Name of model to use
#' @param data List of data for rjags
#' @param inits Function for initializing starting values
#' @param n.adapt Length of adaptation/burnin phase
#' @param params List of parameters to monitor
#' @param n.iter Number of iterations after burnin
#' @param thin Thinning rate (not actually a rate)
#' @param export.list List of objects to make available to clusters for parallel processing
#' 
#' @details
#' blah, blah, blah
runJagsModel <- function(model = "src/Yoichiro-DM-Model.R",  data, inits, export.list, n.chains = 3, n.adapt=1000, params, n.iter = 1000, thin = 3, jags.out = FALSE) {

if(jags.out) {
  CL <- makeCluster(3)
  clusterExport(cl=CL, export.list, envir = environment())
  clusterSetRNGStream(cl=CL, iseed = 2345642)
  
  system.time(out <- clusterEvalQ(CL, {
    library(rjags)
    load.module('glm')
    jm <- jags.model(paste0(model), data=data, inits=inits, n.adapt=n.burn, n.chains=1)
    fm <- jags.samples(jm, params, n.iter = n.it, thin = n.thin)
    return(fm)
  }))
  
  model.out <- mcmc.list(out)
  stopCluster(CL)
} else {
  CL <- makeCluster(3)
  clusterExport(cl=CL, export.list, envir = environment())
  clusterSetRNGStream(cl=CL, iseed = 2345642)
  
  system.time(out <- clusterEvalQ(CL, {
    library(rjags)
    load.module('glm')
    jm <- jags.model(paste0(model), data=Dat, inits=Init, n.adapt=n.burn, n.chains=1)
    fm <- coda.samples(jm, params, n.iter = n.it, thin = n.thin)
    return(as.mcmc(fm))
  }))
  
  model.out <- mcmc.list(out)
  stopCluster(CL)
}

return(model.out)
}



