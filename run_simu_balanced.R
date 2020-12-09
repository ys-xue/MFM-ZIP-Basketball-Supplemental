## read in basis functions
## natural log and normalize

basis_functions <- readRDS("nmf_basis.rds")
nmf_basis <- log(basus_functions)
for (i in 1:ncol(nmf_basis)) {
  nmf_basis[,i] <- (nmf_basis[, i] - mean(nmf_basis[, i])) / sd(nmf_basis[, i])
}

X <- rbind(1, t(nmf_basis))


ZIPcode <- nimbleCode({
  p ~ dunif(0,1)
  lambda ~ dunif(0,10)
  for(i in 1:N){
    y[i] ~ dZIP(lambda,zeroProb = p)
  }
})

dZIP <- nimbleFunction(
  run = function(x = integer(), lambda = double(), zeroProb = double(), 
                 log = logical(0, default = 0)) {
    returnType(double())
    ## First handle non-zero data
    if(x != 0) {
      ## return the log probability if log = TRUE
      if(log) return(dpois(x, lambda, log = TRUE) + log(1-zeroProb))
      ## or the probability if log = FALSE
      else return((1-zeroProb) * dpois(x, lambda, log = FALSE))
    }
    ## From here down we know x is 0
    totalProbZero <- zeroProb + (1-zeroProb) * dpois(0, lambda, log = FALSE)
    if(log) return(log(totalProbZero))
    return(totalProbZero)
  })

rZIP <- nimbleFunction(
  run = function(n = integer(), lambda = double(), zeroProb = double()) {
    returnType(integer())
    isStructuralZero <- rbinom(1, prob = zeroProb, size = 1)
    if(isStructuralZero) return(0)
    return(rpois(1, lambda))
  })

registerDistributions(list(
  dZIP = list(
    BUGSdist = "dZIP(lambda, zeroProb)",
    discrete = TRUE,
    range = c(0, Inf),
    types = c('value = integer()', 'lambda = double()', 'zeroProb = double()')
  )))


piFun2 <- nimbleFunction(
  run = function(r = double(1)) {
    rlength <- length(r)
    rsum <- rep(0, rlength)
    pi <- rep(0, rlength)
    rsum[1] <- pi[1] <- r[1]
    for (i in 2:rlength) {
      rsum[i] <- rsum[i - 1] + r[i]
      if (rsum[i] >= 1) {
        pi[i] <- 1 - rsum[i - 1]
      }
      else {pi[i] <- r[i]}
    }
    for (i in 1:rlength) {
      if (pi[i] < 0) {pi[i] <- 0}
    }
    returnType(double(1))
    return(pi)
  }
)


MainCode <- nimbleCode({
  for (s in 1:S) {
    # S number of  players
    for (j in 1:Q) {
      Y[s, j] ~ dZIP(lambda[s, j], zeroProb = p.ZIP[s])
      lambda[s, j] <- exp(b[s, 1] * x1[j] + b[s, 2] * x2[j] + 
                            b[s, 3] * x3[j] + b[s, 4] * x4[j] +
                            b[s, 5] * x5[j] + b[s, 6] * x6[j])
    }
    b[s, 1:6] <- bm[latent[s], 1:6]
    p.ZIP[s] <- pm[latent[s]]
    latent[s] ~ dcat(pi_beta[1:M])
  }
  
  ### clustered coefficients: MFM ###
  for (i in 1:M) {
    r_beta[i] ~ dexp(rate = lambda_beta)
  }
  pi_beta[1:M] <- piFun2(r = r_beta[1:M])
  
  lambda_beta ~ dlnorm(0, varlog = 1)
  
  for (k in 1:M) {
    bm[k, 1:6] ~ dmnorm(mu_bm[1:6], cov = var_bm[1:6, 1:6])
    pm[k] ~ dunif(0, 1)
  }
})


sim <- function(seed, Y, X, iter, burn_in) {
  set.seed(seed)
  S <- dim(Y)[1]
  Q <- dim(X)[2]
  
  MFMdata <- list(Y = Y,
                  x1 = X[1,], x2 = X[2,], x3 = X[3,], x4 = X[4,], x5 = X[5,],
                  x6 = X[6,]

  
  MFMConsts <- list(S = dim(Y)[1], 
                    Q = Q, 
                    M = 10,
                    mu_bm = rep(0,6),
                    var_bm = diag(1,6))
  
  MFMInits <- list(latent = rep(1,MFMConsts$S),
                   pi_beta = rep(0.1, MFMConsts$M),
                   r_beta = rep(0.2,MFMConsts$M),
                   lambda_beta = 1,
                   pm = rep(0.1,MFMConsts$M),
                   bm = matrix(1,MFMConsts$M,6)) 
  
  
  
  MFM <- nimbleModel(code = MainCode, name = "MFM", constants = MFMConsts,
                     data = MFMdata, inits = MFMInits)
  cMFM <- compileNimble(MFM)
  MFMconf <- configureMCMC(MFM, print = FALSE)
  
  MFMconf$addMonitors(c( "b",
                         "bm",
                         "p.ZIP",
                         "pm",
                         "latent",
                         "pi_beta"
  ))
  
  MFMmcmc <- buildMCMC(MFMconf)
  cMFMmcmc <- compileNimble(MFMmcmc, project = MFM)
  cMFM$setInits(MFMInits)
  mcmc.out <- runMCMC(cMFMmcmc, niter = iter, setSeed = seed, thin = 1)
  pos_mcmc <- as.mcmc(mcmc.out[-c(1:burn_in),])
  return(pos_mcmc)
}


## load simulated data
Y <- readRDS("simulated_balance.rds")

result <- purrr::map(1:100, ~sim(seed = 1, Y = Y[, , .x], X = X, iter = 7000,
                                 burn_in = 2000))
latentColumns <- 512:586

getDahl <- function(latent_iter) {
    membership_matrices <- purrr::map(1:nrow(latent_iter),
                                ~outer(latent_iter[, .x],
                                       latent_iter[, .x], "=="))
    avg_membership <- Reduce("+", membership_matrices) / nrow(latent_iter)
    l2Dist <- purrr::map_dbl(membership_matrices,
                            ~sum((.x - avg_membership)^2))
    finalCluster <- as.numeric(latent_iter[which.min(l2Dist)])
}

latentList <- purrr::map(result, ~getDahl(.x[, latentColumns]))
randValues <- purrr::map(latentList,
                        ~fossil::rand.index(.x, rep(c(1, 2, 3), each = 25)))
