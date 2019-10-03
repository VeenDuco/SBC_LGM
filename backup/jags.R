library(blavaan)

model <- '
I =~ 1 * x1 + 1 * x2 + 1 * x3
S =~ .083 * x1 + .25 * x2 + 1 * x3
I ~ prior("dnorm(35, 0.1)") * 1
S ~ prior("dnorm(-5, 0.2)") * 1

I ~~ I
S ~~ S
I ~~ S

x1 ~~ c(a) * x1
x2 ~~ c(a) * x2
x3 ~~ c(a) * x3

'

X <- data.frame(x1 = prior_y[ , 1, 1],
                x2 = prior_y[ , 2, 1],
                x3 = prior_y[ , 3, 1])

dpri <- dpriors(alpha = "dnorm(35, .1)",
                beta = "dnorm(35, .1)",
                itheta = "dlnorm(0, 1)",
                target = "jags")
fit_blav <- blavaan(model, data = X, dp = dpri,
                    target = "jags", mcmcfile = T, test = "none")
parTable(fit_blav)

load("C:/Users/5507553/Dropbox/Werk/PhD UU/SBC_LGM/lavExport/semjags.rda")
jagtrans$data


library(rjags)
reps <- 255
reps <- 4

tempMatrix <- matrix(NA, ncol = 12, nrow = reps)
for(i in 1:reps){
  print(i)
  jm <- jags.model(file = "lavExport/sem.jag", 
                   data = list(x1 = prior_y[ , 1, i],
                               x2 = prior_y[ , 2, i],
                               x3 = prior_y[ , 3, i],
                               g = rep(1, n_ind),
                               N = n_ind),#jagtrans$data,
                   n.adapt = 10000)
  
  
  thin <- 1
  draws <- L
  
  fit_jags <- coda.samples(jm, 
                           variable.names = 
                             c("lambda[1,1,1]", "lambda[2,1,1]", "lambda[3,1,1]", 
                               "lambda[1,2,1]", 
                               "lambda[2,2,1]", "lambda[3,2,1]", "psistar[1,1,1]", 
                               "psistar[2,2,1]", 
                               "psi[1,2,1]", "theta[1,1,1]", "nu[1,1,1]", 
                               "nu[2,1,1]", "nu[3,1,1]", 
                               "alpha[1,1,1]", "alpha[2,1,1]", "beta[1,1,1]", 
                               "beta[2,2,1]", 
                               "psi[3,3,1]", "psi[1,1,1]", "psi[2,2,1]"),
                           n.iter = draws, thin = thin)
  
  at <- autocorr.diag(fit_jags, lags = c(1))
  at <- c(na.omit(c(at)))
  
  while(max(at) > .1) { 
    thin <- thin * 2
    draws <- draws * 2
    fit_jags <- coda.samples(jm, 
                             variable.names = jagtrans$monitors[-c(11, 12)],
                             n.iter = draws, thin = thin)
    
    at <- autocorr.diag(fit_jags, lags = c(1))
    at <- c(na.omit(c(at)))
  }
  
  
  fit_alpha1 <- fit_jags[, "alpha[1,1,1]"][[1]]
  fit_alpha2 <- fit_jags[, "alpha[2,1,1]"][[1]]
  fit_theta <- fit_jags[, "theta[1,1,1]"][[1]]
  fit_psi11 <- fit_jags[, "psi[1,1,1]"][[1]]
  fit_psi22 <- fit_jags[, "psi[2,2,1]"][[1]]
  fit_psi21 <- fit_jags[, "psi[1,2,1]"][[1]]
  
  
  alpha1_bias <- summary(fit_jags)[[1]]["alpha[1,1,1]", 1] - prior_alpha1[i]
  alpha2_bias <- summary(fit_jags)[[1]]["alpha[2,1,1]", 1] - prior_alpha2[i]
  theta_bias <- summary(fit_jags)[[1]]["theta[1,1,1]", 1] - prior_theta[i]
  psi11_bias <- summary(fit_jags)[[1]]["psi[1,1,1]", 1] - prior_psi11[i]
  psi22_bias <- summary(fit_jags)[[1]]["psi[2,2,1]", 1] - prior_psi22[i]
  psi21_bias <- summary(fit_jags)[[1]]["psi[1,2,1]", 1] - prior_psi21[i]
  
  tempMatrix[i, 1] <- sum(fit_alpha1 < prior_alpha1[i])
  tempMatrix[i, 2] <- sum(fit_alpha2 < prior_alpha2[i])
  tempMatrix[i, 3] <- sum(fit_theta < prior_theta[i])
  tempMatrix[i, 4] <- sum(fit_psi11 < prior_psi11[i])
  tempMatrix[i, 5] <- sum(fit_psi22 < prior_psi22[i])
  tempMatrix[i, 6] <- sum(fit_psi21 < prior_psi21[i])
  
  tempMatrix[i, 7] <- alpha1_bias
  tempMatrix[i, 8] <- alpha2_bias
  tempMatrix[i, 9] <- theta_bias
  tempMatrix[i, 10] <- psi11_bias
  tempMatrix[i, 11] <- psi22_bias
  tempMatrix[i, 12] <- psi21_bias
  
}

  tempMatrix
  