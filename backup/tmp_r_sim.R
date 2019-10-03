# LGM SBC:
library(rstan)
library(foreach)
library(parallel)
library(doParallel)
library(bayesplot)
library(tcltk)
library(dplyr)
library(brms)
library(matrixcalc)

sbc_rank <- function(ranks, reps, L, title = NULL){
  
  rank_df = data.frame(ranks = ranks)
  
  lb <- qbinom(.005, reps, (L + 1)^-1)
  mean <- qbinom(.5, reps, (L + 1)^-1)
  ub <- qbinom(.995, reps, (L + 1)^-1)
  
  ggplot(rank_df, aes(x = ranks))  + 
    geom_segment(aes(x = -0.5, y = mean, xend = L - 0.5, yend = mean), colour = "grey25") + 
    geom_polygon(data = data.frame(x = c(-1, -0.5, -1, L + 1, L + 0.5, L + 1),
                                   y=c(lb, mean, ub, ub, mean, lb)), 
                 aes(x = x, y = y), fill="grey45", color="grey25", alpha=0.5) +
    geom_histogram(binwidth = 1, fill="#A25050",colour="black") + 
    theme_default() + xlab("Rank Statistic") + ylab("") +
    if(!is.null(title)) ggtitle(title) else ggtitle("")
  
}

estimate <- stan_model("stan_models/LGM_estimate_wish.stan")
# detectCores()
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

n_ind <- 10
n_time <- 3
# informative <- TRUE
reps <- 255
bins <- 8
L <- bins - 1

set.seed(1236)
# sigma1 <- rlnorm(reps, 0, 1)
# sigma2 <- rlnorm(reps, 0, 1)
# cov_sig <- runif(reps, -1, 1)
# PSI <- array(NA, dim = c(2, 2, reps))
for(i in 1:reps){
  PSI <- rWishart(reps, 3, diag(2))
  # PSI[ , , i] <- matrix(c(sigma1[i], cov_sig[i], cov_sig[i], sigma2[i]), ncol = 2)
  # print(is.positive.definite(PSI[ , , i]))
}

alpha1 <- rnorm(reps, 35, 10)
alpha2 <- rnorm(reps, -5, 5)
theta <- rlnorm(reps, 0, 1)

tau <- array(NA, dim = c(n_ind, 2, reps))
for(i in 1:reps){
  tau[ , , i] <- mvtnorm::rmvnorm(n = n_ind, mean = c(alpha1[i], alpha2[i]), sigma = PSI[ , , i])
}

ysim <- array(NA, dim = c(n_ind, n_time, reps))

time <- c(0.083, 0.25, 1)
# set.seed(123)
for(i in 1:n_ind){
  for(j in 1:n_time){
    for(k in 1:reps){
      ysim[i, j, k] <- rnorm(1, mean = I(1 * tau[i, 1, k] + time[j] * tau[i, 2, k]),
                              sd = theta[k])
    }
  }
}
# 


# ysim

prior_y <- ysim

prior_alpha1 <- alpha1
prior_alpha2 <- alpha2
prior_theta <- theta
prior_psi11 <- PSI[1, 1, ]
prior_psi22 <- PSI[2, 2, ]
prior_psi21 <- PSI[2, 1, ]

# data.frame(a1 = prior_alpha1,
#            a2 = prior_alpha2,
#            th = prior_theta,
#            p11 = prior_psi11,
#            p22 = prior_psi22,
#            p21 = prior_psi21) %>%
#   ggplot() + 
#   # geom_density(aes(x = a1, fill = "alpha1"), alpha = .3) +
#   # geom_density(aes(x = a2, fill = "alpha2"), alpha = .3) +
#   # geom_density(aes(x = th, fill = "theta"), alpha = .3) +
#   geom_density(aes(x = p11, fill = "psi11"), alpha = .3) +
#   geom_density(aes(x = p22, fill = "psi22"), alpha = .3) +
#   geom_density(aes(x = p21, fill = "psi21"), alpha = .3) + 
#   theme_classic()


set.seed(123)
warm.init <- 1000

finalMatrix  <- foreach(i = 1:reps, .packages = c("rstan", "tcltk"), 
                        .export = c("L", "warm.init", "reps",
                                    "prior_y", 
                                    "prior_alpha1", "prior_alpha2",
                                    "prior_psi11", "prior_psi22", "prior_psi21",
                                    "prior_theta"),
                        .combine = rbind,
                        .verbose = TRUE) %dopar% {
                          
                          if(!exists("pb")) pb <- tkProgressBar("Parallel task", min = 1, max = reps)
                          setTkProgressBar(pb, i)  
                          
                          # w_init <- rnorm(sim.data$D)
                          # w0_init <- rnorm(1)
                          # 
                          # init_fun <- function(n) {
                          #   list(w = w_init,
                          #        w0 = w0_init)
                          # }
                          
                          thin <- 1
                          draws <- L
                          warm <- warm.init
                          fit <- sampling(estimate, chains = 1, warmup = warm, iter = warm + draws, 
                                          data = list(n_ind = n_ind,
                                                      n_time = n_time,
                                                      y = prior_y[, , i],
                                                      x = c(.083, .25, 1)))
                          
                          
                          ## note 2 * L, was L but still signs of autocorrelation.
                          while(min(na.omit(summary(fit)$summary[, "n_eff"])) < 2 * L | max(na.omit(summary(fit)$summary[, "Rhat"])) > 1.02) { 
                            thin <- thin * 2
                            draws <- draws * 2
                            warm <- warm * 1.1
                            fit <- sampling(estimate, chains = 1, warmup = warm, iter = warm + draws, 
                                            data = list(n_ind = n_ind,
                                                        n_time = n_time,
                                                        y = prior_y[, , i],
                                                        x = c(.083, .25, 1)))
                          }
                          
                          # saveRDS(fit, paste0("temp/fit_", i ,".rds"))
                          
                          fit_alpha1 <- rstan::extract(fit, pars = "alpha1")[[1]]
                          fit_alpha2 <- rstan::extract(fit, pars = "alpha2")[[1]]
                          fit_theta <- rstan::extract(fit, pars = "theta")[[1]]
                          fit_psi11 <- rstan::extract(fit, pars = "Psi[1,1]")[[1]]
                          fit_psi22 <- rstan::extract(fit, pars = "Psi[2,2]")[[1]]
                          fit_psi21 <- rstan::extract(fit, pars = "Psi[2,1]")[[1]]
                          
                          
                          fit_alpha1 <- fit_alpha1[(1:L) * thin]
                          fit_alpha2 <- fit_alpha2[(1:L) * thin]
                          fit_theta <- fit_theta[(1:L) * thin]
                          fit_psi11 <- fit_psi11[(1:L) * thin]
                          fit_psi22 <- fit_psi22[(1:L) * thin]
                          fit_psi21 <- fit_psi21[(1:L) * thin]
                          
                          ## bias
                          alpha1_bias <- summary(fit, pars = "alpha1")$summary[1] - prior_alpha1[i]
                          alpha2_bias <- summary(fit, pars = "alpha2")$summary[1] - prior_alpha2[i]
                          theta_bias <- summary(fit, pars = "theta")$summary[1] - prior_theta[i]
                          psi11_bias <- summary(fit, pars = "Psi[1,1]")$summary[1] - prior_psi11[i]
                          psi22_bias <- summary(fit, pars = "Psi[2,2]")$summary[1] - prior_psi22[i]
                          psi21_bias <- summary(fit, pars = "Psi[2,1]")$summary[1] - prior_psi21[i]
                          
                          tempMatrix <- matrix(NA, ncol = 12, nrow = 1)
                          
                          tempMatrix[1, 1] <- sum(fit_alpha1 < prior_alpha1[i])
                          tempMatrix[1, 2] <- sum(fit_alpha2 < prior_alpha2[i])
                          tempMatrix[1, 3] <- sum(fit_theta < prior_theta[i])
                          tempMatrix[1, 4] <- sum(fit_psi11 < prior_psi11[i])
                          tempMatrix[1, 5] <- sum(fit_psi22 < prior_psi22[i])
                          tempMatrix[1, 6] <- sum(fit_psi21 < prior_psi21[i])
                          
                          tempMatrix[1, 7] <- alpha1_bias
                          tempMatrix[1, 8] <- alpha2_bias
                          tempMatrix[1, 9] <- theta_bias
                          tempMatrix[1, 10] <- psi11_bias
                          tempMatrix[1, 11] <- psi22_bias
                          tempMatrix[1, 12] <- psi21_bias
                          
                          tempMatrix
                          
                        }


rank_alpha1 <- finalMatrix[, 1]
rank_alpha2 <- finalMatrix[, 2]
rank_theta <- finalMatrix[, 3]
rank_psi11 <- finalMatrix[, 4]
rank_psi22 <- finalMatrix[, 5]
rank_psi21 <- finalMatrix[, 6]

sbc_alpha1 <- sbc_rank(rank_alpha1, reps = reps, L = L, title = "Alpha 1")
sbc_alpha2 <- sbc_rank(rank_alpha2, reps = reps, L = L, title = "Alpha 2")
sbc_theta <- sbc_rank(rank_theta, reps = reps, L = L, title = "Theta")
sbc_psi11 <- sbc_rank(rank_psi11, reps = reps, L = L, title = "Psi 11")
sbc_psi22 <- sbc_rank(rank_psi22, reps = reps, L = L, title = "Psi 22")
sbc_psi21 <- sbc_rank(rank_psi21, reps = reps, L = L, title = "Psi 21")


bias_alpha1 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X7)) + theme_classic() + xlab(expression(paste("bias ", alpha[1])))

bias_alpha2 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X8)) + theme_classic() + xlab(expression(paste("bias ", alpha[2])))

bias_theta <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X9)) + theme_classic() + xlab(expression(paste("bias ", theta)))

bias_psi11 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X10)) + theme_classic() + xlab(expression(paste("bias ", psi[11])))

bias_psi22 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X11)) + theme_classic() + xlab(expression(paste("bias ", psi[22])))

bias_psi21 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X12)) + theme_classic() + xlab(expression(paste("bias ", psi[21])))



gridExtra::grid.arrange(sbc_alpha1, sbc_alpha2, sbc_theta,
                        sbc_psi11, sbc_psi22, sbc_psi21, ncol = 3)

gridExtra::grid.arrange(bias_alpha1, bias_alpha2, bias_theta,
                        bias_psi11, bias_psi22, bias_psi21, ncol = 3)




pdf(paste0("results/plots_wish_", n_ind, "ind_", n_time, "time_", reps, "reps", ".pdf"))
gridExtra::grid.arrange(sbc_alpha1, sbc_alpha2, sbc_theta,
                        sbc_psi11, sbc_psi22, sbc_psi21, ncol = 3)

gridExtra::grid.arrange(bias_alpha1, bias_alpha2, bias_theta,
                        bias_psi11, bias_psi22, bias_psi21, ncol = 3)

dev.off()

saveRDS(finalMatrix, 
        paste0("results/finalMatrix_wish_", n_ind, "ind_", n_time, "time_", reps, "reps", ".rds"))




######### JAGS
library(blavaan)
library(rjags)

# model <- '
# I =~ 1 * x1 + 1 * x2 + 1 * x3
# S =~ .083 * x1 + .25 * x2 + 1 * x3
# I ~ prior("dnorm(35, 0.1)") * 1
# S ~ prior("dnorm(-5, 0.2)") * 1
# 
# I ~~ I
# S ~~ S
# I ~~ S
# 
# x1 ~~ c(a) * x1
# x2 ~~ c(a) * x2
# x3 ~~ c(a) * x3
# 
# '
# 
# X <- data.frame(x1 = prior_y[ , 1, 1],
#                 x2 = prior_y[ , 2, 1],
#                 x3 = prior_y[ , 3, 1])
# 
# dpri <- dpriors(alpha = "dnorm(35, .1)",
#                 beta = "dnorm(35, .1)",
#                 itheta = "dlnorm(0, 1)",
#                 target = "jags")
# fit_blav <- blavaan(model, data = X, dp = dpri,
#                     target = "jags", mcmcfile = T, test = "none")
# parTable(fit_blav)

# load("C:/Users/5507553/Dropbox/Werk/PhD UU/SBC_LGM/lavExport/semjags.rda")
# jagtrans$data


# reps <- 255
# reps <- 32
set.seed(123)
# tempMatrix column 13 to show how many times we need to run before
# we satisfy rhat < 1.02
tempMatrix <- matrix(NA, ncol = 13, nrow = reps)
for(i in 1:reps){
  print(i)
  jm <- jags.model(file = "jags_models/LGM_estimate_wish.jag", 
                   data = list(x1 = prior_y[ , 1, i],
                               x2 = prior_y[ , 2, i],
                               x3 = prior_y[ , 3, i],
                               g = rep(1, n_ind),
                               N = n_ind),#jagtrans$data,
                   n.adapt = 10000)
  
  
  thin <- 1
  draws <- L
  multiply <- 1
  counter <- 1
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
                           n.iter = draws, thin = 1)
  
  # at <- autocorr.diag(fit_jags, lags = c(1))
  # at <- c(na.omit(c(at)))
  
  
  rhat <- max(data.frame(fit_jags[[1]]) %>%
    apply(., 2, rstan::Rhat))
  # ess_bulk <- min(data.frame(fit_jags[[1]]) %>%
  #   apply(., 2, rstan::ess_bulk))
  # ess_tail <- min(data.frame(fit_jags[[1]]) %>%
  #   apply(., 2, rstan::ess_tail))
  
  
  # while(max(at) > .03) {
  while(rhat > 1.01 & counter < 60){#| ess_bulk < 2 * L | ess_tail < 2 * L) {
    counter <- counter + 1
    multiply <- multiply * 2
    if(multiply < 4096){
      thin <- thin * 2
      draws <- draws * 2
    } 
    if((counter %% 50) == 0) {
      cat("\014") 
      jm <- jags.model(file = "jags_models/LGM_estimate_wish.jag", 
                       data = list(x1 = prior_y[ , 1, i],
                                   x2 = prior_y[ , 2, i],
                                   x3 = prior_y[ , 3, i],
                                   g = rep(1, n_ind),
                                   N = n_ind),#jagtrans$data,
                       n.adapt = 10000 + (counter * 20))
  
    }
    
    
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
                             n.iter = draws, thin = 1)
    
    # at <- autocorr.diag(fit_jags, lags = c(1))
    # at <- c(na.omit(c(at)))
    rhat <- max(data.frame(fit_jags[[1]]) %>%
                  apply(., 2, rstan::Rhat))
    print(paste0("rep: ",i, ", attempt: ",
                 counter, ", rhat: ",
                 round(rhat, 2)))
    # ess_bulk <- min(data.frame(fit_jags[[1]]) %>%
    #                   apply(., 2, rstan::ess_bulk))
    # ess_tail <- min(data.frame(fit_jags[[1]]) %>%
    #                   apply(., 2, rstan::ess_tail))
  }
  
  
  fit_alpha1 <- fit_jags[, "alpha[1,1,1]"][[1]]
  fit_alpha2 <- fit_jags[, "alpha[2,1,1]"][[1]]
  fit_theta <- fit_jags[, "theta[1,1,1]"][[1]]
  fit_psi11 <- fit_jags[, "psi[1,1,1]"][[1]]
  fit_psi22 <- fit_jags[, "psi[2,2,1]"][[1]]
  fit_psi21 <- fit_jags[, "psi[1,2,1]"][[1]]
  
  fit_alpha1 <- fit_alpha1[(1:L) * thin]
  fit_alpha2 <- fit_alpha2[(1:L) * thin]
  fit_theta <- fit_theta[(1:L) * thin]
  fit_psi11 <- fit_psi11[(1:L) * thin]
  fit_psi22 <- fit_psi22[(1:L) * thin]
  fit_psi21 <- fit_psi21[(1:L) * thin]
  
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
  
  
  tempMatrix[i, 13] <- counter
}

# tempMatrix

rank_alpha1 <- tempMatrix[, 1]
rank_alpha2 <- tempMatrix[, 2]
rank_theta <- tempMatrix[, 3]
rank_psi11 <- tempMatrix[, 4]
rank_psi22 <- tempMatrix[, 5]
rank_psi21 <- tempMatrix[, 6]

# ttt <- which(tempMatrix[,13] == 1000)
# 
# sbc_alpha1 <- sbc_rank(rank_alpha1[-ttt], reps = reps - length(ttt), L = L, title = "Alpha 1")
# sbc_alpha2 <- sbc_rank(rank_alpha2[-ttt], reps = reps - length(ttt), L = L, title = "Alpha 2")
# sbc_theta <- sbc_rank(rank_theta[-ttt], reps = reps - length(ttt), L = L, title = "Theta")
# sbc_psi11 <- sbc_rank(rank_psi11[-ttt], reps = reps - length(ttt), L = L, title = "Psi 11")
# sbc_psi22 <- sbc_rank(rank_psi22[-ttt], reps = reps - length(ttt), L = L, title = "Psi 22")
# sbc_psi21 <- sbc_rank(rank_psi21[-ttt], reps = reps - length(ttt), L = L, title = "Psi 21")

sbc_alpha1 <- sbc_rank(rank_alpha1, reps = reps, L = L, title = "Alpha 1")
sbc_alpha2 <- sbc_rank(rank_alpha2, reps = reps, L = L, title = "Alpha 2")
sbc_theta <- sbc_rank(rank_theta, reps = reps, L = L, title = "Theta")
sbc_psi11 <- sbc_rank(rank_psi11, reps = reps, L = L, title = "Psi 11")
sbc_psi22 <- sbc_rank(rank_psi22, reps = reps, L = L, title = "Psi 22")
sbc_psi21 <- sbc_rank(rank_psi21, reps = reps, L = L, title = "Psi 21")


# tt <- tempMatrix[!which(tempMatrix[,13] == 500), ]

bias_alpha1 <- data.frame(tempMatrix) %>%
  ggplot() + geom_histogram(aes(x = X7)) + theme_classic() + xlab(expression(paste("bias ", alpha[1])))

bias_alpha2 <- data.frame(tempMatrix) %>%
  ggplot() + geom_histogram(aes(x = X8)) + theme_classic() + xlab(expression(paste("bias ", alpha[2])))

bias_theta <- data.frame(tempMatrix) %>%
  ggplot() + geom_histogram(aes(x = X9)) + theme_classic() + xlab(expression(paste("bias ", theta)))

bias_psi11 <- data.frame(tempMatrix) %>%
  ggplot() + geom_histogram(aes(x = X10)) + theme_classic() + xlab(expression(paste("bias ", psi[11])))

bias_psi22 <- data.frame(tempMatrix) %>%
  ggplot() + geom_histogram(aes(x = X11)) + theme_classic() + xlab(expression(paste("bias ", psi[22])))

bias_psi21 <- data.frame(tempMatrix) %>%
  ggplot() + geom_histogram(aes(x = X12)) + theme_classic() + xlab(expression(paste("bias ", psi[21])))



gridExtra::grid.arrange(sbc_alpha1, sbc_alpha2, sbc_theta,
                        sbc_psi11, sbc_psi22, sbc_psi21, ncol = 3)

gridExtra::grid.arrange(bias_alpha1, bias_alpha2, bias_theta,
                        bias_psi11, bias_psi22, bias_psi21, ncol = 3)




pdf(paste0("results/jags_plots_wish_", n_ind, "ind_", n_time, "time_", reps, "reps", ".pdf"))
gridExtra::grid.arrange(sbc_alpha1, sbc_alpha2, sbc_theta,
                        sbc_psi11, sbc_psi22, sbc_psi21, ncol = 3)

gridExtra::grid.arrange(bias_alpha1, bias_alpha2, bias_theta,
                        bias_psi11, bias_psi22, bias_psi21, ncol = 3)

dev.off()


saveRDS(tempMatrix, 
        paste0("results/jags_finalMatrix_wish_", n_ind, "ind_", n_time, "time_", reps, "reps", ".rds"))


