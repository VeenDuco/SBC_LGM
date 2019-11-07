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

estimate <- stan_model("stan_models/LGM_estimate_inv_wish.stan")

Ncores <- detectCores(logical = T)
Ncores

if(Sys.info()['sysname'] != "Windows"){
  require("doMC")
  registerDoMC(Ncores - 1)
}else{
  require("doParallel")
  cl <- makeCluster(Ncores - 1)
  registerDoParallel(cl)
  #snow is also an option
}



n_ind <- 10
n_time <- 3
# informative <- TRUE
reps <- 511
bins <- 8
L <- bins - 1

set.seed(1236)
PSI <- array(NA, dim = c(2, 2, reps))
for(i in 1:reps){
  PSI[, , i] <- MCMCpack::riwish(3, diag(2))
  # print(is.positive.definite(round(PSI[ , , i], 10)))
}
# for(i in 1:reps){
#   PSI <- rWishart(reps, 3, diag(2))
#   # PSI[ , , i] <- matrix(c(sigma1[i], cov_sig[i], cov_sig[i], sigma2[i]), ncol = 2)
#   # print(is.positive.definite(PSI[ , , i]))
# }
# 
alpha1 <- rnorm(reps, 35, 10)
alpha2 <- rnorm(reps, -5, 5)
theta <- list()
theta[[1]] <- rlnorm(reps, 0, 1)
theta[[2]] <- rlnorm(reps, 0, 1)
theta[[3]] <- rlnorm(reps, 0, 1)
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
                             sd = theta[[j]][k])
    }
  }
}
# 


# ysim

prior_y <- ysim

prior_alpha1 <- alpha1
prior_alpha2 <- alpha2
prior_theta11 <- theta[[1]]
prior_theta22 <- theta[[2]]
prior_theta33 <- theta[[3]]
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
                        .combine = rbind,
                        .verbose = TRUE) %dopar% {
                          
                          if(!exists("pb")) pb <- tkProgressBar("Parallel task", min = 1, max = reps)
                          setTkProgressBar(pb, i)  
                          
                          thin <- 1
                          draws <- L
                          warm <- warm.init
                          counter <- 1
                          fit <- sampling(estimate, chains = 1, warmup = warm, iter = warm + draws, 
                                          seed = i,
                                          data = list(n_ind = n_ind,
                                                      n_time = n_time,
                                                      y = prior_y[, , i],
                                                      x = c(.083, .25, 1)))
                          
                          
                          ## note 2 * L, was L but still signs of autocorrelation.
                          while(min(na.omit(summary(fit)$summary[, "n_eff"])) < 2 * L | max(na.omit(summary(fit)$summary[, "Rhat"])) > 1.02) { 
                            counter <- counter + 1
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
                          fit_theta11 <- rstan::extract(fit, pars = "theta[1]")[[1]]
                          fit_theta22 <- rstan::extract(fit, pars = "theta[2]")[[1]]
                          fit_theta33 <- rstan::extract(fit, pars = "theta[3]")[[1]]
                          fit_psi11 <- rstan::extract(fit, pars = "Psi[1,1]")[[1]]
                          fit_psi22 <- rstan::extract(fit, pars = "Psi[2,2]")[[1]]
                          fit_psi21 <- rstan::extract(fit, pars = "Psi[2,1]")[[1]]
                          
                          
                          fit_alpha1 <- fit_alpha1[(1:L) * thin]
                          fit_alpha2 <- fit_alpha2[(1:L) * thin]
                          fit_theta11 <- fit_theta11[(1:L) * thin]
                          fit_theta22 <- fit_theta22[(1:L) * thin]
                          fit_theta33 <- fit_theta33[(1:L) * thin]
                          fit_psi11 <- fit_psi11[(1:L) * thin]
                          fit_psi22 <- fit_psi22[(1:L) * thin]
                          fit_psi21 <- fit_psi21[(1:L) * thin]
                          
                          ## bias
                          alpha1_bias <- summary(fit, pars = "alpha1")$summary[1] - prior_alpha1[i]
                          alpha2_bias <- summary(fit, pars = "alpha2")$summary[1] - prior_alpha2[i]
                          theta11_bias <- summary(fit, pars = "theta[1]")$summary[1] - prior_theta11[i]
                          theta22_bias <- summary(fit, pars = "theta[2]")$summary[1] - prior_theta22[i]
                          theta33_bias <- summary(fit, pars = "theta[3]")$summary[1] - prior_theta33[i]
                          psi11_bias <- summary(fit, pars = "Psi[1,1]")$summary[1] - prior_psi11[i]
                          psi22_bias <- summary(fit, pars = "Psi[2,2]")$summary[1] - prior_psi22[i]
                          psi21_bias <- summary(fit, pars = "Psi[2,1]")$summary[1] - prior_psi21[i]
                          
                          tempMatrix <- matrix(NA, ncol = 17, nrow = 1)
                          
                          tempMatrix[1, 1] <- sum(fit_alpha1 < prior_alpha1[i])
                          tempMatrix[1, 2] <- sum(fit_alpha2 < prior_alpha2[i])
                          tempMatrix[1, 3] <- sum(fit_theta11 < prior_theta11[i])
                          tempMatrix[1, 4] <- sum(fit_theta22 < prior_theta22[i])
                          tempMatrix[1, 5] <- sum(fit_theta33 < prior_theta33[i])
                          tempMatrix[1, 6] <- sum(fit_psi11 < prior_psi11[i])
                          tempMatrix[1, 7] <- sum(fit_psi22 < prior_psi22[i])
                          tempMatrix[1, 8] <- sum(fit_psi21 < prior_psi21[i])
                          
                          tempMatrix[1, 9] <- alpha1_bias
                          tempMatrix[1, 10] <- alpha2_bias
                          tempMatrix[1, 11] <- theta11_bias
                          tempMatrix[1, 12] <- theta22_bias
                          tempMatrix[1, 13] <- theta33_bias
                          tempMatrix[1, 14] <- psi11_bias
                          tempMatrix[1, 15] <- psi22_bias
                          tempMatrix[1, 16] <- psi21_bias
                          
                          tempMatrix[1, 17] <- counter
                          
                          tempMatrix
                          
                        }


rank_alpha1 <- finalMatrix[, 1]
rank_alpha2 <- finalMatrix[, 2]
rank_theta11 <- finalMatrix[, 3]
rank_theta22 <- finalMatrix[, 4]
rank_theta33 <- finalMatrix[, 5]
rank_psi11 <- finalMatrix[, 6]
rank_psi22 <- finalMatrix[, 7]
rank_psi21 <- finalMatrix[, 8]

sbc_alpha1 <- sbc_rank(rank_alpha1, reps = reps, L = L, title = "Alpha 1")
sbc_alpha2 <- sbc_rank(rank_alpha2, reps = reps, L = L, title = "Alpha 2")
sbc_theta11 <- sbc_rank(rank_theta11, reps = reps, L = L, title = "Theta 11")
sbc_theta22 <- sbc_rank(rank_theta22, reps = reps, L = L, title = "Theta 22")
sbc_theta33 <- sbc_rank(rank_theta33, reps = reps, L = L, title = "Theta 33")
sbc_psi11 <- sbc_rank(rank_psi11, reps = reps, L = L, title = "Psi 11")
sbc_psi22 <- sbc_rank(rank_psi22, reps = reps, L = L, title = "Psi 22")
sbc_psi21 <- sbc_rank(rank_psi21, reps = reps, L = L, title = "Psi 21")


bias_alpha1 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X9)) + theme_classic() + xlab(expression(paste("bias ", alpha[1])))

bias_alpha2 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X10)) + theme_classic() + xlab(expression(paste("bias ", alpha[2])))

bias_theta11 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X11)) + theme_classic() + xlab(expression(paste("bias ", theta[11])))

bias_theta22 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X12)) + theme_classic() + xlab(expression(paste("bias ", theta[22])))

bias_theta33 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X13)) + theme_classic() + xlab(expression(paste("bias ", theta[33])))

bias_psi11 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X14)) + theme_classic() + xlab(expression(paste("bias ", psi[11])))

bias_psi22 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X15)) + theme_classic() + xlab(expression(paste("bias ", psi[22])))

bias_psi21 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X16)) + theme_classic() + xlab(expression(paste("bias ", psi[21])))



gridExtra::grid.arrange(sbc_alpha1, sbc_alpha2, sbc_theta11,
                        sbc_theta22, sbc_theta33,
                        sbc_psi11, sbc_psi22, sbc_psi21, ncol = 3)

gridExtra::grid.arrange(bias_alpha1, bias_alpha2, bias_theta11,
                        bias_theta22, bias_theta33,
                        bias_psi11, bias_psi22, bias_psi21, ncol = 3)




pdf(paste0("results/plots_inv_wish_", n_ind, "ind_", n_time, "time_", reps, "reps", ".pdf"))
gridExtra::grid.arrange(sbc_alpha1, sbc_alpha2, sbc_theta11,
                        sbc_theta22, sbc_theta33,
                        sbc_psi11, sbc_psi22, sbc_psi21, ncol = 3)

gridExtra::grid.arrange(bias_alpha1, bias_alpha2, bias_theta11,
                        bias_theta22, bias_theta33,
                        bias_psi11, bias_psi22, bias_psi21, ncol = 3)

dev.off()

saveRDS(finalMatrix, 
        paste0("results/finalMatrix_inv_wish_", n_ind, "ind_", n_time, "time_", reps, "reps", ".rds"))




######### JAGS
library(blavaan)
library(rjags)

model <- '
I =~ 1 * x1 + 1 * x2 + 1 * x3
S =~ .083 * x1 + .25 * x2 + 1 * x3
I ~ prior("dnorm(35, 0.01)") * 1
S ~ prior("dnorm(-5, 0.04)") * 1

x1 ~~  x1
x2 ~~  x2
x3 ~~ x3

x1 ~ 0
x2 ~ 0
x3 ~ 0

'

dpri <- dpriors(itheta = "dlnorm(0, 1)[sd]",
                target = "jags")

set.seed(123)
warm.init <- 1000

finalMatrix  <- foreach(i = 1:reps, .packages = c("blavaan", "rjags", "tcltk", "dplyr"), 
                        .combine = rbind,
                        .verbose = TRUE) %dopar% {    
                          
                          thin <- 1
                          draws <- L
                          warm <- warm.init
                          counter <- 1
                          
                          if(!exists("pb")) pb <- tkProgressBar("Parallel task", min = 1, max = reps)
                          setTkProgressBar(pb, i)
                          
                          fit_blav <- bcfa(model, data = data.frame(x1 = prior_y[ , 1, i],
                                                                    x2 = prior_y[ , 2, i],
                                                                    x3 = prior_y[ , 3, i]),
                                           dp = dpri, cp = "fa", save.lvs = TRUE, seed = i,
                                           burnin = warm, adapt = warm, sample = draws, n.chains = 1,
                                           target = "jags", mcmcfile = F, test = "none")
                          
                          
                          rhat <- max(data.frame(blavInspect(fit_blav, what = "mcmc")[[1]]) %>%
                                        apply(., 2, rstan::Rhat))
                          
                          neff <- min(blavInspect(fit_blav, what = "neff"))
                          
                          
                          while(neff < 2 * L | rhat > 1.02) {
                            counter <- counter + 1
                            thin <- thin * 2
                            draws <- draws * 2
                            warm <- round(warm * 1.1)
                            
                            fit_blav <- bcfa(model, data = data.frame(x1 = prior_y[ , 1, i],
                                                                      x2 = prior_y[ , 2, i],
                                                                      x3 = prior_y[ , 3, i]),
                                             dp = dpri, cp = "fa", save.lvs = TRUE,
                                             burnin = warm, adapt = warm, sample = draws, n.chains = 1,
                                             target = "jags", mcmcfile = F, test = "none")
                            
                            
                            rhat <- max(data.frame(blavInspect(fit_blav, what = "mcmc")[[1]]) %>%
                                          apply(., 2, rstan::Rhat))
                            
                            neff <- min(blavInspect(fit_blav, what = "neff"))
                          }
                          
                          
                          # blavInspect(fit_blav, what = "mcmc")[, "alpha[1,1,1]"]
                          
                          fit_alpha1 <- blavInspect(fit_blav, what = "mcmc")[, "alpha[1,1,1]"][[1]]
                          fit_alpha2 <- blavInspect(fit_blav, what = "mcmc")[, "alpha[2,1,1]"][[1]]
                          fit_theta11 <- blavInspect(fit_blav, what = "mcmc")[, "theta[1,1,1]"][[1]]
                          fit_theta22 <- blavInspect(fit_blav, what = "mcmc")[, "theta[2,2,1]"][[1]]
                          fit_theta33 <- blavInspect(fit_blav, what = "mcmc")[, "theta[3,3,1]"][[1]]
                          fit_psi11 <- blavInspect(fit_blav, what = "mcmc")[, "psi[1,1,1]"][[1]]
                          fit_psi22 <- blavInspect(fit_blav, what = "mcmc")[, "psi[2,2,1]"][[1]]
                          fit_psi21 <- blavInspect(fit_blav, what = "mcmc")[, "psi[1,2,1]"][[1]]
                          
                          fit_alpha1 <- fit_alpha1[(1:L) * thin] %>% unname()
                          fit_alpha2 <- fit_alpha2[(1:L) * thin] %>% unname()
                          fit_theta11 <- fit_theta11[(1:L) * thin] %>% unname()
                          fit_theta11 <- sqrt(fit_theta11)
                          fit_theta22 <- fit_theta22[(1:L) * thin] %>% unname()
                          fit_theta22 <- sqrt(fit_theta22)
                          fit_theta33 <- fit_theta33[(1:L) * thin] %>% unname()
                          fit_theta33 <- sqrt(fit_theta33)
                          fit_psi11 <- fit_psi11[(1:L) * thin] %>% unname()
                          fit_psi22 <- fit_psi22[(1:L) * thin] %>% unname()
                          fit_psi21 <- fit_psi21[(1:L) * thin] %>% unname()
                          
                          
                          alpha1_bias <- unname(blavInspect(fit_blav, what = "postmean")[1]) - prior_alpha1[i]
                          alpha2_bias <- unname(blavInspect(fit_blav, what = "postmean")[2]) - prior_alpha2[i]
                          theta11_bias <- mean(fit_theta11) - prior_theta11[i]
                          theta22_bias <- mean(fit_theta22) - prior_theta22[i]
                          theta33_bias <- mean(fit_theta33) - prior_theta33[i]
                          psi11_bias <- unname(blavInspect(fit_blav, what = "postmean")[6]) - prior_psi11[i]
                          psi22_bias <- unname(blavInspect(fit_blav, what = "postmean")[7]) - prior_psi22[i]
                          psi21_bias <- unname(blavInspect(fit_blav, what = "postmean")[8]) - prior_psi21[i]
                          
                          tempMatrix <- matrix(NA, ncol = 17, nrow = 1)
                          
                          tempMatrix[, 1] <- sum(fit_alpha1 < prior_alpha1[i])
                          tempMatrix[, 2] <- sum(fit_alpha2 < prior_alpha2[i])
                          tempMatrix[, 3] <- sum(fit_theta11 < prior_theta11[i])
                          tempMatrix[, 4] <- sum(fit_theta22 < prior_theta22[i])
                          tempMatrix[, 5] <- sum(fit_theta33 < prior_theta33[i])
                          tempMatrix[, 6] <- sum(fit_psi11 < prior_psi11[i])
                          tempMatrix[, 7] <- sum(fit_psi22 < prior_psi22[i])
                          tempMatrix[, 8] <- sum(fit_psi21 < prior_psi21[i])
                          
                          tempMatrix[, 9] <- alpha1_bias
                          tempMatrix[, 10] <- alpha2_bias
                          tempMatrix[, 11] <- theta11_bias
                          tempMatrix[, 12] <- theta22_bias
                          tempMatrix[, 13] <- theta33_bias
                          tempMatrix[, 14] <- psi11_bias
                          tempMatrix[, 15] <- psi22_bias
                          tempMatrix[, 16] <- psi21_bias
                          
                          tempMatrix[, 17] <- counter
                          
                          tempMatrix
                          
                        }

finalMatrix <- tempMatrix

rank_alpha1 <- finalMatrix[, 1]
rank_alpha2 <- finalMatrix[, 2]
rank_theta11 <- finalMatrix[, 3]
rank_theta22 <- finalMatrix[, 4]
rank_theta33 <- finalMatrix[, 5]
rank_psi11 <- finalMatrix[, 6]
rank_psi22 <- finalMatrix[, 7]
rank_psi21 <- finalMatrix[, 8]

sbc_alpha1 <- sbc_rank(rank_alpha1, reps = reps, L = L, title = "Alpha 1")
sbc_alpha2 <- sbc_rank(rank_alpha2, reps = reps, L = L, title = "Alpha 2")
sbc_theta11 <- sbc_rank(rank_theta11, reps = reps, L = L, title = "Theta 11")
sbc_theta22 <- sbc_rank(rank_theta22, reps = reps, L = L, title = "Theta 22")
sbc_theta33 <- sbc_rank(rank_theta33, reps = reps, L = L, title = "Theta 33")
sbc_psi11 <- sbc_rank(rank_psi11, reps = reps, L = L, title = "Psi 11")
sbc_psi22 <- sbc_rank(rank_psi22, reps = reps, L = L, title = "Psi 22")
sbc_psi21 <- sbc_rank(rank_psi21, reps = reps, L = L, title = "Psi 21")


bias_alpha1 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X9)) + theme_classic() + xlab(expression(paste("bias ", alpha[1])))

bias_alpha2 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X10)) + theme_classic() + xlab(expression(paste("bias ", alpha[2])))

bias_theta11 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X11)) + theme_classic() + xlab(expression(paste("bias ", theta[11])))

bias_theta22 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X12)) + theme_classic() + xlab(expression(paste("bias ", theta[22])))

bias_theta33 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X13)) + theme_classic() + xlab(expression(paste("bias ", theta[33])))

bias_psi11 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X14)) + theme_classic() + xlab(expression(paste("bias ", psi[11])))

bias_psi22 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X15)) + theme_classic() + xlab(expression(paste("bias ", psi[22])))

bias_psi21 <- data.frame(finalMatrix) %>%
  ggplot() + geom_histogram(aes(x = X16)) + theme_classic() + xlab(expression(paste("bias ", psi[21])))



gridExtra::grid.arrange(sbc_alpha1, sbc_alpha2, sbc_theta11,
                        sbc_theta22, sbc_theta33,
                        sbc_psi11, sbc_psi22, sbc_psi21, ncol = 3)

gridExtra::grid.arrange(bias_alpha1, bias_alpha2, bias_theta11,
                        bias_theta22, bias_theta33,
                        bias_psi11, bias_psi22, bias_psi21, ncol = 3)




pdf(paste0("results/jags_plots_inv_wish_", n_ind, "ind_", n_time, "time_", reps, "reps", ".pdf"))
gridExtra::grid.arrange(sbc_alpha1, sbc_alpha2, sbc_theta11,
                        sbc_theta22, sbc_theta33,
                        sbc_psi11, sbc_psi22, sbc_psi21, ncol = 3)

gridExtra::grid.arrange(bias_alpha1, bias_alpha2, bias_theta11,
                        bias_theta22, bias_theta33,
                        bias_psi11, bias_psi22, bias_psi21, ncol = 3)

dev.off()

saveRDS(finalMatrix, 
        paste0("results/jags_finalMatrix_inv_wish_", n_ind, "ind_", n_time, "time_", reps, "reps", ".rds"))




