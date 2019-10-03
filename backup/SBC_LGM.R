# LGM SBC:
library(rstan)
library(foreach)
library(parallel)
library(doParallel)
library(bayesplot)
library(tcltk)
library(dplyr)

# detectCores()
no_cores <- detectCores() - 1
cl <- makeCluster(no_cores)
registerDoParallel(cl)

informative <- TRUE
reps <- 20
bins <- 5
L <- bins - 1

if(informative == FALSE){
  # generate <- stan_model("stan_models/LGM_generate_uninformative.stan")
  # saveRDS(generate ,"stan_models/LGM_generate_uninformative.rds")
  generate <- readRDS("stan_models/LGM_generate_uninformative.rds")
  # estimate <- stan_model("stan_models/LGM_estimate_uninformative.stan")
  # saveRDS(estimate ,"stan_models/LGM_estimate_uninformative.rds")
  estimate <- readRDS("stan_models/LGM_estimate_uninformative.rds")
} else {
  # generate <- stan_model("stan_models/LGM_generate.stan")
  # saveRDS(generate ,"stan_models/LGM_generate.rds")
  generate <- readRDS("stan_models/LGM_generate.rds")
  # estimate <- stan_model("stan_models/LGM_estimate.stan")
  # saveRDS(estimate ,"stan_models/LGM_estimate.rds")
  estimate <- readRDS("stan_models/LGM_estimate.rds")
}

n_ind <- 30
n_time <- 3



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

set.seed(123)
warm.init <- 1000




prior_sample <- sampling(generate, data = list(n_ind = n_ind,
                                               n_time = n_time),
                         algorithm = "Fixed_param",
                         seed = 20092019, chains = 1,
                         warmup = 0, iter = reps)


prior_y <- rstan::extract(prior_sample, pars = c("y_rep"))[[1]]
prior_alpha1 <- rstan::extract(prior_sample, pars = c("alpha1"))[[1]]
prior_alpha2 <- rstan::extract(prior_sample, pars = c("alpha2"))[[1]]
prior_theta <- rstan::extract(prior_sample, pars = c("theta"))[[1]]
prior_psi11 <- rstan::extract(prior_sample, pars = c("Psi[1,1]"))[[1]]
prior_psi22 <- rstan::extract(prior_sample, pars = c("Psi[2,2]"))[[1]]
prior_psi21 <- rstan::extract(prior_sample, pars = c("Psi[2,1]"))[[1]]

# prior_psi <- rstan::extract(prior_sample, pars = c("Psi"))[[1]]
# prior_psi11 <- prior_psi[, 1, 1]
# prior_psi22 <- prior_psi[, 2, 2]
# prior_psi21 <- prior_psi[, 2, 1]

# prior_y[1, ,]









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
                                                      y = prior_y[i, , ],
                                                      x = c(.083, .25, 1)))
                          
                          
                          ## note 2 * L, was L but still signs of autocorrelation.
                          while(min(na.omit(summary(fit)$summary[, "n_eff"])) < 2 * L | max(na.omit(summary(fit)$summary[, "Rhat"])) > 1.02) { 
                            thin <- thin * 2
                            draws <- draws * 2
                            warm <- warm * 1.1
                            fit <- sampling(estimate, chains = 1, warmup = warm, iter = warm + draws, 
                                            data = list(n_ind = n_ind,
                                                        n_time = n_time,
                                                        y = prior_y[i, , ],
                                                        x = c(.083, .25, 1)))
                          }
                          
                          saveRDS(fit, paste0("temp/fit_", i ,".rds"))
                          
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


pdf(paste0("results/plots_", n_ind, "ind_", n_time, "time_", reps, "reps", ".pdf"))
gridExtra::grid.arrange(sbc_alpha1, sbc_alpha2, sbc_theta,
                        sbc_psi11, sbc_psi22, sbc_psi21, ncol = 3)

gridExtra::grid.arrange(bias_alpha1, bias_alpha2, bias_theta,
                        bias_psi11, bias_psi22, bias_psi21, ncol = 3)

dev.off()

saveRDS(finalMatrix, 
        paste0("results/finalMatrix_", n_ind, "ind_", n_time, "time_", reps, "reps_",
               "informative_priors_", informative, ".rds"))
