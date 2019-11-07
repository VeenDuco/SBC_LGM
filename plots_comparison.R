library(ggplot2)
library(gridExtra)
library(ggpubr)
library(dplyr)
fm_jags <- readRDS("results/jags_finalMatrix_inv_wish_10ind_3time_511reps.rds")
fm_stan <- readRDS("results/finalMatrix_inv_wish_10ind_3time_511reps.rds")

n_ind <- 10
n_time <- 3
reps <- 511
bins <- 8
L <- bins - 1

set.seed(1236)
PSI <- array(NA, dim = c(2, 2, reps))
for(i in 1:reps){
  PSI[, , i] <- MCMCpack::riwish(3, diag(2))
}
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
for(i in 1:n_ind){
  for(j in 1:n_time){
    for(k in 1:reps){
      ysim[i, j, k] <- rnorm(1, mean = I(1 * tau[i, 1, k] + time[j] * tau[i, 2, k]),
                             sd = theta[[j]][k])
    }
  }
}
prior_y <- ysim
prior_alpha1 <- alpha1
prior_alpha2 <- alpha2
prior_theta11 <- theta[[1]]
prior_theta22 <- theta[[2]]
prior_theta33 <- theta[[3]]
prior_psi11 <- PSI[1, 1, ]
prior_psi22 <- PSI[2, 2, ]
prior_psi21 <- PSI[2, 1, ]


# plot(prior_alpha1, (fm_jags[,9] + prior_alpha1), col = "red", 
#      xlim = c(0, 70), ylim = c(0, 70), 
#      xlab = expression(paste("simulated ", alpha[1])),
#      ylab = "Estimated mean", bty = "n")
# points(prior_alpha1, (fm_stan[,9] + prior_alpha1), col = "blue", pch = 3)

kleur <- c("jags" = c("red"),
           "stan" = c("blue"),
           "J" = 2)
# vorm <- c("jags" = 3,
#           "stan" = 2)

df <- data.frame(a1 = prior_alpha1,
           a2 = prior_alpha2,
           th1 = prior_theta11,
           th2 = prior_theta22,
           th3 = prior_theta33,
           p11 = prior_psi11,
           p22 = prior_psi22,
           p21 = prior_psi21,
           fm_jags[, 9:16], fm_stan[, 9:16]) 

# df_prior <- data.frame(a1 = prior_alpha1,
#                        a2 = prior_alpha2,
#                        th1 = prior_theta11,
#                        th2 = prior_theta22,
#                        th3 = prior_theta33,
#                        p11 = prior_psi11,
#                        p22 = prior_psi22,
#                        p21 = prior_psi21)
# 
# df_jags <- data.frame(fm_jags[, 9:16])
# df_stan <- data.frame(fm_stan[, 9:16])
# 
# bind_rows(df_prior, df_jags, df_stan)

p_a1 <- df %>%
  ggplot() +
  geom_point(aes(x = a1, y = X1 + a1, color = "jags"), 
             size = 3, alpha = .8) + 
  geom_point(aes(x = a1, y = X1.1 + a1, color = "stan"), 
             size = 3, alpha = .5) +
  geom_abline(slope = 1, lty = 2) +
  theme_classic() + xlab(expression(paste("Simulated ", alpha[1]))) +
  ylab("Estimated mean") + scale_color_manual(values = kleur, name = "")

p_a2 <- df %>%
  ggplot() +
  geom_point(aes(x = a2, y = X2 + a2, color = "jags"), 
             size = 3, alpha = .8) + 
  geom_point(aes(x = a2, y = X2.1 + a2, color = "stan"), 
             size = 3, alpha = .5) +
  geom_abline(slope = 1, lty = 2) +
  theme_classic() + xlab(expression(paste("Simulated ", alpha[2]))) +
  ylab("Estimated mean") + scale_color_manual(values = kleur, name = "")


p_th1 <- df %>%
  ggplot() +
  geom_point(aes(x = th1, y = X3 + th1, color = "jags"), 
             size = 3, alpha = .8) + 
  geom_point(aes(x = th1, y = X3.1 + th1, color = "stan"), 
             size = 3, alpha = .5) +
  geom_abline(slope = 1, lty = 2) +
  theme_classic() + xlab(expression(paste("Simulated ", theta[11]))) +
  ylab("Estimated mean") + scale_color_manual(values = kleur, name = "")

p_th2 <- df %>%
  ggplot() +
  geom_point(aes(x = th2, y = X4 + th2, color = "jags"), 
             size = 3, alpha = .8) + 
  geom_point(aes(x = th2, y = X4.1 + th2, color = "stan"), 
             size = 3, alpha = .5) +
  geom_abline(slope = 1, lty = 2) +
  theme_classic() + xlab(expression(paste("Simulated ", theta[22]))) +
  ylab("Estimated mean") + scale_color_manual(values = kleur, name = "")

p_th3 <- df %>%
  ggplot() +
  geom_point(aes(x = th3, y = X5 + th3, color = "jags"), 
             size = 3, alpha = .8) + 
  geom_point(aes(x = th3, y = X5.1 + th3, color = "stan"), 
             size = 3, alpha = .5) +
  geom_abline(slope = 1, lty = 2) +
  theme_classic() + xlab(expression(paste("Simulated ", theta[33]))) +
  ylab("Estimated mean") + scale_color_manual(values = kleur, name = "")


p_p11 <- df %>%
  ggplot() +
  geom_point(aes(x = p11, y = X6 + p11, color = "jags"), 
             size = 3, alpha = .8) + 
  geom_point(aes(x = p11, y = X6.1 + p11, color = "stan"), 
             size = 3, alpha = .5) +
  geom_abline(slope = 1, lty = 2) +
  theme_classic() + xlab(expression(paste("Simulated ", psi[11]))) +
  ylab("Estimated mean") + scale_color_manual(values = kleur, name = "")

p_p22 <- df %>%
  ggplot() +
  geom_point(aes(x = p22, y = X7 + p22, color = "jags"), 
             size = 3, alpha = .8) + 
  geom_point(aes(x = p22, y = X7.1 + p22, color = "stan"), 
             size = 3, alpha = .5) +
  geom_abline(slope = 1, lty = 2) +
  theme_classic() + xlab(expression(paste("Simulated ", psi[22]))) +
  ylab("Estimated mean") + scale_color_manual(values = kleur, name = "")

p_p21 <- df %>%
  ggplot() +
  geom_point(aes(x = p21, y = X8 + p21, color = "jags"), 
             size = 3, alpha = .8) + 
  geom_point(aes(x = p21, y = X8.1 + p21, color = "stan"), 
             size = 3, alpha = .5) +
  geom_abline(slope = 1, lty = 2) +
  theme_classic() + xlab(expression(paste("Simulated ", psi[21]))) +
  ylab("Estimated mean") + scale_color_manual(values = kleur, name = "")

ggarrange(p_a1, p_a2, p_th1, p_th2, p_th3,
          p_p11, p_p22, p_p21,
          ncol = 2, nrow =4,
          common.legend = T, legend = "top")


#   scale_size_manual(values=c(3, 3)) +
#   theme(legend.position="top")
