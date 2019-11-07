library(MASS)
library(rstan)
library(foreach)
library(parallel)
library(doParallel)
library(bayesplot)
library(tcltk)
library(dplyr)
library(brms)
library(matrixcalc)
lasso_reg_pred <- stan_model("regression/stan_models/lasso_reg_pred.stan")
# lasso_reg <- stan_model("regression/stan_models/lasso_reg.stan")
simple_reg_model_pred <- stan_model("regression/stan_models/simple_reg_pred.stan")
# simple_reg_model <- stan_model("regression/stan_models/simple_reg.stan")

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

n_ind <- 4
n_test <- 30
reps <- 1023
bins <- 16
# reps <- 1023
# bins <- 32

source("regression/data_gen.r")

warm.init <- 1000

source("regression/SBC_lasso_pred.r")
# source("regression/SBC_lasso.r")
source("regression/SBC_reg_pred.r")
# source("regression/SBC_reg.r")
finalMatrix_correct <- finalMatrix
simple_reg_model_pred <- stan_model("regression/stan_models/simple_reg_pred_wrong.stan")
source("regression/SBC_reg_pred.r")
finalMatrix_wrong <- finalMatrix

# correct reg pcitures
rank_b0_correct <- finalMatrix_correct[, 1]
rank_b1_correct <- finalMatrix_correct[, 2]
rank_b2_correct <- finalMatrix_correct[, 3]
rank_b3_correct <- finalMatrix_correct[, 4]
rank_b4_correct <- finalMatrix_correct[, 5]
rank_b5_correct <- finalMatrix_correct[, 6]
rank_sigma_correct <- finalMatrix_correct[, 7]

sbc_b0_correct <- sbc_rank(rank_b0_correct, reps = reps, L = L, title = "b0")
sbc_b1_correct <- sbc_rank(rank_b1_correct, reps = reps, L = L, title = "b1" )
sbc_b2_correct <- sbc_rank(rank_b2_correct, reps = reps, L = L, title = "b2")
sbc_b3_correct <- sbc_rank(rank_b3_correct, reps = reps, L = L, title = "b3")
sbc_b4_correct <- sbc_rank(rank_b4_correct, reps = reps, L = L, title = "b4")
sbc_b5_correct <- sbc_rank(rank_b5_correct, reps = reps, L = L, title = "b5")
sbc_sigma_correct <- sbc_rank(rank_sigma_correct, reps = reps, L = L, title = "sigma")


# wrong reg pcitures

rank_b0_wrong <- finalMatrix_wrong[, 1]
rank_b1_wrong <- finalMatrix_wrong[, 2]
rank_b2_wrong <- finalMatrix_wrong[, 3]
rank_b3_wrong <- finalMatrix_wrong[, 4]
rank_b4_wrong <- finalMatrix_wrong[, 5]
rank_b5_wrong <- finalMatrix_wrong[, 6]
rank_sigma_wrong <- finalMatrix_wrong[, 7]

sbc_b0_wrong <- sbc_rank(rank_b0_wrong, reps = reps, L = L, title = "b0")
sbc_b1_wrong <- sbc_rank(rank_b1_wrong, reps = reps, L = L, title = "b1" )
sbc_b2_wrong <- sbc_rank(rank_b2_wrong, reps = reps, L = L, title = "b2")
sbc_b3_wrong <- sbc_rank(rank_b3_wrong, reps = reps, L = L, title = "b3")
sbc_b4_wrong <- sbc_rank(rank_b4_wrong, reps = reps, L = L, title = "b4")
sbc_b5_wrong <- sbc_rank(rank_b5_wrong, reps = reps, L = L, title = "b5")
sbc_sigma_wrong <- sbc_rank(rank_sigma_wrong, reps = reps, L = L, title = "sigma")

# lasso pcitures

rank_b0_lasso <- finalMatrix_lasso[, 1]
rank_b1_lasso <- finalMatrix_lasso[, 2]
rank_b2_lasso <- finalMatrix_lasso[, 3]
rank_b3_lasso <- finalMatrix_lasso[, 4]
rank_b4_lasso <- finalMatrix_lasso[, 5]
rank_b5_lasso <- finalMatrix_lasso[, 6]
rank_sigma_lasso <- finalMatrix_lasso[, 7]

sbc_b0_lasso <- sbc_rank(rank_b0_lasso, reps = reps, L = L, title = "b0")
sbc_b1_lasso <- sbc_rank(rank_b1_lasso, reps = reps, L = L, title = "b1" )
sbc_b2_lasso <- sbc_rank(rank_b2_lasso, reps = reps, L = L, title = "b2")
sbc_b3_lasso <- sbc_rank(rank_b3_lasso, reps = reps, L = L, title = "b3")
sbc_b4_lasso <- sbc_rank(rank_b4_lasso, reps = reps, L = L, title = "b4")
sbc_b5_lasso <- sbc_rank(rank_b5_lasso, reps = reps, L = L, title = "b5")
sbc_sigma_lasso <- sbc_rank(rank_sigma_lasso, reps = reps, L = L, title = "sigma")


# both
gridExtra::grid.arrange(sbc_b0_correct, sbc_b1_correct, sbc_b2_correct, 
                        sbc_b3_correct, sbc_b4_correct,
                        sbc_b5_correct, sbc_sigma_correct, ncol = 3)

gridExtra::grid.arrange(sbc_b0_wrong, sbc_b1_wrong, sbc_b2_wrong, 
                        sbc_b3_wrong, sbc_b4_wrong,
                        sbc_b5_wrong, sbc_sigma_wrong, ncol = 3)

gridExtra::grid.arrange(sbc_b0_lasso, sbc_b1_lasso, sbc_b2_lasso, 
                        sbc_b3_lasso, sbc_b4_lasso,
                        sbc_b5_lasso, sbc_sigma_lasso, ncol = 3)



sum(colMeans(finalMatrix_correct[, 16:(15 + n_test)])^2)
sum(colMeans(finalMatrix_wrong[, 16:(15 + n_test)])^2)
sum(colMeans(finalMatrix_lasso[, 16:(15 + n_test)])^2)

sum(abs(colMeans(finalMatrix_correct)[16:(15 + n_test)]))
sum(abs(colMeans(finalMatrix_wrong)[16:(15 + n_test)]))
sum(abs(colMeans(finalMatrix_lasso)[16:(15 + n_test)]))


