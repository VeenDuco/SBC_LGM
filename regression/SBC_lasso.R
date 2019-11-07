finalMatrix_lasso  <- foreach(i = 1:reps, .packages = c("rstan", "tcltk"), 
                              .combine = rbind,
                              .verbose = TRUE) %dopar% {
                                
                                if(!exists("pb")) pb <- tkProgressBar("Parallel task", min = 1, max = reps)
                                setTkProgressBar(pb, i)  
                                
                                thin <- 1
                                draws <- L
                                warm <- warm.init
                                counter <- 1
                                fit <- sampling(lasso_reg, 
                                                data = list(X_train = X[, , i], N_train = n_ind, p = 5, y_train = ysim[, i]),
                                                warmup = warm, iter = warm + draws, seed = i, chains = 1)
                                
                                ## note 2 * L, was L but still signs of autocorrelation.
                                while(min(na.omit(summary(fit)$summary[, "n_eff"])) < 2 * L | max(na.omit(summary(fit)$summary[, "Rhat"])) > 1.02) { 
                                  counter <- counter + 1
                                  thin <- thin * 2
                                  draws <- draws * 2
                                  warm <- warm * 1.1
                                  fit <- sampling(lasso_reg, 
                                                  data = list(X_train = X[, , i], N_train = n_ind, p = 5, y_train = ysim[, i]),
                                                  warmup = warm, iter = warm + draws, seed = i, chains = 1)
                                }
                                
                                # saveRDS(fit, paste0("temp/fit_", i ,".rds"))
                                
                                fit_b0 <- rstan::extract(fit, pars = "mu")[[1]]
                                fit_b1 <- rstan::extract(fit, pars = "beta[1]")[[1]]
                                fit_b2 <- rstan::extract(fit, pars = "beta[2]")[[1]]
                                fit_b3 <- rstan::extract(fit, pars = "beta[3]")[[1]]
                                fit_b4 <- rstan::extract(fit, pars = "beta[4]")[[1]]
                                fit_b5 <- rstan::extract(fit, pars = "beta[5]")[[1]]
                                fit_sigma <- rstan::extract(fit, pars = "sigma")[[1]]
                                
                                
                                
                                fit_b0 <- fit_b0[(1:L) * thin]
                                fit_b1 <- fit_b1[(1:L) * thin]
                                fit_b2 <- fit_b2[(1:L) * thin]
                                fit_b3 <- fit_b3[(1:L) * thin]
                                fit_b4 <- fit_b4[(1:L) * thin]
                                fit_b5 <- fit_b5[(1:L) * thin]
                                fit_sigma <- fit_sigma[(1:L) * thin]
                                
                                ## bias
                                bias_b0 <- summary(fit, pars = "mu")$summary[1] - b0[i]
                                bias_b1 <- summary(fit, pars = "beta[1]")$summary[1] - b1[i]
                                bias_b2 <- summary(fit, pars = "beta[2]")$summary[1] - b2[i]
                                bias_b3 <- summary(fit, pars = "beta[3]")$summary[1] - b3[i]
                                bias_b4 <- summary(fit, pars = "beta[4]")$summary[1] - b4[i]
                                bias_b5 <- summary(fit, pars = "beta[5]")$summary[1] - b5[i]
                                bias_sigma <- summary(fit, pars = "sigma")$summary[1] - sigma[i]
                                
                                
                                tempMatrix <- matrix(NA, ncol = 15, nrow = 1)
                                
                                tempMatrix[1, 1] <- sum(fit_b0 < b0[i])
                                tempMatrix[1, 2] <- sum(fit_b1 < b1[i])
                                tempMatrix[1, 3] <- sum(fit_b2 < b2[i])
                                tempMatrix[1, 4] <- sum(fit_b3 < b3[i])
                                tempMatrix[1, 5] <- sum(fit_b4 < b4[i])
                                tempMatrix[1, 6] <- sum(fit_b5 < b5[i])
                                tempMatrix[1, 7] <- sum(fit_sigma < sigma[i])
                                
                                tempMatrix[1, 8] <- bias_b0
                                tempMatrix[1, 9] <- bias_b1
                                tempMatrix[1, 10] <- bias_b2
                                tempMatrix[1, 11] <- bias_b3
                                tempMatrix[1, 12] <- bias_b4
                                tempMatrix[1, 13] <- bias_b5
                                tempMatrix[1, 14] <- bias_sigma
                                
                                tempMatrix[1, 15] <- counter
                                
                                tempMatrix
                                
                              }
