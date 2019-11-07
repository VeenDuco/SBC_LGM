
L <- bins - 1

set.seed(1236)

b0 <- rnorm(reps, 1, 1)
b1 <- rnorm(reps, 1, 1)
b2 <- rnorm(reps, 2, 1)
b3 <- rnorm(reps, 3, 1)
b4 <- rnorm(reps, 0, 1)
b5 <- rnorm(reps, 0, 1)
sigma <- rlnorm(reps, 0, 1)

X <- array(NA, dim = c(n_ind, 5, reps))
for(k in 1:reps) X[, , k] <- mvrnorm(n_ind, mu = rep(0, 5), Sigma = diag(5))

ysim <- array(NA, dim = c(n_ind, reps))

for(i in 1:n_ind){
  for(k in 1:reps){
    ysim[i, k] <- rnorm(1, mean = I(1 * b0[k] + 
                                      b1[k] * X[i, 1, k] +
                                      b2[k] * X[i, 2, k] + 
                                      b3[k] * X[i, 3, k] + 
                                      b4[k] * X[i, 4, k] + 
                                      b5[k] * X[i, 5, k]), 
                        sd = sigma[k])
  }
}

X_test <- array(NA, dim = c(n_test, 5, reps))
for(k in 1:reps) X_test[, , k] <- mvrnorm(n_test, mu = rep(0, 5), Sigma = diag(5))

ysim_test <- array(NA, dim = c(n_test, reps))

for(i in 1:n_test){
  for(k in 1:reps){
    ysim_test[i, k] <- rnorm(1, mean = I(1 * b0[k] + 
                                      b1[k] * X_test[i, 1, k] +
                                      b2[k] * X_test[i, 2, k] + 
                                      b3[k] * X_test[i, 3, k] + 
                                      b4[k] * X_test[i, 4, k] + 
                                      b5[k] * X_test[i, 5, k]), 
                        sd = sigma[k])
  }
}
