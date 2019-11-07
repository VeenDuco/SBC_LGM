data {
  int<lower = 0> n_ind;
  vector[n_ind] y;
  matrix[n_ind, 5] X;
}

parameters {
  real b0;
  real b1;
  real b2;
  real b3;
  real b4;
  real b5;
  real sigma;
}

model {
  target += normal_lpdf(b0 | 1, 1);
  target += normal_lpdf(b1 | 1, 1);
  target += normal_lpdf(b2 | 2, 1);
  target += normal_lpdf(b3 | 3, 1);
  target += normal_lpdf(b4 | 0, 1);
  target += normal_lpdf(b5 | 0, 1);
  target += lognormal_lpdf(sigma | 0, 1);
  
  for(i in 1:n_ind){
    target += normal_lpdf(y[i] | b0 + b1 * X[i, 1] + b2 * X[i, 2] + b3 * X[i, 3] + b4 * X[i, 4] + b5 * X[i, 5], sigma);
  }
  
}


