data {
  int<lower = 0> n_ind;
  vector[n_ind] y;
  matrix[n_ind, 5] X;
  //test set
	int N_test; //number of observations test set
	matrix[N_test, 5] X_test; //model matrix test set

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
  target += normal_lpdf(b0 | 1, 10);
  target += normal_lpdf(b1 | 1, 10);
  target += normal_lpdf(b2 | 2, 10);
  target += normal_lpdf(b3 | 3, 10);
  target += normal_lpdf(b4 | 0, 10);
  target += normal_lpdf(b5 | 0, 10);
  target += lognormal_lpdf(sigma | 0, 1);
  
  for(i in 1:n_ind){
    target += normal_lpdf(y[i] | b0 + b1 * X[i, 1] + b2 * X[i, 2] + b3 * X[i, 3] + b4 * X[i, 4] + b5 * X[i, 5], sigma);
  }
  
}


generated quantities{ //predict responses test set
	real y_test[N_test]; //predicted responses
	for(i in 1:N_test){
		y_test[i] = normal_rng(b0 + b1 * X_test[i, 1] + b2 * X_test[i, 2] + b3 * X_test[i, 3] + b4 * X_test[i, 4] + b5 * X_test[i, 5], sigma);
	}
}	
