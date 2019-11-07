data {
  int n_ind;
  int n_time;
  matrix[n_ind, n_time] y;
  vector[n_time] x;
}


parameters {
  // latent variable structure
  real alpha1;
  real alpha2;
  vector[2] tau[n_ind];
  //cholesky_factor_corr[2] Lcorr;
  //vector<lower = 0>[2] sigma;  
  cov_matrix[2] Psi;
  
  // residuals observed variables
  real<lower = 0> theta[3];
  //real<lower = 0> theta22;
  //real<lower = 0> theta33;
  
}

transformed parameters {
  vector[n_ind] Intercept;
  vector[n_ind] Slope;
  //matrix[2,2] Omega;
  
  
  
  for(i in 1:n_ind) {
    Intercept[i] = alpha1 + tau[i, 1];
    Slope[i] = alpha2 + tau[i, 2];
  }
  
  //Omega = multiply_lower_tri_self_transpose(Lcorr);
  //Psi = quad_form_diag(Omega, sigma); 
}

model {
  // zero structure mvt normal latent variables
  vector[2] zero_combined[n_ind];  // zero matrix mvt_normal to get (co)variances latent variables
  vector[2] latent_variables[n_ind];  // response matrix latent variables
  vector[n_ind] zero_i = rep_vector(0, n_ind);
  vector[n_ind] zero_s = rep_vector(0, n_ind);
  matrix[2,2] q;
  q[1,1] = 1;
  q[2,1] = 0;
  q[2,2] = 1;
  q[1,2] = 0;
    
  for (n in 1:n_ind) {
    zero_combined[n] = [zero_i[n], zero_s[n]]';
  }

    // priors
    //target += lognormal_lpdf(sigma[1] | 0, 1);
    //target += lognormal_lpdf(sigma[2] | 0, 1);
    //target += lkj_corr_cholesky_lpdf(Lcorr | 1);
   
    
    target += inv_wishart_lpdf(Psi | 3, q);
    
    
    // target += normal_lpdf(tau_i | 0, psi11);
    target += normal_lpdf(alpha1 | 30, 10);
    target += normal_lpdf(alpha2 | -5, 5);
    target += lognormal_lpdf(theta | 0, 1);
    //target += normal_lpdf(theta22 | 0, 5);
    //target += normal_lpdf(theta33 | 0, 5);
    
    target += multi_normal_lpdf(tau | zero_combined, Psi);
    //latent_variables ~ multi_normal(zero_combined, psi);
    
    

    for(i in 1:n_ind){
    for(j in 1:n_time) {
    target += normal_lpdf(y[i, j] | 1 * Intercept[i] + x[j] * Slope[i], theta[j]); 
    }
    }

    
  }
    
  