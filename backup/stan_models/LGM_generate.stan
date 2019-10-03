data {
  int n_ind;
  int n_time;
}


parameters {
  
}

transformed parameters {
  
}

model {
  
}

generated quantities {
  
  matrix[n_ind, n_time] y_rep;
  // latent variable structure
  real alpha1 = normal_rng(30, 10);
  real alpha2 = normal_rng(-5, 5);
  real theta = lognormal_rng(0, 1);
  //real psi11 = fabs(normal_rng(0, 2.5));
  //real psi22 = fabs(normal_rng(0, 2.5));
  //real psi21 = normal_rng(0, 2.5);
  cholesky_factor_corr[2] Lcorr = lkj_corr_cholesky_rng(2, 1);
  vector[2] sigma;  // prior on the standard deviations
  matrix[2,2] Omega;
  matrix[2,2] Psi;
  
  //cov_matrix[2] psi = wishart_rng(1,);
  //matrix[2, 2] psi;
  vector[2] tau[n_ind];
  //vector[n_ind] tau_i;
  //vector[n_ind] tau_s;
  
  vector[n_ind] Intercept;
  vector[n_ind] Slope;
  
  vector[2] zero_combined[n_ind];  // zero matrix mvt_normal to get (co)variances latent variables
  //vector[2] latent_variables[n_ind];  // response matrix latent variables
  vector[n_ind] zero_i = rep_vector(0, n_ind);
  vector[n_ind] zero_s = rep_vector(0, n_ind);
  
  //psi[1, 1] = psi11;
  //psi[2, 2] = psi22;
  //psi[2, 1] = psi21;
  //psi[1, 2] = psi21;
  sigma[1] = lognormal_rng(0, 1);
  sigma[2] = lognormal_rng(0, 1);
  Omega = multiply_lower_tri_self_transpose(Lcorr);
  Psi = quad_form_diag(Omega, sigma); 
  
  for (n in 1:n_ind) {
    //latent_variables[n] = [tau_i[n], tau_s[n]]';
    zero_combined[n] = [zero_i[n], zero_s[n]]';
  }
  
  
  //latent_variables = multi_normal_rng(zero_combined, Psi);
  tau = multi_normal_rng(zero_combined, Psi);
  
  
  for(i in 1:n_ind) {
    Intercept[i] = alpha1 + tau[i, 1];
    Slope[i] = alpha2 + tau[i, 2];
  }
  
  
  
  for(i in 1:n_ind){
    y_rep[i, 1] = normal_rng( 1 * Intercept[i] + 0.083 * Slope[i], theta);
    y_rep[i, 2] = normal_rng( 1 * Intercept[i] + 0.25 * Slope[i], theta);
    y_rep[i, 3] = normal_rng( 1 * Intercept[i] + 1 * Slope[i], theta);
  }
  
  }


