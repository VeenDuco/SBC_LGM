data{
	int N_train; //number of observations training and validation set
	int p; //number of predictors
	real y_train[N_train]; //response vector
	matrix[N_train, p] X_train; //model matrix
	//test set
}
parameters{
	real mu; //intercept
	real<lower = 0> sigma2; //error variance
	vector[p] beta; // regression parameters
	//hyperparameters prior
	real<lower = 0> lambda; //penalty parameter
}
transformed parameters{
	real<lower = 0> sigma; //error sd
	vector[N_train] linpred; //mean normal model
	sigma = sqrt(sigma2);
	linpred = mu + X_train*beta;
}
model{
 //prior regression coefficients: lasso
	beta ~ double_exponential(0, sigma/lambda);
	lambda ~ cauchy(0, 1);
	
 //priors nuisance parameters: uniform on log(sigma^2) & mu
	target += -2 * log(sigma); 
	
 //likelihood
	y_train ~ normal(linpred, sigma);
}
