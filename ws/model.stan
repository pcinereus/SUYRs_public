data{
  int N; // sample size
  int P; // number of species
  int D; // number of dimensions
  int<lower=0> Y[N,P]; // data matrix of order [N,P]
}
transformed data{
  int<lower=1> M;
  M = D*(P-D)+ D*(D-1)/2;  // number of lower-triangular, non-zero loadings
}
parameters{
  //Parameters
  vector[N] alpha; //Row intercepts
  row_vector[P] b0; // Intercept per species
  vector[M] L_lower; //Lower diagonal loadings
  vector<lower=0>[D] L_diag; //Positive diagonal elements of loadings
  matrix[N,D] FS; //Factor scores, matrix of order [N,D]
  cholesky_factor_corr[D] Rho; //Correlation matrix between factors
  //Hyperparameters
  real<lower=0> sigma_a; //Variance of the row intercepts
  real<lower=0> sigma_b0; //Variance of the species intercepts
  real<lower=0> mu_lower;
  real<lower=0> sigma_lower;
  real<lower=0> eta; //Parameter of LKJ prior for Rho
}
transformed parameters{
  cholesky_factor_cov[P,D] lambda; //Final matrix of laodings
  vector[D] FS_mu; // factor means
  vector<lower=0>[D] FS_sd; // factor standard deviations
  matrix[D,D] Ld; // cholesky decomposition of the covariance matrix between factors
  matrix[N,P] temp; //intermediate predictor
  matrix<lower=0>[N,P] mu; // predicted values
  
  for (m in 1:D) {
    FS_mu[m] = 0; //Mean of factors = 0
    FS_sd[m] = 1;} //Sd of factors = 1
  
  // Correlation matrix of factors
  Ld = diag_matrix(FS_sd) * Rho;
  {
    int idx2; //Index for the lower diagonal loadings
    idx2 = 0;
    
    // Constraints to allow identifiability of loadings
  	 for(i in 1:(D-1)) { for(j in (i+1):(D)){ lambda[i,j] = 0; } } //0 on upper diagonal
  	 for(i in 1:D) lambda[i,i] = L_diag[i]; //Positive values on diagonal
  	 for(j in 1:D) {
  	   for(i in (j+1):P) {
  	     idx2 = idx2+1;
  	     lambda[i,j] = L_lower[idx2];
  	   }
  	 }
  }
  
  // Predictor
  temp = FS * lambda';
  for(n in 1:N) mu[n] = exp(alpha[n] + b0 + temp[n]);
  
}
model{
  // Hyperpriors
  sigma_a ~ gamma(2,0.1); //Row intercepts hyperprior
  sigma_b0 ~ gamma(2,0.1); //Species intercept effects hyperprior
  mu_lower ~ gamma(2,0.1);//Mean of lower loadings
  sigma_lower ~ gamma(2,0.1); //Variance of lower loadings
  eta ~ gamma(2,0.1); //Parameter of the cholesky prior for FS correlation structure
  
  // Priors
  alpha ~ normal(0, sigma_a); //Regularizing prior for row intercepts
  b0 ~ normal(0, sigma_b0); //Regularizing prior for species intercepts
  L_diag ~ normal(0,0.1); //Prior for diagonal elements
  //L_lower ~ normal(mu_lower, sigma_lower); //Regularizing prior for lower loadings
  L_lower ~ normal(0, sigma_lower); //Regularizing prior for lower loadings
  Rho ~ lkj_corr_cholesky(eta); //Regularizing prior for Rho
  
  for(i in 1:N){	
    Y[i,] ~ poisson(mu[i]);
    FS[i,] ~ multi_normal_cholesky(FS_mu, Ld);
  }
}
generated quantities{
  //cov_matrix[P] cov_lambda;
  
  //cov_lambda = lambda * lambda´
}