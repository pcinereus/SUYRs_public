data {
  int<lower=1> N; //Number of samples
  int<lower=1> S; //Number of species
  int<lower=1> D; //Number of latent dimensions

  array[N, S] int Y; //Species matrix
}
transformed data{
  // Number of non-zero lower triangular factor loadings
  // Ensures identifiability of the model - no rotation of factors
  int<lower=1> M; 
  M = D * (S - D) + D * (D - 1) / 2 + D;
}
parameters {
  // Site intercepts
  real a_bar;
  real<lower=0> sigma_a;
  vector[N] a;
  
  // Species intercepts
  real<lower=0> sigma_b0;
  vector[S] b0;
  
  // Factor parameters
  vector[M] L; // lower triangle of species loadings
  real<lower=0> sigma_L; // variance of species loadings
  
  // Latent variables
  matrix[D, N] LV_uncor; // Per-site latent variable

  // NegBin parameters
  real<lower=0> kappa;
}
transformed parameters {
  // Construct factor loading matrix
  matrix[S, D] Lambda_uncor;
  // Constraints to allow identifiability of loadings
  for (i in 1:(D-1)) { 
    for (j in (i+1):(D)){ 
      Lambda_uncor[i,j] = 0; 
    } 
  }
  {
    int index;
    index = 0;
    for (j in 1:D) {
      for (i in j:S) {
        index = index + 1;
        Lambda_uncor[i, j] = L[index];
      } 
    }
  }
}
model {
  // Factor priors
  to_vector(LV_uncor) ~ std_normal();
  L ~ std_normal();

  // Random effect priors
  a ~ std_normal();
  b0 ~ std_normal();
  a_bar ~ std_normal();
  sigma_a ~ exponential(1);
  sigma_b0 ~ exponential(1);
  sigma_L ~ exponential(1);
  
  // Negative Binomial scale parameter
  kappa ~ exponential(1);
  
  array[N] vector[S] mu;
  for (i in 1:N) {
      mu[i,] = a_bar + a[i] * sigma_a + b0 * sigma_b0 + (Lambda_uncor * sigma_L) * LV_uncor[,i];
      Y[i,] ~ neg_binomial_2_log(mu[i, ], kappa);
  }
}
generated quantities {
  // Sign correct factor loadings and factors
  matrix[D, N] LV;
  matrix[S, D] Lambda;
  for(d in 1:D){
    if(Lambda_uncor[d,d] < 0){
      Lambda[,d] = -1 * Lambda_uncor[,d];
      LV[d,] = -1 * LV_uncor[d,];
    } else {
      Lambda[,d] = Lambda_uncor[,d];
      LV[d,] = LV_uncor[d,];
    }
  }
  
  // Calculate species covariance matrix
  matrix[S, S] COV;
  COV = multiply_lower_tri_self_transpose(Lambda);
}
