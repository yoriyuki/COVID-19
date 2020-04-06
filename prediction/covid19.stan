data {
  int<lower=0> T; // Time horizon
  int<lower=0> N; // Number of countries
  int<lower=1> L; // Delay
  matrix<lower=0>[T, N] S; // Cummulative infection
  matrix<lower=0>[T, N] R; // Recovered
}
parameters {
  matrix<lower=0>[N, N] c_alpha; // Reproduction matrix
  matrix<lower=0>[N, N] c_beta;
  matrix<lower=0>[N, N] c;
  vector<lower=0>[N] p_alpha; // Effective population
  vector<lower=0>[N] p_beta;
  vector<lower=0>[N] p;
  real<lower=0> a_alpha; // Recovery rate
  real<lower=0> a_beta;
  real<lower=0, upper=1> a;
  vector<lower=0>[N] sigma_b; //noise on infection rate
  vector<lower=0>[N] sigma_S; //observation error on S
  vector<lower=0>[N] sigma_R; //observation error on R
}  
model {
    p ~ gamma(p_alpha, p_beta);
    a ~ beta(a_alpha, a_beta);
    for (i in 1:N) {
      for (j in 1:N) {
          c[i, j] ~ gamma(c_alpha[i, j], c_beta[i, j]);
      };
    };
    for (t in L+1:T){
      S[t] ~ normal(S[t-1] + (S[t-L] - R[t-L]) * c .* (1 - S[t-L]./to_row_vector(p)), sigma_S);
      R[t] ~ normal(R[t-1] + a * (S[t-1] - R[t-1]), sigma_R);
  }
}