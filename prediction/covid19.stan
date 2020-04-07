data {
  int<lower=0> T; // Time horizon
  int<lower=0> N; // Number of countries
  int<lower=1> L; // Delay
  matrix<lower=0>[T, N] S; // Cummulative infection
  matrix<lower=0>[T, N] R; // Recovered
}
parameters {
  vector<lower=0>[N] c_alpha; // Reproduction matrix
  vector<lower=0>[N] c_beta;
  vector<lower=0>[N] c;
  vector<lower=0>[N] p_alpha; // Effective population
  vector<lower=0>[N] p_beta;
  vector<lower=0>[N] p;
  vector<lower=0>[N] a_alpha; // Recovery rate
  vector<lower=0>[N] a_beta;
  vector<lower=0, upper=1>[N] a;
  vector<lower=0>[N] sigma_S; // noise factor for the cumulative infection
  vector<lower=0>[N] sigma_R; // noise factor for the recorvery
}  
model {
    a ~ beta(a_alpha, a_beta);
    p ~ gamma(p_alpha, p_beta);
    c ~ gamma(c_alpha, c_beta);
    for (t in L+1:T){
      S[t] ~ normal(S[t-1] + (S[t-L] - R[t-L]) .* to_row_vector(c) .* (1 - S[t-L] ./ to_row_vector(p)), to_row_vector(sigma_S));
      R[t] ~ normal(R[t-1] + to_row_vector(a) .* (S[t-1] - R[t-1]), to_row_vector(sigma_R));
  }
}