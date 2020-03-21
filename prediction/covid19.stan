data {
  int<lower=0> T; // Time horizon
  int<lower=0> N; // Number of countries
  int<lower=0> L; // Period of auto regression
  matrix<lower=0>[T, N] S; // Cummulative infection
  matrix<lower=0>[T, N] R; // Recovered
}
parameters {
  matrix<lower=0>[N, N] c_alpha[L]; // Reproduction matrix
  matrix<lower=1>[N, N] c_beta[L];
  matrix<lower=0>[N, N] c[L];
  vector<lower=0>[N] p_alpha; // Effective population
  vector<lower=1>[N] p_beta;
  vector<lower=0>[N] p;
  real<lower=0> a_alpha; // Recovery rate
  real<lower=1> a_beta;
  real<lower=0, upper=1> a;
  real<lower=0> sigma_S; // noise factor for the cumulative infection
  real<lower=0> sigma_R; // noise factor for the recorvery
}  
transformed parameters {

}
model {
    row_vector[N] r;
    a ~ beta(a_alpha, a_beta);
    p ~ gamma(p_alpha, p_beta);
    for (i in 1:N) {
      for (j in 1:N) {
        for (k in 1:L){
          c[k][i, j] ~ gamma(c_alpha[k][i, j], c_beta[k][i, j]);
        }
      };
    };
    for (t in L+1:T){
      r = (S[t-L] - R[t-L]) * c[L];
      for (k in 1:L-1){
        r += ((S[t-k] - R[t-k]) - (S[t-k-1] - R[t-k-1])) * c[k];
      };
      S[t] ~ normal(S[t-1] + r .* (1 - S[t-1]./to_row_vector(p)), sigma_S);
      R[t] ~ normal(R[t-1] + a * (S[t-1] - R[t-1]), sigma_R);
  }
}