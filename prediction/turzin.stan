data {
  int<lower=0> T; // Time horizon
  int<lower=0> P; // Population
  vector<lower=0>[T] S0; // Cummulative infection
  vector<lower=0>[T] R0; // Recovered
  vector<lower=0>[T] D0; // Cummulative death

}
parameters {
  vector<lower=0, upper=P>[T] S; // real cummulative infection
  vector<lower=0, upper=P>[T] R; // real recovered
  real<lower=0> b_alpha; // growth rate
  real<lower=0> b_beta;
  real<lower=0> b;
  real<lower=0> q_alpha; // detection rate
  real<lower=0> q_beta;
  real<lower=0, upper=1> q;
  real<lower=0> a_alpha; // Recovery rate
  real<lower=0> a_beta;
  real<lower=0, upper=1> a;
  real<lower=0> d_alpha; // Death rate
  real<lower=0> d_beta;
  real<lower=0, upper=1> d;
  real<lower=0> sigma_S_alpha;
  real<lower=0> sigma_S_beta;
  real<lower=0> sigma_S; // noise factor for cumulative infection
  real<lower=0> sigma_R_alpha;
  real<lower=0> sigma_R_beta;
  real<lower=0> sigma_R; // noise factor for recorvery
  real<lower=0> sigma_D_alpha;
  real<lower=0> sigma_D_beta;
  real<lower=0> sigma_D; // noise factor for death
  real<lower=0> sigma_S0_alpha;
  real<lower=0> sigma_S0_beta;
  real<lower=0> sigma_S0; // observation error for cumulative infection
  // real<lower=0> sigma_R0; // noise factor for detected recorvery
}  
model {
    a ~ beta(a_alpha, a_beta);
    d ~ beta(d_alpha, d_beta);
    b ~ gamma(b_alpha, b_beta);
    q ~ beta(q_alpha, q_beta);
    sigma_S ~ gamma(sigma_S_alpha, sigma_S_beta);
    sigma_R ~ gamma(sigma_R_alpha, sigma_R_beta);
    sigma_D ~ gamma(sigma_D_alpha, sigma_D_beta);
    sigma_S0 ~ gamma(sigma_S0_alpha, sigma_S0_beta);
    
    for (t in 2:T){
      S[t] ~ normal(S[t-1] + (S[t-1] - R[t-1] - D0[t-1]) * b * (1 - S[t-1]/P), sigma_S);
      R[t] ~ normal(R[t-1] + a * (S[t-1] - R[t-1] - D0[t-1]), sigma_R);
      D0[t] ~ normal(D0[t-1] + d * (S[t-1] - R[t-1] - D0[t-1]), sigma_D);
      S0[t] ~ normal(q * S[t], sigma_S0);
      R0[t] ~ normal(a * S0[t], 1);
  }
}