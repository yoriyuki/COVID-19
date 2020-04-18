functions{
  vector ni(real a, real d, 
            real b0, real b1, real theta_b, real b_date, 
            real b2, real theta_b2, real b2_date, 
            real init_inf, real P, int T){
    real C;
    real R;
    real I;
    real NR;
    real ND;
    real D;
    real b;
    vector[T-1] NI;

    I = init_inf;
    R = 0;
    D = 0;
    C = 0;
    for (t in 1:T-1){
      b = b0 + (b1 - b0) * inv_logit(theta_b * (t - b_date)) + (b2 - b1) * inv_logit(theta_b2 * (t - b2_date)) ;
      NI[t] = I * b * (1 - C/P);
      NR = a * I;
      ND = d * I;
      D = D + ND;
      I = I + NI[t] - NR - ND;
      C = C + NI[t];
      R = R + NR;
    }
    return NI;
  }
}
data {
  int<lower=0> T; // Time horizon
  int<lower=0> P; // Population
  int C0[T]; // Cummulative infection
  int R0[T]; // Recovered
  int D0[T]; // Cummulative death
  real ED; // Emergency Declaration

}
parameters {
  real<lower=0> init_inf;
  real<lower=0> b0;
  real<lower=0> b1;
  real<lower=0> theta_b;
  real<lower=0, upper=T> b_date;
  real<lower=0> b2;
  real<lower=0> theta_b2;
  real<lower=0, upper=T> b2_date;
  real<lower=0, upper=1> q0;
  real<lower=0, upper=1> q1;
  real<lower=0> theta_q;
  real<lower=0, upper=T> q_date;
  real<lower=0, upper=1> q2;
  real<lower=0> theta_q2;
  real<lower=0, upper=T> q2_date;

  real<lower=0, upper=1> a;
  real<lower=0, upper=1> d;
  }  
transformed parameters {
  vector[T-1] NI;
  vector[T-1] q;
   for (t in 1:T-1){
    q[t] = q0 + (q1 - q0) * inv_logit(theta_q * (t - q_date)) + (q2 - q1) * inv_logit(theta_q2 * (t - q2_date));
  }
  NI = ni(a, d, b0, b1, theta_b, b_date, b2, theta_b2, b2_date, init_inf, P, T);
}
model {
    a ~ beta(0.1, 1);
    d ~ beta(0.001, 1);

    b0 ~ gamma(0.1, 1);
    b1 ~ gamma(0.1, 1);
    theta_b ~ gamma(5, 1);
    b_date ~ uniform(0, T);
    b2 ~ gamma(0.1, 1);
    theta_b2 ~ gamma(5, 1);
    b2_date ~ uniform(b_date, T);

    q0 ~ beta(1, 1);
    q1 ~ beta(1, 1);
    theta_q ~ gamma(5, 1);
    q_date ~ uniform(0, T);
    q2 ~ beta(1, 1);
    theta_q2 ~ gamma(5, 1);
    q2_date ~ uniform(q_date, T);

    init_inf ~ gamma(1, 1);
    C0[1] ~ poisson(q[1] * init_inf);
    for (t in 1:T-1){
      C0[t+1] - C0[t] ~ poisson(q[t] * NI[t]);
      D0[t+1] - D0[t] ~ poisson(d * (C0[t] - R0[t] - D0[t]));   
      R0[t+1] - R0[t] ~ poisson(a * (C0[t] - R0[t] - D0[t]));
    }
}
generated quantities {
  vector[T-1] log_lik;
  
  for (t in 1:T-1) {
    log_lik[t] = poisson_lpmf(C0[t+1] - C0[t] | q[t] * NI[t]);
  }
}