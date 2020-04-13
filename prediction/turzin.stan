data {
  int<lower=0> T; // Time horizon
  int<lower=0> P; // Population
  int S0[T]; // Cummulative infection
  int R0[T]; // Recovered
  int D0[T]; // Cummulative death

}
parameters {
  real<lower=0> init_inf_alpha;
  real<lower=0> init_inf_beta;
  real<lower=0> init_inf;
  real<lower=0> b0_alpha; // initial growth rate
  real<lower=0> b0_beta;
  real<lower=0> b0;
  real<lower=0> b1_alpha; // final growth rate
  real<lower=0> b1_beta;
  real<lower=0> b1;
  real<lower=0> theta_b_alpha; // curve steepness
  real<lower=0> theta_b_beta;
  real<lower=0> theta_b;
  real<lower=0> b2_alpha; // final growth rate
  real<lower=0> b2_beta;
  real<lower=0> b2;
  real<lower=0> theta_b2_alpha; // curve steepness
  real<lower=0> theta_b2_beta;
  real<lower=0> theta_b2;
  real<lower=0, upper=T> b_date;
  real<lower=0, upper=T> b2_date;
  real<lower=0> q0_alpha; // initial detection rate
  real<lower=0> q0_beta;
  real<lower=0, upper=1> q0;
  real<lower=0> q1_alpha; // final detection rate
  real<lower=0> q1_beta;
  real<lower=0, upper=1> q1;
  real<lower=0> theta_q_alpha; // curve steepness
  real<lower=0> theta_q_beta;
  real<lower=0> theta_q;
  real<lower=0, upper=T> q_date;
  real<lower=0> q2_alpha; // final detection rate
  real<lower=0> q2_beta;
  real<lower=0, upper=1> q2;
  real<lower=0> theta_q2_alpha; // curve steepness
  real<lower=0> theta_q2_beta;
  real<lower=0> theta_q2;
  real<lower=0, upper=T> q2_date;

  real<lower=0> a_alpha; // Recovery rate
  real<lower=0> a_beta;
  real<lower=0, upper=1> a;
  real<lower=0> d_alpha; // Death rate
  real<lower=0> d_beta;
  real<lower=0, upper=1> d;
  }  
model {
    real S;
    real R;
    real I;
    real NR;
    real ND;
    real D;
    real b; 
    real q;
    a ~ beta(a_alpha, a_beta);
    d ~ beta(d_alpha, d_beta);

    b0 ~ gamma(b0_alpha, b0_beta);
    b1 ~ gamma(b1_alpha, b1_beta);
    b2 ~ gamma(b2_alpha, b2_beta);
    theta_b ~ gamma(theta_b_alpha, theta_b_beta);
    theta_b2 ~ gamma(theta_b2_alpha, theta_b2_beta);
    b_date ~ uniform(0, T);
    b2_date ~ uniform(b_date, T);

    q0 ~ beta(q0_alpha, q0_beta);
    q1 ~ beta(q1_alpha, q1_beta);
    theta_q ~ gamma(theta_q_alpha, theta_q_beta);
    q_date ~ uniform(0, T);
    q2 ~ beta(q2_alpha, q2_beta);
    theta_q2 ~ gamma(theta_q_alpha, theta_q_beta);
    q2_date ~ uniform(q_date, T);

    init_inf ~ gamma(init_inf_alpha, init_inf_beta);
    S = init_inf;
    R = 0;
    D = 0;
    
    for (t in 2:T){
      b = b0 + (b1-b0) * inv_logit(theta_b * (t - b_date)) + (b2 - b1) * inv_logit(theta_b2 * (t - b2_date));
      I = (S - R - D) * b * (1 - S/P);
      q = q0 + (q1-q0) * inv_logit(theta_q * (t - q_date)) + (q2-q1) * inv_logit(theta_q2 * (t - q2_date));
      NR = a * (S - R - D);
      ND = d * (S - R - D);
      D = D + ND;
      S = S + I;
      R = R + NR;
      S0[t] - S0[t-1] ~ poisson(q * I);
      D0[t] - D0[t-1] ~ poisson(ND);      
      R0[t] - R0[t-1] ~ poisson(a * (S0[t-1] - R0[t-1] - D0[t]));
  }
}