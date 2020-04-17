data {
  int<lower=0> T; // Time horizon
  int<lower=0> P; // Population
  int C0[T]; // Cummulative infection
  int R0[T]; // Recovered
  int D0[T]; // Cummulative death

}
parameters {
  real<lower=0> init_inf;
  real<lower=0> b0;
  real<lower=0> b1;
  real<lower=0> theta_b;
  real<lower=0, upper=T> b_date;
  real<lower=0, upper=1> q0;
  real<lower=0, upper=1> q1;
  real<lower=0> theta_q;
  real<lower=0, upper=T> q_date;

  real<lower=0, upper=1> a;
  real<lower=0, upper=1> d;
  }  
model {
    real C;
    real R;
    real I;
    real NI;
    real NR;
    real ND;
    real D;
    real b;
    real bb; 
    real q;
    a ~ beta(0.1, 1);
    d ~ beta(5, 1000);

    b0 ~ gamma(0.1, 1);
    b1 ~ gamma(0.1, 1);
    theta_b ~ gamma(1, 2);
    b_date ~ uniform(0, T);

    q0 ~ beta(1, 1);
    q1 ~ beta(1, 1);
    theta_q ~ gamma(1, 2);
    q_date ~ uniform(0, T);

    init_inf ~ gamma(1, 1);
    I = init_inf;
    R = 0;
    D = 0;
    C = 0;
    
    for (t in 2:T){
      b = b0 + (b1 - b0) * inv_logit(theta_b * (t - b_date));
      NI = I * b * (1 - C/P);
      q = q0 + (q1 - q0) * inv_logit(theta_q * (t - q_date));
      NR = a * I;
      ND = d * I;
      D = D + ND;
      I = I + NI - NR - ND;
      C = C + NI;
      R = R + NR;
      C0[t] - C0[t-1] ~ poisson(q * NI);
      D0[t] - D0[t-1] ~ poisson(d * (C0[t-1] - R0[t-1] - D0[t-1]));   
      R0[t] - R0[t-1] ~ poisson(a * (C0[t-1] - R0[t-1] - D0[t-1]));
  }
}