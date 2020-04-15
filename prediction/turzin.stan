data {
  int<lower=0> T; // Time horizon
  int<lower=0> P; // Population
  int C0[T]; // Cummulative infection
  int R0[T]; // Recovered
  int D0[T]; // Cummulative death

}
parameters {
  real<lower=0, upper=1> p;

  // Constant b and q
  real<lower=0> init_infc;
  real<lower=0> bc;
  real<lower=0, upper=1> qc;
  real<lower=0, upper=1> ac;
  real<lower=0, upper=1> dc;

  // moving b and q
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
    real q;

    real lp1;
    real lp2;

    p ~ beta(1, 2);

    // Fixed 
    ac ~ beta(0.1, 1);
    dc ~ beta(1, 1000);

    bc ~ gamma(0.1, 1);

    qc ~ beta(1, 2);
    
    // Changing b and q
    a ~ beta(0.1, 1);
    d ~ beta(1, 1000);

    b0 ~ gamma(0.1, 1);
    b1 ~ gamma(0.1, 1);
    theta_b ~ gamma(1, 1);
    b_date ~ uniform(0, T);

    q0 ~ beta(1, 2);
    q1 ~ beta(1, 2);
    theta_q ~ exponential(1);
    q_date ~ uniform(0, T);

    lp1 = 0;
    lp2 = 0;

    init_infc ~ gamma(1, 1);
    I = init_infc;
    R = 0;
    D = 0;
    C = 0;
    //constant b and q
    for (t in 2:T){
      NI = I * bc * (1 - C/P);
      NR = ac * I;
      ND = dc * I;
      D = D + ND;
      I = I + NI - NR - ND;
      C = C + NI;
      R = R + NR;

      //loglikellyfood
      lp1 = log_sum_exp(lp1, poisson_lpmf(C0[t]-C0[t-1]| qc * NI));
      lp1 = log_sum_exp(lp1, poisson_lpmf(D0[t]-D0[t-1]| ND));
      lp1 = log_sum_exp(lp1, poisson_lpmf(R0[t]-R0[t-1]| ac * (C0[t-1] - R0[t-1] - D0[t-1])));
    }

    init_inf ~ gamma(1, 1);
    I = init_inf;
    R = 0;
    D = 0;
    C = 0;
    //moving b and q
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
      //loglikellyfood
      lp2 = log_sum_exp(lp2, poisson_lpmf(C0[t]-C0[t-1]| q * NI));
      lp2 = log_sum_exp(lp2, poisson_lpmf(D0[t]-D0[t-1]| ND));
      lp2 = log_sum_exp(lp2, poisson_lpmf(R0[t]-R0[t-1]| a * (C0[t-1] - R0[t-1] - D0[t-1])));
  }
  target += log_sum_exp(log(p) + lp1, log(1-p) + lp2);
}