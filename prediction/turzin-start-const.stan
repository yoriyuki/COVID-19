functions{
  real[] ni_eq(real t,        // time
           real[] y,      // state
           real[] theta,  // parameters
           real[] x_r,    // data (real)
           int[] x_i) {   // data (integer)
  real dydt[2];
  real a, d, b, P;
  r = theta[1];
  b = theta[2];
  P = x_r[1];
  dydt[1] = y[1] * b * y[2] * (1 - y[2]/P) - r*y[1];
  dydt[2] = - y[1];
  return dydt;
}
  vector ni(real a, real d, real b, real init_inf, P, int T){
    real NI[T]
    real inital_state[2]
    real times[T-1];
    real theta[2];
    real x_r[1];
    int x_i[1];
    initial_state[1] = init_inf;
    initial_state[2] = P;
    initial_time = 0;
    for (t in 1:T-1){
      times[t] = t;
    }
    theta[1] = a+d;
    theta[2] = b;
    x_r[1] = P;
    x_i[1] = 0;
    NI = integrate_ode_rk45(ni_eq, real[] initial_state, real initial_time, real[] times, real[] theta, real[] x_r, int[] x_i)
}
data {
  int<lower=0> T; // Time horizon
  int<lower=0> P; // Population
  int C0[T]; // Cummulative infection
  int R0[T]; // Recovered
  int D0[T]; // Cummulative death

}
parameters {
  real<lower=0> init_inf;
  real<lower=0> b;
  real<lower=0, upper=1> q;

  real<lower=0, upper=1> a;
  real<lower=0, upper=1> d;
  }  
  transformed parameters {
  vector[T-1] NI;
  vector[T-1] q;
  NI = ni(a, d, b0, b1, theta_b, b_date, init_inf, P, T);
}    
model{
    a ~ beta(0.1, 1);
    d ~ beta(5, 1000);
    b ~ gamma(0.1, 1);
    q ~ beta(1, 1);
    init_inf ~ gamma(1, 1);
    C0[1] ~ poisson(init_inf);

    for (t in 1:T-1){
      C0[t+1] - C0[t] ~ poisson(q * NI[t]);
      D0[t+1] - D0[t] ~ poisson(d * (C0[t] - R0[t] - D0[t]));   
      R0[t+1] - R0[t] ~ poisson(a * (C0[t] - R0[t] - D0[t]));
    }
}
generated quantities {
  vector[T] log_lik;
  
  for (t in 1:T-1) {
    log_lik[t] = poisson_lpmf(C0[t+1] - C0[t] | q * NI[t]);
  }
}