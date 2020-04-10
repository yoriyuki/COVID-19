data {
    int<lower=0> T; // Time horizon
    int<lower=0> N; // Population
    vector<lower=0>[T] I0; // observed infected 
    vector<lower=0>[T] R0; // observed recovery
    vector<lower=0>[T] D; // Death
}
parameters {
    real<lower=0> l0_alpha; // Initial infection
    real<lower=0> l0_beta;
    real<lower=10> l0;
    real<lower=0> beta0_alpha; // transmission rate at t=0
    real<lower=0> beta0_beta;
    real<lower=0> beta0;
    real<lower=0> beta1_alpha; // transmission rate at t=inf
    real<lower=0> beta1_beta;
    real<lower=0> beta1;
    real<lower=0> b_date; // startig date of intervention
    real<lower=0> theta_alpha; //steepness of intervention logistic curve
    real<lower=0> theta_beta;
    real<lower=0> theta;
    real<lower=0> gamma_alpha; // recovery rate
    real<lower=0> gamma_beta;
    real<lower=0, upper=1> gamma;
    real<lower=0> delta_alpha; // death rate
    real<lower=0> delta_beta;
    real<lower=0, upper=1> delta;
    real<lower=0> q1_alpha; //initial detection rate
    real<lower=0> q1_beta;
    real<lower=0, upper=1> q1;
    real<lower=0> q2_alpha; //final detection rate
    real<lower=0> q2_beta;
    real<lower=0, upper=1> q2;
    real<lower=0> theta_q_alpha; //steepness of logistic curve
    real<lower=0> theta_q_beta;
    real<lower=0> theta_q;
    real<lower=0> q_date;
    real<lower=0> sigma_I; // noise factor for the infection
    real<lower=0> sigma_R; // noise factor for recovery
    real<lower=0> sigma_D; // noise factor for death
    real<lower=0> sigma_I0; // observation error for I0
    real<lower=0> sigma_R0; // observation error for R;
    vector<lower=0>[T] S; // real cummulative infection
    vector<lower=0>[T] I;
    real<lower=0> b;
    real<lower=0, upper=1> q;
    real<lower=0> NI;
    real<lower=0> R;
}  
model {
    l0 ~ gamma(l0_alpha, l0_beta);
    beta0 ~ gamma(beta0_alpha, beta0_beta);
    beta1 ~ gamma(beta1_alpha, beta1_beta);
    b_date ~ uniform(0, T);
    theta ~ gamma(theta_alpha, theta_beta);
    gamma ~ beta(gamma_alpha, gamma_beta);
    delta ~ beta(delta_alpha, delta_beta);
    q1 ~ beta(q1_alpha, q1_beta);
    q2 ~ beta(q2_alpha, q2_beta);
    theta_q ~ gamma(theta_q_alpha, theta_q_beta);
    q_date ~ uniform(0, T);
    S[1] ~ normal(N, 0.1);
    I[1] ~ normal(l0, sigma_I);
    D[1] ~ normal(0, 0.1);
    I0[1] ~ normal(0, 0.1);
    R0[1] ~ normal(0, 0.1);

    for (t in 2:T){
        b ~ normal(beta0 + (beta1 - beta0) * 1 / (1 + exp(-theta * (t - b_date))), 0.01); // transmission rate
        q ~ normal(q1 + (q2 - q1) * 1 / (1 + exp(-theta_q * (t - q_date))), 0.01); // observation rate
        NI ~ normal(b * S[t-1] * I[t-1] / N, sigma_I);
        R ~ normal(gamma * I[t-1], sigma_R);
        D[t] ~ normal(delta * I[t-1], sigma_D);
        S[t]  ~ normal(S[t-1] - NI, 0.1);
        I[t]  ~ normal(I[t-1] + NI - R - D[t], 0.1) ;
        I0[t] ~ normal(q * I[t], sigma_I0);
        R0[t] ~ normal(q * R, sigma_R0);
  }
}