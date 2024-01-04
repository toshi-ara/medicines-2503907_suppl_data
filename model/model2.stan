// model2.stan

// function to avoid 0 and 1
functions{
  real adjust_p(real x){
    return min([max([x, 0.00001]), 0.99999]);
  }
}

data {
  int<lower=1> N;  // number of all data (all drugs, animals, and time points)
  int<lower=1> ND; // number of drugs
  int<lower=1> NG; // number of animals (groups)
  int GROUP[N];    // ID of each animal (group)
  int DRUG1[N];    // ID of each drug (integer)
  int DRUG2[N];    // without Adr: 0, with Adr: 1
  real TIME[N];    // time point
  int<lower=0,upper=6> SCORE[N]; // score by experiments
}

parameters {
  // parameters of each drug (fixed effect)
  real mu0[ND];
  real log_sigma0[ND];
  real<lower=0,upper=1> adr;

  // standard deviation of parameter in each drug (random effect)
  real<lower=0> s_mu0[ND];
  real<lower=0> log_s_sigma0[ND];

  // offset value for each individual
  real d[NG];

  // parameters of each drug and animal
  real mu[ND, NG];
  real<lower=0> sigma[ND, NG];
}

transformed parameters {
  real<lower=0.00001,upper=0.99999> p[N];
  real X;
  real Mean;
  real Std;

  for (i in 1:N) {
    X = 100 - (1 - adr * DRUG2[i]) * TIME[i];
    Mean = mu[DRUG1[i], GROUP[i]];
    Std = sigma[DRUG1[i], GROUP[i]];
    p[i] = adjust_p(1.0 - Phi_approx((X - Mean) / Std));
  }
}

model {
  // parameter of each drug (each animal)
  // with offset value for each animal
  for (i in 1:ND) {
    for (j in 1:NG) {
      mu[i, j] ~ normal(mu0[i] + d[j], s_mu0[i]);
      sigma[i, j] ~ lognormal(log_sigma0[i], log_s_sigma0[i]);
    }
  }

  mu0 ~ cauchy(50, 20);
  s_mu0 ~ cauchy(0, 1);  // halfCauchy
  log_sigma0 ~ normal(2.5, 1);
  d ~ normal(0, 20);

  SCORE ~ binomial(6, p);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = binomial_lpmf(SCORE[n] | 6, p[n]);
  }
}
