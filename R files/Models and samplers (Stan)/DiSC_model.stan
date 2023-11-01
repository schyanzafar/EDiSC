// Stan model for DiSC

data {
  int<lower=2> V;  // num words
  int<lower=2> K;  // num senses
  int<lower=1> G;  // num genres
  int<lower=1> T;  // num time periods
  int<lower=1> D;  // num snippets
  int<lower=1> L;  // snippet length
  
  int<lower=0, upper=V> snippets[D,L];  // snippets, i.e. words. zero indicates position unoccupied
  int<lower=0, upper=L> lengths[D];     // snippet lengths
  int<lower=1, upper=T> times[D];       // snippet time periods
  int<lower=1, upper=G> genres[D];      // snippet genres
  
  real<lower=0, upper=1> alpha_phi;   // prior hyperparameter
  real<lower=0, upper=1> alpha_theta; // prior hyperparameter
  real<lower=0> kappa_phi;            // prior hyperparameter
  real<lower=0> kappa_theta;          // prior hyperparameter
  real<lower=0> kappa_chi;            // prior hyperparameter
}

transformed data {
  real phi_stat_sd;    // std dev of the phi AR(1) stationary distribution
  real theta_stat_sd;  // std dev of the theta AR(1) stationary distribution
  
  phi_stat_sd = sqrt(kappa_phi / (1 - alpha_phi^2));
  theta_stat_sd = sqrt(kappa_theta / (1 - alpha_theta^2));
}

parameters {
  vector[K] phi[G,T];  // sense prevalence parameter
  vector[V] theta[T];  // context-word-time parameter
  vector[V] chi[K];    // context-word-sense parameter
}

transformed parameters {
  simplex[K] phi_tilde[G,T];     // sense prevalence distribution
  simplex[V] psi_tilde[K,T];     // context-word distribution
  vector[K] log_sense_probs[D];  // unnormalised log probabilities for each snippet belonging to each sense

  for (g in 1:G) {
    for (t in 1:T) {
      phi_tilde[g,t] = softmax(phi[g,t]);
    }
  }

  for (k in 1:K) {
    for (t in 1:T) {
      psi_tilde[k,t] = softmax(theta[t] + chi[k]);
    }
  }

  for (d in 1:D) {
    for (k in 1:K) {
      log_sense_probs[d,k] = log(phi_tilde[genres[d],times[d],k])
                             + sum(log(psi_tilde[k,times[d],snippets[d,1:lengths[d]]]));
    } // for k
  } // for d
}

model {
  // prior for chi
  for (k in 1:K) {
    chi[k] ~ normal(0, sqrt(kappa_chi));
  }
  
  // prior for theta
  theta[1] ~ normal(0, theta_stat_sd);
  for (t in 2:T) {
    theta[t] ~ normal(alpha_theta * theta[t-1], sqrt(kappa_theta)); 
  }
  
  // prior for phi
  for (g in 1:G) {
    phi[g,1] ~ normal(0, phi_stat_sd);
    for (t in 2:T) {
      phi[g,t] ~ normal(alpha_phi * phi[g,t-1], sqrt(kappa_phi));
    }
  }
  
  // log likelihood
  for (d in 1:D) {
    target += log_sum_exp(log_sense_probs[d]);
  }
}
