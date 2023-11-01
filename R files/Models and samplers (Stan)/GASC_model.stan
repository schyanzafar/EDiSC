// Stan model for GASC

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
  
  real<lower=0> kappa_psi;        // prior hyperparameter
  real<lower=0> kappa_phi_shape;  // prior hyperparameter
  real<lower=0> kappa_phi_rate;   // prior hyperparameter
}

parameters {
  vector[K] phi[G,T];      // sense prevalence parameter
  vector[V] psi[K,T];      // context-word parameter
  real<lower=0> kappa_phi; // phi variance parameter
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
      psi_tilde[k,t] = softmax(psi[k,t]);
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
  // prior for psi
  for (k in 1:K) {
    //psi[k,1] ~ uniform(-2,2);
    for (t in 2:T) {
      psi[k,t] ~ normal(psi[k,t-1], sqrt(2*kappa_psi));
    }
  }
  
  // prior for phi
  for (g in 1:G) {
    //phi[g,1] ~ uniform(-2,2);
    for (t in 2:T) {
      phi[g,t] ~ normal(phi[g,t-1], sqrt(2*kappa_phi));
    }
  }
  
  // prior for kappa_phi
  kappa_phi ~ inv_gamma(kappa_phi_shape, kappa_phi_rate);
  
  // log likelihood
  for (d in 1:D) {
    target += log_sum_exp(log_sense_probs[d]);
  }
}
