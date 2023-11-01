options(digits = 4)
library(rstan)
rstan_options(auto_write = TRUE)

#Set working directory
# setwd("...")
print(getwd())

#Read snippets data
snippets_dir = getwd()
snippets_name = "snippets.RData"
load(paste(snippets_dir, snippets_name, sep = "/"))

#Specify model path
model_dir  = ".../Stan"
model_name = "DiSC_model.stan"
model_path = paste(model_dir, model_name, sep = "/")


#Set number of senses K
num.senses = 3

#Set coefficient alpha in the model X(t) = alpha * X(t-1) + noise
alpha.phi   = 0.9
alpha.theta = 0.9

#Set variances for phi, chi and theta
kappa.phi   = 0.25
kappa.theta = 0.25
kappa.chi   = 1.25

#Prepare snippets to feed into Stan model
stan_snippets = snippets[,1:snippet.length] #discard extra columns
stan_snippets = t(apply(stan_snippets, 1, function(x){c(x[!is.na(x)], x[is.na(x)])})) #move NAs to end of each snippet
stan_snippets = replace(stan_snippets, is.na(stan_snippets), 0) #replace NAs with zero

#Prepare list of data to feed into Stan
stan_data = list(V = num.words,
                 K = num.senses,
                 G = num.genres,
                 T = num.periods,
                 D = num.snippets,
                 L = snippet.length, 
                 
                 snippets = stan_snippets, 
                 lengths  = snippet.lengths, 
                 times    = snippets$Time, 
                 genres   = snippets$genre,
                 
                 alpha_phi   = alpha.phi,
                 alpha_theta = alpha.theta, 
                 kappa_phi   = kappa.phi,
                 kappa_theta = kappa.theta, 
                 kappa_chi   = kappa.chi)

rm(stan_snippets); gc()


#Save data 
save(stan_data, file = "DiSC_fit_data.RData")
save(snippets, snippets.info, words.used, num.words, num.periods, num.genres,
     num.senses, num.snippets, snippet.length, snippet.lengths,  file = "snippets.RData")


#Fit model using Stan NUTS
seeds = c(100, 200)
target_accept_rate = 0.6
num_iterations = 2000
warmup = 1000
thin = 1

for (seed in seeds) {
  DiSC_fit = stan(file = model_path, model_name = model_name, data = stan_data, verbose = FALSE,
                  pars = c("log_sense_probs", "phi_tilde", "psi_tilde", "theta", "chi"),
                  control = list(metric = "diag_e", max_treedepth = 10, stepsize_jitter = 1,
                                 # adapt_init_buffer = 75, adapt_term_buffer = 50, adapt_window = 25, 
                                 adapt_delta = target_accept_rate), 
                  cores = 1, chains = 1, iter = num_iterations, warmup = warmup, thin = thin, seed = seed)
  
  #Save results
  save(DiSC_fit, file = paste0("DiSC_fit_seed_",seed,".RData"))
}


#Fit model using Stan variational inference
seeds = c(100, 200)
num_samples = 1000
max_iterations = 10^5

compiled_model = stan_model(model_path, model_name)

for (seed in seeds) {
  DiSC_fit = vb(object = compiled_model, data = stan_data, algorithm = "meanfield",
                pars = c("log_sense_probs", "phi_tilde", "psi_tilde", "theta", "chi"),
                seed = seed, output_samples = num_samples, iter = max_iterations,
                tol_rel_obj = 0.0005, eval_elbo = 1000, elbo_samples = 1000, grad_samples = 1)

  #Save results
  save(DiSC_fit, file = paste0("DiSC_fit_seed_",seed,".RData"))
}
