library(tidyr)
library(plyr)
library(MCMCpack)
options(digits = 4)

#Load previously saved data and results
seed = 100
load(paste0("GASC_fit_seed_",seed,".RData"))
load("snippets.RData")
# load("GASC_fit_data.RData")


# PERMUTED ----------------------------------------------------------------

#Extract posterior
GASC_posterior = rstan::extract(GASC_fit, permuted = TRUE)

#Rearrange posterior samples and variable names into familiar forms
num.samples = (GASC_fit@sim$iter - GASC_fit@sim$warmup) / GASC_fit@sim$thin
burn.in = 0

sense.probs.sim = array(data = aperm(exp(GASC_posterior$log_sense_probs), c(1,3,2)),
                        dim = c(num.samples, num.senses, num.snippets),
                        dimnames = list(Sample = 1:num.samples,
                                        Sense  = 1:num.senses, SnippetID = snippets$SnippetID))

psi.tilde.sim = array(data = aperm(GASC_posterior$psi_tilde, c(1,4,2,3)),
                      dim = c(num.samples, num.words, num.senses, num.periods),
                      dimnames = list(Sample = 1:num.samples,
                                      Word = 1:num.words, Sense = 1:num.senses, Time = 1:num.periods))

phi.tilde.sim = array(data = aperm(GASC_posterior$phi_tilde, c(1,4,2,3)),
                      dim = c(num.samples, num.senses, num.genres, num.periods),
                      dimnames = list(Sample = 1:num.samples, Sense = 1:num.senses,
                                      Genre = 1:num.genres, Time = 1:num.periods))

kappa.phi.sim = as.vector(GASC_posterior$kappa_phi)

#Other diagnostics
Run.Time = lubridate::as.duration(rstan::get_elapsed_time(GASC_fit)[2]) #run time for sampling phase only


# UNPERMUTED --------------------------------------------------------------

#Run above code then the following
sense.probs.sim[] = sense.probs.sim[order(GASC_fit@sim$permutation[[1]]),,]
psi.tilde.sim[] = psi.tilde.sim[order(GASC_fit@sim$permutation[[1]]),,,]
phi.tilde.sim[] = phi.tilde.sim[order(GASC_fit@sim$permutation[[1]]),,,]
kappa.phi.sim[] = kappa.phi.sim[order(GASC_fit@sim$permutation[[1]])]


# INC WARMUP --------------------------------------------------------------

#Extract posterior
GASC_posterior = rstan::extract(GASC_fit, permuted = FALSE, inc_warmup = TRUE)[,1,]

#Rearrange posterior samples and variable names into familiar forms
num.samples = (GASC_fit@sim$iter - GASC_fit@sim$warmup) / GASC_fit@sim$thin
burn.in = 0

sense.probs.sim = aperm(array(data = exp(as.matrix(dplyr::select(as.data.frame(GASC_posterior), starts_with("log_sense_probs")))),
                              dim = c(num.samples, num.snippets, num.senses),
                              dimnames = list(Sample = 1:num.samples, SnippetID = snippets$SnippetID, Sense  = 1:num.senses)),
                        c(1,3,2))

psi.tilde.sim = aperm(array(data = as.matrix(dplyr::select(as.data.frame(GASC_posterior), starts_with("psi_tilde"))),
                            dim = c(num.samples, num.senses, num.periods, num.words),
                            dimnames = list(Sample = 1:num.samples, Sense = 1:num.senses, 
                                            Time = 1:num.periods, Word = 1:num.words)), 
                      c(1,4,2,3))

phi.tilde.sim = aperm(array(data = as.matrix(dplyr::select(as.data.frame(GASC_posterior), starts_with("phi_tilde"))),
                            dim = c(num.samples, num.genres, num.periods, num.senses),
                            dimnames = list(Sample = 1:num.samples, Genre = 1:num.genres, 
                                            Time = 1:num.periods, Sense = 1:num.senses)),
                      c(1,4,2,3))

kappa.phi.sim = unlist(dplyr::select(as.data.frame(GASC_posterior), starts_with("kappa_phi")))

#Other diagnostics
Run.Time = lubridate::as.duration(sum(rstan::get_elapsed_time(GASC_fit))) #total run time, including warmup

