library(tidyr)
library(plyr)
library(MCMCpack)
options(digits = 4)

#Load previously saved data and results
seed = 100
load(paste0("EDiSC_fit_seed_",seed,".RData"))
load("snippets.RData")
load("embeddings.RData")
# load("EDiSC_fit_data.RData")


# PERMUTED ----------------------------------------------------------------

#Extract posterior
EDiSC_posterior = rstan::extract(EDiSC_fit, permuted = TRUE)

#Rearrange posterior samples and variable names into familiar forms
num.samples = (EDiSC_fit@sim$iter - EDiSC_fit@sim$warmup) / EDiSC_fit@sim$thin
burn.in = 0

sense.probs.sim = array(data = aperm(exp(EDiSC_posterior$log_sense_probs), c(1,3,2)),
                        dim = c(num.samples, num.senses, num.snippets),
                        dimnames = list(Sample = 1:num.samples,
                                        Sense  = 1:num.senses, SnippetID = snippets$SnippetID))

psi.tilde.sim = array(data = aperm(EDiSC_posterior$psi_tilde, c(1,4,2,3)),
                      dim = c(num.samples, num.words, num.senses, num.periods),
                      dimnames = list(Sample = 1:num.samples,
                                      Word = 1:num.words, Sense = 1:num.senses, Time = 1:num.periods))

phi.tilde.sim = array(data = aperm(EDiSC_posterior$phi_tilde, c(1,4,2,3)),
                      dim = c(num.samples, num.senses, num.genres, num.periods),
                      dimnames = list(Sample = 1:num.samples, Sense = 1:num.senses,
                                      Genre = 1:num.genres, Time = 1:num.periods))

chi.sim = array(data = aperm(EDiSC_posterior$chi, c(1,3,2)),
                dim = c(num.samples, embedding.dim, num.senses),
                dimnames = list(Sample = 1:num.samples, Embedding = 1:embedding.dim, Sense = 1:num.senses))

theta.sim = array(data = aperm(EDiSC_posterior$theta, c(1,3,2)),
                  dim = c(num.samples, embedding.dim, num.periods),
                  dimnames = list(Sample = 1:num.samples, Embedding = 1:embedding.dim, Time = 1:num.periods))

sigma.sim = array(data = EDiSC_posterior$sigma, dim = c(num.samples, num.words), 
                  dimnames = list(Sample = 1:num.samples, Word = 1:num.words))

#Other diagnostics
Run.Time = lubridate::as.duration(rstan::get_elapsed_time(EDiSC_fit)[2]) #run time for sampling phase only


# UNPERMUTED --------------------------------------------------------------

#Run above code then the following
sense.probs.sim[] = sense.probs.sim[order(EDiSC_fit@sim$permutation[[1]]),,]
psi.tilde.sim[] = psi.tilde.sim[order(EDiSC_fit@sim$permutation[[1]]),,,]
phi.tilde.sim[] = phi.tilde.sim[order(EDiSC_fit@sim$permutation[[1]]),,,]
chi.sim[] = chi.sim[order(EDiSC_fit@sim$permutation[[1]]),,]
theta.sim[] = theta.sim[order(EDiSC_fit@sim$permutation[[1]]),,]
sigma.sim[] = sigma.sim[order(EDiSC_fit@sim$permutation[[1]]),]


# INC WARMUP --------------------------------------------------------------

#Extract posterior
EDiSC_posterior = rstan::extract(EDiSC_fit, permuted = FALSE, inc_warmup = FALSE)[,1,]

#Rearrange posterior samples and variable names into familiar forms
num.samples = (EDiSC_fit@sim$iter - EDiSC_fit@sim$warmup) / EDiSC_fit@sim$thin
burn.in = 0

sense.probs.sim = aperm(array(data = exp(as.matrix(dplyr::select(as.data.frame(EDiSC_posterior), starts_with("log_sense_probs")))),
                              dim = c(num.samples, num.snippets, num.senses),
                              dimnames = list(Sample = 1:num.samples, SnippetID = snippets$SnippetID, Sense  = 1:num.senses)),
                        c(1,3,2))

psi.tilde.sim = aperm(array(data = as.matrix(dplyr::select(as.data.frame(EDiSC_posterior), starts_with("psi_tilde"))),
                            dim = c(num.samples, num.senses, num.periods, num.words),
                            dimnames = list(Sample = 1:num.samples, Sense = 1:num.senses, 
                                            Time = 1:num.periods, Word = 1:num.words)), 
                      c(1,4,2,3))

phi.tilde.sim = aperm(array(data = as.matrix(dplyr::select(as.data.frame(EDiSC_posterior), starts_with("phi_tilde"))),
                            dim = c(num.samples, num.genres, num.periods, num.senses),
                            dimnames = list(Sample = 1:num.samples, Genre = 1:num.genres, 
                                            Time = 1:num.periods, Sense = 1:num.senses)),
                      c(1,4,2,3))

chi.sim = aperm(array(data = as.matrix(dplyr::select(as.data.frame(EDiSC_posterior), starts_with("chi"))),
                      dim = c(num.samples, num.senses, embedding.dim),
                      dimnames = list(Sample = 1:num.samples, Sense = 1:num.senses, Embedding = 1:embedding.dim)), 
                c(1,3,2))

theta.sim = aperm(array(data = as.matrix(dplyr::select(as.data.frame(EDiSC_posterior), starts_with("theta"))),
                        dim = c(num.samples, num.periods, embedding.dim),
                        dimnames = list(Sample = 1:num.samples, Time = 1:num.periods, Embedding = 1:embedding.dim)), 
                  c(1,3,2))

sigma.sim = array(data = as.matrix(dplyr::select(as.data.frame(EDiSC_posterior), starts_with("sigma"))),
                  dim = c(num.samples, num.words),
                  dimnames = list(Sample = 1:num.samples, Word = 1:num.words))

#Other diagnostics
Run.Time = lubridate::as.duration(sum(rstan::get_elapsed_time(EDiSC_fit))) #total run time, including warmup

