library(tidyr)
library(plyr)
library(MCMCpack)
options(digits = 4)

#Load previously saved data and results
seed = 100
load(paste0("DiSC_fit_seed_",seed,".RData"))
load("snippets.RData")
# load("DiSC_fit_data.RData")


# PERMUTED ----------------------------------------------------------------

#Extract posterior
DiSC_posterior = rstan::extract(DiSC_fit, permuted = TRUE)

#Rearrange posterior samples and variable names into familiar forms
num.samples = (DiSC_fit@sim$iter - DiSC_fit@sim$warmup) / DiSC_fit@sim$thin
burn.in = 0

sense.probs.sim = array(data = aperm(exp(DiSC_posterior$log_sense_probs), c(1,3,2)),
                        dim = c(num.samples, num.senses, num.snippets),
                        dimnames = list(Sample = 1:num.samples,
                                        Sense  = 1:num.senses, SnippetID = snippets$SnippetID))

psi.tilde.sim = array(data = aperm(DiSC_posterior$psi_tilde, c(1,4,2,3)),
                      dim = c(num.samples, num.words, num.senses, num.periods),
                      dimnames = list(Sample = 1:num.samples,
                                      Word = 1:num.words, Sense = 1:num.senses, Time = 1:num.periods))

phi.tilde.sim = array(data = aperm(DiSC_posterior$phi_tilde, c(1,4,2,3)),
                      dim = c(num.samples, num.senses, num.genres, num.periods),
                      dimnames = list(Sample = 1:num.samples, Sense = 1:num.senses,
                                      Genre = 1:num.genres, Time = 1:num.periods))

chi.sim = array(data = aperm(DiSC_posterior$chi, c(1,3,2)),
                dim = c(num.samples, num.words, num.senses),
                dimnames = list(Sample = 1:num.samples, Word = 1:num.words, Sense = 1:num.senses))

theta.sim = array(data = aperm(DiSC_posterior$theta, c(1,3,2)),
                  dim = c(num.samples, num.words, num.periods),
                  dimnames = list(Sample = 1:num.samples, Word = 1:num.words, Time = 1:num.periods))

#Other diagnostics
Run.Time = lubridate::as.duration(rstan::get_elapsed_time(DiSC_fit)[2]) #run time for sampling phase only


# UNPERMUTED --------------------------------------------------------------

#Run above code then the following
sense.probs.sim[] = sense.probs.sim[order(DiSC_fit@sim$permutation[[1]]),,]
psi.tilde.sim[] = psi.tilde.sim[order(DiSC_fit@sim$permutation[[1]]),,,]
phi.tilde.sim[] = phi.tilde.sim[order(DiSC_fit@sim$permutation[[1]]),,,]
chi.sim[] = chi.sim[order(DiSC_fit@sim$permutation[[1]]),,]
theta.sim[] = theta.sim[order(DiSC_fit@sim$permutation[[1]]),,]


# INC WARMUP --------------------------------------------------------------

#Extract posterior
DiSC_posterior = rstan::extract(DiSC_fit, permuted = FALSE, inc_warmup = TRUE)[,1,]

#Rearrange posterior samples and variable names into familiar forms
num.samples = (DiSC_fit@sim$iter - DiSC_fit@sim$warmup) / DiSC_fit@sim$thin
burn.in = 0

sense.probs.sim = aperm(array(data = exp(as.matrix(dplyr::select(as.data.frame(DiSC_posterior), starts_with("log_sense_probs")))),
                              dim = c(num.samples, num.snippets, num.senses),
                              dimnames = list(Sample = 1:num.samples, SnippetID = snippets$SnippetID, Sense  = 1:num.senses)),
                        c(1,3,2))

psi.tilde.sim = aperm(array(data = as.matrix(dplyr::select(as.data.frame(DiSC_posterior), starts_with("psi_tilde"))),
                            dim = c(num.samples, num.senses, num.periods, num.words),
                            dimnames = list(Sample = 1:num.samples, Sense = 1:num.senses, 
                                            Time = 1:num.periods, Word = 1:num.words)), 
                      c(1,4,2,3))

phi.tilde.sim = aperm(array(data = as.matrix(dplyr::select(as.data.frame(DiSC_posterior), starts_with("phi_tilde"))),
                            dim = c(num.samples, num.genres, num.periods, num.senses),
                            dimnames = list(Sample = 1:num.samples, Genre = 1:num.genres, 
                                            Time = 1:num.periods, Sense = 1:num.senses)),
                      c(1,4,2,3))

chi.sim = aperm(array(data = as.matrix(dplyr::select(as.data.frame(DiSC_posterior), starts_with("chi"))),
                      dim = c(num.samples, num.senses, num.words),
                      dimnames = list(Sample = 1:num.samples, Sense = 1:num.senses, Word = 1:num.words)), 
                c(1,3,2))

theta.sim = aperm(array(data = as.matrix(dplyr::select(as.data.frame(DiSC_posterior), starts_with("theta"))),
                        dim = c(num.samples, num.periods, num.words),
                        dimnames = list(Sample = 1:num.samples, Time = 1:num.periods, Word = 1:num.words)), 
                  c(1,3,2))

#Other diagnostics
Run.Time = lubridate::as.duration(sum(rstan::get_elapsed_time(DiSC_fit))) #total run time, including warmup

