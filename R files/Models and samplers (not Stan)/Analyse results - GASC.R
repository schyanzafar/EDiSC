#############################################################################
# Required libraries ------------------------------------------------------ #
#############################################################################

library(MCMCpack)

options(digits = 4)
#options(scipen = 10)


#############################################################################
# Load previously saved results ------------------------------------------- #
#############################################################################

load("snippets.RData")
load("other.parameters.RData")
# load("initial.z.RData")

b = 1
filename.suffix = paste((b-1)*num.iterations+1, b*num.iterations, sep = "-")

load(paste("psi.tilde.sim", filename.suffix, "RData", sep = "."))
load(paste("phi.tilde.sim", filename.suffix, "RData", sep = "."))
load(paste("kappa.phi.sim", filename.suffix, "RData", sep = "."))
load(paste("sense.probs.sim", filename.suffix, "RData", sep = "."))
load(paste("run.time", filename.suffix, "RData", sep = "."))
load(paste("HMC.parameters.RData", filename.suffix, "RData", sep = "."))



#############################################################################
# Analyse MCMC results ---------------------------------------------------- #
#############################################################################

#Run time
print(Run.Time)

#Acceptance rate
accept.count.phi / num.iterations
accept.count.psi / num.iterations


#Likelihood
log.lik.sim = colSums(log(apply(sense.probs.sim, 1, colSums)))
plot(log.lik.sim, type = "l")
mean(log.lik.sim)

num.samples = num.iterations/N
burn.in = 100
#burn.in = 0
plot(log.lik.sim[(burn.in+1):num.samples], type = "l")
mean(log.lik.sim[(burn.in+1):num.samples])


#Trace plots - phi
Time = 8
Genre = 1
plot(as.mcmc(phi.tilde.sim[(burn.in+1):num.samples,,Genre,Time]))


#Trace plots - psi
Time = 8
Sense = 1
Words = 1:3
plot(as.mcmc(psi.tilde.sim[(burn.in+1):num.samples,Words,Sense,Time]))


#Trace plots - kappa.phi
plot(as.mcmc(kappa.phi.sim))
