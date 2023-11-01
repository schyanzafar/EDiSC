#############################################################################
# Required libraries ------------------------------------------------------ #
#############################################################################

library(tidyr)
library(plyr)
library(MCMCpack)
library(Matrix)

options(digits = 4)
#options(scipen = 10)



#############################################################################
# Specify locations for reading data and writing results ------------------ #
#############################################################################

#Set working directory
# setwd("...")

#Read snippets data
snippets_dir = getwd()
snippets_name = "snippets.RData"
load(paste(snippets_dir, snippets_name, sep = "/"))



#############################################################################
# Set model parameters ---------------------------------------------------- #
#############################################################################

#Set number of senses K
num.senses = 4

#Set coefficient alpha in the model X(t) = alpha * X(t-1) + noise
alpha.phi   = 0.9
alpha.theta = 0.9

#Set variances for phi, chi and theta
kappa.phi   = 0.25
kappa.theta = 0.25
kappa.chi   = 1.25



#############################################################################
# Set MCMC parameters ----------------------------------------------------- #
#############################################################################

#Set seed for initialisation
random.seed = 100

#Specify number of MCMC iterations
num.iterations = 10000

#Sample every Nth iteration
N = 5

#Specify number of times to run the above MCMC iterations
num.runs = 1

#Set tempering parameters
#we will target prior x likelihood ^ inv.temp, where inv.temp increases from 
#the minimum specified to 1 in the given number of MCMC iterations
temper.phi   = TRUE
temper.chi   = TRUE
temper.theta = FALSE
min.inv.temp = 0.1
num.iter.temp = num.iterations/2

#Stop tuning HMC scale parameters after this many MCMC iterations
stop.tuning.after = num.iterations/2 #num.iterations * num.runs

#Set number of leapfrog steps for HMC proposals
num.leapfrog.steps.phi = 5
num.leapfrog.steps.chi = 10
num.leapfrog.steps.theta = 5

#Set HMC tuning parameters - refer to Shaby and Wells (2010)
target.accept.rate.phi = 0.651
target.accept.rate.chi = 0.651
target.accept.rate.theta = 0.651
block.size = 10 #k in the paper

#c0.phi = 1
c1.phi = 0.8
#c0.chi = 1
c1.chi = 0.8
#c0.theta = 1
c1.theta = 0.8

#Set initial scales for HMC proposals - these are the epsilon^2 parameters
phi.prop.scale.init = 2.4^2 / num.senses / num.leapfrog.steps.phi
chi.prop.scale.init = 2.4^2 / (num.words * num.senses) / num.leapfrog.steps.chi
theta.prop.scale.init = 2.4^2 / num.words / num.leapfrog.steps.theta



#############################################################################
# Functions for parameter updates ----------------------------------------- #
#############################################################################

#Function that takes as input a single snippet, phi.g.t and psi.t, 
#and returns a vector of (unnormalised) sense probabilities for that snippet
get.sense.prob = function(snippet, phi.g.t.tilde, psi.t.tilde) {
  apply(psi.t.tilde[snippet,], 2, prod, na.rm = TRUE) * phi.g.t.tilde
}#get.sense.prob


#Function to calculate the log posterior gradient wrt phi for a single time slice and genre
get.phi.g.t.gradient = function(phi.g.t, phi.g.t.tilde, phi.g.t.prior.mean, phi.t.prior.var, 
                                snippet.counts.g.t, sense.probs.g.t, temper = FALSE, inv.temp = 1) {
  
  - (phi.g.t - phi.g.t.prior.mean) / phi.t.prior.var + ifelse(temper, inv.temp, 1) * 
    (colSums(t(sense.probs.g.t) / colSums(sense.probs.g.t)) - snippet.counts.g.t * phi.g.t.tilde)
  
}#get.phi.g.t.gradient


#Function to calculate the log posterior gradient wrt chi for all senses
get.chi.gradient = function(chi, psi.tilde, chi.prior.mean, chi.prior.var, sense.probs, word.counts.in.snippets, 
                            snippet.lengths, temper = FALSE, inv.temp = 1) {
  
  P = t(sense.probs)/colSums(sense.probs)
  
  I.P.psi.tilde = sapply(1:num.senses, function(l){psi.tilde[,l,snippets$Time] %*% (snippet.lengths*P[,l])})
  
  - (chi - chi.prior.mean) / chi.prior.var + 
    ifelse(temper, inv.temp, 1) * (as.matrix(word.counts.in.snippets %*% P) - I.P.psi.tilde)
  
}#get.chi.gradient


#Function to calculate the log posterior gradient wrt theta for a single time slice
get.theta.t.gradient = function(theta.t, psi.t.tilde, theta.t.prior.mean, theta.t.prior.var, sense.probs.t, 
                                word.counts.in.snippets.t, snippet.lengths.t, temper = FALSE, inv.temp = 1) {

  as.numeric(
    - (theta.t - theta.t.prior.mean) / theta.t.prior.var + ifelse(temper, inv.temp, 1) * 
      (rowSums(word.counts.in.snippets.t - t(snippet.lengths.t * ((t(sense.probs.t) / colSums(sense.probs.t)) %*% t(psi.t.tilde)))) )
  )
  
}#get.theta.t.gradient


#Function to propose phi for genre g time t using HMC
#Function adapted from Neal, R.M., 2011. MCMC using Hamiltonian dynamics
propose.new.phi.g.t.HMC = function(num.steps, epsilon, phi.g.t, phi.g.t.tilde, psi.t.tilde, phi.g.t.prior.mean, 
                                   phi.t.prior.var, sense.probs.g.t, snippet.counts.g.t, temper = FALSE, inv.temp = 1){
  
  #Generate initial momentum vector
  p = rnorm(num.senses)
  
  #Evaluate potential and kinetic energies at start of trajectory
  current_U = -(sum(dnorm(phi.g.t, phi.g.t.prior.mean, sqrt(phi.t.prior.var), log = TRUE)) + 
                  ifelse(temper, inv.temp, 1) * sum(log(colSums(sense.probs.g.t))))
  current_K = sum(p^2) / 2
  
  #Make a half step for momentum at the beginning
  grad_U = -get.phi.g.t.gradient(phi.g.t, phi.g.t.tilde, phi.g.t.prior.mean, phi.t.prior.var, 
                                 snippet.counts.g.t, sense.probs.g.t, temper, inv.temp)
  p = p - epsilon * grad_U / 2
  
  #Alternate full steps for position and momentum
  for (s in 1:num.steps) {
    
    #Make a full step for the position
    phi.g.t = phi.g.t + epsilon * p
    exp.phi.g.t = exp(phi.g.t)
    phi.g.t.tilde = exp.phi.g.t / sum(exp.phi.g.t)
    
    sense.probs.g.t[] = apply(snippets[SnippetIDs[[t]][[g]], 1:snippet.length], 1, get.sense.prob, 
                              phi.g.t.tilde, psi.t.tilde)
    
    grad_U = -get.phi.g.t.gradient(phi.g.t, phi.g.t.tilde, phi.g.t.prior.mean, phi.t.prior.var, 
                                   snippet.counts.g.t, sense.probs.g.t, temper, inv.temp)
    
    #Make a full step for the momentum, except at end of trajectory
    if (s != num.steps) p = p - epsilon * grad_U
    
  }#for s
  
  #Make a half step for momentum at the end
  p = p - epsilon * grad_U / 2
  
  #Evaluate potential and kinetic energies at end of trajectory
  proposed_U = -(sum(dnorm(phi.g.t, phi.g.t.prior.mean, sqrt(phi.t.prior.var), log = TRUE)) + 
                   ifelse(temper, inv.temp, 1) * sum(log(colSums(sense.probs.g.t))))
  proposed_K = sum(p^2) / 2
  
  #Compute log MH ratio given current and proposed energies
  MH = current_U - proposed_U + current_K - proposed_K
  
  
  out = list(MH = MH, phi.g.t = phi.g.t, phi.g.t.tilde = phi.g.t.tilde, sense.probs.g.t = sense.probs.g.t)
  
  return(out)
  
}#propose.new.phi.g.t.HMC


#Function to propose chi for all senses using HMC
#Function adapted from Neal, R.M., 2011. MCMC using Hamiltonian dynamics
propose.new.chi.HMC = function(num.steps, epsilon, chi, psi.tilde, chi.prior.mean, chi.prior.var, sense.probs, 
                               word.counts.in.snippets, snippet.lengths, temper = FALSE, inv.temp = 1){

  #Generate initial momentum vector
  p = rnorm(num.words*num.senses)
  
  #Evaluate potential and kinetic energies at start of trajectory
  current_U = -(sum(dnorm(chi, chi.prior.mean, sqrt(chi.prior.var), log = TRUE)) + 
                  ifelse(temper, inv.temp, 1) * sum(log(colSums(sense.probs))))
  current_K = sum(p^2) / 2
  
  #Make a half step for momentum at the beginning
  grad_U = -get.chi.gradient(chi, psi.tilde, chi.prior.mean, chi.prior.var, sense.probs,
                             word.counts.in.snippets, snippet.lengths, temper, inv.temp)
  p = p - epsilon * grad_U / 2
  
  #Alternate full steps for position and momentum
  for (s in 1:num.steps) {
    
    #Make a full step for the position
    chi = chi + epsilon * p
    psi = sapply(1:num.periods, function(t) {sapply(1:num.senses, function(k) {chi[,k]+theta[,t]} )}, simplify = "array")
    psi.tilde = aperm(aaply(exp(psi), 3, function(x) t(t(x)/colSums(x))), c(2,3,1))
    
    for(t in 1:num.periods) {
      for(g in 1:num.genres) {
        sense.probs[,SnippetIDs[[t]][[g]]] = apply(snippets[SnippetIDs[[t]][[g]], 1:snippet.length], 1, get.sense.prob,
                                                   phi.tilde[,g,t], psi.tilde[,,t])
      }#for g
    }#for t

    grad_U = -get.chi.gradient(chi, psi.tilde, chi.prior.mean, chi.prior.var, sense.probs,
                               word.counts.in.snippets, snippet.lengths, temper, inv.temp)
    
    #Make a full step for the momentum, except at end of trajectory
    if (s != num.steps) p = p - epsilon * grad_U
    
  }#for s
  
  #Make a half step for momentum at the end
  p = p - epsilon * grad_U / 2
  
  #Evaluate potential and kinetic energies at end of trajectory
  proposed_U = -(sum(dnorm(chi, chi.prior.mean, sqrt(chi.prior.var), log = TRUE)) + 
                   ifelse(temper, inv.temp, 1) * sum(log(colSums(sense.probs))))
  proposed_K = sum(p^2) / 2
  
  #Compute log MH ratio given current and proposed energies
  MH = current_U - proposed_U + current_K - proposed_K
  
  
  out = list(MH = MH, chi = chi, psi.tilde = psi.tilde, sense.probs = sense.probs)
  
  return(out)
  
}#propose.new.chi.HMC


#Function to propose theta for time t using HMC
#Function adapted from Neal, R.M., 2011. MCMC using Hamiltonian dynamics
propose.new.theta.t.HMC = function(num.steps, epsilon, theta.t, psi.t.tilde, theta.t.prior.mean, theta.t.prior.var, 
                                   sense.probs.t, word.counts.in.snippets.t, snippet.lengths.t, temper = FALSE, inv.temp = 1){
  
  #Generate initial momentum vector
  p = rnorm(num.words)
  
  #Evaluate potential and kinetic energies at start of trajectory
  current_U = -(sum(dnorm(theta.t, theta.t.prior.mean, sqrt(theta.t.prior.var), log = TRUE)) + 
                  ifelse(temper, inv.temp, 1) * sum(log(colSums(sense.probs.t))))
  current_K = sum(p^2) / 2

  #Make a half step for momentum at the beginning
  grad_U = -get.theta.t.gradient(theta.t, psi.t.tilde, theta.t.prior.mean, theta.t.prior.var, sense.probs.t, 
                                 word.counts.in.snippets.t, snippet.lengths.t, temper, inv.temp)
  p = p - epsilon * grad_U / 2
  
  #Alternate full steps for position and momentum
  for (s in 1:num.steps) {
    
    #Make a full step for the position
    theta.t = theta.t + epsilon * p
    psi.t = theta.t + chi
    exp.psi.t = exp(psi.t)
    psi.t.tilde = t(t(exp.psi.t) / colSums(exp.psi.t))
    
    for(g in 1:num.genres) {
      sense.probs.t[,paste(SnippetIDs[[t]][[g]])] = apply(snippets[SnippetIDs[[t]][[g]], 1:snippet.length], 1,
                                                          get.sense.prob, phi.tilde[,g,t], psi.t.tilde)
    }#for g
    
    grad_U = -get.theta.t.gradient(theta.t, psi.t.tilde, theta.t.prior.mean, theta.t.prior.var, sense.probs.t, 
                                   word.counts.in.snippets.t, snippet.lengths.t, temper, inv.temp)
    
    #Make a full step for the momentum, except at end of trajectory
    if (s != num.steps) p = p - epsilon * grad_U
    
  }#for s
  
  #Make a half step for momentum at the end
  p = p - epsilon * grad_U / 2

  #Evaluate potential and kinetic energies at end of trajectory
  proposed_U = -(sum(dnorm(theta.t, theta.t.prior.mean, sqrt(theta.t.prior.var), log = TRUE)) + 
                   ifelse(temper, inv.temp, 1) * sum(log(colSums(sense.probs.t))))
  proposed_K = sum(p^2) / 2
  
  #Compute log MH ratio given current and proposed energies
  MH = current_U - proposed_U + current_K - proposed_K
  
  
  out = list(MH = MH, theta.t = theta.t, psi.t.tilde = psi.t.tilde, sense.probs.t = sense.probs.t)
  
  return(out)
  
}#propose.new.theta.t.HMC



#############################################################################
# Functions for tuning HMC parameters and tempering ----------------------- #
#############################################################################

#Refer to Shaby and Wells (2010) "Exploring an Adaptive Metropolis Algorithm"
update.scale = function(scale, accept.rate, gamma2, target.accept.rate) {
  
  scale * exp(gamma2 * (accept.rate - target.accept.rate))
  
}#update.scale

#Compute the inverse temperature for simulated tempering of the likelihood
get.inv.temp = function(i, min.inv.temp, num.iter.temp) {
  
  #linearly increasing
  # ifelse(i < num.iter.temp, min.inv.temp + (1 - min.inv.temp) * i / num.iter.temp, 1)
  
  #increasing at a decreasing rate
  ifelse(i < num.iter.temp, min.inv.temp + (1 - min.inv.temp) * (i / num.iter.temp) ^ (1/3), 1)
  
}#get.inv.temp



#############################################################################
# Initialisation ---------------------------------------------------------- #
#############################################################################

#Get the counts of snippets in each time period and genre
snippet.counts = table(factor(snippets$genre, levels = 1:num.genres),
                       factor(snippets$Time, levels = 1:num.periods),
                       dnn = c("Genre", "Time"))


#Initialise parameters AT RANDOM

#Initialise z randomly
set.seed(random.seed)
z = sample(1:num.senses, num.snippets, TRUE)
# z = as.numeric(factor(snippets.info$sense.id)) #initialise z at truth

#Initialise phi based on the initial z
sense.counts = table(factor(z, levels = 1:num.senses), factor(snippets$genre, levels = 1:num.genres),
                     factor(snippets$Time, levels = 1:num.periods), dnn = c("Sense", "Genre", "Time"))

phi.tilde = aperm(sapply(1:num.senses, function(k) {(sense.counts[k,,] + 0.01) / (snippet.counts + 0.01*num.senses)}, simplify = "array"), c(3,1,2))
dimnames(phi.tilde) = list(Sense = 1:num.senses, Genre = 1:num.genres, Time = 1:num.periods)

phi = sapply(1:num.periods, function(t) {t(log(t(phi.tilde[,,t]) / phi.tilde[num.senses,,t]))}, simplify = "array")
dimnames(phi) = dimnames(phi.tilde)

#Adjust phi for identifiability
phi = aperm(sapply(1:num.genres, function(g) {t(t(phi[,g,]) - colMeans(phi[,g,]))}, simplify = "array"), c(1,3,2))
dimnames(phi) = dimnames(phi.tilde)


#Initialise psi based on the initial z
snippets.expanded = gather(cbind(snippets, sense = z), key = position, value = word, 1:snippet.length)
word.counts = table(factor(snippets.expanded$word, levels = 1:num.words),
                    factor(snippets.expanded$sense, levels = 1:num.senses),
                    factor(snippets.expanded$Time, levels = 1:num.periods),
                    dnn = c("Word", "Sense", "Time"))

psi.tilde = sapply(1:num.periods, function(t) {sapply(1:num.senses, function(k)
{(word.counts[,k,t]+0.01)/(sum(word.counts[,k,t])+0.01*num.words)}, simplify = "aray" )}, simplify = "array" )
dimnames(psi.tilde) = list(Word = 1:num.words, Sense = 1:num.senses, Time = 1:num.periods)

psi = sapply(1:num.periods, function(t) {t(log(t(psi.tilde[,,t]) / psi.tilde[num.words,,t]))}, simplify = "array")
dimnames(psi) = dimnames(psi.tilde)

rm(snippets.expanded, word.counts, sense.counts); gc()

#Fit a linear regression to estimate theta and chi. (We have VKT equations and V(K+T) variables)
vkt.combos = expand.grid(dimnames(psi)) #word-sense-time combinations
vkt.combos = data.frame(apply(vkt.combos,2,as.numeric)) #convert to numeric
design.matrix.idx = data.frame(row.idx = rep(1:(num.words*num.senses*num.periods),2),
                               col.idx = c((vkt.combos$Sense-1)*num.words+vkt.combos$Word,
                                           num.words*num.senses + (vkt.combos$Time-1)*num.words+vkt.combos$Word))
design.matrix = sparseMatrix(i = design.matrix.idx$row.idx, j = design.matrix.idx$col.idx, x = 1)
design.matrix = design.matrix[,-((num.words*(num.senses-1)+1):(num.words*num.senses))] #remove last sense of chi for identifiability
#design.matrix = cbind(1, design.matrix) #add intercept column
psi.fit = MatrixModels:::lm.fit.sparse(design.matrix, as.vector(psi), method = "cholesky")

chi = array(data = c(psi.fit$coef[1:(num.words*(num.senses-1))], numeric(num.words)), 
            dim = c(num.words, num.senses), dimnames = list(Word = 1:num.words, Sense = 1:num.senses))

theta = array(data = psi.fit$coef[-(1:(num.words*(num.senses-1)))], 
              dim = c(num.words, num.periods),dimnames = list(Word = 1:num.words, Time = 1:num.periods))

#Adjust chi and theta for identifiability
chi.mean = rowMeans(chi)
theta = theta + chi.mean
chi = chi - chi.mean
theta = t(t(theta) - colMeans(theta))
chi = t(t(chi) - colMeans(chi))

#Re-calculate psi and psi.tilde based on initial theta and chi
psi = sapply(1:num.periods, function(t) {sapply(1:num.senses, function(k) {chi[,k]+theta[,t]} )}, simplify = "array")
dimnames(psi) = list(Word = 1:num.words, Sense = 1:num.senses, Time = 1:num.periods)
psi.tilde = aperm(aaply(exp(psi), 3, function(x) t(t(x)/colSums(x))), c(2,3,1))

rm(vkt.combos, design.matrix.idx, design.matrix, psi, psi.fit, chi.mean); gc()


#Set up structures for prior means (to avoid errors later for phi when num.genres = 1, and for better efficiency)
phi.t.prior.mean = array(data = 0, dim = c(num.senses, num.genres), dimnames = list(Sense = 1:num.senses, Genre = 1:num.genres))
chi.prior.mean = 0 #this remains fixed throughout


#Make a list of all Snippet IDs in each time period and genre
SnippetIDs = setNames(replicate(num.periods, list( setNames(vector("list", length = num.genres), 
                                                            paste("Genre", 1:num.genres, sep = "")) )), 
                      paste("Time", 1:num.periods, sep = "") )
for (t in 1:num.periods) {
  for (g in 1:num.genres) {
    SnippetIDs[[t]][[g]] = snippets[with(snippets, Time == t & genre == g), "SnippetID"]
  }#for
}#for

SnippetIDs.by.time = setNames(vector("list", length = num.periods), paste("Time", 1:num.periods, sep = ""))
for (t in 1:num.periods) {
  SnippetIDs.by.time[[t]] = snippets[with(snippets, Time == t), "SnippetID"]
}#for


#Get the counts of words used in each snippet across all time periods and genres
snippets.expanded = gather(snippets, key = position, value = word, 1:snippet.length)[,-c(1,3,4)]
snippets.expanded = snippets.expanded[!is.na(snippets.expanded$word),]
word.counts.in.snippets = sparseMatrix(i = snippets.expanded$word, j = snippets.expanded$SnippetID, 
                                       x = 1, dims = c(num.words, num.snippets),
                                       dimnames = list(Word = 1:num.words, SnippetID = snippets$SnippetID), 
                                       use.last.ij = FALSE) 
rm(snippets.expanded); gc()


#Get initial probabilities for all snippets belonging to each sense
sense.probs = array(dim = c(num.senses, num.snippets),
                    dimnames = list(Sense = 1:num.senses, SnippetID = snippets$SnippetID))
for(t in 1:num.periods) {
  for(g in 1:num.genres) {
    sense.probs[,SnippetIDs[[t]][[g]]] = apply(snippets[SnippetIDs[[t]][[g]], 1:snippet.length], 1, get.sense.prob,
                                               phi.tilde[,g,t], psi.tilde[,,t])
  }#for
}#for


#Set scales for HMC proposals - these are the epsilon^2 parameters
phi.prop.scale = array(phi.prop.scale.init, dim = c(num.genres, num.periods), 
                       dimnames = list(Genre = 1:num.genres, Time = 1:num.periods) )
chi.prop.scale = chi.prop.scale.init
theta.prop.scale = array(theta.prop.scale.init, dim = c(num.periods), 
                         dimnames = list(Time = 1:num.periods))

#Set up arrays to hold acceptance counts - for HMC tuning
phi.jumps = array(0, dim = c(block.size, num.genres, num.periods), dimnames = 
                    list(Iteration = 1:block.size, Genre = 1:num.genres, Time = 1:num.periods))
chi.jumps = array(0, dim = c(block.size), dimnames = list(Iteration = 1:block.size))
theta.jumps = array(0, dim = c(block.size, num.periods), dimnames = 
                    list(Iteration = 1:block.size, Time = 1:num.periods))


#Set up arrays to hold acceptance counts - for acceptance rate calculation
accept.count.phi = array(0, dim = c(num.genres, num.periods), 
                         dimnames = list(Genre = 1:num.genres, Time = 1:num.periods)) #To calculate acceptance rate for phi for each time period and genre
accept.count.chi = 0 #To calculate acceptance rate for chi
accept.count.theta = array(0, dim = c(num.periods), 
                         dimnames = list(Sense = 1:num.periods)) #To calculate acceptance rate for theta for each time period


#We will alternately go forward and backward in time for each successive iteration
time.direction = 1:num.periods


#Set up arrays to hold simulated values
num.samples = num.iterations/N

phi.tilde.sim = array(dim = c(num.samples, num.senses, num.genres, num.periods),
                      dimnames = list(Sample = 1:num.samples, Sense = 1:num.senses, 
                                      Genre = 1:num.genres, Time = 1:num.periods))

psi.tilde.sim = array(dim = c(num.samples, num.words, num.senses, num.periods),
                      dimnames = list(Sample = 1:num.samples, Word = 1:num.words, 
                                      Sense = 1:num.senses, Time = 1:num.periods))

chi.sim = array(dim = c(num.samples, num.words, num.senses), 
                dimnames = list(Sample = 1:num.samples, Word = 1:num.words, Sense = 1:num.senses))

theta.sim = array(dim = c(num.samples, num.words, num.periods), 
                  dimnames = list(Sample = 1:num.samples, Word = 1:num.words, Time = 1:num.periods))

sense.probs.sim = array(dim = c(num.samples, num.senses, num.snippets),
                        dimnames = list(Sample = 1:num.samples, 
                                        Sense  = 1:num.senses, SnippetID = snippets$SnippetID))



#############################################################################
# Run MCMC ---------------------------------------------------------------- #
#############################################################################

progress.bar = txtProgressBar(min = 1 , max = num.iterations, style = 3) #Show progress
print(getwd())

save(snippets, snippets.info, words.used, num.words, num.snippets, snippet.length, 
     snippet.lengths, num.periods, num.genres, num.senses, file = "snippets.RData")
save(z, file = "initial.z.RData")
save(random.seed, num.runs, stop.tuning.after, num.iterations, N, num.samples,
     temper.phi, temper.chi, temper.theta, min.inv.temp, num.iter.temp,
     kappa.phi, kappa.chi, kappa.theta, alpha.phi, alpha.theta, 
     num.leapfrog.steps.phi, num.leapfrog.steps.chi, num.leapfrog.steps.theta, 
     target.accept.rate.phi, target.accept.rate.chi, target.accept.rate.theta,
     phi.prop.scale.init, chi.prop.scale.init, theta.prop.scale.init,
     block.size, c1.phi, c1.chi, c1.theta, #c0.phi, c0.chi, c0.theta,
     file = "other.parameters.RData")

for (b in 1:num.runs) {
  
  Start.Time = Sys.time() #Save start time
  for (i in 1:num.iterations) {
    j = (b-1)*num.iterations + i
    
    inv.temp = get.inv.temp(j, min.inv.temp, num.iter.temp)
    
    # --- chi updates --- #
    
    chi.HMC.proposal = propose.new.chi.HMC(num.leapfrog.steps.chi, sqrt(chi.prop.scale), chi, psi.tilde, chi.prior.mean, 
                                           kappa.chi, sense.probs, word.counts.in.snippets, snippet.lengths, temper.chi, inv.temp)
    
    #Accept or reject
    if(log(runif(1)) < chi.HMC.proposal$MH) {
      accept.count.chi = accept.count.chi + 1
      chi.jumps[i %% block.size + 1] = 1
      sense.probs = chi.HMC.proposal$sense.probs
      chi = chi.HMC.proposal$chi
      psi.tilde = chi.HMC.proposal$psi.tilde
    } else {
      chi.jumps[i %% block.size + 1] = 0
    }#else
    
    #Update HMC parameters
    if (j <= stop.tuning.after) {
      if (j >= block.size) {
        gamma1 = ((j+1)/block.size)^(-c1.chi)
        #gamma2 = c0.chi * gamma1
        accept.rate = sum(chi.jumps)/block.size
        chi.prop.scale = update.scale(chi.prop.scale, accept.rate, gamma1, target.accept.rate.chi)
      }#if
    }#if
    
    
    #Cycle through each time period
    for (t in time.direction) {
      
      # --- theta updates --- #
      
      if (t == 1) {
        theta.t.prior.mean = alpha.theta * theta[,t+1]
        theta.t.prior.var = kappa.theta
      } else {
        if (t == num.periods) {
          theta.t.prior.mean = alpha.theta * theta[,t-1]
          theta.t.prior.var = kappa.theta
        } else {
          theta.t.prior.mean = alpha.theta*(theta[,t-1]+theta[,t+1]) / (1+alpha.theta^2)
          theta.t.prior.var = kappa.theta / (1 + alpha.theta^2)
        }#else
      }#if
      
      theta.t.HMC.proposal = propose.new.theta.t.HMC(num.leapfrog.steps.theta, sqrt(theta.prop.scale[t]), theta[,t], 
                                                     psi.tilde[,,t], theta.t.prior.mean, theta.t.prior.var, 
                                                     sense.probs[,SnippetIDs.by.time[[t]], drop=FALSE], 
                                                     word.counts.in.snippets[,SnippetIDs.by.time[[t]]],
                                                     snippet.lengths[SnippetIDs.by.time[[t]]], temper.theta, inv.temp)
      
      #Accept or reject
      if(log(runif(1)) < theta.t.HMC.proposal$MH) {
        accept.count.theta[t] = accept.count.theta[t] + 1
        theta.jumps[i %% block.size + 1, t] = 1
        sense.probs[,SnippetIDs.by.time[[t]]] = theta.t.HMC.proposal$sense.probs.t
        theta[,t] = theta.t.HMC.proposal$theta.t
        psi.tilde[,,t] = theta.t.HMC.proposal$psi.t.tilde
      } else {
        theta.jumps[i %% block.size + 1, t] = 0
      }#else
      
      #Update HMC parameters
      if (j <= stop.tuning.after) {
        if (j >= block.size) {
          gamma1 = ((j+1)/block.size)^(-c1.theta)
          #gamma2 = c0.theta * gamma1
          accept.rate = sum(theta.jumps[,t])/block.size
          theta.prop.scale[t] = update.scale(theta.prop.scale[t], accept.rate, gamma1, target.accept.rate.theta)
        }#if
      }#if
      
      
      # --- phi updates --- #
      
      if (t == 1) {
        phi.t.prior.mean[] = alpha.phi * phi[,,t+1]
        phi.t.prior.var = kappa.phi
      } else {
        if (t == num.periods) {
          phi.t.prior.mean[] = alpha.phi * phi[,,t-1]
          phi.t.prior.var = kappa.phi
        } else {
          phi.t.prior.mean[] = alpha.phi*(phi[,,t-1]+phi[,,t+1]) / (1+alpha.phi^2)
          phi.t.prior.var = kappa.phi / (1 + alpha.phi^2)
        }#else
      }#if
      
      #Cycle through each genre
      for (g in 1:num.genres) {
        
        phi.g.t.HMC.proposal = propose.new.phi.g.t.HMC(num.leapfrog.steps.phi, sqrt(phi.prop.scale[g,t]), phi[,g,t], 
                                                       phi.tilde[,g,t], psi.tilde[,,t], phi.t.prior.mean[,g], phi.t.prior.var, 
                                                       sense.probs[,SnippetIDs[[t]][[g]], drop=FALSE], snippet.counts[g,t], 
                                                       temper.phi, inv.temp)
        
        #Accept or reject
        if(log(runif(1)) < phi.g.t.HMC.proposal$MH) {
          accept.count.phi[g,t] = accept.count.phi[g,t] + 1
          phi.jumps[i %% block.size + 1, g, t] = 1
          sense.probs[,SnippetIDs[[t]][[g]]] = phi.g.t.HMC.proposal$sense.probs.g.t
          phi[,g,t] = phi.g.t.HMC.proposal$phi.g.t
          phi.tilde[,g,t] = phi.g.t.HMC.proposal$phi.g.t.tilde
        } else {
          phi.jumps[i %% block.size + 1, g, t] = 0
        }#else
        
        #Update HMC parameters
        if (j <= stop.tuning.after) {
          if (j >= block.size) {
            gamma1 = ((j+1)/block.size)^(-c1.phi)
            #gamma2 = c0.phi * gamma1
            accept.rate = sum(phi.jumps[,g,t])/block.size
            phi.prop.scale[g,t] = update.scale(phi.prop.scale[g,t], accept.rate, gamma1, target.accept.rate.phi)
          }#if
        }#if
        
      }#for g
      
    }#for t
    
    
    #Store values
    if (i%%N == 0) {
      sense.probs.sim[i/N,,] = sense.probs
      phi.tilde.sim[i/N,,,] = phi.tilde
      psi.tilde.sim[i/N,,,] = psi.tilde
      theta.sim[i/N,,] = theta
      chi.sim[i/N,,] = chi
    }#if
    
    
    #Reverse time direction for next iteration
    time.direction = rev(time.direction)
    
    #Show progress
    setTxtProgressBar(progress.bar, i)
    
  }#for i
  close(progress.bar)
  End.Time = Sys.time() #Save finish time
  
  #Duration of algorithm
  Run.Time = End.Time - Start.Time
  print(Run.Time)
  
  #Write results to file
  filename.suffix = paste((b-1)*num.iterations+1, b*num.iterations, sep = "-")
  
  save(sense.probs.sim, file = paste("sense.probs.sim", filename.suffix, "RData", sep = "."))
  save(psi.tilde.sim, file = paste("psi.tilde.sim", filename.suffix, "RData", sep = "."))
  save(phi.tilde.sim, file = paste("phi.tilde.sim", filename.suffix, "RData", sep = "."))
  save(chi.sim, file = paste("chi.sim", filename.suffix, "RData", sep = "."))
  save(theta.sim, file = paste("theta.sim", filename.suffix, "RData", sep = "."))
  save(Start.Time, End.Time, Run.Time, file = paste("run.time", filename.suffix, "RData", sep = "."))
  save(phi.prop.scale, chi.prop.scale, theta.prop.scale, accept.count.phi,
       accept.count.chi, accept.count.theta, phi.jumps, chi.jumps, theta.jumps,
       file = paste("HMC.parameters.RData", filename.suffix, "RData", sep = "."))
  
  #Reset accept counts
  if (b != num.runs) {
    accept.count.phi[,] = 0
    accept.count.chi[] = 0
    accept.count.theta[] = 0
  }
  
  #Delete variables not required
  if (b == num.runs) {
    rm(b, i, j, t, g, gamma1, accept.rate, chi.HMC.proposal, theta.t.HMC.proposal, theta.t.prior.mean, 
       theta.t.prior.var, phi.g.t.HMC.proposal, phi.t.prior.var); gc()
  }
  
}#for b
