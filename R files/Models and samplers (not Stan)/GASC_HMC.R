#############################################################################
# Required libraries ------------------------------------------------------ #
#############################################################################

library(tidyr)
library(plyr)
library(MCMCpack)
library(invgamma)
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
num.senses = 2

#Set variance for psi
kappa.psi = 1/10

#Set initial variance for phi
kappa.phi = 1/4

#Prior hyperparameters for kappa.phi (which has inv Gamma distribution)
kappa.phi.shape.prior = 7
kappa.phi.rate.prior  = 3



#############################################################################
# Set MCMC parameters ----------------------------------------------------- #
#############################################################################

#Set seed for initialisation
random.seed = 500

#Specify number of MCMC iterations
num.iterations = 10000

#Sample every Nth iteration
N = 5

#Specify number of times to run the above MCMC iterations
num.runs = 1

#Stop tuning HMC scale parameters after this many MCMC iterations
stop.tuning.after = num.iterations * num.runs #num.iterations/2

#Set parameters for resampling kappa.phi
resample.kappa.phi.every = 50 #resample every this many iterations
resample.kappa.phi.after = 150 #start resampling from this iteration

#Set number of leapfrog steps for HMC proposals
num.leapfrog.steps.phi = 5
num.leapfrog.steps.psi = 10

#Set HMC tuning parameters - refer to Shaby and Wells (2010)
target.accept.rate.phi = 0.651
target.accept.rate.psi = 0.651
block.size = 10 #k in the paper

#c0.phi = 1
c1.phi = 0.8
#c0.psi = 1
c1.psi = 0.8

#Set initial scales for HMC proposals - these are the epsilon^2 parameters
phi.prop.scale.init = 2.4^2 / num.senses / num.leapfrog.steps.phi
psi.prop.scale.init = 2.4^2 / num.words / num.leapfrog.steps.psi



#############################################################################
# Functions for parameter updates ----------------------------------------- #
#############################################################################

#Write a function to update the rate hyperparameter for kappa.phi using current phi
update.kappa.phi.rate.single.genre = function(phi.g) {
  
  differences = t(diff(t(phi.g)))
  
  sum(differences^2)
  
}#update.kappa.phi.rate.single.genre

update.kappa.phi.rate = function(phi) {
  
  kappa.phi.rate.prior + 0.25 * sum(apply(phi, 2, update.kappa.phi.rate.single.genre))
  
}#update.kappa.phi.rate


#Write a function to update kappa.phi given hyperparameters
update.kappa.phi = function(kappa.phi.shape, kappa.phi.rate) {
  
  invgamma::rinvgamma(1, shape = kappa.phi.shape, rate = kappa.phi.rate)
  
}#update.kappa.phi


#Write a function that takes as input a single snippet, phi and psi, 
#and returns a vector of (unnormalised) sense probabilities for that snippet
get.sense.prob = function(snippet, phi.tilde, psi.tilde) {
  apply(psi.tilde[snippet,], 2, prod, na.rm = TRUE) * phi.tilde
}#get.sense.prob

#Use the following function instead when only the probability of belonging to a single sense is required
get.single.sense.prob = function(snippet, phi.k.tilde, psi.k.tilde) {
  prod(psi.k.tilde[snippet], na.rm = TRUE) * phi.k.tilde
}#get.single.sense.prob


#Write a function to calculate the gradient of the log posterior for phi genre g time t
get.phi.gradient = function(phi, phi.tilde, phi.prior.mean, phi.prior.var, num.snippets, sense.probs) {
  
  - (phi - phi.prior.mean) / phi.prior.var - num.snippets * phi.tilde +
    colSums(t(sense.probs) / colSums(sense.probs))
  
}#get.phi.gradient


#Write a function to calculate the gradient  for psi sense k time t
get.psi.k.gradient = function(k, psi.k, psi.k.tilde, psi.prior.mean, psi.prior.var, sense.probs, word.counts.in.snippets, snippet.lengths) {
  
  - (psi.k - psi.prior.mean) / psi.prior.var +
    as.numeric( (word.counts.in.snippets - outer(psi.k.tilde, snippet.lengths)) %*% (sense.probs[k,] / colSums(sense.probs)) )
  
}#get.psi.k.gradient


#Write a function to propose phi for genre g time t using HMC
#Function adapted from Neal, R.M., 2011. MCMC using Hamiltonian dynamics
propose.new.phi.g.t.HMC = function(num.steps, epsilon, phi.g.t, phi.g.t.tilde, phi.g.t.prior.mean, phi.t.prior.var,
                                   sense.probs.g.t, snippet.counts.g.t){
  
  #Generate initial momentum vector
  p = rnorm(num.senses)
  
  #Evaluate potential and kinetic energies at start of trajectory
  current_U = -(sum(dnorm(phi.g.t, phi.g.t.prior.mean, sqrt(phi.t.prior.var), log = TRUE)) + 
                  sum(log(colSums(sense.probs.g.t))))
  current_K = sum(p^2) / 2
  
  #Make a half step for momentum at the beginning
  grad_U = -get.phi.gradient(phi.g.t, phi.g.t.tilde, phi.g.t.prior.mean, phi.t.prior.var, snippet.counts.g.t, 
                             sense.probs.g.t)
  p = p - epsilon * grad_U / 2
  
  #Alternate full steps for position and momentum
  for (s in 1:num.steps) {
    
    #Make a full step for the position
    phi.g.t = phi.g.t + epsilon * p
    exp.phi.g.t = exp(phi.g.t)
    phi.g.t.tilde = exp.phi.g.t / sum(exp.phi.g.t)
    
    sense.probs.g.t[] = apply(snippets[SnippetIDs[[t]][[g]], 1:snippet.length], 1, get.sense.prob, 
                              phi.g.t.tilde, psi.tilde[,,t])
    
    grad_U = -get.phi.gradient(phi.g.t, phi.g.t.tilde, phi.g.t.prior.mean, phi.t.prior.var, snippet.counts.g.t, 
                               sense.probs.g.t)
    
    #Make a full step for the momentum, except at end of trajectory
    if (s != num.steps) p = p - epsilon * grad_U
    
  }#for s
  
  #Make a half step for momentum at the end
  p = p - epsilon * grad_U / 2
  
  #Evaluate potential and kinetic energies at end of trajectory
  proposed_U = -(sum(dnorm(phi.g.t, phi.g.t.prior.mean, sqrt(phi.t.prior.var), log = TRUE)) + 
                   sum(log(colSums(sense.probs.g.t))))
  proposed_K = sum(p^2) / 2
  
  #Compute log MH ratio given current and proposed energies
  MH = current_U - proposed_U + current_K - proposed_K
  
  
  out = list(MH = MH, phi.g.t = phi.g.t, phi.g.t.tilde = phi.g.t.tilde, sense.probs.g.t = sense.probs.g.t)
  
  return(out)
  
}#propose.new.phi.g.t.HMC


#Write a function to propose psi for sense k time t using HMC
#Function adapted from Neal, R.M., 2011. MCMC using Hamiltonian dynamics
propose.new.psi.k.t.HMC = function(num.steps, epsilon, k, psi.k.t, psi.k.t.tilde, psi.k.t.prior.mean, psi.t.prior.var,
                                   sense.probs.t, word.counts.in.snippets.t, snippet.lengths.t){
  
  #Generate initial momentum vector
  p = rnorm(num.words)
  
  #Evaluate potential and kinetic energies at start of trajectory
  current_U = -(sum(dnorm(psi.k.t, psi.k.t.prior.mean, sqrt(psi.t.prior.var), log = TRUE)) + 
                  sum(log(colSums(sense.probs.t))))
  current_K = sum(p^2) / 2
  
  #Make a half step for momentum at the beginning
  grad_U = -get.psi.k.gradient(k, psi.k.t, psi.k.t.tilde, psi.k.t.prior.mean, psi.t.prior.var,
                               sense.probs.t, word.counts.in.snippets.t, snippet.lengths.t)
  p = p - epsilon * grad_U / 2
  
  #Alternate full steps for position and momentum
  for (s in 1:num.steps) {
    
    #Make a full step for the position
    psi.k.t = psi.k.t + epsilon * p
    exp.psi.k.t = exp(psi.k.t)
    psi.k.t.tilde = exp.psi.k.t / sum(exp.psi.k.t)
    
    for(g in 1:num.genres) {
      sense.probs.t[k,paste(SnippetIDs[[t]][[g]])] = apply(snippets[SnippetIDs[[t]][[g]], 1:snippet.length], 1,
                                                           get.single.sense.prob, phi.tilde[k,g,t], psi.k.t.tilde)
    }#for
    
    grad_U = -get.psi.k.gradient(k, psi.k.t, psi.k.t.tilde, psi.k.t.prior.mean, psi.t.prior.var,
                                 sense.probs.t, word.counts.in.snippets.t, snippet.lengths.t)
    
    #Make a full step for the momentum, except at end of trajectory
    if (s != num.steps) p = p - epsilon * grad_U
    
  }#for s
  
  #Make a half step for momentum at the end
  p = p - epsilon * grad_U / 2
  
  #Evaluate potential and kinetic energies at end of trajectory
  proposed_U = -(sum(dnorm(psi.k.t, psi.k.t.prior.mean, sqrt(psi.t.prior.var), log = TRUE)) + 
                   sum(log(colSums(sense.probs.t))))
  proposed_K = sum(p^2) / 2
  
  #Compute log MH ratio given current and proposed energies
  MH = current_U - proposed_U + current_K - proposed_K
  
  
  out = list(MH = MH, psi.k.t = psi.k.t, psi.k.t.tilde = psi.k.t.tilde, sense.probs.t = sense.probs.t)
  
  return(out)
  
}#propose.new.psi.k.t.HMC



#############################################################################
# Function for tuning HMC parameters -------------------------------------- #
#############################################################################

# Refer to Shaby and Wells (2010) "Exploring an Adaptive Metropolis Algorithm"

update.scale = function(scale, accept.rate, gamma2, target.accept.rate) {
  
  scale * exp(gamma2 * (accept.rate - target.accept.rate))
  
}#update.scale



#############################################################################
# Initialisation ---------------------------------------------------------- #
#############################################################################

#Get the counts of snippets in each time period and genre
snippet.counts = table(factor(snippets$genre, levels = 1:num.genres),
                       factor(snippets$Time, levels = 1:num.periods),
                       dnn = c("Genre", "Time"))

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

#Adjust psi for identifiability
psi = sapply(1:num.periods, function(t) {t(t(psi[,,t]) - colMeans(psi[,,t]))}, simplify = "array")
dimnames(psi) = dimnames(psi.tilde)


#Posterior shape hyperparameter for kappa.phi - this remains fixed throughout
kappa.phi.shape = kappa.phi.shape.prior + 0.5 * num.senses * num.genres * (num.periods-1)


#Set up structures for phi and psi prior means (to avoid errors later when num.genres = 1)
phi.prior.mean = array(data = 0, dim = c(num.senses, num.genres), dimnames = list(Sense = 1:num.senses, Genre = 1:num.genres))
psi.prior.mean = array(data = 0, dim = c(num.words, num.senses), dimnames = list(Word = 1:num.words, Sense = 1:num.senses))


#Make a list of all Snippet IDs in each time period and genre
SnippetIDs = setNames(replicate(num.periods, list( setNames(vector("list", length = num.genres), 
                                                            paste("Genre", 1:num.genres, sep = "")) )), 
                      paste("Time", 1:num.periods, sep = "") )
for (t in 1:num.periods) {
  for (g in 1:num.genres) {
    SnippetIDs[[t]][[g]] = snippets[with(snippets, Time == t & genre == g), "SnippetID"]
  }#for
}#for

SnippetIDs.by.time = lapply(SnippetIDs, unlist)


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
psi.prop.scale = array(psi.prop.scale.init, dim = c(num.senses, num.periods), 
                       dimnames = list(Sense = 1:num.senses, Time = 1:num.periods))

#Set up arrays to hold acceptance counts - for HMC tuning
phi.jumps = array(0, dim = c(block.size, num.genres, num.periods), 
                  dimnames = list(Iteration = 1:block.size, Genre = 1:num.genres, Time = 1:num.periods))

psi.jumps = array(0, dim = c(block.size, num.senses, num.periods), dimnames = 
                    list(Iteration = 1:block.size, Sense = 1:num.senses, Time = 1:num.periods))


#Set up arrays to hold acceptance counts - for acceptance rate calculation
accept.count.phi = array(0, dim = c(num.genres, num.periods), 
                         dimnames = list(Genre = 1:num.genres, Time = 1:num.periods)) #To calculate acceptance rate for phi for each time period and genre

accept.count.psi = array(0, dim = c(num.senses, num.periods), 
                         dimnames = list(Sense = 1:num.senses, Time = 1:num.periods)) #To calculate acceptance rate for psi for each time period and sense


#We will alternately go forward and backward in time for each successive iteration
time.direction = 1:num.periods
#And also alternately go forward and backward over senses for each successive iteration for psi
sense.direction = 1:num.senses


#Set up arrays to hold simulated values
num.samples = num.iterations/N

phi.tilde.sim = array(dim = c(num.samples, num.senses, num.genres, num.periods),
                      dimnames = list(Sample = 1:num.samples, Sense = 1:num.senses, 
                                      Genre = 1:num.genres, Time = 1:num.periods))

psi.tilde.sim = array(dim = c(num.samples, num.words, num.senses, num.periods),
                      dimnames = list(Sample = 1:num.samples, 
                                      Word = 1:num.words, Sense = 1:num.senses, Time = 1:num.periods))

kappa.phi.sim = numeric(num.samples)

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
     kappa.phi, kappa.psi, num.leapfrog.steps.phi, num.leapfrog.steps.psi,
     target.accept.rate.phi, target.accept.rate.psi, phi.prop.scale.init, psi.prop.scale.init,
     block.size, c1.phi, c1.psi, #c0.phi, c0.psi,
     file = "other.parameters.RData")

for (b in 1:num.runs) {
  
  Start.Time = Sys.time() #Save start time
  for (i in 1:num.iterations) {
    j = (b-1)*num.iterations + i
    
    #Cycle through each time period
    for (t in time.direction) {
      
      # --- psi updates --- #
      
      if (t == 1) {
        psi.prior.mean[] = psi[,,t+1]
        psi.prior.var = 2*kappa.psi
      } else {
        if (t == num.periods) {
          psi.prior.mean[] = psi[,,t-1]
          psi.prior.var = 2*kappa.psi
        } else {
          psi.prior.mean[] = (psi[,,t-1] + psi[,,t+1]) / 2
          psi.prior.var = kappa.psi
        }#else
      }#if
      
      #Cycle through each sense in turn
      for (k in sense.direction) {
        
        psi.k.t.HMC.proposal = propose.new.psi.k.t.HMC(num.leapfrog.steps.psi, sqrt(psi.prop.scale[k,t]), k, 
                                                       psi[,k,t], psi.tilde[,k,t], psi.prior.mean[,k], psi.prior.var, 
                                                       sense.probs[,SnippetIDs.by.time[[t]], drop=FALSE], 
                                                       word.counts.in.snippets[,SnippetIDs.by.time[[t]]], 
                                                       snippet.lengths[SnippetIDs.by.time[[t]]])
        
        #Accept or reject
        if(log(runif(1)) < psi.k.t.HMC.proposal$MH) {
          accept.count.psi[k,t] = accept.count.psi[k,t] + 1
          psi.jumps[i %% block.size + 1, k, t] = 1
          sense.probs[,SnippetIDs.by.time[[t]]] = psi.k.t.HMC.proposal$sense.probs.t
          psi[,k,t] = psi.k.t.HMC.proposal$psi.k.t 
          psi.tilde[,k,t] = psi.k.t.HMC.proposal$psi.k.t.tilde
        } else {
          psi.jumps[i %% block.size + 1, k, t] = 0
        }#else
        
        #Update HMC parameters
        if (j <= stop.tuning.after) {
          if (j >= block.size) {
            gamma1 = ((j+1)/block.size)^(-c1.psi)
            #gamma2 = c0.psi * gamma1
            accept.rate = sum(psi.jumps[,k,t])/block.size
            psi.prop.scale[k,t] = update.scale(psi.prop.scale[k,t], accept.rate, gamma1, target.accept.rate.psi)
          }#if
        }#if
        
      }#for k
      
      
      # --- phi updates --- #
      
      if (t == 1) {
        phi.prior.mean[] = phi[,,t+1]
        phi.prior.var = 2*kappa.phi
      } else {
        if (t == num.periods) {
          phi.prior.mean[] = phi[,,t-1]
          phi.prior.var = 2*kappa.phi
        } else {
          phi.prior.mean[] = (phi[,,t-1] + phi[,,t+1]) / 2
          phi.prior.var = kappa.phi
        }#else
      }#if
      
      #Cycle through each genre
      for (g in 1:num.genres) {
        
        phi.g.t.HMC.proposal = propose.new.phi.g.t.HMC(num.leapfrog.steps.phi, sqrt(phi.prop.scale[g,t]), phi[,g,t], 
                                                       phi.tilde[,g,t], phi.prior.mean[,g], phi.prior.var, 
                                                       sense.probs[,SnippetIDs[[t]][[g]], drop=FALSE], snippet.counts[g,t])
        
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
    
    
    # --- kappa.phi update --- #
    
    if( (j >= resample.kappa.phi.after) & ((j %% resample.kappa.phi.every) == 0) ) { #Check whether we need to resample
      kappa.phi.rate = update.kappa.phi.rate(phi)
      kappa.phi = update.kappa.phi(kappa.phi.shape, kappa.phi.rate)
    }#if
    
    
    #Store values
    if (i%%N == 0) {
      psi.tilde.sim[i/N,,,] = psi.tilde
      phi.tilde.sim[i/N,,,] = phi.tilde
      sense.probs.sim[i/N,,] = sense.probs
      kappa.phi.sim[i/N] = kappa.phi
    }#if
    
    
    #Reverse sense direction for next iteration
    sense.direction = rev(sense.direction)
    
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

  save(psi.tilde.sim, file = paste("psi.tilde.sim", filename.suffix, "RData", sep = "."))
  save(phi.tilde.sim, file = paste("phi.tilde.sim", filename.suffix, "RData", sep = "."))
  save(kappa.phi.sim, file = paste("kappa.phi.sim", filename.suffix, "RData", sep = "."))
  save(sense.probs.sim, file = paste("sense.probs.sim", filename.suffix, "RData", sep = "."))
  save(Start.Time, End.Time, Run.Time, file = paste("run.time", filename.suffix, "RData", sep = "."))
  save(phi.prop.scale, psi.prop.scale, accept.count.phi, accept.count.psi, phi.jumps, psi.jumps, 
       file = paste("HMC.parameters.RData", filename.suffix, "RData", sep = "."))
  
  #Reset accpet counts
  if (b != num.runs) {
    accept.count.phi[,] = 0
    accept.count.psi[,] = 0
  }
  
  #Delete variables not required
  if (b == num.runs) {
    rm(b, i, j, t, g, gamma1, accept.rate, psi.k.t.HMC.proposal, psi.prior.var,
       phi.g.t.HMC.proposal, phi.t.prior.var); gc()
  }
  
}#for b
