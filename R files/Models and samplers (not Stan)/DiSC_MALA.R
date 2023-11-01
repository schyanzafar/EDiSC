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
# Parameter updates ------------------------------------------------------- #
#############################################################################

#Write a function that takes as input a single snippet, phi and psi, 
#and returns a vector of (unnormalised) sense probabilities for that snippet
get.sense.prob = function(snippet, phi.tilde, psi.tilde) {
  apply(psi.tilde[snippet,], 2, prod, na.rm = TRUE) * phi.tilde
}#get.sense.prob


#Write a function to calculate the MALA gradient for phi updates for a single time slice
get.phi.MALA.gradient = function(phi, phi.tilde, phi.prior.mean, phi.prior.var, num.snippets, sense.probs) {
  
  grad = - (phi - phi.prior.mean) / phi.prior.var - num.snippets * phi.tilde +
    colSums(t(sense.probs) / colSums(sense.probs))
  
  grad
  
}#get.phi.MALA.gradient


#Write a function to calculate the MALA gradient for chi updates for all senses
get.chi.MALA.gradient = function(chi, psi.tilde, chi.prior.mean, chi.prior.var, sense.probs, word.counts.in.snippets, snippet.lengths) {

  P = t(sense.probs)/colSums(sense.probs)

  I.P.psi.tilde = sapply(1:num.senses, function(l){psi.tilde[,l,snippets$Time] %*% (snippet.lengths*P[,l])})

  grad = - (chi - chi.prior.mean) / chi.prior.var +
    as.matrix(word.counts.in.snippets %*% P) - I.P.psi.tilde

  grad

}#get.chi.MALA.gradient


#Write a function to calculate the MALA gradient for theta updates for a single time slice
get.theta.MALA.gradient = function(theta, psi.tilde, theta.prior.mean, theta.prior.var, sense.probs, word.counts.in.snippets, snippet.lengths) {
  
  grad = as.numeric(
    - (theta - theta.prior.mean) / theta.prior.var +
      rowSums(word.counts.in.snippets - t(snippet.lengths * ((t(sense.probs) / colSums(sense.probs)) %*% t(psi.tilde))))
  )
  
  grad
  
}#get.theta.MALA.gradient


#Write a function to generate MALA proposals for phi, chi and theta given given other variables
propose.new.MALA = function(X, size, gradient, prop.scale) {
  
  X + rnorm(size, prop.scale/2 * gradient, (prop.scale)^0.5)
  
}#propose.new.MALA


#Write a function to calcualte the log MH ratio given current and proposed phi and theta
get.mh.ratio = function(X, X.new, X.prior.mean, X.prior.var, sense.probs, sense.probs.new, gradient, gradient.new, prop.scale) {
  
  sum(dnorm(X.new, X.prior.mean, sqrt(X.prior.var), log = TRUE) - 
        dnorm(X, X.prior.mean, sqrt(X.prior.var), log = TRUE)) +
    
    sum(dnorm(X, X.new + prop.scale/2 * gradient.new, sqrt(prop.scale), log = TRUE) -
          dnorm(X.new, X + prop.scale/2 * gradient, sqrt(prop.scale), log = TRUE)) +
    
    sum(log(colSums(sense.probs.new)) - 
          log(colSums(sense.probs)))
  
}#get.mh.ratio



#############################################################################
# Tuning MALA parameters -------------------------------------------------- #
#############################################################################

# Refer to Shaby and Wells (2010) "Exploring an Adaptive Metropolis Algorithm"

update.scale = function(scale, accept.rate, gamma2, target.accept.rate) {
  
  scale * exp(gamma2 * (accept.rate - target.accept.rate))
  
}#update.scale



#############################################################################
# MCMC setup -------------------------------------------------------------- #
#############################################################################

# num.senses = 3
# num.genres = 1

# #Initialise parameters AT ZERO
# phi = array(data = 0, dim = c(num.senses, num.genres, num.periods),
#             dimnames = list(Sense = 1:num.senses, Genre = 1:num.genres, Time = 1:num.periods))
# phi.tilde = aperm(aaply(exp(phi), 2, function(x) t(t(x)/colSums(x)), .drop = FALSE), c(2,1,3))
# 
# chi = array(data = 0, dim = c(num.words, num.senses),
#             dimnames = list(Word = 1:num.words, Sense = 1:num.senses))
# 
# theta = array(data = 0, dim = c(num.words, num.periods),
#             dimnames = list(Word = 1:num.words, Time = 1:num.periods))

#Get the counts of snippets in each time period and genre
snippet.counts = table(factor(snippets$genre, levels = 1:num.genres),
                       factor(snippets$Time, levels = 1:num.periods),
                       dnn = c("Genre", "Time"))


#Initialise parameters AT RANDOM

#Initialise z randomly
random.seed = 500
set.seed(random.seed)
z = sample(1:num.senses, num.snippets, TRUE)
# z = as.numeric(factor(snippets.info$sense.id))

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

rm(vkt.combos, design.matrix.idx, design.matrix, psi.fit, chi.mean); gc()


#Set coefficient alpha in the model X(t) = alpha * X(t-1) + noise
alpha.phi = 0.9
alpha.theta = 0.9


#Set variances for phi, chi and theta
kappa.phi = 0.25
kappa.theta = 0.25
kappa.chi = 1.25


#Set up structures for prior means (to avoid errors later for phi when num.genres = 1)
phi.prior.mean = array(data = 0, dim = c(num.senses, num.genres), dimnames = list(Sense = 1:num.senses, Genre = 1:num.genres))
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


#Set scales for MALA proposals - these are the sigma^2 parameters
phi.MALA.prop.scale = array(2.4^2 / num.senses, dim = c(num.genres, num.periods), 
                            dimnames = list(Genre = 1:num.genres, Time = 1:num.periods) )

chi.MALA.prop.scale = 2.4^2 / (num.words * num.senses)

theta.MALA.prop.scale = array(2.4^2 / num.words, dim = c(num.periods), dimnames = list(Time = 1:num.periods))


#Set MALA tuning parameters - refer to Shaby and Wells (2010)
target.accept.rate.phi = 0.574
target.accept.rate.chi = 0.574
target.accept.rate.theta = 0.574
block.size = 10 #k in the paper

#c0.phi = 1
c1.phi = 0.8

#c0.chi = 1
c1.chi = 0.8

#c0.theta = 1
c1.theta = 0.8


#Set up arrays to hold acceptance counts - for MALA tuning
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
                         dimnames = list(Time = 1:num.periods)) #To calculate acceptance rate for theta for each time period


#We will alternately go forward and backward in time for each successive iteration
time.direction = 1:num.periods


#Specify number of MCMC iterations
num.iterations = 5000

#Sample every Nth iteration
N = 1
num.samples = num.iterations/N

#Specify number of times to run the above MCMC iterations
num.runs = 1

#Stop tuning MALA scale parameters after this many runs
stop.MALA.tuning.after = 99


#Set up arrays to hold simulated values
phi.tilde.sim = array(dim = c(num.samples, num.senses, num.genres, num.periods),
                      dimnames = list(Sample = 1:num.samples, Sense = 1:num.senses, 
                                      Genre = 1:num.genres, Time = 1:num.periods))

psi.tilde.sim = array(dim = c(num.samples, num.words, num.senses, num.periods),
                      dimnames = list(Sample = 1:num.samples, 
                                      Word = 1:num.words, Sense = 1:num.senses, Time = 1:num.periods))

chi.sim = array(dim = c(num.samples, num.words, num.senses), 
                dimnames = list(Sample = 1:num.samples, Word = 1:num.words, Sense = 1:num.senses))

theta.sim = array(dim = c(num.samples, num.words, num.periods), 
                  dimnames = list(Sample = 1:num.samples, Word = 1:num.words, Time = 1:num.periods))

sense.probs.sim = array(dim = c(num.samples, num.senses, num.snippets),
                        dimnames = list(Sample = 1:num.samples, 
                                        Sense  = 1:num.senses, SnippetID = snippets$SnippetID))



### ------------------------------ Run MCMC ------------------------------ ###

progress.bar = txtProgressBar(min = 1 , max = num.iterations, style = 3) #Show progress

save(snippets, snippets.info, words.used, num.words, num.snippets, snippet.length, snippet.lengths, num.periods, num.genres, num.senses, file = "snippets.RData")
save(z, file = "initial.z.RData")
save(num.iterations, N, kappa.phi, kappa.chi, kappa.theta, alpha.phi, alpha.theta, random.seed, file = "other.parameters.RData")

for (b in 1:num.runs) {
  
  Start.Time = Sys.time() #Save start time
  for (i in 1:num.iterations) {
    
    # --- chi updates --- #

    #Get chi gradient
    chi.gradient = get.chi.MALA.gradient(chi, psi.tilde, chi.prior.mean, kappa.chi, sense.probs,
                                             word.counts.in.snippets, snippet.lengths)

    #Propose new chi
    chi.new = propose.new.MALA(chi, num.words*num.senses, chi.gradient, chi.MALA.prop.scale)
    
    #Calculate new psi
    psi.new = sapply(1:num.periods, function(t) {sapply(1:num.senses, function(k) {
      chi.new[,k]+theta[,t]} )}, simplify = "array")
    psi.tilde.new = aperm(aaply(exp(psi.new), 3, function(x) t(t(x)/colSums(x))), c(2,3,1))

    #Get probabilities for all snippets belonging to each sense
    sense.probs.new = sense.probs
    for(t in 1:num.periods) {
      for(g in 1:num.genres) {
        sense.probs.new[,SnippetIDs[[t]][[g]]] = apply(snippets[SnippetIDs[[t]][[g]], 1:snippet.length], 1, get.sense.prob,
                                                   phi.tilde[,g,t], psi.tilde.new[,,t])
      }#for
    }#for

    #Get new chi gradient
    chi.gradient.new = get.chi.MALA.gradient(chi.new, psi.tilde.new, chi.prior.mean, kappa.chi, sense.probs.new,
                                         word.counts.in.snippets, snippet.lengths)

    #Calculate MH ratio
    MH = get.mh.ratio(chi, chi.new, chi.prior.mean, kappa.chi, sense.probs, sense.probs.new, chi.gradient, chi.gradient.new,
                      chi.MALA.prop.scale)

    #Accept or reject
    if(log(runif(1)) < MH) {
      accept.count.chi = accept.count.chi + 1
      chi.jumps[i %% block.size + 1] = 1
      sense.probs = sense.probs.new
      chi = chi.new
      psi.tilde = psi.tilde.new
    } else {
      chi.jumps[i %% block.size + 1] = 0
    }#else

    #Update MALA parameters
    if (b <= stop.MALA.tuning.after) {
      j = (b-1)*num.iterations + i
      if (j >= block.size) {
        gamma1 = ((j+1)/block.size)^(-c1.chi)
        #gamma2 = c0.chi * gamma1

        accept.rate = sum(chi.jumps)/block.size
        chi.MALA.prop.scale = update.scale(chi.MALA.prop.scale, accept.rate, gamma1, target.accept.rate.chi)
      }#if
    }#if

    
    #Cycle through each time period
    for (t in time.direction) {

      # --- theta updates --- #

      sense.probs.t = sense.probs[,SnippetIDs.by.time[[t]], drop=FALSE]

      if (t == 1) {
        theta.prior.mean = alpha.theta * theta[,t+1]
        theta.prior.var = kappa.theta
      } else {
        if (t == num.periods) {
          theta.prior.mean = alpha.theta * theta[,t-1]
          theta.prior.var = kappa.theta
        } else {
          theta.prior.mean = alpha.theta*(theta[,t-1]+theta[,t+1]) / (1+alpha.theta^2)
          theta.prior.var = kappa.theta / (1 + alpha.theta^2)
        }#else
      }#if

      #Get theta gradient
      theta.t.gradient = get.theta.MALA.gradient(theta[,t], psi.tilde[,,t], theta.prior.mean, theta.prior.var, 
                                                 sense.probs.t, word.counts.in.snippets[,SnippetIDs.by.time[[t]]],
                                                 snippet.lengths[SnippetIDs.by.time[[t]]])
      
      #Propose new theta
      theta.t.new = propose.new.MALA(theta[,t], num.words, theta.t.gradient, theta.MALA.prop.scale[t])
      
      #Calculate new psi
      psi.t.new = theta.t.new + chi
      exp.psi.t.new = exp(psi.t.new)
      psi.t.tilde.new = t(t(exp.psi.t.new) / colSums(exp.psi.t.new))

      #Get probabilities for all snippets belonging to each sense
      sense.probs.t.new = sense.probs.t
      for(g in 1:num.genres) {
        sense.probs.t.new[,paste(SnippetIDs[[t]][[g]])] = apply(snippets[SnippetIDs[[t]][[g]], 1:snippet.length], 1,
                                                                get.sense.prob, phi.tilde[,g,t], psi.t.tilde.new)
      }#for

      #Get new theta gradient
      theta.t.gradient.new = get.theta.MALA.gradient(theta.t.new, psi.t.tilde.new, theta.prior.mean, theta.prior.var,
                                                     sense.probs.t.new, word.counts.in.snippets[,SnippetIDs.by.time[[t]]],
                                                     snippet.lengths[SnippetIDs.by.time[[t]]])

      #Calculate MH ratio
      MH = get.mh.ratio(theta[,t], theta.t.new, theta.prior.mean, theta.prior.var, 
                        sense.probs.t, sense.probs.t.new, theta.t.gradient, theta.t.gradient.new, theta.MALA.prop.scale[t])

      #Accept or reject
      if(log(runif(1)) < MH) {
        accept.count.theta[t] = accept.count.theta[t] + 1
        theta.jumps[i %% block.size + 1, t] = 1
        sense.probs[,SnippetIDs.by.time[[t]]] = sense.probs.t.new
        theta[,t] = theta.t.new
        psi.tilde[,,t] = psi.t.tilde.new
      } else {
        theta.jumps[i %% block.size + 1, t] = 0
      }#else

      #Update MALA parameters
      if (b <= stop.MALA.tuning.after) {
        j = (b-1)*num.iterations + i
        if (j >= block.size) {
          gamma1 = ((j+1)/block.size)^(-c1.theta)
          #gamma2 = c0.theta * gamma1

          accept.rate = sum(theta.jumps[,t])/block.size
          theta.MALA.prop.scale[t] = update.scale(theta.MALA.prop.scale[t], accept.rate, gamma1, target.accept.rate.theta)
        }#if
      }#if


      # --- phi updates --- #

      if (t == 1) {
        phi.prior.mean[] = alpha.phi * phi[,,t+1]
        phi.prior.var = kappa.phi
      } else {
        if (t == num.periods) {
          phi.prior.mean[] = alpha.phi * phi[,,t-1]
          phi.prior.var = kappa.phi
        } else {
          phi.prior.mean[] = alpha.phi*(phi[,,t-1]+phi[,,t+1]) / (1+alpha.phi^2)
          phi.prior.var = kappa.phi / (1 + alpha.phi^2)
        }#else
      }#if

      #Cycle through each genre
      for (g in 1:num.genres) {

        sense.probs.g.t = sense.probs[,SnippetIDs[[t]][[g]], drop=FALSE]

        #Get phi gradient
        phi.g.t.gradient = get.phi.MALA.gradient(phi[,g,t], phi.tilde[,g,t], phi.prior.mean[,g], phi.prior.var, 
                                                 snippet.counts[g,t], sense.probs.g.t)

        #Propose new phi
        phi.g.t.new = propose.new.MALA(phi[,g,t], num.senses, phi.g.t.gradient, phi.MALA.prop.scale[g,t])
        exp.phi.g.t.new = exp(phi.g.t.new)
        phi.g.t.tilde.new = exp.phi.g.t.new / sum(exp.phi.g.t.new)

        #Get probabilities for all snippets belonging to each sense
        sense.probs.g.t.new = sense.probs.g.t
        sense.probs.g.t.new[] = apply(snippets[SnippetIDs[[t]][[g]], 1:snippet.length], 1, get.sense.prob, 
                                      phi.g.t.tilde.new, psi.tilde[,,t])

        #Get new phi gradient
        phi.g.t.gradient.new = get.phi.MALA.gradient(phi.g.t.new, phi.g.t.tilde.new, phi.prior.mean[,g], 
                                                     phi.prior.var, snippet.counts[g,t], sense.probs.g.t.new)

        #Calculate MH ratio
        MH = get.mh.ratio(phi[,g,t], phi.g.t.new, phi.prior.mean[,g], phi.prior.var, 
                          sense.probs.g.t, sense.probs.g.t.new, phi.g.t.gradient, phi.g.t.gradient.new, 
                          phi.MALA.prop.scale[g,t])

        #Accept or reject
        if(log(runif(1)) < MH) {
          accept.count.phi[g,t] = accept.count.phi[g,t] + 1
          phi.jumps[i %% block.size + 1, g, t] = 1
          sense.probs[,SnippetIDs[[t]][[g]]] = sense.probs.g.t.new
          phi[,g,t] = phi.g.t.new
          phi.tilde[,g,t] = phi.g.t.tilde.new
        } else {
          phi.jumps[i %% block.size + 1, g, t] = 0
        }#else

        #Update MALA parameters
        if (b <= stop.MALA.tuning.after) {
          j = (b-1)*num.iterations + i
          if (j >= block.size) {
            gamma1 = ((j+1)/block.size)^(-c1.phi)
            #gamma2 = c0.phi * gamma1

            accept.rate = sum(phi.jumps[,g,t])/block.size
            phi.MALA.prop.scale[g,t] = update.scale(phi.MALA.prop.scale[g,t], accept.rate, gamma1, target.accept.rate.phi)
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
    
    #Show pro,gress
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
  save(chi.sim, file = paste("chi.sim", filename.suffix, "RData", sep = "."))
  save(theta.sim, file = paste("theta.sim", filename.suffix, "RData", sep = "."))
  save(sense.probs.sim, file = paste("sense.probs.sim", filename.suffix, "RData", sep = "."))
  save(Start.Time, End.Time, Run.Time, file = paste("run.time", filename.suffix, "RData", sep = "."))
  save(phi.MALA.prop.scale, chi.MALA.prop.scale, theta.MALA.prop.scale, accept.count.phi, accept.count.chi,
       accept.count.theta, phi.jumps, chi.jumps, theta.jumps,
       file = paste("MALA.parameters.RData", filename.suffix, "RData", sep = "."))
  
  #Reset accpet counts
  if (b != num.runs) {
    accept.count.phi[,] = 0
    accept.count.chi[] = 0
    accept.count.theta[] = 0
  }
  
  #Delete variables not required
  if (b == num.runs) {
    rm(chi.gradient, chi.new, psi.new, psi.tilde.new, sense.probs.new, chi.gradient.new, MH, gamma1, 
       accept.rate, sense.probs.t, sense.probs.t.new, theta.prior.mean, theta.prior.var, theta.t.gradient, 
       theta.t.gradient.new, theta.t.new, psi.t.new, exp.psi.t.new, psi.t.tilde.new, phi.prior.var, 
       sense.probs.g.t, sense.probs.g.t.new, phi.g.t.new, exp.phi.g.t.new, phi.g.t.tilde.new, phi.g.t.gradient, 
       phi.g.t.gradient.new); gc()
  }
  
}#for b



#############################################################################
# MCMC results ------------------------------------------------------------ #
#############################################################################

#Acceptance rate
accept.count.phi / num.iterations
accept.count.chi / num.iterations
accept.count.theta / num.iterations


#Likelihood
log.lik.sim = colSums(log(apply(sense.probs.sim, 1, colSums)))
plot(log.lik.sim, type = "l")
mean(log.lik.sim)

num.samples = num.iterations/N
burn.in = 2000
#burn.in = 0
plot(log.lik.sim[(burn.in+1):num.samples], type = "l")
mean(log.lik.sim[(burn.in+1):num.samples])


#Trace plots - phi
Time = 8
Genre = 1
plot(as.mcmc(phi.tilde.sim[,,Genre,Time]))
plot(as.mcmc(phi.tilde.sim[(burn.in+1):num.samples,,Genre,Time]))


#Trace plots - psi
Time = 8
Sense = 1
Words = 1:3
plot(as.mcmc(psi.tilde.sim[(burn.in+1):num.samples,Words,Sense,Time]))

plot(as.mcmc(chi.sim[(burn.in+1):num.samples,Words,Sense]))
plot(as.mcmc(theta.sim[(burn.in+1):num.samples,Words,Time]))


# #Subtract means from all columns of chi.sim to get mean zero and add to all columns of theta.sim for identifiability
# chi.sim.means = apply(chi.sim, 1:2, mean)
# theta.sim.adj = sapply(1:num.periods, function(t) {theta.sim[,,t] + chi.sim.means}, simplify = "array")
# chi.sim.adj = sapply(1:num.senses, function(k) {chi.sim[,,k] - chi.sim.means}, simplify = "array")
# 
# # #Set column means of chi.sim.adj and theta.sim.adj to zero for identifiability
# # theta.sim.adj = sapply(1:num.periods, function(t) {theta.sim.adj[,,t] - rowMeans(theta.sim.adj[,,t])}, simplify = "array")
# # chi.sim.adj = sapply(1:num.senses, function(k) {chi.sim.adj[,,k] - rowMeans(chi.sim.adj[,,k])}, simplify = "array")
# 
# plot(as.mcmc(chi.sim.adj[(burn.in+1):num.samples,Words,Sense]))
# plot(as.mcmc(theta.sim.adj[(burn.in+1):num.samples,Words,Time]))



#############################################################################
# Load results ------------------------------------------------------------ #
#############################################################################

load("snippets.RData")
load("other.parameters.RData")
load("initial.z.RData")

b = 1
filename.suffix = paste((b-1)*num.iterations+1, b*num.iterations, sep = "-")

load(paste("psi.tilde.sim", filename.suffix, "RData", sep = "."))
load(paste("phi.tilde.sim", filename.suffix, "RData", sep = "."))
load(paste("chi.sim", filename.suffix, "RData", sep = "."))
load(paste("theta.sim", filename.suffix, "RData", sep = "."))
load(paste("sense.probs.sim", filename.suffix, "RData", sep = "."))
load(paste("run.time", filename.suffix, "RData", sep = "."))
load(paste("MALA.parameters.RData", filename.suffix, "RData", sep = "."))
