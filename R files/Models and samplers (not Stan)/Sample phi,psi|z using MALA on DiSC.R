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

#Write a function to calculate the MALA gradient for phi updates for a single time slice
get.phi.MALA.gradient = function(phi, phi.tilde, phi.prior.mean, phi.prior.var, counts, counts.total) {
  
  - (phi - phi.prior.mean) / phi.prior.var + counts - counts.total * phi.tilde
  
}#get.phi.MALA.gradient


#Write a function to calculate the MALA gradient for chi updates for all senses
get.chi.MALA.gradient = function(chi, psi.tilde, chi.prior.mean, chi.prior.var, counts, counts.total) {
  
  - (chi - chi.prior.mean) / chi.prior.var + counts - t(apply(psi.tilde, 1, function(x){rowSums(x * counts.total)}))
  
}#get.chi.MALA.gradient


#Write a function to calculate the MALA gradient for theta updates for a single time slice
get.theta.MALA.gradient = function(theta, psi.tilde, theta.prior.mean, theta.prior.var, counts, counts.total) {
  
  - (theta - theta.prior.mean) / theta.prior.var + counts - as.numeric(psi.tilde %*% counts.total)
  
}#get.theta.MALA.gradient


#Write a function to generate MALA proposals for phi, chi and theta given given other variables
propose.new.MALA = function(X, size, gradient, prop.scale) {
  
  X + rnorm(size, prop.scale/2 * gradient, (prop.scale)^0.5)
  
}#propose.new.MALA


#Write a function to calcualte the log MH ratio given current and proposed phi, chi and theta
get.mh.ratio = function(X, X.new, X.prior.mean, X.prior.var, log.lik, log.lik.new, gradient, gradient.new, prop.scale) {
  
  sum(dnorm(X.new, X.prior.mean, sqrt(X.prior.var), log = TRUE) - 
        dnorm(X, X.prior.mean, sqrt(X.prior.var), log = TRUE)) +
    
    sum(dnorm(X, X.new + prop.scale/2 * gradient.new, sqrt(prop.scale), log = TRUE) -
          dnorm(X.new, X + prop.scale/2 * gradient, sqrt(prop.scale), log = TRUE)) +
    
    log.lik.new - log.lik
  
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

#Index of non-empty snippets -- Greek data only
# idx = which(snippet.lengths > 0)
# idx = which(snippet.lengths > 0 & snippets.info$disambiguation == 1) #of type collocates
sense.ids = levels(factor(snippets.info$sense.id))
idx = which(snippets.info$sense.id %in% sense.ids[1:3])
names(idx) = idx 
rm(sense.ids)

# #Index of snippets with either the riverbank or the institution bank true sense -- bank data only
# sense.ids = levels(snippets.info$sense.id)
# idx = which(snippets.info$sense.id %in% sense.ids[1:2])
# names(idx) = idx #required for consistency with the Greek data
# rm(sense.ids)


#Initialise parameters based on true z
z = as.numeric(factor(snippets.info$sense.id))

#Get the counts of snippets in each time period and genre
num.senses = length(unique(z[idx]))
snippet.counts = table(factor(snippets$genre[idx], levels = 1:num.genres),
                       factor(snippets$Time[idx], levels = 1:num.periods),
                       dnn = c("Genre", "Time"))

#Initialise phi based on the initial z
sense.counts = table(factor(z[idx], levels = 1:num.senses), factor(snippets$genre[idx], levels = 1:num.genres), 
                     factor(snippets$Time[idx], levels = 1:num.periods), dnn = c("Sense", "Genre", "Time"))

phi.tilde = aperm(sapply(1:num.senses, function(k) {(sense.counts[k,,] + 0.01) / (snippet.counts + 0.01*num.senses)}, simplify = "array"), c(3,1,2))
dimnames(phi.tilde) = list(Sense = 1:num.senses, Genre = 1:num.genres, Time = 1:num.periods)

phi = sapply(1:num.periods, function(t) {t(log(t(phi.tilde[,,t]) / phi.tilde[num.senses,,t]))}, simplify = "array")
dimnames(phi) = dimnames(phi.tilde)

#Adjust phi for identifiability
phi = aperm(sapply(1:num.genres, function(g) {t(t(phi[,g,]) - colMeans(phi[,g,]))}, simplify = "array"), c(1,3,2))
dimnames(phi) = dimnames(phi.tilde)


#Initialise psi based on the initial z
snippets.expanded = gather(cbind(snippets[idx,], sense = z[idx]), key = position, value = word, 1:snippet.length)
# snippets.expanded = pivot_longer(cbind(snippets[idx,], sense = z[idx]), 1:snippet.length, names_to = "position", values_to = "word")
word.counts = table(factor(snippets.expanded$word, levels = 1:num.words),
                    factor(snippets.expanded$sense, levels = 1:num.senses),
                    factor(snippets.expanded$Time, levels = 1:num.periods),
                    dnn = c("Word", "Sense", "Time"))

psi.tilde = sapply(1:num.periods, function(t) {sapply(1:num.senses, function(k)
  {(word.counts[,k,t]+0.01)/(sum(word.counts[,k,t])+0.01*num.words)}, simplify = "aray" )}, simplify = "array" )
dimnames(psi.tilde) = list(Word = 1:num.words, Sense = 1:num.senses, Time = 1:num.periods)

psi = sapply(1:num.periods, function(t) {t(log(t(psi.tilde[,,t]) / psi.tilde[num.words,,t]))}, simplify = "array")
dimnames(psi) = dimnames(psi.tilde)

rm(snippets.expanded); gc()

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


#Memoise word count totals
word.counts.KT = apply(word.counts, 2:3, sum)
word.counts.VK = apply(word.counts, 1:2, sum)
word.counts.VT = apply(word.counts, c(1,3), sum)


#Set coefficient alpha in the model X(t) = alpha * X(t-1) + noise
alpha.phi = 0.9
alpha.theta = 0.9


#Set (initial) variances for phi, chi and theta
kappa.phi = 0.25
kappa.theta = 0.25
kappa.chi = 1.25


#Set up structures for prior means (to avoid errors later for phi when num.genres = 1)
phi.prior.mean = array(data = 0, dim = c(num.senses, num.genres), dimnames = list(Sense = 1:num.senses, Genre = 1:num.genres))
chi.prior.mean = 0 #this remains fixed throughout


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
                         dimnames = list(Sense = 1:num.periods)) #To calculate acceptance rate for theta for each time period


#We will alternately go forward and backward in time for each successive iteration
time.direction = 1:num.periods


#Specify number of MCMC iterations
num.iterations = 2000

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


### ------------------------------ Run MCMC ------------------------------ ###

progress.bar = txtProgressBar(min = 1 , max = num.iterations, style = 3) #Show progress

random.seed = 500
set.seed(random.seed)

save(snippets, snippets.info, words.used, num.words, num.snippets, snippet.length, snippet.lengths, num.periods, num.genres, num.senses, idx, file = "snippets.RData")
save(z, file = "true.z.RData")
save(num.iterations, N, kappa.phi, kappa.chi, kappa.theta, alpha.phi, alpha.theta, random.seed, file = "other.parameters.RData")

for (b in 1:num.runs) {
  
  Start.Time = Sys.time() #Save start time
  for (i in 1:num.iterations) {
    j = (b-1)*num.iterations + i
    
    # --- chi updates --- #

    #Get chi gradient
    chi.gradient = get.chi.MALA.gradient(chi, psi.tilde, chi.prior.mean, kappa.chi, word.counts.VK, word.counts.KT)
    chi.log.lik = sum(word.counts * log(psi.tilde))

    #Propose new chi
    chi.new = propose.new.MALA(chi, num.words*num.senses, chi.gradient, chi.MALA.prop.scale)
    
    #Calculate new psi
    psi.new = sapply(1:num.periods, function(t) {sapply(1:num.senses, function(k) {
      chi.new[,k]+theta[,t]} )}, simplify = "array")
    psi.tilde.new = aperm(aaply(exp(psi.new), 3, function(x) t(t(x)/colSums(x))), c(2,3,1))

    #Get new chi gradient
    chi.gradient.new = get.chi.MALA.gradient(chi.new, psi.tilde.new, chi.prior.mean, kappa.chi, word.counts.VK, word.counts.KT)
    chi.log.lik.new = sum(word.counts * log(psi.tilde.new))

    #Calculate MH ratio
    MH = get.mh.ratio(chi, chi.new, chi.prior.mean, kappa.chi, chi.log.lik, chi.log.lik.new, chi.gradient, chi.gradient.new,
                      chi.MALA.prop.scale)

    #Accept or reject
    if(log(runif(1)) < MH) {
      accept.count.chi = accept.count.chi + 1
      chi.jumps[i %% block.size + 1] = 1
      chi = chi.new
      psi.tilde = psi.tilde.new
    } else {
      chi.jumps[i %% block.size + 1] = 0
    }#else

    #Update MALA parameters
    if (b <= stop.MALA.tuning.after) {
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
                                                 word.counts.VT[,t], word.counts.KT[,t])
      theta.log.lik = sum(word.counts[,,t] * log(psi.tilde[,,t]))
      
      #Propose new theta
      theta.t.new = propose.new.MALA(theta[,t], num.words, theta.t.gradient, theta.MALA.prop.scale[t])
      
      #Calculate new psi
      psi.t.new = theta.t.new + chi
      exp.psi.t.new = exp(psi.t.new)
      psi.t.tilde.new = t(t(exp.psi.t.new) / colSums(exp.psi.t.new))

      #Get new theta gradient
      theta.t.gradient.new = get.theta.MALA.gradient(theta.t.new, psi.t.tilde.new, theta.prior.mean, theta.prior.var,
                                                     word.counts.VT[,t], word.counts.KT[,t])
      theta.log.lik.new = sum(word.counts[,,t] * log(psi.t.tilde.new))

      #Calculate MH ratio
      MH = get.mh.ratio(theta[,t], theta.t.new, theta.prior.mean, theta.prior.var, theta.log.lik, theta.log.lik.new, 
                        theta.t.gradient, theta.t.gradient.new, theta.MALA.prop.scale[t])

      #Accept or reject
      if(log(runif(1)) < MH) {
        accept.count.theta[t] = accept.count.theta[t] + 1
        theta.jumps[i %% block.size + 1, t] = 1
        theta[,t] = theta.t.new
        psi.tilde[,,t] = psi.t.tilde.new
      } else {
        theta.jumps[i %% block.size + 1, t] = 0
      }#else

      #Update MALA parameters
      if (b <= stop.MALA.tuning.after) {
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

        #Get phi gradient
        phi.g.t.gradient = get.phi.MALA.gradient(phi[,g,t], phi.tilde[,g,t], phi.prior.mean[,g], phi.prior.var, 
                                                 sense.counts[,g,t], snippet.counts[g,t])
        phi.log.lik = sum(sense.counts[,g,t] * log(phi.tilde[,g,t]))

        #Propose new phi
        phi.g.t.new = propose.new.MALA(phi[,g,t], num.senses, phi.g.t.gradient, phi.MALA.prop.scale[g,t])
        exp.phi.g.t.new = exp(phi.g.t.new)
        phi.g.t.tilde.new = exp.phi.g.t.new / sum(exp.phi.g.t.new)

        #Get new phi gradient
        phi.g.t.gradient.new = get.phi.MALA.gradient(phi.g.t.new, phi.g.t.tilde.new, phi.prior.mean[,g], 
                                                     phi.prior.var, sense.counts[,g,t], snippet.counts[g,t])
        phi.log.lik.new = sum(sense.counts[,g,t] * log(phi.g.t.tilde.new))

        #Calculate MH ratio
        MH = get.mh.ratio(phi[,g,t], phi.g.t.new, phi.prior.mean[,g], phi.prior.var, phi.log.lik, phi.log.lik.new, 
                          phi.g.t.gradient, phi.g.t.gradient.new, phi.MALA.prop.scale[g,t])

        #Accept or reject
        if(log(runif(1)) < MH) {
          accept.count.phi[g,t] = accept.count.phi[g,t] + 1
          phi.jumps[i %% block.size + 1, g, t] = 1
          phi[,g,t] = phi.g.t.new
          phi.tilde[,g,t] = phi.g.t.tilde.new
        } else {
          phi.jumps[i %% block.size + 1, g, t] = 0
        }#else

        #Update MALA parameters
        if (b <= stop.MALA.tuning.after) {
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
    rm(chi.gradient, chi.new, psi.new, psi.tilde.new, chi.gradient.new, MH, gamma1, accept.rate, 
       theta.prior.mean, theta.prior.var, theta.t.gradient, theta.t.gradient.new, theta.t.new, 
       psi.t.new, exp.psi.t.new, psi.t.tilde.new, phi.prior.var, phi.g.t.new, exp.phi.g.t.new, 
       phi.g.t.tilde.new, phi.g.t.gradient, phi.g.t.gradient.new); gc()
  }
  
}#for b



#############################################################################
# MCMC results ------------------------------------------------------------ #
#############################################################################

#Acceptance rate
accept.count.phi / num.iterations
accept.count.chi / num.iterations
accept.count.theta / num.iterations


#Trace plots - phi
Time = 8
Genre = 2
plot(as.mcmc(phi.tilde.sim[,,Genre,Time]))


#Trace plots - psi
Time = 8
Sense = 1
Words = 1:3
plot(as.mcmc(psi.tilde.sim[,Words,Sense,Time]))

plot(as.mcmc(chi.sim[,Words,Sense]))
plot(as.mcmc(theta.sim[,Words,Time]))


# #Subtract means from all columns of chi.sim to get mean zero and add to all columns of theta.sim for identifiability
# chi.sim.means = apply(chi.sim, 1:2, mean)
# theta.sim.adj = sapply(1:num.periods, function(t) {theta.sim[,,t] + chi.sim.means}, simplify = "array")
# chi.sim.adj = sapply(1:num.senses, function(k) {chi.sim[,,k] - chi.sim.means}, simplify = "array")
# 
# # #Set column means of chi.sim.adj and theta.sim.adj to zero for identifiability
# # theta.sim.adj = sapply(1:num.periods, function(t) {theta.sim.adj[,,t] - rowMeans(theta.sim.adj[,,t])}, simplify = "array")
# # chi.sim.adj = sapply(1:num.senses, function(k) {chi.sim.adj[,,k] - rowMeans(chi.sim.adj[,,k])}, simplify = "array")
# 
# plot(as.mcmc(chi.sim.adj[-(1:burn.in),Words,Sense]))
# plot(as.mcmc(theta.sim.adj[-(1:burn.in),Words,Time]))



#############################################################################
# Load results ------------------------------------------------------------ #
#############################################################################

load("snippets.RData")
load("other.parameters.RData")
load("true.z.RData")

b = 1
filename.suffix = paste((b-1)*num.iterations+1, b*num.iterations, sep = "-")

load(paste("psi.tilde.sim", filename.suffix, "RData", sep = "."))
load(paste("phi.tilde.sim", filename.suffix, "RData", sep = "."))
load(paste("chi.sim", filename.suffix, "RData", sep = "."))
load(paste("theta.sim", filename.suffix, "RData", sep = "."))
load(paste("run.time", filename.suffix, "RData", sep = "."))
load(paste("MALA.parameters.RData", filename.suffix, "RData", sep = "."))
