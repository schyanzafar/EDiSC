#Load all the packages and functions for DiSC (MALA sampling), then run the following code

Vs = c(250, 500, 1000, 2000, 3000) #choices for num.words
Ds = c(250, 500, 1000, 2500, 5000, 7500, 10000) #choices for num.snippets
runs = c(1,2,3)

for (num.words in Vs) {
  
  for (num.snippets in Ds) {
    
    print((paste("V = ", num.words, ", D = ", num.snippets, sep = "")), quote = FALSE)
    
    load(paste("snippets_V=",num.words,"_D=",num.snippets,".RData", sep=""))
    
    for (run in runs) {
      
      print((paste("run", run)), quote = FALSE)
      
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
      num.iterations = 500
      
      #Sample every Nth iteration
      N = 1
      
      #Specify number of times to run the above MCMC iterations
      num.runs = 1
      
      #Stop tuning MALA scale parameters after this many runs
      stop.MALA.tuning.after = 99
      
      
      #Set up arrays to hold simulated values
      phi.tilde.sim = array(dim = c(num.iterations/N, num.senses, num.genres, num.periods),
                            dimnames = list(Iteration = 1:(num.iterations/N), Sense = 1:num.senses, 
                                            Genre = 1:num.genres, Time = 1:num.periods))
      
      psi.tilde.sim = array(dim = c(num.iterations/N, num.words, num.senses, num.periods),
                            dimnames = list(Iteration = 1:(num.iterations/N), 
                                            Word = 1:num.words, Sense = 1:num.senses, Time = 1:num.periods))
      
      chi.sim = array(dim = c(num.iterations/N, num.words, num.senses), 
                      dimnames = list(Iteration = 1:(num.iterations/N), Word = 1:num.words, Sense = 1:num.senses))
      
      theta.sim = array(dim = c(num.iterations/N, num.words, num.periods), 
                        dimnames = list(Iteration = 1:(num.iterations/N), Word = 1:num.words, Time = 1:num.periods))
      
      sense.probs.sim = array(dim = c(num.iterations/N, num.senses, num.snippets),
                              dimnames = list(Iteration = 1:(num.iterations/N), 
                                              Sense  = 1:num.senses, SnippetID = snippets$SnippetID))
      
      
      
      ### ------------------------------ Run MCMC ------------------------------ ###
      
      progress.bar = txtProgressBar(min = 1 , max = num.iterations, style = 3) #Show progress
      
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
        
        save(Start.Time, End.Time, Run.Time, file = paste("run.time_DiSC_sims=",num.iterations,"_V=",num.words,
                                                          "_D=",num.snippets,"_run",run,".RData", sep=""))
        
        #Reset accept counts
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
      
    }#for run
    
  }#for num.snippets
  
}#for num.words
