#Load all the packages and functions for EDiSC (MALA sampling), then run the following code

Vs = c(250, 500, 1000, 2000, 3000) #choices for num.words
Ds = c(250, 500, 1000, 2500, 5000, 7500, 10000) #choices for num.snippets
Ms = c(25, 200)
runs = c(1,2,3)

for (embedding.dim in Ms) {
  
  load(paste("/data/gnateater/sczafar/DPhil work (local)/Embeddings - COHA/COHA.word_vectors.min_count_10.dim_",embedding.dim,".RData", sep = ""))
  
  for (num.words in Vs) {
    
    for (num.snippets in Ds) {
      
      print((paste("M = ", embedding.dim, ", V = ", num.words, ", D = ", num.snippets, sep = "")), quote = FALSE)
      
      load(paste("snippets_V=",num.words,"_D=",num.snippets,".RData", sep=""))
      
      
      #word embeddings (rho)
      rho = word_vectors[paste(words.used),] 
      embedding.dim = ncol(rho)
      dimnames(rho) = list(Word = 1:num.words, Embedding = 1:embedding.dim)
      
      gc()
      
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
        
        
        #Set correction term (sigma) = 0, and initialise psi and xi based on the initial z
        sigma = numeric(num.words) 
        
        snippets.expanded = gather(cbind(snippets, sense = z), key = position, value = word, 1:snippet.length)
        
        word.counts = table(factor(snippets.expanded$word, levels = 1:num.words),
                            factor(snippets.expanded$sense, levels = 1:num.senses),
                            factor(snippets.expanded$Time, levels = 1:num.periods),
                            dnn = c("Word", "Sense", "Time"))
        
        psi.tilde = sapply(1:num.periods, function(t) {sapply(1:num.senses, function(k)
        {(word.counts[,k,t]+0.01)/(sum(word.counts[,k,t])+0.01*num.words)}, simplify = "aray" )}, simplify = "array" )
        dimnames(psi.tilde) = list(Word = 1:num.words, Sense = 1:num.senses, Time = 1:num.periods)
        
        psi = sapply(1:num.periods, function(t) {t(log(t(psi.tilde[,,t]) / psi.tilde[num.words,,t]))}, simplify = "array")
        psi = sapply(1:num.periods, function(t) {t(t(psi[,,t]) - colMeans(psi[,,t]))}, simplify = "array") #adjust psi for identifiability
        
        rho.svd = svd(rho)
        xi = sapply(1:num.periods, function(t) {rho.svd$v %*% diag(1/rho.svd$d) %*% t(rho.svd$u) %*% psi[,,t]}, simplify = "array") #best-fitting xi
        dimnames(xi) = list(Embedding = 1:embedding.dim, Sense = 1:num.senses, Time = 1:num.periods)
        
        #Fit a linear regression to estimate theta and chi. (We have MKT equations and M(K+T) variables)
        mkt.combos = expand.grid(dimnames(xi)) #embedding-sense-time combinations
        mkt.combos = data.frame(apply(mkt.combos,2,as.numeric)) #convert to numeric
        design.matrix.idx = data.frame(row.idx = rep(1:(embedding.dim*num.senses*num.periods),2),
                                       col.idx = c((mkt.combos$Sense-1)*embedding.dim+mkt.combos$Embedding,
                                                   embedding.dim*num.senses + (mkt.combos$Time-1)*embedding.dim+mkt.combos$Embedding))
        design.matrix = sparseMatrix(i = design.matrix.idx$row.idx, j = design.matrix.idx$col.idx, x = 1)
        design.matrix = design.matrix[,-((embedding.dim*(num.senses-1)+1):(embedding.dim*num.senses))] #remove last sense of chi for identifiability
        #design.matrix = cbind(1, design.matrix) #add intercept column
        xi.fit = MatrixModels:::lm.fit.sparse(design.matrix, as.vector(xi), method = "cholesky")
        
        chi = array(data = c(xi.fit$coef[1:(embedding.dim*(num.senses-1))], numeric(embedding.dim)), 
                    dim = c(embedding.dim, num.senses), dimnames = list(Embedding = 1:embedding.dim, Sense = 1:num.senses))
        
        theta = array(data = xi.fit$coef[-(1:(embedding.dim*(num.senses-1)))], 
                      dim = c(embedding.dim, num.periods),dimnames = list(Embedding = 1:embedding.dim, Time = 1:num.periods))
        
        #Adjust chi and theta for identifiability
        chi.mean = rowMeans(chi)
        theta = theta + chi.mean
        chi = chi - chi.mean
        theta = t(t(theta) - colMeans(theta))
        chi = t(t(chi) - colMeans(chi))
        
        #Re-calculate xi, psi and psi.tilde based on initial theta and chi
        xi = sapply(1:num.periods, function(t) {sapply(1:num.senses, function(k) {chi[,k]+theta[,t]} )}, simplify = "array")
        dimnames(xi) = list(Embedding = 1:embedding.dim, Sense = 1:num.senses, Time = 1:num.periods)
        psi = sapply(1:num.periods, function(t) {rho %*% xi[,,t]}, simplify = "array")
        dimnames(psi) = dimnames(psi.tilde)
        psi.tilde = aperm(aaply(exp(psi), 3, function(x) t(t(x)/colSums(x))), c(2,3,1))
        
        rm(sense.counts, snippets.expanded, word.counts, rho.svd, mkt.combos, design.matrix.idx, design.matrix, xi.fit, chi.mean); gc()
        
        
        #Set coefficient alpha in the model X(t) = alpha * X(t-1) + noise
        alpha.phi = 0.9
        alpha.theta = 0.9
        
        
        #Set variances for phi, chi, theta and sigma
        rho.squared.diffs = unlist(lapply(1:(nrow(rho)-1), function(x){as.numeric(colSums((rho[x,] - t(rho[-(1:x),,drop=FALSE]))^2))}))
        median.rho.squared.diff = median(rho.squared.diffs)
        kappa.phi = 0.25
        kappa.theta = 0.5 / median.rho.squared.diff
        kappa.chi = 2.5 / median.rho.squared.diff
        kappa.sigma = 0.25
        
        rm(rho.squared.diffs); gc()
        
        
        #Set up structures for prior means (to avoid errors later for phi when num.genres = 1, and for better efficiency)
        phi.prior.mean = array(data = 0, dim = c(num.senses, num.genres), dimnames = list(Sense = 1:num.senses, Genre = 1:num.genres))
        chi.prior.mean = 0 #this remains fixed throughout
        sigma.prior.mean = 0 #this remains fixed throughout
        
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
        
        # #Get the total counts for each word in each time period
        # word.counts.by.time = sapply(1:num.periods, function(t){rowSums(word.counts.in.snippets[,SnippetIDs.by.time[[t]]])}) 
        # dimnames(word.counts.by.time) = list(Word = 1:num.words, Time = 1:num.periods)
        
        #Get the total counts for each word across all snippets
        word.counts.totals = rowSums(word.counts.in.snippets)
        
        #Lemmatise fixed values that will be called repeatedly
        rho.times.word.counts = t(rho) %*% word.counts.in.snippets
        
        
        #Get initial probabilities for all snippets belonging to each sense
        sense.probs = array(dim = c(num.senses, num.snippets),
                            dimnames = list(Sense = 1:num.senses, SnippetID = snippets$SnippetID))
        for(t in 1:num.periods) {
          for(g in 1:num.genres) {
            sense.probs[,SnippetIDs[[t]][[g]]] = apply(snippets[SnippetIDs[[t]][[g]], 1:snippet.length], 1, get.sense.prob,
                                                       phi.tilde[,g,t], psi.tilde[,,t])
          }#for
        }#for
        
        
        #Set scales (variances) for MALA proposals
        phi.MALA.prop.scale = array(2.4^2 / num.senses, dim = c(num.genres, num.periods), 
                                    dimnames = list(Genre = 1:num.genres, Time = 1:num.periods) )
        
        chi.MALA.prop.scale = 2.4^2 / (embedding.dim * num.senses)^2
        
        theta.MALA.prop.scale = array(2.4^2 / embedding.dim^2, dim = c(num.periods), dimnames = list(Time = 1:num.periods))
        
        sigma.MALA.prop.scale = 2.4^2 / num.words
        
        
        #Set MALA tuning parameters - refer to Shaby and Wells (2010)
        target.accept.rate.phi = 0.574
        target.accept.rate.chi = 0.574
        target.accept.rate.theta = 0.574
        target.accept.rate.sigma = 0.574
        block.size = 10 #k in the paper
        
        #c0.phi = 1
        c1.phi = 0.8
        
        #c0.chi = 1
        c1.chi = 0.8
        
        #c0.theta = 1
        c1.theta = 0.8
        
        #c0.sigma = 1
        c1.sigma = 0.8
        
        
        #Set up arrays to hold acceptance counts - for MALA tuning
        phi.jumps = array(0, dim = c(block.size, num.genres, num.periods), dimnames = 
                            list(Iteration = 1:block.size, Genre = 1:num.genres, Time = 1:num.periods))
        
        chi.jumps = array(0, dim = c(block.size), dimnames = list(Iteration = 1:block.size))
        
        theta.jumps = array(0, dim = c(block.size, num.periods), dimnames = 
                              list(Iteration = 1:block.size, Time = 1:num.periods))
        
        sigma.jumps = array(0, dim = c(block.size), dimnames = list(Iteration = 1:block.size))
        
        
        #Set up arrays to hold acceptance counts - for acceptance rate calculation
        accept.count.phi = array(0, dim = c(num.genres, num.periods), 
                                 dimnames = list(Genre = 1:num.genres, Time = 1:num.periods)) #To calculate acceptance rate for phi for each time period and genre
        
        accept.count.chi = 0 #To calculate acceptance rate for chi
        
        accept.count.theta = array(0, dim = c(num.periods), 
                                   dimnames = list(Sense = 1:num.periods)) #To calculate acceptance rate for theta for each time period
        
        accept.count.sigma = 0 #To calculate acceptance rate for sigma
        
        
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
        
        chi.sim = array(dim = c(num.iterations/N, embedding.dim, num.senses), 
                        dimnames = list(Iteration = 1:(num.iterations/N), Embedding = 1:embedding.dim, Sense = 1:num.senses))
        
        theta.sim = array(dim = c(num.iterations/N, embedding.dim, num.periods), 
                          dimnames = list(Iteration = 1:(num.iterations/N), Embedding = 1:embedding.dim, Time = 1:num.periods))
        
        sigma.sim = array(dim = c(num.iterations/N, num.words), 
                          dimnames = list(Iteration = 1:(num.iterations/N), Word = 1:num.words))
        
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
                                                 rho.times.word.counts, snippet.lengths)
            
            #Propose new chi
            chi.new = propose.new.MALA(chi, embedding.dim*num.senses, chi.gradient, chi.MALA.prop.scale)
            
            #Calculate new xi, psi and psi.tilde
            xi.new = sapply(1:num.periods, function(t) {sapply(1:num.senses, function(k) {
              chi.new[,k]+theta[,t]} )}, simplify = "array")
            psi.new = sapply(1:num.periods, function(t) {rho %*% xi.new[,,t] + sigma}, simplify = "array")
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
                                                     rho.times.word.counts, snippet.lengths)
            
            #Calculate MH ratio
            MH = get.mh.ratio(chi, chi.new, chi.prior.mean, kappa.chi, sense.probs, sense.probs.new, chi.gradient, chi.gradient.new,
                              chi.MALA.prop.scale)
            
            #Accept or reject
            if(log(runif(1)) < MH) {
              accept.count.chi = accept.count.chi + 1
              chi.jumps[i %% block.size + 1] = 1
              sense.probs = sense.probs.new
              chi = chi.new
              xi = xi.new
              psi = psi.new
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
                                                         sense.probs.t, rho.times.word.counts[,SnippetIDs.by.time[[t]]],
                                                         snippet.lengths[SnippetIDs.by.time[[t]]])
              
              #Propose new theta
              theta.t.new = propose.new.MALA(theta[,t], embedding.dim, theta.t.gradient, theta.MALA.prop.scale[t])
              
              #Calculate new xi, psi and psi.tilde
              xi.t.new = theta.t.new + chi
              psi.t.new = rho %*% xi.t.new + sigma
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
                                                             sense.probs.t.new, rho.times.word.counts[,SnippetIDs.by.time[[t]]],
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
                xi[,,t] = xi.t.new
                psi[,,t] = psi.t.new
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
            
            
            # --- sigma updates --- #
            
            #Get sigma gradient
            sigma.gradient = get.sigma.MALA.gradient(sigma, psi.tilde, sigma.prior.mean, kappa.sigma, sense.probs,
                                                     word.counts.totals, snippet.lengths)
            
            #Propose new sigma
            sigma.new = propose.new.MALA(sigma, num.words, sigma.gradient, sigma.MALA.prop.scale)
            
            #Calculate new psi and psi.tilde
            psi.new = sapply(1:num.periods, function(t) {psi[,,t] - sigma + sigma.new}, simplify = "array")
            psi.tilde.new = aperm(aaply(exp(psi.new), 3, function(x) t(t(x)/colSums(x))), c(2,3,1))
            
            #Get probabilities for all snippets belonging to each sense
            sense.probs.new = sense.probs
            for(t in 1:num.periods) {
              for(g in 1:num.genres) {
                sense.probs.new[,SnippetIDs[[t]][[g]]] = apply(snippets[SnippetIDs[[t]][[g]], 1:snippet.length], 1, get.sense.prob,
                                                               phi.tilde[,g,t], psi.tilde.new[,,t])
              }#for
            }#for
            
            #Get new sigma gradient
            sigma.gradient.new = get.sigma.MALA.gradient(sigma.new, psi.tilde.new, sigma.prior.mean, kappa.sigma, sense.probs.new,
                                                         word.counts.totals, snippet.lengths)
            
            #Calculate MH ratio
            MH = get.mh.ratio(sigma, sigma.new, sigma.prior.mean, kappa.sigma, sense.probs, sense.probs.new, 
                              sigma.gradient, sigma.gradient.new, sigma.MALA.prop.scale)
            
            #Accept or reject
            if(log(runif(1)) < MH) {
              accept.count.sigma = accept.count.sigma + 1
              sigma.jumps[i %% block.size + 1] = 1
              sense.probs = sense.probs.new
              sigma = sigma.new
              psi = psi.new
              psi.tilde = psi.tilde.new
            } else {
              sigma.jumps[i %% block.size + 1] = 0
            }#else
            
            #Update MALA parameters
            if (b <= stop.MALA.tuning.after) {
              j = (b-1)*num.iterations + i
              if (j >= block.size) {
                gamma1 = ((j+1)/block.size)^(-c1.sigma)
                #gamma2 = c0.sigma * gamma1
                
                accept.rate = sum(sigma.jumps)/block.size
                sigma.MALA.prop.scale = update.scale(sigma.MALA.prop.scale, accept.rate, gamma1, target.accept.rate.sigma)
              }#if
            }#if
            
            
            #Store values
            if (i%%N == 0) {
              sense.probs.sim[i/N,,] = sense.probs
              phi.tilde.sim[i/N,,,] = phi.tilde
              psi.tilde.sim[i/N,,,] = psi.tilde
              theta.sim[i/N,,] = theta
              chi.sim[i/N,,] = chi
              sigma.sim[i/N,] = sigma
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
          
          save(Start.Time, End.Time, Run.Time, file = paste("run.time_EDiSC_M=",embedding.dim,"_sims=",num.iterations,
                                                            "_V=",num.words,"_D=",num.snippets,"_run",run,".RData", sep=""))
          
          #Reset accpet counts
          if (b != num.runs) {
            accept.count.phi[,] = 0
            accept.count.chi[] = 0
            accept.count.theta[] = 0
            accept.count.sigma[] = 0
          }
          
          #Delete variables not required
          if (b == num.runs) {
            rm(chi.gradient, chi.new, xi.new, psi.new, psi.tilde.new, sense.probs.new, chi.gradient.new, MH, gamma1, 
               accept.rate, sense.probs.t, sense.probs.t.new, theta.prior.mean, theta.prior.var, theta.t.gradient, 
               theta.t.gradient.new, theta.t.new, xi.t.new, psi.t.new, exp.psi.t.new, psi.t.tilde.new, phi.prior.var, 
               sense.probs.g.t, sense.probs.g.t.new, phi.g.t.new, exp.phi.g.t.new, phi.g.t.tilde.new, phi.g.t.gradient, 
               phi.g.t.gradient.new, sigma.gradient, sigma.new, sigma.gradient.new); gc()
          }
          
        }#for b
        
      }#for run
      
    }#for num.snippets
    
  }#for num.words
  
}#for embedding.dim


