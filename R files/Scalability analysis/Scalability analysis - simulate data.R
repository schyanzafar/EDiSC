options(digits = 4)

#############################################################################
# Get vocabulary ---------------------------------------------------------- #
#############################################################################

#Load COHA data
load(".../Embeddings/COHA.tokens.RData")

library(text2vec)
it = itoken(tokens)
vocab = create_vocabulary(it)
vocab = vocab[-which(vocab$term=="."),] #remove stopwords (fullstops) from vocabulary
min_count = 10L #remove words appearing fewer than 10 times
vocab = prune_vocabulary(vocab, term_count_min = min_count) 

#retain only the top 90% of the vocabulary
vocab = vocab[order(vocab$term_count, decreasing = TRUE),]
cutoff = Position(function(x) x >= 0.9, cumsum(vocab$term_count)/sum(vocab$term_count))
vocab = vocab$term[1:cutoff]

save(vocab, file = "vocab.RData")
load("vocab.RData")

rm(tokens, it, cutoff, min_count); gc()



#############################################################################
# Set parameters for synthetic data simulation ---------------------------- #
#############################################################################

Vs = c(250, 500, 1000, 2000, 3000, 5000) #choices for num.words
Ds = c(250, 500, 1000, 2500, 5000, 7500, 10000) #choices for num.snippets

#Set parameters to mimic the 'bank' data
num.genres = 1
num.senses = 2
num.periods = 10
snippet.length = 14
num.words.per.snippet = 4 #fixed, to control variation in run times due to randomness



#############################################################################
# Simulate data ----------------------------------------------------------- #
#############################################################################

for (num.words in Vs) {
  set.seed(100)
  words.used = sample(vocab, num.words)
  
  #all senses have equal probability - irrelevant for checking computational speed
  true.phi.tilde = array(1/num.senses, dim = c(num.senses, num.genres, num.periods),
                         dimnames = list(Sense = 1:num.senses, Genre = 1:num.genres, Time = 1:num.periods))
  
  #all words have equal probability - irrelevant for checking computational speed
  true.psi.tilde = array(1/num.words, dim = c(num.words, num.senses, num.periods),
                         dimnames = list(Word = 1:num.words, Sense = 1:num.senses, Time = 1:num.periods))
  
  for (num.snippets in Ds) {
    
    #Set up data frame and vectors to hold data
    snippets = matrix(nrow = num.snippets, ncol = snippet.length + 3)
    colnames(snippets) = c(paste("V", 1:snippet.length, sep = ""), "Time", "SnippetID", "genre")
    snippets = as.data.frame(snippets, stringsAsFactors=FALSE)
    snippets$Time = rep(1:num.periods)
    snippets$SnippetID = 1:num.snippets
    snippets$genre = 1
    
    true.z = numeric(num.snippets)
    snippet.lengths = rep(num.words.per.snippet, num.snippets)
    
    #Simulate data
    set.seed(500)
    
    for (d in 1:num.snippets) {
      #Sample sense
      true.z[d] = sample(1:num.senses, size = 1, prob = true.phi.tilde[,snippets$genre[d],snippets$Time[d]])
      
      #Sample words
      snippets[d, sample(1:snippet.length, num.words.per.snippet)] = 
        sample(1:num.words, size = num.words.per.snippet, replace = FALSE, prob = true.psi.tilde[, true.z[d], snippets$Time[d]])
      
    }#for d
    
    #Write data to file
    save(snippets, num.words, num.snippets, snippet.length, snippet.lengths, num.periods,
         num.genres, num.senses, words.used, true.z, true.phi.tilde, true.psi.tilde,
         file = paste("snippets_V=",num.words,"_D=",num.snippets,".RData", sep=""))
    
  }#for num.snippets
  
}#for num.words



