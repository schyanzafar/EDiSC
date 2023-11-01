#############################################################################
# Setup ------------------------------------------------------------------- #
#############################################################################

# required packages
library(foreach)
library(Matrix)
library(text2vec)

# parallel setup
library(doParallel)
num.cores = 10
cl = makeCluster(num.cores, type = "FORK")
registerDoParallel(cl)
# stopCluster(cl)

setwd(".../Greek data")
load("words.RData")
lemmas.unique = read.csv("lemmas.unique.with.SW.tags.csv", header = TRUE, colClasses = "character")



#############################################################################
# Tokens ------------------------------------------------------------------ #
#############################################################################

# Create iterator over tokens
# replace stopwords with fullstops and collect lemmaIDs
stopword.lemmaIDs = lemmas.unique$lemmaID[lemmas.unique$SW == TRUE]
stopword.lemmaIDs = c(stopword.lemmaIDs, "unknown")
tokens = foreach(i = 1:length(words)) %do% {
  ifelse(words[[i]]$lemmaID %in% stopword.lemmaIDs, ".",words[[i]]$lemmaID)
}

rm(lemmas.unique, stopword.lemmaIDs, words, i); gc()

# Save tokens for later use
setwd(".../Embeddings")
save(tokens, file = "Greek_data.tokens.RData")
load("Greek_data.tokens.RData")



#############################################################################
# GloVe embeddings -------------------------------------------------------- #
#############################################################################

# Following code copied and adapted from https://cran.r-project.org/web/packages/text2vec/vignettes/glove.html on 1 Jun 2022

# Create vocabulary. Terms will be lemmas
it = itoken(tokens)
vocab = create_vocabulary(it)
vocab = vocab[-which(vocab$term=="."),] #remove stopwords (fullstops) from vocabulary
min_count = 10L #remove words appearing fewer than this many times
vocab = prune_vocabulary(vocab, term_count_min = min_count) 

# Use our filtered vocabulary
vectorizer = vocab_vectorizer(vocab)

# term co-occurence matric (TCM)
context.window = 10L
tcm = create_tcm(it, vectorizer, skip_grams_window = context.window, skip_grams_window_context = "symmetric") 

# Save TCM for later use
setwd(".../Embeddings")
save(tcm, file = "Greek_data.tcm.min_count_10.RData")
load("Greek_data.tcm.min_count_10.RData")

# GloVe embeddings
set.seed(100) #for master process
clusterSetRNGStream(cl, 500) #set seed for slave processes
embedding.dim = 100 #dimension of embedding vectors 
glove = GlobalVectors$new(rank = embedding.dim, x_max = 100, learning_rate = 0.05)
wv_main = glove$fit_transform(tcm, n_iter = 50, convergence_tol = 0.01, n_threads = num.cores)

# word vectors
wv_context = glove$components
word_vectors = wv_main + t(wv_context)

# Save word vectors for later use
setwd(".../Embeddings")
save(word_vectors, file = "Greek_data.word_vectors.min_count_10.dim_100.RData")
load("Greek_data.word_vectors.min_count_10.dim_100.RData")
