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

setwd(".../COHA data processed")
load("texts.combined.RData")

#replace stopwords (blanks) with fullstops for convenience
texts.combined$lemma = ifelse(texts.combined$lemma=="",".",texts.combined$lemma) 



#############################################################################
# Tokens ------------------------------------------------------------------ #
#############################################################################

# Create iterator over tokens
textIDs = unique(texts.combined$textID); textIDs = textIDs[!is.na(textIDs)]
text.word.counts = table(texts.combined$textID)[paste(textIDs)]
text.indices.start = match(textIDs, texts.combined$textID)
text.indices.end = text.indices.start + as.integer(text.word.counts) - 1L

tokens = foreach(i = 1:length(textIDs)) %do% {
  paste(texts.combined$lemma[text.indices.start[i]:text.indices.end[i]])
}

rm(texts.combined, textIDs, text.indices.start, text.indices.end, text.word.counts, i); gc()

# Save tokens for later use
setwd(".../Embeddings")
save(tokens, file = "COHA.tokens.RData")
load("COHA.tokens.RData")



#############################################################################
# GloVe embeddings -------------------------------------------------------- #
#############################################################################

# Following code copied and adapted from https://cran.r-project.org/web/packages/text2vec/vignettes/glove.html on 1 Jun 2022

# Create vocabulary. Terms will be lemmas
it = itoken(tokens)
vocab = create_vocabulary(it)
vocab = vocab[-which(vocab$term=="."),] #remove stopwords (fullstops) from vocabulary
min_count = 10L #remove words appearing fewer than 10 times
vocab = prune_vocabulary(vocab, term_count_min = min_count) 

# Use our filtered vocabulary
vectorizer = vocab_vectorizer(vocab)

# term co-occurence matric (TCM)
context.window = 10L
tcm = create_tcm(it, vectorizer, skip_grams_window = context.window, skip_grams_window_context = "symmetric") 

# Save TCM for later use
setwd(".../Embeddings")
save(tcm, file = "COHA.tcm.min_count_10.RData")
load("COHA.tcm.min_count_10.RData")

# GloVe embeddings
set.seed(100) #for master process
clusterSetRNGStream(cl, 500) #set seed for slave processes
embedding.dim = 200 #dimension of embedding vectors 
glove = GlobalVectors$new(rank = embedding.dim, x_max = 100, learning_rate = 0.05)
wv_main = glove$fit_transform(tcm, n_iter = 50, convergence_tol = 0.01, n_threads = num.cores)

# word vectors
wv_context = glove$components
word_vectors = wv_main + t(wv_context)

# Save word vectors for later use
setwd(".../Embeddings")
save(word_vectors, file = "COHA.word_vectors.min_count_10.dim_200.RData")
load("COHA.word_vectors.min_count_10.dim_200.RData")
