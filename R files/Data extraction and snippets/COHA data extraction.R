#############################################################################
# Required libraries ------------------------------------------------------ #
#############################################################################

library(stringi)
library(stringr)
library(stopwords)
library(foreach)
library(plyr)
library(tidyr)

# The COHA data needs to be purchased from https://www.english-corpora.org/coha/



#############################################################################
# Identify stopwords in lexicon based on PoS ------------------------------ #
#############################################################################

#Read dictionary of all words used in the corpus
lexicon = read.table("lexicon.txt", header = TRUE, colClasses = "character", sep="\t", quote = "", fill = TRUE, skipNul = TRUE)
lexicon = lexicon[-1,]
row.names(lexicon) = lexicon$wID


#Identify stopwords based on PoS tags - refer to http://ucrel.lancs.ac.uk/claws7tags.html
#Everything except the following parts of speech (and their variants) is a stopword or function word: JJ, RR, VV
#plus NN except for NNA, NNB
pos.tags = noquote(sort(unique(lexicon$PoS)))
allowed.pos = c("jj", "nn", "rr", "vv")
allowed.pos.except = c("nna", "nnb")
pos.allowed.tags = as.logical(sapply(1:length(pos.tags), function(i) max(str_detect(pos.tags[i], allowed.pos)))
                              - sapply(1:length(pos.tags), function(i) max(str_detect(pos.tags[i], allowed.pos.except))))
pos.SW.tags = data.frame(posID = 1:length(pos.tags), pos = pos.tags, SW = !pos.allowed.tags)

#Add PoS SW tags to lexicon
lexicon$posID = match(lexicon$PoS, pos.tags)
lexicon$SW.PoS = pos.SW.tags[lexicon$posID, 3]

#Tag any words containing numbers as stopwords
lexicon$SW.num = str_detect(lexicon$word, "[:digit:]")

#Make an overall stopword column
lexicon$SW = as.logical(pmax(lexicon$SW.PoS, lexicon$SW.num))

#Remove unwanted variables/data
rm(pos.tags, allowed.pos, allowed.pos.except, pos.allowed.tags, pos.SW.tags); gc()



#############################################################################
# Clean data -------------------------------------------------------------- #
#############################################################################

#Fill in missing lemmas with the relevant words instead and remove any whitespace
lexicon$lemma.noblanks = ifelse(lexicon$lemma == "", lexicon$word, lexicon$lemma)
lexicon$lemma.noblanks = ifelse(lexicon$lemma.noblanks == "", lexicon$wordCS, lexicon$lemma.noblanks)
lexicon$lemma.noblanks = stri_enc_toutf8(lexicon$lemma.noblanks, TRUE) #to avoid an error in the next line
lexicon$lemma.noblanks = str_trim(lexicon$lemma.noblanks)


#Identify any "words" that are just punctuation
lexicon$punct = str_detect(lexicon$wordCS, "[:alnum:]", negate = TRUE)
non.punct.idx = which(!lexicon$punct) #index of rows that are not /just/ punctuation
lexicon$SW = as.logical(pmax(lexicon$SW, lexicon$punct)) #tag these as stopwords


#Remove any punctuation/symbols such as hyphens from the start of words
punct.idx = which(str_detect(lexicon$lemma.noblanks[non.punct.idx], "^[[:alnum:]]", negate = TRUE))
lemmas.with.punct = lexicon$lemma.noblanks[non.punct.idx][punct.idx]
punct.idx.sub = 1:length(lemmas.with.punct)

counter = length(punct.idx.sub)
while (counter > 0) {
  lemmas.with.punct[punct.idx.sub] = str_sub(lemmas.with.punct[punct.idx.sub], 2, -1)
  punct.idx.sub = which(str_detect(lemmas.with.punct, "^[[:alnum:]]", negate = TRUE))
  counter = length(punct.idx.sub)
  print(paste(counter, "remaining"))
}
lexicon$lemma.noblanks[non.punct.idx][punct.idx] = lemmas.with.punct


#Remove any punctuation/symbols such as hyphens from the end of words
punct.idx = which(str_detect(lexicon$lemma.noblanks[non.punct.idx], "[[:alnum:]]$", negate = TRUE))
lemmas.with.punct = lexicon$lemma.noblanks[non.punct.idx][punct.idx]
punct.idx.sub = 1:length(lemmas.with.punct)

counter = length(punct.idx.sub)
while (counter > 0) {
  lemmas.with.punct[punct.idx.sub] = str_sub(lemmas.with.punct[punct.idx.sub], 1, -2)
  punct.idx.sub = which(str_detect(lemmas.with.punct, "[[:alnum:]]$", negate = TRUE))
  counter = length(punct.idx.sub)
  print(paste(counter, "remaining"))
}
lexicon$lemma.noblanks[non.punct.idx][punct.idx] = lemmas.with.punct


#Tag any blanks (including any new blanks created from above procedure) as stopwords
lexicon$SW.blank = (lexicon$lemma.noblanks == "")
lexicon$SW = as.logical(pmax(lexicon$SW, lexicon$SW.blank))
lexicon = lexicon[,c(1:8,12,9:11)]


#Remove unwanted variables/data
rm(non.punct.idx, punct.idx, punct.idx.sub, lemmas.with.punct, counter); gc()



#############################################################################
# Identify common English stopwords --------------------------------------- #
#############################################################################

#sources containing English language in the 'stopwords' package
eng.sources = c(1,2,4,5,7) 

#collect all stopwords
stopwords.all = foreach(l=eng.sources, .combine = c) %do% {
  stopwords(language = "en", source = stopwords_getsources()[l])
}
stopwords.all = sort(unique(stopwords.all))

#Add stopword tags to lexicon
lexicon$SW.eng = (lexicon$lemma.noblanks %in% stopwords.all)
lexicon$SW = as.logical(pmax(lexicon$SW, lexicon$SW.eng))
lexicon = lexicon[,c(1:9,13,10:12)]


#Remove unwanted variables/data
rm(l, eng.sources, stopwords.all); gc()



#############################################################################
# Assign lemma IDs -------------------------------------------------------- #
#############################################################################

non.SW.idx = which(!lexicon$SW) #index of rows that are not stopwords
lemmas.unique = sort(unique(lexicon$lemma.noblanks[non.SW.idx]))
lemmas.unique = data.frame(lemmaID = 1:length(lemmas.unique), lemma = lemmas.unique, stringsAsFactors = FALSE)

lexicon$lemmaID = ""
lexicon$lemmaID[non.SW.idx] = match(lexicon$lemma.noblanks[non.SW.idx], lemmas.unique$lemma)
lexicon$lemma.noSWs = ifelse(lexicon$SW,"",lexicon$lemma.noblanks)


#remove unnecessary columns
lexicon = lexicon[,-(7:10)]
lexicon = lexicon[,c(1:8,10:11,9)]

#Remove unwanted variables/data
rm(non.SW.idx); gc()

#Save files that will be needed later
# save(lexicon, file = "lexicon.RData")
# write.csv(lemmas.unique, "lemmas.unique.csv", row.names = FALSE)

#Remove unwanted variables/data
rm(non.SW.idx, lemmas.unique); gc()

#Load files
load("lexicon.RData")
lemmas.unique = read.csv("lemmas.unique.csv", header = TRUE, colClasses = "character")



#############################################################################
# Process database files -------------------------------------------------- #
#############################################################################

# Load each database file and perform the following actions:
# add lemmas (SWs treated as blanks), remove punctuation, save file as an R object
filenames = list.files(getwd())
n = length(filenames)

for (k in 1:n) {
  
  print(paste("File", k, "of", n, "..."))
  
  current.file = read.table(filenames[k], col.names = c("textID", "ID", "wordID"))
  current.file = cbind(current.file,
                       word = lexicon[paste(current.file$wordID), "wordCS"],
                       lemma = lexicon[paste(current.file$wordID), "lemma.noSWs"],
                       punct = lexicon[paste(current.file$wordID), "punct"])
  current.file = current.file[!current.file$punct, -c(2,3,6)]
  
  save(current.file, file = paste("processed_", filenames[k], ".RData", sep = ""))
  
  print("Done!")
  
}#for k

#Remove unwanted variables/data
rm(current.file, lexicon, filenames, k, n); gc()


# Manual step: move processed files to a new location and change directory to this location
# Load all processed files and stack them above each other into a single object
filenames = list.files(getwd())
n = length(filenames)

for (k in 1:n) {
  
  print(paste("File", k, "of", n, "..."))
  
  load(paste(filenames[k]))
  
  if (k == 1) {
    texts.combined = current.file
  } else {
    texts.combined = rbind(texts.combined, current.file, make.row.names = FALSE)
  }
  
  print("Done!")
  
}#for k

#Remove unwanted variables/data
rm(current.file, filenames, k, n); gc()


#Save data to file
# save(texts.combined, file = "texts.combined.RData")

#Load data
load("texts.combined.RData")



#############################################################################
# Generate snippets ------------------------------------------------------- #
#############################################################################

#Function to get the indices for a snippet, given the idx of the target word
get.indices = function(idx, snippet.radius) {
  (idx - snippet.radius):(idx + snippet.radius)
}#get.indices


#Function to replace the indices with NA if outside text
filter.indices = function(idx, start, end) {
  idx = ifelse(idx<start,NA,idx) 
  idx = ifelse(idx>end,  NA,idx) 
  return(idx)
}#filter.indices

filter.indices.all = function(indices, start, end) {
  apply(indices, 2, filter.indices, start=start, end=end)
}#filter.indices.all


#Function to generate snippets for a given lemma and snippet size, including the target word itself
generate.snippets = function(target.lemma, snippet.size, words = FALSE) {
  
  #snippet.size must be an odd integer - need to add a check for this
  snippet.radius = (snippet.size - 1) / 2
  
  #indices of positions where target word occurs
  word.instances = which(texts.combined$lemma %in% target.lemma)
  
  #list of texts containing the target word
  textIDs = texts.combined$textID[word.instances]
  textIDs.unique = unique(textIDs)
  textIDs.index = match(textIDs, textIDs.unique)
  
  #word counts, start and end indices for each unique text
  text.word.counts = table(texts.combined$textID)[paste(textIDs.unique)]
  text.indices.start = match(textIDs.unique, texts.combined$textID)
  text.indices.end = text.indices.start + as.numeric(text.word.counts) - 1
  
  #date and genre of each unique text
  doc.data.idx = match(textIDs.unique, document.data$textID)
  date = document.data$year[doc.data.idx]
  genre = document.data$genre[doc.data.idx]
  
  #combind all unique text info into a single object
  text.info = cbind(textIDs.unique, text.indices.start, text.indices.end, date, genre)
  
  #row indices in texts.combined for the lemmas used in the snippets
  lemma.indices = t(sapply(word.instances, get.indices, snippet.radius=snippet.radius))
  
  #append text info to above
  lemma.indices.with.text.info = cbind(lemma.indices, text.info[textIDs.index,])
  
  #replace with NA if lemma index outside text
  lemma.indices.filtered = filter.indices.all(lemma.indices,
                                              as.numeric(lemma.indices.with.text.info[,"text.indices.start"]),
                                              as.numeric(lemma.indices.with.text.info[,"text.indices.end"]))
  
  #replace lemma indices with lemmas (blanks indicate stopwords)
  snippets.lemmas = cbind(matrix(texts.combined$lemma[lemma.indices.filtered], ncol = snippet.size),
                          textID = lemma.indices.with.text.info[,"textIDs.unique"],
                          date = lemma.indices.with.text.info[,"date"],
                          genre = lemma.indices.with.text.info[,"genre"])
  
  colnames(snippets.lemmas)[1:snippet.size] = c(paste("V", 1:snippet.radius, sep = ""),
                                                "target.word", paste("V", (snippet.radius+1):(snippet.size-1), sep = ""))
  
  snippets.lemmas = as.data.frame(snippets.lemmas, stringsAsFactors = FALSE)
  
  #if (words == TRUE) then generate the snippets with actual words rather than lemmas 
  if (words) {
    snippets.words = cbind(matrix(texts.combined$word[lemma.indices.filtered], ncol = snippet.size),
                           snippets.lemmas[,-(1:snippet.size)])
    
    colnames(snippets.words)[1:snippet.size] = c(paste("V", 1:snippet.radius, sep = ""),
                                                 "target.word", paste("V", (snippet.radius+1):(snippet.size-1), sep = ""))
    
    snippets.words = as.data.frame(snippets.words, stringsAsFactors = FALSE)
    
    out = list(lemmas = snippets.lemmas, words = snippets.words)
    
  } else {
    
    out = snippets.lemmas
    
  }#else
  
  return(out)
  
}#generate.snippets



#############################################################################
# Quick start notes ------------------------------------------------------- #
#############################################################################

# In order to generate snippets, need to load three items:

# 1. 
load("texts.combined.RData")

# 2. 
library(openxlsx)
document.data = read.xlsx("sources_coha.xlsx")
document.data = data.frame(textID = document.data$textID, year = document.data$year, genre = document.data$genre, stringsAsFactors = FALSE)

# 3. 
# Load functions defined in the previous section 


# Set snippet length and generate snippets
snippet.length = 14 #not counting target word
target.lemma = "bank"
snippets = generate.snippets(target.lemma = target.lemma, snippet.size = snippet.length + 1, TRUE)
snippets.words = snippets$words
snippets = snippets$lemmas


#Save snippets
# write.csv(snippets, "snippets.bank.length14.csv", row.names = FALSE)
# write.csv(snippets.words, "snippets.bank.length14.words.csv", row.names = FALSE)

#Load snippets
snippets = read.csv("snippets.bank.length14.csv", header = TRUE, colClasses = "character")
snippets.words = read.csv("snippets.bank.length14.words.csv", header = TRUE, colClasses = "character")



#############################################################################
# Subsample snippets (optional) ------------------------------------------- #
#############################################################################

# Discard "empty" snippets, i.e. snippets consisting only of stopwords
snippet.length  = ncol(snippets) - 4 #Not counting target word
empty.snippets = which(rowSums(snippets[,1:(snippet.length+1)] == "") == snippet.length)
snippets = snippets[-empty.snippets,]
snippets.words = snippets.words[-empty.snippets,]


# If too many snippets, we can subsample a max number of snippets for each time period and genre 
# (prioritising snippets with fewer stopwords)

snippet.lengths = snippet.length - as.numeric(rowSums(snippets[,1:(snippet.length+1)] == "", na.rm = TRUE))

#Convert date to numeric variable
snippets$date = as.numeric(snippets$date)

#Label time periods
break.points = seq(1810, 2010, 20)
num.break.points = length(break.points)
period.labels = paste(break.points[1:(num.break.points-1)], break.points[2:num.break.points], sep = "-")
snippets$Period = cut(snippets$date, breaks = break.points, labels = period.labels, right = FALSE, include.lowest = TRUE)
#table(snippets$Period)
snippets$Time = cut(snippets$date, breaks = break.points, labels = FALSE, right = FALSE, include.lowest = TRUE)
#table(snippets$Time)

#Get counts of snippets in each time period and genre
snippet.counts = table(snippets$genre, snippets$Time, dnn = c("Genre", "Time"))


#Sample snippets
set.seed(100)
max.snippets.per.block = 100

sample.indices = NULL
for (g in seq_along(dimnames(snippet.counts)$Genre)) {
  for (t in as.numeric(dimnames(snippet.counts)$Time)) {
    snippet.indices = which((snippets$genre == dimnames(snippet.counts)$Genre[g]) & (snippets$Time == t))
    if (length(snippet.indices) == 0) next
    sample.indices = c(sample.indices, sample(snippet.indices, min(max.snippets.per.block, snippet.counts[g,t]), 
                                              prob = snippet.lengths[snippet.indices]^2))
  }#for t
}#for g
sample.indices = sort(sample.indices)

snippets.all = snippets
snippets.not.sampled = snippets[-sample.indices,]
snippets = snippets[sample.indices,]
snippets.words = snippets.words[sample.indices,]



#############################################################################
# Identify informative/uninformative words -------------------------------- #
#############################################################################

# Identify informative/uninformative words based on their frequency in the not sampled data
# Mark words representing the top M proportion as Imp
M = 0.7

#Calculate word frequencies
snippet.length  = ncol(snippets) - 6 #Not counting target word
empirical.word.freqs = sort(table(factor(stack(snippets.not.sampled[,1:(snippet.length+1)])$values, levels = unique(stack(snippets.not.sampled[,1:(snippet.length+1)])$values))), decreasing = TRUE)
empirical.word.freqs = empirical.word.freqs[-which(names(empirical.word.freqs) == "")] #exclude stopwords (blanks)

#Identify uninformative words
m = Position(function(x) x >= M, cumsum(empirical.word.freqs)/sum(empirical.word.freqs))
words.sampled = unique(stack(snippets[,1:(snippet.length+1)])$values)
Imp = c("", names(empirical.word.freqs[1:m]))
Unimp = words.sampled[which(!(words.sampled %in% Imp))]

#Replace with "Unimp"
for (i in 1:(snippet.length+1)) {
  snippets[snippets[,i] %in% Unimp, i] = "Unimp"
}

# Discard "empty" snippets, i.e. snippets consisting only of Unimp
empty.snippets = which(rowSums(snippets[,1:(snippet.length+1)] == "" | snippets[,1:(snippet.length+1)] == "Unimp") == snippet.length)
snippets = snippets[-empty.snippets,]
snippets.words = snippets.words[-empty.snippets,]



#############################################################################
# Prepare file for hand annotation ---------------------------------------- #
#############################################################################

snippets.words = unite(snippets.words, sentence, 1:(snippet.length+1), sep = " ", remove = FALSE)
snippets.words = snippets.words[,c(2:ncol(snippets.words),1)]

snippets.words$sense.id = ""

library(openxlsx)
write.xlsx(snippets.words, "snippets.bank.length14.words.annotated.xlsx")

#Hand annotate in Excel and read into R again
snippets.words = read.xlsx("snippets.bank.length14.words.annotated.xlsx")
snippets.words$sense.id = as.factor(snippets.words$sense.id)

sense.ids = levels(snippets.words$sense.id)
sense.ids = sense.ids[c(6,2,7,5,3,4,1,8)] #Re-order so that senses 1 and 2 are river-bank and institution bank
snippets.words$sense.id = factor(snippets.words$sense.id, levels = sense.ids)



#############################################################################
# Re-arrange snippet data ------------------------------------------------- #
#############################################################################

#remove target word itself
snippets = snippets[,-which(colnames(snippets) == "target.word")]

#Convert date to numeric variable
snippets$date = as.numeric(snippets$date)

#Label time periods
break.points = seq(1810, 2010, 20)
num.break.points = length(break.points)
period.labels = paste(break.points[1:(num.break.points-1)], break.points[2:num.break.points], sep = "-")
snippets$Period = cut(snippets$date, breaks = break.points, labels = period.labels, right = FALSE, include.lowest = TRUE)
#table(snippets$Period)
snippets$Time = cut(snippets$date, breaks = break.points, labels = FALSE, right = FALSE, include.lowest = TRUE)
#table(snippets$Time)

#Create unique ID for each snippet
snippets$SnippetID = 1:nrow(snippets)
snippets.words$SnippetID = snippets$SnippetID
rownames(snippets) = snippets$SnippetID
rownames(snippets.words) = snippets.words$SnippetID

#We only need the words, time period and ID for each snippet for the analysis
#Other data can be recovered using the snippet ID
#Therefore, split the data into two parts: snippets and snippets.info
nc = ncol(snippets)
snippet.length = nc-6
snippets.info = snippets[,c(nc,(nc-5):(nc-1))]
snippets = snippets[,c(1:snippet.length, (nc-1):nc)]
snippets.info$sense.id = snippets.words$sense.id


#Extract fixed parameters from snippet data
num.snippets    = nrow(snippets)
snippet.length  = ncol(snippets) - 2 #Not counting target word
snippet.lengths = snippet.length - as.numeric(rowSums(snippets[,1:snippet.length] == "Unimp", na.rm = TRUE) + 
                                                rowSums(snippets[,1:snippet.length] == "", na.rm = TRUE)) #Not counting Unimp or SWs (blanks)
num.periods     = max(snippets$Time) #Number of time periods
words.used      = sort(unique(stack(snippets[,1:snippet.length])$values)) #Vector of all unique words used across all snippets
words.used      = words.used[-which(words.used %in% c("Unimp",""))] #Ignore Unimp and SWs (blanks)
num.words       = length(words.used) #Size of vocabulary
#Imp[which(!(Imp %in% words.used))] #Imp words not used in the sampled snippet vocabulary. Can ignore if small number; 
				     #otherwise the sample is probably not very representative of the actual snippet population
#snippet.counts = table(factor(snippets$Time, levels = 1:num.periods)) #Count of snippets for each time period

#Replace each word in each snippet with its position in the vector words.used
#This will be its unique ID and will correspond to a row in the psi matrix
#This will save time when looking up an entry for a word in the psi matrix
#If we want the actual word, it can be recovered by looking up the word ID in vector words.used
snippets.word.index = apply(snippets[,1:snippet.length], 2, match, table = words.used)
snippets[,1:snippet.length] = snippets.word.index

#Add genre information
num.genres = 1
snippets$genre = 1

# genres = c("FIC", "Non-FIC")
# num.genres = 2
# snippets$genre = (snippets.info$genre != "FIC") + 1 #bank: FIC = 1, Non-FIC = 2

# genres = c("NEWS/MAG", "FIC/NF")
# num.genres = 2
# snippets$genre = (snippets.info$genre == "NEWS" | snippets.info$genre == "MAG") + 1 #bank: FIC/NF = 1, NEWS/MAG = 2

# genres = levels(factor(snippets.info$genre))
# num.genres = length(genres)
# snippets$genre = as.numeric(factor(snippets.info$genre))



#Save snippets
# save(snippets, snippets.info, words.used, num.words, num.snippets, snippet.length, snippet.lengths, num.periods, num.genres, file = "bank_snippets.RData")
# save(snippets, snippets.words, snippets.info, words.used, num.words, num.snippets, snippet.length, snippet.lengths, num.periods, num.genres, file = "bank_snippets.words.RData")

#Load snippets
load("bank_snippets.RData")
load("bank_snippets.words.RData")



#############################################################################
# Discard variables/data not needed --------------------------------------- #
#############################################################################

rm(texts.combined, document.data)
rm(g, t, empty.snippets, snippet.counts, max.snippets.per.block, sample.indices, snippet.indices, snippets.all, snippets.not.sampled)
rm(M, m, i, words.sampled, empirical.word.freqs, Imp, Unimp)
rm(nc, break.points, num.break.points, period.labels, snippets.word.index, sense.ids, idx, genres)
rm(generate.snippets, get.indices, filter.indices, filter.indices.all)
gc()
