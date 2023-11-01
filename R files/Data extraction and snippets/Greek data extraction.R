#############################################################################
# Load packages ----------------------------------------------------------- #
#############################################################################

library(XML)
library(methods)
library(plyr)



#############################################################################
# Data extraction - document level ---------------------------------------- #
#############################################################################

# Obtain files from https://figshare.com/articles/dataset/The_Diorisis_Ancient_Greek_Corpus/6187256
# extract and set working directory to this location
setwd(".../Diorisis")
filenames = list.files(getwd())


# Create vectors to store document level metadata

n = length(filenames)

AuthorID = character(n)
Author   = character(n)
ID       = character(n)
Title    = character(n)
date     = character(n)
genre    = character(n)
subgenre = character(n)

# Cycle through each document and extract metadata

progress.bar = txtProgressBar(min = 1, max = n) #Show progress

for (k in 1:n) {
  
  xmlobject = xmlParse(filenames[k])
  
  AuthorID[k] = xpathSApply(xmlobject, "/TEI.2/teiHeader/fileDesc/titleStmt/tlgAuthor", xmlValue, recursive = FALSE)
  Author[k]   = xpathSApply(xmlobject, "/TEI.2/teiHeader/fileDesc/titleStmt/author", xmlValue, recursive = FALSE)
  ID[k]       = xpathSApply(xmlobject, "/TEI.2/teiHeader/fileDesc/titleStmt/tlgId", xmlValue, recursive = FALSE)
  Title[k]    = xpathSApply(xmlobject, "/TEI.2/teiHeader/fileDesc/titleStmt/title", xmlValue, recursive = FALSE)
  date[k]     = xpathSApply(xmlobject, "/TEI.2/teiHeader/profileDesc/creation/date", xmlValue, recursive = FALSE)
  genre[k]    = xpathSApply(xmlobject, "/TEI.2/teiHeader/xenoData/genre", xmlValue, recursive = FALSE)
  subgenre[k] = xpathSApply(xmlobject, "/TEI.2/teiHeader/xenoData/subgenre", xmlValue, recursive = FALSE)
  
  setTxtProgressBar(progress.bar, k) #Show progress
  
}

close(progress.bar)

# filenames[165]
# filenames[222]
# File "Demosthenes (0014) - kata meidiou peri tou kondulou (021).xml" doesn't work on Windows PC due to Greek letters in the filename. Filename has to be changed to read in the data.


# Create unique ID for each document
DocID = paste("Doc", AuthorID, ID, sep = ".") 

# Doc.0086.029 and Doc.0060.001 are not unique, so manually rename these
DocID[DocID == "Doc.0086.029"] = c("Doc.0086.029a", "Doc.0086.029b")
DocID[DocID == "Doc.0060.001"] = c("Doc.0060.001a", "Doc.0060.001b", "Doc.0060.001c")

# Save document level data in a data frame
document.data = data.frame(cbind(DocID, Filename = filenames, AuthorID, Author, ID, Title, date, genre, subgenre))


# Save data to file
# write.csv(document.data, file = "Document_data.csv", row.names = FALSE)


# Read data
document.data = read.csv("Document_data.csv", header = TRUE, colClasses = "character")



#############################################################################
# Data extraction - word level -------------------------------------------- #
#############################################################################

# Create a list of empty (named) objects to store data for each document
n = nrow(document.data)
words = vector("list", n)
names(words) = document.data$DocID

# Cycle through each file, extract words and assign to an element in the list
for (k in 1:n) {
  
  xmlobject = xmlParse(document.data$Filename[k])
  
  node.set.word = getNodeSet(xmlobject, "/TEI.2/text/body/sentence/word") #Set of all "word" nodes
  
  num.words = xmlSize(node.set.word) #number of words in document
  
  #Set up data frame to store words, lemma IDs and part of speech tags, plus sentence id and location
  df = data.frame(form    = character(num.words), 
                  lemmaID = character(num.words),
                  lemma   = character(num.words), 
                  POS     = character(num.words),
                  sentence.id       = character(num.words),
                  sentence.location = character(num.words),
                  stringsAsFactors = FALSE
  )
  
  #Cycle through each word in the document and extract word-level data
  
  progress.bar = txtProgressBar(min = 1, max = num.words) #Show progress
  
  print(paste("Document", k, "of", n))
  
  for (w in 1:num.words) {
    
    df$form[w]    = xmlAttrs(node.set.word[[w]])["form"] #word form
    df$lemmaID[w] = xmlAttrs(node.set.word[[w]][["lemma"]])["id"]    #lemma id
    df$lemma[w]   = xmlAttrs(node.set.word[[w]][["lemma"]])["entry"] #lemma entry
    df$POS[w]     = xmlAttrs(node.set.word[[w]][["lemma"]])["POS"]   #part of speech
    
    df$sentence.id[w]       = xmlAttrs(xmlParent(node.set.word[[w]]))["id"] #sentence id
    df$sentence.location[w] = xmlAttrs(xmlParent(node.set.word[[w]]))["location"] #sentence location
    
    setTxtProgressBar(progress.bar, w) #Show progress
    
  }
  
  close(progress.bar)
  
  write.csv(df, paste(document.data$DocID[k], "csv", sep = "."), row.names = FALSE) #Write data to file
  
  words[[k]] = df #Save data frame in list
  
}

# Save list to file
# save(words, file = "words.RData")


# Read data
load("words.RData")



#############################################################################
# Create dictionary of lemmas and their IDs ------------------------------- #
#############################################################################

lemmas = rbind.fill(words) #Stack all data frames one above the other

lemmas.unique = unique(lemmas[,-c(1,5,6)]) #remove word forms, sentence id/location and duplicates
lemmas.unique = lemmas.unique[order(lemmas.unique$lemma),] #sort alphabetically

# Save data to file
# write.csv(lemmas, "lemmas.csv", row.names = FALSE)
# write.csv(lemmas.unique, "lemmas.unique.csv", row.names = FALSE)

# Read data
lemmas.unique = read.csv("lemmas.unique.csv", header = TRUE, colClasses = "character")



#############################################################################
# Identify stopwords based on PoS ----------------------------------------- #
#############################################################################

# This section works on the assumption that all articles, particles, conjunctions,
# interjections, prepositions and pronouns are stopwords or function words

lemmas.unique = read.csv("lemmas.unique.csv", header = TRUE, colClasses = "character")
lemmas.unique$SW = with(lemmas.unique, ifelse(POS == "adjective" | POS == "noun" | 
                                                POS == "adverb" | POS == "verb" | 
                                                POS == "proper", FALSE, TRUE))



#############################################################################
# Identify further stopwords ---------------------------------------------- #
#############################################################################

# Identify further stopwords based on stopword lists from Ancient Greek
# stopwords.list = unlist(read.table("Stopwords.txt", colClasses = "character")) # NOT USED IN THIS VERSION

library(stopwords)
# stopwords(language = "el", source = "misc")
# stopwords(language = "el", source = "stopwords-iso")

stopwords.all = c(#stopwords.list,
                  # stopwords(language = "el", source = "misc"), # NOT USED IN THIS VERSION
                  stopwords(language = "el", source = "stopwords-iso"))

lemmas.unique$SW = as.logical(pmax(lemmas.unique$SW, lemmas.unique$lemma %in% stopwords.all))


rm(stopwords.list, stopwords.all)


# Save data to file
#write.csv(lemmas.unique, "lemmas.unique.with.SW.tags.csv", row.names = FALSE)

# Read data
lemmas.unique = read.csv("lemmas.unique.with.SW.tags.csv", header = TRUE, colClasses = "character")



#############################################################################
# Generate snippets ------------------------------------------------------- #
#############################################################################

# # Function to generate snippets for a given lemma ID and snippet size
# 
# generate.snippets = function(target.lemmaID, snippet.size) {
# 
#   #snippet.size must be an odd integer - need to add a check for this
# 
#   snippet.radius = (snippet.size - 1) / 2
# 
#   #Create list of empty (named) objects to store data for each document
#   n = length(words)
#   snippets = vector("list", n)
#   names(snippets) = names(words)
# 
#   #Cycle through each document and generate snippets of given size around each instance of target word
# 
#   progress.bar = txtProgressBar(min = 1, max = n) #Show progress
# 
#   for (k in 1:n) {
# 
#     word.instances = which(words[[k]]$lemmaID %in% target.lemmaID) #indices of positions where target word occurs
# 
#     l = length(word.instances) #number of times target word occurs
# 
#     if (l == 0) next #target word does not appear in the document in this case
# 
#     doc.length = nrow(words[[k]]) #total number of words in the document
# 
#     #Cycle through each instance of target word and generate snippets
#     for (i in 1:l) { 
#       
#       index = (word.instances[i] - snippet.radius):(word.instances[i] + snippet.radius)
#       
#       #Deal with cases where the target word occurs at the start or end of a document
#       index = ifelse(index<1,NA,index)
#       index = ifelse(index>doc.length,NA,index)
# 
#       #Save snippets as a data frame
#       snippets[[k]] = as.data.frame(rbind(snippets[[k]],
#                                           unlist(c(words[[k]]$lemmaID[index],
#                                                                   words[[k]][word.instances[i], c("sentence.id", "sentence.location")]))),
#                                     stringsAsFactors = FALSE)
#     }
# 
#     #Append document ID, date and genre information to each snippet
#     doc.info = subset(document.data, DocID == names(snippets)[k], drop = TRUE)
#     snippets[[k]] = cbind(snippets[[k]], doc.info, stringsAsFactors = FALSE)
# 
#     setTxtProgressBar(progress.bar, k) #Show progress
# 
#   }
#   close(progress.bar)
# 
#   #Stack all snippets one above the other
#   out = rbind.fill(snippets)
# 
#   return(out)
# 
# }#generate.snippets



#############################################################################
# Generate snippets (with target word itself excluded) -------------------- #
#############################################################################

# This function is the same as that above, except that the target word itself 
# is excluded from each snippet. This allows us to focus just on the words around
# the target word.

# Function to generate snippets for a given lemma ID and snippet size

generate.snippets = function(target.lemmaID, snippet.size) {
  
  #snippet.size must be an odd integer - need to add a check for this
  
  snippet.radius = (snippet.size - 1) / 2
  
  #Create list of empty (named) objects to store data for each document
  n = length(words)
  snippets = vector("list", n)
  names(snippets) = names(words)
  
  #Cycle through each document and generate snippets of given size around each instance of target word
  
  progress.bar = txtProgressBar(min = 1, max = n) #Show progress
  
  for (k in 1:n) {
    
    word.instances = which(words[[k]]$lemmaID %in% target.lemmaID) #indices of positions where target word occurs
    
    l = length(word.instances) #number of times target word occurs
    
    if (l == 0) next #target word does not appear in the document in this case
    
    doc.length = nrow(words[[k]]) #total number of words in the document
    
    #Cycle through each instance of target word and generate snippets
    for (i in 1:l) { 
      
      index = (word.instances[i] - snippet.radius):(word.instances[i] + snippet.radius)
      index = index[-(snippet.radius + 1)] #remove target word itself
      
      #Deal with cases where the target word occurs at the start or end of a document
      index = ifelse(index<1,NA,index)
      index = ifelse(index>doc.length,NA,index)
      
      #Save snippets as a data frame
      snippets[[k]] = as.data.frame(rbind(snippets[[k]], 
                                          unlist(c(words[[k]]$lemmaID[index], 
                                                   words[[k]][word.instances[i], c("sentence.id", "sentence.location")]))), 
                                    stringsAsFactors = FALSE)
    }
    
    #Append document ID, date and genre information to each snippet
    doc.info = subset(document.data, DocID == names(snippets)[k], drop = TRUE)
    snippets[[k]] = cbind(snippets[[k]], doc.info, stringsAsFactors = FALSE)
    
    setTxtProgressBar(progress.bar, k) #Show progress
    
  }
  close(progress.bar)
  
  #Stack all snippets one above the other
  out = rbind.fill(snippets) 
  
  return(out)
  
}#generate.snippets



#############################################################################
# Quick start notes ------------------------------------------------------- #
#############################################################################

# In order to generate snippets, need to load three items:

# 1.
load("words.RData")

# 2. 
document.data = read.csv("Document_data.csv", header = TRUE, colClasses = "character")

# 3. 
# Function generate.snippets() - see previous section

# Set snippet length and generate snippets
snippet.length = 14 #not counting target word
snippets = generate.snippets(target.lemmaID = "59339", snippet.size = snippet.length + 1) #kosmos "κόσμος" lemmaID 59339: 2603 instances
#snippets = generate.snippets(target.lemmaID = "15281", snippet.size = snippet.length + 1) #harmonia "ἁρμονία" lemmaID 15281: 655 instances
#snippets = generate.snippets(target.lemmaID = "69419", snippet.size = snippet.length + 1) #mus μῦς lemmaID 69419: 214 instances



#############################################################################
# Append expert-annotated data -------------------------------------------- #
#############################################################################

# Obtain files from https://figshare.com/articles/dataset/Ancient_Greek_semantic_annotation_datasets/7882940
# extract and set working directory to this location

# Read expert-annotated data
senses_data = read.table("senses_kosmos.txt", sep="\t", header = TRUE, colClasses = "character")
#senses_data = read.table("senses_harmonia.txt", sep="\t", header = TRUE, colClasses = "character")
#senses_data = read.table("senses_mus.txt", sep="\t", header = TRUE, colClasses = "character")

#Merge the two datasets based on specified common columns
snippets.annotated = merge(senses_data, snippets, sort = FALSE,
                           by.x = c("tlg.author", "tlg.work", "sentence.id", "sentence.location", "date", 
                                    "genre", "subgenre", "author", "work"), 
                           by.y = c("AuthorID", "ID", "sentence.id", "sentence.location", "date", 
                                    "genre", "subgenre", "Author", "Title"))

#Re-arrange columns and discard the datasets not required
snippets = snippets.annotated[,c(16:ncol(snippets.annotated), 1:15)]

rm(snippets.annotated, senses_data, document.data, generate.snippets); gc()



#############################################################################
# Identify stopwords and infrequent words in snippets --------------------- #
#############################################################################

lemmas.unique = read.csv("lemmas.unique.with.SW.tags.csv", header = TRUE, colClasses = "character")

# Identify stopwords replace these with value "SW"
stopword.lemmaIDs = lemmas.unique$lemmaID[lemmas.unique$SW == TRUE]
for (i in 1:snippet.length) {
  snippets[snippets[,i] %in% stopword.lemmaIDs, i] = "SW"
}


# Identify words that appear less than the specified number of times in the entire corpus and replace these with value "infreq"

#Calculate word frequencies
words = rbind.fill(words)
empirical.word.freqs = table(factor(words$lemmaID, levels = unique(stack(snippets[,1:snippet.length])$values)))
empirical.word.freqs = sort(empirical.word.freqs, decreasing = TRUE)
empirical.word.freqs = empirical.word.freqs[-length(empirical.word.freqs)] #discard SW

#Get a plot of vocabulary size vs the cutoff min frequency
cutoff.freq = 0:50
vocab.size = sapply(cutoff.freq, function(m){sum(empirical.word.freqs>m)})
plot(cutoff.freq, vocab.size, type = "l", lwd = 2, xlab = "cutoff frequency", ylab = "vocabulary size")

#specify min count
min.count = 10
abline(v = min.count, h = vocab.size[min.count+1], col = "red", lty = 2)

#Identify infrequent words
infreqs = names(which(empirical.word.freqs < min.count))

#Replace with "infreq"
for (i in 1:snippet.length) {
  snippets[snippets[,i] %in% infreqs, i] = "infreq"
}

# Write data to file
write.csv(snippets, "snippets_kosmos.csv", row.names = FALSE)
#write.csv(snippets, "snippets_harmonia.csv", row.names = FALSE)
#write.csv(snippets, "snippets_mus.csv", row.names = FALSE)

rm(i, words, stopword.lemmaIDs, empirical.word.freqs, cutoff.freq, vocab.size, min.count, infreqs); gc()



#############################################################################
# Re-arrange snippet data ------------------------------------------------- #
#############################################################################

# START HERE IF SNIPPETS FILE IS ALREADY SAVED #

#Load snippets (excluding target word)
snippets = read.csv("snippets_kosmos.csv", header = TRUE, colClasses = "character")
# snippets = read.csv("snippets_harmonia.csv", header = TRUE, colClasses = "character")
# snippets = read.csv("snippets_mus.csv", header = TRUE, colClasses = "character")

#Convert date to numeric variable
snippets$date = as.numeric(snippets$date)

#Label centuries
snippets$Century = cut(snippets$date, breaks = seq(-800, 400, 100), labels = c(paste(8:1, "BC"), paste(1:4, "AD")), right = FALSE, include.lowest = TRUE)
#table(snippets$Century)

#Label time periods (NOTE: may need to adjust the endpoints manually)
min.date = floor(min(snippets$date) / 100) * 100
max.date = ceiling(max(snippets$date) / 100) * 100
snippets$Time = cut(snippets$date, breaks = seq(min.date, max.date, 100), labels = FALSE, right = FALSE, include.lowest = TRUE)
#table(snippets$Time)

#Create unique ID for each snippet
snippets$SnippetID = 1:nrow(snippets)
rownames(snippets) = snippets$SnippetID

#We only need the words, time period and ID for each snippet for the analysis
#Other data can be recovered using the snippet ID
#Therefore, split the data into two parts: snippets and snippets.info
nc = ncol(snippets)
snippet.length = nc-20
snippets.info = snippets[,c(nc,(nc-19):(nc-1))]
snippets = snippets[,c(1:snippet.length, (nc-1):nc)]

#Extract fixed parameters from snippet data
num.snippets    = nrow(snippets)
num.periods     = max(snippets$Time) #Number of time periods
words.used      = unique(stack(snippets[,1:snippet.length])$values) #Vector of all unique words used across all snippets
words.used      = words.used[-which(words.used %in% c("unknown", "SW", "infreq"))] #Ignore infrequent words, stopwords and unknowns
num.words       = length(words.used) #Size of vocabulary
snippet.length  = ncol(snippets) - 2 #Not counting target word
snippet.lengths = apply(snippets[,1:snippet.length], 1, function(i) {sum(!(i %in% c("unknown", "SW", "infreq")))}) #Not counting infrequent words, stopwords and unknowns
#snippet.counts = table(factor(snippets$Time, levels = 1:num.periods)) #Count of snippets for each time period

#Replace each word in each snippet with its position in the vector words.used
#This will be its unique ID and will correspond to a row in the psi matrix
#This will save time when looking up an entry for a word in the psi matrix
#If we want the actual word, it can be recovered by looking up the word ID in vector words.used
snippets.word.index = apply(snippets[,1:snippet.length], 2, match, table = words.used)
snippets[,1:snippet.length] = snippets.word.index

#Add genre information
num.genres = 2
snippets$genre = (snippets.info$genre != "Narrative") + 1 #kosmos: Narrative = 1, Non-narrative = 2
#snippets$genre = (snippets.info$genre != "Technical") + 1 #harmonia and mus: Technical = 1, Non-technical = 2


#Save snippets
save(snippets, snippets.info, words.used, num.words, num.snippets, snippet.length, snippet.lengths, num.periods, num.genres, file = "kosmos_snippets.RData")
# save(snippets, snippets.info, words.used, num.words, num.snippets, snippet.length, snippet.lengths, num.periods, num.genres, file = "harmonia_snippets.RData")
# save(snippets, snippets.info, words.used, num.words, num.snippets, snippet.length, snippet.lengths, num.periods, num.genres, file = "mus_snippets.RData")

#Load snippets
load("kosmos_snippets.RData")
# load("harmonia_snippets.RData")
# load("mus_snippets.RData")

#Load lemmas - may be needed when analysing results
lemmas.unique = read.csv("lemmas.unique.csv", header = TRUE, colClasses = "character")


rm(nc, max.date, min.date, snippets.word.index); gc()
