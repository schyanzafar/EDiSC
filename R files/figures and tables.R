#############################################################################
# Required libraries ------------------------------------------------------ #
#############################################################################

library(tidyr)
library(plyr)
library(MCMCpack)

options(digits = 4)



#############################################################################
# Top words --------------------------------------------------------------- #
#############################################################################

m = 10 

psi.tilde.post.mean = sapply(1:num.periods, function(t) {sapply(1:num.senses, function(k) {apply(psi.tilde.sim[(burn.in+1):num.samples,,k,t], 2, mean)})}, simplify = "array" )
dimnames(psi.tilde.post.mean) = list(Word = 1:num.words, Sense = 1:num.senses, Time = 1:num.periods)

# lemmas.unique = read.csv("lemmas.unique.csv", header = TRUE, colClasses = "character")

# # top m words with highest posterior mean for each sense separately for each time period
# top.m.words = sapply(1:num.periods, function(t) {sapply(1:num.senses, function(k) {order(psi.tilde.post.mean[,k,t], decreasing = TRUE)[1:m]})}, simplify = "array")
# dimnames(top.m.words) = list(Rank = 1:m, Sense = 1:num.senses, Time = 1:num.periods)
# 
# top.m.LemmaIDs = aaply(top.m.words, 1:3, function(x) {words.used[x]})
# #noquote(aperm(top.m.LemmaIDs,c(2,1,3))) #Gives the top words for "bank"
# top.m.original = aaply(top.m.LemmaIDs, 1:3, function(x) {lemmas.unique[match(x, lemmas.unique[,1]),2]})
# noquote(aperm(top.m.original,c(2,1,3)))

# top m words with highest posterior mean for each sense across all time periods
psi.tilde.post.mean.marg.over.time = apply(psi.tilde.post.mean, 1:2, sum) / num.periods
top.m.words.marg.over.time = sapply(1:num.senses, function(k) {order(psi.tilde.post.mean.marg.over.time[,k], decreasing = TRUE)[1:m]})
dimnames(top.m.words.marg.over.time) = list(Rank = 1:m, Sense = 1:num.senses)

top.m.marg.over.time.LemmaIDs = aaply(top.m.words.marg.over.time, 1:2, function(x) {words.used[x]})
#noquote(t(top.m.marg.over.time.LemmaIDs)) #Gives the top words for "bank"
top.m.marg.over.time.original = aaply(top.m.marg.over.time.LemmaIDs, 1:2, function(x) {lemmas.unique[match(x, lemmas.unique[,1]),2]})
noquote(t(top.m.marg.over.time.original))


# Top words - expert annotation
#idx = which(snippet.lengths > 0) #Index of non-empty snippets
idx = which(snippet.lengths > 0 & snippets.info$disambiguation == 1) #Index of non-empty snippets of type collocates
z = as.numeric(factor(snippets.info$sense.id))
snippets.expanded = gather(cbind(snippets[idx,], sense = z[idx]), key = position, value = word, 1:snippet.length)
word.counts = table(factor(snippets.expanded$word, levels = 1:num.words),
                    factor(snippets.expanded$sense, levels = 1:num.senses),
                    dnn = c("Word", "Sense"))

psi.tilde.marg.over.time = sapply(1:num.senses, function(k) {(word.counts[,k])/(sum(word.counts[,k]))}, simplify = "aray")
dimnames(psi.tilde.marg.over.time) = list(Word = 1:num.words, Sense = 1:num.senses)

#Top m words with highest empirical psi.tilde for each sense across all time periods
m = 10 
top.m.words.marg.over.time = sapply(1:num.senses, function(k) {order(psi.tilde.marg.over.time[,k], decreasing = TRUE)[1:m]})
dimnames(top.m.words.marg.over.time) = list(Rank = 1:m, Sense = 1:num.senses)

top.m.marg.over.time.LemmaIDs = aaply(top.m.words.marg.over.time, 1:2, function(x) {words.used[x]})
top.m.marg.over.time.original = aaply(top.m.marg.over.time.LemmaIDs, 1:2, function(x) {lemmas.unique[match(x, lemmas.unique[,1]),2]})
noquote(t(top.m.marg.over.time.original))



#############################################################################
# Brier scores - when model senses equal true senses ---------------------- #
#############################################################################

desired.labelling = 1:num.senses

#Posterior normalised sense probabilities
sense.probs.sim.normalised = sapply(1:num.snippets, function(d){
  sense.probs.sim[,,d]/rowSums(sense.probs.sim[,,d])}, simplify = "array")
sense.probs.sim.normalised[] = sense.probs.sim.normalised[,desired.labelling,]

#Posterior mean sense probabilities
sense.probs.post.mean = apply(sense.probs.sim.normalised[(burn.in+1):num.samples,,], 2:3, mean)

#Index of non-empty snippets
# idx = which(snippet.lengths > 0)
idx = which(snippet.lengths > 0 & snippets.info$disambiguation == 1) #of type collocates

# #Index of snippets with either the riverbank or the institution bank true sense -- run the following 2 lines for bank only
# sense.ids = levels(snippets.info$sense.id)
# idx = which(snippets.info$sense.id %in% sense.ids[1:2])
# names(idx) = idx #required for consistency with the Greek data

#Brier score
(BS = measures::multiclass.Brier(t(sense.probs.post.mean[,idx]), as.numeric(factor(snippets.info$sense.id[idx]))))



#############################################################################
# Brier scores - Greek target words with more than 3 senses --------------- #
#############################################################################

#Define sense order corresponding to
#kosmos -- decoration, order, world, noise
#mus -- mouse, muscle, mussel, noise
#harmonia -- abstract, concrete, musical, noise
desired.labelling = list(c(1), c(2), c(3,4), c())

#Normalise posterior sense probabilities over all senses (incl noise)
sense.probs.sim.normalised = sapply(1:num.snippets, function(d){
  sense.probs.sim[,,d]/rowSums(sense.probs.sim[,,d])}, simplify = "array")

#Combine senses into the 3 true senses
sense.probs.sim.normalised.3 = sense.probs.sim.normalised[,1:3,]
sense.probs.sim.normalised.3[,1,] = apply(sense.probs.sim.normalised[,desired.labelling[[1]],,drop=FALSE], 3, rowSums)
sense.probs.sim.normalised.3[,2,] = apply(sense.probs.sim.normalised[,desired.labelling[[2]],,drop=FALSE], 3, rowSums)
sense.probs.sim.normalised.3[,3,] = apply(sense.probs.sim.normalised[,desired.labelling[[3]],,drop=FALSE], 3, rowSums)
sense.probs.sim.normalised = sapply(1:num.snippets, function(d){
  sense.probs.sim.normalised.3[,,d] / rowSums(sense.probs.sim.normalised.3[,,d])}, simplify = "array")
rm(sense.probs.sim.normalised.3); gc()

#Posterior mean sense probabilities
sense.probs.post.mean = apply(sense.probs.sim.normalised[(burn.in+1):num.samples,,], 2:3, mean)

#Index of non-empty snippets
# idx = which(snippet.lengths > 0)
idx = which(snippet.lengths > 0 & snippets.info$disambiguation == 1) #of type collocates

#Brier score
(BS = measures::multiclass.Brier(t(sense.probs.post.mean[,idx]), as.numeric(factor(snippets.info$sense.id[idx]))))



#############################################################################
# Brier scores - "bank" with more than 2 senses --------------------------- #
#############################################################################

#Define sense order corresponding to: riverbank, institution, noise
desired.labelling = list(c(1), c(2,3), c())

#Normalise posterior sense probabilities over all senses (incl noise)
sense.probs.sim.normalised = sapply(1:num.snippets, function(d){
  sense.probs.sim[,,d]/rowSums(sense.probs.sim[,,d])}, simplify = "array")

#Combine senses into the 3 true senses
sense.probs.sim.normalised.2 = sense.probs.sim.normalised[,1:2,]
sense.probs.sim.normalised.2[,1,] = apply(sense.probs.sim.normalised[,desired.labelling[[1]],,drop=FALSE], 3, rowSums)
sense.probs.sim.normalised.2[,2,] = apply(sense.probs.sim.normalised[,desired.labelling[[2]],,drop=FALSE], 3, rowSums)
sense.probs.sim.normalised = sapply(1:num.snippets, function(d){
  sense.probs.sim.normalised.2[,,d] / rowSums(sense.probs.sim.normalised.2[,,d])}, simplify = "array")
rm(sense.probs.sim.normalised.2); gc()

#Posterior mean sense probabilities
sense.probs.post.mean = apply(sense.probs.sim.normalised[(burn.in+1):num.samples,,], 2:3, mean)

#Index of snippets with either the riverbank or the institution bank true sense
sense.ids = levels(snippets.info$sense.id)
idx = which(snippets.info$sense.id %in% sense.ids[1:2])
names(idx) = idx #required for consistency with the Greek data

#Brier score
(BS = measures::multiclass.Brier(t(sense.probs.post.mean[,idx]), as.numeric(factor(snippets.info$sense.id[idx]))))



#############################################################################
# Brier scores - confidence intervals ------------------------------------- #
#############################################################################

alpha = 0.05
conf.level = 1-alpha/2

#Divide MCMC output into chunks, compute BS on each chunk, then compute CI around mean BS
num.chunks = 10
mcmc.blocks = split((burn.in+1):num.samples, sort(((burn.in+1):num.samples) %% num.chunks))

sense.probs.post.means = lapply(mcmc.blocks, function(i){apply(sense.probs.sim.normalised[paste(i),,], 2:3, mean)})

brier.scores = sapply(1:num.chunks, function(k){
  measures::multiclass.Brier(t(sense.probs.post.means[[k]][,idx]), as.numeric(factor(snippets.info$sense.id[idx])))})

BS.sd = sd(brier.scores)

dec.places = 4
noquote( paste(round(BS - qt(conf.level, num.chunks-1) * BS.sd / sqrt(num.chunks), dec.places), 
               round(BS + qt(conf.level, num.chunks-1) * BS.sd / sqrt(num.chunks), dec.places), sep = " - ") )
plot(brier.scores)



#############################################################################
# WAIC -------------------------------------------------------------------- #
#############################################################################

library(LaplacesDemon)
log.lik.by.snippet.sim = log(apply(sense.probs.sim, 1, colSums))
WAIC(log.lik.by.snippet.sim[,(burn.in+1):num.samples])

#alternative with diagnostics
library(loo)
log.lik.by.snippet.sim = log(apply(sense.probs.sim, 1, colSums))
waic(t(log.lik.by.snippet.sim[,(burn.in+1):num.samples]))
loo(t(log.lik.by.snippet.sim[,(burn.in+1):num.samples]))



#############################################################################
# BS vs WAIC -------------------------------------------------------------- #
#############################################################################

embedding.dims = c(50, 100, 200, 300)

#Hard-code results from runs:

BS.bank   = c(0.1391, 0.1388, 0.1333, 0.1431)
WAIC.bank = c(154693, 154459, 154154, 154069)

BS.kosmos   = c(0.3488, 0.3286, 0.3316, 0.3263)
WAIC.kosmos = c(137255, 136856, 136490, 136721)

BS.mus   = c(0.13460, 0.09342, 0.09940, 0.1005)
WAIC.mus = c(  19434,   19417,   19450,  19519)

#Plot just the WAIC
dev.off()
par(mfrow = c(1,3), oma = c(0.25,0,0,0), mar=c(2.5,2.5,2,1), cex.lab = 1, cex.axis = 1, cex.main = 1, mgp = c(1.5,0.5,0))
plot(embedding.dims, WAIC.bank, type = "o", pch = 18, col = "black",
     lty = 1, main = "bank", xlab = "embedding dimension", ylab = "WAIC")
plot(embedding.dims, WAIC.kosmos, type = "o", pch = 18, col = "black",
     lty = 1, main = "kosmos", xlab = "embedding dimension", ylab = "WAIC")
plot(embedding.dims, WAIC.mus, type = "o", pch = 18, col = "black",
     lty = 1, main = "mus", xlab = "embedding dimension", ylab = "WAIC")


#Plot WAIC and BS on the same graphs
dev.off()
# par(mfrow = c(1,3), oma = c(0.25,0,0,0.25), mar=c(2.5,2.5,2,2.5), cex.lab = 1, cex.axis = 1, cex.main = 1, mgp = c(1.5,0.5,0)) #Use this if legend is not required
par(mfrow = c(1,3), oma = c(2,0,0,0.25), mar=c(2.5,2.5,2,2.5), cex.lab = 1, cex.axis = 1, cex.main = 1, mgp = c(1.5,0.5,0))

plot(embedding.dims, WAIC.bank, type = "o", pch = 18, col = "black",  
     lty = 1, main = "bank", xlab = "embedding dimension", ylab = "WAIC")
par(new = TRUE)
plot(embedding.dims, BS.bank, type = "o", pch = 18, col = "blue",  
     lty = 2, axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, at = pretty(range(BS.bank)), col = "blue", col.axis = "blue")
mtext("BS", side=4, line=1.5, cex = 0.66, col = "blue")

plot(embedding.dims, WAIC.kosmos, type = "o", pch = 18, col = "black",  
     lty = 1, main = "kosmos", xlab = "embedding dimension", ylab = "WAIC")
par(new = TRUE)
plot(embedding.dims, BS.kosmos, type = "o", pch = 18, col = "blue",  
     lty = 2, axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, at = pretty(range(BS.kosmos)), col = "blue", col.axis = "blue")
mtext("BS", side=4, line=1.5, cex = 0.66, col = "blue")

plot(embedding.dims, WAIC.mus, type = "o", pch = 18, col = "black",  
     lty = 1, main = "mus", xlab = "embedding dimension", ylab = "WAIC")
par(new = TRUE)
plot(embedding.dims, BS.mus, type = "o", pch = 18, col = "blue",  
     lty = 2, axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, at = pretty(range(BS.mus)), col = "blue", col.axis = "blue")
mtext("BS", side=4, line=1.5, cex = 0.66, col = "blue")

par(mfrow = c(1,1), oma = c(0,0,0,0), mar = c(0.5,0,0,0), new = TRUE)
plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE) #overlay emply plot
legend("bottom", horiz = TRUE, bty = "n", legend = c("WAIC", "BS"),
       col = c("black", "blue"), pch = 18, lty = c(1,2), cex = 0.7, x.intersp = 0.7)

#Width = 9, Height = 3



#############################################################################
# Phi error bars - bank --------------------------------------------------- #
#############################################################################

# bank error bars - DiSC/EDiSC models - run model with 1 genre

conf.level = 0.95

# First load the phi.tilde results from either the DiSC or the EDiSC model and set burn-in

#Define sense order order corresponding to: riverbank, institution, noise
desired.labelling = list(c(1), c(2), c())

#Combine phi.tilde.sim into the 3 true senses
phi.tilde.sim.normalised = phi.tilde.sim[,1:2,,,drop=FALSE]
phi.tilde.sim.normalised[,1,,] = apply(phi.tilde.sim[,desired.labelling[[1]],,,drop=FALSE], c(3,4), rowSums)
phi.tilde.sim.normalised[,2,,] = apply(phi.tilde.sim[,desired.labelling[[2]],,,drop=FALSE], c(3,4), rowSums)
phi.tilde.sim.normalised = sapply(1:num.periods, function(t) {sapply(1:num.genres, function(g) {
  phi.tilde.sim.normalised[,,g,t] / rowSums(phi.tilde.sim[,,g,t]) }, simplify = "array") }, simplify = "array")


# Get the HPD intervals and posterior means for phi|W
phi.tilde.given.W.conf = array(dim = c(num.senses, 2, num.genres, num.periods),
                               dimnames = list(Sense = 1:2, Limit = c("lower","upper"),
                                               Genre = 1:num.genres, Time = 1:num.periods))
for (t in 1:num.periods){
  for (g in 1:num.genres){
    phi.tilde.given.W.conf[,,g,t] = HPDinterval( as.mcmc(phi.tilde.sim.normalised[(burn.in+1):num.samples,,g,t]), conf.level )
  }#for g
}#for t


#Phi.tilde posterior means
phi.tilde.given.W.mean = sapply(1:num.periods, function(t) {sapply(1:num.genres, function(g) {apply(phi.tilde.sim.normalised[(burn.in+1):num.samples,,g,t], 2, mean)}) }, simplify = "array")
dimnames(phi.tilde.given.W.mean) = list(Sense = 1:2, Genre = 1:num.genres, Time = 1:num.periods)


# save HPD intervals and means as relevant for DiSC or EDiSC, then calculate it again for the other
phi.tilde.given.W.conf.DiSC = phi.tilde.given.W.conf
phi.tilde.given.W.mean.DiSC = phi.tilde.given.W.mean
rm(phi.tilde.given.W.conf, phi.tilde.given.W.mean)

phi.tilde.given.W.conf.EDiSC = phi.tilde.given.W.conf
phi.tilde.given.W.mean.EDiSC = phi.tilde.given.W.mean
rm(phi.tilde.given.W.conf, phi.tilde.given.W.mean)


# Now load phi.tilde.sim from the phi|z run, and get the HPD intervals and posterior means
phi.tilde.given.z.conf = array(dim = c(num.senses, 2, num.genres, num.periods),
                               dimnames = list(Sense = 1:num.senses, Limit = c("lower","upper"),
                                               Genre = 1:num.genres, Time = 1:num.periods))
for (t in 1:num.periods){
  for (g in 1:num.genres){
    phi.tilde.given.z.conf[,,g,t] = HPDinterval( as.mcmc(phi.tilde.sim[,,g,t]), conf.level ) #burn-in should be zero
  }#for g
}#for t

#Phi.tilde posterior means
phi.tilde.given.z.mean = sapply(1:num.periods, function(t) {sapply(1:num.genres, function(g) {apply(phi.tilde.sim[,,g,t], 2, mean)}) }, simplify = "array")
dimnames(phi.tilde.given.z.mean) = list(Sense = 1:num.senses, Genre = 1:num.genres, Time = 1:num.periods)


#Index of snippets with either the riverbank or the institution bank true sense
sense.ids = levels(snippets.info$sense.id)
idx = which(snippets.info$sense.id %in% sense.ids[1:2])
names(idx) = idx #required for consistency with the Greek data
rm(sense.ids)

# expert annotation
z = as.numeric(factor(snippets.info$sense.id[idx]))

# empirical phi
sense.counts = table(factor(z, levels = 1:num.senses), factor(snippets$genre[idx], levels = 1:num.genres), 
                     factor(snippets$Time[idx], levels = 1:num.periods), dnn = c("Sense", "Genre", "Time"))

snippet.counts = table(factor(snippets$genre[idx], levels = 1:num.genres),
                       factor(snippets$Time[idx], levels = 1:num.periods),
                       dnn = c("Genre", "Time"))

phi.tilde = aperm(sapply(1:num.senses, function(k) {sense.counts[k,,] / snippet.counts}, simplify = "array"), c(3,1,2))
dimnames(phi.tilde) = list(Sense = 1:num.senses, Genre = 1:num.genres, Time = 1:num.periods)
phi.tilde[is.nan(phi.tilde)] = 0

#barplots
phi.all.genres = phi.tilde[,1,]
periods = levels(factor(snippets.info$Period))
dimnames(phi.all.genres) = list(Sense = c("river bank","institution bank"), Time = substr(periods, 1, 4))

#set up plot window
dev.off()
par(oma= c(2,0,0,0), mar=c(3,2.5,1.5,0), cex.lab = 0.8, cex.axis = 0.8, cex.main = 0.8, mgp = c(1.75,0,0))
# par(oma= c(2,0,0,0), mar=c(2.5,2.5,1.5,0), cex.lab = 0.8, cex.axis = 0.8, cex.main = 0.8, mgp = c(1,0,0))
disp = 0.3 #displace DiSC error bars this much to the left and EDiSC error bars this much to the right

#bars for all genres
plot.bar = barplot(phi.all.genres, beside = TRUE, space = c(0,0.5), col = c("green3","Orange2"), 
                   xlab = "Time period (snippet count)", yaxt = 'n', main = "All genres", ylim = c(0,1))
axis(2, mgp = c(1.5,0.5,0)); title(ylab = "Sense prevalence", mgp = c(1.5,0.5,0))
axis(1, at = apply(plot.bar, 2, mean), labels = paste0("(",snippet.counts[1,],")"), tick = FALSE, mgp = c(1.5,0.75,0), )

#errors for DiSC
arrows(x0 = as.vector(plot.bar)-disp, 
       y0 = as.vector(phi.tilde.given.W.conf.DiSC[,1,1,]), y1 = as.vector(phi.tilde.given.W.conf.DiSC[,2,1,]),
       code = 3, angle = 90, length = 0.015, col = "red", lty = 1, lwd = 2)
points(as.vector(plot.bar)-disp, as.vector(phi.tilde.given.W.mean.DiSC[,1,]), pch = 16, cex = 0.75, col = "red")

#errors for EDiSC
arrows(x0 = as.vector(plot.bar)+disp, 
       y0 = as.vector(phi.tilde.given.W.conf.EDiSC[,1,1,]), y1 = as.vector(phi.tilde.given.W.conf.EDiSC[,2,1,]),
       code = 3, angle = 90, length = 0.015, col = "blue", lty = 1, lwd = 2)
points(as.vector(plot.bar)+disp, as.vector(phi.tilde.given.W.mean.EDiSC[,1,]), pch = 16, cex = 0.75, col = "blue")

#errors for phi|z
arrows(x0 = as.vector(plot.bar),
       y0 = as.vector(phi.tilde.given.z.conf[,1,1,]), y1 = as.vector(phi.tilde.given.z.conf[,2,1,]),
       code = 3, angle = 90, length = 0.015, col = "black", lty = 2, lwd = 1)
points(as.vector(plot.bar), as.vector(phi.tilde.given.z.mean[,1,]), pch = 16, cex = 0.75, col = "black")

#legend
par(mfrow = c(1,1), oma = c(0,0,0,0), mar = c(0.5,0,0,0), new = TRUE)
plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE) #overlay emply plot
legend("bottom", legend = c(rownames(phi.all.genres), expression(paste("DiSC ",tilde(phi),"|W")), 
                            expression(paste("DiSC ",tilde(phi),"|z")), expression(paste("EDiSC ",tilde(phi),"|W"))), 
       fill = c("green3","Orange2",NA,NA,NA), col = c("green3","Orange2","red","black","blue"), 
       border = c("black","black",NA,NA,NA),  horiz = TRUE, bty = "n", lty = c(NA,NA,1,2,1), lwd = c(NA,NA,2,1,2),
       cex = 0.8, x.intersp = c(-1.25,-1.25,0.75,0.75,0.75), text.width = c(0.12,0.135,0.065,0.06,0.065))

#Width = 7, Height = 3.5



#############################################################################
# Phi error bars - mus ---------------------------------------------------- #
#############################################################################

conf.level = 0.95


# First load the phi.tilde results from either the DiSC or the EDiSC model and set burn-in

#Define sense order order corresponding to: mouse, muscle, mussel, noise
desired.labelling = list(c(1), c(2), c(3), c())

#Combine phi.tilde.sim into the 3 true senses
phi.tilde.sim.normalised = phi.tilde.sim[,1:3,,,drop=FALSE]
phi.tilde.sim.normalised[,1,,] = apply(phi.tilde.sim[,desired.labelling[[1]],,,drop=FALSE], c(3,4), rowSums)
phi.tilde.sim.normalised[,2,,] = apply(phi.tilde.sim[,desired.labelling[[2]],,,drop=FALSE], c(3,4), rowSums)
phi.tilde.sim.normalised[,3,,] = apply(phi.tilde.sim[,desired.labelling[[3]],,,drop=FALSE], c(3,4), rowSums)
phi.tilde.sim.normalised = sapply(1:num.periods, function(t) {sapply(1:num.genres, function(g) {
  phi.tilde.sim.normalised[,,g,t] / rowSums(phi.tilde.sim[,,g,t]) }, simplify = "array") }, simplify = "array")


# Get the HPD intervals and posterior means for phi|W
phi.tilde.given.W.conf = array(dim = c(3, 2, num.genres, num.periods),
                               dimnames = list(Sense = 1:3, Limit = c("lower","upper"),
                                               Genre = 1:num.genres, Time = 1:num.periods))
for (t in 1:num.periods){
  for (g in 1:num.genres){
    phi.tilde.given.W.conf[,,g,t] = HPDinterval( as.mcmc(phi.tilde.sim.normalised[(burn.in+1):num.samples,,g,t]), conf.level )
  }#for g
}#for t

#Phi.tilde posterior means
phi.tilde.given.W.mean = sapply(1:num.periods, function(t) {sapply(1:num.genres, function(g) {apply(phi.tilde.sim.normalised[(burn.in+1):num.samples,,g,t], 2, mean)}) }, simplify = "array")
dimnames(phi.tilde.given.W.mean) = list(Sense = 1:3, Genre = 1:num.genres, Time = 1:num.periods)


# save HPD intervals and means as relevant for DiSC or EDiSC, then calculate it again for the other
phi.tilde.given.W.conf.DiSC = phi.tilde.given.W.conf
phi.tilde.given.W.mean.DiSC = phi.tilde.given.W.mean
rm(phi.tilde.given.W.conf, phi.tilde.given.W.mean)

phi.tilde.given.W.conf.EDiSC = phi.tilde.given.W.conf
phi.tilde.given.W.mean.EDiSC = phi.tilde.given.W.mean
rm(phi.tilde.given.W.conf, phi.tilde.given.W.mean)


# Now load snippets and phi.tilde.sim from the phi|z run, and get the HPD intervals and posterior means
phi.tilde.given.z.conf = array(dim = c(num.senses, 2, num.genres, num.periods),
                               dimnames = list(Sense = 1:num.senses, Limit = c("lower","upper"),
                                               Genre = 1:num.genres, Time = 1:num.periods))
for (t in 1:num.periods){
  for (g in 1:num.genres){
    phi.tilde.given.z.conf[,,g,t] = HPDinterval( as.mcmc(phi.tilde.sim[,,g,t]), conf.level ) #burn-in should be zero
  }#for g
}#for t

#Phi.tilde posterior means
phi.tilde.given.z.mean = sapply(1:num.periods, function(t) {sapply(1:num.genres, function(g) {apply(phi.tilde.sim[,,g,t], 2, mean)}) }, simplify = "array")
dimnames(phi.tilde.given.z.mean) = list(Sense = 1:num.senses, Genre = 1:num.genres, Time = 1:num.periods)


#Index of non-empty snippets
# idx = which(snippet.lengths > 0)
# idx = which(snippet.lengths > 0 & snippets.info$disambiguation == 1) #of type collocates
sense.ids = levels(factor(snippets.info$sense.id))
idx = which(snippets.info$sense.id %in% sense.ids[1:3])
names(idx) = idx 
rm(sense.ids)

# expert annotation
z = as.numeric(factor(snippets.info$sense.id[idx]))

# empirical phi
sense.counts = table(factor(z, levels = 1:num.senses), factor(snippets$genre[idx], levels = 1:num.genres), 
                     factor(snippets$Time[idx], levels = 1:num.periods), dnn = c("Sense", "Genre", "Time"))

snippet.counts = table(factor(snippets$genre[idx], levels = 1:num.genres),
                       factor(snippets$Time[idx], levels = 1:num.periods),
                       dnn = c("Genre", "Time"))

phi.tilde = aperm(sapply(1:num.senses, function(k) {sense.counts[k,,] / snippet.counts}, simplify = "array"), c(3,1,2))
dimnames(phi.tilde) = list(Sense = 1:num.senses, Genre = 1:num.genres, Time = 1:num.periods)
phi.tilde[is.nan(phi.tilde)] = 0


#barplots
phi.technical = phi.tilde[,1,]
phi.non.technical = phi.tilde[,2,]
# dimnames(phi.technical) = list(Sense = c("Mouse","Muscle","Mussel"), Time = c((-5):(-1),1:4))
# dimnames(phi.non.technical) = list(Sense = c("Mouse","Muscle","Mussel"), Time = c((-5):(-1),1:4))
dimnames(phi.technical) = list(Sense = c("Mouse","Muscle","Mussel"), Time = c(paste(5:1, "BC"), paste(1:4, "AD")))
dimnames(phi.non.technical) = list(Sense = c("Mouse","Muscle","Mussel"), Time = c(paste(5:1, "BC"), paste(1:4, "AD")))


#set up plot window
dev.off()
par(mfrow = c(2,1), oma= c(2,0,0,0), mar=c(3,2.5,1.5,0), cex.lab = 0.8, cex.axis = 0.8, cex.main = 0.8, mgp = c(1.75,0,0))
disp = 0.3 #displace DiSC error bars this much to the left and EDiSC error bars this much to the right


#technical genre
plot.tech = barplot(phi.technical, beside = TRUE, space = c(0,0.5), col = c("green3","Orange2","Purple1"), 
                    xlab = "Century (snippet count)", yaxt = 'n', main = "Technical", ylim = c(0,1))
axis(2, mgp = c(1.5,0.5,0)); title(ylab = "Sense prevalence", mgp = c(1.5,0.5,0))
axis(1, at = plot.tech[2,], labels = paste0("(",snippet.counts[1,],")"), tick = FALSE, mgp = c(1.5,0.75,0), )

#errors for DiSC
arrows(x0 = as.vector(plot.tech)-disp, 
       y0 = as.vector(phi.tilde.given.W.conf.DiSC[,1,1,]), y1 = as.vector(phi.tilde.given.W.conf.DiSC[,2,1,]),
       code = 3, angle = 90, length = 0.015, col = "red", lty = 1, lwd = 2)
points(as.vector(plot.tech)-disp, as.vector(phi.tilde.given.W.mean.DiSC[,1,]), pch = 16, cex = 0.75, col = "red")

#errors for EDiSC
arrows(x0 = as.vector(plot.tech)+disp, 
       y0 = as.vector(phi.tilde.given.W.conf.EDiSC[,1,1,]), y1 = as.vector(phi.tilde.given.W.conf.EDiSC[,2,1,]),
       code = 3, angle = 90, length = 0.015, col = "blue", lty = 1, lwd = 2)
points(as.vector(plot.tech)+disp, as.vector(phi.tilde.given.W.mean.EDiSC[,1,]), pch = 16, cex = 0.75, col = "blue")

#errors for phi|z
arrows(x0 = as.vector(plot.tech),
       y0 = as.vector(phi.tilde.given.z.conf[,1,1,]), y1 = as.vector(phi.tilde.given.z.conf[,2,1,]),
       code = 3, angle = 90, length = 0.015, col = "black", lty = 2, lwd = 1)
points(as.vector(plot.tech), as.vector(phi.tilde.given.z.mean[,1,]), pch = 16, cex = 0.75, col = "black")


#non-technical genre
plot.non.tech = barplot(phi.non.technical, beside = TRUE, space = c(0,0.5), col = c("green3","Orange2","Purple1"), 
                        xlab = "Century (snippet count)", yaxt = 'n', main = "Non-technical", ylim = c(0,1))
axis(2, mgp = c(1.5,0.5,0)); title(ylab = "Sense prevalence", mgp = c(1.5,0.5,0))
axis(1, at = plot.non.tech[2,], labels = paste0("(",snippet.counts[2,],")"), tick = FALSE, mgp = c(1.5,0.75,0), )

#errors for DiSC
arrows(x0 = as.vector(plot.non.tech)-disp, 
       y0 = as.vector(phi.tilde.given.W.conf.DiSC[,1,2,]), y1 = as.vector(phi.tilde.given.W.conf.DiSC[,2,2,]),
       code = 3, angle = 90, length = 0.015, col = "red", lty = 1, lwd = 2)
points(as.vector(plot.non.tech)-disp, as.vector(phi.tilde.given.W.mean.DiSC[,2,]), pch = 16, cex = 0.75, col = "red")

#errors for EDiSC
arrows(x0 = as.vector(plot.non.tech)+disp, 
       y0 = as.vector(phi.tilde.given.W.conf.EDiSC[,1,2,]), y1 = as.vector(phi.tilde.given.W.conf.EDiSC[,2,2,]),
       code = 3, angle = 90, length = 0.015, col = "blue", lty = 1, lwd = 2)
points(as.vector(plot.non.tech)+disp, as.vector(phi.tilde.given.W.mean.EDiSC[,2,]), pch = 16, cex = 0.75, col = "blue")

#errors for phi|z
arrows(x0 = as.vector(plot.non.tech),
       y0 = as.vector(phi.tilde.given.z.conf[,1,2,]), y1 = as.vector(phi.tilde.given.z.conf[,2,2,]),
       code = 3, angle = 90, length = 0.015, col = "black", lty = 2, lwd = 1)
points(as.vector(plot.non.tech), as.vector(phi.tilde.given.z.mean[,2,]), pch = 16, cex = 0.75, col = "black")


#legend
par(mfrow = c(1,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE) #overlay emply plot
legend(x=0.05,y=0.01, legend = c(rownames(phi.technical), expression(paste("DiSC ",tilde(phi),"|W")),
                                 expression(paste("DiSC ",tilde(phi),"|z")), expression(paste("EDiSC ",tilde(phi),"|W"))),
       fill = c("green3","Orange2","Purple1",NA,NA,NA), col = c("green3","Orange2","Purple1","red","black","blue"),
       border = c("black","black","black",NA,NA,NA), lty = c(NA,NA,NA,1,2,1), lwd = c(NA,NA,NA,2,1,2), bty = "n", horiz = TRUE,
       cex = 0.9, x.intersp = c(-1.75,-1.75,-1.75,0.25,0.25,0.25), text.width = c(0.065,0.065,0.065,0.065,0.06,0.065))

#Width = 7, Height = 6



#############################################################################
# Phi error bars - harmonia ----------------------------------------------- #
#############################################################################

conf.level = 0.95


# First load the phi.tilde results from either the DiSC or the EDiSC model and set burn-in

#Define sense order order corresponding to: abstract, concrete, musical, noise
desired.labelling = list(c(1), c(2), c(3,4), c())

#Combine phi.tilde.sim into the 3 true senses
phi.tilde.sim.normalised = phi.tilde.sim[,1:3,,,drop=FALSE]
phi.tilde.sim.normalised[,1,,] = apply(phi.tilde.sim[,desired.labelling[[1]],,,drop=FALSE], c(3,4), rowSums)
phi.tilde.sim.normalised[,2,,] = apply(phi.tilde.sim[,desired.labelling[[2]],,,drop=FALSE], c(3,4), rowSums)
phi.tilde.sim.normalised[,3,,] = apply(phi.tilde.sim[,desired.labelling[[3]],,,drop=FALSE], c(3,4), rowSums)
phi.tilde.sim.normalised = sapply(1:num.periods, function(t) {sapply(1:num.genres, function(g) {
  phi.tilde.sim.normalised[,,g,t] / rowSums(phi.tilde.sim[,,g,t]) }, simplify = "array") }, simplify = "array")


# Get the HPD intervals and posterior means for phi|W
phi.tilde.given.W.conf = array(dim = c(3, 2, num.genres, num.periods),
                               dimnames = list(Sense = 1:3, Limit = c("lower","upper"),
                                               Genre = 1:num.genres, Time = 1:num.periods))
for (t in 1:num.periods){
  for (g in 1:num.genres){
    phi.tilde.given.W.conf[,,g,t] = HPDinterval( as.mcmc(phi.tilde.sim.normalised[(burn.in+1):num.samples,,g,t]), conf.level )
  }#for g
}#for t

#Phi.tilde posterior means
phi.tilde.given.W.mean = sapply(1:num.periods, function(t) {sapply(1:num.genres, function(g) {apply(phi.tilde.sim.normalised[(burn.in+1):num.samples,,g,t], 2, mean)}) }, simplify = "array")
dimnames(phi.tilde.given.W.mean) = list(Sense = 1:3, Genre = 1:num.genres, Time = 1:num.periods)


# save HPD intervals and means as relevant for DiSC or EDiSC, then calculate it again for the other
phi.tilde.given.W.conf.DiSC = phi.tilde.given.W.conf
phi.tilde.given.W.mean.DiSC = phi.tilde.given.W.mean
rm(phi.tilde.given.W.conf, phi.tilde.given.W.mean)

phi.tilde.given.W.conf.EDiSC = phi.tilde.given.W.conf
phi.tilde.given.W.mean.EDiSC = phi.tilde.given.W.mean
rm(phi.tilde.given.W.conf, phi.tilde.given.W.mean)


# Now load snippets and phi.tilde.sim from the phi|z run, and get the HPD intervals and posterior means
phi.tilde.given.z.conf = array(dim = c(num.senses, 2, num.genres, num.periods),
                               dimnames = list(Sense = 1:num.senses, Limit = c("lower","upper"),
                                               Genre = 1:num.genres, Time = 1:num.periods))
for (t in 1:num.periods){
  for (g in 1:num.genres){
    phi.tilde.given.z.conf[,,g,t] = HPDinterval( as.mcmc(phi.tilde.sim[,,g,t]), conf.level ) #burn-in should be zero
  }#for g
}#for t

#Phi.tilde posterior means
phi.tilde.given.z.mean = sapply(1:num.periods, function(t) {sapply(1:num.genres, function(g) {apply(phi.tilde.sim[,,g,t], 2, mean)}) }, simplify = "array")
dimnames(phi.tilde.given.z.mean) = list(Sense = 1:num.senses, Genre = 1:num.genres, Time = 1:num.periods)


#Index of non-empty snippets
# idx = which(snippet.lengths > 0)
# idx = which(snippet.lengths > 0 & snippets.info$disambiguation == 1) #of type collocates
sense.ids = levels(factor(snippets.info$sense.id))
idx = which(snippets.info$sense.id %in% sense.ids[1:3])
names(idx) = idx 
rm(sense.ids)

# expert annotation
z = as.numeric(factor(snippets.info$sense.id[idx]))

# empirical phi
sense.counts = table(factor(z, levels = 1:num.senses), factor(snippets$genre[idx], levels = 1:num.genres), 
                     factor(snippets$Time[idx], levels = 1:num.periods), dnn = c("Sense", "Genre", "Time"))

snippet.counts = table(factor(snippets$genre[idx], levels = 1:num.genres),
                       factor(snippets$Time[idx], levels = 1:num.periods),
                       dnn = c("Genre", "Time"))

phi.tilde = aperm(sapply(1:num.senses, function(k) {sense.counts[k,,] / snippet.counts}, simplify = "array"), c(3,1,2))
dimnames(phi.tilde) = list(Sense = 1:num.senses, Genre = 1:num.genres, Time = 1:num.periods)
phi.tilde[is.nan(phi.tilde)] = 0


#barplots
phi.technical = phi.tilde[,1,]
phi.non.technical = phi.tilde[,2,]
# dimnames(phi.technical) = list(Sense = c("Abstract","Concrete","Musical"), Time = c((-8):(-1),1:4))
# dimnames(phi.non.technical) = list(Sense = c("Abstract","Concrete","Musical"), Time = c((-8):(-1),1:4))
dimnames(phi.technical) = list(Sense = c("Abstract","Concrete","Musical"), Time = c(paste(8:1, "BC"), paste(1:4, "AD")))
dimnames(phi.non.technical) = list(Sense = c("Abstract","Concrete","Musical"), Time = c(paste(8:1, "BC"), paste(1:4, "AD")))


#set up plot window
dev.off()
par(mfrow = c(2,1), oma= c(2,0,0,0), mar=c(3,2.5,1.5,0), cex.lab = 0.8, cex.axis = 0.8, cex.main = 0.8, mgp = c(1.75,0,0))
disp = 0.3 #displace DiSC error bars this much to the left and EDiSC error bars this much to the right


#technical genre
plot.tech = barplot(phi.technical, beside = TRUE, space = c(0,0.5), col = c("green3","Orange2","Purple1"), 
                    xlab = "Century (snippet count)", yaxt = 'n', main = "Technical", ylim = c(0,1))
axis(2, mgp = c(1.5,0.5,0)); title(ylab = "Sense prevalence", mgp = c(1.5,0.5,0))
axis(1, at = plot.tech[2,], labels = paste0("(",snippet.counts[1,],")"), tick = FALSE, mgp = c(1.5,0.75,0), )

#errors for DiSC
arrows(x0 = as.vector(plot.tech)-disp, 
       y0 = as.vector(phi.tilde.given.W.conf.DiSC[,1,1,]), y1 = as.vector(phi.tilde.given.W.conf.DiSC[,2,1,]),
       code = 3, angle = 90, length = 0.015, col = "red", lty = 1, lwd = 2)
points(as.vector(plot.tech)-disp, as.vector(phi.tilde.given.W.mean.DiSC[,1,]), pch = 16, cex = 0.75, col = "red")

#errors for EDiSC
arrows(x0 = as.vector(plot.tech)+disp, 
       y0 = as.vector(phi.tilde.given.W.conf.EDiSC[,1,1,]), y1 = as.vector(phi.tilde.given.W.conf.EDiSC[,2,1,]),
       code = 3, angle = 90, length = 0.015, col = "blue", lty = 1, lwd = 2)
points(as.vector(plot.tech)+disp, as.vector(phi.tilde.given.W.mean.EDiSC[,1,]), pch = 16, cex = 0.75, col = "blue")

#errors for phi|z
arrows(x0 = as.vector(plot.tech),
       y0 = as.vector(phi.tilde.given.z.conf[,1,1,]), y1 = as.vector(phi.tilde.given.z.conf[,2,1,]),
       code = 3, angle = 90, length = 0.015, col = "black", lty = 2, lwd = 1)
points(as.vector(plot.tech), as.vector(phi.tilde.given.z.mean[,1,]), pch = 16, cex = 0.75, col = "black")


#non-technical genre
plot.non.tech = barplot(phi.non.technical, beside = TRUE, space = c(0,0.5), col = c("green3","Orange2","Purple1"), 
                        xlab = "Century (snippet count)", yaxt = 'n', main = "Non-technical", ylim = c(0,1))
axis(2, mgp = c(1.5,0.5,0)); title(ylab = "Sense prevalence", mgp = c(1.5,0.5,0))
axis(1, at = plot.non.tech[2,], labels = paste0("(",snippet.counts[2,],")"), tick = FALSE, mgp = c(1.5,0.75,0), )

#errors for DiSC
arrows(x0 = as.vector(plot.non.tech)-disp, 
       y0 = as.vector(phi.tilde.given.W.conf.DiSC[,1,2,]), y1 = as.vector(phi.tilde.given.W.conf.DiSC[,2,2,]),
       code = 3, angle = 90, length = 0.015, col = "red", lty = 1, lwd = 2)
points(as.vector(plot.non.tech)-disp, as.vector(phi.tilde.given.W.mean.DiSC[,2,]), pch = 16, cex = 0.75, col = "red")

#errors for EDiSC
arrows(x0 = as.vector(plot.non.tech)+disp, 
       y0 = as.vector(phi.tilde.given.W.conf.EDiSC[,1,2,]), y1 = as.vector(phi.tilde.given.W.conf.EDiSC[,2,2,]),
       code = 3, angle = 90, length = 0.015, col = "blue", lty = 1, lwd = 2)
points(as.vector(plot.non.tech)+disp, as.vector(phi.tilde.given.W.mean.EDiSC[,2,]), pch = 16, cex = 0.75, col = "blue")

#errors for phi|z
arrows(x0 = as.vector(plot.non.tech),
       y0 = as.vector(phi.tilde.given.z.conf[,1,2,]), y1 = as.vector(phi.tilde.given.z.conf[,2,2,]),
       code = 3, angle = 90, length = 0.015, col = "black", lty = 2, lwd = 1)
points(as.vector(plot.non.tech), as.vector(phi.tilde.given.z.mean[,2,]), pch = 16, cex = 0.75, col = "black")


#legend
par(mfrow = c(1,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE) #overlay emply plot
legend(x=0.05,y=0.01, legend = c(rownames(phi.technical), expression(paste("DiSC ",tilde(phi),"|W")),
                                 expression(paste("DiSC ",tilde(phi),"|z")), expression(paste("EDiSC ",tilde(phi),"|W"))),
       fill = c("green3","Orange2","Purple1",NA,NA,NA), col = c("green3","Orange2","Purple1","red","black","blue"),
       border = c("black","black","black",NA,NA,NA), lty = c(NA,NA,NA,1,2,1), lwd = c(NA,NA,NA,2,1,2), bty = "n", horiz = TRUE,
       cex = 0.9, x.intersp = c(-1.75,-1.75,-1.75,0.25,0.25,0.25), text.width = c(0.07,0.075,0.07,0.065,0.06,0.065))

#Width = 7, Height = 6



#############################################################################
# Phi error bars - kosmos ------------------------------------------------- #
#############################################################################

conf.level = 0.95

# First load the phi.tilde results from either the DiSC or the EDiSC model and set burn-in

#Define sense order order corresponding to: decoration, order, world, noise
desired.labelling = list(c(2), c(1), c(4,3), c())

#Combine phi.tilde.sim into the 3 true senses
phi.tilde.sim.normalised = phi.tilde.sim[,1:3,,,drop=FALSE]
phi.tilde.sim.normalised[,1,,] = apply(phi.tilde.sim[,desired.labelling[[1]],,,drop=FALSE], c(3,4), rowSums)
phi.tilde.sim.normalised[,2,,] = apply(phi.tilde.sim[,desired.labelling[[2]],,,drop=FALSE], c(3,4), rowSums)
phi.tilde.sim.normalised[,3,,] = apply(phi.tilde.sim[,desired.labelling[[3]],,,drop=FALSE], c(3,4), rowSums)
phi.tilde.sim.normalised = sapply(1:num.periods, function(t) {sapply(1:num.genres, function(g) {
  phi.tilde.sim.normalised[,,g,t] / rowSums(phi.tilde.sim[,,g,t]) }, simplify = "array") }, simplify = "array")


# Get the HPD intervals and posterior means for phi|W
phi.tilde.given.W.conf = array(dim = c(3, 2, num.genres, num.periods),
                               dimnames = list(Sense = 1:3, Limit = c("lower","upper"),
                                               Genre = 1:num.genres, Time = 1:num.periods))
for (t in 1:num.periods){
  for (g in 1:num.genres){
    phi.tilde.given.W.conf[,,g,t] = HPDinterval( as.mcmc(phi.tilde.sim.normalised[(burn.in+1):num.samples,,g,t]), conf.level )
  }#for g
}#for t

#Phi.tilde posterior means
phi.tilde.given.W.mean = sapply(1:num.periods, function(t) {sapply(1:num.genres, function(g) {apply(phi.tilde.sim.normalised[(burn.in+1):num.samples,,g,t], 2, mean)}) }, simplify = "array")
dimnames(phi.tilde.given.W.mean) = list(Sense = 1:3, Genre = 1:num.genres, Time = 1:num.periods)


# save HPD intervals and means as relevant for DiSC or EDiSC, then calculate it again for the other
phi.tilde.given.W.conf.DiSC = phi.tilde.given.W.conf
phi.tilde.given.W.mean.DiSC = phi.tilde.given.W.mean
rm(phi.tilde.given.W.conf, phi.tilde.given.W.mean)

phi.tilde.given.W.conf.EDiSC = phi.tilde.given.W.conf
phi.tilde.given.W.mean.EDiSC = phi.tilde.given.W.mean
rm(phi.tilde.given.W.conf, phi.tilde.given.W.mean)


# Now load snippets and phi.tilde.sim from the phi|z run, and get the HPD intervals and posterior means
phi.tilde.given.z.conf = array(dim = c(num.senses, 2, num.genres, num.periods),
                               dimnames = list(Sense = 1:num.senses, Limit = c("lower","upper"),
                                               Genre = 1:num.genres, Time = 1:num.periods))
for (t in 1:num.periods){
  for (g in 1:num.genres){
    phi.tilde.given.z.conf[,,g,t] = HPDinterval( as.mcmc(phi.tilde.sim[,,g,t]), conf.level ) #burn-in should be zero
  }#for g
}#for t

#Phi.tilde posterior means
phi.tilde.given.z.mean = sapply(1:num.periods, function(t) {sapply(1:num.genres, function(g) {apply(phi.tilde.sim[,,g,t], 2, mean)}) }, simplify = "array")
dimnames(phi.tilde.given.z.mean) = list(Sense = 1:num.senses, Genre = 1:num.genres, Time = 1:num.periods)


# #Index of non-empty snippets
# idx = which(snippet.lengths > 0)
# idx = which(snippet.lengths > 0 & snippets.info$disambiguation == 1) #of type collocates
sense.ids = levels(factor(snippets.info$sense.id))
idx = which(snippets.info$sense.id %in% sense.ids[1:3])
names(idx) = idx 
rm(sense.ids)

# expert annotation
z = as.numeric(factor(snippets.info$sense.id[idx]))

# empirical phi
sense.counts = table(factor(z, levels = 1:num.senses), factor(snippets$genre[idx], levels = 1:num.genres), 
                     factor(snippets$Time[idx], levels = 1:num.periods), dnn = c("Sense", "Genre", "Time"))

snippet.counts = table(factor(snippets$genre[idx], levels = 1:num.genres),
                       factor(snippets$Time[idx], levels = 1:num.periods),
                       dnn = c("Genre", "Time"))

phi.tilde = aperm(sapply(1:num.senses, function(k) {sense.counts[k,,] / snippet.counts}, simplify = "array"), c(3,1,2))
dimnames(phi.tilde) = list(Sense = 1:num.senses, Genre = 1:num.genres, Time = 1:num.periods)
phi.tilde[is.nan(phi.tilde)] = 0


#barplots
phi.narrative = phi.tilde[,1,]
phi.non.narrative = phi.tilde[,2,]
# dimnames(phi.narrative) = list(Sense = c("Decoration","Order","World"), Time = c((-7):(-1),1:2))
# dimnames(phi.non.narrative) = list(Sense = c("Decoration","Order","World"), Time = c((-7):(-1),1:2))
dimnames(phi.narrative) = list(Sense = c("Decoration","Order","World"), Time = c(paste(7:1, "BC"), paste(1:2, "AD")))
dimnames(phi.non.narrative) = list(Sense = c("Decoration","Order","World"), Time = c(paste(7:1, "BC"), paste(1:2, "AD")))


#set up plot window
dev.off()
par(mfrow = c(2,1), oma= c(2,0,0,0), mar=c(3,2.5,1.5,0), cex.lab = 0.8, cex.axis = 0.8, cex.main = 0.8, mgp = c(1.75,0,0))
disp = 0.3 #displace DiSC error bars this much to the left and EDiSC error bars this much to the right


#narrative genre
plot.narr = barplot(phi.narrative, beside = TRUE, space = c(0,0.5), col = c("green3","Orange2","Purple1"), 
                    xlab = "Century (snippet count)", yaxt = 'n', main = "Narrative", ylim = c(0,1))
axis(2, mgp = c(1.5,0.5,0)); title(ylab = "Sense prevalence", mgp = c(1.5,0.5,0))
axis(1, at = plot.narr[2,], labels = paste0("(",snippet.counts[1,],")"), tick = FALSE, mgp = c(1.5,0.75,0), )

#errors for DiSC
arrows(x0 = as.vector(plot.narr)-disp, 
       y0 = as.vector(phi.tilde.given.W.conf.DiSC[,1,1,]), y1 = as.vector(phi.tilde.given.W.conf.DiSC[,2,1,]),
       code = 3, angle = 90, length = 0.015, col = "red", lty = 1, lwd = 2)
points(as.vector(plot.narr)-disp, as.vector(phi.tilde.given.W.mean.DiSC[,1,]), pch = 16, cex = 0.75, col = "red")

#errors for EDiSC
arrows(x0 = as.vector(plot.narr)+disp, 
       y0 = as.vector(phi.tilde.given.W.conf.EDiSC[,1,1,]), y1 = as.vector(phi.tilde.given.W.conf.EDiSC[,2,1,]),
       code = 3, angle = 90, length = 0.015, col = "blue", lty = 1, lwd = 2)
points(as.vector(plot.narr)+disp, as.vector(phi.tilde.given.W.mean.EDiSC[,1,]), pch = 16, cex = 0.75, col = "blue")

#errors for phi|z
arrows(x0 = as.vector(plot.narr),
       y0 = as.vector(phi.tilde.given.z.conf[,1,1,]), y1 = as.vector(phi.tilde.given.z.conf[,2,1,]),
       code = 3, angle = 90, length = 0.015, col = "black", lty = 2, lwd = 1)
points(as.vector(plot.narr), as.vector(phi.tilde.given.z.mean[,1,]), pch = 16, cex = 0.75, col = "black")


#non-narrative genre
plot.non.narr = barplot(phi.non.narrative, beside = TRUE, space = c(0,0.5), col = c("green3","Orange2","Purple1"), 
                        xlab = "Century (snippet count)", yaxt = 'n', main = "Non-narrative", ylim = c(0,1))
axis(2, mgp = c(1.5,0.5,0)); title(ylab = "Sense prevalence", mgp = c(1.5,0.5,0))
axis(1, at = plot.non.narr[2,], labels = paste0("(",snippet.counts[2,],")"), tick = FALSE, mgp = c(1.5,0.75,0), )

#errors for DiSC
arrows(x0 = as.vector(plot.non.narr)-disp, 
       y0 = as.vector(phi.tilde.given.W.conf.DiSC[,1,2,]), y1 = as.vector(phi.tilde.given.W.conf.DiSC[,2,2,]),
       code = 3, angle = 90, length = 0.015, col = "red", lty = 1, lwd = 2)
points(as.vector(plot.non.narr)-disp, as.vector(phi.tilde.given.W.mean.DiSC[,2,]), pch = 16, cex = 0.75, col = "red")

#errors for EDiSC
arrows(x0 = as.vector(plot.non.narr)+disp, 
       y0 = as.vector(phi.tilde.given.W.conf.EDiSC[,1,2,]), y1 = as.vector(phi.tilde.given.W.conf.EDiSC[,2,2,]),
       code = 3, angle = 90, length = 0.015, col = "blue", lty = 1, lwd = 2)
points(as.vector(plot.non.narr)+disp, as.vector(phi.tilde.given.W.mean.EDiSC[,2,]), pch = 16, cex = 0.75, col = "blue")

#errors for phi|z
arrows(x0 = as.vector(plot.non.narr),
       y0 = as.vector(phi.tilde.given.z.conf[,1,2,]), y1 = as.vector(phi.tilde.given.z.conf[,2,2,]),
       code = 3, angle = 90, length = 0.015, col = "black", lty = 2, lwd = 1)
points(as.vector(plot.non.narr), as.vector(phi.tilde.given.z.mean[,2,]), pch = 16, cex = 0.75, col = "black")


#legend
par(mfrow = c(1,1), oma = c(0,0,0,0), mar = c(0,0,0,0), new = TRUE)
plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE) #overlay emply plot
legend(x=0.05,y=0.01, legend = c(rownames(phi.narrative), expression(paste("DiSC ",tilde(phi),"|W")),
                            expression(paste("DiSC ",tilde(phi),"|z")), expression(paste("EDiSC ",tilde(phi),"|W"))),
       fill = c("green3","Orange2","Purple1",NA,NA,NA), col = c("green3","Orange2","Purple1","red","black","blue"),
       border = c("black","black","black",NA,NA,NA), lty = c(NA,NA,NA,1,2,1), lwd = c(NA,NA,NA,2,1,2), bty = "n", horiz = TRUE,
       cex = 0.9, x.intersp = c(-1.75,-1.75,-1.75,0.25,0.25,0.25), text.width = c(0.1,0.05,0.05,0.065,0.06,0.065))

#Width = 7, Height = 6



#############################################################################
# Metastable states plot for paper ---------------------------------------- #
#############################################################################

#Stan NUTS on EDiSC (M=300) for "mus"; seed 300; 105k sims, 5k warmup

desired.labelling = c(2,1,3)

#Posterior normalised sense probabilities
sense.probs.sim.normalised = sapply(1:num.snippets, function(d){
  sense.probs.sim[,,d]/rowSums(sense.probs.sim[,,d])}, simplify = "array")
sense.probs.sim.normalised[] = sense.probs.sim.normalised[,desired.labelling,]

#Index of non-empty snippets
# idx = which(snippet.lengths > 0)
idx = which(snippet.lengths > 0 & snippets.info$disambiguation == 1) #of type collocates

#Split mcmc into chunks
num.chunks = 100
mcmc.blocks = split(1:num.samples, sort((1:num.samples) %% num.chunks))

#Compute Brier scores on each chunk
sense.probs.post.means = lapply(mcmc.blocks, function(i){apply(sense.probs.sim.normalised[paste(i),,], 2:3, mean)})
brier.scores = sapply(1:num.chunks, function(k){
  measures::multiclass.Brier(t(sense.probs.post.means[[k]][,idx]), as.numeric(factor(snippets.info$sense.id[idx])))})

#Plot Brier scores
dev.off()
par(oma = c(0,0,0,0), mar=c(2.75,2.75,1.5,0.75), cex.lab = 0.8, cex.axis = 0.8, cex.main = 0.8, mgp = c(1.5,0.5,0))
plot(brier.scores, xlab = "MCMC chunk", ylab = "Brier score")

#Trace plots - phi
dev.off()
par(mfrow = c(1,3), oma = c(0,0,0,0), mar=c(2.75,1.75,1.5,0.75), cex.lab = 1.2, cex.axis = 1.2, cex.main = 1.2, mgp = c(1.5,0.5,0))
Time = 8; Genre = 2
plot(phi.tilde.sim[1:num.samples, desired.labelling[1], Genre, Time], type = "l", 
     xlab = "Iterations", ylab = "", ylim = c(0,1), main = "Sense 1")
plot(phi.tilde.sim[1:num.samples, desired.labelling[2], Genre, Time], type = "l", 
     xlab = "Iterations", ylab = "", ylim = c(0,1), main = "Sense 2")
plot(phi.tilde.sim[1:num.samples, desired.labelling[3], Genre, Time], type = "l", 
     xlab = "Iterations", ylab = "", ylim = c(0,1), main = "Sense 3")

#Width = 12, Height = 3



#############################################################################
# ESS per unit time ------------------------------------------------------- #
#############################################################################

#phi.tilde
phi.tilde.ess = sapply(1:num.periods, function(t) {
  sapply(1:num.genres, function(g) effectiveSize(as.mcmc(phi.tilde.sim[(burn.in+1):num.samples,,g,t])))})/
  as.double(Run.Time * (1 - burn.in/num.samples), units = "hours")
median(phi.tilde.ess)
quantile(phi.tilde.ess, probs = c(0.25,0.75))

#psi.tilde
psi.tilde.post.mean = sapply(1:num.periods, function(t) {
  sapply(1:num.senses, function(k) {apply(psi.tilde.sim[(burn.in+1):num.samples,,k,t], 2, mean)})}, simplify = "array" )
dimnames(psi.tilde.post.mean) = list(Word = 1:num.words, Sense = 1:num.senses, Time = 1:num.periods)

m = 20
psi.tilde.post.mean.marg.over.time = apply(psi.tilde.post.mean, 1:2, sum) / num.periods
top.m.words.marg.over.time = sapply(1:num.senses, function(k) {order(psi.tilde.post.mean.marg.over.time[,k], decreasing = TRUE)[1:m]})
dimnames(top.m.words.marg.over.time) = list(Rank = 1:m, Sense = 1:num.senses)

Words1 = top.m.words.marg.over.time[,1]
Words2 = top.m.words.marg.over.time[,2]

psi.tilde.ess = array(c(
  sapply(1:num.periods, function(t) {
    effectiveSize(as.mcmc(psi.tilde.sim[(burn.in+1):num.samples, Words1, 1, t])) / 
      as.double(Run.Time * (1 - burn.in/num.samples), units = "hours")
  }),
  sapply(1:num.periods, function(t) {
    effectiveSize(as.mcmc(psi.tilde.sim[(burn.in+1):num.samples, Words2, 2, t])) / 
      as.double(Run.Time * (1 - burn.in/num.samples), units = "hours")
  })
), dim = c(m, num.periods, 2))

median(psi.tilde.ess)
quantile(psi.tilde.ess, probs = c(0.25,0.75))

