#############################################################################
# Scalability analysis ---------------------------------------------------- #
#############################################################################

library(tidyr)
library(plyr)

Ms = c(25, 200)
Ns = c(500)
Vs = c(250, 500, 1000, 2000, 3000)
Ds = c(250, 500, 1000, 2500, 5000, 7500, 10000)
runs = c(1,2,3)

#Read EDiSC run times
run.times.EDiSC = data.frame(embedding.dim = numeric(), sims = numeric(), num.words = numeric(), 
                             num.snippets = numeric(), run = numeric(), run.time = numeric())
for (M in Ms) {
  for (N in Ns) {
    for (V in Vs) {
      for (D in Ds) {
        for (run in runs) {
          if (file.exists(paste("run.time_EDiSC_M=",M,"_sims=",N,"_V=",V,"_D=",D,"_run",run,".RData", sep=""))) {
            load(paste("run.time_EDiSC_M=",M,"_sims=",N,"_V=",V,"_D=",D,"_run",run,".RData", sep=""))
            run.times.EDiSC[nrow(run.times.EDiSC)+1,] = c(M,N,V,D,run, as.double(Run.Time, units = "secs"))
          }#if
        }#for run
      }#for D
    }#for V
  }#for N
}#for M

#Read DiSC run times
run.times.DiSC = data.frame(sims = numeric(), num.words = numeric(), num.snippets = numeric(), 
                            run = numeric(), run.time = numeric())
for (N in Ns) {
  for (V in Vs) {
    for (D in Ds) {
      for (run in runs) {
        if (file.exists(paste("run.time_DiSC_sims=",N,"_V=",V,"_D=",D,"_run",run,".RData", sep=""))) {
          load(paste("run.time_DiSC_sims=",N,"_V=",V,"_D=",D,"_run",run,".RData", sep=""))
          run.times.DiSC[nrow(run.times.DiSC)+1,] = c(N,V,D,run, as.double(Run.Time, units = "secs"))
        }#if
      }#for run
    }#for D
  }#for V
}#for N

#Re-arrange data
run.times.EDiSC = pivot_wider(run.times.EDiSC, names_from = run, values_from = run.time, names_prefix = "run")
run.times.DiSC = pivot_wider(run.times.DiSC, names_from = run, values_from = run.time, names_prefix = "run")

run.times.EDiSC$mean.run.time = rowMeans(run.times.EDiSC[,5:7])
run.times.DiSC$mean.run.time = rowMeans(run.times.DiSC[,4:6])


#Run time vs number of snippets, for a range of V and M
par(mfrow = c(1,3), oma = c(2,0,0,0), mar=c(2.5,2.5,2,1), cex.lab = 1, cex.axis = 1, cex.main = 1, mgp = c(1.5,0.5,0))

plot(NULL, xlim = c(0,max(Ds)), ylim = c(0,max(run.times.EDiSC$mean.run.time, run.times.DiSC$mean.run.time)), 
     main = "DiSC", xlab = "snippets", ylab = "run time (secs)")
for (i in 1:length(Vs)) {
  lines(Ds, run.times.DiSC$mean.run.time[run.times.DiSC$num.words==Vs[i]], type = "o", pch = 18, col = i)
}#for i

M = 25 
plot(NULL, xlim = c(0,max(Ds)), ylim = c(0,max(run.times.EDiSC$mean.run.time, run.times.DiSC$mean.run.time)), 
     main = paste0("EDiSC (M = ",M,")"), xlab = "snippets", ylab = "run time (secs)")
for (i in 1:length(Vs)) {
  lines(Ds, run.times.EDiSC$mean.run.time[run.times.EDiSC$embedding.dim==M & run.times.EDiSC$num.words==Vs[i]],
        type = "o", pch = 18, col = i,  lty = 1)  
}#for i

M = 200
plot(NULL, xlim = c(0,max(Ds)), ylim = c(0,max(run.times.EDiSC$mean.run.time, run.times.DiSC$mean.run.time)), 
     main = paste0("EDiSC (M = ",M,")"), xlab = "snippets", ylab = "run time (secs)")
for (i in 1:length(Vs)) {
  lines(Ds, run.times.EDiSC$mean.run.time[run.times.EDiSC$embedding.dim==M & run.times.EDiSC$num.words==Vs[i]],
        type = "o", pch = 18, col = i,  lty = 1)  
}#for i

par(mfrow = c(1,1), oma = c(0,0,0,0), mar = c(0.5,0,0,0), new = TRUE)
plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE) #overlay emply plot
legend("bottom", horiz = TRUE, bty = "n", legend = paste0("V = ", rev(Vs)), 
       col = length(Vs):1, pch = 18, lty = c(1,1,1), cex = 0.7, x.intersp = 0.7)

#Width = 8, Height = 3


#Regression parameters
M = 25
EDiSC.run.times = lm(mean.run.time ~ num.words * num.snippets, data = run.times.EDiSC[run.times.EDiSC$embedding.dim==M,])
summary(EDiSC.run.times)

M = 200
EDiSC.run.times = lm(mean.run.time ~ num.words * num.snippets, data = run.times.EDiSC[run.times.EDiSC$embedding.dim==M,])
summary(EDiSC.run.times)

DiSC.run.times = lm(mean.run.time ~ num.words * num.snippets, data = run.times.DiSC)
summary(DiSC.run.times)

