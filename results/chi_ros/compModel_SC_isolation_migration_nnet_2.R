#!/usr/bin/env Rscript

library(nnet)
source("cv4abc.R")
nsim = 10 * 100000

ss = c(2, 4, 6, 8, 10, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38)

target = read.table("/scratch/ul/monthly/croux/ABCh/simulations/largePrior_v4/ABCmonolocus/analysis/summary_stats_lengthScaled_2.txt", h=T)[, ss]

isolation = matrix(NA, ncol = length(ss), nrow = nsim)
migration = matrix(NA, ncol = length(ss), nrow = nsim)

path = "/scratch/ul/monthly/croux/ABCh/simulations/largePrior_v4/ABCmonolocus"
for(i in 1:10){
	# isolation
	tmp = matrix(scan(paste(path, "/isolation/rep_", i-1, "/ABCstat.txt", sep=""), skip=2), byrow=T, ncol=42)[, ss] 
	isolation[(i*100000-99999):((i)*100000), ] = tmp

	# migration
	tmp = matrix(scan(paste(path, "/migration/rep_", i-1, "/ABCstat.txt", sep=""), skip=2), byrow=T, ncol=42)[, ss] 
	migration[(i*100000-99999):((i)*100000), ] = tmp
}

for(i in 1:length(ss)){
	isolation[which(isolation[,i]=="NaN"),i]=mean(isolation[,i], na.rm=T)
	migration[which(migration[,i]=="NaN"),i]=mean(migration[,i], na.rm=T)
	
}
# SI AM IM SC
setwd(path)
x=c(rep(1:2, each=nrow(isolation)))
res_isolation_migration = model_selection_abc_nnet(target=target, x=x, sumstat=rbind(isolation, migration), tol=1000/(2*nrow(isolation)), noweight=F, rejmethod=F, nb.nnet=10, size.nnet=8, output="OBS_isolation_migration_2.txt")


