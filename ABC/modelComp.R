#!/usr/bin/env Rscript
for(arg in commandArgs()){
        arg=strsplit(arg, "=")
        if(arg[[1]][1]=="input"){species=arg[[1]][2]}
}

library(nnet)
setwd("/Home/croux/work/tools/nnetABC")
source("cv4abc.R")
setwd(paste("/scratch/cluster/monthly/croux/metaanalyse/", species, sep=""))
target=as.numeric(read.table("target.txt",skip=2,h=F))


