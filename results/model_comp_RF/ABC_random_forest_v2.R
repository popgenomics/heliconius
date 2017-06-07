library(abcrf)
setwd("/scratch/ul/monthly/croux/ABCh/simulations/largePrior_v4")

SI = AMhomo = AMhetero = IMhomo = IMhetero = SChomo = SChetero = migInchetero = migDechetero = NULL

for(i in 1:10){
	SI_stats = read.table(paste("SI_Nhetero_", i, "/ABCstat.txt", sep=""), skip=1, h=T)[, -1]
	SI_param = read.table(paste("SI_Nhetero_", i, "/priorfile", sep=""), skip=1, h=F)[, -1]
	tmp = cbind(SI_stats, SI_param)
	SI = rbind(SI, tmp)

	AM_stats = read.table(paste("AM_Mhomo_Nhetero_", i, "/ABCstat.txt", sep=""), skip=1, h=T)[, -1]
	AM_param = read.table(paste("AM_Mhomo_Nhetero_", i, "/priorfile", sep=""), skip=1, h=F)[, -1]
	tmp = cbind(AM_stats, AM_param)
	AMhomo = rbind(AMhomo, tmp)

	AM_stats = read.table(paste("AM_Mhetero_Nhetero_", i, "/ABCstat.txt", sep=""), skip=1, h=T)[, -1]
	AM_param = read.table(paste("AM_Mhetero_Nhetero_", i, "/priorfile", sep=""), skip=1, h=F)[, -1]
	tmp = cbind(AM_stats, AM_param)
	AMhetero = rbind(AMhetero, tmp)

#	PAM_stats = read.table(paste("PAM_Mhetero_Nhetero_", i, "/ABCstat.txt", sep=""), skip=1, h=T)[, -1]
#	PAM_param = read.table(paste("PAM_Mhetero_Nhetero_", i, "/priorfile", sep=""), skip=1, h=F)[, -1]
#	tmp = cbind(PAM_stats, PAM_param)
#	PAM = rbind(PAM, tmp)

	IM_stats = read.table(paste("IM_Mhomo_Nhetero_", i, "/ABCstat.txt", sep=""), skip=1, h=T)[, -1]
	IM_param = read.table(paste("IM_Mhomo_Nhetero_", i, "/priorfile", sep=""), skip=1, h=F)[, -1]
	tmp = cbind(IM_stats, IM_param)
	IMhomo = rbind(IMhomo, tmp)

	IM_stats = read.table(paste("IM_Mhetero_Nhetero_", i, "/ABCstat.txt", sep=""), skip=1, h=T)[, -1]
	IM_param = read.table(paste("IM_Mhetero_Nhetero_", i, "/priorfile", sep=""), skip=1, h=F)[, -1]
	tmp = cbind(IM_stats, IM_param)
	IMhetero = rbind(IMhetero, tmp)
	
	SC_stats = read.table(paste("SC_Mhomo_Nhetero_", i, "/ABCstat.txt", sep=""), skip=1, h=T)[, -1]
	SC_param = read.table(paste("SC_Mhomo_Nhetero_", i, "/priorfile", sep=""), skip=1, h=F)[, -1]
	tmp = cbind(SC_stats, SC_param)
	SChomo = rbind(SChomo, tmp)

	SC_stats = read.table(paste("SC_Mhetero_Nhetero_", i, "/ABCstat.txt", sep=""), skip=1, h=T)[, -1]
	SC_param = read.table(paste("SC_Mhetero_Nhetero_", i, "/priorfile", sep=""), skip=1, h=F)[, -1]
	tmp = cbind(SC_stats, SC_param)
	SChetero = rbind(SChetero, tmp)
	
	migInc_stats = read.table(paste("migInc_Mhetero_Nhetero_", i, "/ABCstat.txt", sep=""), skip=1, h=T)[, -1]
	migInc_param = read.table(paste("migInc_Mhetero_Nhetero_", i, "/priorfile", sep=""), skip=1, h=F)[, -1]
	tmp = cbind(migInc_stats, migInc_param)
	migInchetero = rbind(migInchetero, tmp)

	migDec_stats = read.table(paste("migDec_Mhetero_Nhetero_", i, "/ABCstat.txt", sep=""), skip=1, h=T)[, -1]
	migDec_param = read.table(paste("migDec_Mhetero_Nhetero_", i, "/priorfile", sep=""), skip=1, h=F)[, -1]
	tmp = cbind(migDec_stats, migDec_param)
	migDechetero = rbind(migDechetero, tmp)

#	PSC_stats = read.table(paste("PSC_Mhetero_Nhetero_", i, "/ABCstat.txt", sep=""), skip=1, h=T)[, -1]
#	PSC_param = read.table(paste("PSC_Mhetero_Nhetero_", i, "/priorfile", sep=""), skip=1, h=F)[, -1]
#	tmp = cbind(PSC_stats, PSC_param)
#	PSC = rbind(PSC, tmp)
}

SI = na.omit(SI)
AMhomo = na.omit(AMhomo)
AMhetero = na.omit(AMhetero)
IMhomo = na.omit(IMhomo)
IMhetero = na.omit(IMhetero)
SChomo = na.omit(SChomo)
SChetero = na.omit(SChetero)
migInchetero = na.omit(migInchetero)
migDechetero = na.omit(migDechetero)

observation = read.table("/scratch/ul/monthly/croux/ABCh/simulations/observed_datasets.txt", h=T)[, -1]

# intra AM
modIndex1 = as.factor(c(rep("AMhomo", nrow(AMhomo)), rep("AMhetero", nrow(AMhetero))))
data1 = data.frame(modIndex1, rbind(AMhomo[, 1:41], AMhetero[, 1:41]))
mod_rf1 = abcrf(modIndex1~., data1, ntree = 10000, paral = T, ncores=10)
mod_rf1
predict_AMhomo_AMhetero = predict(mod_rf1, observation, data1, ntree = 10000, paral = T, ncores=10)

predict_AMhomo_AMhetero

y=read.table("../observed_datasets.txt",h=T)[,1]
res1 = cbind(as.character(y), predict_AMhomo_AMhetero$vote, predict_AMhomo_AMhetero$allocation, predict_AMhomo_AMhetero$post.prob)
colnames(res1) = c("pairs", "AMhomo", "AMhetero", "vote", "post_prob")

write.table(res1, "res_randomForest_AM_homo_hetero.txt", col.names=T, sep="\t", quote=F)

# intra IM
modIndex2 = as.factor(c(rep("IMhomo", nrow(IMhomo)), rep("IMhetero", nrow(IMhetero))))
data2 = data.frame(modIndex2, rbind(IMhomo[, 1:41], IMhetero[, 1:41]))
mod_rf2 = abcrf(modIndex2~., data2, ntree = 10000, paral = T, ncores=10)
mod_rf2
predict_IMhomo_IMhetero = predict(mod_rf2, observation, data2, ntree = 10000, paral = T, ncores=10)

predict_IMhomo_IMhetero

y=read.table("../observed_datasets.txt",h=T)[,1]
res2 = cbind(as.character(y), predict_IMhomo_IMhetero$vote, predict_IMhomo_IMhetero$allocation, predict_IMhomo_IMhetero$post.prob)
colnames(res2) = c("pairs", "IMhomo", "IMhetero", "vote", "post_prob")

write.table(res2, "res_randomForest_IM_homo_hetero.txt", col.names=T, sep="\t", quote=F)

# intra SC
modIndex3 = as.factor(c(rep("SChomo", nrow(SChomo)), rep("SChetero", nrow(SChetero))))
data3 = data.frame(modIndex3, rbind(SChomo[, 1:41], SChetero[, 1:41]))
mod_rf3 = abcrf(modIndex3~., data3, ntree = 10000, paral = T, ncores=10)
mod_rf3
predict_SChomo_SChetero = predict(mod_rf3, observation, data3, ntree = 10000, paral = T, ncores=10)

predict_SChomo_SChetero

y=read.table("../observed_datasets.txt",h=T)[,1]
res3 = cbind(as.character(y), predict_SChomo_SChetero$vote, predict_SChomo_SChetero$allocation, predict_SChomo_SChetero$post.prob)
colnames(res3) = c("pairs", "SChomo", "SChetero", "vote", "post_prob")

write.table(res3, "res_randomForest_SC_homo_hetero.txt", col.names=T, sep="\t", quote=F)

# SI vs AM, IM, SC, migInc, migDec
# SI vs AM
modIndex4 = as.factor(c(rep("SI", nrow(SI)), rep("AMhetero", nrow(AMhetero))))
data4 = data.frame(modIndex4, rbind(SI[, 1:41], AMhetero[, 1:41]))
mod_rf4 = abcrf(modIndex4~., data4, ntree = 10000, paral = T, ncores=10)
mod_rf4
predict_SI_AMhetero = predict(mod_rf4, observation, data4, ntree = 10000, paral = T, ncores=10)

predict_SI_AMhetero

y=read.table("../observed_datasets.txt",h=T)[,1]
res4 = cbind(as.character(y), predict_SI_AMhetero$vote, predict_SI_AMhetero$allocation, predict_SI_AMhetero$post.prob)
colnames(res4) = c("pairs", "SI", "AMhetero", "vote", "post_prob")

write.table(res4, "res_randomForest_SI_AMhetero.txt", col.names=T, sep="\t", quote=F)

# SI vs IM
modIndex5 = as.factor(c(rep("SI", nrow(SI)), rep("IMhetero", nrow(IMhetero))))
data5 = data.frame(modIndex5, rbind(SI[, 1:41], IMhetero[, 1:41]))
mod_rf5 = abcrf(modIndex5~., data5, ntree = 10000, paral = T, ncores=10)
mod_rf5
predict_SI_IMhetero = predict(mod_rf5, observation, data5, ntree = 10000, paral = T, ncores=10)

predict_SI_IMhetero

y=read.table("../observed_datasets.txt",h=T)[,1]
res5 = cbind(as.character(y), predict_SI_IMhetero$vote, predict_SI_IMhetero$allocation, predict_SI_IMhetero$post.prob)
colnames(res5) = c("pairs", "SI", "IMhetero", "vote", "post_prob")

write.table(res5, "res_randomForest_SI_IMhetero.txt", col.names=T, sep="\t", quote=F)

# SI vs SC
modIndex6 = as.factor(c(rep("SI", nrow(SI)), rep("SChetero", nrow(SChetero))))
data6 = data.frame(modIndex6, rbind(SI[, 1:41], SChetero[, 1:41]))
mod_rf6 = abcrf(modIndex6~., data6, ntree = 10000, paral = T, ncores=10)
mod_rf6
predict_SI_SChetero = predict(mod_rf6, observation, data6, ntree = 10000, paral = T, ncores=10)

predict_SI_SChetero

y=read.table("../observed_datasets.txt",h=T)[,1]
res6 = cbind(as.character(y), predict_SI_SChetero$vote, predict_SI_SChetero$allocation, predict_SI_SChetero$post.prob)
colnames(res6) = c("pairs", "SI", "SChetero", "vote", "post_prob")

write.table(res6, "res_randomForest_SI_SChetero.txt", col.names=T, sep="\t", quote=F)

# SI vs migInc
modIndex7 = as.factor(c(rep("SI", nrow(SI)), rep("migInchetero", nrow(migInchetero))))
data7 = data.frame(modIndex7, rbind(SI[, 1:41], migInchetero[, 1:41]))
mod_rf7 = abcrf(modIndex7~., data7, ntree = 10000, paral = T, ncores=10)
mod_rf7
predict_SI_migInchetero = predict(mod_rf7, observation, data7, ntree = 10000, paral = T, ncores=10)

predict_SI_migInchetero

y=read.table("../observed_datasets.txt",h=T)[,1]
res7 = cbind(as.character(y), predict_SI_migInchetero$vote, predict_SI_migInchetero$allocation, predict_SI_migInchetero$post.prob)
colnames(res7) = c("pairs", "SI", "migInchetero", "vote", "post_prob")

write.table(res7, "res_randomForest_SI_migInchetero.txt", col.names=T, sep="\t", quote=F)

# SI vs migDec
modIndex8 = as.factor(c(rep("SI", nrow(SI)), rep("migDechetero", nrow(migDechetero))))
data8 = data.frame(modIndex8, rbind(SI[, 1:41], migDechetero[, 1:41]))
mod_rf8 = abcrf(modIndex8~., data8, ntree = 10000, paral = T, ncores=10)
mod_rf8
predict_SI_migDechetero = predict(mod_rf8, observation, data8, ntree = 10000, paral = T, ncores=10)

predict_SI_migDechetero

y=read.table("../observed_datasets.txt",h=T)[,1]
res8 = cbind(as.character(y), predict_SI_migDechetero$vote, predict_SI_migDechetero$allocation, predict_SI_migDechetero$post.prob)
colnames(res8) = c("pairs", "SI", "migDechetero", "vote", "post_prob")

write.table(res8, "res_randomForest_SI_migDechetero.txt", col.names=T, sep="\t", quote=F)


# AM vs IM, SC, migInc, migDec
# AM vs IM
modIndex9 = as.factor(c(rep("AMhetero", nrow(AMhetero)), rep("IMhetero", nrow(IMhetero))))
data9 = data.frame(modIndex9, rbind(AMhetero[, 1:41], IMhetero[, 1:41]))
mod_rf9 = abcrf(modIndex9~., data9, ntree = 10000, paral = T, ncores=10)
mod_rf9
predict_AMhetero_IMhetero = predict(mod_rf9, observation, data9, ntree = 10000, paral = T, ncores=10)

predict_AMhetero_IMhetero

y=read.table("../observed_datasets.txt",h=T)[,1]
res9 = cbind(as.character(y), predict_AMhetero_IMhetero$vote, predict_AMhetero_IMhetero$allocation, predict_AMhetero_IMhetero$post.prob)
colnames(res9) = c("pairs", "AMhetero", "IMhetero", "vote", "post_prob")

write.table(res9, "res_randomForest_AMhetero_IMhetero.txt", col.names=T, sep="\t", quote=F)

# AM vs SC
modIndex10 = as.factor(c(rep("AMhetero", nrow(AMhetero)), rep("SChetero", nrow(SChetero))))
data10 = data.frame(modIndex10, rbind(AMhetero[, 1:41], SChetero[, 1:41]))
mod_rf10 = abcrf(modIndex10~., data10, ntree = 10000, paral = T, ncores=10)
mod_rf10
predict_AMhetero_SChetero = predict(mod_rf10, observation, data10, ntree = 10000, paral = T, ncores=10)

predict_AMhetero_SChetero

y=read.table("../observed_datasets.txt",h=T)[,1]
res10 = cbind(as.character(y), predict_AMhetero_SChetero$vote, predict_AMhetero_SChetero$allocation, predict_AMhetero_SChetero$post.prob)
colnames(res10) = c("pairs", "AMhetero", "SChetero", "vote", "post_prob")

write.table(res10, "res_randomForest_AMhetero_SChetero.txt", col.names=T, sep="\t", quote=F)

# AM vs migInc
modIndex11 = as.factor(c(rep("AMhetero", nrow(AMhetero)), rep("migInchetero", nrow(migInchetero))))
data11 = data.frame(modIndex11, rbind(AMhetero[, 1:41], migInchetero[, 1:41]))
mod_rf11 = abcrf(modIndex11~., data11, ntree = 10000, paral = T, ncores=10)
mod_rf11
predict_AMhetero_migInchetero = predict(mod_rf11, observation, data11, ntree = 10000, paral = T, ncores=10)

predict_AMhetero_migInchetero

y=read.table("../observed_datasets.txt",h=T)[,1]
res11 = cbind(as.character(y), predict_AMhetero_migInchetero$vote, predict_AMhetero_migInchetero$allocation, predict_AMhetero_migInchetero$post.prob)
colnames(res11) = c("pairs", "AMhetero", "migInchetero", "vote", "post_prob")

write.table(res11, "res_randomForest_AMhetero_migInchetero.txt", col.names=T, sep="\t", quote=F)

# AM vs migDec
modIndex12 = as.factor(c(rep("AMhetero", nrow(AMhetero)), rep("migDechetero", nrow(migDechetero))))
data12 = data.frame(modIndex12, rbind(AMhetero[, 1:41], migDechetero[, 1:41]))
mod_rf12 = abcrf(modIndex12~., data12, ntree = 10000, paral = T, ncores=10)
mod_rf12
predict_AMhetero_migDechetero = predict(mod_rf12, observation, data12, ntree = 10000, paral = T, ncores=10)

predict_AMhetero_migDechetero

y=read.table("../observed_datasets.txt",h=T)[,1]
res12 = cbind(as.character(y), predict_AMhetero_migDechetero$vote, predict_AMhetero_migDechetero$allocation, predict_AMhetero_migDechetero$post.prob)
colnames(res12) = c("pairs", "AMhetero", "migDechetero", "vote", "post_prob")

write.table(res12, "res_randomForest_AMhetero_migDechetero.txt", col.names=T, sep="\t", quote=F)

# IM vs SC, migInc, migDec
# IM vs SC
modIndex13 = as.factor(c(rep("IMhetero", nrow(IMhetero)), rep("SChetero", nrow(SChetero))))
data13 = data.frame(modIndex13, rbind(IMhetero[, 1:41], SChetero[, 1:41]))
mod_rf13 = abcrf(modIndex13~., data13, ntree = 10000, paral = T, ncores=10)
mod_rf13
predict_IMhetero_SChetero = predict(mod_rf13, observation, data13, ntree = 10000, paral = T, ncores=10)

predict_IMhetero_SChetero

y=read.table("../observed_datasets.txt",h=T)[,1]
res13 = cbind(as.character(y), predict_IMhetero_SChetero$vote, predict_IMhetero_SChetero$allocation, predict_IMhetero_SChetero$post.prob)
colnames(res13) = c("pairs", "IMhetero", "SChetero", "vote", "post_prob")

write.table(res13, "res_randomForest_IMhetero_SChetero.txt", col.names=T, sep="\t", quote=F)

# IM vs migInc 
modIndex14 = as.factor(c(rep("IMhetero", nrow(IMhetero)), rep("migInchetero", nrow(migInchetero))))
data14 = data.frame(modIndex14, rbind(IMhetero[, 1:41], migInchetero[, 1:41]))
mod_rf14 = abcrf(modIndex14~., data14, ntree = 10000, paral = T, ncores=10)
mod_rf14
predict_IMhetero_migInchetero = predict(mod_rf14, observation, data14, ntree = 10000, paral = T, ncores=10)

predict_IMhetero_migInchetero

y=read.table("../observed_datasets.txt",h=T)[,1]
res14 = cbind(as.character(y), predict_IMhetero_migInchetero$vote, predict_IMhetero_migInchetero$allocation, predict_IMhetero_migInchetero$post.prob)
colnames(res14) = c("pairs", "IMhetero", "migInchetero", "vote", "post_prob")

write.table(res14, "res_randomForest_IMhetero_migInchetero.txt", col.names=T, sep="\t", quote=F)

# IM vs migDec 
modIndex15 = as.factor(c(rep("IMhetero", nrow(IMhetero)), rep("migDechetero", nrow(migDechetero))))
data15 = data.frame(modIndex15, rbind(IMhetero[, 1:41], migDechetero[, 1:41]))
mod_rf15 = abcrf(modIndex15~., data15, ntree = 10000, paral = T, ncores=10)
mod_rf15
predict_IMhetero_migDechetero = predict(mod_rf15, observation, data15, ntree = 10000, paral = T, ncores=10)

predict_IMhetero_migDechetero

y=read.table("../observed_datasets.txt",h=T)[,1]
res15 = cbind(as.character(y), predict_IMhetero_migDechetero$vote, predict_IMhetero_migDechetero$allocation, predict_IMhetero_migDechetero$post.prob)
colnames(res15) = c("pairs", "IMhetero", "migDechetero", "vote", "post_prob")

write.table(res15, "res_randomForest_IMhetero_migDechetero.txt", col.names=T, sep="\t", quote=F)


# SC vs migInc, migDec
# SC vs migInc 
modIndex16 = as.factor(c(rep("SChetero", nrow(SChetero)), rep("migInchetero", nrow(migInchetero))))
data16 = data.frame(modIndex16, rbind(SChetero[, 1:41], migInchetero[, 1:41]))
mod_rf16 = abcrf(modIndex16~., data16, ntree = 10000, paral = T, ncores=10)
mod_rf16
predict_SChetero_migInchetero = predict(mod_rf16, observation, data16, ntree = 10000, paral = T, ncores=10)

predict_SChetero_migInchetero

y=read.table("../observed_datasets.txt",h=T)[,1]
res16 = cbind(as.character(y), predict_SChetero_migInchetero$vote, predict_SChetero_migInchetero$allocation, predict_SChetero_migInchetero$post.prob)
colnames(res16) = c("pairs", "SChetero", "migInchetero", "vote", "post_prob")

write.table(res16, "res_randomForest_SChetero_migInchetero.txt", col.names=T, sep="\t", quote=F)

# SC vs migDec 
modIndex17 = as.factor(c(rep("SChetero", nrow(SChetero)), rep("migDechetero", nrow(migDechetero))))
data17 = data.frame(modIndex17, rbind(SChetero[, 1:41], migDechetero[, 1:41]))
mod_rf17 = abcrf(modIndex17~., data17, ntree = 10000, paral = T, ncores=10)
mod_rf17
predict_SChetero_migDechetero = predict(mod_rf17, observation, data17, ntree = 10000, paral = T, ncores=10)

predict_SChetero_migDechetero

y=read.table("../observed_datasets.txt",h=T)[,1]
res17 = cbind(as.character(y), predict_SChetero_migDechetero$vote, predict_SChetero_migDechetero$allocation, predict_SChetero_migDechetero$post.prob)
colnames(res17) = c("pairs", "SChetero", "migDechetero", "vote", "post_prob")

write.table(res17, "res_randomForest_SChetero_migDechetero.txt", col.names=T, sep="\t", quote=F)


# migInc vs migDec 
modIndex18 = as.factor(c(rep("migInchetero", nrow(migInchetero)), rep("migDechetero", nrow(migDechetero))))
data18 = data.frame(modIndex18, rbind(migInchetero[, 1:41], migDechetero[, 1:41]))
mod_rf18 = abcrf(modIndex18~., data18, ntree = 10000, paral = T, ncores=10)
mod_rf18
predict_migInchetero_migDechetero = predict(mod_rf18, observation, data18, ntree = 10000, paral = T, ncores=10)

predict_migInchetero_migDechetero

y=read.table("../observed_datasets.txt",h=T)[,1]
res18 = cbind(as.character(y), predict_migInchetero_migDechetero$vote, predict_migInchetero_migDechetero$allocation, predict_migInchetero_migDechetero$post.prob)
colnames(res18) = c("pairs", "migInchetero", "migDechetero", "vote", "post_prob")

write.table(res18, "res_randomForest_migInchetero_migDechetero.txt", col.names=T, sep="\t", quote=F)


#
output = NULL
x = predict_AMhomo_AMhetero
	allocation = as.character(x$allocation)
	vote = x$vote
	prob = x$post.prob
	output = cbind(output, allocation, vote, prob)

x = predict_IMhomo_IMhetero
	allocation = as.character(x$allocation)
	vote = x$vote
	prob = x$post.prob
	output = cbind(output, allocation, vote, prob)

x = predict_SChomo_SChetero
	allocation = as.character(x$allocation)
	vote = x$vote
	prob = x$post.prob
	output = cbind(output, allocation, vote, prob)

x = predict_SI_AMhetero
	allocation = as.character(x$allocation)
	vote = x$vote
	prob = x$post.prob
	output = cbind(output, allocation, vote, prob)

x = predict_SI_IMhetero
	allocation = as.character(x$allocation)
	vote = x$vote
	prob = x$post.prob
	output = cbind(output, allocation, vote, prob)

x = predict_SI_SChetero
	allocation = as.character(x$allocation)
	vote = x$vote
	prob = x$post.prob
	output = cbind(output, allocation, vote, prob)

x = predict_SI_migInchetero
	allocation = as.character(x$allocation)
	vote = x$vote
	prob = x$post.prob
	output = cbind(output, allocation, vote, prob)

x = predict_SI_migDechetero
	allocation = as.character(x$allocation)
	vote = x$vote
	prob = x$post.prob
	output = cbind(output, allocation, vote, prob)

x = predict_AMhetero_IMhetero
	allocation = as.character(x$allocation)
	vote = x$vote
	prob = x$post.prob
	output = cbind(output, allocation, vote, prob)

x = predict_AMhetero_SChetero
	allocation = as.character(x$allocation)
	vote = x$vote
	prob = x$post.prob
	output = cbind(output, allocation, vote, prob)

x = predict_AMhetero_migInchetero
	allocation = as.character(x$allocation)
	vote = x$vote
	prob = x$post.prob
	output = cbind(output, allocation, vote, prob)

x = predict_AMhetero_migDechetero
	allocation = as.character(x$allocation)
	vote = x$vote
	prob = x$post.prob
	output = cbind(output, allocation, vote, prob)

x = predict_IMhetero_SChetero
	allocation = as.character(x$allocation)
	vote = x$vote
	prob = x$post.prob
	output = cbind(output, allocation, vote, prob)

x = predict_IMhetero_migInchetero
	allocation = as.character(x$allocation)
	vote = x$vote
	prob = x$post.prob
	output = cbind(output, allocation, vote, prob)

x = predict_IMhetero_migDechetero
	allocation = as.character(x$allocation)
	vote = x$vote
	prob = x$post.prob
	output = cbind(output, allocation, vote, prob)

x = predict_SChetero_migInchetero
	allocation = as.character(x$allocation)
	vote = x$vote
	prob = x$post.prob
	output = cbind(output, allocation, vote, prob)

x = predict_SChetero_migDechetero
	allocation = as.character(x$allocation)
	vote = x$vote
	prob = x$post.prob
	output = cbind(output, allocation, vote, prob)

x = predict_migInchetero_migDechetero
	allocation = as.character(x$allocation)
	vote = x$vote
	prob = x$post.prob
	output = cbind(output, allocation, vote, prob)

# confusion matrix:
for(i in 1:18){
	a = get(paste("mod_rf",i, sep=""))
	a = a$model.rf$confusion.matrix
	write.table(a, paste("confusion_matrix_", i, sep=""), col.names=T, row.names=T)
}


