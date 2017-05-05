# AM_Mhetero_Nhetero; AM_Mhetero_Nhomo; AM_Mhomo_Nhetero; AM_Mhomo_Nhomo
# IM_Mhetero_Nhetero; IM_Mhetero_Nhomo; IM_Mhomo_Nhetero; IM_Mhomo_Nhomo
# SC_Mhetero_Nhetero; SC_Mhetero_Nhomo; SC_Mhomo_Nhetero; SC_Mhomo_Nhomo
# SI_Nhetero; SI_Nhomo

library(abcrf)
#pair = "chi_ros"
pair = "ros_txn"
mypath = paste("/home/croux/Documents/heliconius/inputABC/", pair, "/simulations/", sep="")
ss = -1

# Data
# AM
AMHH = read.table(paste(mypath, "AM_Mhetero_Nhetero/ABCstat.txt", sep=""), h=T)[,ss]
nStats = ncol(AMHH)
AMHh = read.table(paste(mypath, "AM_Mhetero_Nhomo/ABCstat.txt", sep=""), h=T)[,ss]
AMhH = read.table(paste(mypath, "AM_Mhomo_Nhetero/ABCstat.txt", sep=""), h=T)[,ss]
AMhh = read.table(paste(mypath, "AM_Mhomo_Nhomo/ABCstat.txt", sep=""), h=T)[,ss]
AMHH_par = read.table(paste(mypath, "AM_Mhetero_Nhetero/priorfile", sep=""), h=T)
AMHh_par = read.table(paste(mypath, "AM_Mhetero_Nhomo/priorfile", sep=""), h=T)
AMhH_par = read.table(paste(mypath, "AM_Mhomo_Nhetero/priorfile", sep=""), h=T)
AMhh_par = read.table(paste(mypath, "AM_Mhomo_Nhomo/priorfile", sep=""), h=T)
AMHH = cbind(AMHH, AMHH_par); AMHh = cbind(AMHh, AMHh_par); AMhH = cbind(AMhH, AMhH_par); AMhh = cbind(AMhh, AMhh_par)

# IM
IMHH = read.table(paste(mypath, "IM_Mhetero_Nhetero/ABCstat.txt", sep=""), h=T)[,ss]
IMHh = read.table(paste(mypath, "IM_Mhetero_Nhomo/ABCstat.txt", sep=""), h=T)[,ss]
IMhH = read.table(paste(mypath, "IM_Mhomo_Nhetero/ABCstat.txt", sep=""), h=T)[,ss]
IMhh = read.table(paste(mypath, "IM_Mhomo_Nhomo/ABCstat.txt", sep=""), h=T)[,ss]
IMHH_par = read.table(paste(mypath, "IM_Mhetero_Nhetero/priorfile", sep=""), h=T)
IMHh_par = read.table(paste(mypath, "IM_Mhetero_Nhomo/priorfile", sep=""), h=T)
IMhH_par = read.table(paste(mypath, "IM_Mhomo_Nhetero/priorfile", sep=""), h=T)
IMhh_par = read.table(paste(mypath, "IM_Mhomo_Nhomo/priorfile", sep=""), h=T)
IMHH = cbind(IMHH, IMHH_par); IMHh = cbind(IMHh, IMHh_par); IMhH = cbind(IMhH, IMhH_par); IMhh = cbind(IMhh, IMhh_par)

# SC
SCHH = read.table(paste(mypath, "SC_Mhetero_Nhetero/ABCstat.txt", sep=""), h=T)[,ss]
SCHh = read.table(paste(mypath, "SC_Mhetero_Nhomo/ABCstat.txt", sep=""), h=T)[,ss]
SChH = read.table(paste(mypath, "SC_Mhomo_Nhetero/ABCstat.txt", sep=""), h=T)[,ss]
SChh = read.table(paste(mypath, "SC_Mhomo_Nhomo/ABCstat.txt", sep=""), h=T)[,ss]
SCHH_par = read.table(paste(mypath, "SC_Mhetero_Nhetero/priorfile", sep=""), h=T)
SCHh_par = read.table(paste(mypath, "SC_Mhetero_Nhomo/priorfile", sep=""), h=T)
SChH_par = read.table(paste(mypath, "SC_Mhomo_Nhetero/priorfile", sep=""), h=T)
SChh_par = read.table(paste(mypath, "SC_Mhomo_Nhomo/priorfile", sep=""), h=T)
SCHH = cbind(SCHH, SCHH_par); SCHh = cbind(SCHh, SCHh_par); SChH = cbind(SChH, SChH_par); SChh = cbind(SChh, SChh_par)

# SI
SIH = read.table(paste(mypath, "SI_Nhetero/ABCstat.txt", sep=""), h=T)[,ss]
SIh = read.table(paste(mypath, "SI_Nhomo/ABCstat.txt", sep=""), h=T)[,ss]
SIH_par = read.table(paste(mypath, "SI_Nhetero/priorfile", sep=""), h=T)
SIh_par = read.table(paste(mypath, "SI_Nhomo/priorfile", sep=""), h=T)
SIH = cbind(SIH, SIH_par); SIh = cbind(SIh, SIh_par)

AMHH = na.omit(AMHH); AMHh = na.omit(AMHh); AMhH = na.omit(AMhH); AMhh = na.omit(AMhh)
IMHH = na.omit(IMHH); IMHh = na.omit(IMHh); IMhH = na.omit(IMhH); IMhh = na.omit(IMhh)
SCHH = na.omit(SCHH); SCHh = na.omit(SCHh); SChH = na.omit(SChH); SChh = na.omit(SChh)
SIH = na.omit(SIH); SIh = na.omit(SIh)

# analyse
observation = read.table(paste("/home/croux/Documents/heliconius/inputABC/", pair, "/observation/target.txt", sep=""), h=T)[ss]

# SI AM IM SC
modIndex = as.factor(c(rep("SIH", nrow(SIH)), rep("AMHH", nrow(AMHH)), rep("IMHH", nrow(IMHH)), rep("SCHH", nrow(SCHH))))
sumsta = rbind(SIH[, 1:nStats], AMHH[, 1:nStats], IMHH[, 1:nStats], SCHH[, 1:nStats])
data1 = data.frame(modIndex, sumsta)

mod_SIAMIMSC = abcrf(modIndex~., data1, ntree = 2000, paral = T)
predict_SIAMIMSC = predict(mod_SIAMIMSC, rbind(observation, observation, observation), data1, ntree = 2000, paral = T)
predict_SIAMIMSC

# best model
#  selected model votes model1 votes model2 votes model3 votes model4 post.proba
#1           SCHH        0.021        0.411       0.4475       0.1205  0.8083417

# "N1"
# "N1past"
# "tauN1past"
# "N2"
# "N2past"
# "tauN2past"
# "Na"
# "Tsplit"
# "Tsc"

# N1
N1 = SCHH$N1
data2 = data.frame(N1, SCHH[, 1:nStats])
mod_N1 = regAbcrf(N1~., data2, ntree=2000, paral = T)
predictN1 = predict(mod_N1, rbind(observation, observation), data2, ntree=2000, paral = T)
predictN1

# N2
N2 = SCHH$N2
data2 = data.frame(N2, SCHH[, 1:nStats])
mod_N2 = regAbcrf(N2~., data2, ntree=2000, paral = T)
predictN2 = predict(mod_N2, rbind(observation, observation), data2, ntree=2000, paral = T)
predictN2

# N1past
N1past = SCHH$N1past
data2 = data.frame(N1past, SCHH[, 1:nStats])
mod_N1past = regAbcrf(N1past~., data2, ntree=2000, paral = T)
predictN1past = predict(mod_N1past, rbind(observation, observation), data2, ntree=2000, paral = T)
predictN1past

# N2past
N2past = SCHH$N2past
data2 = data.frame(N2past, SCHH[, 1:nStats])
mod_N2past = regAbcrf(N2past~., data2, ntree=2000, paral = T)
predictN2past = predict(mod_N2past, rbind(observation, observation), data2, ntree=2000, paral = T)
predictN2past

# tauN1past
tauN1past = SCHH$tauN1past
data2 = data.frame(tauN1past, SCHH[, 1:nStats])
mod_tauN1past = regAbcrf(tauN1past~., data2, ntree=2000, paral = T)
predicttauN1past = predict(mod_tauN1past, rbind(observation, observation), data2, ntree=2000, paral = T)
predicttauN1past

# tauN2past
tauN2past = SCHH$tauN2past
data2 = data.frame(tauN2past, SCHH[, 1:nStats])
mod_tauN2past = regAbcrf(tauN2past~., data2, ntree=2000, paral = T)
predicttauN2past = predict(mod_tauN2past, rbind(observation, observation), data2, ntree=2000, paral = T)
predicttauN2past

# Tsplit
Tsplit = SCHH$Tsplit
data2 = data.frame(Tsplit, SCHH[, 1:nStats])
mod_Tsplit = regAbcrf(Tsplit~., data2, ntree=2000, paral = T)
predictTsplit = predict(mod_Tsplit, rbind(observation, observation), data2, ntree=2000, paral = T)
predictTsplit

# Tsc
Tsc = SCHH$Tsc
data2 = data.frame(Tsc, SCHH[, 1:nStats])
mod_Tsc = regAbcrf(Tsc~., data2, ntree=2000, paral = T)
predictTsc = predict(mod_Tsc, rbind(observation, observation), data2, ntree=2000, paral = T)
predictTsc

