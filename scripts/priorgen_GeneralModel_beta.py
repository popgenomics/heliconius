#!/software/bin/python2.7
# -*- coding: utf-8 -*-
from numpy.random import uniform
from numpy.random import beta
from numpy.random import binomial
from random import shuffle
import sys
n1, n1past, n2, n2past, nA, tau, M1, M2, shape1, shape2 = [], [], [], [], [], [], [], [], [], []
help="\n\t\033[32mExternal required library: numpy \033[1;m(sudo apt-get install python-numpy)\n\t\
priorgen.py generates prior distributions for multiple multilocus simulations under 14 different models of speciation. The output can be used from the stdout by ms (Hudson 2002), msnsam (Ross-Ibarra 2008) and msms (Ewing and Hermisson 2010) using the 'tbs' feature.\n\t\
It requires one input file containing six lines: \n\t\t\
L1=description line, non read by priorgen.py\n\t\t\
L2=a vector with the lengths (L) for each of the surveyed locus\n\t\t\
L3=a vector with the number of sampled individuals (nspA) for each locus, for the first population\n\t\t\
L4=a vector with the number of sampled individuals (nspB) for each locus, for the second population\n\t\t\
L5=a vector with the populational mutation rates theta(i)=4.N.µ.L(i) for each locus 'i'\n\t\t\
L6=a vector with the populational recombination rates rho(i)=4.N.r.L(i) for each locus 'i'\n\n\t\
Values print in the stdout are used by ms-like coalescent simulators, values written in a file are the multilocus parameters useful for an ABC analysis\n\n\t\t\
parameters: name of the output file name. Ex \033[1;35mparameters=listOfParameters.txt\033[1;m\n\t\t\
n1: prior for N1 (the effective population size of the first population). Ex \033[1;35mn1=0 n1=10\033[1;m\n\t\t\
n1past: prior for a scalar of N1 (>1 = recent contraction, <1 = recent expansion). The effective population size of the first population between 'tau' and 'tauN1past' = N1 * n1past). Ex \033[1;35mn1past=0 n1past=2\033[1;m\n\n\t\t\
n2: prior for N2 (the effective population size of the second population). Ex\033[1;35m n2=0 n2=10\033[1;m\n\t\t\
n2past: prior for a scalar of N2 (>1 = recent contraction, <1 = recent expansion). The effective population size of the first population between 'tau' and 'tauN1past' = N2 * n1past). Ex \033[1;35mn2past=0 n2past=2\033[1;m\n\n\t\t\
nA: prior for NA (the effective population size of the ancetral population). Ex\033[1;35m nA=0 nA=10\033[1;m\n\t\t\
tau: prior for Tsplit (the time of speciation). Ex\033[1;35m tau=0 tau=3\033[1;m\n\t\t\
M1 (M2): prior for migration rate 4.N1.m1 (4.N2.m2) into the first (second) population. Ex\033[1;35m M1=0 M1=4 M2=0 M2=4\033[1;m\n\t\t\
shape1: prior for the first shape parameter of the Beta distribution. Ex\033[1;35m shape1=0 shape1=10\033[1;m\n\t\t\
shape2: prior for the second shape parameter of the Beta distribution. Ex\033[1;35m shape2=0 shape2=50\033[1;m\n\t\t\
model: migInc and migDec models always assume ancestral homogeneous migration, ongoing migration can be heterogenous.\n\t\t\t\033[1;35m=migInc\033[1;m (migration increases with time; past migration at lower rates than ongoing migration)\n\t\t\t\033[1;35m=SC\033[1;m (no migration in the past), \n\t\t\t\033[1;35m=migDec\033[1;m (past migration is greater than ongoing migration); \n\t\t\t\033[1;35m=AM\033[1;m (no ongoing migration), \n\t\t\t\033[1;35m=IM\033[1;m (same migration over time; even for heteroM models); \n\t\t\t\033[1;35m=SI\033[1;m (no migration between the time of split and present time)\n\t\t\
Nvariation: \033[1;35m=homo\033[1;m (shared values of Ne throughout genome for N1, N2 and Nanc) or \033[1;35m=hetero\033[1;m (variation of Ne throughout genome for N1, N2 and Nanc)\n\t\t\
Mvariation: \033[1;35m=homo\033[1;m (shared values of M throughout genome) or \033[1;35m=hetero\033[1;m (variation of M throughout genome)\n\t\t\
nreps: number of multilocus simulations. Ex \033[1;35mnreps=1000\033[1;m\n\t\t\
symMig: \033[1;35m=sym\033[1;m (M1=M2) or \033[1;35m=asym\033[1;m (M1 and M2 are independently chosen)\n\n\t\
Ex:\n\t\
\033[1;32m./priorgen_GeneralModel_beta.py bpfile=bpfile n1=0 n1=1 n1past=0.1 n1past=0.1 n2=1 n2=2 n2past=2 n2past=2 nA=2 nA=3 tau=3 tau=4 M1=4 M1=5 M2=5 M2=6 shape1=0 shape1=10 shape2=0 shape2=100 nreps=4 Nvariation=homo Mvariation=hetero symMig=asym parameters=priorfile.txt model=migDec\033[1;m\n\n\t\
msnsam tbs 20000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -en tbs 1 tbs -n 2 tbs -en tbs 2 tbs -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs\t\n\t\
camille.roux.1983@gmail.com\n\t\
05/06/2017\n"

for i in sys.argv:
	if("help" in i):
		print(help)
		sys.exit(0)


listOfParam = ["bpfile", "parameters", "n1", "n1past", "n2", "n2past", "nA", "tau", "M1", "M2", "shape1", "shape2", "model", "nreps", "Nvariation", "Mvariation", "symMig"]

if len(sys.argv) != 28:
	print(help)
	sys.exit(0)

for i in sys.argv:
	if "=" in i:
		i=i.split("=")
		if i[0] not in listOfParam:
			print('\033[1;35m \n\t\t{0} is not an allowed parameters like {1}\033[1;m\n\n'.format(i[0], ", ".join(listOfParam)))
			sys.exit()
		if(i[0]=="bpfile"):
			bpfile=i[1]
		if(i[0]=="parameters"):
			outputParameters=i[1]
		if(i[0]=="n1"):
			n1.append(float(i[1]))
		if(i[0]=="n1past"):
			n1past.append(float(i[1]))
		if(i[0]=="n2"):
			n2.append(float(i[1]))
		if(i[0]=="n2past"):
			n2past.append(float(i[1]))
		if(i[0]=="nA"):
			nA.append(float(i[1]))
		if(i[0]=="tau"):
			tau.append(float(i[1]))
		if(i[0]=="M1"):
			M1.append(float(i[1]))
		if(i[0]=="M2"):
			M2.append(float(i[1]))
		if(i[0]=="shape1"):
			shape1.append(float(i[1]))
		if(i[0]=="shape2"):
			shape2.append(float(i[1]))
		if(i[0]=="model"):
			model=i[1]
		if(i[0]=="nreps"):
			nreps=int(i[1])
		if(i[0]=="Nvariation"):
			Nvariation=i[1]
		if(i[0]=="Mvariation"):
			Mvariation=i[1]
		if(i[0]=="symMig"):
			sym=i[1]

listOfModels = ["migInc", "SC", "migDec", "AM", "IM", "SI"]
if model not in listOfModels:
	print('\033[1;35m \n\t\tmodel has to be in ["migInc", "SC", "migDec", "AM", "IM", "SI"]\033[1;m\n\n')
	sys.exit(0)


def binomBeta(nlocus, shape1, shape2, scalar):
	neutre=[0]
	hetero=[1]
#	nNeutre=int(uniform(0, nlocus))
	nNeutre=0
	nHetero=nlocus-nNeutre
	status=nNeutre*neutre+nHetero*hetero
	shuffle(status)
	values=[]
	for i in status:
		if i==0:
			values.append(scalar)
		if i==1:
			values.append(scalar*beta(shape1, shape2))
	res={}
	res["values"]=values
	res["nNeutre"]=nNeutre
	return(res)


def heteroBeta(nlocus, shape1, shape2, scalar):
	values=[]
	for i in range(nlocus):
		values.append(scalar*beta(shape1, shape2))
	res={}
	res["values"] = values
	res["max"] = max(values)
	res["min"] = min(values)
	return(res)


infile=open(bpfile, "r")
tmp=infile.readline()	#skip the header
L=infile.readline().strip().replace(" ", "\t")
L=L.split("\t")
nspA=infile.readline().strip().replace(" ", "\t")
nspA=nspA.split("\t")
nspB=infile.readline().strip().replace(" ", "\t")
nspB=nspB.split("\t")
theta=infile.readline().strip().replace(" ", "\t")
theta=theta.split("\t")
rho=infile.readline().strip().replace(" ", "\t")
rho=rho.split("\t")

nlocus=int(len(L))

res = "N1\tN1past\ttauN1past\tN2\tN2past\ttauN2past\tNa\tTsplit\tTsmall"
if sym == "asym":
	if(Nvariation=="homo"):
		if(Mvariation=="homo"):
			res += "\tM1pres\tM1past\tM2pres\tM2past\n"
		if(Mvariation=="hetero"):
			res += "\tM1pres\tM2pres\tshape1M1\tshape2M1\tshape1M2\tshape2M2\tM1past\tM2past\n"
	if(Nvariation=="hetero"):
		if(Mvariation=="homo"):
			res += "\tshape1Ne\tshape2Ne\tM1pres\tM2pres\tM1past\tM2past\n"
		if(Mvariation=="hetero"):
			res += "\tshape1Ne\tshape2Ne\tM1pres\tM2pres\tshape1M1\tshape2M1\tshape1M2\tshape2M2\tM1past\tM2past\n"
if sym == "sym":
	if(Nvariation=="homo"):
		if(Mvariation=="homo"):
			res += "\tM1pres\tM1past\n"
		if(Mvariation=="hetero"):
			res += "\tM1pres\tshape1M1\tshape2M1\tM1past\n"
	if(Nvariation=="hetero"):
		if(Mvariation=="homo"):
			res += "\tshape1Ne\tshape2Ne\tM1pres\tM1past\n"
		if(Mvariation=="hetero"):
			res += "\tshape1Ne\tshape2Ne\tM1pres\tshape1M1\tshape2M1\tM1past\n"




for i in range(nreps):
	n1prior = uniform(n1[0], n1[1]) # N1 = n1prior * N0 (see msdoc.pdf)
	n1pastprior = uniform(n1past[0], n1past[1]) # N1past = n1pastprior * N1
	n2prior = uniform(n2[0], n2[1])
	n2pastprior = uniform(n2past[0], n2past[1])
	nAprior = uniform(nA[0], nA[1])
	Tsplit = uniform(tau[0], tau[1])
	Tsmall = uniform(min(tau), Tsplit)
	tauN1past = uniform(min(tau), Tsplit)
	tauN2past = uniform(min(tau), Tsplit)
	M1prior = uniform(M1[0], M1[1])
	M2prior = uniform(M2[0], M2[1])
	shape1mig1 = uniform(shape1[0], shape1[1]*2)
	shape2mig1 = uniform(shape2[0], shape2[1]*2)
	shape1mig2 = uniform(shape1[0], shape1[1]*2) # tmp: fois 2 pour un prior du shape_mig plus large que pour le shape_Ne
	shape2mig2 = uniform(shape2[0], shape2[1]*2)
	shape1Ne = uniform(shape1[0], shape1[1])
	shape2Ne = uniform(shape2[0], shape2[1])
	TsplitGenomic = [Tsplit]*nlocus
	TsmallGenomic = [Tsmall]*nlocus
	
	if(Nvariation=="hetero"):
		nPriorGenomic = heteroBeta(nlocus=nlocus, shape1=shape1Ne, shape2=shape2Ne, scalar=1)
		n1priorGenomic = nPriorGenomic["values"]*n1prior
		n2priorGenomic = nPriorGenomic["values"]*n2prior
		nApriorGenomic = nPriorGenomic["values"]*nAprior
#		n1priorGenomic=heteroBeta(nlocus=nlocus, shape1=shape1Ne, shape2=shape2Ne, scalar=n1prior)
#		n2priorGenomic=heteroBeta(nlocus=nlocus, shape1=shape1Ne, shape2=shape2Ne, scalar=n2prior)
#		nApriorGenomic=heteroBeta(nlocus=nlocus, shape1=shape1Ne, shape2=shape2Ne, scalar=nAprior)
	
	if(Nvariation=="homo"):
#		n1priorGenomic={}
#		n2priorGenomic={}
#		nApriorGenomic={}
#		n1priorGenomic["values"]=[n1prior]*nlocus
#		n2priorGenomic["values"]=[n2prior]*nlocus
#		nApriorGenomic["values"]=[nAprior]*nlocus
		n1priorGenomic = [n1prior]*nlocus
		n2priorGenomic = [n2prior]*nlocus
		nApriorGenomic = [nAprior]*nlocus
	
	if(Mvariation=="hetero"):
		M1priorGenomic=heteroBeta(nlocus=nlocus, shape1=shape1mig1, shape2=shape2mig1, scalar=M1prior)
		if(sym=="sym"):
			M2priorGenomic={}
#			M2priorGenomic=M1priorGenomic["values"]
			M2priorGenomic=M1priorGenomic.copy()	# Sophie Galina
			M2prior=M1prior
		if(sym=="asym"):
			M2priorGenomic=heteroBeta(nlocus=nlocus, shape1=shape1mig2, shape2=shape2mig2, scalar=M2prior)
	
	if(Mvariation=="homo"):
		M1priorGenomic={}
		M1priorGenomic["values"]=[M1prior]*nlocus
		M1priorGenomic["min"]=M1prior
		M1priorGenomic["max"]=M1prior
		if(sym=="sym"):
			M2priorGenomic={}
#			M2priorGenomic["values"]=M1priorGenomic["values"]
			M2priorGenomic=M1priorGenomic.copy()	# Sophie Galina
			M2prior=M1prior
		if(sym=="asym"):
			M2priorGenomic={}
			M2priorGenomic["values"]=[M2prior]*nlocus

	if(model == "migInc"): # increased migration with time + transition from homoM to heteroM
		M1past = [M1priorGenomic["min"]]*nlocus
		M2past = [M2priorGenomic["min"]]*nlocus
		M1current =  M1priorGenomic["values"]
		M2current =  M2priorGenomic["values"]
		M1past_print = M1priorGenomic["min"]
		M2past_print = M2priorGenomic["min"]
		M1current_print = M1prior
		M2current_print = M2prior
	
	if(model == "SC"): # special case of migInc where past migration == 0
		M1past = [0]*nlocus
		M2past = [0]*nlocus
		M1current = M1priorGenomic["values"]
		M2current = M2priorGenomic["values"]
		M1past_print = 0
		M2past_print = 0
		M1current_print = M1prior
		M2current_print = M2prior
	
	if(model == "migDec"): # decreased migration with time + transition from homoM to heteroM
		M1past = [M1priorGenomic["max"]]*nlocus
		M2past = [M2priorGenomic["max"]]*nlocus
		M1current = M1priorGenomic["values"]
		M2current = M2priorGenomic["values"]
		M1past_print = M1priorGenomic["max"]
		M2past_print = M2priorGenomic["max"]
		M1current_print = M1prior
		M2current_print = M2prior
	
	if(model == "AM"): # special case of migDec where current migration == 0
		M1past = M1priorGenomic["values"]
		M2past = M2priorGenomic["values"]
		M1current = [0]*nlocus
		M2current = [0]*nlocus
		M1past_print = M1prior
		M2past_print = M2prior
		M1current_print = 0
		M2current_print = 0
	
	if(model == "IM"):
		M1past = M1priorGenomic["values"]
		M2past = M2priorGenomic["values"]
		M1current = M1priorGenomic["values"]
		M2current = M2priorGenomic["values"]
		M1past_print = M1prior
		M2past_print = M2prior
		M1current_print = M1prior
		M2current_print = M2prior
	
	if(model == "SI"):
		M1past = [0] * nlocus
		M2past = [0] * nlocus
		M1current = [0] * nlocus
		M2current = [0] * nlocus
		M1past_print = 0
		M2past_print = 0 
		M1current_print = 0
		M2current_print = 0
	
	res+="{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}".format(n1prior, n1pastprior, tauN1past, n2prior, n2pastprior, tauN2past, nAprior, Tsplit, Tsmall)
	
	if sym == "asym":
		if(Nvariation=="homo"):
			if(Mvariation=="homo"):
				res+="\t{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\n".format(M1current_print, M2current_print, M1past_print, M2past_print)
			if(Mvariation=="hetero"):
				res+="\t{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\n".format(M1current_print, M2current_print, shape1mig1, shape2mig1, shape1mig2, shape2mig2, M1past_print, M2past_print)
		if(Nvariation=="hetero"):
			if(Mvariation=="homo"):
				res+="\t{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\n".format(shape1Ne, shape2Ne, M1current_print, M2current_print, M1past_print, M2past_print)
			if(Mvariation=="hetero"):
				res+="\t{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\n".format(shape1Ne, shape2Ne, M1current_print, M2current_print, shape1mig1, shape2mig1, shape1mig2, shape2mig2, M1past_print, M2past_print)
	if sym == "sym":
		if(Nvariation=="homo"):
			if(Mvariation=="homo"):
				res+="\t{0:.5f}\t{1:.5f}\n".format(M1current_print, M1past_print)
			if(Mvariation=="hetero"):
				res+="\t{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\n".format(M1current_print, shape1mig1, shape2mig1, M1past_print)
		if(Nvariation=="hetero"):
			if(Mvariation=="homo"):
				res+="\t{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\n".format(shape1Ne, shape2Ne, M1current_print, M1past_print)
			if(Mvariation=="hetero"):
				res+="\t{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\n".format(shape1Ne, shape2Ne, M1current_print, shape1mig1, shape2mig1, M1past_print)
	
	for loc in range(nlocus):
		#msnsam tbs 200 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -en tbs 1 tbs -n 2 tbs -en tbs 2 tbs -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs # like "AM"
		print("{0} {1:.5f} {2:.5f} {3} {4} {5} {6:.5f} {7:.5f} {8:.5f} {9:.5f} {10:.5f} {11:.5f} {12:.5f} {13:.5f} {14:.5f} {15:.5f} {16:.5f} {17:.5f} {18:.5f} {19:.5f}".format(int(nspA[loc])+int(nspB[loc]), float(theta[loc]), float(rho[loc]), int(L[loc]), int(nspA[loc]), int(nspB[loc]), M1current[loc], M2current[loc], float(n1priorGenomic[loc]), tauN1past, n1pastprior*float(n1priorGenomic[loc]), float(n2priorGenomic[loc]), tauN2past, n2pastprior*float(n2priorGenomic[loc]), float(TsmallGenomic[loc]), float(M1past[loc]), float(M2past[loc]), float(TsplitGenomic[loc]), float(TsplitGenomic[loc]), float(nApriorGenomic[loc])))

outputfile=open(outputParameters, "w")
outputfile.write(res)
outputfile.close()

