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
model: \033[1;35m=SI\033[1;m (Strict Isolation), \033[1;35m=IM\033[1;m (Isolation with Migration), \033[1;35m=AM\033[1;m (Ancient Migration) or \033[1;35m=SC\033[1;m (Secondary Contact)\n\t\t\
Nvariation: \033[1;35m=homo\033[1;m (shared values of Ne throughout genome for N1, N2 and Nanc) or \033[1;35m=hetero\033[1;m (variation of Ne throughout genome for N1, N2 and Nanc)\n\t\t\
Mvariation: \033[1;35m=homo\033[1;m (shared values of M throughout genome) or \033[1;35m=hetero\033[1;m (variation of M throughout genome)\n\t\t\
nreps: number of multilocus simulations. Ex \033[1;35mnreps=1000\033[1;m\n\t\t\
symMig: \033[1;35m=sym\033[1;m (M1=M2) or \033[1;35m=asym\033[1;m (M1 and M2 are independently chosen)\n\n\t\
Ex:\n\t\
\033[1;32m./priorgen.py bpfile=bpfile_test.txt n1=0 n1=1 n1past=0.1 n1past=2 n2=1 n2=2 n2past=1 n2past=1 nA=2 nA=3 tau=3 tau=4\n\t\tM1=4 M1=5 M2=5 M2=6 shape1=0 shape1=10 shape2=0 shape2=100 model=IM\n\t\tnreps=2 Nvariation=homo Mvariation=hetero symMig=asym parameters=output.txt\033[1;m\n\n\t\
More details about the models in \033[33mhttp://onlinelibrary.wiley.com/doi/10.1111/jeb.12425/abstract\033[1;m and \033[33mhttp://mbe.oxfordjournals.org/content/30/7/1574\033[1;m\n\n\t\
msnsam tbs 20000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -en tbs 1 tbs -n 2 tbs -en tbs 2 tbs -ej tbs 2 1 -eN tbs tbs\t#for 'SI'\n\t\
msnsam tbs 20000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -en tbs 1 tbs -n 2 tbs -en tbs 2 tbs -ej tbs 2 1 -eN tbs tbs\t#for 'IM'\n\t\
msnsam tbs 20000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -en tbs 1 tbs -n 2 tbs -en tbs 2 tbs -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs\t#for 'AM'\n\t\
msnsam tbs 20000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -en tbs 1 tbs -n 2 tbs -en tbs 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs\t#for 'SC'\n\t\
msnsam tbs 20000 -t tbs -r tbs tbs -I 2 tbs tbs 0 -ema 0 2 0 tbs tbs 0 -ema tbs 2 0 tbs tbs 0 -ema tbs 2 0 tbs tbs 0 -ema tbs 2 0 tbs tbs 0 -n 1 tbs -en tbs 1 tbs -n 2 tbs -en tbs 2 tbs -ej tbs 2 1 -eN tbs tbs\t#for 'PAM or PSC'\n\n\t\
camille.roux.1983@gmail.com\n\t\
19/05/2017\n"
for i in sys.argv:
	if("help" in i):
		print(help)
		sys.exit(0)

if len(sys.argv)<=1:
	print(help)
	sys.exit(0)

for i in sys.argv:
	if "=" in i:
		i=i.split("=")
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

if(model=="SI"):
	if(Nvariation=="homo"):
		res="N1\tN1past\ttauN1past\tN2\tN2past\ttauN2past\tNa\tTsplit\n"
	if(Nvariation=="hetero"):
		res="N1\tN1past\ttauN1past\tN2\tN2past\ttauN2past\tNa\tTsplit\tshape1Ne\tshape2Ne\tpropNtrlNe1\tpropNtrlNe2\n"

if(model=="IM"):
	if(Nvariation=="homo"):
		if(Mvariation=="homo"):
			res="N1\tN1past\ttauN1past\tN2\tN2past\ttauN2past\tNa\tTsplit\tM1\tM2\n"
		if(Mvariation=="hetero"):
			res="N1\tN1past\ttauN1past\tN2\tN2past\ttauN2past\tNa\tTsplit\tM1\tM2\tshape1M1\tshape2M1\tshape1M2\tshape2M2\tpropNtrlM1\tpropNtrlM2\n"
	if(Nvariation=="hetero"):
		if(Mvariation=="homo"):
			res="N1\tN1past\ttauN1past\tN2\tN2past\ttauN2past\tNa\tTsplit\tshape1Ne\tshape2Ne\tpropNtrlNe1\tpropNtrlNe2\tM1\tM2\n"
		if(Mvariation=="hetero"):
			res="N1\tN1past\ttauN1past\tN2\tN2past\ttauN2past\tNa\tTsplit\tshape1Ne\tshape2Ne\tpropNtrlNe1\tpropNtrlNe2\tM1\tM2\tshape1M1\tshape2M1\tshape1M2\tshape2M2\tpropNtrlM1\tpropNtrlM2\n"

if(model=="AM" or model=="SC"):
	if "AM" in model:
		tau_name = "Tam"
	if "SC" in model:
		tau_name = "Tsc"
	if(Nvariation=="homo"):
		if(Mvariation=="homo"):
			res="N1\tN1past\ttauN1past\tN2\tN2past\ttauN2past\tNa\tTsplit\t{0}\tM1\tM2\n".format(tau_name)
		if(Mvariation=="hetero"):
			res="N1\tN1past\ttauN1past\tN2\tN2past\ttauN2past\tNa\tTsplit\t{0}\tM1\tM2\tshape1M1\tshape2M1\tshape1M2\tshape2M2\tpropNtrlM1\tpropNtrlM2\n".format(tau_name)
	if(Nvariation=="hetero"):
		if(Mvariation=="homo"):
			res="N1\tN1past\ttauN1past\tN2\tN2past\ttauN2past\tNa\tTsplit\t{0}\tshape1Ne\tshape2Ne\tpropNtrlNe1\tpropNtrlNe2\tM1\tM2\n".format(tau_name)
		if(Mvariation=="hetero"):
			res="N1\tN1past\ttauN1past\tN2\tN2past\ttauN2past\tNa\tTsplit\t{0}\tshape1Ne\tshape2Ne\tpropNtrlNe1\tpropNtrlNe2\tM1\tM2\tshape1M1\tshape2M1\tshape1M2\tshape2M2\tpropNtrlM1\tpropNtrlM2\n".format(tau_name)

if(model=="PSC" or model=="PAM"):
	if(Nvariation=="homo"):
		if(Mvariation=="homo"):
			res="N1\tN1past\ttauN1past\tN2\tN2past\ttauN2past\tNa\tTsplit\tM1\tM2\n"
		if(Mvariation=="hetero"):
			res="N1\tN1past\ttauN1past\tN2\tN2past\ttauN2past\tNa\tTsplit\tM1\tM2\tshape1M1\tshape2M1\tshape1M2\tshape2M2\tpropNtrlM1\tpropNtrlM2\n"
	if(Nvariation=="hetero"):
		if(Mvariation=="homo"):
			res="N1\tN1past\ttauN1past\tN2\tN2past\ttauN2past\tNa\tTsplit\tshape1Ne\tshape2Ne\tpropNtrlNe1\tpropNtrlNe2\tM1\tM2\n"
		if(Mvariation=="hetero"):
			res="N1\tN1past\ttauN1past\tN2\tN2past\ttauN2past\tNa\tTsplit\tshape1Ne\tshape2Ne\tpropNtrlNe1\tpropNtrlNe2\tM1\tM2\tshape1M1\tshape2M1\tshape1M2\tshape2M2\tpropNtrlM1\tpropNtrlM2\n"


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
		n1priorGenomic=binomBeta(nlocus=nlocus, shape1=shape1Ne, shape2=shape2Ne, scalar=n1prior)
		n2priorGenomic=binomBeta(nlocus=nlocus, shape1=shape1Ne, shape2=shape2Ne, scalar=n2prior)
		nApriorGenomic=binomBeta(nlocus=nlocus, shape1=shape1Ne, shape2=shape2Ne, scalar=nAprior)
	
	if(Nvariation=="homo"):
		n1priorGenomic={}
		n2priorGenomic={}
		nApriorGenomic={}
		n1priorGenomic["values"]=[n1prior]*nlocus
		n2priorGenomic["values"]=[n2prior]*nlocus
		nApriorGenomic["values"]=[nAprior]*nlocus
	
	if(Mvariation=="hetero"):
		M1priorGenomic=binomBeta(nlocus=nlocus, shape1=shape1mig1, shape2=shape2mig1, scalar=M1prior)
		if(sym=="sym"):
			M2priorGenomic={}
#			M2priorGenomic=M1priorGenomic["values"]
			M2priorGenomic=M1priorGenomic.copy()	# Sophie Galina
			M2prior=M1prior
		if(sym=="asym"):
			M2priorGenomic=binomBeta(nlocus=nlocus, shape1=shape1mig2, shape2=shape2mig2, scalar=M2prior)
	
	if(Mvariation=="homo"):
		M1priorGenomic={}
		M1priorGenomic["values"]=[M1prior]*nlocus
		if(sym=="sym"):
			M2priorGenomic={}
#			M2priorGenomic["values"]=M1priorGenomic["values"]
			M2priorGenomic=M1priorGenomic.copy()	# Sophie Galina
			M2prior=M1prior
		if(sym=="asym"):
			M2priorGenomic={}
			M2priorGenomic["values"]=[M2prior]*nlocus
	
	res+="{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}".format(n1prior, n1pastprior, tauN1past, n2prior, n2pastprior, tauN2past, nAprior, Tsplit)
	if(model=="SI"):
		if(Nvariation=="homo"):
			res+=" \n"
		if(Nvariation=="hetero"):
			res+=" \t{0:.5f}\t{1:.5f}\t{2}\t{3}\n".format(shape1Ne, shape2Ne, n1priorGenomic["nNeutre"], n2priorGenomic["nNeutre"])
	
	if(model=="IM"):
		if(Nvariation=="homo"):
			if(Mvariation=="homo"):
				res+=" \t{0:.5f}\t{1:.5f}\n".format(M1prior, M2prior)
			if(Mvariation=="hetero"):
				res+=" \t{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6}\t{7}\n".format(M1prior, M2prior, shape1mig1, shape2mig1, shape1mig2, shape2mig2, M1priorGenomic["nNeutre"], M2priorGenomic["nNeutre"])
		if(Nvariation=="hetero"):
			if(Mvariation=="homo"):
				res+=" \t{0:.5f}\t{1:.5f}\t{2}\t{3}\t{4:.5f}\t{5:.5f}\n".format(shape1Ne, shape2Ne, n1priorGenomic["nNeutre"], n2priorGenomic["nNeutre"], M1prior, M2prior)
			if(Mvariation=="hetero"):
				res+=" \t{0:.5f}\t{1:.5f}\t{2}\t{3}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10}\t{11}\n".format(shape1Ne, shape2Ne, n1priorGenomic["nNeutre"], n2priorGenomic["nNeutre"], M1prior, M2prior, shape1mig1, shape2mig1, shape1mig2, shape2mig2, M1priorGenomic["nNeutre"], M2priorGenomic["nNeutre"])
	
	if(model=="AM" or model=="SC"):
		if(Nvariation=="homo"):
			if(Mvariation=="homo"):
				res+=" \t{2:.5f}\t{0:.5f}\t{1:.5f}\n".format(M1prior, M2prior, Tsmall)
			if(Mvariation=="hetero"):
				res+=" \t{8:.5f}\t{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6}\t{7}\n".format(M1prior, M2prior, shape1mig1, shape2mig1, shape1mig2, shape2mig2, M1priorGenomic["nNeutre"], M2priorGenomic["nNeutre"], Tsmall)
		if(Nvariation=="hetero"):
			if(Mvariation=="homo"):
				res+=" \t{6:.5f}\t{0:.5f}\t{1:.5f}\t{2}\t{3}\t{4:.5f}\t{5:.5f}\n".format(shape1Ne, shape2Ne, n1priorGenomic["nNeutre"], n2priorGenomic["nNeutre"], M1prior, M2prior, Tsmall)
			if(Mvariation=="hetero"):
				res+=" \t{12:.5f}\t{0:.5f}\t{1:.5f}\t{2}\t{3}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10}\t{11}\n".format(shape1Ne, shape2Ne, n1priorGenomic["nNeutre"], n2priorGenomic["nNeutre"], M1prior, M2prior, shape1mig1, shape2mig1, shape1mig2, shape2mig2, M1priorGenomic["nNeutre"], M2priorGenomic["nNeutre"], Tsmall)
	
	if(model=="PAM" or model=="PSC"):
		if(Nvariation=="homo"):
			if(Mvariation=="homo"):
				res+=" \t{2:.5f}\t{0:.5f}\t{1:.5f}\n".format(M1prior, M2prior, Tsmall)
			if(Mvariation=="hetero"):
				res+=" \t{8:.5f}\t{0:.5f}\t{1:.5f}\t{2:.5f}\t{3:.5f}\t{4:.5f}\t{5:.5f}\t{6}\t{7}\n".format(M1prior, M2prior, shape1mig1, shape2mig1, shape1mig2, shape2mig2, M1priorGenomic["nNeutre"], M2priorGenomic["nNeutre"], Tsmall)
		if(Nvariation=="hetero"):
			if(Mvariation=="homo"):
				res+=" \t{0:.5f}\t{1:.5f}\t{2}\t{3}\t{4:.5f}\t{5:.5f}\n".format(shape1Ne, shape2Ne, n1priorGenomic["nNeutre"], n2priorGenomic["nNeutre"], M1prior, M2prior)
			if(Mvariation=="hetero"):
				res+=" \t{0:.5f}\t{1:.5f}\t{2}\t{3}\t{4:.5f}\t{5:.5f}\t{6:.5f}\t{7:.5f}\t{8:.5f}\t{9:.5f}\t{10}\t{11}\n".format(shape1Ne, shape2Ne, n1priorGenomic["nNeutre"], n2priorGenomic["nNeutre"], M1prior, M2prior, shape1mig1, shape2mig1, shape1mig2, shape2mig2, M1priorGenomic["nNeutre"], M2priorGenomic["nNeutre"])
	
	for loc in range(nlocus):
		cout=""
		if(model=="SI"):
			#msnsam tbs 200 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -en tbs 1 tbs -n 2 tbs -en tbs 2 tbs -ej tbs 2 1 -eN tbs tbs # for 'SI'
			cout+="{0} {1:.5f} {2:.5f} {3} {4} {5} {6} {7} {8:.5f} {9:.5f} {10:.5f} {11:.5f} {12:.5f} {13:.5f} {14:.5f} {15:.5f} {16:.5f}".format(int(nspA[loc])+int(nspB[loc]), float(theta[loc]), float(rho[loc]), int(L[loc]), int(nspA[loc]), int(nspB[loc]), 0, 0, float(n1priorGenomic["values"][loc]), tauN1past, n1pastprior*float(n1priorGenomic["values"][loc]), float(n2priorGenomic["values"][loc]), tauN2past, n2pastprior*float(n2priorGenomic["values"][loc]), float(TsplitGenomic[loc]), float(TsplitGenomic[loc]), float(nApriorGenomic["values"][loc]))
		
		if(model=="IM"):
			#msnsam tbs 200 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -en tbs 1 tbs -n 2 tbs -en tbs 2 tbs -ej tbs 2 1 -eN tbs tbs # for 'IM'
                        cout+="{0} {1:.5f} {2:.5f} {3} {4} {5} {6:.5f} {7:.5f} {8:.5f} {9:.5f} {10:.5f} {11:.5f} {12:.5f} {13:.5f} {14:.5f} {15:.5f} {16:.5f}".format(int(nspA[loc])+int(nspB[loc]), float(theta[loc]), float(rho[loc]), int(L[loc]), int(nspA[loc]), int(nspB[loc]), float(M1priorGenomic["values"][loc]), float(M2priorGenomic["values"][loc]), float(n1priorGenomic["values"][loc]), tauN1past, n1pastprior*float(n1priorGenomic["values"][loc]), float(n2priorGenomic["values"][loc]), tauN2past, n2pastprior*float(n2priorGenomic["values"][loc]), float(TsplitGenomic[loc]), float(TsplitGenomic[loc]), float(nApriorGenomic["values"][loc]))
		
		if(model=="AM"):
			#msnsam tbs 200 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 0 -m 2 1 0 -n 1 tbs -en tbs 1 tbs -n 2 tbs -en tbs 2 tbs -ema tbs 2 0 tbs tbs 0 -ej tbs 2 1 -eN tbs tbs # for "AM"
                        cout+="{0} {1:.5f} {2:.5f} {3} {4} {5} {6:.5f} {7:.5f} {8:.5f} {9:.5f} {10:.5f} {11:.5f} {12:.5f} {13:.5f} {14:.5f} {15:.5f} {16:.5f} {17:.5f}".format(int(nspA[loc])+int(nspB[loc]), float(theta[loc]), float(rho[loc]), int(L[loc]), int(nspA[loc]), int(nspB[loc]), float(n1priorGenomic["values"][loc]), tauN1past, n1pastprior*float(n1priorGenomic["values"][loc]), float(n2priorGenomic["values"][loc]), tauN2past, n2pastprior*float(n2priorGenomic["values"][loc]), float(TsmallGenomic[loc]), float(M1priorGenomic["values"][loc]), float(M2priorGenomic["values"][loc]), float(TsplitGenomic[loc]), TsplitGenomic[loc], float(nApriorGenomic["values"][loc]))
		
		if(model=="SC"):
			#msnsam tbs 200 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -en tbs 1 tbs -n 2 tbs -en tbs 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs # for 'SC'
                        cout+="{0} {1:.5f} {2:.5f} {3} {4} {5} {6:.5f} {7:.5f} {8:.5f} {9:.5f} {10:.5f} {11:.5f} {12:.5f} {13:.5f} {14:.5f} {15:.5f} {16:.5f} {17:.5f}".format(int(nspA[loc])+int(nspB[loc]), float(theta[loc]), float(rho[loc]), int(L[loc]), int(nspA[loc]), int(nspB[loc]), float(M1priorGenomic["values"][loc]), float(M2priorGenomic["values"][loc]), float(n1priorGenomic["values"][loc]), tauN1past, n1pastprior*float(n1priorGenomic["values"][loc]), float(n2priorGenomic["values"][loc]), tauN2past, n2pastprior*float(n2priorGenomic["values"][loc]), float(TsmallGenomic[loc]), float(TsplitGenomic[loc]), float(TsplitGenomic[loc]), float(nApriorGenomic["values"][loc]))
		
		if(model=="PAM"):
			#msnsam tbs 200 -t tbs -r tbs tbs -I 2 tbs tbs 0 -ema 0 2 0 tbs tbs 0 -ema tbs 2 0 tbs tbs 0 -ema tbs 2 0 tbs tbs 0 -ema tbs 2 0 tbs tbs 0 -n 1 tbs -en tbs 1 tbs -n 2 tbs -en tbs 2 tbs -ej tbs 2 1 -eN tbs tbs # for 'PAM'
			Tsplit_tmp = float(TsplitGenomic[loc])
			M1_tmp = float(M1priorGenomic["values"][loc])
			M2_tmp = float(M2priorGenomic["values"][loc])
			cout+="{0} {1:.5f} {2:.5f} {3:.5f} {4:.5f} {5:.5f} {6:.5f} {7:.5f} {8:.5f} {9:.5f} {10:.5f} {11:.5f} {12:.5f} {13:.5f} {14:.5f} {15:.5f} {16:.5f} {17:.5f} {18:.5f} {19:.5f} {20:.5f} {21:.5f} {22:.5f} {23:.5f} {24:.5f} {25:.5f}".format(int(nspA[loc])+int(nspB[loc]), float(theta[loc]), float(rho[loc]), int(L[loc]), int(nspA[loc]), int(nspB[loc]), 0, 0, Tsplit_tmp/4, M1_tmp, M2_tmp, Tsplit_tmp/2, 0, 0, Tsplit_tmp*3/4, M1_tmp, M2_tmp, float(n1priorGenomic["values"][loc]), tauN1past, n1pastprior*float(n1priorGenomic["values"][loc]), float(n2priorGenomic["values"][loc]), tauN2past, n2pastprior*float(n2priorGenomic["values"][loc]), Tsplit_tmp, Tsplit_tmp, float(nApriorGenomic["values"][loc]))
		
		if(model=="PSC"):
			#msnsam tbs 200 -t tbs -r tbs tbs -I 2 tbs tbs 0 -ema 0 2 0 tbs tbs 0 -ema tbs 2 0 tbs tbs 0 -ema tbs 2 0 tbs tbs 0 -ema tbs 2 0 tbs tbs 0 -n 1 tbs -en tbs 1 tbs -n 2 tbs -en tbs 2 tbs -ej tbs 2 1 -eN tbs tbs # for 'PSC'
			Tsplit_tmp = float(TsplitGenomic[loc])
			M1_tmp = float(M1priorGenomic["values"][loc])
			M2_tmp = float(M2priorGenomic["values"][loc])
			cout+="{0} {1:.5f} {2:.5f} {3:.5f} {4:.5f} {5:.5f} {6:.5f} {7:.5f} {8:.5f} {9:.5f} {10:.5f} {11:.5f} {12:.5f} {13:.5f} {14:.5f} {15:.5f} {16:.5f} {17:.5f} {18:.5f} {19:.5f} {20:.5f} {21:.5f} {22:.5f} {23:.5f} {24:.5f} {25:.5f}".format(int(nspA[loc])+int(nspB[loc]), float(theta[loc]), float(rho[loc]), int(L[loc]), int(nspA[loc]), int(nspB[loc]), M1_tmp, M2_tmp, Tsplit_tmp/4, 0, 0, Tsplit_tmp/2, M1_tmp, M2_tmp, Tsplit_tmp*3/4, 0, 0, float(n1priorGenomic["values"][loc]), tauN1past, n1pastprior*float(n1priorGenomic["values"][loc]), float(n2priorGenomic["values"][loc]), tauN2past, n2pastprior*float(n2priorGenomic["values"][loc]), Tsplit_tmp, Tsplit_tmp, float(nApriorGenomic["values"][loc]))
		print(cout)
outputfile=open(outputParameters, "w")
outputfile.write(res)
outputfile.close()

