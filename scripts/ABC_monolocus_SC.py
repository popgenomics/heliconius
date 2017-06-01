#!/usr/bin/python
from numpy.random import beta
from numpy.random import choice
import os

posteriorFile = "Posterior_v4_SC_Mhetero_Nhetero_Beta_1"
#N1, N1past, tauN1past, N2, N2past, tauN2past, Na, Tsplit, Tsc, shape1Ne, shape2Ne, M1, M2, shape1M1, shape2M1, shape1M2, shape2M2
nRep = 100000
nsamA = 20
nsamB = 20
L = 175
mu = 3e-9
N = 100000
theta = 4*100000*mu*L
rho = 4*100000*mu*L

posterior = {}
for i in range(17):
	posterior[i] = []


infile = open(posteriorFile, "r")
for i in infile:
	i = i.strip().split(" ")
	for j in range(17):
		posterior[j].append(float(i[j]))
infile.close()


os.system("mkdir migration")

spinput = "\n{0}\n{1}\n{2}\n{3}\n{4}\nmyfifo".format(1, nsamA, nsamB, L, nRep)

for sub_rep in range(10):
	os.system("mkdir migration/rep_{0}".format(sub_rep))
	outfile = open("migration/rep_{0}/spinput.txt".format(sub_rep), "w")
	outfile.write(spinput)
	outfile.close()
	
	outfile = open("migration/rep_{0}/prior_monolocus_SC_migration_{0}.txt".format(sub_rep), "w")
	nPosterior = len(posterior[0])
	for i in range(nRep):
		setOfParam = choice(nPosterior)
		N1 = posterior[0][setOfParam] * beta(posterior[9][setOfParam], posterior[10][setOfParam])
		N1past = N1 * posterior[1][setOfParam]
		tauN1past = posterior[2][setOfParam]
		
		N2 = posterior[3][setOfParam] * beta(posterior[9][setOfParam], posterior[10][setOfParam])
		N2past = N2 * posterior[4][setOfParam]
		tauN2past = posterior[5][setOfParam]
		
		Na = posterior[6][setOfParam] * beta(posterior[9][setOfParam], posterior[10][setOfParam])
		
		Tsplit = posterior[7][setOfParam]
		
		Tsc = posterior[8][setOfParam]
		
		M1 = posterior[11][setOfParam] * beta(posterior[13][setOfParam], posterior[14][setOfParam])
		M2 = posterior[12][setOfParam] * beta(posterior[15][setOfParam], posterior[16][setOfParam])
		
		# tbs nRep -t tbs -r tbs tbs -I 2 tbs tbs -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -en tbs 1 tbs -en tbs 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs
		res = "\t".join([ str(k) for k in [nsamA+nsamB, theta, rho, L, nsamA, nsamB, M1, M2, N1, N2, tauN1past, N1past, tauN2past, N2past, Tsc, Tsplit, Tsplit, Na] ]) + "\n"
		outfile.write(res)
	outfile.close()


os.system("mkdir isolation")
for sub_rep in range(10):
	os.system("mkdir isolation/rep_{0}".format(sub_rep))
	outfile = open("isolation/rep_{0}/spinput.txt".format(sub_rep), "w")
	outfile.write(spinput)
	outfile.close()
	
	outfile = open("isolation/rep_{0}/prior_monolocus_SC_isolation_{0}.txt".format(sub_rep), "w")
	nPosterior = len(posterior[0])
	for i in range(nRep):
		setOfParam = choice(nPosterior)
		N1 = posterior[0][setOfParam] * beta(posterior[9][setOfParam], posterior[10][setOfParam])
		N1past = N1 * posterior[1][setOfParam]
		tauN1past = posterior[2][setOfParam]
		
		N2 = posterior[3][setOfParam] * beta(posterior[9][setOfParam], posterior[10][setOfParam])
		N2past = N2 * posterior[4][setOfParam]
		tauN2past = posterior[5][setOfParam]
		
		Na = posterior[6][setOfParam] * beta(posterior[9][setOfParam], posterior[10][setOfParam])
		
		Tsplit = posterior[7][setOfParam]
		
		Tsc = posterior[8][setOfParam]
		
		M1 = 0
		M2 = 0
		
		# tbs nRep -t tbs -r tbs tbs -I 2 tbs tbs -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -en tbs 1 tbs -en tbs 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs
		res = "\t".join([ str(k) for k in [nsamA+nsamB, theta, rho, L, nsamA, nsamB, M1, M2, N1, N2, tauN1past, N1past, tauN2past, N2past, Tsc, Tsplit, Tsplit, Na] ]) + "\n"
		outfile.write(res)
	outfile.close()


print('bsub -q normal <<<"mknod myfifo p; mscalc < myfifo & msnsam tbs {0} -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -en tbs 1 tbs -n 2 tbs -en tbs 2 tbs -eM tbs 0 -ej tbs 2 1 -eN tbs tbs >myfifo"'.format(nRep))



