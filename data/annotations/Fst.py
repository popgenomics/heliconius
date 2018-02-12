#!/gepv/home2/croux/bin/pypy

import sys

pop = ['ama', 'chi', 'flo', 'mal', 'ros', 'txn', 'vul', 'zel']
statNames = ""
nStats = 0

# write the header for the 8 Hs
for i in pop:
	statNames += "HS_{0}\t".format(i)
	nStats += 1

# write the header with the 28 FST
for i in range(len(pop)-1):
	for j in range(i+1, len(pop), 1):
		statNames += "FST_{0}_{1}\t".format(pop[i], pop[j])
		nStats += 1

# write the header with the categories of sites for the 28 pairs (SxA, SxB, Ss, Sf)
for i in range(len(pop)-1):
	for j in range(i+1, len(pop), 1):
		statNames += "site_{0}_{1}\t".format(pop[i], pop[j])
		nStats += 1
statNames = statNames.strip()


def heteroZ(nRef, nAlt):
	# computes the heterozygosity Hs
	nInd = nRef + nAlt
	num = nRef*nAlt*1.0
	denom = nInd * (nInd - 1) / 2.0
	return(num/denom)


def sites(nRefA, nAltA, nRefB, nAltB, site):
	# determines the category of site
	if 0 in [nRefA, nAltA]: # if monomorphic in A
		if 0 in [nRefB, nAltB]: # if monomorphic in B
			if nRefA==0: # if Alt allele fixed in A
				if nRefB==0: # if Alt allele fixed in B
					site.append("same")
				else: # if Ref allele fixed in B
					site.append("Sf")
			else: # if Ref allele fixed in A
				if nAltB==0: # if Ref allele fixed in B
					site.append("same")
				else: # if Alt allele fixed in B
					site.append("Sf")
		else: # if polymorphic in B
			site.append("SxB")
	else: # if polymorphic in A
		if 0 in [nRefB, nAltB]: # if monomorphic in B
			site.append("SxA")
		else: # if polymorphic in B
			site.append("Ss")
		

def treatLine(x, pop, header, nStats):
	seuil = 10 # maximum number of N tolerated at a position within a species
	x = x.split("\t")
	if x[6] == "N" or x[6] == "multiallelic":
			res = "NA\t" * nStats
			res = res.strip() + "\n"
			return(res)
	else:
		if x[7] == ".":
			res = "0\t" * 8
			res += "-9\t" * 28
			res += "same\t" * 28
			res = res.strip() + "\n"
			return(res)
		else:
			Fst = []
			HS = []
			site = []
			for popA in pop: # loop to compute Hs within populations
				index_nN = header.index("nN_{0}".format(popA))
				index_nRef = header.index("nRefA_{0}".format(popA))
				index_nAlt = header.index("nAltA_{0}".format(popA))
				nN = int(x[index_nN])
				nRef = int(x[index_nRef])
				nAlt = int(x[index_nAlt])
				if nN<seuil: # if not too much 'N' in the alignment
					HS.append(heteroZ(nRef, nAlt))
				else: # if too much 'N'
					HS.append('NA')
			for i in range(len(pop)-1): # loop to compute Fst and Sites among pairs of populations: i = population A
				for j in range(i+1, len(pop), 1): # j = population B
					index_nN_A = header.index("nN_{0}".format(pop[i]))
					index_nRef_A = header.index("nRefA_{0}".format(pop[i]))
					index_nAlt_A = header.index("nAltA_{0}".format(pop[i]))
					nN_A = int(x[index_nN_A])
					nRef_A = int(x[index_nRef_A])
					nAlt_A = int(x[index_nAlt_A])
					
					index_nN_B = header.index("nN_{0}".format(pop[j]))
					index_nRef_B = header.index("nRefA_{0}".format(pop[j]))
					index_nAlt_B = header.index("nAltA_{0}".format(pop[j]))
					nN_B = int(x[index_nN_B])
					nRef_B = int(x[index_nRef_B])
					nAlt_B = int(x[index_nAlt_B])
					
					if nN_A<seuil and nN_B<seuil:
						nRef_tot = nRef_A + nRef_B
						nAlt_tot = nAlt_A + nAlt_B
						HS_mean = (HS[i] + HS[j])/2.0
						Htot = heteroZ(nRef_tot, nAlt_tot)
						sites(nRef_A, nAlt_A, nRef_B, nAlt_B, site)
						if Htot != 0:
							Fst.append(1-HS_mean/Htot)
						else:
							Fst.append(-9)
					else:
						Fst.append('NA')
						site.append('NA')
			res = "\t".join([ str(i) for i in HS+Fst+site ]) + "\n"
			return(res)


infile = open("{0}_freq.txt".format(sys.argv[1]), "r")

tmp = infile.readline().strip()
header = tmp.split("\t")

res = "{0}\t{1}\n".format(tmp, statNames)
outfile = open("{0}_Fst.txt".format(sys.argv[1]), "w")
outfile.write(res)
outfile.close()

outfile = open("{0}_Fst.txt".format(sys.argv[1]), "a")
for i in infile:
	i = i.strip()
	res = i + "\t" + treatLine(i, pop, header, nStats) 
	outfile.write(res)
infile.close()
outfile.close()


