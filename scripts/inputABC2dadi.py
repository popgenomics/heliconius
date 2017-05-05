#!/usr/bin/pypy
import os
from os import listdir
from os.path import isfile, join
import sys

# check the arguments
if len(sys.argv) != 2:
        print("\n\tproduces a table with allele frequencies for DaDi. Needs '*.ms' and '*_info.txt' files produced by fasta2ABC.py")
        print("\n\033[1;33m\tExample: ./inputABC2dadi.py ama_txn\033[0m\n")
        print("\t\targ1 =\tname of the directory containin *.ms, *_info.txt and *.fas files")
        sys.exit("\n\033[1;31m ERROR: 1 arguments is required: {0} provided\033[0m\n".format(len(sys.argv)))


pair = sys.argv[1]
#pair = "ama_chi"

#mypath = '/home/croux/Documents/ABCheliconius/pairs/' + pair
#mypath = '/home/croux/Documents/heliconius/inputABC/' + pair
mypath = '/home/croux/Documents/data_heliconius/inputABC/' + pair

test = os.path.isdir(mypath)


if test == False:
        sys.exit("\n\033[1;31m ERROR: {0} is not found\033[0m\n".format(mypath))


files = [ f for f in listdir(mypath) if isfile(join(mypath, f)) and '.ms' in f ]
listOfGenes = [ i.split(".")[0] for i in files ]


outfileName = "/home/croux/Documents/data_heliconius/inputDADI/data_dadi_{0}.txt".format(pair)
outfile = open(outfileName, "w")
res = "Ing\tOut\tAllele1\t{0}\t{1}\tAllele2\t{0}\t{1}\tGene\tPosition\n".format(pair.split("_")[0], pair.split("_")[1])
outfile.write(res)
outfile.close()

res = "" 
for i in listOfGenes:
	#print(i)
        # get nSamA, nSamB and L
        infile = open("{0}/{1}_info.txt".format(mypath, i), "r")
        header = infile.readline().strip().split("\t")
        info = infile.readline().strip().split("\t")
        infile.close()
        nsamA = int(info[header.index("nsamA")])
        nsamB = int(info[header.index("nsamB")])
	# get ms output file
        infile = open("{0}/{1}.ms".format(mypath, i), "r")
	line = infile.readline() # line 1, *.ms
	line = infile.readline() # line 2, *.ms
	line = infile.readline() # line 3, *.ms
	line = infile.readline() # line 4, *.ms
	line = infile.readline() # line 5, *.ms
	nSegSites = int(line.split(":")[1]) # get the number of SegSites
	line = infile.readline() # line 6, *.ms
	speciesA, speciesB = {}, {}
	cnt=0
	for line in infile:
		cnt += 1
		if cnt <= nsamA:
			speciesA[cnt] = line.strip()
		else:
			speciesB[cnt-nsamA] = line.strip()
#	print("locus {0}: {1} speciesA and {2} speciesB, {3} segsites".format(i, len(speciesA), len(speciesB), nSegSites))
	infile.close()
	for snp in range(nSegSites):
		nA_spA, nG_spA, nA_spB, nG_spB = 0, 0, 0, 0
		for ind in speciesA:
			allele = speciesA[ind][snp]
			if allele == '0':
				nA_spA += 1
			if allele == '1':
				nG_spA += 1
		for ind in speciesB:
			allele = speciesB[ind][snp]
			if allele == '0':
				nA_spB += 1
			if allele == '1':
				nG_spB += 1
		res = "-A-\t---\tA\t{0}\t{1}\tG\t{2}\t{3}\t{4}\t{5}\n".format(nA_spA, nA_spB, nG_spA, nG_spB, i, snp+1)
		outfile = open(outfileName, "a")
		outfile.write(res)
		outfile.close()
 
