#!/usr/bin/python
from Bio.SeqIO import parse
from Bio.Seq import Seq
import sys
from os import listdir
from os.path import isfile, join
from numpy.random import choice
from numpy import ceil


if len(sys.argv) != 6:
	print("\n\tproduces bpfile, spinput.txt and locus.ms")
	print("\033[1;33m\tExample: ./make_ABC_inputFiles.py ama_txn 100 900 100000 0.000000003\033[0m")
	sys.exit("\n\033[1;31m\tNeeds 5 arguments\033[0m\n")


pair = sys.argv[1]
Lmin = int(sys.argv[2])
nLoci = int(sys.argv[3])
Nindiv = int(sys.argv[4])
mu = float(sys.argv[5])

#pair = "ama_txn"
#Lmin = 100
#nLoci = 900
#Nindiv = 100000
#mu = 0.00000002

mypath = '/home/croux/Documents/heliconius/inputABC/' + pair

files = [ f for f in listdir(mypath) if isfile(join(mypath, f)) and '.fas' in f ]


listOfGenes = [ i.split(".")[0] for i in files ]


genes = {}
for i in listOfGenes:
	if i not in genes:
		genes[i] = {}
	infile = open("{0}/{1}_info.txt".format(mypath, i), "r")
	header = infile.readline().strip().split("\t")
	info = infile.readline().strip().split("\t")
	infile.close()
	for j in range(len(info)):
		genes[i][header[j]] = info[j]


candidateLoci = {}
for i in genes:
	nStopCodon = 0
	infile = parse("{0}/{1}.fas".format(mypath, i), "fasta")
	for j in infile:
		nStopCodon += j.seq.translate().count("*")
	infile.close()
	if int(genes[i]['nSynSegSite']) != 0 and float(genes[i]['Lsyno']) >= Lmin and nStopCodon == 0:
		candidateLoci[i] = genes[i]
		infile = open("{0}/{1}.ms".format(mypath, i), "r")
		ms = ""
		cnt = 0
		for k in infile:
			cnt += 1
			if cnt > 2:
				ms += k
		infile.close()
		candidateLoci[i]['ms'] = ms


sampledLoci = choice(len(candidateLoci), size = nLoci, replace = False)


bpfile_header = "#\t{0}\t{1}\t{2}\t{3}\n".format(pair.split("_")[0], pair.split("_")[1], Nindiv, mu)
bpfile_line1 = "" # L
bpfile_line2 = "" # nsamA
bpfile_line3 = "" # nsamB
bpfile_line4 = "" # theta
bpfile_line5 = "" # rho

spinput = "\n{0}\n".format(nLoci)

loci_ms = "./msnsam tbs 20 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -ej tbs 1 2 -eN tbs tbs\n3579 27011 59243\n"

for i in sampledLoci:
	name = candidateLoci.keys()[i]
	loci_ms += candidateLoci[name]['ms']
	L = float(candidateLoci[name]['Lsyno'])
	nsamA = candidateLoci[name]['nsamA']
	nsamB = candidateLoci[name]['nsamB']
	theta = 4 * Nindiv * mu * float(L)
	bpfile_line1 += "{0}\t".format(int(ceil(L)))
	bpfile_line2 += "{0}\t".format(nsamA)
	bpfile_line3 += "{0}\t".format(nsamB)
	bpfile_line4 += "{0}\t".format(theta)
	bpfile_line5 += "{0}\t".format(theta)
	spinput += "{0}\n{1}\n{2}\n".format(nsamA, nsamB, int(ceil(float(L))))
spinput	+= "1\nloci.ms\n"

bpfile = bpfile_header + bpfile_line1 + "\n" + bpfile_line2 + "\n" + bpfile_line3 + "\n" + bpfile_line4 + "\n" + bpfile_line5 + "\n"
 
outfile = open("loci.ms", "w")
outfile.write(loci_ms)
outfile.close()


outfile = open("bpfile", "w")
outfile.write(bpfile)
outfile.close()


outfile = open("spinput.txt", "w")
outfile.write(spinput)
outfile.close()

