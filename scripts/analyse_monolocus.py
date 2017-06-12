#!/usr/bin/python
import os
from os import listdir
from os.path import isfile, join
import sys
from numpy import ceil
from Bio.SeqIO import parse
from Bio.Seq import Seq

# check the arguments
if len(sys.argv) != 2:
	print("\n\tproduces a table with summary statistics for all loci")
	print("\n\033[1;33m\tExample: ./analyse_monolocus.py ama_txn\033[0m\n")
	print("\t\targ1 =\tname of the directory containin *.ms, *_info.txt and *.fas files")
	sys.exit("\n\033[1;31m ERROR: 1 arguments is required: {0} provided\033[0m\n".format(len(sys.argv)))


pair = sys.argv[1]


#mypath = '/home/croux/Documents/ABCheliconius/pairs/' + pair
#mypath = '/home/croux/Documents/heliconius/inputABC/' + pair
mypath = '/home/croux/Documents/data_heliconius/inputABC/' + pair

test = os.path.isdir(mypath)


if test == False:
	sys.exit("\n\033[1;31m ERROR: {0} is not found\033[0m\n".format(mypath))


files = [ f for f in listdir(mypath) if isfile(join(mypath, f)) and '.ms' in f ]
listOfGenes = [ i.split(".")[0] for i in files ]

res = ""
for i in listOfGenes:
	# get nSamA, nSamB and L
	infile = open("{0}/{1}_info.txt".format(mypath, i), "r")
	header = infile.readline().strip().split("\t")
	info = infile.readline().strip().split("\t")
	infile.close()
	nsamA = int(info[header.index("nsamA")])
	nsamB = int(info[header.index("nsamB")])
	Lsyno = int(ceil(float(info[header.index("Lsyno")])))
	# get L syno
	nStopCodon = 0
	infile = parse("{0}/{1}.fas".format(mypath, i), "fasta")
	for j in infile:
		nStopCodon += j.seq.translate().count("*")
	infile.close()
	# get summary statistics
	spinput = "\n1\n{0}\n{1}\n{2}\n1\n{3}.ms\n".format(nsamA, nsamB, Lsyno, mypath + "/" + i)
	outfile = open("spinput.txt", "w")
	outfile.write(spinput)
	outfile.close()
	os.system("mscalc") 
	infile = open("ABCstat.txt", "r")
	header = infile.readline()
	header = infile.readline().strip()
	stats = infile.readline().strip().split("\t")
	infile.close()
	stats[0] = i
	stats = "\t".join(stats)
	if res == "":
		res = header + "\t" + "nSamA\tnSamB\tLsyno\n"
	res += "{0}\t{1}\t{2}\t{3}\n".format(stats, nsamA, nsamB, Lsyno)

outfile = open("summary_stats.txt", "w")
outfile.write(res)
outfile.close()

