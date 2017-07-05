#!/usr/bin/python

import os.path
import sys
from Bio.Seq import Seq
from Bio.SeqIO import parse

if len(sys.argv) != 2:
	sys.exit("error: needs one argument = geneName (ex: Hmel200001_1)")


geneName = sys.argv[1]

txtFileName = "{0}.txt".format(geneName)
fastaFileName = "PHASED_{0}.fas".format(geneName)

if os.path.exists(txtFileName) == False:
	sys.exit("the file {0} is not found\n".format(txtFileName))

if os.path.exists(fastaFileName) == False:
	sys.exit("the file {0} is not found\n".format(fastaFileName))

infile = open(txtFileName, "r")
strand = infile.readline().strip().split("\t")[4]
infile.close()

if strand == "+":
	sys.exit()

if strand == "-":
	res = ""
	infile = parse(fastaFileName, "fasta")
	for i in infile:
		res += ">{0}\n".format(i.id)
		res += "{0}\n".format(i.seq.reverse_complement())
	infile.close()
	
	outfile = open(fastaFileName, "w")
	outfile.write(res)
	outfile.close()
	
