#!/usr/local/Cluster-Apps/python/2.7.5/bin/python

import sys
from numpy.random import choice
import os.path


# check the arguments
if len(sys.argv) != 3:
	print("\n\tsample N SNPs from SNP data table")
	print("\n\033[1;33m\tExample: ./snp_sampler.py cleaned_output_ama_chi.txt 1000000\033[0m\n")
	print("\t\targ1 =\tname of the table containing SNP data")
	print("\t\targ2 =\tnumber of SNPs to sample")
	sys.exit("\n\033[1;31m ERROR: 2 arguments are required: {0} provided\033[0m\n".format(len(sys.argv)-1))


# get the params
infileName = sys.argv[1]
nSNPs = int(sys.argv[2])


# test file
test = os.path.exists(infileName)
if test==False:
	sys.exit("\n\033[1;31m ERROR: {0} is not found\033[0m\n".format(infileName))


nLines=0
infile = open(infileName, "r")
for i in infile:
	nLines += 1
nLines -= 1 # remove the header
infile.close()


# start the sampling
list_of_SNPs = choice(nLines, nSNPs, replace=False) 

outfileName = "subSampled_" + infileName
outfile = open(outfileName, "w")


infile = open(infileName, "r")
header = infile.readline()
outfile.write(header)
cnt = -1
for i in infile:
	cnt += 1
	if cnt in list_of_SNPs:
		outfile.write(i)
infile.close()
outfile.close()


