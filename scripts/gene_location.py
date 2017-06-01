#!/usr/bin/python
import os.path
import sys

infileName = "Hmel2.gff"

testFile = os.path.exists(infileName)

if testFile == False:
	print("\n\treads through the Hmel2.gff file and gives the location of coding genes.\n\tAll genes are counted but only those with CDS are outputed")
	sys.exit("\n\033[1;31m \tERROR: Hmel2.gff is missing\033[0m\n")

infile = open(infileName, "r")

outfile = open("gene_location.txt", "w")
outfile.write("geneName\tcontig\tstart\tend\n")

list_of_contigs = []
test = 0
cnt = 0

for i in infile:
	if i[0] != "#":
		tmp = i.strip().split("\t")
		contig = tmp[0]
		feature = tmp[2]
		
		if contig not in list_of_contigs:
			list_of_contigs.append(contig)
			cnt = 0 
		
		if feature == "gene":
			test = 0
			start = tmp[3]
			end = tmp[4]
			cnt += 1
		
		if feature == "CDS" and test == 0:
			test = 1
			res = "{0}\t{1}\t{2}\t{3}\n".format(contig + "_" + str(cnt), contig, start, end)
			outfile.write(res)

outfile.close()
infile.close()

