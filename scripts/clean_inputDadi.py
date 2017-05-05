#!/usr/bin/pypy

# script used to clean the SNP table input for dadi (from vcf2dadi.py)
import sys

fileName = sys.argv[1]

infile = open(fileName, "r")

header = infile.readline()
lenHeader = len(header.strip().split("\t"))

outfile = open("cleaned_" + fileName, "w")
outfile.write(header)
outfile.close()

outfile = open("cleaned_" + fileName, "a")

res = header
for i in infile:
	j = i.strip().split("\t")
	if len(j)==lenHeader:
		if int(j[3])+int(j[4]) > 0:
			if int(j[6])+int(j[7]) > 0:
				outfile.write(i)

infile.close()
outfile.close()
