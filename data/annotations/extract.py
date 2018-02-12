#!/usr/bin/python

import sys

gffFile = sys.argv[1]
contig_target = sys.argv[2]

infile = open(gffFile, "r")

genome = {}

geneName = 0

for i in infile:
	if i[0] == "#":
		if contig_target in i:
			i = i.strip().split(" ")
			genome[i[1]] = {}
			genome[i[1]]["position"] = range(int(i[3]))
			genome[i[1]]["feature_strandFOR"] = ["ig"] * int(i[3])
			genome[i[1]]["geneName_strandFOR"] = ["ig"] * int(i[3])
			genome[i[1]]["feature_strandREV"] = ["ig"] * int(i[3])
			genome[i[1]]["geneName_strandREV"] = ["ig"] * int(i[3])
	else:
		i = i.strip().split("\t")
		contig = i[0]
		if i[2] != "exon":
			feature = i[2]
			if feature == "gene" or "pseudo" in feature:
				geneName += 1
			start = int(i[3])-1
			end = int(i[4])
			strand = i[6]
			if contig == contig_target:
				if strand == "+":
					genome[contig]["feature_strandFOR"][start:end] = [feature]*(end-start)
					genome[contig]["geneName_strandFOR"][start:end] = ["gene_{0}".format(geneName)]*(end-start)
				if strand == "-":
					genome[contig]["feature_strandREV"][start:end] = [feature]*(end-start)
					genome[contig]["geneName_strandREV"][start:end] = ["gene_{0}".format(geneName)]*(end-start)
infile.close()


res = "contig\tposition\tfeature_strandFOR\tgeneName_strandFOR\tfeature_strandREV\tgeneName_strandREV\n"
for i in genome[contig_target]["position"]:
	res += "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(contig_target, i, genome[contig_target]["feature_strandFOR"][i], genome[contig_target]["geneName_strandFOR"][i], genome[contig_target]["feature_strandREV"][i], genome[contig_target]["geneName_strandREV"][i])

outfile = open("annotations_{0}.txt".format(contig_target), "w")
outfile.write(res)
outfile.close()

print("{0} done".format(contig_target))

