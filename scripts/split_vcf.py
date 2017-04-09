#!/usr/bin/pypy
#script made to split vcf in shorter files for future multi-thread analysis

import sys

vcfFileName = "/home/roux/Documents/ABCheliconius/test_ama10.Hmel2.bwa.default.HC.vcf"
species = "mel"

group_previous_locus = -1

vcfFileName = sys.argv[1] 
species = sys.argv[2]

nContig = 0
cnt = 0
res = {} 

infile = open(vcfFileName, 'r')
for i in infile:
	if "#" not in i:
		contig = i.strip().split("\t")[0]
#		print("locus {0} du groupe {1}: {2}".format(contig, res[contig], cntTMP))
		group_current_locus = res[contig]
		if group_previous_locus == -1:
			tmpVCF = header + i
			group_previous_locus = group_current_locus
			continue
		if group_previous_locus != -1:
			if group_current_locus == group_previous_locus:
				tmpVCF += i
				group_previous_locus = group_current_locus
				continue
			if group_current_locus != group_previous_locus:
				outfile = open("{0}_{1}_subVCF.vcf".format(species, group_previous_locus), "w")
				outfile.write(tmpVCF)
				outfile.close()
				tmpVCF = header
				group_previous_locus = group_current_locus
				continue
			continue
		continue
	if "#CHROM" in i:
		header = i
		continue
	if "##contig" in i:
		i = i.strip()
		i = i.split("=")[2]
		i = i.split(",")[0]
		nContig += 1
		if i not in res:
			res[i] = cnt
			continue
		if nContig%50 == 0:
			cnt += 1
			continue
infile.close()


outfile = open("{0}_{1}_subVCF.vcf".format(species, group_current_locus), "w")
outfile.write(tmpVCF)
outfile.close()


res2 = ""
for i in res:
	res2 += "{0}\t{1}".format(i, res[i])

outfile = open("{0}_correspondance_btw_contigs_subVCF.txt", "w")
outfile.write(res2)
outfile.close()

