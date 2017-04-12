#!/usr/bin/python

import sys

if len(sys.argv) > 1:
	print("\n\tparse_polydNdS.py read the output of Kevin R. Thornton's polydNdS (https://github.com/molpopgen/analysis)\n")
	print("\033[1;33m\tExample: polydNdS -i txn_Hmel201001_4.fas -N | parse_polydNdS.py \033[0m")
	sys.exit("\n\033[1;31m\tReads stdin, doesn't need arguments\033[0m\n")

cnt = -1

test = ""
for line in sys.stdin:
	cnt += 1
	line = line.strip()
	if cnt == 0:
		locus = line
		if "/" in locus:
			locus = locus.split("/")[-1:][0]
		continue
	if "Mean # of replacement sites" in line:
		nNonSyn = line.replace(" ", "").split("=")[1]
		continue
	if "Mean # of synonymous sites" in line:
		nSyn = line.replace(" ", "").split("=")[1]
		continue
	if "The # of third positions" in line:
		nThird = line.replace(" ", "").split("=")[1]
		continue
	if "Mean fourfold degenerate sites" in line:
		nFourFold = line.split(" ")[-1:]
		continue
	if line=="Polymorphism in entire rgion:":
		test = "allRegion"
		continue
	if line=="Replacement polymorphisms in coding region:":
		test = "nonSyn"
		continue
	if line=="Silent polymorphisms in coding region:":
		test = "syn"
		continue
	if line=="Polymorphisms at third positions":
		test = "third"
		continue
	if line=="All fourfold-degenerate positions:":
		test = "fourFold"
		continue
	if test == "allRegion":
		if "Num_Sites" in line:
			nSites = line.split("\t")[-1:]
			continue
		if "ThetaW" in line:
			theta_total = line.split("\t")[-1:]
			continue
		if "ThetaPi" in line:
			pi_total = line.split("\t")[-1:]
			test=""
			continue
		continue
	if test == "nonSyn":
		if "Num_Sites" in line:
			nSites_N = line.split("\t")[-1:]
			continue
		if "ThetaW" in line:
			theta_N = line.split("\t")[-1:]
			continue
		if "ThetaPi" in line:
			pi_N = line.split("\t")[-1:]
			test=""
			continue
		continue
	if test == "syn":
		if "Num_Sites" in line:
			nSites_S = line.split("\t")[-1:]
			continue
		if "ThetaW" in line:
			theta_S = line.split("\t")[-1:]
			continue
		if "ThetaPi" in line:
			pi_S = line.split("\t")[-1:]
			test=""
			continue
		continue
	if test == "third":
		if "Num_Sites" in line:
			nSites_third = line.split("\t")[-1:]
			continue
		if "ThetaW" in line:
			theta_third = line.split("\t")[-1:]
			continue
		if "ThetaPi" in line:
			pi_third = line.split("\t")[-1:]
			test=""
			continue
		continue
	if test == "fourFold":
		if "Num_Sites" in line:
			nSites_fourF = line.split("\t")[-1:]
			continue
		if "ThetaW" in line:
			theta_fourF = line.split("\t")[-1:]
			continue
		if "ThetaPi" in line:
			pi_fourF = line.split("\t")[-1:]
			test=""
			continue
		continue

header = "locus\tnSites_tot\tnSites_N\tnSites_S\tnSites_third\tnSites_fourF\tthetaW_total\tpi_total\tthetaW_N\tpi_N\tthetaW_S\tpi_S\tthetaW_third\tpi_third\tthetaW_fourF\tpi_fourF\n"
res = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\t{12}\t{13}\t{14}\t{15}".format(locus, nSites[0], nNonSyn, nSyn, nThird, nFourFold[0], theta_total[0], pi_total[0], theta_N[0], pi_N[0], theta_S[0], pi_S[0], theta_third[0], pi_third[0], theta_fourF[0], pi_fourF[0])

print(header+res)
#Polymorphisms at third positions
#Segsites        9
#Mutations       9
#Singletons      3
#Num_Sites       604
#ThetaW/site     0.00420004
#ThetaPi/site    0.00402579
#
#All fourfold-degenerate positions:
#Segsites        6
#Mutations       6
#Singletons      3
#Num_Sites       382.8
#ThetaW/site     0.00441802
#ThetaPi/site    0.00288731
#
#All silent polymorphisms (coding + noncoding):
#Segsites        9
#Mutations       9
#Singletons      3
#Num_Sites       469.078
#ThetaW/site     0.00540811
#ThetaPi/site    0.00518374
#

