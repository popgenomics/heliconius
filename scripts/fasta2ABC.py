#!/usr/bin/python
import sys

# nN = number of non-synonymous sites in the codon i: example for CGG -> nN = 2/3 + 3/3 + 0/3
# nS = number of synonymous sites in the codon i: example for CGG -> n> = 1/3 + 0/3 + 3/3
codonTable = {'AAA': {'aa': 'K', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AAC': {'aa': 'N', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AAG': {'aa': 'K', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AAT': {'aa': 'N', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'ACA': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'ACC': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'ACG': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'ACT': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'AGA': {'aa': 'R', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'AGC': {'aa': 'S', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AGG': {'aa': 'R', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'AGT': {'aa': 'S', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'ATA': {'aa': 'I', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'ATC': {'aa': 'I', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'ATG': {'aa': 'M', 'nN': 3.0, 'nS': 0.0}, 'ATT': {'aa': 'I', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'CAA': {'aa': 'Q', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CAC': {'aa': 'H', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CAG': {'aa': 'Q', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CAT': {'aa': 'H', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CCA': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CCC': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CCG': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CCT': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CGA': {'aa': 'R', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CGC': {'aa': 'R', 'nN': 2.0, 'nS': 1.0}, 'CGG': {'aa': 'R', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CGT': {'aa': 'R', 'nN': 2.0, 'nS': 1.0}, 'CTA': {'aa': 'L', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CTC': {'aa': 'L', 'nN': 2.0, 'nS': 1.0}, 'CTG': {'aa': 'L', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CTT': {'aa': 'L', 'nN': 2.0, 'nS': 1.0}, 'GAA': {'aa': 'E', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GAC': {'aa': 'D', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GAG': {'aa': 'E', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GAT': {'aa': 'D', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GCA': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GCC': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GCG': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GCT': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GGA': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GGC': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GGG': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GGT': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GTA': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'GTC': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'GTG': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'GTT': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'TAC': {'aa': 'Y', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TAT': {'aa': 'Y', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TCA': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TCC': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TCG': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TCT': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TGC': {'aa': 'C', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TGG': {'aa': 'W', 'nN': 3.0, 'nS': 0.0}, 'TGT': {'aa': 'C', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TTA': {'aa': 'L', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'TTC': {'aa': 'F', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TTG': {'aa': 'L', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'TTT': {'aa': 'F', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}}


def fasta2dic(f):
	fasta = open(f).readlines()
	seqName = [x.split(" ")[0].rstrip().replace('>','') for x in fasta if x[0] == '>']
	seq = ''.join([x.rstrip() if x[0]!='>' else '@' for x in fasta])[1:].split('@')
	res = {}
	for i in range(len(seq)):
		res[seqName[i]] = seq[i]
	return (res)


seqA = sys.argv[1] # alignment of sequences from species A
seqB = sys.argv[2] # alignment of sequences from species B
threshold_N = float(sys.argv[3]) # if an allele has %N > threshold_N --> sequence is rejected
nMin = int(sys.argv[4]) # minimum number of individuals within a species


alignA = fasta2dic(seqA)
alignB = fasta2dic(seqB)


L = 1.0 * len(alignA[alignA.keys()[0]])


interspe = {}
nA = 0
nB = 0
for i in alignA:
	propN = alignA[i].count("N")/L
	if propN < threshold_N:
		nA += 1
		interspe[i] = alignA[i]


for i in alignB:
	propN = alignB[i].count("N")/L
	if propN < threshold_N:
		nB += 1
		interspe[i] = alignB[i]


if nA < nMin or nB < nMin:
	exit(0)


