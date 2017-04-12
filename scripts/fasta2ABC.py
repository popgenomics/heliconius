#!/usr/bin/pypy
import sys

def coloredSeq(seq):
	# print sequences with the standard color code
	seq = seq.replace("A", '\x1b[5;30;41m' + 'A' + '\x1b[0m')
	seq = seq.replace("T", '\x1b[5;30;44m' + 'T' + '\x1b[0m')
	seq = seq.replace("G", '\x1b[6;30;43m' + 'G' + '\x1b[0m')
	seq = seq.replace("C", '\x1b[5;30;42m' + 'C' + '\x1b[0m')
	return(seq)


def trunc2triplets(size):
	# trunc a value to its closest and smaller multiple of 3
	size = int(size)
	for i in range(2):
		if size%3 != 0:
			size -= 1
	return(size)


# nN = number of non-synonymous sites in the codon i: example for CGG -> nN = 2/3 + 3/3 + 0/3
# nS = number of synonymous sites in the codon i: example for CGG -> n> = 1/3 + 0/3 + 3/3
codonTable = {'AAA': {'aa': 'K', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AAC': {'aa': 'N', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AAG': {'aa': 'K', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AAT': {'aa': 'N', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'ACA': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'ACC': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'ACG': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'ACT': {'aa': 'T', 'nN': 2.0, 'nS': 1.0}, 'AGA': {'aa': 'R', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'AGC': {'aa': 'S', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'AGG': {'aa': 'R', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'AGT': {'aa': 'S', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'ATA': {'aa': 'I', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'ATC': {'aa': 'I', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'ATG': {'aa': 'M', 'nN': 3.0, 'nS': 0.0}, 'ATT': {'aa': 'I', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'CAA': {'aa': 'Q', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CAC': {'aa': 'H', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CAG': {'aa': 'Q', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CAT': {'aa': 'H', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'CCA': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CCC': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CCG': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CCT': {'aa': 'P', 'nN': 2.0, 'nS': 1.0}, 'CGA': {'aa': 'R', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CGC': {'aa': 'R', 'nN': 2.0, 'nS': 1.0}, 'CGG': {'aa': 'R', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CGT': {'aa': 'R', 'nN': 2.0, 'nS': 1.0}, 'CTA': {'aa': 'L', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CTC': {'aa': 'L', 'nN': 2.0, 'nS': 1.0}, 'CTG': {'aa': 'L', 'nN': 1.6666666666666667, 'nS': 1.3333333333333333}, 'CTT': {'aa': 'L', 'nN': 2.0, 'nS': 1.0}, 'GAA': {'aa': 'E', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GAC': {'aa': 'D', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GAG': {'aa': 'E', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GAT': {'aa': 'D', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'GCA': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GCC': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GCG': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GCT': {'aa': 'A', 'nN': 2.0, 'nS': 1.0}, 'GGA': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GGC': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GGG': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GGT': {'aa': 'G', 'nN': 2.0, 'nS': 1.0}, 'GTA': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'GTC': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'GTG': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'GTT': {'aa': 'V', 'nN': 2.0, 'nS': 1.0}, 'TAC': {'aa': 'Y', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TAT': {'aa': 'Y', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TCA': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TCC': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TCG': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TCT': {'aa': 'S', 'nN': 2.0, 'nS': 1.0}, 'TGC': {'aa': 'C', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TGG': {'aa': 'W', 'nN': 3.0, 'nS': 0.0}, 'TGT': {'aa': 'C', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TTA': {'aa': 'L', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'TTC': {'aa': 'F', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}, 'TTG': {'aa': 'L', 'nN': 2.3333333333333335, 'nS': 0.6666666666666666}, 'TTT': {'aa': 'F', 'nN': 2.6666666666666665, 'nS': 0.3333333333333333}}


def fasta2dic(fastaFile):
	fasta = open(fastaFile).readlines()
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


L = len(alignA[alignA.keys()[0]]) - 3 # remove the last 3 bases to excluse final stop codon
L = trunc2triplets(L) # convert the remaining length into a multiple of 3


interspe = []
nA = 0
nB = 0
for i in alignA:
	seq = alignA[i][0:L]
	propN = seq.count("N")/L
	if propN < threshold_N:
		nA += 1
		interspe.append(seq)


for i in alignB:
	seq = alignB[i][0:L]
	propN = seq.count("N")/L
	if propN < threshold_N:
		nB += 1
		interspe.append(seq)


if nA < nMin:
	if nB >= nMin:
		sys.exit("Species A has less than {0} sequences. Species B is ok".format(nMin))
	if nB < nMin:
		sys.exit("Species A and B have less than {0} sequences".format(nMin))
if nB < nMin:
	sys.exit("Species B have less than {0} sequences. Species A is ok". format(nMin))


#for i in interspe:
#	print(">{0}\n{1}".format(i, coloredSeq(interspe[i][0:120])))



nSites = 0 # total number of synonymous sites within the sequence, computed using codonTable
positions = [] # list of synonymous polymorphic positions
msStyle = [] # contains the msStyle format
for ind in range(nA):
	msStyle.append([])
for ind in range(nB):
	msStyle.append([])

# loop over codons:
for pos in range(L)[::3]:
	alignmentOfCodons = [] # set of codons in the alignment, starting at the position 'pos1'
	# loop over individuals:
	# get all codons in the alignment
	for ind in range(nA + nB):
		pos1 = interspe[ind][pos]
		pos2 = interspe[ind][pos + 1]
		pos3 = interspe[ind][pos + 2]
		base = pos1 + pos2 + pos3 
		alignmentOfCodons.append(base)
	polyMcodons = list(set(alignmentOfCodons)) # list of codons found in the alignment
	nCodons = 0
	nCodons = len(polyMcodons)
	testN = False # False if no codon with 'N'; True if a 'N' is found in at least one codon for one individual
	testStopCodon = False # False if no stop codon was found; True if a stop codon was found
	for i in polyMcodons: # loop to test for some 'N'
		if 'N' in i:
			testN = True
		if i not in codonTable:
			testStopCodon = True
	# if: 1) a maximum of 2 polymorphic codons, and, 2) no codon with 'N', and, 3) all codons effectively code for an amino acid
	if nCodons <= 2 and testN==False and testStopCodon==False: 
		nSites_pos = 0.0
		for i in alignmentOfCodons:
			nSites_pos += codonTable[i]['nS']
		nSites += nSites_pos/len(alignmentOfCodons)
		# if two codons --> there is a polymorphism
		if nCodons == 2:
			alignmentOfAminoAcids = []
			for i in alignmentOfCodons:
				alignmentOfAminoAcids.append(codonTable[i]['aa'])
			setOfAminoAcids = list(set(alignmentOfAminoAcids))
			# if two codons but one amino acids --> synonymous polymorphism
			if len(setOfAminoAcids) == 1: 
				positions.append(pos)
				ancestralAllele = polyMcodons[0] # in absence of outgroup --> the ancestral allele is the first in the alignement
				derivedAllele = polyMcodons[1] # without outgroup --> the derived allele is the one who is not the first...
				for i in range(nA + nB):
					if alignmentOfCodons[i] == ancestralAllele:
						msStyle[i].append('0')
					if alignmentOfCodons[i] == derivedAllele:
						msStyle[i].append('1')


print("# of synonymous positions = {0}".format(nSites))
print("\t".join( [ str(i) for i in positions ] ))